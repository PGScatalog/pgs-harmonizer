import pandas as pd
import re
import json
import gzip

remap_header = {
    'PGS ID' : 'pgs_id',
    'PGS Name' : 'pgs_name',
    'Reported Trait' : 'trait_reported',
    'Original Genome Build' : 'genome_build',
    'Number of Variants' : 'variants_number',
    'PGP ID' : 'pgp_id',
    'Citation' : 'citation',
    'LICENSE' : 'pgs_license'
} # ToDo remove once Scoring File headers are fixed


chromosomes = [str(x) for x in range(1,23)] + ['X', 'Y', 'MT']
acceptable_alleles = re.compile('[ACGT]*$')  # alleles that can be reverse complemented

def reversecomplement(x):
    if acceptable_alleles.match(x):
        return ''.join([{'A':'T','C':'G','G':'C','T':'A'}[B] for B in x][::-1])
    else:
        return None

def conv2int(n):
    try:
        return int(n)
    except:
        return n

def read_scorefile(loc_scorefile):
    """Loads PGS Catalog Scoring file and parses the header into a dictionary"""
    if loc_scorefile.endswith('.gz'):
        f = gzip.open(loc_scorefile,'rt')
    else:
        f = open(loc_scorefile, 'rt')

    header = {}
    for line in f:
        line = line.strip()
        if line.startswith('#'):
            if '=' in line:
                line = line[1:].split('=')
                field, val = [x.strip() for x in line]
                header[remap_header[field]] = val # ToDo change once Scoring File headers are fixed
    f.close()

    if header['genome_build'] == 'NR':
        header['genome_build'] = None

    df_scoring = pd.read_table(loc_scorefile, float_precision='round_trip', comment='#',
                               dtype = {'chr_position' : 'object',
                                        'chr_name' : 'str'})

    # Make sure certain columns maintain specific datatypes
    if 'reference_allele' in df_scoring.columns:
        df_scoring = df_scoring.rename(columns={"reference_allele": "other_allele"})
    if 'chr_name' in df_scoring.columns:
        df_scoring['chr_name'] = df_scoring['chr_name'].astype('str') # Asssert character in case there are only chr numbers
    if 'chr_position' in df_scoring.columns:
        df_scoring.loc[df_scoring['chr_position'].isnull() == False, 'chr_position'] = [conv2int(x) for x in df_scoring.loc[df_scoring['chr_position'].isnull() == False, 'chr_position']]

    return(header, df_scoring)


def DetermineHarmonizationCode(hm_matchesVCF, hm_isPalindromic, hm_isFlipped,alleles = [], ):
    hm_coding = {
        (True, False, False): 5,
        (True, True, False): 4,
        (True, False, True): -4,
        (False, False, False): -5
    }
    hm_code = hm_coding.get((hm_matchesVCF, hm_isPalindromic, hm_isFlipped), -1)
    if hm_code == 4:
        if len(alleles) > 1:
            rc_alleles = [reversecomplement(x) for x in alleles]
            for A in alleles:
                if A in rc_alleles:
                    hm_code = 3
    return hm_code

class Harmonizer:
    """Class to select and harmonize variant locations in a PGS Scoring file."""
    def __init__(self, cols, ensureOtherAllele=False, returnVariantID = False):
        '''Used to select the columns and ordering of the output PGS Scoring file'''
        self.cols_previous = cols
        self.hm_fields = ['rsID', 'chr_name', 'chr_position', 'effect_allele', 'other_allele']

        if returnVariantID is True:
            self.cols_order = ['variant_id', 'chr_name', 'chr_position']
        else:
            self.cols_order = ['chr_name', 'chr_position']

        # Check if the rsID will be added
        if 'rsID' in self.cols_previous:
            self.cols_order.append('rsID')

        # Check if there is a reference allele will be added
        self.cols_order.append('effect_allele')
        if ('other_allele' in self.cols_previous) or (ensureOtherAllele is True):
            self.cols_order.append('other_allele')

        # Check which other columns need to be added
        for c in self.cols_previous:
            if c not in self.cols_order:
                self.cols_order.append(c)

        self.cols_order += ['hm_code', 'hm_info']

    def format_line(self, v, hm, hm_source, build, rsid=None,vcfid=None, fixflips=True):
        """Method that takes harmonized variant location and compares it with the old information.
        Outputs any changes to an hm_info dictionary"""
        if type(hm) == tuple:
            hm = list(hm)
        v = dict(v)
        v['variant_id'] = vcfid

        hm_info = {'hm_source' : hm_source}
        if (hm[2] is None) or (hm[2] < 0):
            if hm[2] is None:
                v['hm_code'] = -1 # Unable to map the variant
            else:
                v['hm_code'] = hm[2] # Variant does not match ENSEMBL Variation

            if (v['hm_code'] == -4) and (fixflips is True):
                hm_info['original_build'] = build
                hm_info['fixedStrandFlip'] = True
                for c in self.hm_fields:
                    f = v.get(c)
                    if f is not None:
                        if c in ['effect_allele', 'other_allele']:
                            flipped_allele = reversecomplement(f)
                            hm_info['original_{}'.format(c)] = f
                            v[c] = flipped_allele
                        elif c == 'rsID':
                            if rsid is not None:  # mapped by rsID
                                if rsid != v['rsID']:
                                    hm_info['previous_rsID'] = v['rsID']
                                    v['rsID'] = rsid
                        else:
                            hm_info[c] = f
                            v[c] = ''
                v['chr_name'] = hm[0]
                v['chr_position'] = hm[1]

            else:
                if (v['hm_code'] == -4) and (fixflips is False):
                    hm_info['fixedStrandFlip'] = False
                # Fill in the original mapping in 'hm_info' and remove the harmonized variant information from the column
                hm_info['original_build'] = build
                for c in self.hm_fields:
                    f = v.get(c)
                    if f:
                        hm_info[c] = f
                        v[c] = ''
                if ('variant_id' in self.cols_order) and (vcfid is not None):
                    hm_info['variant_id'] = vcfid
                    v['variant_id'] = ''

                # If the variant was lifted but the allele's don't match add the lifted over positions
                if hm[0] is not None:
                    hm_info['hm_chr'] = hm[0]
                if hm[1] is not None:
                    hm_info['hm_chr_position'] = hm[1]
        else:
            v['chr_name'] = hm[0]
            v['chr_position'] = hm[1]
            v['hm_code'] = hm[2]
            if rsid is not None: # mapped by rsID
                if rsid != v['rsID']:
                    hm_info['previous_rsID'] = v['rsID']
                    v['rsID'] = rsid

        # Create output
        o = []
        for c in self.cols_order:
            if c == 'hm_info':
                if len(hm_info) > 0:
                    o.append(json.dumps(hm_info))
                else:
                    o.append('')
            elif c in v:
                val = str(v[c]).replace('None', '')
                if val == 'nan':
                    val = ''
                o.append(val)
            else:
                o.append('')
        return o




