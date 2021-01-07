import pandas as pd
import re
import json
import gzip

remap_header = {
    'PGS ID' : 'pgs_id',
    'Reported Trait' : 'trait_reported',
    'Original Genome Build' : 'genome_build',
    'Number of Variants' : 'variants_number',
    'PGP ID' : 'pgp_id',
    'Citation' : 'citation'
} # ToDo remove once Scoring File headers are fixed


chromosomes = [str(x) for x in range(1,23)] + ['X', 'Y', 'MT']
acceptable_alleles = re.compile('[ACGT]*$') #alleles that can be reverse complemented

def reversecomplement(x):
    if acceptable_alleles.match(x):
        return ''.join([{'A':'T','C':'G','G':'C','T':'A'}[B] for B in x][::-1])
    else:
        return None

#reversecomplement = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A'}[B] for B in x][::-1])

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

    df_scoring = pd.read_table(loc_scorefile, float_precision='round_trip', comment='#')

    # Make sure certain columns maintin specific datatypes
    if 'chr_name' in df_scoring.columns:
        df_scoring['chr_name'] = df_scoring['chr_name'].astype('str') # Asssert character in case there are only chr numbers
    if 'chr_position' in df_scoring.columns:
        df_scoring['chr_position'] = df_scoring['chr_position'].astype('int')  # Asssert int

    return(header, df_scoring)

class Harmonizer:
    """Class to select and harmonize variant locations in a PGS Scoring file."""
    def __init__(self, cols):
        '''Used to select the columns and ordering of the output PGS Scoring file'''
        self.cols_previous = cols
        self.hm_fields = ['rsID', 'chr_name', 'chr_position', 'effect_allele', 'reference_allele']

        self.cols_order = ['chr_name', 'chr_position']

        # Check if the rsID will be added
        if 'rsID' in self.cols_previous:
            self.cols_order.append('rsID')

        # Check if there is a reference allele will be added
        self.cols_order.append('effect_allele')
        if 'reference_allele' in self.cols_previous:
            self.cols_order.append('reference_allele')

        # Check which other columns need to be added
        for c in self.cols_previous:
            if c not in self.cols_order:
                self.cols_order.append(c)

        self.cols_order += ['hm_code', 'hm_info']

    def format_line(self, v, hm, build, rsid=None):
        """Method that takes harmonized variant location and compares it with the old information.
        Outputs any changes to an hm_info dictionary"""
        if type(hm) == tuple:
            hm = list(hm)
        v = dict(v)

        hm_info = {}
        if (hm[2] is None) or (hm[2] < 0):
            if hm[2] is None:
                v['hm_code'] = '-1' # Unable to map the variant
            else:
                v['hm_code'] = hm[2] # Variant does not match ENSEMBL Variation

            # Fill in the original mapping in 'hm_info' and remove the harmonized variant information from the column
            hm_info['original_build'] = build
            for c in self.hm_fields:
                f = v.get(c)
                if f:
                    hm_info[c] = f
                    v[c] = ''

            # If the variant was lifted but the allele's don't match add the lifted over positions
            if hm[0] != None:
                hm_info['hm_chr'] = hm[0]
            if hm[1] != None:
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
                o.append(str(v[c]).replace('nan', ''))
            else:
                o.append('')
        return o




