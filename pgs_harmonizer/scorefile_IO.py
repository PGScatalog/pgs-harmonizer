import pandas as pd

remap_header = {
    'PGS ID' : 'pgs_id',
    'Reported Trait' : 'trait_reported',
    'Original Genome Build' : 'genome_build',
    'Number of Variants' : 'variants_number',
    'PGP ID' : 'pgp_id',
    'Citation' : 'citation'
} # ToDo remove once Scoring File headers are fixed

def read_scorefile(loc_scorefile):
    if loc_scorefile.endswith('.gz'):
        import gzip
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

    return(header, df_scoring)

class WriteHarmonized:
    '''Class to choose which columns will be added, and in what order to the output.'''
    def __init__(self, cols):
        '''Used to select the columns and ordering of the output file'''
        self.cols_previous = cols
        self.hm_fields = ['rsID', 'chr_name', 'effect_allele', 'reference_allele']

        self.cols_order = ['chr_name', 'chr_position']

        # Check if the rsID will be added
        if 'rsID' in self.cols_previous:
            self.cols_order.append('rsID')

        # Check if there is a reference allele will be added
        if 'reference_allele' in self.cols_previous:
            self.cols_order += ['effect_allele', 'reference_allele']
        else:
            self.cols_order.append('effect_allele')

        # Check which other columns need to be added
        for c in self.cols_previous:
            if c not in self.cols_order:
                self.cols_order.append(c)

        self.cols_order += ['hm_code', 'hm_info']

    def format_line(self, v, hm, build):
        v = dict(v)

        hm_info = {}
        if hm[2] is None:
            hm[2] = '-'
            hm_info['original_build'] = build
            for c in self.hm_fields:
                f = v.get(c)
                if f:
                    hm_info[c] = f
                    v[c] = '' # insert blank, or '-' ?
        else:
            v['chr_name'] = hm[0]
            v['chr_position'] = hm[1]
            v['hm_code'] = hm[2]
            if hm[2] == 1: #mapped by rsID
                if hm[3] != v['rsID']:
                    hm_info['previous_rsID'] = v['rsID']
                    v['rsID'] = hm[3]

        # Create output
        o = []
        for c in self.cols_order:
            if c =='hm_info':
                if len(hm_info) > 0:
                    o.append(repr(hm_info))
                else:
                    o.append('')
            elif c in v:
                o.append(str(v[c]))
            else:
                o.append('')
        return o




