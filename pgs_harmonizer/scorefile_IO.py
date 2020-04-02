import pandas as pd


remap_header = {
    'PGS ID' : 'pgs_id',
    'Reported Trait' : 'trait_reported',
    'Original Genome Build' : 'genome_build',
    'Number of Variants' : 'variants_number',
    'PGP ID' : 'pgp_id',
    'Citation' : 'citation'
} # ToDo remove once Scoring File headers are fixed


def read_ScoreFile(loc_scorefile):
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
