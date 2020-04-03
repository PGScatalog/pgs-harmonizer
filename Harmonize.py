from pgs_harmonizer.ensembl_tools import ensembl_post
from pgs_harmonizer.scorefile_IO import read_ScoreFile
from pgs_harmonizer.liftover_tools import liftover

#Globals
chromosomes = [str(x) for x in range(1,23)] + ['X', 'Y']


target_build = 'GRCh38'

loc_scorefile = '../pgs_ScoringFiles/PGS000035.txt.gz'

print('Reading Score File')
header, df_scoring  = read_ScoreFile(loc_scorefile)
print('PGS ID: {} | Build: {}'.format(header['pgs_id'], header['genome_build']))
print('Number of variants = {} (lines of score file)'.format(header['variants_number']))

# Sorting out the genome build
mappable = False
if header['genome_build'] == None:
    if 'rsID' in df_scoring.columns:
        if all([x.startswith('rs') for x in df_scoring['rsID']]):
            mappable = True #By rsID
    else:
        print('Log: Need to guess the source ge')
    # ToDo guess the genome build
else:
    mappable = True

if mappable:
    print('MAPPABLE')

# Getting Chromosome positions
df_scoring['hm_chr'] = None
df_scoring['hm_pos'] = int()
df_scoring['hm_flag'] = None

if header['genome_build'] == target_build:
    print('Using author-reported variant annotations')
    df_scoring['hm_chr'] = df_scoring['chr_name']
    df_scoring['hm_pos'] = df_scoring['chr_position']
    df_scoring['hm_flag'] = 'Author-reported'
else:
    if 'chr_name' or 'chr_position' not in df_scoring.columns:
        # First try and map everything by rsID (if available)
        if 'rsID' in df_scoring.columns:
            print('Retrieving rsID mappings from ENSEMBL API')
            mapping_ensembl = ensembl_post(list(df_scoring['rsID']), target_build) #retireve the SNP info from ENSEMBL

            print('Joining rsID mappings')
            map_count = 0
            for i, v in df_scoring.iterrows():
                v_map = mapping_ensembl.get(v['rsID'])
                if v_map is not None:
                    df_scoring.loc[i,['hm_chr', 'hm_pos', 'hm_flag']] = v_map.select_canonical_data(chromosomes)
                    map_count += 1
            print('{} coordinates obtained by rsID mapping'.format(map_count))

    # Check if everything has a chr/pos
    if None in df_scoring['hm_chr'].values:
        # Check if it can be lifted over
        if 'chr_name' and 'chr_position' in df_scoring.columns:
            build_map = liftover(header['genome_build'], target_build) # Get the chain file
            print('Starting Liftover: {} -> {} ({})'.format(header['genome_build'], target_build, build_map.chain_name))
            map_count = 0
            for i, v in df_scoring.iterrows():
                if v['hm_chr'] == None:
                    df_scoring.loc[i,['hm_chr', 'hm_pos', 'hm_flag']] = build_map.lift(v['chr_name'], v['chr_position'])
                    map_count += 1
                print('{} coordinates obtained by liftover'.format(map_count))

