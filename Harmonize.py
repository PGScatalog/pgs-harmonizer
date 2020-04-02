from pgs_harmonizer.ensembl_tools import ensembl_post
from pgs_harmonizer.scorefile_IO import read_ScoreFile
import

#Globals
chromosomes = [str(x) for x in range(1,23)] + ['X', 'Y']


build_target = 'GRCh38'

loc_scorefile = '../pgs_ScoringFiles/PGS000001.txt.gz'

header, df_scoring  = read_ScoreFile(loc_scorefile)

# Sorting out the genome build
mappable = False
if header['genome_build'] == None:
    if 'rsID' in df_scoring.columns:
        if all([x.startswith('rs') for x in df_scoring['rsID']]):
            mappable = True #By rsID
    # ToDo guess the genome build
else:
    mappable = True

# Getting the Chromosome positions
if header['genome_build'] != build_target:
    df_scoring['hm_chr'] = ''
    df_scoring['hm_pos'] = int()

    # First try and map everything by rsID (if available)
    if 'rsID' in df_scoring.columns:
        mapping_ensembl = ensembl_post(list(df_scoring['rsID']), build_target) #retireve the SNP info from ENSEMBL
        for i, v in df_scoring.iterrows():
            v_map = mapping_ensembl.get(v['rsID'])
            if v_map is not None:
                df_scoring.loc[i,['hm_chr', 'hm_pos']] = v_map.select_canonical_data(chromosomes)

    # Check if everything has a chr/pos




