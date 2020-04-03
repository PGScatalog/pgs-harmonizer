from pgs_harmonizer.ensembl_tools import ensembl_post
from pgs_harmonizer.scorefile_IO import read_ScoreFile
from pgs_harmonizer.liftover_tools import liftover

#Globals
chromosomes = [str(x) for x in range(1,23)] + ['X', 'Y']


target_build = 'GRCh38'

loc_scorefile = '../pgs_ScoringFiles/PGS000111.txt.gz'

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
df_scoring['hm_chr'] = None
df_scoring['hm_pos'] = int()
df_scoring['hm_flag'] = None

if 'chr_name' or 'chr_position' not in df_scoring.columns:
    # First try and map everything by rsID (if available)
    if 'rsID' in df_scoring.columns:
        mapping_ensembl = ensembl_post(list(df_scoring['rsID']), target_build) #retireve the SNP info from ENSEMBL
        for i, v in df_scoring.iterrows():
            v_map = mapping_ensembl.get(v['rsID'])
            if v_map is not None:
                df_scoring.loc[i,['hm_chr', 'hm_pos', 'hm_flag']] = v_map.select_canonical_data(chromosomes)

    # Check if everything has a chr/pos
    if None in df_scoring['hm_chr'].values:
        # Check if it can be lifted over
        if (header['genome_build'] != target_build) and ('chr_name' and 'chr_position' in df_scoring.columns):
            build_map = liftover(header['genome_build'], target_build) # Get the chain file
            for i, v in df_scoring.iterrows():
                if v['hm_chr'] == None:
                    df_scoring.loc[i,['hm_chr', 'hm_pos', 'hm_flag']] = build_map.lift(v['chr_name'], v['chr_position'])




