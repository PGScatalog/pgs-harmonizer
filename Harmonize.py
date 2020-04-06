from pgs_harmonizer.ensembl_tools import ensembl_post
from pgs_harmonizer.scorefile_IO import read_ScoreFile
from pgs_harmonizer.liftover_tools import liftover
from collections import Counter
import gzip

#Globals
chromosomes = [str(x) for x in range(1,23)] + ['X', 'Y']

target_build = 'GRCh38'

loc_scorefile = '../pgs_ScoringFiles/PGS000027.txt.gz'

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
        print('Log: Need to guess the source genome build')
    # ToDo guess the genome build
else:
    mappable = True

if mappable:
    print('Mapping -> {}'.format(target_build))
    if header['genome_build'] == target_build:
        print('Using author-reported variant annotations')
        df_scoring['hm_chr'] = df_scoring['chr_name']
        df_scoring['hm_pos'] = df_scoring['chr_position']
        df_scoring['hm_flag'] = 'Author-reported'
        df_scoring.to_csv('{}_hm.txt.gz'.format(header['pgs_id']), compression='gzip', sep='\t', index=False)
    else:

        if 'rsID' in df_scoring.columns:
            print('Retrieving rsID mappings from ENSEMBL API')
            mapping_ensembl = ensembl_post(list(df_scoring['rsID']), target_build) #retireve the SNP info from ENSEMBL
        else:
            mapping_ensembl = None

        if 'chr_name' and 'chr_position' in df_scoring.columns:
            build_map = liftover(header['genome_build'], target_build) # Get the chain file
            print('Retrieved Liftover chain: {} -> {} ({})'.format(header['genome_build'], target_build, build_map.chain_name))

        print('Starting Mapping')
        mapped_rsID = 0
        mapped_lift = 0
        with gzip.open('./hm_coords/{}_hm.txt.gz'.format(header['pgs_id']), 'wt') as hm_out:
            cols = list(df_scoring.columns) + ['hm_chr', 'hm_pos', 'hm_flag']
            hm_out.write('\t'.join(cols) + '\n')
            for i, v in df_scoring.iterrows():

                hm = [None, None, None]

                # Try to map by rsID
                if mapping_ensembl:
                    v_map = mapping_ensembl.get(v['rsID'])
                else:
                    v_map = None
                if v_map:
                    hm = list(v_map.select_canonical_data(chromosomes))
                    mapped_rsID += 1
                elif 'chr_name' and 'chr_position' in df_scoring.columns:
                    hm = list(build_map.lift(v['chr_name'], v['chr_position'])) # Mapping by liftover
                    mapped_lift += 1

                hm_out.write('\t'.join(map(str, list(v) + hm)) + '\n')

                if i%100000 == 0:
                    print('Mapped {} / {} lines'.format(i, df_scoring.shape[0]))
print('Mapping complete')
print('{} lines mapped by rsID'.format(mapped_rsID))
print('{} lines mapped by liftover'.format(mapped_lift))

# ##### Strand counting
# strand_counter = Counter()
# for i, v in df_scoring.iterrows():
