import glob
from pgs_harmonizer.harmonize import *
from pgs_harmonizer.ensembl_tools import clean_rsIDs

allIDs = set()
rs_count = 0

for loc_scorefile in glob.glob('../pgs_ScoringFiles/*.txt.gz'):
    # Read Score File
    header, df_scoring = read_scorefile(loc_scorefile)
    if 'rsID' in df_scoring.columns:
        tomap_rsIDs = clean_rsIDs(list(df_scoring['rsID']))
        rs_count += len(tomap_rsIDs)
        allIDs.update(tomap_rsIDs)
        # Write list of rsIDs that need mapping via ENSEMBL
        print('{}: {} rsIDs'.format(header['pgs_id'], len(tomap_rsIDs)))
        with open('EnsemblMappings/variants/{}.txt'.format(header['pgs_id']), 'w') as outf:
            outf.write('\n'.join(tomap_rsIDs))

print('All PGS: {} unique rsIDs (count={})'.format(len(allIDs), rs_count))
with open('EnsemblMappings/variants/UNION.txt', 'w') as outf:
    allIDs = list(allIDs)
    allIDs.sort()
    outf.write('\n'.join(allIDs))