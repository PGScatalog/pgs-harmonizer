from pgs_harmonizer.api_tools import ensembl_post
import pandas as pd

loc_scorefile = '../pgs_ScoringFiles/PGS000001.txt.gz'
df_scoring = pd.read_table(loc_scorefile, float_precision='round_trip', comment = '#')

#Get mappings for 38
mapping_ensembl = ensembl_post(list(df_scoring['rsID']))