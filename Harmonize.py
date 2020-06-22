import argparse
from pgs_harmonizer.ensembl_tools import ensembl_post, clean_rsIDs, parse_var2location
from pgs_harmonizer.scorefile_IO import read_scorefile, WriteHarmonized
from pgs_harmonizer.liftover_tools import liftover
import gzip
import os

parser = argparse.ArgumentParser(description='Harmonize a PGS Catalog Scoring file (PGS######.txt.gz) to a specific genome build.')
parser.add_argument("-id", dest="pgs_id",
                    help="PGS Catalog Score ID", metavar="PGSID", required=True)
parser.add_argument("-build", dest="target_build",
                    help="Target genome build choices = ['GRCh37', GRCh38']", metavar="GENOMEBUILD",
                    choices=['GRCh37', 'GRCh38'], required=True)
parser.add_argument("-loc_scorefiles", dest="loc_scorefiles",
                    help="Root directory where the PGS files are located, otherwise assumed to be in: ../pgs_ScoringFiles/", metavar="DIR",
                    default='../pgs_ScoringFiles/')
parser.add_argument('--var2location', action='store_true')


args = parser.parse_args()

# Globals
chromosomes = [str(x) for x in range(1,23)] + ['X', 'Y', 'MT']

# Inputs
if 'loc_scorefiles' in args:
    if not args.loc_scorefiles.endswith('/'):
        args.loc_scorefiles += '/'
    loc_scorefile = args.loc_scorefiles + args.pgs_id + '.txt.gz'
else:
    loc_scorefile = '../pgs_ScoringFiles/{}.txt.gz'.format(args.pgs_id)

print('Reading Score File')
header, df_scoring  = read_scorefile(loc_scorefile)
print('PGS ID: {} | Build: {}'.format(header['pgs_id'], header['genome_build']))
print('Number of variants (score file lines) = {}'.format(header['variants_number']))

# Sorting out the genome build
mappable = False
if header['genome_build'] == None:
    if 'rsID' in df_scoring.columns:
        mappable = True
    else:
        print('Log: Need to guess the source genome build')
    # ToDo guess the genome build
else:
    mappable = True

if mappable:
    print('Mapping -> {}'.format(args.target_build))
    if header['genome_build'] == args.target_build:
        print('Using author-reported variant annotations')
        df_scoring['hm_chr'] = df_scoring['chr_name']
        df_scoring['hm_pos'] = df_scoring['chr_position']
        df_scoring['hm_code'] = '0' # Author-reported
        df_scoring.to_csv('{}_hm.txt.gz'.format(header['pgs_id']), compression='gzip', sep='\t', index=False)
    else:
        if 'rsID' in df_scoring.columns:
            tomap_rsIDs = clean_rsIDs(list(df_scoring['rsID']))
            if args.var2location:
                with open('EnsemblMappings/variants/{}.txt'.format(header['pgs_id']), 'w') as outf:
                    outf.write('\n'.join(tomap_rsIDs))
                #ToDo add command to run var2location.pl on the EBI cluster (using local ENSEMBL mirror)
                loc_mapping = 'EnsemblMappings/{}/{}.out'.format(args.target_build, header['pgs_id'])
                if os.path.isfile(loc_mapping):
                    print('Retrieving rsID mappings from ENSEMBL Mirror (var2location.pl)')
                    mapping_ensembl = parse_var2location(loc_mapping)
                else:
                    sys.exit('Error: No rsID mappings from ENSEMBL Mirror (var2location.pl)')
            else:
                print('Retrieving rsID mappings from ENSEMBL API')
                mapping_ensembl = ensembl_post(tomap_rsIDs, args.target_build) #retireve the SNP info from ENSEMBL
        else:
            mapping_ensembl = None

        if 'chr_name' and 'chr_position' in df_scoring.columns:
            build_map = liftover(header['genome_build'], args.target_build) # Get the chain file
            print('Retrieved Liftover chain: {} -> {} ({})'.format(header['genome_build'], target_build, build_map.chain_name))

        print('Starting Mapping')
        mapped_rsID = 0
        mapped_lift = 0
        mapped_nope = 0
        with gzip.open('./hm_coords/{}_hm{}.txt.gz'.format(header['pgs_id'], args.target_build), 'wt') as hm_out:
            hm_formatter = WriteHarmonized(df_scoring.columns)
            hm_out.write('\t'.join(hm_formatter.cols_order) + '\n')
            for i, v in df_scoring.iterrows():
                if mapping_ensembl:
                    v_map = mapping_ensembl.get(v['rsID'])
                else:
                    v_map = None
                if v_map:
                    hm = v_map.select_canonical_data(chromosomes)
                    mapped_rsID += 1
                elif 'chr_name' and 'chr_position' in df_scoring.columns:
                    hm = build_map.lift(v['chr_name'], v['chr_position']) # Mapping by liftover
                    mapped_lift += 1
                else:
                    hm = (None, None, None)
                    mapped_nope += 1

                hm_out.write('\t'.join(hm_formatter.format_line(v,hm,header['genome_build'])) + '\n')

                if (i % 100000 == 0) and (i != 0):
                    print('Mapped {} / {} lines'.format(i, df_scoring.shape[0]))
print('Mapping complete')
if mapped_rsID > 0:
    print('{} lines mapped by rsID'.format(mapped_rsID))
if mapped_lift > 0:
    print('{} lines mapped by liftover'.format(mapped_lift))