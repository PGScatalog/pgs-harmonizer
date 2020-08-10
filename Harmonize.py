import argparse
from pgs_harmonizer.ensembl_tools import ensembl_post, clean_rsIDs, parse_var2location
from pgs_harmonizer.harmonize import *
from pgs_harmonizer.liftover_tools import liftover
from pgs_harmonizer.vcf_tools import *
import os, sys, gzip
from collections import Counter

# Inputs

parser = argparse.ArgumentParser(description='Harmonize a PGS Catalog Scoring file (PGS######.txt.gz) to a specific genome build.')
parser.add_argument("-id", dest="pgs_id",
                    help="PGS Catalog Score ID", metavar="PGSID", required=True)
parser.add_argument("-build", dest="target_build",
                    help="Target genome build choices = ['GRCh37', GRCh38']", metavar="GENOMEBUILD",
                    choices=['GRCh37', 'GRCh38'], required=True)
parser.add_argument("-loc_scorefiles", dest="loc_scorefiles",
                    help="Root directory where the PGS files are located, otherwise assumed to be in: ../pgs_ScoringFiles/", metavar="DIR",
                    default='../pgs_ScoringFiles/')
parser.add_argument("-source_build", dest="source_build",
                    help="Source genome build [overwrites the scoring file header information]", metavar="BUILD",
                    default=None)
parser.add_argument('--var2location', action='store_true')
args = parser.parse_args()

# Locations of relevant files
if 'loc_scorefiles' in args:
    if not args.loc_scorefiles.endswith('/'):
        args.loc_scorefiles += '/'
    loc_scorefile = args.loc_scorefiles + args.pgs_id + '.txt.gz'
else:
    loc_scorefile = '../pgs_ScoringFiles/{}.txt.gz'.format(args.pgs_id)

# Read Score File
print('Reading Score File')
header, df_scoring  = read_scorefile(loc_scorefile)

source_build = header['genome_build']
if args.source_build is not None:
    source_build = args.source_build

print('PGS ID: {} | Build: {}'.format(header['pgs_id'], source_build))
print('Number of variants (score file lines) = {}'.format(header['variants_number']))
# ToDo - print columns available for mapping

# Sorting out the genome build
mappable = False
if source_build == None:
    if 'rsID' in df_scoring.columns:
        mappable = True
    else:
        print('Log: Need to guess the source genome build')
        # ToDo guess the genome build
else:
    mappable = True

if mappable:
    # Load Liftover Chains
    if source_build is not None:
        build_map = liftover(source_build, args.target_build)  # Get the chain file

    # Source ENSEMBL DB/API variant mappings if required
    print('Mapping -> {}'.format(args.target_build))
    if 'rsID' in df_scoring.columns:
        tomap_rsIDs = clean_rsIDs(list(df_scoring['rsID']))
        if args.var2location:
            # Write list of rsIDs that need mapping via ENSEMBL
            with open('EnsemblMappings/variants/{}.txt'.format(header['pgs_id']), 'w') as outf:
                outf.write('\n'.join(tomap_rsIDs))

            # ToDo add command to run var2location.pl on the EBI cluster (using local ENSEMBL mirror)

            # Load ENSEMBL mappings
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

    print('Starting Mapping')
    mapped_counter = Counter()
    #with gzip.open('./hm_coords/{}_hm{}.txt.gz'.format(header['pgs_id'], args.target_build), 'wt') as hm_out:
    with open('./hm_coords/{}_hm{}.txt'.format(header['pgs_id'], args.target_build), 'w') as hm_out:
        # Create harmonization/output formatting object
        hm_formatter = Harmonizer(df_scoring.columns)
        hm_out.write('\t'.join(hm_formatter.cols_order) + '\n')

        #Loop through variants
        for i, v in df_scoring.iterrows():
            current_rsID = None

            hm = (None, None, None)
            if mapping_ensembl:
                v_map = mapping_ensembl.get(v['rsID'])
                if v_map is not None:
                    hm = list(v_map.select_canonical_data(chromosomes))
                    current_rsID = v_map.id
                    mapped_counter['mapped_rsID'] += 1
                    ## Check if this  matches the ENSEMBL annotation
                    ## DEPRECATED so that everything is evaluated in the same way - also ENSEMBL allele strings
                    ## are encoded different than the VCF line (e.g encoding of indels)
                    # if 'reference_allele' in v:
                    #     hm[:3] = v_map.check_alleles(ref=v['reference_allele'], eff=v['effect_allele'])
                    # else:
                    #     hm[:3] = v_map.check_alleles(eff=v['effect_allele'])
            elif 'chr_name' and 'chr_position' in df_scoring.columns:
                if build_map.chain:
                    hm = list(build_map.lift(v['chr_name'], v['chr_position'])) # Mapping by liftover
                    mapped_counter['mapped_lift'] += 1
            if all([x == None for x in hm]):
                mapped_counter['mapped_unable'] += 1

            # Check the VCF for variant
            # ToDo Revise harmonization codes
            if hm[2] is not None:
                v_records = VCFResult(chr=hm[0], pos=hm[1], build=args.target_build)
                if 'reference_allele' in v:
                    hm[2] = v_records.check_alleles(ref=v['reference_allele'], eff=v['effect_allele'])
                else:
                    hm[2] = v_records.check_alleles(eff=v['effect_allele'])

            hm_out.write('\t'.join(hm_formatter.format_line(v, hm, source_build, rsid=current_rsID)) + '\n')

            if (i % 250000 == 0) and (i != 0):
                print('Mapped {} / {} lines'.format(i, df_scoring.shape[0]))

    if mapped_counter['mapped_rsID'] > 0:
        print('{} lines mapped by rsID'.format(mapped_counter['mapped_rsID']))
    if mapped_counter['mapped_lift'] > 0:
        print('{} lines mapped by liftover'.format(mapped_counter['mapped_lift']))
    if mapped_counter['mapped_author'] > 0:
        print('{} lines used author-reported mappings'.format(mapped_counter['mapped_author']))
    print('Mapping complete')
