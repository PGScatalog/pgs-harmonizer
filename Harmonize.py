import argparse, os, sys, gzip
from collections import Counter
from pgs_harmonizer.harmonize import *
from pgs_harmonizer.ensembl_tools import ensembl_post, clean_rsIDs, parse_var2location
from pgs_harmonizer.liftover_tools import liftover, map_release
from pgs_harmonizer.variantlookup_tools import VCFs, load_cohortvariants

# Inputs
parser = argparse.ArgumentParser(
    description='Harmonize a PGS Catalog Scoring file (PGS######.txt.gz) to a specific genome build.')
parser.add_argument("-id", dest="pgs_id", help="PGS Catalog Score ID", metavar="PGS######", required=True)
parser.add_argument("-build", dest="target_build",
                    help="Target genome build choices: 'GRCh37'or GRCh38'", metavar="GRCh##",
                    choices=['GRCh37', 'GRCh38'], required=True)
parser.add_argument("-loc_scorefiles", dest="loc_scorefiles",
                    help="Root directory where the PGS files are located, otherwise assumed to be in: ../pgs_ScoringFiles/",
                    metavar="DIR",
                    default='../pgs_ScoringFiles/', required=False)
parser.add_argument("-source_build", dest="source_build",
                    help="Source genome build [overwrites information in the scoring file header]",
                    metavar="GENOMEBUILD",
                    default=None, required=False)
parser.add_argument('--var2location',
                    help='Uses the annotations from the var2location.pl script (ENSEMBL SQL connection)',
                    action='store_true', required=False)
parser.add_argument('--addReferenceAllele',
                    help='Adds a reference_allele(s) column for PGS that only have a recorded effect_allele',
                    action='store_true', required=False)
parser.add_argument('--ignore_rsid', help='Ignores rsID mappings and harmonizes variants using only liftover',
                    action='store_true', required=False)
parser.add_argument('--gzip', help='Writes gzipped harmonized output',
                    action='store_true', required=False)
args = parser.parse_args()

## I/O File focations
# Scoring file location
if 'loc_scorefiles' in args:
    if not args.loc_scorefiles.endswith('/'):
        args.loc_scorefiles += '/'
    loc_scorefile = args.loc_scorefiles + args.pgs_id + '.txt.gz'
else:
    loc_scorefile = '../pgs_ScoringFiles/{}.txt.gz'.format(args.pgs_id)
# Define output location
loc_hm_out = './hm_coords/{}_hm{}.txt'.format(args.pgs_id, args.target_build)
if args.gzip is True:
    loc_hm_out += '.gz'


# Read Score File
print('Reading Score File')
header, df_scoring = read_scorefile(loc_scorefile)


# Get consistent source build (e.g. NCBI/GRC)
source_build = header['genome_build']
source_build_mapped = None
if args.source_build is not None:
    source_build = args.source_build

if source_build in map_release.values():
    for grc, hg in map_release.items():
        if hg == source_build:
            source_build_mapped = grc
    print('PGS ID: {} | Build: {}/{}'.format(header['pgs_id'], source_build, source_build_mapped))
elif source_build in map_release.keys():
    source_build_mapped = source_build

    print('PGS ID: {} | Build: {}'.format(header['pgs_id'], source_build))
print('Number of variants (score file lines) = {}'.format(header['variants_number']))
# ToDo - print columns available for mapping

# Sorting out the genome build
mappable = False
if source_build_mapped is None:
    if 'rsID' in df_scoring.columns:
        mappable = True
    else:
        sys.exit(
            '{} CAN NOT BE HARMONIZED: Need to specify the source genome build (-source_build)'.format(header['pgs_id']))
        # ToDo possibly implement a method to guess the genome build using GRCh37/38 VCFs
else:
    mappable = True

if mappable is False:
    sys.exit('{} CAN NOT BE HARMONIZED: Insufficient variant or genome build data for mapping'.format(header['pgs_id']))

if args.target_build == source_build_mapped:
    print('Harmonizing -> {}'.format(args.target_build))
else:
    print('Re-Mapping/Liftingl + Harmonizing -> {}'.format(args.target_build))

# Load Liftover Chains
if source_build is not None:
    build_map = liftover(source_build, args.target_build)  # Get the chain file

# Source ENSEMBL DB/API variant mappings if required
mapping_ensembl = None
if 'rsID' in df_scoring.columns and args.ignore_rsid is False:
    tomap_rsIDs = clean_rsIDs(list(df_scoring['rsID']))
    if args.var2location:
        # Write list of rsIDs that need mapping via ENSEMBL
        with open('EnsemblMappings/variants/{}.txt'.format(header['pgs_id']), 'w') as outf:
            outf.write('\n'.join(tomap_rsIDs))

        # ToDo add command to run var2location.pl on the EBI cluster (using local ENSEMBL mirror)

        # Load ENSEMBL mappings
        loc_mapping = 'EnsemblMappings/{}/{}.out'.format(args.target_build, header['pgs_id'])
        loc_mapping_UNION = 'EnsemblMappings/{}/UNION.out'.format(args.target_build)
        if os.path.isfile(loc_mapping):
            print('Retrieving rsID mappings from ENSEMBL Mirror (var2location.pl:{}.out)'.format(header['pgs_id']))
            mapping_ensembl = parse_var2location(loc_mapping)
        elif os.path.isfile(loc_mapping_UNION):
            print('Retrieving rsID mappings from ENSEMBL Mirror (var2location.pl:UNION.out)')
            mapping_ensembl = parse_var2location(loc_mapping_UNION, tomap_rsIDs)
        else:
            sys.exit('Error: No rsID mappings from ENSEMBL Mirror (var2location.pl)')
    else:
        print('Retrieving rsID mappings from ENSEMBL API')
        mapping_ensembl = ensembl_post(tomap_rsIDs, args.target_build)  # Retrieve the SNP info from ENSEMBL

# Load Variant References (VCF & Cohort)
vcfs_targetbuild = VCFs(build=args.target_build)  # ENSEMBL VCF


print('Starting Mapping')
mapped_counter = Counter()
counter_hmcodes = Counter()

# Create harmonization/output formatting objects
if args.gzip is True:
    hm_out = gzip.open(loc_hm_out, 'wt')
else:
    hm_out = open(loc_hm_out, 'w')

# Check if extra columns (reference_allele) will need to be added
if args.addReferenceAllele is True and ('reference_allele' in df_scoring.columns is False):
    hm_formatter = Harmonizer(list(df_scoring.columns) + ['reference_allele'])
else:
    hm_formatter = Harmonizer(df_scoring.columns)
hm_out.write('\t'.join(hm_formatter.cols_order) + '\n')

#Loop through variants
for i, v in df_scoring.iterrows():
    # Variant harmonization information
    current_rsID = None

    hm_chr = None
    hm_pos = None

    hm_source = None  # { 0 : 'Author-reported', 1 : 'ENSEMBL Variation', 2 : 'liftover' }
    hm_alleles = []
    hm_matchesVCF = False
    hm_isPalindromic = False
    hm_isFlipped = False
    hm_OtherAllele = None

    hm_code = None  # Should be derived from the above information

    if mapping_ensembl and v['rsID'] in mapping_ensembl:
        v_map = mapping_ensembl.get(v['rsID'])
        if v_map is not None:
            hm_chr, hm_pos, hm_alleles = list(v_map.select_canonical_data(chromosomes))
            hm_source = 1
            current_rsID = v_map.id
            if 'reference_allele' not in v and args.addReferenceAllele:
                # Add reference allele(s) based on ENSMEBL - provided the variant matches the allele list
                hm_OtherAllele = v_map.infer_reference_allele(v['effect_allele'])
                v['reference_allele'] = hm_OtherAllele
            mapped_counter['mapped_rsID'] += 1
    elif 'chr_name' and 'chr_position' in df_scoring.columns:
        if source_build_mapped == args.target_build:
            hm_chr = v['chr_name']
            hm_pos = v['chr_position']
            hm_source = 0 # Author-reported
        elif build_map.chain:
            hm_chr, hm_pos, hm_liftcode = list(build_map.lift(v['chr_name'], v['chr_position']))  # Mapping by liftover
            mapped_counter['mapped_lift'] += 1
    if all([x == None for x in [hm_chr, hm_pos]]):
        mapped_counter['mapped_unable'] += 1

    # Check the VCF for variant
    # ToDo Revise harmonization codes
    if hm_source is not None:
        v_records = vcfs_targetbuild.vcf_lookup(chromosome=hm_chr, position=hm_pos)
        if ('reference_allele' in v) and (pd.isnull(v['reference_allele']) is False) and (('/' in v['reference_allele']) is False):
            hm_matchesVCF, hm_isPalindromic, hm_isFlipped = v_records.check_alleles(eff=v['effect_allele'],
                                                                                    ref=v['reference_allele'])
        else:
            hm_matchesVCF, hm_isPalindromic, hm_isFlipped = v_records.check_alleles(eff=v['effect_allele'])

        # If the variant does not work revert to author-reported if possible
        # if hm[2] < 0 and source_build_mapped == args.target_build:
        #     if 'chr_name' and 'chr_position' in df_scoring.columns:
        #         hm = [v['chr_name'], v['chr_position'], 0]  # author-reported

    hm_code = DetermineHarmonizationCode(hm_matchesVCF, hm_isPalindromic, hm_isFlipped, hm_source)
    hm = [hm_chr, hm_pos, hm_code]

    # Harmonize and write to file
    v_hm = hm_formatter.format_line(v, hm, hm_source, source_build, rsid=current_rsID)
    counter_hmcodes[v_hm[-2]] += 1
    hm_out.write('\t'.join(v_hm) + '\n')

    if (i % 250000 == 0) and (i != 0):
        print('Mapped {} / {} lines'.format(i, df_scoring.shape[0]))

hm_out.close()

if mapped_counter['mapped_rsID'] > 0:
    print('{} lines mapped by rsID'.format(mapped_counter['mapped_rsID']))
if mapped_counter['mapped_lift'] > 0:
    print('{} lines mapped by liftover'.format(mapped_counter['mapped_lift']))
if mapped_counter['mapped_author'] > 0:
    print('{} lines used author-reported mappings'.format(mapped_counter['mapped_author']))
print(counter_hmcodes)
print('Mapping complete')
