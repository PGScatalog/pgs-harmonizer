import argparse, os, sys, gzip
from tqdm import tqdm
from pgs_harmonizer.harmonize import *
from datetime import datetime

# Inputs
parser = argparse.ArgumentParser(
    description='Harmonize a PGS Catalog Scoring file (PGS######.txt.gz) to a specific genome build.')
subparsers = parser.add_subparsers(help='Harmonization Commands help', dest='HmAction')

# Sub-parser for HmPOS
parser_POS = subparsers.add_parser('HmPOS', help='HmPOS - Harmonizing position information (adding/updating chr/pos information)')
parser_POS.add_argument(dest="pgs_id", help="PGS Catalog Score ID", metavar="PGS######", type=str)
parser_POS.add_argument(dest="target_build", help="Target genome build choices: 'GRCh37'or GRCh38'", metavar="GRCh3#",
                        choices=['GRCh37', 'GRCh38'])
parser_POS.add_argument("-loc_files", dest="loc_scorefiles",
                        help="Root directory where the PGS files are located, otherwise assumed to be in: ../pgs_ScoringFiles/",
                        metavar="DIR",
                        default='../pgs_ScoringFiles/', required=False)
parser_POS.add_argument("-source_build", dest="source_build",
                        help="Source genome build [overwrites information in the scoring file header]",
                        metavar="GENOMEBUILD",
                        default=None, required=False)
parser_POS.add_argument("-loc_hmoutput", dest="loc_outputs",
                        help="Directory where the harmonization output will be saved (default: PGS_HmPOS/)",
                        metavar="DIR",
                        default='./PGS_HmPOS/', required=False)
parser_POS.add_argument('--var2location',
                        help='Uses the annotations from the var2location.pl script (ENSEMBL SQL connection)',
                        action='store_true', required=False)
parser_POS.add_argument('--silent_tqdm', help='Disables tqdm progress bar',
                        action='store_true', required=False)
parser_POS.add_argument('--ignore_rsid', help='Ignores rsID mappings and harmonizes variants using only liftover',
                        action='store_true', required=False)
parser_POS.add_argument('--gzip', help='Writes gzipped harmonized output',
                        action='store_true', required=False)

# Sub-parser for HmVCF
parser_VCF = subparsers.add_parser('HmVCF', help='HmVCF - Checking positional information and/or adding other_alleles')
parser_VCF.add_argument(dest="pgs_id", help="PGS Catalog Score ID", metavar="PGS######", type=str)
parser_VCF.add_argument("-loc_files", dest="loc_scorefiles",
                        help="Root directory where the PGS files are located, otherwise assumed to be in: PGS_HmPOS/",
                        metavar="DIR",
                        default='./PGS_HmPOS/', required=False)
parser_VCF.add_argument(dest="target_build",
                        help="Target genome build choices: 'GRCh37'or GRCh38'",
                        metavar="GRCh3#",
                        choices=['GRCh37', 'GRCh38'])
parser_VCF.add_argument("-cohort_vcf", dest="cohort_name",
                        help="Cohort VCF: Used to check if a variant is present in the genotyped/imputed variants for "
                             "a cohort and add other allele when the information from ENSEMBL is ambiguous "
                             "(multiple potential alleles)",
                        metavar="COHORT",
                        default=None, required=False)
parser_VCF.add_argument('--addOtherAllele',
                        help='Adds a other_allele(s) column for PGS that only have a recorded effect_allele',
                        action='store_true', required=False)
parser_VCF.add_argument('--addVariantID',
                        help='Returns a column with the ID from the VCF corresponding to the match variant/allele(s)',
                        action='store_true', required=False)
parser_VCF.add_argument('--author_reported',
                        help='Replaces unmappable variants (hm_code = -5) with the author-reported code (hm_code = 0)',
                        action='store_true', required=False)
parser_VCF.add_argument('--skip_strandflips',
                        help='This flag will stop the harmonizing from trying to correct strand flips',
                        action='store_true', required=False)
parser_VCF.add_argument('--silent_tqdm', help='Disables tqdm progress bar',
                        action='store_true', required=False)
parser_VCF.add_argument('--gzip', help='Writes gzipped harmonized output',
                        action='store_true', required=False)

args = parser.parse_args()


def variant_HmPOS(v, rsIDmaps=None, liftchain=None, isSameBuild=False, inferOtherAllele=False):
    hm_source = ''  # {'Author-reported', 'ENSEMBL Variation', 'liftover' }
    hm_rsID = ''
    hm_chr = ''
    hm_pos = ''
    hm_inferOtherAllele = None  # Field to capture the inferred other/reference allele

    if rsIDmaps and (v['rsID'] in rsIDmaps):
        v_map = rsIDmaps.get(v['rsID'])
    else:
        v_map = None

    if v_map is not None:
        hm_chr, hm_pos, hm_alleles = list(v_map.select_canonical_data(chromosomes))
        hm_source = 'ENSEMBL'
        hm_rsID = v_map.id
        if inferOtherAllele:
            if 'other_allele' in v:
                if pd.isnull(v.get('other_allele')) is True:
                    hm_inferOtherAllele = v_map.infer_OtherAllele(v['effect_allele'])
            else:
                hm_inferOtherAllele = v_map.infer_OtherAllele(v['effect_allele'])  # Based on the rsID
    elif 'chr_name' and 'chr_position' in v:
        if isSameBuild:
            hm_chr = v['chr_name']
            hm_pos = v['chr_position']
            hm_source = 'Author-reported'  # Author-reported
        elif (liftchain is not None) and (liftchain.chain is not None):
            if (pd.isnull(v['chr_name']) is False) and (pd.isnull(v['chr_position']) is False):
                hm_chr, hm_pos, hm_liftover_multimaps = list(liftchain.lift(v['chr_name'], v['chr_position']))  # liftover
                hm_source = 'liftover'
    if all([x is None for x in [hm_chr, hm_pos]]):
        hm_source = 'Unknown'

    if hm_pos != '':
        hm_pos = str(hm_pos)

    if inferOtherAllele:
        if hm_inferOtherAllele is None:
            hm_inferOtherAllele = ''
        return pd.Series([hm_source, hm_rsID, hm_chr, hm_pos, hm_inferOtherAllele])
    else:
        return pd.Series([hm_source, hm_rsID, hm_chr, hm_pos])


def run_HmPOS(args):
    # Module-specifc imports
    from pgs_harmonizer.ensembl_tools import ensembl_post, clean_rsIDs, parse_var2location
    from pgs_harmonizer.liftover_tools import liftover, map_release

    ## Set I/O File focations
    # Scoring file location
    if 'loc_scorefiles' in args:
        if not args.loc_scorefiles.endswith('/'):
            args.loc_scorefiles += '/'
        loc_scorefile = args.loc_scorefiles + args.pgs_id + '.txt.gz'
    else:
        loc_scorefile = '../pgs_ScoringFiles/{}.txt.gz'.format(args.pgs_id)
    # Define output location
    ofolder = args.loc_outputs
    if ofolder.endswith('/'):
        ofolder = ofolder[:-1]
    if os.path.isdir(ofolder) is False:
        os.mkdir(ofolder)
    loc_hm_out = '{}/{}_hmPOS{}.txt'.format(ofolder, args.pgs_id, args.target_build)
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
    else:
        print('PGS ID: {} | Build: {}'.format(header['pgs_id'], source_build))
    print('Number of variants (score file lines) = {}'.format(header['variants_number']))
    # ToDo - print columns available for mapping

    # Sorting out the genome build
    mappable = False
    tf_unmappable2authorreported = False
    if source_build_mapped is None:
        if 'rsID' in df_scoring.columns:
            mappable = True
        else:
            sys.exit(
                '{} CAN NOT BE HARMONIZED: Need to specify the source genome build (-source_build)'.format(
                    header['pgs_id']))
            # ToDo possibly implement a method to guess the genome build using GRCh37/38 VCFs
    else:
        mappable = True

    if mappable is False:
        sys.exit(
            '{} CAN NOT BE HARMONIZED: Insufficient variant or genome build data for mapping'.format(header['pgs_id']))

    if args.target_build == source_build_mapped:
        print('Harmonizing -> {}'.format(args.target_build))
        isSameBuild = True
    else:
        print('Re-Mapping/Lifting + Harmonizing -> {}'.format(args.target_build))
        isSameBuild = False

    # Load Liftover Chains
    if source_build is not None:
        build_map = liftover(source_build, args.target_build)  # Get the chain file
    else:
        build_map = None

    # Source ENSEMBL DB/API variant mappings if required
    mapping_ensembl = None
    if 'rsID' in df_scoring.columns and args.ignore_rsid is False:
        tomap_rsIDs = clean_rsIDs(list(df_scoring['rsID']))
        if args.var2location:
            # Write list of rsIDs that need mapping via ENSEMBL
            with open('EnsemblMappings/variants/{}.txt'.format(header['pgs_id']), 'w') as outf:
                outf.write('\n'.join(tomap_rsIDs))

            # ToDo add command to run var2location.pl on the EBI cluster (using local ENSEMBL mirror) or loop through ENSEMBL VCF for mappings

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
    # Run Variant Harmonization
    tqdm.pandas()
    # df_scoring[['hm_source', 'hm_rsID', 'hm_chr', 'hm_pos']] = df_scoring.progress_apply(variant_HmPOS,
    #                                                                                      axis=1,
    #                                                                                      rsIDmaps=mapping_ensembl,
    #                                                                                      liftchain=build_map,
    #                                                                                      isSameBuild=isSameBuild,
    #                                                                                      inferOtherAllele=False)
    df_scoring[['hm_source', 'hm_rsID', 'hm_chr', 'hm_pos', 'hm_inferOtherAllele']] = df_scoring.progress_apply(variant_HmPOS,
                                                                                                                axis=1,
                                                                                                                rsIDmaps=mapping_ensembl,
                                                                                                                liftchain=build_map,
                                                                                                                isSameBuild=isSameBuild,
                                                                                                                inferOtherAllele=True)
    print('Mapped {} -> {}'.format(dict(df_scoring['hm_source'].value_counts()), loc_hm_out))
    # Append information to header:
    header.update({'HmPOS_build': args.target_build,
                   'HmPOS_date': str(datetime.date(datetime.now()))}) # ToDo Consider adding information about the ENSEMBL build?

    # Output Data
    if args.gzip is True:
        hm_out = gzip.open(loc_hm_out, 'wt')
    else:
        hm_out = open(loc_hm_out, 'w')
    hm_out.write('\n'.join(create_scoringfileheader(header)))
    hm_out.write('\n')
    df_scoring.to_csv(hm_out, mode='a', index=False, sep='\t')# Write output using pandas
    hm_out.close() # Close file


def variant_HmVCF(v, vcfs_targetbuild, CohortVCF=None, returnOtherAllele=True):
    hm_source = v['hm_source']  # {'Author-reported', 'ENSEMBL Variation', 'liftover' }
    hm_matchesVCF = False  # T/F whether the variant is consistent with the VCF/Variant Lookup
    hm_isPalindromic = False  # T/F whether the alleles are consistent with being palindromic
    hm_isFlipped = False  # T/F whether the alleles are consistent with the negative strand (from VCF)
    hm_vid = None
    hm_code = None  # Derived from the above True/False information

    # Sort out non-effect/other allele
    other_allele = None
    if 'other_allele' in v:
        if pd.isnull(v['other_allele']) is False:
            other_allele = v['other_allele']
    hm_inferOtherAllele = None  # Field to capture the inferred other/reference allele
    if 'hm_inferOtherAllele' in v:
        if pd.isnull(v['hm_inferOtherAllele']) is False:
            hm_inferOtherAllele = v['hm_inferOtherAllele']

    # Check/select allales
    if pd.isnull(v['hm_source']) is False:
        v_records = vcfs_targetbuild.vcf_lookup(chromosome=v['hm_chr'], position=v['hm_pos'], rsid=v['hm_rsID'])
        if CohortVCF is not None:
            hm_source += '+{}'.format(CohortVCF)

        if other_allele is None:
            if returnOtherAllele is True:
                other_allele, hm_TF, hm_vid, hm_code = v_records.infer_OtherAllele(eff=v['effect_allele'],
                                                                                          oa_ensembl=hm_inferOtherAllele)
                hm_matchesVCF, hm_isPalindromic, hm_isFlipped = hm_TF
            else:
                hm_TF, hm_vid = v_records.check_alleles(eff=v['effect_allele'])
                hm_matchesVCF, hm_isPalindromic, hm_isFlipped = hm_TF
        else:
            hm_TF, hm_vid = v_records.check_alleles(eff=v['effect_allele'],
                                                    oa=other_allele)
            hm_matchesVCF, hm_isPalindromic, hm_isFlipped = hm_TF

    if hm_code is None:
        if 'other_allele' in v:
            hm_code = DetermineHarmonizationCode(hm_matchesVCF, hm_isPalindromic, hm_isFlipped,
                                                 alleles=[v['effect_allele'], v['other_allele']])
        else:
            hm_code = DetermineHarmonizationCode(hm_matchesVCF, hm_isPalindromic, hm_isFlipped,
                                                 alleles=[v['effect_allele']])

    # ToDo handle INDEL lookups in VCFs (e.g. ENSEMBL) better
    # ToDo (use allele frequency to resolve ambiguous variants hm_code=3)

    #Output
    if returnOtherAllele is True:
        return pd.Series([hm_source, hm_vid, hm_code, other_allele])
    else:
        return pd.Series([hm_source, hm_vid, hm_code])


def run_HmVCF(args):
    from pgs_harmonizer.variantlookup_tools import VCFs

    ## Set I/O File focations
    # Scoring file location
    if 'loc_scorefiles' in args:
        if not args.loc_scorefiles.endswith('/'):
            args.loc_scorefiles += '/'
        loc_scorefile = args.loc_scorefiles + '{}_hmPOS{}.txt.gz'.format(args.pgs_id, args.target_build)
    else:
        loc_scorefile = 'PGS_HmPOS/{}_hmPOS{}.txt.gz'.format(args.pgs_id, args.target_build)
    try:
        print('Reading Score File')
        header, df_scoring = read_scorefile(loc_scorefile)
    except:
        print('There was an error opening the file!')
        raise IOError

    # Load Variant References (VCF & Cohort)
    usingCohortVCF = None
    if args.cohort_name is not None:
        vcfs_targetbuild = VCFs(build=args.target_build, cohort_name=args.cohort_name)
        usingCohortVCF = args.cohort_name
        args.addOtherAllele = True
    else:
        vcfs_targetbuild = VCFs(build=args.target_build)  # ENSEMBL VCF
    if (vcfs_targetbuild.VCF is None) and (len(vcfs_targetbuild.by_chr) == 0):
        print('ERROR: Could not find the VCF')
        raise IOError

    # Apply the harmonization to the df
    tqdm.pandas()
    if args.addOtherAllele is True:
        df_scoring[['hm_source', 'hm_vid', 'hm_code', 'other_allele']] = df_scoring.progress_apply(variant_HmVCF, axis=1,
                                                                                                   vcfs_targetbuild=vcfs_targetbuild,
                                                                                                   CohortVCF=usingCohortVCF,
                                                                                                   returnOtherAllele=True)
    else:
        df_scoring[['hm_source', 'hm_vid', 'hm_code']] = df_scoring.progress_apply(variant_HmVCF, axis=1,
                                                                                   vcfs_targetbuild=vcfs_targetbuild,
                                                                                   CohortVCF=usingCohortVCF,
                                                                                   returnOtherAllele=False)
    # Post-Harmonization Fixes
    if args.skip_strandflips is False:
        df_scoring = FixStrandFlips(df_scoring) # also returns new column 'hm_fixedStrandFlip'

    # ToDo unmappable2authorreported
    # if args.author_reported is True:
    #     df_scoring = unmappable2authorreported(df_scoring)

    # Summarize Hm_Codes
    print(df_scoring['hm_code'].value_counts())

    # Format Output
    header['HmVCF_date'] = str(datetime.date(datetime.now()))
    if usingCohortVCF is not None:
        header['HmVCF_ref'] = usingCohortVCF
    else:
        header['HmVCF_ref'] = 'Ensembl Variation / dbSNP'#ToDo Consider adding information about the ENSEMBL build?

    # ToDo Write Output
    print(df_scoring.head())
    print(df_scoring.tail())
    hm_formatter = Harmonizer(df_scoring.columns, returnVariantID=args.addVariantID)
    hm_df = df_scoring.progress_apply(hm_formatter.format_line, axis=1,
                                      original_build=header['genome_build'])
    hm_df.columns = hm_formatter.cols_order
    hm_df.to_csv('TestHMVCF.csv', index=False)
    return


if __name__ == "__main__":
    if args.HmAction == 'HmPOS':
        run_HmPOS(args)
    elif args.HmAction == 'HmVCF':
        run_HmVCF(args)
    else:
        print('Not a valid method, try running: `python Harmonize.py -h` for valid options and more details')