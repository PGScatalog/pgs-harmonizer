import argparse
import os
import ftplib


vcf_chrs = [str(x) for x in range(1,23)] + ['X', 'Y', 'MT']

parser = argparse.ArgumentParser(
    description='Download VCFs from ENSEMBL.')
parser.add_argument(dest="genome_build",
                    help="Target genome build choices: GRCh37, GRCh38, or both (GRCh3*)", metavar="GRCh3#",
                    choices=['GRCh37', 'GRCh38', 'GRCh3*'])
parser.add_argument(dest="loc_outputs",
                    help="Directory where the VCF files will be saved (default: PGS_HmPOS/)",
                    metavar="OUTDIR")
parser.add_argument('-chrom', dest='single_chrom',
                    help='Option to download a single chrom, otherwise all chromosomes will be downloaded',
                    metavar="CHR", choices=vcf_chrs, required=False)
args = parser.parse_args()

# Check/create structure to hold VCF files
if os.path.isdir(args.loc_outputs) is False:
    os.mkdir(args.loc_outputs)

if args.genome_build == 'GRCh3*':
    vcf_builds = ['GRCh37', 'GRCh38']
else:
    vcf_builds = [args.genome_build]

for build in vcf_builds:
    if os.path.isdir('map/vcf_ref/{}'.format(build)) is False:
        os.mkdir('map/vcf_ref/{}'.format(build))
        os.mkdir('map/vcf_ref/{}/cohort_ref'.format(build))

# Check which chromosomes to download
if args.single_chrom is not None:
    dl_chrs = [args.single_chrom]
else:
    dl_chrs = vcf_chrs

# Download VCF files from ENSEMBL FTP
ensembl_ftp = ftplib.FTP("ftp.ensembl.org")
ensembl_ftp.login()

for build in vcf_builds:
    vcf_loc = 'pub/current_variation/vcf/homo_sapiens/'
    if build == 'GRCh37':
        vcf_loc = 'pub/grch37/current/variation/vcf/homo_sapiens/'

    for chr in dl_chrs:
        print('Downloading (vcf):', build, 'chr{}'.format(chr))
        try:
            # Download vcf
            fn_vcf = 'homo_sapiens-chr{}.vcf.gz'.format(chr)
            vcf_loc_chr = vcf_loc + fn_vcf
            with open("./map/vcf_ref/{}/{}".format(build, fn_vcf), 'wb') as lf:
                ensembl_ftp.retrbinary('RETR {}'.format(vcf_loc_chr), lf.write)
            # Download index .csi / .tbi
            fn_vcf = 'homo_sapiens-chr{}.vcf.gz.csi'.format(chr)
            vcf_loc_chr = vcf_loc + fn_vcf
            with open("./map/vcf_ref/{}/{}".format(build, fn_vcf), 'wb') as lf:
                ensembl_ftp.retrbinary('RETR {}'.format(vcf_loc_chr), lf.write)
        except:
            print('! ERROR downloading:', build, 'chr{}'.format(chr))

    # ToDo: Download variant_location_37/38 from PGS Catalog FTP (to be distributed/updated each release)

ensembl_ftp.quit()
