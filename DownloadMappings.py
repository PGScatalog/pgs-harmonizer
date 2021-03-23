import os
import ftplib

# Check/create structure to hold VCF files
if os.path.isdir('map/') is False:
    os.mkdir('map/')
if os.path.isdir('map/vcf_ref/') is False:
    os.mkdir('map/vcf_ref/')

vcf_builds = ['GRCh37', 'GRCh38']
vcf_chrs = list(range(1,23)) + ['X', 'Y', 'MT']

# Download VCF files from ENSEMBL FTP
ensembl_ftp = ftplib.FTP("ftp.ensembl.org")
ensembl_ftp.login()

for build in vcf_builds:
    vcf_loc = 'pub/current_variation/vcf/homo_sapiens/'
    if build == 'GRCh37':
        vcf_loc = 'pub/grch37/current/variation/vcf/homo_sapiens/'

    for chr in vcf_chrs:
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

ensembl_ftp.quit()
