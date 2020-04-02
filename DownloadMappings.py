import ftplib

download_chains = False

download_vcf = True
vcf_builds = ['GRCh37', 'GRCh38']
vcf_chrs = ['1', '2', '6']

ensembl_ftp = ftplib.FTP("ftp.ensembl.org")
ensembl_ftp.login()

if download_chains:
    list_chain = ensembl_ftp.nlst('pub/current_assembly_chain/homo_sapiens/*chain.gz')
    for fn_chain in list_chain:
        with open("./map/liftover_chain/{}".format(fn_chain.split('/')[-1]), 'wb') as lf:
            ensembl_ftp.retrbinary('RETR {}'.format(fn_chain), lf.write)

if download_vcf:
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
                if build == 'GRCh37':
                    fn_vcf = 'homo_sapiens-chr{}.vcf.gz.tbi'.format(chr)
                vcf_loc_chr = vcf_loc + fn_vcf
                with open("./map/vcf_ref/{}/{}".format(build, fn_vcf), 'wb') as lf:
                    ensembl_ftp.retrbinary('RETR {}'.format(vcf_loc_chr), lf.write)
            except:
                print('! ERROR downloading:', build, 'chr{}'.format(chr))

ensembl_ftp.quit()
