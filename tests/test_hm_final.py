import os.path, re
import json
import pandas as pd
import unittest
import gzip

class TestHmFinal(unittest.TestCase):

    hm_cols = {
        'chr_name': 'chr_name',
        'chr_pos': 'chr_position',
        'var_id': 'variant_id',
        'effect_allele': 'effect_allele',
        'other_allele': 'other_allele',
        'effect_weight': 'effect_weight',
        'hm_code': 'hm_code',
        'hm_info': 'hm_info'
    }
    cols = hm_cols.values()

    valid_data = {
        'allele': ['A','T','G','C','N','-'],
        'build': ['GRCh37','GRCh38'],
        'chr': ['1', '2', '3', '4', '5', '6', '7', '8',
                 '9', '10', '11', '12', '13', '14', '15', '16',
                 '17', '18', '19', '20', '21', '22',
                 'X', 'x', 'Y', 'y', 'XY', 'xy', 'MT', 'Mt', 'mt'],
        'code': ['5','4','3','1','0','-1','-4','-5'],
        'match': ['True', 'False'],
        'source': ['ENSEMBL','Author-reported']
    }

    hm_matched_header = '#Hm_variants_number_matched'
    hm_unmapped_header = '#Hm_variants_number_unmapped'
    hm_build_header = '#Hm_genome_build'
    hm_date = '#Hm_creation_date'

    header = [
        '##HARMONIZATION DETAILS',
        hm_build_header,
        hm_date,
        '#Hm_reference_source',
        hm_matched_header,
        hm_unmapped_header
    ]


    def __init__(self, file):
        self.hm_file = file
        self.hm_filename = os.path.basename(file)
        self.count_hm_matched = 0
        self.count_hm_unmapped = 0
        self.df_hm = pd.read_table(file, comment='#', engine = 'python')
        super().__init__()


    def test_header(self):
        count_header=0
        with gzip.open(self.hm_file , 'rb') as f:
            for f_line in f: # read file line by line
                line = f_line.decode()
                if line.startswith('#'):
                    for header_key in self.header:
                        if line.startswith(header_key):
                           count_header += 1
                        if line.startswith(self.hm_build_header):
                            valid_builds =  self.valid_data['build']
                            build = line.strip().split('=')[1]
                            self.assertTrue(build in valid_builds, msg=f"The HmPOS_build ({build}) is not in the list of valid builds ({valid_builds})")
                        if line.startswith(self.hm_date):
                            hm_vcf_date = line.strip().split('=')[1]
                            val = re.match('^\d{4}\-\d{2}\-\d{2}',hm_vcf_date)
                            self.assertTrue(bool(val))
                        if line.startswith(self.hm_matched_header):
                            self.count_hm_matched = int(line.strip().split('=')[1])
                        if line.startswith(self.hm_unmapped_header):
                            self.count_hm_unmapped = int(line.strip().split('=')[1])
                        
                else:
                    break
        self.assertEqual(count_header,len(self.header), msg=f'{self.hm_filename} - The number of Harmonization header information is different from the number expected!')


    def test_columns(self):
        for col in self.cols:
            self.assertTrue(col in self.df_hm.columns, msg=f'{self.hm_filename} - The column {col} is missing!')


    def test_data(self):
        chr_name_col = self.hm_cols['chr_name']
        df_chr_name = self.df_hm[chr_name_col]
        
        valid_chr = self.valid_data['chr']
        valid_allele = self.valid_data['allele']

        # Compare number of matched variants with 'HmVCF_n_matched'
        i_chr_notnull = (df_chr_name.isnull() == False)
        chr_notnull_entries = self.df_hm.loc[i_chr_notnull, chr_name_col]
        self.assertEquals(len(chr_notnull_entries), self.count_hm_matched)

        # Compare number of unmapped variants 'HmVCF_n_unmapped'
        i_chr_null = (df_chr_name.isnull() == True)
        chr_null_entries = self.df_hm.loc[i_chr_null, chr_name_col]
        self.assertEquals(len(chr_null_entries), self.count_hm_unmapped)

        # Test effect_allele
        hm_effect_alleles = set()
        df_hm_effect_allele = self.df_hm[self.hm_cols['effect_allele']]
        for effect_allele in dict(df_hm_effect_allele.value_counts()).keys():
            if not pd.isna(effect_allele):
                hm_effect_alleles.add(effect_allele)
        if hm_effect_alleles:
            for effect_allele in hm_effect_alleles:
                self.assertTrue(str(effect_allele) in valid_allele, msg=f"effect_allele '{effect_allele}' is not in the list of valid alleles: {valid_allele}!")
        
        # Test other_allele
        hm_other_alleles = set()
        df_hm_other_allele = self.df_hm[self.hm_cols['other_allele']]
        for other_allele in dict(df_hm_other_allele.value_counts()).keys():
            if not pd.isna(other_allele):
                hm_other_alleles.add(other_allele)
        if hm_other_alleles:
            for other_allele in hm_other_alleles:
                self.assertTrue(str(other_allele) in valid_allele, msg=f"other_allele '{other_allele}' is not in the list of valid alleles: {valid_allele}!")

        # Test effect_weight1
        df_effect_weight = self.df_hm[self.hm_cols['effect_weight']]
        self.assertEqual(df_effect_weight.dtypes, 'float64')

        # Test hm_code
        df_hm_code = self.df_hm[self.hm_cols['hm_code']]
        valid_codes = self.valid_data['code']
        for hm_code in dict(df_hm_code.value_counts()).keys():
            self.assertTrue(str(hm_code) in valid_codes, msg=f"hm_code '{hm_code}' is not in the list of valid codes: {valid_codes}!")

        # Test hm_info
        df_hm_info = self.df_hm[self.hm_cols['hm_info']]
        for hm_info in dict(df_hm_info.value_counts()).keys(): # Remove duplicates
            info = json.loads(hm_info)
            self.assertTrue('hm_source' in info.keys())
            # self.assertTrue('hm_match_chr' in info.keys())
            # self.assertTrue('hm_match_pos' in info.keys())


    def test_file(self):
        self.test_header()
        self.test_columns()
        self.test_data()