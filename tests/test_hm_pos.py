import os.path, re
import pandas as pd
import unittest
import gzip

class TestHmPOS(unittest.TestCase):


    hm_cols = {
        'chr_name': 'chr_name',
        'chr_pos': 'chr_position',
        'effect_allele': 'effect_allele',
        'other_allele': 'other_allele',
        'effect_weight': 'effect_weight',
        'source': 'hm_source',
        'hm_rsID': 'hm_rsID',
        'hm_chr': 'hm_chr',
        'hm_pos': 'hm_pos',
        'hm_inferOtherAllele': 'hm_inferOtherAllele',
        'match_chr': 'hm_match_chr',
        'match_pos': 'hm_match_pos'
    }
    cols = hm_cols.values()

    valid_data = {
        'allele': ['A','T','G','C','N','-'],
        'chr': ['1', '2', '3', '4', '5', '6', '7', '8',
                 '9', '10', '11', '12', '13', '14', '15', '16',
                 '17', '18', '19', '20', '21', '22',
                 'X', 'x', 'Y', 'y', 'XY', 'xy', 'MT', 'Mt', 'mt'],
        'match': ['True', 'False'],
        'source': ['ENSEMBL','Author-reported']
    }

    hm_date = '#HmPOS_date'

    header = [
        '##HARMONIZATION DETAILS',
        '#HmPOS_build',
        hm_date,
        '#HmPOS_match_chr',
        '#HmPOS_match_pos'
    ]

   
    def __init__(self, file):
        self.hm_file = file
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
                        if line.startswith(self.hm_date):
                            hm_pos_date = line.strip().split('=')[1]
                            val = re.match('^\d{4}\-\d{2}\-\d{2}',hm_pos_date)
                            self.assertTrue(bool(val))
                else:
                    break
        self.assertEqual(count_header,len(self.header), msg='The number of Harmonization header information is different from the number expected!')


    def test_columns(self):
        for col in self.cols:
            self.assertTrue(col in self.df_hm.columns, msg=f'The column {col} is missing!')


    def test_data(self):
        
        valid_chr = self.valid_data['chr']
        valid_allele = self.valid_data['allele']

        # Test chr_name
        df_hm_chr_name = self.df_hm[self.hm_cols['chr_name']]
        for chr_name in dict(df_hm_chr_name.value_counts()).keys():
            self.assertTrue(str(chr_name) in valid_chr, msg=f"chr_name '{chr_name}' is not in the list of valid chromosomes: {valid_chr}!")

        # Test chr_position
        df_hm_chr_pos = self.df_hm[self.hm_cols['chr_pos']]
        self.assertEqual(df_hm_chr_pos.dtypes, 'int64')

        # Test effect_allele
        df_effect_allele = self.df_hm[self.hm_cols['effect_allele']]
        for effect_allele in dict(df_effect_allele.value_counts()).keys():
            self.assertTrue(str(effect_allele) in valid_allele, msg=f"effect_allele '{effect_allele}' is not in the list of valid alleles: {valid_allele}!")
        
        # Test other_allele
        df_other_allele = self.df_hm[self.hm_cols['other_allele']]
        for other_allele in dict(df_other_allele.value_counts()).keys():
            self.assertTrue(str(other_allele) in valid_allele, msg=f"other_allele '{other_allele}' is not in the list of valid alleles: {valid_allele}!")

        # Test effect_weight
        df_effect_weight = self.df_hm[self.hm_cols['effect_weight']]
        self.assertEqual(df_effect_weight.dtypes, 'float64')
        
        # Test hm_source
        df_hm_source = self.df_hm[self.hm_cols['source']]
        valid_sources = self.valid_data['source']
        for hm_source in dict(df_hm_source.value_counts()).keys():
            self.assertTrue(str(hm_source) in valid_sources, msg=f"hm_source '{hm_source}' is not in the list of valid sources: {valid_sources}!")

        # Test hm_rsID
        rsIDs = set()
        df_hm_rsid = self.df_hm[self.hm_cols['hm_rsID']]
        for hm_rsid in dict(df_hm_rsid.value_counts()).keys():
            if not pd.isna(hm_rsid):
                rsIDs.add(hm_rsid)
        if rsIDs:
            for rsID in rsIDs:
                self.assertTrue(rsID.startswith('rs'))

        # Test hm_chr
        df_hm_chr = self.df_hm[self.hm_cols['hm_chr']]
        for hm_chr in dict(df_hm_chr.value_counts()).keys():
            self.assertTrue(str(hm_chr) in valid_chr, msg=f"hm_chr '{hm_chr}' is not in the list of valid chromosomes: {valid_chr}!")

        # Test hm_pos
        df_hm_pos = self.df_hm[self.hm_cols['hm_pos']]
        self.assertEqual(df_hm_pos.dtypes, 'int64')

        # hm_inferOtherAllele
        df_hm_infer = self.df_hm[self.hm_cols['hm_inferOtherAllele']]
        self.assertEqual(df_hm_infer.dtypes, 'float64')

        # hm_match_chr
        df_hm_match_chr = self.df_hm[self.hm_cols['match_chr']]
        valid_match = self.valid_data['match']
        for hm_match_chr in dict(df_hm_match_chr.value_counts()).keys():
            self.assertTrue(str(hm_match_chr) in valid_match, msg=f"hm_match_chr '{hm_match_chr}' is not in the list of valid values: {valid_match}!")

        # hm_match_pos
        df_hm_match_pos = self.df_hm[self.hm_cols['match_pos']]
        valid_match = self.valid_data['match']
        for hm_match_pos in dict(df_hm_match_pos.value_counts()).keys():
            self.assertTrue(str(hm_match_pos) in valid_match, msg=f"hm_match_pos '{hm_match_pos}' is not in the list of valid values: {valid_match}!")


    def test_file(self):
        self.test_header()
        self.test_columns()
        self.test_data()