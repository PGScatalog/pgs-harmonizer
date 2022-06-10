import sys, os.path
from argparse import Namespace
import unittest
from test_hm_pos import TestHmPOS
from test_hm_vcf import TestHmVCF
from test_hm_final import TestHmFinal
sys_paths = ['../','./']
sys.path.extend(sys_paths)
from Harmonize import run_HmPOS,run_HmVCF
from pgs_harmonizer.finalise_harmonized_file import FinaliseHarmonizedScoringFiles

class TestArgs:
    
    def __init__(self, data):
        for name, value in data.items():
            setattr(self, name, self._wrap(value))

    def _wrap(self, value):
        return TestArgs(value) if isinstance(value, dict) else value


class TestHm(unittest.TestCase):

    current_dir = os.path.dirname(os.path.abspath(__file__))
    test_dir  = current_dir+'/data/'
    hmPOS_dir = test_dir+'hmPOS'
    hmVCF_dir = test_dir+'hmVCF'
    hmFinal_dir = test_dir+'hmFinal'
    hm_knowledge_base = test_dir+'kb/harmonization_kb.db'

    genebuilds = [37] #[37,38]
    pgs_ids = ['PGS000015','PGS000065']

    def generate_HmPOS_files(self):
        self.create_hm_directory(self.hmPOS_dir)

        for gb in self.genebuilds:
            hmPOS_gb_dir = f'{self.hmPOS_dir}/{gb}'
            target_assembly = f'GRCh{gb}'
            self.create_hm_directory(hmPOS_gb_dir)

            for pgs_id in self.pgs_ids:
                hm_args = Namespace(
                    pgs_id=pgs_id,
                    loc_scorefiles=self.test_dir,
                    loc_outputs=hmPOS_gb_dir,
                    var2location=self.test_dir+'kb',
                    target_build=target_assembly,
                    source_build=None,
                    useAPI=False,
                    catchmissingAPI=False,
                    silent_tqdm=False,
                    ignore_rsid=False,
                    gzip=True
                )
                run_HmPOS(hm_args)


    def generate_HmVCF_files(self):
        self.create_hm_directory(self.hmVCF_dir)

        for gb in self.genebuilds:
            hmPOS_gb_dir = f'{self.hmPOS_dir}/{gb}'
            hmVCF_gb_dir = f'{self.hmVCF_dir}/{gb}'

            target_assembly = f'GRCh{gb}'
            self.create_hm_directory(hmVCF_gb_dir)

            for pgs_id in self.pgs_ids:
                    hm_args = Namespace(
                        pgs_id=pgs_id,
                        loc_scorefiles=hmPOS_gb_dir,
                        loc_outputs=hmVCF_gb_dir,
                        loc_vcfref=f'{self.test_dir}/vcfs/{target_assembly}',
                        target_build=target_assembly,
                        cohort_name=None,
                        addOtherAllele=False,
                        addVariantID=False,
                        author_reported=False,
                        skip_strandflips=False,
                        split_unmappable=False,
                        keep_duplicates=False,
                        gzip=True,
                        silent_tqdm=False
                    )
                    run_HmVCF(hm_args)


    def generate_HmFinal_files(self):
        self.create_hm_directory(self.hmFinal_dir)

        for gb in self.genebuilds:
            hmVCF_gb_dir = f'{self.hmVCF_dir}/{gb}'
            hmFinal_gb_dir = f'{self.hmFinal_dir}/{gb}'

            self.create_hm_directory(hmFinal_gb_dir)

            for pgs_id in self.pgs_ids:
                finalise_harmonized_file = FinaliseHarmonizedScoringFiles(pgs_id,hmVCF_gb_dir,hmFinal_gb_dir,gb,self.hm_knowledge_base)
                finalise_harmonized_file.create_entry_in_knowledge_base()
                finalise_harmonized_file.update_header_comments()


    def create_hm_directory(self,path):
        """
        Creates Hm directory for a given PGS
        > Parameters:
            - path: path of the directory
        """
        # Create directory if it doesn't exist
        if not os.path.isdir(path): 
            try:
                os.mkdir(path, 0o755)
            except OSError:
                print (f'Creation of the directory {path} failed')
                exit()

    

if __name__ == "__main__":
    test_hm = TestHm()
    print("- Test files creation")
    test_hm.generate_HmPOS_files()
    test_hm.generate_HmVCF_files()
    test_hm.generate_HmFinal_files()

    print("- Test files formats and contents")
    for gb in test_hm.genebuilds:
        # Test HmPOS
        hmPOS_gb_dir = f'{test_hm.hmPOS_dir}/{gb}'
        for pgs_id in test_hm.pgs_ids:
            # print(f"# {pgs_id} - {gb}")
            hmPOS_file = f'{hmPOS_gb_dir}/{pgs_id}_hmPOS_GRCh{gb}.txt.gz'
            test_hmPOS = TestHmPOS(hmPOS_file)
            test_hmPOS.test_file()
        
        # Test HmVCF
        hmVCF_gb_dir = f'{test_hm.hmVCF_dir}/{gb}'
        for pgs_id in test_hm.pgs_ids:
            # print(f"# {pgs_id} - {gb}")
            hmVCF_file = f'{hmVCF_gb_dir}/{pgs_id}_hmVCF_GRCh{gb}.txt.gz'
            test_hmVCF = TestHmVCF(hmVCF_file)
            test_hmVCF.test_file()

        # Test HmFinal
        hmFinal_gb_dir = f'{test_hm.hmFinal_dir}/{gb}'
        for pgs_id in test_hm.pgs_ids:
            
            hmFinal_file = f'{hmFinal_gb_dir}/{pgs_id}_h{gb}_v1.txt.gz'
            test_hmFinal = TestHmFinal(hmFinal_file)
            test_hmFinal.test_file()