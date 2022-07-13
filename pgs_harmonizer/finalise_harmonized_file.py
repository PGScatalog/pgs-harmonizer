import os
import sqlite3
import argparse
import gzip
from datetime import date


class FinaliseHarmonizedScoringFiles:

    headers_to_replace = [
        ('#HmPOS_build','#Hm_genome_build'),
        ('#HmVCF_ref','#Hm_reference_source'),
        ('#HmVCF_date', '#Hm_creation_date'),
        ('#HmVCF_n_matched','#Hm_variants_number_matched'),
        ('#HmVCF_n_unmapped','#Hm_variants_number_unmapped')
    ]

    def __init__(self, score_id, hmvcf_dir, staged_dir, genebuild, sqlite_file):

        if not os.path.exists(hmvcf_dir):
            print(f'Error: The path to the HmVCF directory can\'t be found ({hmvcf_dir}).')
            exit(1)
        if not os.path.exists(staged_dir):
            print(f'Error: The path to the staged directory ({staged_dir}).')
            exit(1)
        if not os.path.isfile(sqlite_file):
            print(f'Error: The path to the SQLite file can\'t be found ({sqlite_file}).')
            exit(1)

        self.score_id = score_id
        self.hmvcf_dir = hmvcf_dir
        self.staged_dir = staged_dir
        self.genebuild = genebuild
        self.sqlite_file = sqlite_file


    def get_header_information(self):
        with gzip.open(self.harmonized_filename , 'rb') as f:
            for f_line in f: # read file line by line
                line = f_line.decode()
                if line.startswith('#'):
                    if line.startswith('#variants_number'):
                        self.variant_number = line.strip().split('=')[1]
                    elif line.startswith('#HmVCF_n_matched'):
                        self.matched_number = line.strip().split('=')[1]
                    elif line.startswith('#HmVCF_n_unmapped'):
                        self.unmapped_number = line.strip().split('=')[1]
                    elif line.startswith('#HmVCF_ref'):
                        self.source = line.strip().split('=')[1]
                else:
                    break


    def create_entry_in_knowledge_base(self):
        self.version = 1
        comment = 'First version of the harmonized file'
        creation_date = date.today()

        sqlite_connection = sqlite3.connect(self.sqlite_file)
        sqlite_cursor = sqlite_connection.cursor()
        sql_query = f'SELECT max(version) from harmonization_version_control WHERE pgs_id=? and genebuild=?;'
        sqlite_cursor.execute(sql_query,[self.score_id,self.genebuild])
        sqlite_result = sqlite_cursor.fetchone()
        
        if len(sqlite_result):
            if (sqlite_result[0]):
                self.version = int(sqlite_result[0]) + 1
                comment = 'Updated version of the harmonized file'

        # Build file names
        self.harmonized_filename = f'{self.hmvcf_dir}/{self.score_id}_hmVCF_GRCh{self.genebuild}.txt.gz'
        self.staged_filename = f'{self.staged_dir}/{self.score_id}_h{self.genebuild}_v{self.version}.txt.gz'

        # Get header information
        self.get_header_information()

        try:
            sql_query2 = 'INSERT INTO harmonization_version_control (pgs_id,version,genebuild,date,source,variant_number,matched_number,unmapped_number,comment) VALUES (?,?,?,?,?,?,?,?,?);'
            sqlite_cursor.execute(sql_query2,[self.score_id,self.version,self.genebuild,creation_date,self.source,self.variant_number,self.matched_number,self.unmapped_number,comment])
            sqlite_connection.commit()
            sqlite_cursor.close()
            sqlite_connection.close()
        except sqlite3.Error as e:
            print(f"Failed to insert record of {self.score_id} [{self.genebuild}] into sqlite table: {e}")
            exit(1)   


    def update_header_comments(self):
        lines = ''
        count_lines = 0
        max_lines = 200
        with gzip.open(self.harmonized_filename, 'rb') as f_in:
            with gzip.open(self.staged_filename, 'wb') as f_out:
                for f_line in f_in:
                    line = f_line.decode()
                    # Update header
                    if line.startswith('#'):
                        if line.startswith('##HARMONIZATION DETAILS') and self.version:
                            line += f'#Hm_file_version={self.version}\n'
                        elif line.startswith('#HmPOS_date') or line.startswith('#HmPOS_match_'):
                            line = ''
                        else:
                            for vcf_header,new_header in self.headers_to_replace:
                                if line.startswith(vcf_header):
                                    line = f'{line.replace(vcf_header,new_header)}'
                                    break
                    # Any lines
                    if '\r' in line:
                        line = line.replace('\r','\n')
                    elif not '\n' in line and line != '':
                        line += '\n'

                    # Add lines
                    lines += line
                    count_lines += 1
                    if count_lines == max_lines:
                        f_out.write(lines.encode())
                        lines = ''
                        count_lines = 0
                if count_lines:
                    f_out.write(lines.encode())    


######################################################################

def main():
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--score_id", help='List of scores IDs', required=True, metavar='PGS_IDS')
    argparser.add_argument("--input_dir", help='Directory hosting the harmonized scoring files', required=True, metavar='VCF_DIR')
    argparser.add_argument("--staged_dir", help='Path to the variants output file', required=True, metavar='STG_DIR')
    argparser.add_argument("--genebuild", help='Directory hosting the scoring files', required=True, metavar='GENEBUILD')
    argparser.add_argument("--sqlite_file", help='Path to the SQLlite file containing the variants with coordinates already assigned', required=True, metavar='SQLITE_FILE')

    args = argparser.parse_args()

    finalise_harmonized_file = FinaliseHarmonizedScoringFiles(args.score_id,args.input_dir,args.staged_dir,args.genebuild,args.sqlite_file)
    finalise_harmonized_file.create_entry_in_knowledge_base()
    finalise_harmonized_file.update_header_comments()


if __name__== "__main__":
    main()