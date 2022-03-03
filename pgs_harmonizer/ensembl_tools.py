import requests
import time
from requests.adapters import HTTPAdapter
from requests.exceptions import ConnectionError
import pandas as pd
import sqlite3
from pgs_harmonizer.harmonize import reversecomplement


class VariationResult:
    """Class to parse the 'mapping 'information from ENSEMBL Variation"""
    def __init__(self, id, json_result):
        self.id = id
        self.json_result = json_result

        self.chrom = None
        self.bp = None
        self.alleles = None
        self.alleles_rc = None

    def select_canonical_data(self, chromosomes):
        """To identify the best mapping (adapted from GWAS Catalog)"""

        mapped_data = self.json_result['mappings']

        chrom = []
        bp = []
        alleles = []
        for mapping in mapped_data:
            if mapping['seq_region_name'] in chromosomes:
                # print(mapping['seq_region_name'], mapping['start'])
                bp.append(mapping['start'])
                chrom.append(mapping['seq_region_name'])
                alleles.append(mapping['allele_string'].split('/'))
        if (len(bp) == 1) or (len(bp) > 1 and all_same(bp)):
            self.chrom = chrom[0]
            self.bp = bp[0]
            self.alleles = alleles[0]
            self.alleles_rc = [reversecomplement(x) for x in self.alleles]

        return self.chrom, self.bp, self.alleles

    def check_alleles(self, ref = None, eff = None):
        """Check if the original scoring file's alleles match the ENSEMBL rsID mapping allele string
        (NB: The allele string encoding of INDELs is different than the VCF allele notation)"""
        hm_consistent = []
        hm_revcomp = []

        # Check Effect Allele
        if eff in self.alleles:
            hm_consistent.append('eff')
        elif eff in self.alleles_rc:
            hm_revcomp.append('eff')

        # Check reference allele
        if ref is not None:
            if ref in self.alleles:
                hm_consistent.append('ref')
            elif ref in self.alleles_rc:
                hm_revcomp.append('ref')

        # isPalindromic
        isPalindromic = False
        for allele in self.alleles:
            if allele in self.alleles_rc:
                isPalindromic = True

        # Check the alleles
        # return hm_matchesVCF, hm_isPalindromic, hm_isFlipped
        if ref is None:
            # Just check effect allele
            if 'eff' in hm_revcomp:
                return True, isPalindromic, True
            elif 'eff' not in hm_consistent:
                return False, isPalindromic, False
        else:
            # Check both alleles
            if 'eff' and 'ref' in hm_consistent:
                return True, isPalindromic, False
            elif 'eff' and 'ref' in hm_revcomp:
                return True, isPalindromic, True
            else:
                return False, isPalindromic, False

    def infer_OtherAllele(self, eff):
        """Try to infer the reference_allele. Report all possible reference alleles '/'-delimited"""
        try:
            oa = None
            if eff in self.alleles:
                REF = self.alleles[0]
                ALT = self.alleles[1:]

                if eff in ALT:
                    oa = REF
                elif len(ALT) == 1:
                    oa = ALT[0]
                else:
                    oa = '/'.join(ALT)
        except:
            oa = None
        return oa

    def synonyms(self):
        return self.json_result.get('synonyms')


def all_same(items):
    return all(x == items[0] for x in items)


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def ensembl_post(rsid_list, build = 'GRCh38'):
    """Retrieve rsID info from ENSEMBL Variation API"""
    # Assign API URL/settings
    valid_build = ['GRCh37', 'GRCh38']
    if build not in valid_build:
        raise ValueError("results: genome build must be one of {}".format(valid_build))

    headers = {"Content-Type": "application/json", "Accept": "application/json"}

    ensembl_adapter = HTTPAdapter(max_retries=3) # ToDo add handler to handle retry-after header
    session = requests.Session()

    url = "https://rest.ensembl.org"
    if build == 'GRCh37':
        url = "https://grch37.rest.ensembl.org"

    session.mount(url, ensembl_adapter)

    results = {}
    # Loop throught the rsID list and add the results to a dictionary
    for c_ids in chunks(rsid_list, 200):
        payload = {'ids': c_ids }
        try:
            r = session.post(url + '/variation/homo_sapiens', headers=headers, json=payload)
            if 'Retry-After' in r.headers:
                retry = r.headers['Retry-After']
                time.sleep(float(retry))  # pragma: no cover
                r = session.post(url + '/variation/homo_sapiens', headers=headers, json=payload)
            else:
                for i,j in r.json().items():
                    v = VariationResult(i, j) #Class object
                    results[i] = v
                    for syn in v.synonyms():
                        results[syn] = v
        except ConnectionError as ce:
            print(ce)
    return results

def clean_rsIDs(raw_rslist):
    """Takes a list of values, removes anything that doesn't look like an rsID and splits any variants that
    are haplotypes, combinations, or interactions"""
    cln_rslist = set()
    for x in raw_rslist:
        if type(x) is str and x.startswith('rs'):
            if '_x_' in x:
                x = [y.strip() for y in x.split('_x_')]
            elif ';' in x:
                x = [y.strip() for y in x.split(';')]
            elif ',' in x:
                x = [y.strip() for y in x.split(',')]
            else:
                cln_rslist.add(x)

            if type(x) == list:
                for i in x:
                    if i.startswith('rs'):
                        cln_rslist.add(i)
    return(list(cln_rslist))

def parse_var2location(loc_var2location_db, rsIDs = None, catchAPI=True):
    """Reads results of var2location DB into the same class as the ENSEMBL API results"""
    results = {}

    sqlite_connection = sqlite3.connect(loc_var2location_db)

    # Too Slow
    # if type(rsIDs) == list:
    #     print('reading filtered')
    #     var2location = pd.read_sql_query("SELECT * FROM variant_coords WHERE varname IN {}".format(tuple(rsIDs)),
    #                                      sqlite_connection)
    # else:
    #     print('reading all')
    #     var2location = pd.read_sql_query("SELECT * FROM variant_coords", sqlite_connection)
    #
    # var2location.columns = ['query_rsid', 'current_rsid', 'seq_region_name', 'start', 'end', 'allele_string']
    #
    # print('Looping')
    # for query_rsid, maps in var2location.groupby('query_rsid'):
    #     mapped_ids = list(set(maps['current_rsid']))
    #     if len(mapped_ids) == 1:
    #         results[query_rsid] = VariationResult(mapped_ids[0], {'name': mapped_ids[0],
    #                                                      'mappings': maps[['seq_region_name', 'start', 'end', 'allele_string']].to_dict('records')})
    #     else:
    #         print('WARNING')

    sqlite_cursor = sqlite_connection.cursor()
    if type(rsIDs) == list:
        qcursor = sqlite_cursor.execute("SELECT vlist.varname, vlist.current_varname, vc.chr, vc.start, vc.end, vc.alleles  FROM variant vlist LEFT JOIN variant_coords vc ON vlist.varname = vc.varname WHERE vlist.varname IN {}".format(tuple(rsIDs)))
    else:
        qcursor = sqlite_cursor.execute("SELECT vlist.varname, vlist.current_varname, vc.chr, vc.start, vc.end, vc.alleles  FROM variant vlist LEFT JOIN variant_coords vc ON vlist.varname = vc.varname")

    for line in qcursor:
        query_rsid, current_rsid, seq_region_name, start, end, allele_string = line
        if query_rsid in results:
            results[query_rsid].json_result['mappings'].append({'seq_region_name': seq_region_name,
                                                                               'start': start,
                                                                               'end': end,
                                                                               'allele_string': allele_string})
        else:
            results[query_rsid] = VariationResult(current_rsid, {'name': current_rsid,
                                                                 'mappings': [{'seq_region_name': seq_region_name,
                                                                               'start': start,
                                                                               'end': end,
                                                                               'allele_string': allele_string}]})

    # Cleanup w/ ENSEMBL API
    if catchAPI is True:
        id_set = set(rsIDs)
        missing = id_set.difference(results.keys())
        if len(missing) > 0:
            print('! Attempting to get {} missing mappings via ENSEMBL'.format(len(missing)))
            mapping_ensembl = ensembl_post(list(missing))
            if len(mapping_ensembl) > 0:
                results.update(mapping_ensembl)

    return results
