import requests
import time
from requests.adapters import HTTPAdapter
from requests.exceptions import ConnectionError
import pandas as pd
from pgs_harmonizer.harmonize import reversecomplement


class VariationResult:
    """Class to parse the 'mapping 'information from ENSEMBL Variation"""
    def __init__(self, id, json_result):
        self.id = id
        self.json_result = json_result

        self.chrom = None
        self.bp = None
        self.hm_code = None
        self.alleles = None

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
            self.hm_code = 1
            self.alleles = alleles[0]

        return self.chrom, self.bp, self.hm_code, self.alleles

    def check_alleles(self, ref = None, eff = None):
        """Check if the original scoring file's alleles match the ENSEMBL rsID mapping allele string
        (NB: The allele string encoding of INDELs is different than the VCF allele notation)"""
        hm_consistent = []
        hm_revcomp = []

        # Check Effect Allele
        if eff in self.alleles:
            hm_consistent.append('eff')
        elif eff in map(reversecomplement, self.alleles):
            hm_revcomp.append('eff')

        # Check reference allele
        if ref is not None:
            if ref in self.alleles:
                hm_consistent.append('ref')
            elif ref in map(reversecomplement, self.alleles):
                hm_revcomp.append('ref')

        # Check the alleles
        if ref is None: # Just check effect allele
            if 'eff' in hm_revcomp:
                self.hm_code = -2
            elif 'eff' not in hm_consistent:
                self.hm_code = -3
        else: # Check both alleles
            if 'eff' and 'ref' in hm_consistent:
                self.hm_code = 1
            elif 'eff' and 'ref' in hm_revcomp:
                self.hm_code = -2
            else:
                self.hm_code = -3

        return self.chrom, self.bp, self.hm_code

    def synonyms(self):
        return self.json_result['synonyms']


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

def parse_var2location(loc_var2location_results):
    """Reads results of var2location.pl mapping into the same class as the ENSEMBL API results"""
    mappings = pd.read_csv(loc_var2location_results, sep = '\t',
                           names=['query_rsid', 'mapped_rsid', 'allele_string', 'seq_region_name', 'start', 'end'])
    mappings['seq_region_name'] = mappings['seq_region_name'].astype('str') # Asssert character in case there are only chr numbers
    results = {}

    # Loop through results and parse to Variant
    for query_rsid, maps in mappings.groupby('query_rsid'):
        q_json = {'name': maps.iloc[0,1], 'mappings': []}
        syn = list(set(maps['query_rsid']).union(maps['mapped_rsid']))
        for m in maps.iterrows():
            q_json['mappings'].append(dict(m[1]))
        v = VariationResult(maps.iloc[0,1], q_json)
        for s in syn:
            results[s] = v

    return results