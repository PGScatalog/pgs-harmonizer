import requests
import time
from requests.adapters import HTTPAdapter
from requests.exceptions import ConnectionError


class VariantResult:
    def __init__(self, id, json_result):
        self.id = id
        self.json_result = json_result

    def select_canonical_data(self, chromosomes):
        '''To identify the best mapping (adapted from GWAS Catalog)'''

        mapped_data = self.json_result['mappings']

        chrom = []
        bp = []
        for mapping in mapped_data:
            if mapping['seq_region_name'] in chromosomes:
                # print(mapping['seq_region_name'], mapping['start'])
                bp.append(mapping['start'])
                chrom.append(mapping['seq_region_name'])
        if len(bp) == 1:
            return chrom[0], bp[0], 'Mapped (rsID)'
        elif len(bp) > 1 and all_same(bp):
            return chrom[0], bp[0], 'Mapped (rsID)'  # "AMBIGUOUS"
        else:
            return None, None, None  # to catch those where they only map to a patch

    def synonyms(self):
        return self.json_result['synonyms']


def all_same(items):
    return all(x == items[0] for x in items)

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def ensembl_post(rsid_list, build = 'GRCh38'):
    '''Retrieve rsID info from ENSEMBL Variation API'''

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
    else:
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
                    v = VariantResult(i, j) #Class object
                    results[i] = v
                    for syn in v.synonyms():
                        results[syn] = v
        except ConnectionError as ce:
            print(ce)
    return results

