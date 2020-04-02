import requests
from requests.adapters import HTTPAdapter
from requests.exceptions import ConnectionError


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

    ensembl_adapter = HTTPAdapter(max_retries=3)
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
            results.update(r.json())
        except ConnectionError as ce:
            print(ce)
    return(results)