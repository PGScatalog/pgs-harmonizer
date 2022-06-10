import requests
import argparse

def run_rest_api(num_from,num_to,server):
    """"
    Method using the PGS Catalog REST API to retrieve all the PGS IDs 
    within the range of numeric IDs provided
    > Parameters:
        - num_from: Numeric range start ID
        - num_to: Numeric range end ID
        - server: Server URL
    > Return type: array
    """
    if not server.endswith('/'):
        server += '/'
    rest_full_url = server+'rest/prod/range/'
    if num_from and num_to:
        rest_full_url += f'?from={num_from}&to={num_to}'
    try:
        response = requests.get(rest_full_url)
        results = response.json()
    except requests.exceptions.RequestException as e:  # This is the correct syntax
        raise SystemExit(e)
    return results



################################################################################

def main():
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--num_from", help='Numeric ID starting the range', required=True, metavar='FROM')
    argparser.add_argument("--num_to", help='Numeric ID ending the rangee', required=True, metavar='TO')
    argparser.add_argument("--output", help='Output file', required=True, metavar='OUT')
    argparser.add_argument("--rest_server", help='URL to REST API server', required=True, metavar='REST')

    args = argparser.parse_args()

    pgs_ids_list = run_rest_api(args.num_from,args.num_to,args.rest_server)

    output_file = open(args.output,'w')
    for pgs_id in pgs_ids_list:
        output_file.write(f"{pgs_id}\n")
    output_file.close()

if __name__ == '__main__':
    main()