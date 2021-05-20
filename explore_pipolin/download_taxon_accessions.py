import os
import sys
from typing import Mapping

import click

import requests

from explore_pipolin.common import CONTEXT_SETTINGS


def request_taxon_data(taxon_id: int):
    url = "https://www.ebi.ac.uk/ena/portal/api/search"
    params = {'query': f'tax_tree({taxon_id})',
              'result': 'assembly',
              'format': 'json',
              'limit': '0',
              'offset': '0'}
    print(f'Requesting data for the taxon {taxon_id}...')
    req = requests.get(url=url, params=params)

    if not req.ok:
        req.raise_for_status()
        sys.exit()
    print('Finished requesting data!')
    return req.json()


def create_acc_version_dict(json_data) -> Mapping[str, str]:
    acc_version_dict = {}
    for assembly in json_data:
        acc = assembly['accession']
        version = assembly['version']
        if acc in acc_version_dict:
            if acc_version_dict[acc] < version:
                acc_version_dict[acc] = version
        else:
            acc_version_dict[acc] = version
    return acc_version_dict


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('taxon-id', type=int)
@click.option('--out-file', type=click.Path())
def main(taxon_id, out_file):
    """
    Given NCBI taxon ID, a list of assembly accessions, registered under the given taxon tree,
    will be downloaded from the ENA database.
    """
    data = request_taxon_data(taxon_id)
    print(f'Number of downloaded assembly accessions: {len(data)}')

    print('NOTE: Only the last version of each assembly accession will be saved!')
    acc_version_dict = create_acc_version_dict(data)
    print(f'Number of assembly accessions after filtering: {len(acc_version_dict)}')

    out_file = out_file if out_file else os.path.join(os.getcwd(), f'{taxon_id}.accessions.txt')

    with open(out_file, 'w') as ouf:
        for acc, version in acc_version_dict.items():
            print('.'.join([acc, version]), file=ouf)
    print(f'Finished! The results can be found in {out_file}')


if __name__ == '__main__':
    main()
