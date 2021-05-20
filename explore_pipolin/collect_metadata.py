import os
import sys
from typing import Tuple
import xml.etree.ElementTree as ElementTree


import click
import requests

from explore_pipolin.common import CONTEXT_SETTINGS


def retrieve_ena_xml(acc) -> str:
    server = "https://www.ebi.ac.uk/ena/browser/api/xml/" + acc
    r = requests.get(server)

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    return r.text


def extract_metadata(ena_xml: str) -> Tuple[str, str, str, str, str, str]:
    asmbl_level, g_repr, tax_id, sci_name, strain, asmbl_url = '-', '-', '-', '-', '-', '-'
    root = ElementTree.fromstring(ena_xml)

    asmbl_level_node = root.find('ASSEMBLY/ASSEMBLY_LEVEL')
    if asmbl_level_node is not None:
        asmbl_level = asmbl_level_node.text

    g_repr_node = root.find('ASSEMBLY/GENOME_REPRESENTATION')
    if g_repr_node is not None:
        g_repr = g_repr_node.text

    tax_id_node = root.find('ASSEMBLY/TAXON/TAXON_ID')
    if tax_id_node is not None:
        tax_id = tax_id_node.text

    sci_name_node = root.find('ASSEMBLY/TAXON/SCIENTIFIC_NAME')
    if sci_name_node is not None:
        sci_name = sci_name_node.text

    strain_node = root.find('ASSEMBLY/TAXON/STRAIN')
    if strain_node is not None:
        strain = strain_node.text

    for link in root.findall('ASSEMBLY/ASSEMBLY_LINKS/ASSEMBLY_LINK/URL_LINK[LABEL="WGS_SET_FASTA"]/URL'):
        asmbl_url = link.text

    return asmbl_level, g_repr, tax_id, sci_name, strain, asmbl_url


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('accessions', type=click.Path(exists=True))
@click.option('--out-file', type=click.Path())
def main(accessions, out_file):
    """
    ACCESSIONS is a file with accession ids (e.g., found_pipolins.txt) for which
    the metadata will be downloaded and extracted from the ENA database.
    """
    default_out_file = os.path.join(
        os.getcwd(), os.path.splitext(os.path.basename(accessions))[0] + '.metadata.txt'
    )
    out_file = out_file if out_file else default_out_file

    with open(accessions) as inf, open(out_file, 'w') as ouf:
        print('assembly_accession', 'assembly_level', 'genome_representation', 'taxon_id',
              'scientific_name', 'strain', 'assembly_url', sep='\t', file=ouf)
        for i, line in enumerate(inf):
            acc = line.strip()
            print(i + 1, acc, sep='\t')

            ena_xml = retrieve_ena_xml(acc)
            asmbl_level, g_repr, tax_id, sci_name, strain, asmbl_url = extract_metadata(ena_xml)
            print(acc, asmbl_level, g_repr, tax_id, sci_name, strain, asmbl_url,
                  sep='\t', file=ouf, flush=True)

    print(f'Finished! The results can be found in {out_file}')


if __name__ == '__main__':
    main()
