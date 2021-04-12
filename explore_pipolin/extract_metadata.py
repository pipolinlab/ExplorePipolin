from typing import Tuple
import xml.etree.ElementTree as ElementTree


import click


from explore_pipolin.common import CONTEXT_SETTINGS


def _extract_metadata(acc, ena_xml) -> Tuple[str, str, str, str, str]:
    asmbl_level, g_repr, tax_id, sci_name, asmbl_url = '-', '-', '-', '-', '-'
    for event, elem in ElementTree.iterparse(ena_xml, events=('start',)):
        if event == 'start' and elem.tag == 'ASSEMBLY' and elem.attrib['accession'] == acc:

            asmbl_level_node = elem.find('ASSEMBLY_LEVEL')
            if asmbl_level_node is not None:
                asmbl_level = asmbl_level_node.text

            g_repr_node = elem.find('GENOME_REPRESENTATION')
            if g_repr_node is not None:
                g_repr = g_repr_node.text

            tax_id_node = elem.find('TAXON/TAXON_ID')
            if tax_id_node is not None:
                tax_id = tax_id_node.text

            sci_name_node = elem.find('TAXON/SCIENTIFIC_NAME')
            if sci_name_node is not None:
                sci_name = sci_name_node.text

            for link in elem.findall('ASSEMBLY_LINKS/ASSEMBLY_LINK/URL_LINK[LABEL="WGS_SET_FASTA"]/URL'):
                asmbl_url = link.text

    return asmbl_level, g_repr, tax_id, sci_name, asmbl_url


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('accessions', type=click.Path(exists=True))
@click.argument('ena-xml', type=click.Path(exists=True))
@click.argument('output-file', type=click.Path())
def extract_metadata(accessions, ena_xml, output_file):
    with open(accessions) as inf, open(output_file, 'w') as ouf:
        print('assembly_accession', 'assembly_level', 'genome_representation', 'taxon_id',
              'scientific_name', 'assembly_url', sep='\t', file=ouf)
        for i, line in enumerate(inf):
            acc = line.strip()
            print(i + 1, acc, sep='\t')
            asmbl_level, g_repr, tax_id, sci_name, asmbl_url = _extract_metadata(acc, ena_xml)
            print(acc, asmbl_level, g_repr, tax_id, sci_name, asmbl_url, sep='\t', file=ouf)


if __name__ == '__main__':
    extract_metadata()
