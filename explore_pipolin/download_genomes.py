import os

import click

from explore_pipolin.common import CONTEXT_SETTINGS
from explore_pipolin.collect_metadata import retrieve_ena_xml, extract_metadata
from explore_pipolin.massive_screening import download_genome


def yield_acc(found_pipolins_file):
    with open(found_pipolins_file) as inf:
        for line in inf:
            yield line.strip()


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('accessions', type=click.Path(exists=True))
@click.argument('out-dir', type=click.Path())
def main(accessions, out_dir):
    """
    ACCESSIONS is a file with accession ids (e.g., found_pipolins.txt) for which
    the genomes will be downloaded from the ENA database.
    """
    os.makedirs(out_dir, exist_ok=True)

    for i, acc in enumerate(yield_acc(accessions)):
        print(i + 1, acc, sep='\t')

        ena_xml = retrieve_ena_xml(acc)
        _, _, _, _, _, url = extract_metadata(ena_xml)

        genome = os.path.join(out_dir, acc + '.fasta')
        download_genome(url, genome)


if __name__ == '__main__':
    main()
