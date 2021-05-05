import os

import click

from explore_pipolin.common import CONTEXT_SETTINGS
from explore_pipolin.extract_metadata_ena_xml import retrieve_ena_xml, extract_metadata
from explore_pipolin.massive_screening_ena_xml import retrieve_fasta_file, unzip_fasta_file


def yield_acc(found_pipolins_file):
    with open(found_pipolins_file) as inf:
        for line in inf:
            yield line.strip()


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('found-pipolins-file', type=click.Path(exists=True))
@click.argument('out-dir', type=click.Path())
def download_pipolin_genomes(found_pipolins_file, out_dir):
    os.makedirs(out_dir, exist_ok=True)

    for i, acc in enumerate(yield_acc(found_pipolins_file)):
        print(i + 1, acc, sep='\t')

        ena_xml = retrieve_ena_xml(acc)
        _, _, _, _, _, url = extract_metadata(ena_xml)

        file_path = os.path.join(out_dir, acc + '.fasta.gz')
        retrieve_fasta_file(url, file_path)
        unzip_fasta_file(file_path)


if __name__ == '__main__':
    download_pipolin_genomes()
