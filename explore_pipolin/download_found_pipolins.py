import os

import click

from explore_pipolin.common import CONTEXT_SETTINGS
from explore_pipolin.extract_metadata_ena_xml import _retrieve_ena_xml, _extract_metadata
from explore_pipolin.massive_screening_ena_xml import retrieve_fasta_file, unzip_fasta_file


def yield_acc_for_pipolin(found_pipolins_file):
    with open(found_pipolins_file) as inf:
        for line in inf:
            yield line.strip()


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('found-pipolins', type=click.Path(exists=True))
@click.argument('out-dir', type=click.Path())
def download_found_pipolins(found_pipolins, out_dir):
    os.makedirs(out_dir, exist_ok=True)

    for i, acc in enumerate(yield_acc_for_pipolin(found_pipolins)):
        print(i + 1, acc, sep='\t')

        ena_xml = _retrieve_ena_xml(acc)
        _, _, _, _, _, url = _extract_metadata(ena_xml)

        file_path = os.path.join(out_dir, acc + '.fasta.gz')
        retrieve_fasta_file(url, file_path)
        unzip_fasta_file(file_path)


if __name__ == '__main__':
    download_found_pipolins()
