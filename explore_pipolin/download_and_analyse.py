import gzip
import os
import shutil

import urllib.request
from urllib.error import URLError
import xml.etree.ElementTree as XMLParser

import click
from Bio import SeqIO

from explore_pipolin.common import CONTEXT_SETTINGS
from explore_pipolin.utilities.io import create_seqio_records_dict


def yield_acc_and_url(ena_xml_file):
    for event, elem in XMLParser.iterparse(ena_xml_file, events=('start',)):
        accession, url = None, None
        if event == 'start' and elem.tag == 'ASSEMBLY':
            accession = elem.attrib['accession']

            if elem.find('ASSEMBLY_LINKS') is not None:
                for link in elem.find('ASSEMBLY_LINKS').findall('ASSEMBLY_LINK'):
                    url_link = link.find('URL_LINK')
                    if url_link is not None:
                        if url_link.find('LABEL') is not None and url_link.find('LABEL').text == 'WGS_SET_FASTA':
                            if url_link.find('URL') is not None:
                                url = url_link.find('URL').text

        if accession is not None and url is not None:
            yield accession, url


def retrieve_fasta_file(fasta_url, file_path):
    urllib.request.urlretrieve(fasta_url, file_path)


def unzip_fasta_file(file_path):
    new_file = os.path.splitext(file_path)[0]
    with gzip.open(file_path, 'rb') as inf, open(new_file, 'wb') as ouf:
        shutil.copyfileobj(inf, ouf)
    os.remove(file_path)
    return new_file


def change_fasta_identifiers(genome_file):
    """This is valid only for identifiers of assemblies from the ENA database"""
    genome_dict = create_seqio_records_dict(genome_file, 'fasta')
    new_dict = {}
    for key, value in genome_dict.items():
        new_id = key.split(sep='|')[-1]
        value.id = new_id
        new_dict[key] = value
    with open(genome_file, 'w') as ouf:
        SeqIO.write(new_dict.values(), ouf, 'fasta')


def _is_analysed(acc, output_dir):
    try:
        with open(os.path.join(output_dir, 'analysed.txt')) as inf:
            for line in inf:
                if line.strip() == acc:
                    return True
        return False
    except FileNotFoundError:
        return False


def _update_analysed(acc, output_dir):
    with open(os.path.join(output_dir, 'analysed.txt'), 'a') as ouf:
        print(acc, file=ouf)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('ena-xml-file', type=click.Path(exists=True))
@click.argument('output-dir', type=click.Path())
@click.option('--force', is_flag=True,
              help='If the output directory already exists, erase its content before starting the analysis.')
def download_and_analyse(ena_xml_file, output_dir, force):
    if force:
        os.remove(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    for acc, url in yield_acc_and_url(ena_xml_file):
        print(acc, url)

        if _is_analysed(acc, output_dir):
            continue

        file_path = os.path.join(output_dir, acc + '.fasta.gz')

        try:
            retrieve_fasta_file(url, file_path)
        except URLError:
            if not _is_analysed(acc, output_dir):
                _update_analysed(acc, output_dir)
            continue

        file_to_analyse = unzip_fasta_file(file_path)
        change_fasta_identifiers(file_to_analyse)
        os.system(f'explore_pipolin --out-dir {output_dir} --add-colours {file_to_analyse}')
        os.remove(file_to_analyse)

        _update_analysed(acc, output_dir)


if __name__ == '__main__':
    download_and_analyse()
