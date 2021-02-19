import gzip
import os
import shutil
import sys

import requests
import urllib.request
import json
import xml.etree.ElementTree as XMLParser

import click
from Bio import SeqIO

from explore_pipolin.common import CONTEXT_SETTINGS
from explore_pipolin.utilities.io import create_seqio_records_dict


def yield_accessions(json_file):
    with open(json_file) as inf:
        content = inf.readline()
    for entry in json.loads(content):
        yield entry['assembly_accession']


def retrieve_ena_data(acc):
    server = "https://www.ebi.ac.uk/ena/browser/api/xml/" + acc
    r = requests.get(server)

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    return r.text


def get_fasta_url(ena_data):
    root = XMLParser.fromstring(ena_data)
    for ass in root.findall('ASSEMBLY'):
        for link in ass.find('ASSEMBLY_LINKS').findall('ASSEMBLY_LINK'):
            url_link = link.find('URL_LINK')
            if url_link.find('LABEL').text == 'WGS_SET_FASTA':
                return url_link.find('URL').text


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
@click.argument('metadata-file', type=click.Path(exists=True))
@click.argument('output-dir', type=click.Path())
@click.option('--force', is_flag=True,
              help='If the output directory already exists, erase its content before starting the analysis.')
def download_and_analyse(metadata_file, output_dir, force):
    if force:
        os.remove(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    for acc in yield_accessions(metadata_file):
        print(acc)

        if _is_analysed(acc, output_dir):
            continue

        ena_data = retrieve_ena_data(acc)
        fasta_url = get_fasta_url(ena_data)
        if fasta_url is not None:
            file_path = os.path.join(output_dir, acc + '.fasta.gz')
            retrieve_fasta_file(fasta_url, file_path)
            file_to_analyse = unzip_fasta_file(file_path)
            change_fasta_identifiers(file_to_analyse)
            os.system(f'explore_pipolin --out-dir {output_dir} --add-colours {file_to_analyse}')
            os.remove(file_to_analyse)

        _update_analysed(acc, output_dir)


if __name__ == '__main__':
    download_and_analyse()
