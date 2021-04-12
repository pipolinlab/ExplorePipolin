import gzip
import os
import shutil

import urllib.request
from urllib.error import URLError
import xml.etree.ElementTree as ElementTree

import click
from Bio import SeqIO

from explore_pipolin.common import CONTEXT_SETTINGS
from explore_pipolin.utilities.io import create_seqio_records_dict


def yield_acc_and_url(ena_xml_file):
    for event, elem in ElementTree.iterparse(ena_xml_file, events=('start',)):
        asmbl_acc, asmbl_url = None, None
        if event == 'start' and elem.tag == 'ASSEMBLY':
            asmbl_acc = elem.attrib['accession']

            for link in elem.findall('ASSEMBLY_LINKS/ASSEMBLY_LINK/URL_LINK[LABEL="WGS_SET_FASTA"]/URL'):
                asmbl_url = link.text

        if asmbl_acc is not None and asmbl_url is not None:
            yield asmbl_acc, asmbl_url


def retrieve_fasta_file(fasta_url, file_path):
    urllib.request.urlretrieve(fasta_url, file_path)


def unzip_fasta_file(file_path):
    new_file = os.path.splitext(file_path)[0]
    with gzip.open(file_path, 'rb') as inf, open(new_file, 'wb') as ouf:
        shutil.copyfileobj(inf, ouf)
    os.remove(file_path)
    return new_file


class SkipAssemblyError(Exception):
    pass


def change_fasta_identifiers(genome_file):
    """This is valid only for identifiers of assemblies from the ENA database"""
    genome_dict = create_seqio_records_dict(genome_file, 'fasta')
    new_dict = {}
    for key, value in genome_dict.items():
        # TODO: here we need to exclude metagenome assemblies somehow. CHANGE THIS!
        if 'Enterobacter' not in value.description:
            raise SkipAssemblyError()

        new_id = key.split(sep='|')[-1]
        if len(new_id) > 14:
            new_id = new_id[-14:]
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


def update_if_not_analysed(acc, output_dir):
    if not _is_analysed(acc, output_dir):
        _update_analysed(acc, output_dir)


def _download_and_analyse(output_dir, acc, url) -> None:
    file_path = os.path.join(output_dir, acc + '.fasta.gz')

    try:
        retrieve_fasta_file(url, file_path)
    except URLError:
        update_if_not_analysed(acc, output_dir)
        return

    file_to_analyse = unzip_fasta_file(file_path)

    try:
        change_fasta_identifiers(file_to_analyse)
    except SkipAssemblyError:
        os.remove(file_to_analyse)
        update_if_not_analysed(acc, output_dir)
        return

    os.system(f'explore_pipolin --out-dir {output_dir} --add-colours {file_to_analyse}')
    os.remove(file_to_analyse)
    update_if_not_analysed(acc, output_dir)


def clean_if_not_found(acc, output_dir):
    log_path = os.path.join(output_dir, 'logs', acc + '.log')
    if os.path.exists(log_path):
        with open(log_path) as inf:
            not_found = False
            if 'No piPolBs were found!' in inf.read():
                not_found = True
                try:
                    os.remove(os.path.join(output_dir, 'pipolbs', acc + '.faa'))
                    os.remove(os.path.join(output_dir, 'pipolbs', acc + '.tbl'))
                except OSError:
                    pass
        if not_found:
            os.remove(os.path.join(output_dir, 'logs', acc + '.log'))


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

        _download_and_analyse(output_dir, acc, url)
        clean_if_not_found(acc, output_dir)


if __name__ == '__main__':
    download_and_analyse()
