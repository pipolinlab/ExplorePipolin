import gzip
import os
import shutil
import subprocess
import asyncio

import urllib.request
from typing import Tuple, Iterator
from urllib.error import URLError
import xml.etree.ElementTree as ElementTree

import click

from explore_pipolin.common import CONTEXT_SETTINGS


def yield_acc_and_url(ena_xml) -> Iterator[Tuple[str, str]]:
    for event, elem in ElementTree.iterparse(ena_xml, events=('start',)):
        asmbl_acc, asmbl_url = None, None
        if event == 'start' and elem.tag == 'ASSEMBLY':
            asmbl_acc = elem.attrib['accession']

            for link in elem.findall('ASSEMBLY_LINKS/ASSEMBLY_LINK/URL_LINK[LABEL="WGS_SET_FASTA"]/URL'):
                asmbl_url = link.text

            if (asmbl_acc is not None) and (asmbl_url is not None):
                yield asmbl_acc, asmbl_url


def retrieve_fasta_file(fasta_url, file_path):
    urllib.request.urlretrieve(fasta_url, file_path)


def unzip_fasta_file(file_path):
    new_file = os.path.splitext(file_path)[0]
    with gzip.open(file_path, 'rb') as inf, open(new_file, 'wb') as ouf:
        shutil.copyfileobj(inf, ouf)
    os.remove(file_path)
    return new_file


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


async def _download_and_analyse(output_dir, acc, url) -> None:
    file_path = os.path.join(output_dir, acc + '.fasta.gz')

    try:
        retrieve_fasta_file(url, file_path)
    except URLError:
        update_if_not_analysed(acc, output_dir)
        return

    file_to_analyse = unzip_fasta_file(file_path)

    proc = await asyncio.subprocess.create_subprocess_shell(
        f'explore_pipolin --out-dir {output_dir} --no-annotation {file_to_analyse}',
        stdout=subprocess.DEVNULL)
    await proc.wait()
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
            else:
                with open(os.path.join(output_dir, 'found_pipolins.txt'), 'a') as ouf:
                    print(acc, file=ouf)
        if not_found:
            os.remove(os.path.join(output_dir, 'logs', acc + '.log'))


async def download_and_analyse(acc: str, url: str, out_dir, sem: asyncio.BoundedSemaphore):
    print(acc, url)
    if not _is_analysed(acc, out_dir):
        await _download_and_analyse(out_dir, acc, url)
        clean_if_not_found(acc, out_dir)
    sem.release()


async def download_and_analyse_all(ena_xml, out_dir: str, p: int):
    sem = asyncio.BoundedSemaphore(p)
    tasks = []
    for acc, url in yield_acc_and_url(ena_xml):
        await sem.acquire()
        tasks.append(asyncio.create_task(download_and_analyse(acc, url, out_dir, sem)))
        tasks = [task for task in tasks if not task.done()]
    await asyncio.gather(*tasks)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('ena-xml', type=click.Path(exists=True))
@click.argument('out-dir', type=click.Path())
@click.option('-p', type=int, default=1, show_default=True, help='Number of processes to run.')
def massive_screening(ena_xml, out_dir, p):
    """
    ENA_XML is a file downloaded from ENA database after a search of genome assemblies
    for an organism of interest.
    An accession of each analysed genome assembly will be written to the analysed.txt file.
    When the process is interrupted and rerun again, these accessions will be skipped.
    Accessions of pipolin-harboring genomes will be written to the found_pipolins.txt file.
    """
    os.makedirs(out_dir, exist_ok=True)

    asyncio.run(download_and_analyse_all(ena_xml, out_dir, p))


if __name__ == '__main__':
    massive_screening()
