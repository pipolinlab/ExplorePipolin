import datetime
import gzip
import logging
import os
import shutil
import subprocess
import asyncio
import sys
import tempfile

import urllib.request
from typing import Tuple, Iterator, Optional
from urllib.error import URLError
import xml.etree.ElementTree as ElementTree

import click

from explore_pipolin.common import CONTEXT_SETTINGS


def yield_acc_and_url(ena_xml: str) -> Iterator[Tuple[str, Optional[str]]]:
    for event, elem in ElementTree.iterparse(ena_xml, events=('start',)):
        asmbl_acc, asmbl_url = None, None
        if event == 'start' and elem.tag == 'ASSEMBLY':
            asmbl_acc = elem.attrib['accession']

            for link in elem.findall('ASSEMBLY_LINKS/ASSEMBLY_LINK/URL_LINK[LABEL="WGS_SET_FASTA"]/URL'):
                asmbl_url = link.text

            if asmbl_acc is not None:
                yield asmbl_acc, asmbl_url


def unzip_genome(genome_zip: str, new_genome: str) -> None:
    with gzip.open(genome_zip, 'rb') as inf, open(new_genome, 'wb') as ouf:
        shutil.copyfileobj(inf, ouf)


_CHECKED_LIST = 'checked.txt'
_FOUND_PIPOLINS_LIST = 'found_pipolins.txt'


def _is_in_list(acc: str, list_file: str) -> bool:
    try:
        with open(list_file) as inf:
            for line in inf:
                if line.strip() == acc:
                    return True
        return False
    except FileNotFoundError:
        return False


def _update_list(acc: str, list_file: str) -> None:
    with open(list_file, 'a') as ouf:
        print(acc, file=ouf)


def _is_found(acc: str, out_dir: str) -> bool:
    log_path = os.path.join(out_dir, acc, acc + '.log')
    if not os.path.exists(log_path):
        raise AssertionError('Log should exist!')
    with open(log_path) as inf:
        return 'No piPolBs were found!' not in inf.read()


def _clean_all(acc: str, out_dir: str) -> None:
    shutil.rmtree(os.path.join(out_dir, acc))


async def analyse_genome(genome: str, out_dir: str) -> None:
    proc = await asyncio.subprocess.create_subprocess_shell(
        f'explore_pipolin --out-dir {out_dir} --no-annotation {genome}',
        stdout=subprocess.DEVNULL)
    await proc.wait()


async def download_and_analyse(acc: str, url: Optional[str], out_dir: str, sem: asyncio.BoundedSemaphore) -> None:
    try:
        if _is_in_list(acc, os.path.join(out_dir, _CHECKED_LIST)):
            return
        await do_download_and_analyse(acc, url, out_dir)
    except URLError:
        logging.info(f'Broken URL for {acc}. Skip.')
        _update_list(acc, os.path.join(out_dir, _CHECKED_LIST))
    except Exception as e:
        logging.exception(str(e))
        raise
    finally:
        sem.release()


async def do_download_and_analyse(acc: str, url: Optional[str], out_dir: str) -> None:
    if url is not None:
        with tempfile.TemporaryDirectory() as tmp:
            genome = os.path.join(tmp, acc + '.fasta')
            download_genome(url, genome)
            await analyse_genome(genome, tmp)
            logging.info(f'Finished analysis for {acc}')

            if _is_found(acc, tmp):
                logging.info(f'Pipolin(s) were found for {acc}')
                # it might be added in the list and copied to out_dir previously,
                # right before a Keyboard Interrupt signal.
                # For that reason, we check first if it's in the list and set dirs_exist_ok.
                if not _is_in_list(acc, os.path.join(out_dir, _FOUND_PIPOLINS_LIST)):
                    _update_list(acc, os.path.join(out_dir, _FOUND_PIPOLINS_LIST))
                shutil.copytree(os.path.join(tmp, acc), os.path.join(out_dir, acc), dirs_exist_ok=True)

    _update_list(acc, os.path.join(out_dir, _CHECKED_LIST))


def download_genome(url: str, genome: str) -> None:
    with tempfile.NamedTemporaryFile() as genome_zip:
        download_genome_to_path(url, genome_zip.name)
        unzip_genome(genome_zip.name, genome)


def download_genome_to_path(url: str, genome_path: str) -> None:
    urllib.request.urlretrieve(url, genome_path)


async def download_and_analyse_all(ena_xml: str, out_dir: str, p: int) -> None:
    sem = asyncio.BoundedSemaphore(p)
    tasks = []

    for acc, url in yield_acc_and_url(ena_xml):

        await sem.acquire()

        logging.info(f'{acc}, {url}')
        tasks.append(asyncio.create_task(download_and_analyse(acc, url, out_dir, sem)))

        if any(task.done() and task.exception() is not None for task in tasks):
            break

        tasks = [task for task in tasks if not task.done()]

    await asyncio.gather(*tasks)


def set_logging_to_file_and_stdout(log_file: str):
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    logger.addHandler(logging.FileHandler(log_file))
    logger.addHandler(logging.StreamHandler(sys.stdout))


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

    log_file = datetime.datetime.now().strftime('%H%M%S') + '.log'
    set_logging_to_file_and_stdout(os.path.join(out_dir, log_file))

    asyncio.run(download_and_analyse_all(ena_xml, out_dir, p))


if __name__ == '__main__':
    massive_screening()
