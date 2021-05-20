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
from concurrent.futures import ThreadPoolExecutor
from typing import Optional
from urllib.error import URLError

import click

from explore_pipolin.common import CONTEXT_SETTINGS
from explore_pipolin.collect_metadata import retrieve_ena_xml, extract_metadata


def yield_accession(accessions):
    with open(accessions) as inf:
        for line in inf:
            yield line.strip()


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


async def download_genome_to_path(url: str, genome_path: str) -> None:
    await asyncio.get_running_loop().run_in_executor(None, urllib.request.urlretrieve, url, genome_path)


def unzip_genome(genome_zip: str, new_genome: str) -> None:
    with gzip.open(genome_zip, 'rb') as inf, open(new_genome, 'wb') as ouf:
        shutil.copyfileobj(inf, ouf)


async def download_genome(url: str, genome: str) -> None:
    with tempfile.NamedTemporaryFile() as genome_zip:
        for i in range(5):
            try:
                await asyncio.wait_for(download_genome_to_path(url, genome_zip.name), timeout=5)
                break
            except (asyncio.TimeoutError, URLError) as e:
                if i == 4:
                    if e == URLError:
                        raise
                    else:
                        raise AssertionError(
                            f'Cannot download genome for {os.path.splitext(os.path.basename(genome))[0]}. '
                            f'Check the link: {url} of check your Internet connection!'
                        )
                else:
                    logging.info('Timeout! Trying again...')
                    continue

        unzip_genome(genome_zip.name, genome)


async def analyse_genome(genome: str, out_dir: str) -> None:
    proc = await asyncio.subprocess.create_subprocess_shell(
        f'explore_pipolin --out-dir {out_dir} --no-annotation {genome}',
        stdout=subprocess.DEVNULL)
    await proc.wait()


async def do_download_and_analyse(acc: str, url: Optional[str], out_dir: str) -> None:
    if url is not None:
        with tempfile.TemporaryDirectory() as tmp:
            genome = os.path.join(tmp, acc + '.fasta')
            await download_genome(url, genome)
            logging.info(f'Started analysing {acc}')
            await analyse_genome(genome, tmp)
            logging.info(f'Finished analysing {acc}')

            if _is_found(acc, tmp):
                logging.info(f'Pipolin(s) were found for {acc}')
                # it might be added in the list and copied to out_dir previously,
                # right before a Keyboard Interrupt signal.
                # For that reason, we check first if it's in the list and set dirs_exist_ok.
                if not _is_in_list(acc, os.path.join(out_dir, _FOUND_PIPOLINS_LIST)):
                    _update_list(acc, os.path.join(out_dir, _FOUND_PIPOLINS_LIST))
                shutil.copytree(os.path.join(tmp, acc), os.path.join(out_dir, acc), dirs_exist_ok=True)

    _update_list(acc, os.path.join(out_dir, _CHECKED_LIST))


async def download_and_analyse(acc: str, url: Optional[str], out_dir: str, sem: asyncio.BoundedSemaphore) -> None:
    try:
        await do_download_and_analyse(acc, url, out_dir)
    except URLError:
        logging.info(f'Broken URL for {acc}. Skip.')
        _update_list(acc, os.path.join(out_dir, _CHECKED_LIST))
    except Exception as e:
        logging.exception(str(e))
        raise
    finally:
        sem.release()


async def download_and_analyse_all(accessions: str, out_dir: str, p: int) -> None:
    sem = asyncio.BoundedSemaphore(p)
    tasks = []

    for acc in yield_accession(accessions):

        if _is_in_list(acc, os.path.join(out_dir, _CHECKED_LIST)):
            continue

        ena_xml = retrieve_ena_xml(acc)
        _, _, _, _, _, url = extract_metadata(ena_xml)
        if url == '-':
            logging.info(f'No URL for {acc}. Skip.')
            _update_list(acc, os.path.join(out_dir, _CHECKED_LIST))
            continue

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
@click.argument('accessions', type=click.Path(exists=True))
@click.argument('out-dir', type=click.Path())
@click.option('-p', type=int, default=1, show_default=True, help='Number of processes to run.')
def main(accessions, out_dir, p):
    """
    ACCESSIONS is a file with assembly accession IDs.
    An accession of each checked genome assembly will be written to the checked.txt file.
    When the process is interrupted and rerun again, these accessions will be skipped.
    Accessions of pipolin-harboring genomes will be written to the found_pipolins.txt file.
    Stdout will be saved to the log file in the output directory.
    """
    os.makedirs(out_dir, exist_ok=True)

    log_file = datetime.datetime.now().strftime('%H%M%S') + '.log'
    set_logging_to_file_and_stdout(os.path.join(out_dir, log_file))

    asyncio.get_event_loop().set_default_executor(ThreadPoolExecutor(p))

    asyncio.run(download_and_analyse_all(accessions, out_dir, p))


if __name__ == '__main__':
    main()
