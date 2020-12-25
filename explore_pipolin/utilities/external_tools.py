import logging
import os
import subprocess
from time import sleep
from typing import List

import pkg_resources
from Bio import SeqIO

from explore_pipolin.common import Window

_PIPOLB_HMM_PROFILE = pkg_resources.resource_filename('explore_pipolin', 'data/pipolb_expanded_definitive.hmm')


def check_blast():
    command_args = ['blastn', '-version']
    _try_command_except_fatal(command_args, 'blastn')


def check_aragorn():
    command_args = ['aragorn', '-h']
    _try_command_except_fatal(command_args, 'aragorn')


def check_prokka():
    command_args = ['prokka', '--version']
    _try_command_except_fatal(command_args, 'prokka')


def check_prodigal():
    command_args = ['prodigal', '-h']
    _try_command_except_fatal(command_args, 'prodigal')


def check_hmmsearch():
    command_args = ['hmmsearch', '-h']
    _try_command_except_fatal(command_args, 'hmmsearch')


def _try_command_except_fatal(command_args: List[str], exec_name: str):
    try:
        subprocess.run(command_args, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except FileNotFoundError:
        logging.fatal(f'Cannot call {exec_name} executable! Is it installed and added to your PATH?')
        exit(1)


def run_prokka(genome_id, pipolins_dir, proteins, prokka_results_dir):
    subprocess.run(['prokka', '--outdir', prokka_results_dir,
                    '--prefix', genome_id,
                    # TODO: number of CPUs is hardcoded. To pass it as an argument?
                    '--rawproduct', '--cdsrnaolap', '--cpus', '4',
                    '--rfam', '--proteins', proteins, '--force',
                    '--locustag', genome_id,
                    os.path.join(pipolins_dir, genome_id + '.fa')],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def blast_for_repeats(windows: List[Window], repeats_dir):
    genome_id = windows[0].pipolbs[0].genome.id
    for i in range(len(windows)):
        with open(os.path.join(repeats_dir, genome_id + f'_{i}.fmt5'), 'w') as ouf:
            subprocess.run(['blastn', '-query', os.path.join(repeats_dir, genome_id + f'_{i}.left'),
                            '-subject', os.path.join(repeats_dir, genome_id + f'_{i}.right'),
                            '-outfmt', '5', '-perc_identity', '90', '-word_size', '6',
                            '-strand', 'plus'], stdout=ouf)


def run_aragorn(genome_file, output_file):
    with open(output_file, 'w') as ouf:
        subprocess.run(['aragorn', '-l', '-w', genome_file], stdout=ouf)


def blastn_against_ref_att(genome_file, ref_att, output_file):
    with open(output_file, 'w') as ouf:
        subprocess.run(['blastn', '-query', ref_att, '-subject', genome_file, '-evalue', '0.1', '-outfmt', '5'],
                       stdout=ouf)


def run_prodigal(genome_file: str, output_file: str):
    mode = 'single' if _is_long_enough(genome_file) else 'meta'
    subprocess.run(['prodigal', '-c', '-m', '-q', '-p', mode, '-a', output_file, '-i', genome_file],
                   stdout=subprocess.DEVNULL)


def _is_long_enough(genome_file):
    records = SeqIO.parse(handle=genome_file, format='fasta')
    length = 0
    for record in records:
        length += len(record.seq)
    return True if length >= 100000 else False


def run_hmmsearch(proteins_file: str, output_file: str):
    subprocess.run(['hmmsearch', '--tblout', output_file, '-E', '0.01', _PIPOLB_HMM_PROFILE, proteins_file],
                   stdout=subprocess.DEVNULL)


class EmptyResult(Exception):
    pass


class NoAssembly(Exception):
    pass


def subprocess_with_retries(*args, **kwargs):
    num_retries = 5
    sleep_time = 0.5
    for i in range(num_retries):
        proc = subprocess.run(*args, **kwargs)

        try:
            proc.check_returncode()
            return proc
        except subprocess.CalledProcessError:
            if proc.stdout is not None and 'Empty result - nothing to do' in proc.stdout:
                raise EmptyResult
            if proc.stdout is not None and 'Query failed on MegaLink server' in proc.stdout:
                raise NoAssembly
            print('FAILED to retrieve the data! Retrying ...')
            sleep(sleep_time)
            sleep_time = sleep_time * 2
            continue

    print('FAILED!!! Maximum number of retries is exceeded.')
