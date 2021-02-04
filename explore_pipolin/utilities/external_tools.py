import logging
import os
import subprocess
from time import sleep
from typing import List

import pkg_resources
from Bio import SeqIO


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


def run_prodigal(genome_file: str, output_file: str):
    mode = 'single' if _is_long_enough(genome_file) else 'meta'
    subprocess.check_call(['prodigal', '-c', '-m', '-q', '-p', mode, '-a', output_file, '-i', genome_file],
                          stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)


def _is_long_enough(genome_file) -> bool:
    records = SeqIO.parse(handle=genome_file, format='fasta')
    length = 0
    for record in records:
        length += len(record.seq)
    return length >= 100000


_PIPOLB_HMM_PROFILE = pkg_resources.resource_filename('explore_pipolin', 'data/pipolb_expanded_definitive.hmm')


def run_hmmsearch(proteins_file: str, output_file: str):
    subprocess.check_call(['hmmsearch', '--tblout', output_file, '-E', '1e-50', _PIPOLB_HMM_PROFILE, proteins_file],
                          stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)


_REF_ATT = pkg_resources.resource_filename('explore_pipolin', 'data/attL.fa')


def blastn_against_ref_att(genome_file, output_file):
    with open(output_file, 'w') as ouf:
        subprocess.check_call(['blastn', '-query', _REF_ATT, '-subject', genome_file, '-evalue', '0.1', '-outfmt', '5'],
                              stdout=ouf, stderr=subprocess.PIPE)


def blast_for_repeats(genome_id: str, repeats_dir: str):
    dir_content = os.listdir(repeats_dir)
    left_sorted = sorted([f for f in dir_content if f.startswith(genome_id) and f.endswith('.left')])
    right_sorted = sorted([f for f in dir_content if f.startswith(genome_id) and f.endswith('.right')])
    if len(left_sorted) != len(right_sorted):
        raise AssertionError('Number of .left files is not equal to the number of .right files!')

    for lr, rr in zip(left_sorted, right_sorted):
        l_rep_index = _extract_index_from_filename(lr)
        r_rep_index = _extract_index_from_filename(rr)
        if l_rep_index != r_rep_index:
            raise AssertionError(f'Wrong pair for file {lr}: {rr}')

        with open(os.path.join(repeats_dir, genome_id + f'_{l_rep_index}.fmt5'), 'w') as ouf:
            subprocess.check_call(['blastn', '-query', os.path.join(repeats_dir, lr),
                                   '-subject', os.path.join(repeats_dir, rr),
                                   '-outfmt', '5', '-perc_identity', '85', '-word_size', '6',
                                   '-strand', 'plus'], stdout=ouf, stderr=subprocess.PIPE)


def _extract_index_from_filename(filename: str) -> int:
    return int(os.path.splitext(filename)[0].split(sep='_')[-1])


def run_aragorn(genome_file, output_file):
    with open(output_file, 'w') as ouf:
        p = subprocess.run(['aragorn', '-l', '-w', genome_file], stdout=ouf, stderr=subprocess.PIPE)
        # If no file, still returns 0 exitcode. Therefore, checking for the message.
        if b'Could not open ' in p.stderr:
            raise FileNotFoundError(genome_file)


_PROTEINS = pkg_resources.resource_filename('explore_pipolin', '/data/HHpred_proteins.faa')


def run_prokka(input_file, prokka_results_dir):
    prefix = os.path.splitext(os.path.basename(input_file))[0]
    subprocess.check_call(['prokka', '--outdir', prokka_results_dir, '--prefix', prefix,
                          # TODO: number of CPUs is hardcoded. To pass it as an argument?
                           '--rawproduct', '--cdsrnaolap', '--cpus', '4',
                           '--rfam', '--proteins', _PROTEINS, '--force', input_file],
                          stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)


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
            if proc.stdout is not None and b'Empty result - nothing to do' in proc.stdout:
                raise EmptyResult
            if proc.stdout is not None and b'Query failed on MegaLink server' in proc.stdout:
                raise NoAssembly
            print('FAILED to retrieve the data! Retrying ...')
            sleep(sleep_time)
            sleep_time = sleep_time * 2
            continue

    print('FAILED!!! Maximum number of retries is exceeded.')
