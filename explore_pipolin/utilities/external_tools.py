import logging
import os
import subprocess
from time import sleep
from typing import List

import pkg_resources
from Bio import SeqIO


def check_external_dependencies():
    check_blast()
    check_aragorn()
    check_prokka()
    check_prodigal()
    check_hmmsearch()


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
    p = subprocess.Popen(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    std_out, std_err = p.communicate()
    if p.returncode != 0:
        print(std_err.decode()) if std_err != b'' else print(std_out.decode())
        logging.fatal(f'Cannot check {exec_name} executable!')
        exit(1)


def run_prodigal(genome_file: str, output_file: str):
    mode = 'single' if _is_long_enough(genome_file) else 'meta'
    command = ['prodigal', '-c', '-m', '-q', '-p', mode, '-a', output_file, '-i', genome_file]
    p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    std_out, std_err = p.communicate()
    if p.returncode != 0:
        print(std_err.decode()) if std_err != b'' else print(std_out.decode())
        raise subprocess.CalledProcessError(p.returncode, ' '.join(command))


def _is_long_enough(genome_file) -> bool:
    records = SeqIO.parse(handle=genome_file, format='fasta')
    length = 0
    for record in records:
        length += len(record.seq)
    return length >= 100000


_PIPOLB_HMM_PROFILE = pkg_resources.resource_filename('explore_pipolin', 'data/pipolb_expanded_definitive.hmm')


def run_hmmsearch(proteins_file: str, pipolb_hmm_profile, output_file: str):
    if pipolb_hmm_profile is None:
        pipolb_hmm_profile = _PIPOLB_HMM_PROFILE
    command = ['hmmsearch', '--tblout', output_file, '-E', '1e-50', pipolb_hmm_profile, proteins_file]
    p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    std_out, std_err = p.communicate()
    if p.returncode != 0:
        print(std_err.decode()) if std_err != b'' else print(std_out.decode())
        raise subprocess.CalledProcessError(p.returncode, ' '.join(command))


_REF_ATT = pkg_resources.resource_filename('explore_pipolin', 'data/attL.fa')


def blastn_against_ref_att(genome_file, output_file, ref_att):
    if ref_att is None:
        ref_att = _REF_ATT
    with open(output_file, 'w') as ouf:
        command = ['blastn', '-query', ref_att, '-subject', genome_file, '-evalue', '0.1', '-outfmt', '5']
        p = subprocess.Popen(command, stdout=ouf, stderr=subprocess.PIPE)
        _, std_err = p.communicate()
        if p.returncode != 0:
            if std_err != b'':
                print(std_err.decode())
            raise subprocess.CalledProcessError(p.returncode, ' '.join(command))


def blast_for_repeats(genome_id: str, repeats_dir: str, perc_identity):
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
            command = ['blastn', '-query', os.path.join(repeats_dir, lr),
                       '-subject', os.path.join(repeats_dir, rr),
                       '-outfmt', '5', '-perc_identity', str(perc_identity),
                       '-word_size', '6', '-strand', 'plus']
            p = subprocess.Popen(command, stdout=ouf, stderr=subprocess.PIPE)
            _, std_err = p.communicate()
            if p.returncode != 0:
                if std_err != b'':
                    print(std_err)
                raise subprocess.CalledProcessError(p.returncode, ' '.join(command))


def _extract_index_from_filename(filename: str) -> int:
    return int(os.path.splitext(filename)[0].split(sep='_')[-1])


def run_aragorn(genome_file, output_file):
    with open(output_file, 'w') as ouf:
        p = subprocess.run(['aragorn', '-l', '-w', genome_file], stdout=ouf, stderr=subprocess.PIPE)
        # If no file, still returns 0 exitcode. Therefore, checking for the message.
        if b'Could not open ' in p.stderr:
            raise FileNotFoundError(genome_file)


_PROTEINS = pkg_resources.resource_filename('explore_pipolin', '/data/HHpred_proteins.faa')


def run_prokka(input_file, prokka_results_dir, proteins, cpus):
    if proteins is None:
        proteins = _PROTEINS
    prefix = os.path.splitext(os.path.basename(input_file))[0]
    command = ['prokka', '--outdir', prokka_results_dir, '--prefix', prefix,
               '--rawproduct', '--cdsrnaolap', '--cpus', str(cpus),
               '--rfam', '--proteins', proteins, '--force', input_file]
    p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    std_out, std_err = p.communicate()
    if p.returncode != 0:
        print(std_err.decode()) if std_err != b'' else print(std_out.decode())
        raise subprocess.CalledProcessError(p.returncode, ' '.join(command))


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
