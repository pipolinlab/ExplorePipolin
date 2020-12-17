import logging
import os
import subprocess
from time import sleep


def check_blast():
    try:
        subprocess.run(['blastn', '-version'], stdout=subprocess.DEVNULL)
        subprocess.run(['tblastn', '-version'], stdout=subprocess.DEVNULL)
    except FileNotFoundError:
        logging.fatal('Cannot find blastn and tblastn executables in your PATH! Is BLAST+ installed?')
        exit(1)


def check_aragorn():
    try:
        subprocess.run(['aragorn', '-h'], stdout=subprocess.DEVNULL)
    except FileNotFoundError:
        logging.fatal('Cannot find aragorn executable in your PATH! Is ARAGORN installed?')
        exit(1)


def check_prokka():
    try:
        subprocess.run(['prokka', '--version'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except FileNotFoundError:
        logging.fatal('Cannot find prokka executable in your PATH! Is Prokka installed?')
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


def blast_for_repeats(genome_id, repeats_dir):
    with open(os.path.join(repeats_dir, genome_id + '.fmt5'), 'w') as ouf:
        subprocess.run(['blastn', '-query', os.path.join(repeats_dir, genome_id + '.left'),
                        '-subject', os.path.join(repeats_dir, genome_id + '.right'),
                        '-outfmt', '5', '-perc_identity', '90', '-word_size', '6',
                        '-strand', 'plus'], stdout=ouf)


def run_aragorn(genome_file, output_file):
    with open(output_file, 'w') as ouf:
        subprocess.run(['aragorn', '-l', '-w', genome_file], stdout=ouf)


def tblastn_against_ref_pipolb(genome_file, ref_pipolb, output_file):
    with open(output_file, 'w') as ouf:
        subprocess.run(['tblastn', '-query', ref_pipolb, '-subject', genome_file, '-evalue', '0.1', '-outfmt', '5'],
                       stdout=ouf)


def blastn_against_ref_att(genome_file, ref_att, output_file):
    with open(output_file, 'w') as ouf:
        subprocess.run(['blastn', '-query', ref_att, '-subject', genome_file, '-evalue', '0.1', '-outfmt', '5'],
                       stdout=ouf)


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
