#!/usr/bin/evn python3
# -*- encoding: utf-8 -*-

import click
import os
import subprocess
from utilities import CONTEXT_SETTINGS
from utilities import check_dir

ABRICATE_PATH = '/home/liubov/repos/abricate/bin/abricate'


def run_abricate(databases, genomes, in_dir, out_dir):
    check_dir(out_dir)
    for database in databases:
        database_dir = os.path.join(out_dir, database)
        os.mkdir(database_dir)
        for genome in genomes:
            with open(os.path.join(database_dir, f'{genome[:-3]}.tab', ), 'w') as ouf:
                subprocess.run([ABRICATE_PATH, '--db', f'{database}', '--threads', '4',
                                '--nopath', f'{os.path.join(in_dir, genome)}'], stdout=ouf)


def get_abricate_databases():
    abricate_list = subprocess.run([ABRICATE_PATH, '--list'],
                                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.decode('utf-8')
    databases = []
    for line in abricate_list.split('\n'):
        databases.append(line.split('\t')[0])
    databases.remove('DATABASE')
    databases.remove('')
    databases.remove('RIP-DB')   # run later
    databases.remove('REL-DB')   # run later
    return databases


def prepare_summaries(databases, out_dir):
    check_dir(os.path.join(out_dir, 'summaries'))
    for db in databases:
        database_path = os.path.join(out_dir, db)
        db_files = os.listdir(database_path)
        with open(os.path.join(out_dir, 'summaries', f'{db}.tab'), 'w') as ouf:
            subprocess.run(f'{ABRICATE_PATH} --nopath --summary '
                           f'{" ".join([os.path.join(database_path, file) for file in db_files])}',
                           shell=True, stdout=ouf)


def parse_summaries(databases, out_dir):
    summaries_data = {}
    for db in databases:
        summaries_data[db] = {}
        with open(os.path.join(out_dir, 'summaries', f'{db}.tab')) as inf:
            for line in inf:
                if line[0] != '#':
                    genome, number = line.strip().split(sep='\t')[:2]
                    summaries_data[db][genome] = number

    return summaries_data


def save_summaries_data_to_csv(genomes, databases, summaries_data, out_dir):
    with open(os.path.join(out_dir, 'short_summary.csv'), 'w') as ouf:
        print(f'FILE,{",".join(databases)}', file=ouf)
        for genome in genomes:
            numbers = ','.join([summaries_data[db][f'{genome[:-3]}.tab'] for db in databases])
            print(f'{genome},{numbers}', file=ouf)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('in-dir', type=click.Path(exists=True))
@click.argument('out-dir', type=click.Path())
def abricate_analysis(in_dir, out_dir):
    """
    TODO
    """
    databases = get_abricate_databases()
    genomes = os.listdir(in_dir)

    # run_abricate(databases, genomes, in_dir, out_dir)
    prepare_summaries(databases, out_dir)
    summaries_data = parse_summaries(databases, out_dir)
    save_summaries_data_to_csv(genomes, databases, summaries_data, out_dir)


if __name__ == '__main__':
    abricate_analysis()
