#!/usr/bin/evn python3
# -*- encoding: utf-8 -*-

import click
import os
import subprocess
from utilities import CONTEXT_SETTINGS
from utilities import check_dir


def run_abricate(databases, genomes, in_dir, out_dir):
    check_dir(out_dir)
    for database in databases:
        database_dir = os.path.join(out_dir, database)
        os.mkdir(database_dir)
        for genome in genomes:
            with open(os.path.join(database_dir, f'{genome[:-3]}.tab', ), 'w') as ouf:
                subprocess.run(['/home/liubov/repos/abricate/bin/abricate', '--db', f'{database}', '--threads', '4',
                                '--nopath', f'{os.path.join(in_dir, genome)}'], stdout=ouf)


def get_abricate_databases():
    abricate_list = subprocess.run(['/home/liubov/repos/abricate/bin/abricate', '--list'],
                                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.decode('utf-8')
    databases = []
    for line in abricate_list.split('\n'):
        databases.append(line.split('\t')[0])
    databases.remove('DATABASE')
    databases.remove('')
    databases.remove('RIP-DB')   # TODO: some problem with protein DBs
    databases.remove('REL-DB')
    return databases


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('in-dir', type=click.Path(exists=True))
@click.argument('out-dir', type=click.Path())
def abricate_analysis(in_dir, out_dir):
    """
    TODO
    """
    databases = get_abricate_databases()
    genomes = os.listdir(in_dir)

    run_abricate(databases, genomes, in_dir, out_dir)


if __name__ == '__main__':
    abricate_analysis()
