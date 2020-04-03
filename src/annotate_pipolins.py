#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os
import click
from prefect import task
import subprocess
from utilities import CONTEXT_SETTINGS


# @task
def annotate_pipolins(pipolins, proteins, out_dir):
    os.makedirs(os.path.join(out_dir, 'prokka'), exist_ok=True)
    for i_p, pipolin in enumerate(pipolins):
        print(i_p)
        genome_id = os.path.basename(pipolin.split(sep="-")[0])
        subprocess.run(['/home/liubov/repos/prokka/bin/prokka', '--outdir', os.path.join(out_dir, 'prokka'),
                        '--prefix', f'{genome_id}', '--rawproduct', '--cdsrnaolap', '--cpus', '4',
                        '--rfam', '--proteins', proteins, '--force',
                        '--locustag', f'{genome_id}', pipolin])

    return os.path.join(out_dir, 'prokka')


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('pipolins', nargs=-1, type=click.Path(exists=True))
@click.argument('proteins', type=click.Path(exists=True))
@click.argument('out-dir')
def main(pipolins, proteins, out_dir):
    """
    TODO
    ~23 min and 1116 files for 93 genomes.
    """
    annotate_pipolins(pipolins, proteins, out_dir)


if __name__ == '__main__':
    main()
