#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os
import click
import subprocess
from utilities import CONTEXT_SETTINGS


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('pipolins-dir', type=click.Path(exists=True))
@click.argument('proteins', type=click.Path(exists=True))
@click.argument('out-dir')
def annotate_pipolins(pipolins_dir, proteins, out_dir):
    """
    TODO
    """
    pipolins = os.listdir(pipolins_dir)
    for i_p, pipolin in enumerate(pipolins):
        print(i_p)
        genome_id = pipolin.split(sep="-")[0]
        subprocess.run(['/home/liubov/repos/prokka/bin/prokka', '--outdir', out_dir,
                    '--prefix', f'{genome_id}', '--rawproduct', '--cdsrnaolap', '--cpus', '4',
                    '--rfam', '--proteins', proteins, '--force',
                    '--locustag', f'{genome_id}', os.path.join(pipolins_dir, pipolin)])


if __name__ == '__main__':
    annotate_pipolins()
