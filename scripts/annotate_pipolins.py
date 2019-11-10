#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os
import click
import subprocess
from utilities import CONTEXT_SETTINGS   # TODO: fix this!


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('pipolins-dir', type=click.Path(exists=True))
def annotate_pipolins(pipolins_dir):
    """
    TODO
    """
    os.chdir(os.path.dirname(pipolins_dir))
    if not os.path.exists('prokka_annotations'):
        os.mkdir('prokka_annotations')

    pipolins = os.listdir(pipolins_dir)
    for i_p, pipolin in enumerate(pipolins):
        print(i_p)
        genome_id = pipolin.split(sep="-")[0]
        subprocess.run(['/home/liubov/repos/prokka/bin/prokka', '--outdir', 'prokka_annotations',
                    '--prefix', f'{genome_id}', '--rfam', '--rawproduct', '--cdsrnaolap', '--cpus', '4',
                    '--locustag', f'{genome_id}', f'./pipolin_regions/{genome_id}_selected.fa'])


if __name__ == '__main__':
    annotate_pipolins()