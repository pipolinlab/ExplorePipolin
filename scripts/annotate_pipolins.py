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

    pipolins = os.listdir(pipolins_dir)
    for i_p, pipolin in enumerate(pipolins):
        print(i_p)
        genome_id = pipolin.split(sep="-")[0]
        subprocess.run(['/home/liubov/repos/prokka/bin/prokka', '--outdir', 'prokka_annotations',
                    '--prefix', f'{genome_id}', '--rawproduct', '--cdsrnaolap', '--cpus', '4',   # '--rfam'
                    '--locustag', f'{genome_id}', '--force', f'./pipolin_regions/{genome_id}-pipolin.fa'])


if __name__ == '__main__':
    annotate_pipolins()