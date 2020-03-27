#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os
import click
import subprocess
from utilities import CONTEXT_SETTINGS


def annotate_pipolins(out_dir, pipolins_dir, proteins):
    pipolins = os.listdir(pipolins_dir)
    for i_p, pipolin in enumerate(pipolins):
        print(i_p)
        genome_id = pipolin.split(sep="-")[0]
        subprocess.run(['/home/liubov/repos/prokka/bin/prokka', '--outdir', out_dir,
                        '--prefix', f'{genome_id}', '--rawproduct', '--cdsrnaolap', '--cpus', '4',
                        '--rfam', '--proteins', proteins, '--force',
                        '--locustag', f'{genome_id}', os.path.join(pipolins_dir, pipolin)])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('pipolins-dir', type=click.Path(exists=True))
@click.argument('proteins', type=click.Path(exists=True))
@click.argument('out-dir')
def main(pipolins_dir, proteins, out_dir):
    """
    TODO
    ~23 min and 1116 files for 93 genomes.
    """
    annotate_pipolins(out_dir, pipolins_dir, proteins)


if __name__ == '__main__':
    main()
