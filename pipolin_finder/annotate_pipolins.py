#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os
import click
from prefect import task
import subprocess
from pipolin_finder.utilities import CONTEXT_SETTINGS


@task
def annotate_pipolins(gquery, pipolins_dir, proteins, root_dir):
    prokka_dir = os.path.join(root_dir, 'prokka')
    os.makedirs(prokka_dir, exist_ok=True)

    subprocess.run(['prokka', '--outdir', prokka_dir,
                    '--prefix', gquery.gquery_id,
                    '--rawproduct', '--cdsrnaolap', '--cpus', '4',
                    '--rfam', '--proteins', proteins, '--force',
                    '--locustag', gquery.gquery_id,
                    os.path.join(pipolins_dir, gquery.gquery_id + '.fa')])

    return prokka_dir


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
