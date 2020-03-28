#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import click
import os
from prefect import task
import subprocess
from utilities import CONTEXT_SETTINGS


@task
def predict_atts_with_hmmer(hmm, pipolins_dir, out_dir):
    pipolin_regions = os.listdir(pipolins_dir)
    os.makedirs(os.path.join(out_dir, 'att_hmmer'), exist_ok=True)
    for pipolin in pipolin_regions:
        genome_id = pipolin.split(sep="-")[0]
        path = os.path.join(pipolins_dir, pipolin)
        subprocess.run(['nhmmscan', '--tblout', f'{os.path.join(out_dir, "att_hmmer", genome_id)}-atts.tbl',
                        '-E', '1e-10', hmm, path])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('hmm', type=click.Path(exists=True))
@click.argument('pipolins-dir', type=click.Path(exists=True))
@click.argument('out-dir', type=click.Path())
def main(hmm, pipolins_dir, out_dir):
    """
    TODO
    """
    predict_atts_with_hmmer(hmm, pipolins_dir, out_dir)


if __name__ == '__main__':
    main()
