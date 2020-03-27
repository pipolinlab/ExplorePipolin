#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import click
import os
import subprocess
from utilities import CONTEXT_SETTINGS


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('hmm', type=click.Path(exists=True))
@click.argument('pipolins-dir', type=click.Path(exists=True))
@click.argument('out-dir', type=click.Path())
def predict_atts_with_hmmer(hmm, pipolins_dir, out_dir):
    """
    TODO
    """
    pipolin_regions = os.listdir(pipolins_dir)
    os.makedirs(out_dir, exist_ok=True)
    os.chdir(out_dir)
    for pipolin in pipolin_regions:
        genome_id = pipolin.split(sep="-")[0]
        path = os.path.join(pipolins_dir, pipolin)
        subprocess.run(['nhmmscan', '--tblout', f'{genome_id}-atts.tbl', '-E', '1e-10', hmm, path])


if __name__ == '__main__':
    predict_atts_with_hmmer()
