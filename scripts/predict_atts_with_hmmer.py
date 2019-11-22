#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import click
import os
import subprocess
from utilities import CONTEXT_SETTINGS


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('hmm', type=click.Path(exists=True))
@click.argument('in-dir', type=click.Path(exists=True))
@click.argument('out-dir', type=click.Path(exists=True))
def predict_atts_with_hmmer(hmm, in_dir, out_dir):
    """
    TODO
    """
    genomes = os.listdir(in_dir)
    os.chdir(out_dir)
    for genome in genomes:
        path = os.path.join(in_dir, genome)
        subprocess.run(['nhmmscan', '--tblout', f'{genome[:-3]}_atts.tbl', '-E', '1e-10', hmm, path])


if __name__ == '__main__':
    predict_atts_with_hmmer()
