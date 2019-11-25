#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import click
import os
import subprocess
from utilities import CONTEXT_SETTINGS
from utilities import check_dir


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('hmm', type=click.Path(exists=True))
@click.argument('in-dir', type=click.Path(exists=True))
@click.argument('out-dir')
def predict_atts_with_hmmer(hmm, in_dir, out_dir):
    """
    TODO
    """
    pipolin_regions = os.listdir(in_dir)
    check_dir(out_dir)
    os.chdir(out_dir)
    for pipolin in pipolin_regions:
        genome_id = pipolin.split(sep="-")[0]
        path = os.path.join(in_dir, pipolin)
        subprocess.run(['nhmmscan', '--tblout', f'{genome_id}-atts.tbl', '-E', '1e-10', hmm, path])


if __name__ == '__main__':
    predict_atts_with_hmmer()
