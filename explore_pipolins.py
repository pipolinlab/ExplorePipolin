#!/usr/bin/env -S PYTHONPATH=${PWD}/src python3
# -*- encoding: utf-8 -*-

import os
import click
from prefect import Flow, Parameter
from utilities import CONTEXT_SETTINGS
from identify_pipolins_roughly import identify_pipolins_roughly
from analyse_pipolin_orientation import analyse_pipolin_orientation

REF_POLB = './data/pi-polB.fa'
REF_ATT = './data/attL.fa'
REF_TRNA = './data/tRNA.fa'


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('genomes', nargs=-1, type=click.Path(exists=True))
@click.option('--out-dir', type=click.Path())
def explore_pipolins(genomes, out_dir):

    with Flow('MAIN') as flow:
        identify_pipolins_roughly(genomes, out_dir, REF_POLB, REF_ATT, REF_TRNA)

        pol_blast_dir = os.path.join(out_dir, 'polb_blast')
        att_blast_dir = os.path.join(out_dir, 'att_blast')
        trna_blast_dir = os.path.join(out_dir, 'trna_blast')
        shelve = os.path.join(out_dir, 'shelve.db')
        analyse_pipolin_orientation(pol_blast_dir, att_blast_dir, trna_blast_dir, shelve)

    state = flow.run()
    assert state.is_successful()


if __name__ == '__main__':
    explore_pipolins()
