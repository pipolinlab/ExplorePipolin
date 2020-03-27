#!/usr/bin/env -S PYTHONPATH=${PWD}/src python3
# -*- encoding: utf-8 -*-

import click
from prefect import Flow, Parameter, task
from utilities import CONTEXT_SETTINGS
from identify_pipolins_roughly import identify_pipolins_roughly
from analyse_pipolin_orientation import analyse_pipolin_orientation

REF_POLB = './data/pi-polB.fa'
REF_ATT = './data/attL.fa'
REF_TRNA = './data/tRNA.fa'


def get_flow():
    with Flow('MAIN') as flow:
        genomes = Parameter('genomes')
        out_dir = Parameter('out_dir')

        identify_pipolins_roughly(genomes, out_dir, REF_POLB, REF_ATT, REF_TRNA)
        analyse_pipolin_orientation(out_dir)

    return flow


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('genomes', nargs=-1, type=click.Path(exists=True))
@click.option('--out-dir', type=click.Path())
def explore_pipolins(genomes, out_dir):
    """
    TODO
    """
    state = get_flow().run(genomes=genomes, out_dir=out_dir)
    assert state.is_successful()


if __name__ == '__main__':
    explore_pipolins()
