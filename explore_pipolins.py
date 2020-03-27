#!/usr/bin/env -S PYTHONPATH=${PWD}/src python3
# -*- encoding: utf-8 -*-

import click
from prefect import Flow
from utilities import CONTEXT_SETTINGS
from download_genomes_ncbi import download_genomes_ncbi


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('metadata-file', type=click.Path(exists=True))
@click.argument('out-dir', type=click.Path())
def explore_pipolins(metadata_file, out_dir):

    with Flow('MAIN') as flow:
        download_genomes_ncbi(metadata_file, out_dir)

    state = flow.run()
    assert state.is_successful()


if __name__ == '__main__':
    explore_pipolins()
