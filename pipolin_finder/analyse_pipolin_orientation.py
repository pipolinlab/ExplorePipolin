#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import click
from prefect import task
from pipolin_finder.utilities import CONTEXT_SETTINGS


@task
def analyse_pipolin_orientation(gquery):
    gquery.is_single_target_trna_per_contig()
    for contig in gquery.contigs:
        gquery.set_contig_orientation(contig)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('in-dir', type=click.Path(exists=True))
def main(in_dir):
    """
    TODO
    """
    analyse_pipolin_orientation(in_dir)


if __name__ == '__main__':
    main()
