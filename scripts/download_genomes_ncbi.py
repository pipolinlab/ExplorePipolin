#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os
import click
import subprocess
from utilities import CONTEXT_SETTINGS


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('metadata-file', type=click.Path(exists=True))
@click.argument('out-dir', type=click.Path())
def download_genomes_ncbi(metadata_file, out_dir):
    """
     and the genome sequences in FASTA format will be downloaded from NCBI
    to the OUT_DIR.
    NOTE: requires "ncbi-entrez-direct" package!
    """
    os.makedirs(out_dir, exist_ok=True)


if __name__ == '__main__':
    download_genomes_ncbi()
