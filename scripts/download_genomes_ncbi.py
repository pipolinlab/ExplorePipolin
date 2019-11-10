#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os
import click
from Bio import SearchIO
from scripts.utilities import CONTEXT_SETTINGS
from scripts.utilities import ncbi_acc_download


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('blast-tab', type=click.Path(exists=True))
def download_genomes_ncbi(blast_tab):
    """
    BLAST_TAB is a tabular file that was generated when searching E. coli genomes against "reference" pi-polB
    from E. coli 3-373-03_S1_C2 (NZ_JNMI01000006.1). E. coli genome ids will be extracted from the file and the
    genome sequences in the fasta format will be downloaded from NCBI (see ./genomes)
    """
    blast_result = SearchIO.read(blast_tab, 'blast-tab', comments=True)
    ids = [i.id for i in blast_result]   # ids are all unique here
    ids.append(blast_result.id)

    # download the genomes
    os.chdir(os.path.dirname(blast_tab))
    if not os.path.exists('genomes'):
        os.mkdir('genomes')
    os.chdir('./genomes')
    ncbi_acc_download(ids)


if __name__ == '__main__':
    download_genomes_ncbi()