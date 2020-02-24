#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os
import click
from Bio import SearchIO
from utilities import CONTEXT_SETTINGS
from utilities import check_dir
from utilities import ncbi_acc_download


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('blast-tab', type=click.Path(exists=True))
@click.argument('out-dir')
def download_genomes_ncbi(blast_tab, out_dir):
    """
    BLAST_TAB is a tabular file (hittable or outfmt=7) that was generated when searching E. coli genomes
    against "reference" pi-polB from E. coli 3-373-03_S1_C2 (NZ_JNMI01000006.1). E. coli genome accessions
    will be extracted from the file and the genome sequences in FASTA format will be downloaded from NCBI
    to the OUT_DIR.
    """
    blast_result = SearchIO.read(blast_tab, 'blast-tab', comments=True)
    ids = [i.id for i in blast_result]   # ids are all unique here
    ids.append(blast_result.id)   # NZ_JNMI01000006.1

    # TODO: it downloads contigs instead of the whole genome sequences in some cases.
    # SOLUTION: for each id in ids it is possible to fetch additional information:
    #  efetch -format native -mode xml -id NZ_JNMI01000006.1 -db sequences > output.xml
    # and then parse xml for the presence of the assembly id. After this it is possible
    # to download the whole genome with the same command:
    #  efetch -format fasta -id <assembly_id> -db sequences > output.fa
    # The efetch command is from ncbi-entrez-direct package.

    check_dir(out_dir)
    os.chdir(out_dir)
    ncbi_acc_download(ids)


if __name__ == '__main__':
    download_genomes_ncbi()
