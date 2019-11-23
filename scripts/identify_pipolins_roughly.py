#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import click
import os
from utilities import CONTEXT_SETTINGS
from utilities import blast_genomes_against_seq
from utilities import check_dir
from utilities import read_blasttab


class Feature:
    def __init__(self, start, end, node):
        self.start = start
        self.end = end
        self.node = node


class Pipolin:
    def __init__(self, strain_id):
        self.strain_id = strain_id

    atts = []
    polbs = []

    def add_att(self, att):
        self.atts.append(att)

    def add_polb(self, polb):
        self.polbs.append(polb)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('ref-polb', type=click.Path(exists=True))
@click.argument('ref-att', type=click.Path(exists=True))
@click.argument('genomes-dir', type=click.Path(exists=True))
@click.argument('out-dir')
def identify_pipolins_roughly(ref_polb, ref_att, genomes_dir, out_dir):
    """
    GENOMES_DIR contains each genome in a separate FASTA file (strain_id.fa).
    If there are several contigs in the genome, each contig should have unique name.
    If OUT_DIR exists, it should be empty.
    """
    check_dir(out_dir)
    polbs_blast_path = os.path.join(out_dir, 'polb_blast')
    atts_blast_path = os.path.join(out_dir, 'att_blast')
    blast_genomes_against_seq(genomes_dir, ref_polb, polbs_blast_path)
    blast_genomes_against_seq(genomes_dir, ref_att, atts_blast_path)

    pipolins = {}
    for file in os.listdir(genomes_dir):
        pipolins[file[:-3]] = Pipolin(file[:-3])

    for pipolin in pipolins:
        polbs = read_blasttab(os.path.join(polbs_blast_path, f'{pipolin}_fmt7.txt'))
        atts = read_blasttab(os.path.join(atts_blast_path, f'{pipolin}_fmt7.txt'))

        for entry in polbs:
            for hit in entry:
                pipolins[pipolin].add_polb(Feature(hit.hit_start, hit.hit_end, entry.id))

        for entry in atts:
            for hit in entry:
                pipolins[pipolin].add_att(Feature(hit.hit_start, hit.hit_end, entry.id))
        print(pipolins[pipolin])


if __name__ == '__main__':
    identify_pipolins_roughly()
