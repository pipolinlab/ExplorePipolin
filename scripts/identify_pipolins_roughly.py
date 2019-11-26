#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os
import click
from utilities import CONTEXT_SETTINGS
from utilities import blast_genomes_against_seq
from utilities import check_dir
from utilities import read_blasttab
from utilities import Feature, Pipolin
from utilities import save_pipolins_to_shelve


def create_pipolins(genomes_dir, polbs_blast_path, atts_blast_path):
    pipolins = []

    for genome in os.listdir(genomes_dir):
        pipolins.append(Pipolin(strain_id=genome[:-3]))
    for i_p, pipolin in enumerate(pipolins):
        polbs = read_blasttab(os.path.join(polbs_blast_path, f'{pipolin.strain_id}_fmt7.txt'))
        atts = read_blasttab(os.path.join(atts_blast_path, f'{pipolin.strain_id}_fmt7.txt'))

        for entry in polbs:
            for hit in entry:
                print(hit.hit_frame)
                polb = Feature(start=hit.hit_start, end=hit.hit_end, node=entry.id)
                pipolins[i_p].polymerases.append(polb)

        for entry in atts:
            for hit in entry:
                att = Feature(start=hit.hit_start, end=hit.hit_end, node=entry.id)
                pipolins[i_p].atts.append(att)

    return pipolins


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

    pipolins = create_pipolins(genomes_dir, polbs_blast_path, atts_blast_path)
    out_file = os.path.join(out_dir, 'shelve.db')
    save_pipolins_to_shelve(out_file, pipolins)


if __name__ == '__main__':
    identify_pipolins_roughly()
