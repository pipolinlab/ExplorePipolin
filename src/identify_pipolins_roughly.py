#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os
import click
from prefect import task
from utilities import CONTEXT_SETTINGS
from utilities import blast_seqs_against_seq
from utilities import Feature, Pipolin
from utilities import save_to_shelve
from utilities import read_blastxml


def feature_from_blasthit(hit, id):
    start = hit.hit_start if hit.hit_frame == 1 else hit.hit_end
    end = hit.hit_end if hit.hit_frame == 1 else hit.hit_start
    return Feature(start=start, end=end, node=id)


def create_pipolins(genomes, polbs_blast_path, atts_blast_path):
    pipolins = []

    for genome in genomes:
        pipolins.append(Pipolin(strain_id=os.path.basename(genome)[:-3]))
    for i_p, pipolin in enumerate(pipolins):
        polbs = read_blastxml(os.path.join(polbs_blast_path, f'{pipolin.strain_id}-fmt5.txt'))
        atts = read_blastxml(os.path.join(atts_blast_path, f'{pipolin.strain_id}-fmt5.txt'))

        for entry in polbs:
            for hit in entry:
                polb_feature = feature_from_blasthit(hit, entry.id)
                pipolins[i_p].polymerases.append(polb_feature)

        for entry in atts:
            for hit in entry:
                att_feature = feature_from_blasthit(hit, entry.id)
                pipolins[i_p].atts.append(att_feature)

    return pipolins


@task
def identify_pipolins_roughly(genomes, out_dir, ref_polb, ref_att, ref_trna):
    os.makedirs(out_dir, exist_ok=True)
    polbs_blast_path = os.path.join(out_dir, 'polb_blast')
    atts_blast_path = os.path.join(out_dir, 'att_blast')
    trna_blast_path = os.path.join(out_dir, 'trna_blast')
    print('>>> Running BLAST against piPolB...')
    blast_seqs_against_seq(genomes, ref_polb, polbs_blast_path)
    print('>>> Running BLAST against AttL...')
    blast_seqs_against_seq(genomes, ref_att, atts_blast_path)
    print('>>> Running BLAST against tRNA...')
    blast_seqs_against_seq(genomes, ref_trna, trna_blast_path)
    print('>>> Creating "pipolins" shelve object...')
    pipolins = create_pipolins(genomes, polbs_blast_path, atts_blast_path)
    out_file = os.path.join(out_dir, 'shelve.db')
    save_to_shelve(out_file, pipolins, 'pipolins')
    print('DONE!')


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('ref-polb', type=click.Path(exists=True))
@click.argument('ref-att', type=click.Path(exists=True))
@click.argument('ref-trna', type=click.Path(exists=True))
@click.argument('genomes', nargs=-1, type=click.Path(exists=True))
@click.argument('out-dir')
def main(ref_polb, ref_att, ref_trna, genomes, out_dir):
    """
    GENOMES_DIR contains each genome in a separate FASTA file (strain_id.fa).
    If there are several contigs in the genome, each contig should have unique name.
    If OUT_DIR exists, it should be empty.
    """
    identify_pipolins_roughly(genomes, out_dir, ref_polb, ref_att, ref_trna)


if __name__ == '__main__':
    main()
