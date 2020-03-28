#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os
import click
from prefect import task
from utilities import CONTEXT_SETTINGS
from utilities import blast_genomes_against_seq
from utilities import Feature, Pipolin
from utilities import save_to_shelve
from utilities import read_blastxml
from utilities import run_aragorn
from utilities import find_repeats


def feature_from_blasthit(hit, id):
    start = hit.hit_start if hit.hit_frame == 1 else hit.hit_end
    end = hit.hit_end if hit.hit_frame == 1 else hit.hit_start
    return Feature(start=start, end=end, frame=hit.hit_frame, node=id)


def create_pipolins(genomes, polbs_blast_path, atts_blast_path):
    pipolins = []

    for genome in genomes:
        pipolins.append(Pipolin(strain_id=os.path.basename(genome)[:-3]))
    for i_p, pipolin in enumerate(pipolins):
        polbs = read_blastxml(os.path.join(polbs_blast_path, f'{pipolin.strain_id}-fmt5.txt'))
        for entry in polbs:
            for hit in entry:
                polb_feature = feature_from_blasthit(hit, entry.id)
                pipolins[i_p].polymerases.append(polb_feature)

        # atts = read_blastxml(os.path.join(atts_blast_path, f'{pipolin.strain_id}-fmt5.txt'))
        # for entry in atts:
        #     for hit in entry:
        #         att_feature = feature_from_blasthit(hit, entry.id)
        #         pipolins[i_p].atts.append(att_feature)

    return pipolins


@task
def run_blast_against_polb(genomes, out_dir, ref_polb):
    polbs_blast_path = os.path.join(out_dir, 'polb_blast')
    blast_genomes_against_seq(genomes, ref_polb, polbs_blast_path)
    return polbs_blast_path


@task
def detect_trnas(genomes, out_dir):
    aragorn_results = os.path.join(out_dir, 'aragorn_results')
    run_aragorn(genomes, aragorn_results)
    return aragorn_results


@task
def find_att_repeats(genomes, out_dir):
    att_repeats = os.path.join(out_dir, 'att_repeats')
    find_repeats(genomes, att_repeats)
    return att_repeats


@task   # TODO: replace with de-novo atts search finction!
def run_blast_against_att(genomes, out_dir, ref_att):
    atts_blast_path = os.path.join(out_dir, 'att_blast')
    blast_genomes_against_seq(genomes, ref_att, atts_blast_path)
    return atts_blast_path


@task   # TODO: replace with aragorn search function!
def run_blast_against_trna(genomes, out_dir, ref_trna):
    trna_blast_path = os.path.join(out_dir, 'trna_blast')
    blast_genomes_against_seq(genomes, ref_trna, trna_blast_path)
    return trna_blast_path


@task
def identify_pipolins_roughly(genomes, out_dir, polbs_blast, atts_blast):
    pipolins = create_pipolins(genomes, polbs_blast, atts_blast)
    out_file = os.path.join(out_dir, 'shelve.db')
    save_to_shelve(out_file, pipolins, 'pipolins')
    return 'pipolins'


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
