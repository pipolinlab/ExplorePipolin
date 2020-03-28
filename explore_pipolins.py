#!/usr/bin/env -S PYTHONPATH=${PWD}/src python3
# -*- encoding: utf-8 -*-

import os
import click
from prefect import task, Flow, Parameter
from utilities import CONTEXT_SETTINGS
from identify_pipolins_roughly import run_blast_against_att, run_blast_against_polb, run_blast_against_trna
from identify_pipolins_roughly import identify_pipolins_roughly
from analyse_pipolin_orientation import analyse_pipolin_orientation
from extract_pipolin_regions import extract_pipolin_regions
from annotate_pipolins import annotate_pipolins
from predict_atts_with_hmmer import predict_atts_with_hmmer
from store_new_att_bounds import store_new_att_bounds
from include_atts_into_annotation import include_atts_into_annotation

REF_POLB = './data/pi-polB.fa'
REF_ATT = './data/attL.fa'
REF_TRNA = './data/tRNA.fa'
PROTEINS = './data/HHpred_proteins.faa'
ATT_HMM = './data/att.hmm'


def get_flow():
    with Flow('MAIN') as flow:
        genomes = Parameter('genomes')
        out_dir = Parameter('out_dir')

        polbs_blast = run_blast_against_polb(genomes, out_dir, REF_POLB)
        atts_blast = run_blast_against_att(genomes, out_dir, REF_ATT)
        trna_blast = run_blast_against_trna(genomes, out_dir, REF_TRNA)
        pipolins = identify_pipolins_roughly(genomes, out_dir, polbs_blast, atts_blast)
        orientations = analyse_pipolin_orientation(out_dir, polbs_blast, atts_blast, trna_blast)
        rough_pipolins = extract_pipolin_regions(genomes, out_dir, pipolins, orientations, long=False)
        prokka = annotate_pipolins(rough_pipolins, PROTEINS, out_dir)
        att_hmmer = predict_atts_with_hmmer(ATT_HMM, rough_pipolins, out_dir)
        short_pipolins = store_new_att_bounds(out_dir, 'short-pipolins', att_hmmer)
        prokka_atts = include_atts_into_annotation(out_dir, short_pipolins, prokka)

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
