#!/usr/bin/env -S PYTHONPATH=${PWD}/src python3
# -*- encoding: utf-8 -*-

import click
from prefect import Flow, Parameter
from utilities import CONTEXT_SETTINGS
from identify_pipolins_roughly import create_gqueries
from identify_pipolins_roughly import run_blast_against_ref
from identify_pipolins_roughly import add_features_from_blast
from identify_pipolins_roughly import detect_trnas
from identify_pipolins_roughly import add_features_from_aragorn
from identify_pipolins_roughly import find_atts_denovo
from analyse_pipolin_orientation import analyse_pipolin_orientation
from extract_pipolin_regions import extract_pipolin_regions
from annotate_pipolins import annotate_pipolins
from predict_atts_with_hmmer import predict_atts_with_hmmer
from store_new_att_bounds import store_new_att_bounds
from include_atts_into_annotation import include_atts_into_annotation
from scaffold_gapped_pipolins import scaffold_gapped_pipolins
from easyfig_add_colours import easyfig_add_colours

REF_POLB = './data/pi-polB.fa'
REF_ATT = './data/attL.fa'
PROTEINS = './data/HHpred_proteins.faa'
ATT_HMM = './data/att.hmm'


def get_flow():
    with Flow('MAIN') as flow:
        genomes = Parameter('genomes')
        out_dir = Parameter('out_dir')

        gqueries = create_gqueries(genomes=genomes)

        polbs_blast = run_blast_against_ref(genomes=genomes, root_dir=out_dir,
                                            reference=REF_POLB, dir_name='polb_blast')
        add_features_from_blast(gqueries=gqueries, blast_dir=polbs_blast, feature_type='polymerases')

        atts_blast = run_blast_against_ref(genomes=genomes, root_dir=out_dir,
                                           reference=REF_ATT, dir_name='att_blast')
        add_features_from_blast(gqueries=gqueries, blast_dir=atts_blast, feature_type='atts')

        aragorn_results = detect_trnas(genomes=genomes, root_dir=out_dir)
        add_features_from_aragorn(gqueries=gqueries, aragorn_dir=aragorn_results)

        find_atts_denovo(genomes=genomes, gqueries=gqueries, root_dir=out_dir,
                         upstream_tasks=[add_features_from_blast, add_features_from_aragorn])

        # orientations = analyse_pipolin_orientation(out_dir, polbs_blast_dir, atts_blast, trna_blast)
        # rough_pipolins = extract_pipolin_regions(genomes, out_dir, gqueries, orientations, long=False)
        # prokka = annotate_pipolins(rough_pipolins, PROTEINS, out_dir)
        # att_hmmer = predict_atts_with_hmmer(ATT_HMM, rough_pipolins, out_dir)
        # short_pipolins = store_new_att_bounds(out_dir, 'short-gqueries', att_hmmer)
        # prokka_atts = include_atts_into_annotation(out_dir, short_pipolins, prokka)
        # prokka_atts_scaffolded = scaffold_gapped_pipolins(prokka_atts, out_dir, long=False)
        # easyfig_add_colours(prokka_atts_scaffolded, abricate_dir=None)

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
