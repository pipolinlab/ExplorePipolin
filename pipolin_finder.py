#!/usr/bin/env -S PYTHONPATH=${PWD}/src python3
# -*- encoding: utf-8 -*-

import click
from prefect import Flow, Parameter, unmapped
# from prefect.tasks.core.constants import Constant
from utilities import CONTEXT_SETTINGS
from identify_pipolins_roughly import create_gquery
from identify_pipolins_roughly import run_blast_against_ref
from identify_pipolins_roughly import add_features_from_blast
from identify_pipolins_roughly import detect_trnas_with_aragorn
from identify_pipolins_roughly import add_features_from_aragorn
from identify_pipolins_roughly import find_atts_denovo
from identify_pipolins_roughly import add_features_atts_denovo
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

        gquery = create_gquery.map(genome=genomes)

        polbs_blast = run_blast_against_ref.map(genome=genomes, root_dir=unmapped(out_dir),
                                                reference=unmapped(REF_POLB), dir_name=unmapped('polb_blast'))
        t1 = add_features_from_blast.map(gquery=gquery, blast_dir=polbs_blast, feature_type=unmapped('polymerases'))

        atts_blast = run_blast_against_ref.map(genome=genomes, root_dir=unmapped(out_dir),
                                               reference=unmapped(REF_ATT), dir_name=unmapped('att_blast'))
        t2 = add_features_from_blast.map(gquery=gquery, blast_dir=atts_blast, feature_type=unmapped('atts'))

        aragorn_results = detect_trnas_with_aragorn.map(genome=genomes, root_dir=unmapped(out_dir))
        t3 = add_features_from_aragorn.map(gquery=gquery, aragorn_dir=aragorn_results)

        atts_denovo = find_atts_denovo.map(genome=genomes, gquery=gquery, root_dir=unmapped(out_dir),
                                           upstream_tasks=[t1, t2, t3])
        t4 = add_features_atts_denovo.map(gquery=gquery, atts_denovo_dir=atts_denovo)
        # analyse_pipolin_orientation(gquery=gquery, upstream_tasks=[t4])
        # rough_pipolins = extract_pipolin_regions(genome, out_dir, gquery, orientations, long=False)
        # prokka = annotate_pipolins(rough_pipolins, PROTEINS, out_dir)
        # att_hmmer = predict_atts_with_hmmer(ATT_HMM, rough_pipolins, out_dir)
        # short_pipolins = store_new_att_bounds(out_dir, 'short-gquery', att_hmmer)
        # prokka_atts = include_atts_into_annotation(out_dir, short_pipolins, prokka)
        # prokka_atts_scaffolded = scaffold_gapped_pipolins(prokka_atts, out_dir, long=False)
        # easyfig_add_colours(prokka_atts_scaffolded, abricate_dir=None)

    return flow


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('genomes', type=click.Path(exists=True), nargs=-1)
@click.option('--out-dir', type=click.Path())
def explore_pipolins(genomes, out_dir):
    """
    TODO
    """
    state = get_flow().run(genomes=genomes, out_dir=out_dir)
    assert state.is_successful()


if __name__ == '__main__':
    explore_pipolins()
