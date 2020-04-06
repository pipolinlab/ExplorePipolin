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
from include_atts_into_annotation import include_atts_into_annotation
from scaffold_gapped_pipolins import is_scaffolding_required
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
        t3 = add_features_from_aragorn.map(gquery=gquery, aragorn_dir=aragorn_results, upstream_tasks=[t1, t2])

        atts_denovo = find_atts_denovo.map(genome=genomes, gquery=gquery, root_dir=unmapped(out_dir),
                                           upstream_tasks=[t3])
        t4 = add_features_atts_denovo.map(gquery=gquery, atts_denovo_dir=atts_denovo)

        t5 = analyse_pipolin_orientation.map(gquery=gquery, upstream_tasks=[t4])
        t6 = is_scaffolding_required.map(gquery=gquery, upstream_tasks=[t5])

        pipolin_sequences = extract_pipolin_regions.map(genome=genomes, gquery=gquery,
                                                        root_dir=unmapped(out_dir), upstream_tasks=[t6])
        prokka = annotate_pipolins.map(gquery=gquery, pipolins_dir=pipolin_sequences,
                                       proteins=unmapped(PROTEINS), root_dir=unmapped(out_dir))
        prokka_atts = include_atts_into_annotation.map(gquery=gquery, prokka_dir=prokka,
                                                       root_dir=unmapped(out_dir))
        easyfig_add_colours.map(gquery=gquery, in_dir=prokka_atts, abricate_dir=unmapped(None))

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
