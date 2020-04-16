#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import click
from prefect import task
from utilities import CONTEXT_SETTINGS


@task
def analyse_pipolin_orientation(gquery):
    target_trnas = gquery.target_trnas

    targeted_contigs = [trna.contig.contig_id for trna in target_trnas]
    if len(target_trnas) != len(targeted_contigs):
        raise AssertionError("I am expecting a single tRNA to overlap with a single att (per contig)!")

    for target_trna in target_trnas:
        targeted_atts = gquery.get_features_of_contig(contig_id=target_trna.contig.contig_id, feature_type='atts')
        atts_frames = [att.frame for att in targeted_atts]
        if len(set(atts_frames)) != 1:
            raise AssertionError('ATTs are expected to be in the same frame, as they are direct repeats!')
        if set(atts_frames).pop() == target_trna.frame:
            raise AssertionError('ATT and tRNA are expected to be on the different strands!')
        target_trna.contig.contig_orientation = - target_trna.frame

    for contig in gquery.contigs:
        contig_atts = gquery.get_features_of_contig(contig_id=contig.contig_id, feature_type='atts')
        if len(contig_atts) != 0:
            atts_frames = [att.frame for att in contig_atts]
            if len(set(atts_frames)) != 1:
                raise AssertionError('ATTs are expected to be located in the same frame, as they are direct repeats!')
            contig.contig_orientation = contig_atts[0].frame

    for contig in gquery.contigs:
        contig_polbs = gquery.get_features_of_contig(contig_id=contig.contig_id, feature_type='polbs')
        if len(contig_polbs) != 0:
            contig_polb_frames = [polb.frame for polb in contig_polbs]
            if len(set(contig_polb_frames)) != 1:   # an ambiguous case
                continue
            are_atts = gquery.get_features_of_contig(contig_id=contig.contig_id, feature_type='atts')
            are_trnas = gquery.get_features_of_contig(contig_id=contig.contig_id, feature_type='trnas')
            if len(are_atts) == 0 and len(are_trnas) == 0:
                contig.contig_orientation = contig_polbs[0].frame


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('in-dir', type=click.Path(exists=True))
def main(in_dir):
    """
    TODO
    """
    analyse_pipolin_orientation(in_dir)


if __name__ == '__main__':
    main()
