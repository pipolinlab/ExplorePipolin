#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import click
from prefect import task
from utilities import CONTEXT_SETTINGS
from identify_pipolins_roughly import is_overlapping


@task
def analyse_pipolin_orientation(gquery):
    # TODO: do we expect a single pipolin per genome?
    # TODO: how to be with atts_denovo if they are?
    anchor_trnas = []
    for trna in gquery.trnas:
        for att in gquery.atts:
            if att.contig.contig_id == trna.contig.contig_id:
                if is_overlapping(range1=(att.start, att.end), range2=(trna.start, trna.end)):
                    anchor_trnas.append(trna)

    anchor_contigs = [trna.contig.contig_id for trna in anchor_trnas]
    if len(anchor_trnas) != len(anchor_contigs):
        raise AssertionError("I am expecting a single tRNA to overlap with a single att (per contig)!")

    for anchor_trna in anchor_trnas:
        anchor_atts = gquery.get_features_of_contig(contig_id=anchor_trna.contig.contig_id, feature_type='atts')
        anchor_att_frames = [att.frame for att in anchor_atts]
        if len(set(anchor_att_frames)) != 1:
            raise AssertionError('ATTs are expected to be located in the same frame, as they are direct repeats!')
        if set(anchor_att_frames).pop() == anchor_trna.frame:
            raise AssertionError('ATT and tRNA should be on the different strands!')
        anchor_trna.contig.contig_orientation = - anchor_trna.frame

    for contig in gquery.contigs:
        contig_atts = gquery.get_features_of_contig(contig_id=contig.contig_id, feature_type='atts')
        if len(contig_atts) != 0:
            contig_att_frames = [att.frame for att in contig_atts]
            if len(set(contig_att_frames)) != 1:
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
