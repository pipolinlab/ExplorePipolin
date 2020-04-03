#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import click
import os
from prefect import task
from utilities import CONTEXT_SETTINGS
from utilities import save_to_shelve
from identify_pipolins_roughly import is_overlapping


@task
def analyse_pipolin_orientation(gquery):
    # TODO: do we expect a single pipolin per genome?
    anchor_trna = []
    for trna in gquery.trnas:
        for att in gquery.atts:
            if att.contig.contig_id == trna.contig.contig_id:
                if is_overlapping(range1=(att.start, att.end), range2=(trna.start, trna.end)):
                    anchor_trna.append(trna)
    if len(anchor_trna) > 1:
        raise AssertionError("I am expecting a single tRNA to overlap with a single att (per contig)!")

    # TODO: how to be with atts_denovo if they are?

    if len(anchor_trna) == 1:
        anchor_atts = gquery.get_features_of_contig(contig_id=anchor_trna[0].contig.contig_id, feature_type='atts')
        anchor_att_frames = [att.frame for att in anchor_atts]
        if len(set(anchor_att_frames)) != 1:
            raise AssertionError('ATTs are expected to be located in the same frame, as they are direct repeats!')
        if set(anchor_att_frames).pop() == anchor_trna[0].frame:
            raise AssertionError('ATT and tRNA should be on the different strands!')
        anchor_trna[0].contig.contig_orientation = - anchor_trna[0].frame
    else:
        raise AssertionError(f'The anchor tRNA was not found for {gquery.gquery_id}!!! ???')

    for contig in gquery.contigs:
        contig_atts = gquery.get_features_of_contig(contig_id=contig.contig_id, feature_type='atts')
        if len(contig_atts) != 0:
            contig_att_frames = [att.frame for att in contig_atts]
            if len(set(contig_att_frames)) != 1:
                raise AssertionError('ATTs are expected to be located in the same frame, as they are direct repeats!')
            contig.contig_orientation = contig_atts[0].frame

    for contig in gquery.contigs:
        contig_polbs = gquery.get_features_of_contig(contig_id=contig.contig_id, feature_type='polymerases')
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
