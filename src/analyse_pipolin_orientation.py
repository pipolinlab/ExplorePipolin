#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import click
import os
from prefect import task
from utilities import CONTEXT_SETTINGS
from utilities import read_blastxml
from utilities import save_to_shelve


@task
def analyse_pipolin_orientation(in_dir, polbs_blast, atts_blast, trna_blast):
    orientations = {}
    strains = [file.split(sep='-')[0] for file in os.listdir(atts_blast)]
    for strain in strains:
        orientations[strain] = {}
    for strain in strains:
        atts = read_blastxml(os.path.join(atts_blast, f'{strain}-fmt5.txt'))
        trnas = read_blastxml(os.path.join(trna_blast, f'{strain}-fmt5.txt'))
        polbs = read_blastxml(os.path.join(polbs_blast, f'{strain}-fmt5.txt'))

        for entry in trnas:
            for hit in entry:
                if abs(hit.hit_start - hit.hit_end) >= 85:
                    trna_frame = hit.hit_frame
                    orientations[strain][entry.id] = -trna_frame

        for entry in atts:
            att_frames = []
            for hit in entry:
                att_frames.append(hit.hit_frame)

            if len(set(att_frames)) != 1:
                raise AssertionError('ATTs are expected to be located in the same frame, as they are direct repeats!')
            if entry.id in orientations[strain]:
                if att_frames[0] != orientations[strain][entry.id]:
                    raise AssertionError('ATT and tRNA should be on the different strands!')
            else:
                orientations[strain][entry.id] = att_frames[0]  # 1 or -1

        for entry in polbs:
            polb_frames = []
            for hit in entry:
                polb_frames.append(hit.hit_frame)

            if len(set(polb_frames)) != 1:  # an ambiguous case
                continue
            if entry.id not in orientations[strain]:  # the case, when no orientation info from tRNA or att
                orientations[strain][entry.id] = polb_frames[0]  # 1 or -1
    save_to_shelve(os.path.join(in_dir, 'shelve.db'), orientations, 'orientations')
    return 'orientations'


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('in-dir', type=click.Path(exists=True))
def main(in_dir):
    """
    TODO
    """
    analyse_pipolin_orientation(in_dir)


if __name__ == '__main__':
    main()
