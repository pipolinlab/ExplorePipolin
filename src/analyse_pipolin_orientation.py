#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import click
import os
from utilities import CONTEXT_SETTINGS
from utilities import read_blastxml
from utilities import save_to_shelve


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('att-blast-dir', type=click.Path(exists=True))
@click.argument('trna-blast-dir', type=click.Path(exists=True))
@click.argument('polb-blast-dir', type=click.Path(exists=True))
@click.argument('shelve-file', type=click.Path(exists=True))
def analyse_pipolin_orientation(att_blast_dir, trna_blast_dir, polb_blast_dir, shelve_file):
    """
    TODO
    """
    orientations = {}

    strains = [file.split(sep='-')[0] for file in os.listdir(att_blast_dir)]
    for strain in strains:
        orientations[strain] = {}

    for strain in strains:
        atts = read_blastxml(os.path.join(att_blast_dir, f'{strain}-fmt5.txt'))
        trnas = read_blastxml(os.path.join(trna_blast_dir, f'{strain}-fmt5.txt'))
        polbs = read_blastxml(os.path.join(polb_blast_dir, f'{strain}-fmt5.txt'))

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
                orientations[strain][entry.id] = att_frames[0]   # 1 or -1

        for entry in polbs:
            polb_frames = []
            for hit in entry:
                polb_frames.append(hit.hit_frame)

            if len(set(polb_frames)) != 1:   # an ambiguous case
                continue
            if entry.id not in orientations[strain]:   # the case, when no orientation info from tRNA or att
                orientations[strain][entry.id] = polb_frames[0]   # 1 or -1

    save_to_shelve(shelve_file, orientations, 'orientations')


if __name__ == '__main__':
    analyse_pipolin_orientation()
