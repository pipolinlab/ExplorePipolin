#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import click
import os
from utilities import CONTEXT_SETTINGS
from utilities import read_from_shelve
from utilities import read_blastxml
from utilities import save_to_shelve


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('shelve-file', type=click.Path(exists=True))
@click.argument('att-blast-dir', type=click.Path(exists=True))
def analyse_pipolin_orientation(shelve_file, att_blast_dir):
    """
    TODO
    """
    pipolins = read_from_shelve(shelve_file, 'pipolins')

    orientations = {}

    for i_p, pipolin in enumerate(pipolins):
        atts = read_blastxml(os.path.join(att_blast_dir, f'{pipolin.strain_id}_fmt5.txt'))
        for entry in atts:
            atts_info = []
            for hit in entry:
                atts_info.append([hit.hit_frame, hit.evalue])

            if len(atts_info) > 1:
                orientations[pipolin.strain_id] = {entry.id: atts_info[0][0]}   # 1 or -1
            else:
                orientations[pipolin.strain_id] = {entry.id: 1}   # we don't know, leave as is

    save_to_shelve(shelve_file, orientations, 'orientations')


if __name__ == '__main__':
    analyse_pipolin_orientation()
