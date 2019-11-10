#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os
import click
from collections import namedtuple
from scripts.utilities import CONTEXT_SETTINGS   # TODO fix this!
from scripts.utilities import get_hit_positions_by_id   # TODO fix this!
from scripts.utilities import blast_genomes_against_seq   # TODO fix this!


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('blast-tab', type=click.Path(exists=True))
@click.argument('ref-att', type=click.Path(exists=True))
@click.argument('genomes-dir', type=click.Path(exists=True))
def identify_pipolins(blast_tab, ref_att, genomes_dir):
    """
    The script takes a BLAST_TAB tabular file to extract imprecise coordinates of pi-polB for
    each genome. Then it runs blast with each genome in GENOMES_DIR against ref_att.fa to
    roughly estimate pipolin region boundaries. It creates a pipolins.csv file as an output.
    """
    fields = ('id', 'polB1_s', 'polB1_e', 'polB2_s', 'polB2_e',
              'att1_s', 'att1_e', 'att2_s', 'att2_e', 'att3_s', 'att3_e')
    # create a namedtuple object to store the info
    Pipolin = namedtuple('Pipolin', fields, defaults=(None, ) * len(fields))
    pipolins = []

    for id, positions in get_hit_positions_by_id(blast_tab).items():
        pipolins.append(Pipolin(id=id, polB1_s=positions[0][0], polB1_e=positions[0][1]))
        if len(positions) == 2:
            pipolins[-1] = pipolins[-1]._replace(polB2_s=positions[1][0], polB2_e=positions[1][1])

    pipolins.append(Pipolin(id='NZ_JNMI01000006.1', polB1_s=80192, polB1_e=82792))

    blast_genomes_against_seq(genomes_dir, ref_att)
    for i_p, pipolin in enumerate(pipolins):
        atts = get_hit_positions_by_id(f'./att_blast/{pipolin.id}_hits.txt')[pipolin.id]
        # check that pi-polB is within the atts
        if pipolins[i_p].polB1_s > atts[0][0] or \
                (pipolins[i_p].polB2_s and pipolins[i_p].polB2_s > atts[0][0]):
            if pipolins[i_p].polB1_e < atts[-1][1] or \
                    (pipolins[i_p].polB2_e and pipolins[i_p].polB2_e < atts[-1][1]):
                pipolins[i_p] = pipolins[i_p]._replace(att1_s=atts[0][0], att1_e=atts[0][1],
                                                       att2_s=atts[1][0], att2_e=atts[1][1])
                if len(atts) == 3:
                    pipolins[i_p] = pipolins[i_p]._replace(att3_s=atts[2][0], att3_e=atts[2][1])

    # store the data as a csv file
    os.chdir(os.path.dirname(ref_att))
    with open('pipolins.csv', 'w') as ouf:
        print(','.join(Pipolin._fields), file=ouf)
        for i_p, _ in enumerate(pipolins):
            words = ['None' if i is None else str(i) for i in pipolins[i]]
            print(','.join(words), file=ouf)


if __name__ == '__main__':
    identify_pipolins()