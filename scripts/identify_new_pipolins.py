#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os
import click
from collections import namedtuple
from utilities import CONTEXT_SETTINGS
from utilities import get_hit_positions_by_id
from utilities import blast_genomes_against_seq
from utilities import save_as_csv


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('ref-polb', type=click.Path(exists=True))
@click.argument('ref-att', type=click.Path(exists=True))
@click.argument('contigs-dir', type=click.Path(exists=True))
def identify_new_pipolins(ref_polb, ref_att, contigs_dir):
    """
    TODO
    """
    fields = ('id', 'polB1_s', 'polB1_e', 'polB1_node',
              'polB2_s', 'polB2_e', 'polB2_node',
              'att1_s', 'att1_e', 'att1_node',
              'att2_s', 'att2_e', 'att2_node',
              'att3_s', 'att3_e', 'att3_node')
    # create a namedtuple object to store the info
    Pipolin = namedtuple('Pipolin', fields, defaults=(None, ) * len(fields))
    pipolins = []

    os.chdir(os.path.dirname(ref_polb))
    blast_genomes_against_seq(contigs_dir, ref_polb, 'polb_blast')
    blast_genomes_against_seq(contigs_dir, ref_att, 'att_blast')

    for file in os.listdir(contigs_dir):
        pipolins.append(Pipolin(id=file[:-3]))

    for i_p, pipolin in enumerate(pipolins):
        polbs = get_hit_positions_by_id(f'./polb_blast/{pipolin.id}_hits.txt')
        atts = get_hit_positions_by_id(f'./att_blast/{pipolin.id}_hits.txt')

        common_nodes = [node for node in atts if node in polbs]
        polbs_nodes_ordered = common_nodes.copy()
        # without list it returns a generator. Don't know why :-///
        polbs_nodes_ordered.extend(list(node for node in polbs if node not in common_nodes))
        atts_nodes_ordered = common_nodes.copy()
        atts_nodes_ordered.extend(list(node for node in atts if node not in common_nodes))

        polbs_regions = []
        for node in polbs_nodes_ordered:
            polbs_regions.extend(polbs[node])
        atts_regions = []
        for node in atts_nodes_ordered:
            atts_regions.extend(atts[node])

        # TODO: nodes are left empty!!!
        pipolins[i_p] = pipolins[i_p]._replace(polB1_s=polbs_regions[0][0], polB1_e=polbs_regions[0][1],
                                               att1_s=atts_regions[0][0], att1_e=atts_regions[0][1])
        if len(polbs_regions) == 2:
            pipolins[i_p] = pipolins[i_p]._replace(polB2_s=polbs_regions[1][0], polB2_e=polbs_regions[1][1])
        if len(atts_regions) > 1:
            pipolins[i_p] = pipolins[i_p]._replace(att2_s=atts_regions[1][0], att2_e=atts_regions[1][1])
        if len(atts_regions) > 2:
            pipolins[i_p] = pipolins[i_p]._replace(att2_s=atts_regions[2][0], att2_e=atts_regions[2][1])

    # store the data as a csv file
    save_as_csv(Pipolin, pipolins)


if __name__ == '__main__':
    identify_new_pipolins()