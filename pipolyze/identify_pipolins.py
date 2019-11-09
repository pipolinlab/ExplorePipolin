#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os
import click
import subprocess
from collections import namedtuple
from Bio import SearchIO
from utilities import CONTEXT_SETTINGS   # TODO fix this!


def extract_subject_boundaries(blast_tab):
    blast_results = SearchIO.read(blast_tab, 'blast-tab', comments=True)
    boundaries = {}
    for res in blast_results:
        boundaries[res.id] = []
        for hit in res:
            boundaries[res.id].append(sorted([hit.hit_start, hit.hit_end]))
    return boundaries


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('blast-tab', type=click.Path(exists=True))
@click.argument('reference', type=click.Path(exists=True))
@click.argument('genomes-dir', type=click.Path(exists=True))
def identify_pipolins(blast_tab, reference, genomes_dir):
    """
    The script takes a BLAST_TAB tabular file to extract imprecise coordinates of pi-polB for
    each genome. Then it runs blast with each genome in GENOMES_DIR against ref_att.fa to
    roughly estimate pipolin region boundaries. It creates a pipolins.csv file as an output.
    """
    # create a namedtuple object to store the info
    fields = ('id', 'polB1_s', 'polB1_e', 'polB2_s', 'polB2_e',
              'att1_s', 'att1_e', 'att2_s', 'att2_e', 'att3_s', 'att3_e')
    PipolinEntries = namedtuple('PolB_Genomes', fields, defaults=(None, ) * len(fields))
    pipolin_entries = []

    for id, regions in extract_subject_boundaries(blast_tab).items():
        regions.sort()
        pipolin_entries.append(PipolinEntries(id=id, polB1_s=regions[0][0], polB1_e=regions[0][1]))
        if len(regions) == 2:
            pipolin_entries[-1] = pipolin_entries[-1]._replace(polB2_s=regions[1][0], polB2_e=regions[1][1])

    # Add polB boundaries for the reference pipolin
    pipolin_entries.append(PipolinEntries(id='NZ_JNMI01000006.1', polB1_s=80192, polB1_e=82792))

    # BLAST each genome against ref_att.fa
    os.chdir(os.path.dirname(reference))
    if not os.path.exists('att_blast'):
        os.mkdir('att_blast')
    for entry in pipolin_entries:
        with open(f'./att_blast/{entry.id}_hits.txt', 'w') as ouf:
            subprocess.run(['blastn', '-query', reference,
                            '-subject', f'{os.path.join(genomes_dir, entry.id)}.fa',
                            '-outfmt', '7'], stdout=ouf)
    # Parse blast output
    for i, entry in enumerate(pipolin_entries):
        atts = extract_subject_boundaries(f'./att_blast/{entry.id}_hits.txt')[entry.id]
        atts.sort()
        # check that pi-polB is within the atts
        if pipolin_entries[i].polB1_s > atts[0][0] or \
                (pipolin_entries[i].polB2_s and pipolin_entries[i].polB2_s > atts[0][0]):
            if pipolin_entries[i].polB1_e < atts[-1][1] or \
                    (pipolin_entries[i].polB2_e and pipolin_entries[i].polB2_e < atts[-1][1]):
                pipolin_entries[i] = pipolin_entries[i]._replace(att1_s=atts[0][0], att1_e=atts[0][1],
                                                                 att2_s=atts[1][0], att2_e=atts[1][1])
                if len(atts) == 3:
                    pipolin_entries[i] = pipolin_entries[i]._replace(att3_s=atts[2][0], att3_e=atts[2][1])

    # store the data as a csv file
    with open('pipolins.csv', 'w') as ouf:
        print(','.join(PipolinEntries._fields), file=ouf)
        for i, _ in enumerate(pipolin_entries):
            words = ['None' if i is None else str(i) for i in pipolin_entries[i]]
            print(','.join(words), file=ouf)


if __name__ == '__main__':
    identify_pipolins()