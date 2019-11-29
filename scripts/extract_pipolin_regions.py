#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os
import click
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from utilities import CONTEXT_SETTINGS
from utilities import Feature, Pipolin
from utilities import check_dir
from utilities import read_pipolins_from_shelve


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('shelve-file', type=click.Path(exists=True))
@click.argument('genomes-dir', type=click.Path(exists=True))
@click.argument('out-dir')
@click.option('--long', is_flag=True, help='If long, pipolin regions are identified by leftmost and rightmost '
                                           'atts. If it is not specified, the closest atts to pi-polB will determine '
                                           'the pipolin bounds.')
def extract_pipolin_regions(shelve_file, genomes_dir, out_dir, long):
    """
    This script retrieves the information about pipolins from SHELVE_FILE
    and creates FASTA files with pipolin regions for their further annotation.
    ~23 min and 1116 files for 93 genomes.
    """
    pipolins = read_pipolins_from_shelve(shelve_file)

    genomes = {}
    for file in os.listdir(genomes_dir):
        genomes[file[:-3]] = SeqIO.to_dict(SeqIO.parse(os.path.join(genomes_dir, file), 'fasta'))

    check_dir(out_dir)
    for pipolin in pipolins:
        if pipolin.is_complete_genome():
            bounds = pipolin.get_pipolin_bounds(long)
            sequence = genomes[pipolin.strain_id][pipolin.strain_id].seq[bounds[0]:bounds[1]]
            records = SeqRecord(seq=sequence, id=pipolin.strain_id, description=f'{len(sequence)}')
        else:
            contings_bounds = pipolin.get_contigs_with_bounds(long)
            records = []
            for node, bounds in contings_bounds.items():
                sequence = genomes[pipolin.strain_id][node].seq[bounds[0]:bounds[1]]
                # TODO: fix LOCUS name, too long for contigs!!! (see include_atts_into_annotation.py)
                new_node_name = f'{node.split(sep="_")[0]}_{node.split(sep="_")[1]}'   # TODO: use node instead
                records.append(SeqRecord(seq=sequence, id=new_node_name, description=f'{len(sequence)}'))

        with open(os.path.join(out_dir, f'{pipolin.strain_id}-pipolin.fa'), 'w') as ouf:
            SeqIO.write(records, ouf, 'fasta')


if __name__ == '__main__':
    extract_pipolin_regions()
