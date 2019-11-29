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
            length_by_contig = {}
            for contig in genomes[pipolin.strain_id].values():
                length_by_contig[contig.id] = len(contig.seq)
            contigs_bounds = pipolin.get_contigs_with_bounds(length_by_contig, long)
            records = []
            for node, bounds in contigs_bounds.items():
                sequence = genomes[pipolin.strain_id][node].seq[bounds[0]:bounds[1]]
                records.append(SeqRecord(seq=sequence, id=node, description=f'{len(sequence)}'))

        with open(os.path.join(out_dir, f'{pipolin.strain_id}-pipolin.fa'), 'w') as ouf:
            SeqIO.write(records, ouf, 'fasta')


if __name__ == '__main__':
    extract_pipolin_regions()
