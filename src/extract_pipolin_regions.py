#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os
import click
from prefect import task
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from utilities import CONTEXT_SETTINGS
from utilities import Feature, GQuery
from utilities import read_from_shelve
from utilities import read_seqio_records


@task
def extract_pipolin_regions(genome, gquery: GQuery, root_dir):
    # TODO
    if gquery.pipolin_fragment is None:
        pass
    genome_dict = read_seqio_records(file=genome, file_format='fasta')

    pipolins_dir = os.path.join(root_dir, 'pipolin_sequences')
    os.makedirs(pipolins_dir, exist_ok=True)

    with open(os.path.join(pipolins_dir, gquery.gquery_id + '.fa'), 'w') as ouf:
        seq = genome_dict[gquery.pipolin_fragment.contig.contig_id]
        SeqIO.write(seq[gquery.pipolin_fragment.start:gquery.pipolin_fragment.end], ouf, 'fasta')

    # for pipolin in pipolins:
    #     if pipolin.is_single_contig():
    #         bounds = pipolin.get_pipolin_bounds(long)
    #         sequence = genomes[pipolin.strain_id][pipolin.strain_id].seq[bounds[0]:bounds[1]]
    #         orientation = orientations[pipolin.strain_id][pipolin.strain_id]
    #         sequence = sequence.reverse_complement() if orientation == -1 else sequence
    #         records = SeqRecord(seq=sequence, id=pipolin.strain_id, description=f'{len(sequence)}')
    #     else:
    #         length_by_contig = {}
    #         for contig in genomes[pipolin.strain_id].values():
    #             length_by_contig[contig.id] = len(contig.seq)
    #         contigs_bounds = pipolin.get_contigs_with_bounds(length_by_contig, long)
    #         records = []
    #         for node, bounds in contigs_bounds.items():
    #             sequence = genomes[pipolin.strain_id][node].seq[bounds[0]:bounds[1]]
    #             orientation = orientations[pipolin.strain_id][node] if node in orientations[pipolin.strain_id] else 1
    #             sequence = sequence.reverse_complement() if orientation == -1 else sequence
    #             records.append(SeqRecord(seq=sequence, id=node, description=f'{len(sequence)}'))
    #
    #     with open(os.path.join(shelve_in_dir, 'rough_pipolins', f'{pipolin.strain_id}-pipolin.fa'), 'w') as ouf:
    #         SeqIO.write(records, ouf, 'fasta')

    return pipolins_dir


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('genomes', type=click.Path(exists=True))
@click.argument('shelve-in-dir', type=click.Path(exists=True))
@click.argument('out-dir', type=click.Path())
@click.option('--long', is_flag=True, help='If long, pipolin regions are identified by leftmost and rightmost '
                                           'atts. If it is not specified, the closest atts to pi-polB will determine '
                                           'the pipolin bounds.')
def main(genomes, shelve_in_dir, out_dir, long):
    """
    This script retrieves the information about pipolins from SHELVE_FILE
    and creates FASTA files with pipolin regions for their further annotation.
    """
    extract_pipolin_regions(genomes, shelve_in_dir, out_dir, long)


if __name__ == '__main__':
    main()
