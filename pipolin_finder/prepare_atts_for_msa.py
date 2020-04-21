#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os
import click
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from pipolin_finder.utilities import CONTEXT_SETTINGS
from pipolin_finder.utilities import Feature, GQuery
from pipolin_finder.utilities import read_from_shelve


def read_fasta_records(genomes_dir):
    genomes = {}
    for file in os.listdir(genomes_dir):
        genomes[file[:-3]] = SeqIO.to_dict(SeqIO.parse(os.path.join(genomes_dir, file), 'fasta'))
    return genomes


def return_att_seq(att, ref_seq):
    if att.start < att.end:
        subseq = ref_seq[att.start - 20:att.end + 20]
    else:
        subseq = ref_seq[att.end - 20:att.start + 20].reverse_complement()
    return subseq


def extract_all_atts(genomes, pipolins):
    att_records = []
    for pipolin in pipolins:
        atts = pipolin.atts
        if pipolin.is_single_contig():
            ref_seq = genomes[pipolin.strain_id][pipolin.strain_id].seq
            att_records.extend(SeqRecord(seq=return_att_seq(att, ref_seq),
                                         id=f'{pipolin.strain_id}_{i + 1}',
                                         description='') for i, att in enumerate(atts))
        else:
            for i, att in enumerate(atts):
                ref_seq = genomes[pipolin.strain_id][att.node].seq
                att_records.append(SeqRecord(seq=return_att_seq(att, ref_seq),
                                             id=f'{pipolin.strain_id}_{i + 1}',
                                             description=''))
    return att_records


def extract_specific_atts(genomes, pipolins, type):
    att_records = []
    for pipolin in pipolins:
        if pipolin.is_single_contig():
            atts = sorted(pipolin.atts, key=lambda x: x.start)
            ref_seq = genomes[pipolin.strain_id][pipolin.strain_id].seq
            if type == 'left':
                att_records.append(SeqRecord(seq=return_att_seq(atts[0], ref_seq),
                                             id=f'{pipolin.strain_id}_attL',
                                             description=''))
            elif type == 'right':
                att_records.append(SeqRecord(seq=return_att_seq(atts[-1], ref_seq),
                                             id=f'{pipolin.strain_id}_attR',
                                             description=''))
    return att_records


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('shelve-file', type=click.Path(exists=True))
@click.argument('genomes-dir', type=click.Path(exists=True))
@click.option('--type', required=True)
@click.argument('out-file')
def prepare_atts_for_msa(shelve_file, genomes_dir, type, out_file):
    """
    The scripts extracts the information about atts from the SHELVE_FILE and generates
    a FASTA file with att sequences. The argument --type is choice out of [all, left, right].
    """
    pipolins = read_from_shelve(shelve_file, 'pipolins')
    genomes = read_fasta_records(genomes_dir)

    if type == 'all':
        att_records = extract_all_atts(genomes, pipolins)
    elif type == 'left' or type == 'right':
        att_records = extract_specific_atts(genomes, pipolins, type)

    print(f'The total number of atts is {len(att_records)}')
    lengths = [len(record.seq) - 40 for record in att_records]
    print(f'> Maximum att length is {max(lengths)}')
    print(f'> Minimum att length is {min(lengths)}')

    with open(out_file, 'w') as ouf:
        SeqIO.write(att_records, ouf, 'fasta')


if __name__ == '__main__':
    prepare_atts_for_msa()
