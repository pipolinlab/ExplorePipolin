#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os
import click
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from utilities import CONTEXT_SETTINGS
from utilities import Feature, Pipolin
from utilities import read_from_shelve


def return_att_seq(att, ref_seq):
    if att.start < att.end:
        subseq = ref_seq[att.start - 20:att.end + 20]
    else:
        subseq = ref_seq[att.end - 20:att.start + 20].reverse_complement()
    return subseq


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('shelve-file', type=click.Path(exists=True))
@click.argument('genomes-dir', type=click.Path(exists=True))
@click.argument('out-file')
def prepare_atts_for_msa(shelve_file, genomes_dir, out_file):
    """
    The scripts extracts the information about atts from the SHELVE_FILE and generates
    a FASTA file with att sequences.
    """
    pipolins = read_from_shelve(shelve_file, 'pipolins')

    genomes = {}
    for file in os.listdir(genomes_dir):
        genomes[file[:-3]] = SeqIO.to_dict(SeqIO.parse(os.path.join(genomes_dir, file), 'fasta'))

    att_records = []
    for pipolin in pipolins:
        atts = pipolin.atts
        if pipolin.is_complete_genome():
            ref_seq = genomes[pipolin.strain_id][pipolin.strain_id].seq
            att_records.extend(SeqRecord(seq=return_att_seq(att, ref_seq),
                                         id=f'{pipolin.strain_id}_{i+1}',
                                         description='') for i, att in enumerate(atts))
        else:
            for i, att in enumerate(atts):
                ref_seq = genomes[pipolin.strain_id][att.node].seq
                att_records.append(SeqRecord(seq=return_att_seq(att, ref_seq),
                                             id=f'{pipolin.strain_id}_{i+1}',
                                             description=''))

    print(f'The total number of atts is {len(att_records)}')
    lengths = [len(record.seq) - 40 for record in att_records]
    print(f'> Maximum att length is {max(lengths)}')
    print(f'> Minimum att length is {min(lengths)}')

    with open(out_file, 'w') as ouf:
        SeqIO.write(att_records, ouf, 'fasta')


if __name__ == '__main__':
    prepare_atts_for_msa()
