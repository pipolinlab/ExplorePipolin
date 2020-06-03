import os
from collections import defaultdict
from typing import MutableMapping, MutableSequence

from BCBio import GFF
from Bio import SeqIO, SearchIO

from explore_pipolin.utilities import Orientation

SeqIORecords = MutableMapping[str, SeqIO.SeqRecord]


def read_seqio_records(file, file_format) -> SeqIORecords:
    file_formats = ['genbank', 'fasta']
    if file_format not in file_formats:
        raise AssertionError('Only genbank or fasta formats are acceptable!')

    records = SeqIO.to_dict(SeqIO.parse(file, format=file_format))
    return records


def read_gff_records(file) -> SeqIORecords:
    gff_records = {}
    with open(file) as inf:
        for entry in GFF.parse(inf):
            gff_records[entry.id] = entry

    return gff_records


def write_genbank_records(gb_records: SeqIORecords, out_dir, gquery):
    records = [record for record in gb_records.values()]
    with open(os.path.join(out_dir, f'{gquery.gquery_id}.gbk'), 'w') as ouf:
        SeqIO.write(records, ouf, 'genbank')


def write_gff_records(gff_records: SeqIORecords, out_dir, gquery):
    records = [record for record in gff_records.values()]
    with open(os.path.join(out_dir, f'{gquery.gquery_id}.gff'), 'w') as ouf:
        GFF.write(records, ouf)
        print('##FASTA', file=ouf)
        SeqIO.write(records, ouf, format='fasta')


def read_aragorn_batch(aragorn_batch) -> MutableMapping[str, MutableSequence]:
    entries = defaultdict(list)
    with open(aragorn_batch) as inf:
        for line in inf:
            if line[0] == '>':
                entry = line.strip().split(sep=' ')[0][1:]
            else:
                hit = line.split(sep='\t')
                if len(hit) > 1:
                    coordinates = hit[0].split(sep=' ')[-1]
                    if coordinates[0] == 'c':
                        start, end = (int(i) for i in coordinates[2:-1].split(sep=','))
                        entries[entry].append((start, end, Orientation.REVERSE))
                    else:
                        start, end = (int(i) for i in coordinates[1:-1].split(sep=','))
                        entries[entry].append((start, end, Orientation.FORWARD))

    return entries


def read_blastxml(blast_xml):
    return SearchIO.read(blast_xml, 'blast-xml')
