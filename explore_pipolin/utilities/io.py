from typing import MutableMapping

from BCBio import GFF
from Bio import SeqIO, SearchIO

SeqIORecords = MutableMapping[str, SeqIO.SeqRecord]


def create_seqio_records_dict(file, file_format) -> SeqIORecords:
    file_formats = ['genbank', 'fasta']
    if file_format not in file_formats:
        raise AssertionError(f'Unknown file format: {file_format}! Only genbank or fasta formats are acceptable.')

    records = SeqIO.to_dict(SeqIO.parse(file, format=file_format))
    if len(records) == 0:
        raise AssertionError(f'Empty file {file}!')
    return records


def read_gff_records(file) -> SeqIORecords:
    gff_records = {}
    with open(file) as inf:
        for entry in GFF.parse(inf):
            gff_records[entry.id] = entry

    return gff_records


def write_genbank_records(gb_records: SeqIORecords, output_file: str):
    records = [record for record in gb_records.values()]
    with open(output_file, 'w') as ouf:
        SeqIO.write(records, ouf, 'genbank')


def write_gff_records(gff_records: SeqIORecords, output_file: str):
    records = [record for record in gff_records.values()]
    with open(output_file, 'w') as ouf:
        GFF.write(records, ouf)
        print('##FASTA', file=ouf)
        SeqIO.write(records, ouf, format='fasta')


def read_blastxml(blast_xml):
    return SearchIO.read(blast_xml, 'blast-xml')
