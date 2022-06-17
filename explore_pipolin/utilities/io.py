from typing import MutableMapping

from BCBio import GFF
from Bio import SeqIO, SearchIO

SeqIORecords = MutableMapping[str, SeqIO.SeqRecord]


def create_seqio_records_dict(file, file_format) -> SeqIORecords:
    file_formats = ['genbank', 'fasta']
    if file_format not in file_formats:
        raise AssertionError(f'Unknown file format: {file_format}! Only genbank or fasta formats are acceptable.')

    records_dict = SeqIO.to_dict(SeqIO.parse(file, format=file_format))
    if len(records_dict) == 0:
        raise AssertionError(f'Empty file {file}!')
    return records_dict


def read_gff_records(file) -> SeqIORecords:
    gff_records = {}
    with open(file) as inf:
        for entry in GFF.parse(inf):
            gff_records[entry.id] = entry

    return gff_records


def write_seqio_records(records_dict: SeqIORecords, output_file: str, file_format) -> None:
    file_formats = ['genbank', 'fasta']
    if file_format not in file_formats:
        raise AssertionError(f'Unknown file format: {file_format}! Only genbank or fasta formats are acceptable.')

    records = [record for record in records_dict.values()]
    with open(output_file, 'w') as ouf:
        SeqIO.write(records, ouf, file_format)


def write_gff_records(gff_records: SeqIORecords, output_file: str):
    records = [record for record in gff_records.values()]
    with open(output_file, 'w') as ouf:
        GFF.write(records, ouf)
        print('##FASTA', file=ouf)
        SeqIO.write(records, ouf, format='fasta')


def read_blastxml(blast_xml):
    return SearchIO.read(blast_xml, 'blast-xml')


PRODUCTS_TO_COLOUR = {'default': '255 250 240'}


def read_colours(colors_tsv: str) -> None:
    with open(colors_tsv) as inf:
        for line in inf:
            if line[0] == '#':
                continue
            values = line.strip().split(sep='\t')
            if len(values) == 0:   # skip empty lines
                continue
            elif len(values) != 3:
                raise AssertionError(f'{len(values)} columns in the line {line.strip()}.\n'
                                     f'3 columns are expected.')
            else:
                PRODUCTS_TO_COLOUR[values[0]] = values[2]
