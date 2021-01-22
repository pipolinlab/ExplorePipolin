from collections import defaultdict
from typing import MutableMapping, MutableSequence, Sequence, Tuple

from BCBio import GFF
from Bio import SeqIO, SearchIO

from explore_pipolin.common import Strand, Contig, ContigID

SeqIORecords = MutableMapping[str, SeqIO.SeqRecord]


def read_genome_contigs_from_file(genome_file: str) -> MutableSequence[Contig]:
    genome_dict = create_seqio_records_dict(file=genome_file, file_format='fasta')
    contigs = []
    for key, value in genome_dict.items():
        contigs.append(Contig(contig_id=ContigID(key), contig_length=len(value.seq)))
    return contigs


def create_seqio_records_dict(file, file_format) -> SeqIORecords:
    file_formats = ['genbank', 'fasta']
    if file_format not in file_formats:
        raise AssertionError(f'Unknown file format: {file_format}! Only genbank or fasta formats are acceptable.')

    records = SeqIO.to_dict(SeqIO.parse(file, format=file_format))
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
                        entries[entry].append((start, end, Strand.REVERSE))
                    else:
                        start, end = (int(i) for i in coordinates[1:-1].split(sep=','))
                        entries[entry].append((start, end, Strand.FORWARD))

    return entries


def read_blastxml(blast_xml):
    return SearchIO.read(blast_xml, 'blast-xml')


def create_pipolb_entries(hmmsearch_table: str) -> Sequence[Tuple[str, int, int, int]]:
    entries = []
    with open(hmmsearch_table) as inf:
        content = list(SearchIO.parse(inf, 'hmmer3-tab'))
        if len(content) != 1:
            raise AssertionError(f'More than a single query in {hmmsearch_table}! Should be only one.')
        for hit in content[0]:
            name = '_'.join(hit.id.split(sep='_')[:-1])
            description = hit.description.split(sep=' ')
            entries.append((name, int(description[1]), int(description[3]), int(description[5])))

    return entries
