import os
from collections import defaultdict
from typing import MutableMapping, MutableSequence, Sequence

from BCBio import GFF
from Bio import SeqIO, SearchIO

from explore_pipolin.common import Genome, Orientation, define_genome_id, Repeat, Contig

SeqIORecords = MutableMapping[str, SeqIO.SeqRecord]


def read_genome_from_file(genome_file: str) -> Genome:
    genome = Genome(genome_id=define_genome_id(genome_file))
    genome_dict = read_seqio_records(file=genome_file, file_format='fasta')
    for key, value in genome_dict.items():
        contig = Contig(contig_id=key, contig_length=len(value.seq))
        genome.contigs.append(contig)
    return genome


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
            gff_records[entry.genome_id] = entry

    return gff_records


def write_genbank_records(gb_records: SeqIORecords, out_dir, genome: Genome):
    records = [record for record in gb_records.values()]
    with open(os.path.join(out_dir, f'{genome.genome_id}.gbk'), 'w') as ouf:
        SeqIO.write(records, ouf, 'genbank')


def write_gff_records(gff_records: SeqIORecords, out_dir, genome: Genome):
    records = [record for record in gff_records.values()]
    with open(os.path.join(out_dir, f'{genome.genome_id}.gff'), 'w') as ouf:
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


def save_left_right_subsequences(genome, left_window, right_window, repeats_dir):
    genome_seq = SeqIO.read(handle=genome, format='fasta')

    left_seq = genome_seq[left_window[0]:left_window[1]]
    right_seq = genome_seq[right_window[0]:right_window[1]]
    SeqIO.write(sequences=left_seq, format='fasta',
                handle=os.path.join(repeats_dir, define_genome_id(genome=genome) + '.left'))
    SeqIO.write(sequences=right_seq, format='fasta',
                handle=os.path.join(repeats_dir, define_genome_id(genome=genome) + '.right'))


def write_repeats(gquery, repeats: Sequence[Repeat], repeats_dir):
    with open(os.path.join(repeats_dir, gquery.gquery_id + '.repeats'), 'w') as ouf:
        polbs_locations = sorted([(x.start, x.end) for x in gquery.pipolbs], key=lambda x: x[0])
        print('left_repeat', 'right_repeat', 'length', 'polbs',
              'd_to_the_left', 'd_to_the_right', 'sequence', sep='\t', file=ouf)
        for repeat in repeats:
            print(repeat.left, repeat.right, repeat.left[1] - repeat.left[0],
                  ','.join([str(i) for i in polbs_locations]),
                  polbs_locations[0][0] - repeat.left[1], repeat.right[0] - polbs_locations[-1][-1], repeat.seq,
                  sep='\t', file=ouf)


def write_atts_denovo(atts_denovo, gquery, repeats_dir):
    with open(os.path.join(repeats_dir, gquery.gquery_id + '.atts'), 'w') as ouf:
        print('attL_start', 'attL_end', 'attR_start', 'attR_end', sep='\t', file=ouf)
        for att_pair in atts_denovo:
            print(att_pair[0][0], att_pair[0][1], att_pair[1][0], att_pair[1][1], sep='\t', file=ouf)
