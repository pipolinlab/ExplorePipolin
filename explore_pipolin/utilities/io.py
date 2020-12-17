import os
from collections import defaultdict
from typing import MutableMapping, MutableSequence, Sequence

from BCBio import GFF
from Bio import SeqIO, SearchIO

from explore_pipolin.common import Genome, Orientation, RepeatPair, Contig, Feature, FeatureType

SeqIORecords = MutableMapping[str, SeqIO.SeqRecord]


def read_genome_contigs_from_file(genome_file: str) -> MutableSequence[Contig]:
    genome_dict = read_seqio_records(file=genome_file, file_format='fasta')
    contigs = []
    for key, value in genome_dict.items():
        contigs.append(Contig(contig_id=key, contig_length=len(value.seq)))
    return contigs


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


def save_left_right_subsequences(left_window: Feature, right_window: Feature, repeats_dir: str):
    genome_seq = SeqIO.read(handle=left_window.genome.genome_file, format='fasta')

    left_seq = genome_seq[left_window.start:left_window.end]
    right_seq = genome_seq[right_window.start:right_window.end]
    SeqIO.write(sequences=left_seq, format='fasta',
                handle=os.path.join(repeats_dir, left_window.genome.genome_id + '.left'))
    SeqIO.write(sequences=right_seq, format='fasta',
                handle=os.path.join(repeats_dir, left_window.genome.genome_id + '.right'))


def write_repeats(genome: Genome, repeats: Sequence[RepeatPair], out_dir: str):
    with open(os.path.join(out_dir, genome.genome_id + '.repeats'), 'w') as ouf:
        polbs_locations = sorted([(x.start, x.end) for x in genome.features.get_features(
            FeatureType.PIPOLB)], key=lambda x: x[0])
        print('left_rep_range', 'right_rep_range', 'length', 'polbs',
              'd_to_the_left', 'd_to_the_right', 'left_rep_seq', 'right_rep_seq', sep='\t', file=ouf)
        for repeat in repeats:
            print((repeat.left.start, repeat.left.end), (repeat.right.start, repeat.right.end),
                  repeat.left.end - repeat.left.start, ','.join([str(i) for i in polbs_locations]),
                  polbs_locations[0][0] - repeat.left.end, repeat.right.start - polbs_locations[-1][-1],
                  repeat.left_seq, repeat.right_seq, sep='\t', file=ouf)


def write_atts_denovo(atts_denovo: Sequence[Feature], genome: Genome, atts_denovo_dir: str):
    with open(os.path.join(atts_denovo_dir, genome.genome_id + '.atts_denovo'), 'w') as ouf:
        print('att_start', 'att_end', sep='\t', file=ouf)
        for att in atts_denovo:
            print(att.start, att.end, sep='\t', file=ouf)
