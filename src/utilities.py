import csv
import os
import shelve
from collections import defaultdict
from enum import Enum, auto
from typing import Sequence, MutableSequence, Mapping, MutableMapping
from itertools import groupby
import subprocess

from BCBio import GFF
from Bio import SearchIO, SeqIO

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


class Orientation(Enum):
    FORWARD = auto()
    REVERSE = auto()

    @staticmethod
    def orientation_from_blast(hit_frame):
        if hit_frame == 1:
            return Orientation.FORWARD
        elif hit_frame == -1:
            return Orientation.REVERSE

    def to_pm_one_encoding(self):
        return +1 if self.FORWARD else -1

    def __neg__(self):
        if self is self.FORWARD:
            return self.REVERSE
        else:
            return self.FORWARD


class Contig:
    def __init__(self, contig_id, contig_length, orientation=Orientation.FORWARD):
        self.contig_id: str = contig_id
        self.contig_length: int = contig_length
        self.contig_orientation: Orientation = orientation


class Feature:
    def __init__(self, start, end, frame, contig: Contig):
        self.start: int = start
        self.end: int = end
        self.frame: Orientation = frame
        self.contig: Contig = contig


class PipolinFragment:
    def __init__(self, contig, start, end):
        self.contig: Contig = contig
        self.start: int = start
        self.end: int = end
        self.atts: MutableSequence[Feature] = []


class GQuery:
    def __init__(self, gquery_id):
        self.gquery_id: str = gquery_id
        self.contigs: MutableSequence[Contig] = []
        self.polbs: MutableSequence[Feature] = []
        self.atts: MutableSequence[Feature] = []
        self.trnas: MutableSequence[Feature] = []
        self.target_trnas: MutableSequence[Feature] = []
        self.denovo_atts: MutableSequence[Feature] = []
        self.pipolin_fragments: MutableSequence[PipolinFragment] = []

    def get_contig_by_id(self, contig_id) -> Contig:
        for contig in self.contigs:
            if contig.contig_id == contig_id:
                return contig

    def get_features_of_contig(self, contig_id, feature_type: str) -> MutableSequence[Feature]:
        features_to_return = []
        features = self.get_features_by_type(feature_type=feature_type)
        for feature in features:
            if feature.contig.contig_id == contig_id:
                features_to_return.append(feature)

        return features_to_return

    def get_features_by_type(self, feature_type: str) -> MutableSequence[Feature]:
        if feature_type == 'polbs':
            return self.polbs
        elif feature_type == 'atts':
            return self.atts
        elif feature_type == 'trnas':
            return self.trnas
        else:
            raise AssertionError('Feature type can be: "polbs", "atts" or "trnas"!')

    def feature_from_blasthit(self, hit, contig_id):
        return Feature(start=hit.hit_start, end=hit.hit_end,
                       frame=Orientation.orientation_from_blast(hit.hit_frame),
                       contig=self.get_contig_by_id(contig_id))

    def define_target_trnas(self):
        for trna in self.trnas:
            for att in self.atts:
                if att.contig.contig_id == trna.contig.contig_id:
                    if self._is_overlapping(range1=(att.start, att.end), range2=(trna.start, trna.end)):
                        self.target_trnas.append(trna)

    def is_single_contig(self):
        return len(self.contigs) == 1

    def is_on_the_same_contig(self):
        target_contigs = []
        target_contigs.extend(i.contig.contig_id for i in self.polbs)
        target_contigs.extend(i.contig.contig_id for i in self.atts)
        target_contigs.extend(i.contig.contig_id for i in self.target_trnas)
        return len(set(target_contigs)) == 1

    def get_left_right_windows(self):
        # `find_atts_denovo()`
        polymerases = sorted((i for i in self.polbs), key=lambda p: p.start)

        if polymerases[-1].start - polymerases[0].start > 10000:   # TODO: is it small/big enough?
            raise AssertionError(f'You have several piPolBs per genome and they are too far from each other: '
                                 f'within the region ({polymerases[0].start}, {polymerases[-1].end}). It might be, '
                                 f'that you have two or more pipolins per genome, but we are expecting only one.')

        if self.is_single_contig():
            length = self.contigs[0].contig_length
            left_edge = polymerases[0].start - 100000
            left_window = (left_edge if left_edge >= 0 else 0, polymerases[0].start)
            right_edge = polymerases[-1].end + 100000
            right_window = (polymerases[-1].end, right_edge if right_edge <= length else length)

            return left_window, right_window
        else:
            raise AssertionError('This method is only for complete genomes!')

    def is_att_denovo(self, left_repeat, right_repeat):
        # `find_atts_denovo()`
        if self._is_overlapping_att(left_repeat=left_repeat):
            return False
        return self._is_overlapping_trna(left_repeat=left_repeat, right_repeat=right_repeat)

    def _is_overlapping_att(self, left_repeat):
        # `find_atts_denovo()`
        for att in self.atts:
            if self._is_overlapping(left_repeat, (att.start, att.end)):
                return True

        return False

    def _is_overlapping_trna(self, left_repeat, right_repeat):
        # `find_atts_denovo()`
        for trna in self.trnas:
            trna_range = (trna.start, trna.end)
            if self._is_overlapping(left_repeat, trna_range) or self._is_overlapping(right_repeat, trna_range):
                return True

        return False

    def get_pipolin_bounds(self, long):
        # When scaffolding is not required
        polymerases = sorted((i for i in self.polbs), key=lambda p: p.start)
        atts = sorted((i for i in self.atts), key=lambda p: p.start)

        if not self._is_polymerase_inside(atts, polymerases):
            raise AssertionError('The polymerase are not within att bounds!')

        if len(atts) > 3:
            raise AssertionError('There are more than 3 atts!')

        if long:
            if len(atts) == 3:
                return atts[0].start - 50, atts[2].end + 50
            else:
                return atts[0].start - 50, atts[1].end + 50
        else:
            if polymerases[-1].end < atts[1].end:
                return atts[0].start - 50, atts[1].end + 50
            else:
                return atts[1].start - 50, atts[2].end + 50

    @staticmethod
    def _is_overlapping(range1, range2):
        max_start = max(range1[0], range2[0])
        min_end = min(range1[1], range2[1])
        return max_start <= min_end

    @staticmethod
    def _is_polymerase_inside(atts, polymerases):
        return atts[0].start < polymerases[0].start and polymerases[-1].end < atts[-1].end

    # OLD FUNCTIONS!!!
    @staticmethod
    def _dict_by_node_normalized(features):
        return {node: sorted(list(ps), key=lambda p: p.start) for node, ps
                in groupby((i.normalize() for i in features), key=lambda x: x.node)}

    @classmethod
    def _add_dangling_atts(cls, atts, things_to_return, length_by_contig):
        for node, features in atts.items():
            if node not in things_to_return:
                things_to_return[node] = cls._get_dangling_feature(features, node, length_by_contig)

    @classmethod
    def _get_dangling_feature(cls, features, node, length_by_contig):
        pos_left = features[0].start - 50000
        left = pos_left if pos_left > 0 else 0
        pos_right = features[-1].end + 50000
        right = pos_right if pos_right < length_by_contig[node] else -1
        return left, right

    @classmethod
    def _get_pipolin_two_atts(cls, atts, polymerases, long):
        if cls._is_polymerase_inside(atts, polymerases):
            if long:
                return atts[0].start - 50, atts[2].end + 50
            else:
                return atts[0].start - 50, atts[1].end + 50
        else:
            raise AssertionError('The polymerases are not within att bounds!')

    @classmethod
    def _get_pipolin_single_att(cls, att: Feature, polymerases, length_by_contig):
        # TODO: there is an error in my assumptions!
        # When there is an att and a polymerase on one contig, I cannot just cut the sequence
        # from the outer side of att +/-50, as there can be an additional upstream att, but
        # on a separate contig.
        if att.end < polymerases[0].start:
            left = att.start - 50
            pos_right = polymerases[-1].end + 50000
            right = pos_right if pos_right < length_by_contig[att.node] else -1
        else:
            pos_left = polymerases[0].start - 50000
            left = pos_left if pos_left > 0 else 0
            right = att.end + 50
        return left, right

    def get_contigs_with_bounds(self, length_by_contig, long):
        # TODO: check boarders when +/-50 nt !
        # TODO: check indexes carefully: 0-based and 1-based !!!
        polymerases = self._dict_by_node_normalized(self.polbs)
        atts = self._dict_by_node_normalized(self.atts)

        things_to_return = {}
        for node, features in polymerases.items():
            self.add_pipolin_for_node(atts, features, node, things_to_return, length_by_contig, long)

        return things_to_return

    def add_pipolin_for_node(self, atts, features, node, things_to_return, length_by_contig, long):
        if node in atts:
            if long:
                if len(atts[node]) == 3:
                    things_to_return[node] = self._get_pipolin_two_atts(atts[node], features, long)
            if len(atts[node]) == 2:
                things_to_return[node] = self._get_pipolin_two_atts(atts[node], features, long=False)
            elif len(atts[node]) == 1:
                things_to_return[node] = self._get_pipolin_single_att(atts[node][0], features, length_by_contig)
                self._add_dangling_atts(atts, things_to_return, length_by_contig)
            else:
                raise AssertionError('More than two atts on one contig!')
        else:
            things_to_return[node] = self._get_dangling_feature(features, node, length_by_contig)
            self._add_dangling_atts(atts, things_to_return, length_by_contig)


def ncbi_acc_download(acc_ids):
    for acc_id in acc_ids:
        print(f'Downloading {acc_id}...')
        # TODO: fix local path!
        subprocess.run(['/home/liubov/.local/bin/ncbi-acc-download', '-F', 'fasta', acc_id])


def read_blastxml(blast_xml):
    return SearchIO.read(blast_xml, 'blast-xml')


def save_to_shelve(out_file, object_to_save, key_name):
    shelve_db = shelve.open(os.path.splitext(out_file)[0])
    shelve_db[key_name] = object_to_save
    shelve_db.close()


def read_from_shelve(shelve_file, key_name):
    shelve_db = shelve.open(os.path.splitext(shelve_file)[0])
    shelve_object = shelve_db[key_name]
    shelve_db.close()
    return shelve_object


# def save_as_csv(Pipolin, pipolins, out_file):
#     # TODO: rewrite this!
#     with open(out_file, 'w') as ouf:
#         print(','.join(Pipolin._fields), file=ouf)
#         for i_p, _ in enumerate(pipolins):
#             words = ['None' if i is None else str(i) for i in pipolins[i_p]]
#             print(','.join(words), file=ouf)


def get_roary_groups(roary_dir) -> Mapping[str, Mapping[str, Sequence[str]]]:
    roary_groups = {}
    with open(os.path.join(roary_dir, 'gene_presence_absence.csv')) as csv_file:
        reader = csv.reader(csv_file, delimiter=',')
        header = next(reader)
        for entry in reader:
            group_name = entry[0]
            genes = {genome: genes.split('\t') if genes else [] for genome, genes in zip(header[14:],entry[14:])}
            roary_groups[group_name] = genes
    return roary_groups


def define_gquery_id(genome):
    return os.path.splitext(os.path.basename(genome))[0]


def blast_genome_against_seq(genome, seq, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    with open(os.path.join(output_dir, define_gquery_id(genome=genome) + '.fmt5'), 'w') as ouf:
        subprocess.run(['blastn', '-query', seq, '-subject', genome, '-outfmt', '5'], stdout=ouf)


SeqIORecords = MutableMapping[str, SeqIO.SeqRecord]


def read_seqio_records(file, file_format) -> SeqIORecords:
    file_formats = ['genbank', 'fasta']
    if file_format not in file_formats:
        raise AssertionError('Only genbank or fasta formats are acceptable!')

    records = SeqIO.to_dict(SeqIO.parse(file, format=file_format))
    return records


def run_aragorn(genome, aragorn_results):
    os.makedirs(aragorn_results, exist_ok=True)
    with open(os.path.join(aragorn_results, define_gquery_id(genome=genome) + '.batch'), 'w') as ouf:
        subprocess.run(['aragorn', '-l', '-w', genome], stdout=ouf)


def write_genbank_records(gb_records: SeqIORecords, out_dir, gquery):
    records = [record for record in gb_records.values()]
    with open(os.path.join(out_dir, f'{gquery.gquery_id}.gbk'), 'w') as ouf:
        SeqIO.write(records, ouf, 'genbank')


def write_gff_records(in_records, out_dir, gquery):
    records = [record for record in in_records.values()]
    with open(os.path.join(out_dir, f'{gquery.gquery_id}.gff'), 'w') as ouf:
        GFF.write(records, ouf)
        print('##FASTA', file=ouf)
        SeqIO.write(records, ouf, format='fasta')


def write_fna_records(gb_records, out_dir):
    for key, value in gb_records.items():
        if len(value) > 1:
            raise AssertionError('The only record is expected!')
        record = list(value.values()).pop()
        record.id = key
        with open(os.path.join(out_dir, f'{key}.fa'), 'w') as ouf:
            SeqIO.write(record, ouf, format='fasta')


def read_from_prokka_dir(prokka_dir, ext):
    files = []
    for file in os.listdir(prokka_dir):
        if file.endswith(ext):
            files.append(os.path.join(prokka_dir, file))
    records = {}
    for file in files:
        records.update(SeqIO.to_dict(SeqIO.parse(file, 'fasta')))
    return records


def read_aragorn_batch(aragorn_batch):
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


def save_left_right_subsequences(genome, left_window, right_window, repeats_dir):
    os.makedirs(repeats_dir, exist_ok=True)
    genome_seq = SeqIO.read(handle=genome, format='fasta')

    left_seq = genome_seq[left_window[0]:left_window[1]]
    right_seq = genome_seq[right_window[0]:right_window[1]]
    SeqIO.write(sequences=left_seq, format='fasta',
                handle=os.path.join(repeats_dir, define_gquery_id(genome=genome) + '.left'))
    SeqIO.write(sequences=right_seq, format='fasta',
                handle=os.path.join(repeats_dir, define_gquery_id(genome=genome) + '.right'))


def blast_for_identical(gquery_id, repeats_dir):
    with open(os.path.join(repeats_dir, gquery_id + '.fmt5'), 'w') as ouf:
        subprocess.run(['blastn', '-query', os.path.join(repeats_dir, gquery_id + '.left'),
                        '-subject', os.path.join(repeats_dir, gquery_id + '.right'),
                        '-outfmt', '5', '-perc_identity', '100', '-word_size', '12',   # TODO: try smaller word_size!
                        '-strand', 'plus'], stdout=ouf)


def extract_repeats(file):
    repeats = read_blastxml(file)
    left_repeats = []
    right_repeats = []
    for entry in repeats:
        for hit in entry:
            left_repeats.append((hit.query_start, hit.query_end))
            right_repeats.append((hit.hit_start, hit.hit_end))

    return left_repeats, right_repeats


def set_proper_location(seq_range, shift):
    return seq_range[0] + shift, seq_range[1] + shift
