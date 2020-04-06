import csv
import os
import shelve
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
        return Orientation.FORWARD if hit_frame == 1 else Orientation.REVERSE

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
        self.hallmarks: MutableSequence[Feature] = []


class GQuery:
    def __init__(self, gquery_id):
        self.gquery_id: str = gquery_id
        self.contigs: MutableSequence[Contig] = []
        self.polymerases: MutableSequence[Feature] = []
        self.atts: MutableSequence[Feature] = []
        self.trnas: MutableSequence[Feature] = []
        self.target_trnas: MutableSequence[Feature] = []
        self.atts_denovo: MutableSequence[Feature] = []
        self.pipolin_fragments: MutableSequence[PipolinFragment] = []

    def get_features_of_contig(self, contig_id, feature_type: str) -> MutableSequence[Feature]:
        features_to_return = []
        features = self.get_features_by_type(feature_type=feature_type)
        for feature in features:
            if feature.contig.contig_id == contig_id:
                features_to_return.append(feature)

        return features_to_return

    def get_contig_by_id(self, contig_id) -> Contig:
        for contig in self.contigs:
            if contig.contig_id == contig_id:
                return contig

    def get_features_by_type(self, feature_type: str) -> MutableSequence[Feature]:
        if feature_type == 'polymerases':
            return self.polymerases
        elif feature_type == 'atts':
            return self.atts
        elif feature_type == 'trnas':
            return self.trnas
        else:
            raise AssertionError('Feature type can be: "polymerases", "atts" or "trnas"!')

    def is_single_contig(self):
        return len(self.contigs) == 1

    def is_on_the_same_contig(self):
        target_contigs = []
        target_contigs.extend(i.contig.contig_id for i in self.polymerases)
        target_contigs.extend(i.contig.contig_id for i in self.atts)
        target_contigs.extend(i.contig.contig_id for i in self.target_trnas)
        return len(set(target_contigs)) == 1

    def get_left_right_windows(self):
        # `find_atts_denovo()`
        # TODO: how to be with several polymerases?
        #  * some check for polymerase integrity?
        #  * some check of how far polymerases are from each other?
        polymerases = sorted((i for i in self.polymerases), key=lambda p: p.start)
        if self.is_single_contig():
            contig_length = self.contigs[0].contig_length
            left_edge = polymerases[0].start - 100000
            left_window = (left_edge if left_edge >= 0 else 0, polymerases[0].start)
            right_edge = polymerases[-1].end + 100000
            right_window = (polymerases[-1].end, right_edge if right_edge <= contig_length else contig_length)

            return left_window, right_window
        else:
            raise AssertionError('This method is only for complete genomes!')

    def get_pipolin_bounds(self, long):
        # Scaffolding is not required
        polymerases = sorted((i for i in self.polymerases), key=lambda p: p.start)
        atts = sorted((i for i in self.atts), key=lambda p: p.start)

        if not self._is_polymerase_inside(atts, polymerases):
            raise AssertionError('The polymerases are not within att bounds!')

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
    def _is_polymerase_inside(atts, polymerases):
        return atts[0].start < polymerases[0].start and polymerases[-1].end < atts[-1].end

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
        polymerases = self._dict_by_node_normalized(self.polymerases)
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


def blast_genome_against_seq(genome, seq, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    with open(os.path.join(output_dir, f'{os.path.splitext(os.path.basename(genome))[0]}.fmt5'), 'w') as ouf:
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
    with open(os.path.join(aragorn_results, f'{os.path.splitext(os.path.basename(genome))[0]}.batch'), 'w') as ouf:
        subprocess.run(['aragorn', '-l', '-w', genome], stdout=ouf)


def write_genbank_records(gb_records, out_dir):
    for key, value in gb_records.items():
        records = [record for record in value.values()]
        with open(os.path.join(out_dir, f'{key}.gbk'), 'w') as ouf:
            SeqIO.write(records, ouf, 'genbank')


def write_gff_records(in_records, out_dir):
    for key, value in in_records.items():
        records = [record for record in value.values()]
        with open(os.path.join(out_dir, f'{key}.gff'), 'w') as ouf:
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


def feature_from_blasthit(hit, gquery, contig_id):
    return Feature(start=hit.hit_start, end=hit.hit_end,
                   frame=Orientation.orientation_from_blast(hit.hit_frame),
                   contig=gquery.get_contig_by_id(contig_id))
