import csv
import os
from collections import defaultdict
from enum import Enum, auto
from typing import Sequence, MutableSequence, Mapping, MutableMapping, Set
from itertools import groupby
import subprocess
from random import randrange

from BCBio import GFF
from Bio import SearchIO, SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

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
        self.target_trnas_denovo: MutableSequence[Feature] = []
        self.pipolin_fragments: MutableSequence[PipolinFragment] = []

    def get_contig_by_id(self, contig_id) -> Contig:
        for contig in self.contigs:
            if contig.contig_id == contig_id:
                return contig

    # TODO: replace it with _dict_by_contig_normalized()!
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
        elif feature_type == 'target_trnas':
            return self.target_trnas
        else:
            raise AssertionError('Feature type can be: "polbs", "atts" or "target_trnas"!')

    def feature_from_blasthit(self, hit, contig_id):
        return Feature(start=hit.hit_start, end=hit.hit_end,
                       frame=Orientation.orientation_from_blast(hit.hit_frame),
                       contig=self.get_contig_by_id(contig_id))

    # `add_features_from_aragorn` and `add_features_atts_denovo`
    def define_target_trna(self, att: Feature):
        for trna in self.trnas:
            if att.contig.contig_id == trna.contig.contig_id:
                if self._is_overlapping(range1=(att.start, att.end), range2=(trna.start, trna.end)):
                    return trna

    def is_single_contig(self):
        return len(self.contigs) == 1

    def is_on_the_same_contig(self):
        target_contigs = []
        target_contigs.extend(i.contig.contig_id for i in self.polbs)
        target_contigs.extend(i.contig.contig_id for i in self.atts)
        target_contigs.extend(i.contig.contig_id for i in self.target_trnas)
        return len(set(target_contigs)) == 1

    # `find_atts_denovo`
    def get_left_right_windows(self):
        polymerases = sorted((i for i in self.polbs), key=lambda p: p.start)
        atts = sorted((i for i in self.atts), key=lambda p: p.start)

        if len(atts) != 0:
            if not self._is_polymerase_inside(atts=atts, polymerases=polymerases):
                raise AssertionError('The piPolBs are not within att bounds!')
        else:
            if polymerases[-1].start - polymerases[0].start > 10000:   # TODO: is it small/big enough?
                raise AssertionError(f'You have several piPolBs per genome and they are too far from each other: '
                                     f'within the region ({polymerases[0].start}, {polymerases[-1].end}). It might be, '
                                     f'that you have two or more pipolins per genome, but we are expecting only one.')

        length = self.contigs[0].contig_length
        left_edge = polymerases[0].start - 100000
        left_window = (left_edge if left_edge >= 0 else 0, polymerases[0].start)
        right_edge = polymerases[-1].end + 100000
        right_window = (polymerases[-1].end, right_edge if right_edge <= length else length)

        return left_window, right_window

    def is_att_denovo(self, left_repeat, right_repeat):
        if self._is_overlapping_att(left_repeat=left_repeat):
            return False
        return self._is_overlapping_trna(left_repeat=left_repeat, right_repeat=right_repeat)

    def _is_overlapping_att(self, left_repeat):
        for att in self.atts:
            if self._is_overlapping(left_repeat, (att.start, att.end)):
                return True

        return False

    def _is_overlapping_trna(self, left_repeat, right_repeat):
        for trna in self.trnas:
            trna_range = (trna.start, trna.end)
            if self._is_overlapping(left_repeat, trna_range) or self._is_overlapping(right_repeat, trna_range):
                return True

        return False

    # `analyse_pipolin_orientation`
    def is_single_target_trna_per_contig(self):
        # there was one case with two target trnas per genome, although usually only one
        targeted_contigs = [trna.contig.contig_id for trna in self.target_trnas]
        if len(self.target_trnas) != len(targeted_contigs):
            raise AssertionError("We are expecting a single tRNA to overlap with a single att per contig!")

    # `analyse_pipolin_orientation`
    def set_contig_orientation(self, contig: Contig):
        target_trnas = self.get_features_of_contig(contig_id=contig.contig_id, feature_type='target_trnas')
        atts = self.get_features_of_contig(contig_id=contig.contig_id, feature_type='atts')
        atts_frames = [att.frame for att in atts]
        polbs = self.get_features_of_contig(contig_id=contig.contig_id, feature_type='polbs')
        polbs_frames = [polb.frame for polb in polbs]

        if len(target_trnas) != 0:
            if len(set(atts_frames)) != 1:
                raise AssertionError('ATTs are expected to be in the same frame, as they are direct repeats!')
            if set(atts_frames).pop() == target_trnas[0].frame:
                raise AssertionError('ATT and tRNA are expected to be on the different strands!')
            contig.contig_orientation = - target_trnas[0].frame

        elif len(atts) != 0:
            if len(set(atts_frames)) != 1:
                raise AssertionError('ATTs are expected to be in the same frame, as they are direct repeats!')
            contig.contig_orientation = atts[0].frame

        if len(polbs) != 0:
            if len(set(polbs_frames)) != 1:  # an ambiguous case
                return
            contig.contig_orientation = polbs[0].frame

    # scaffolding is not required
    def get_pipolin_bounds(self):
        polymerases = sorted((i for i in self.polbs), key=lambda p: p.start)
        atts = sorted((i for i in self.atts), key=lambda p: p.start)

        length = polymerases[0].contig.contig_length
        left_edge = atts[0].start - 50 if atts[0].start < polymerases[0].start else polymerases[0].start - 50
        right_edge = atts[-1].end + 50 if atts[-1].end > polymerases[-1].end else polymerases[-1].end + 50
        return left_edge if left_edge >= 0 else 0, right_edge if right_edge <= length else length

    # scaffolding is required
    def try_creating_single_record(self):
        polbs_dict = self._dict_by_contig_normalized(self.polbs)
        atts_dict = self._dict_by_contig_normalized(self.atts)
        target_trnas_dict = self._dict_by_contig_normalized(self.target_trnas)

        unchangeable_contigs = self._get_unchangeable_contigs()

        if len(unchangeable_contigs) == 1:
            unchangeable_contig = unchangeable_contigs[0]
            print(f'The unchangeable contig is {unchangeable_contig.contig_id}!')

            att_only_contigs = self._order_att_only_contigs()
            if len(att_only_contigs) == 1:
                direction = self._get_direction_of_unchangeable(
                    polbs_sorted=polbs_dict[unchangeable_contig.contig_id],
                    atts_sorted=atts_dict[unchangeable_contig.contig_id])

                orphan_atts = atts_dict[att_only_contigs[0].contig_id]

                if direction == 'none':
                    contig_length = unchangeable_contig.contig_length
                    left_edge = atts_dict[unchangeable_contig.contig_id][0].start - 50
                    right_edge = atts_dict[unchangeable_contig.contig_id][-1].end + 50
                    pipolin = PipolinFragment(contig=self.get_contig_by_id(unchangeable_contig.contig_id),
                                              start=left_edge if left_edge >= 0 else 0,
                                              end=right_edge if right_edge <= contig_length else contig_length)
                    pipolin.atts.extend(atts_dict[unchangeable_contig.contig_id])
                    self.pipolin_fragments.append(pipolin)
                    # TODO: if att_only_contig has a target_trna, it could be added on the right
                elif direction == 'right':
                    left_fragment = self._create_unchangeable_fragment(atts_dict, unchangeable_contig, direction)
                    self.pipolin_fragments.append(left_fragment)

                    right_fragment = self._create_att_contig_fragment(direction=direction, atts=orphan_atts)
                    self.pipolin_fragments.append(right_fragment)
                elif direction == 'left':
                    left_fragment = self._create_att_contig_fragment(direction=direction, atts=orphan_atts)
                    self.pipolin_fragments.append(left_fragment)

                    right_fragment = self._create_unchangeable_fragment(atts_dict, unchangeable_contig, direction)
                    self.pipolin_fragments.append(right_fragment)
            elif len(att_only_contigs) == 2:
                # TODO: the order can be also [middle_fragment, left_fragment, right_fragment]
                middle_fragment = PipolinFragment(contig=self.get_contig_by_id(unchangeable_contig.contig_id),
                                                  start=0, end=unchangeable_contig.contig_length)
                middle_fragment.atts.extend(atts_dict[unchangeable_contig.contig_id])

                left_atts = atts_dict[att_only_contigs[0].contig_id]
                left_fragment = self._create_att_contig_fragment(direction='right', atts=left_atts)

                right_atts = atts_dict[att_only_contigs[1].contig_id]
                right_fragment = self._create_att_contig_fragment(direction='left', atts=right_atts)

                self.pipolin_fragments.extend([left_fragment, middle_fragment, right_fragment])
        elif len(unchangeable_contigs) == 2:
            target_trnas_contig0 = target_trnas_dict[unchangeable_contigs[0].contig_id]
            target_trnas_contig1 = target_trnas_dict[unchangeable_contigs[1].contig_id]

            if len(target_trnas_contig0) != 0:
                right_contig = unchangeable_contigs[0]
                left_contig = unchangeable_contigs[1]
            elif len(target_trnas_contig1) != 0:
                right_contig = unchangeable_contigs[1]
                left_contig = unchangeable_contigs[0]
            else:
                raise NotImplementedError

            left_direction = self._get_direction_of_unchangeable(polbs_sorted=polbs_dict[left_contig.contig_id],
                                                                 atts_sorted=atts_dict[left_contig.contig_id])
            left_fragment = self._create_unchangeable_fragment(atts_dict, left_contig, left_direction)

            right_direction = self._get_direction_of_unchangeable(polbs_sorted=polbs_dict[right_contig.contig_id],
                                                                  atts_sorted=atts_dict[right_contig.contig_id])
            right_fragment = self._create_unchangeable_fragment(atts_dict, right_contig, right_direction)

            self.pipolin_fragments.extend([left_fragment, right_fragment])
        else:
            self._try_finish_separate(atts_dict)

    def _create_unchangeable_fragment(self, atts_dict, contig, direction):
        if direction == 'none':
            raise NotImplementedError

        if direction == 'right':
            edge = atts_dict[contig.contig_id][0].start - 50
            fragment = PipolinFragment(contig=self.get_contig_by_id(contig.contig_id),
                                       start=edge if edge >= 0 else 0, end=contig.contig_length)
        else:
            edge = atts_dict[contig.contig_id][-1].end + 50
            fragment = PipolinFragment(contig=self.get_contig_by_id(contig.contig_id), start=0,
                                       end=edge if edge <= contig.contig_length else contig.contig_length)

        fragment.atts.extend(atts_dict[contig.contig_id])
        return fragment

    @staticmethod
    def _create_att_contig_fragment(direction: str, atts: MutableSequence[Feature]):
        contig_length = atts[0].contig.contig_length
        if direction == 'left':
            if atts[0].contig.contig_orientation == Orientation.FORWARD:
                left_edge = 0
                right_edge = atts[-1].end + 50
            else:
                left_edge = atts[0].start - 50
                right_edge = atts[0].contig.contig_length
        else:
            if atts[0].contig.contig_orientation == Orientation.FORWARD:
                left_edge = atts[0].start - 50
                right_edge = contig_length
            else:
                left_edge = 0
                right_edge = atts[-1].end + 50

        fragment = PipolinFragment(contig=atts[0].contig, start=left_edge if left_edge >= 0 else 0,
                                   end=right_edge if right_edge <= contig_length else contig_length)
        fragment.atts.extend(atts)
        return fragment

    def _try_finish_separate(self, atts_dict):
        polbs_fragments = self._create_polbs_fragments()

        if len(self.target_trnas) > 1 or len(self.target_trnas) == 0:
            raise NotImplementedError

        right_fragment = self._create_att_contig_fragment(direction='left',
                                                          atts=atts_dict[self.target_trnas[0].contig.contig_id])

        atts_contigs = set([i.contig for i in self.atts])
        if len(atts_contigs) == 2:
            print('The single record can be created!!!\n')

            atts_contigs = list(atts_contigs)
            if atts_contigs[0].contig_id == self.target_trnas[0].contig.contig_id:
                left_contig = atts_contigs[1]
            else:
                left_contig = atts_contigs[0]

            left_fragment = self._create_att_contig_fragment(direction='right',
                                                             atts=atts_dict[left_contig.contig_id])
        else:
            raise NotImplementedError

        self.pipolin_fragments.append(left_fragment)
        self.pipolin_fragments.extend(polbs_fragments)
        self.pipolin_fragments.append(right_fragment)

    def _create_polbs_fragments(self) -> MutableSequence[PipolinFragment]:
        polbs_contigs = set([i.contig for i in self.polbs])
        if len(polbs_contigs) == 2:
            print('Two "polb only contigs" were found!')
            polbs_contigs = list(polbs_contigs)
            if polbs_contigs[0].contig_length + polbs_contigs[1].contig_length < 15000:
                polb0_fragment = PipolinFragment(contig=self.get_contig_by_id(polbs_contigs[0].contig_id), start=0,
                                                 end=polbs_contigs[0].contig_length)
                polb1_fragment = PipolinFragment(contig=self.get_contig_by_id(polbs_contigs[1].contig_id), start=0,
                                                 end=polbs_contigs[1].contig_length)

                polb0_length = sum((i.end - i.start) for i in self.get_features_of_contig(
                    contig_id=polbs_contigs[0].contig_id,
                    feature_type='polbs'))
                polb1_length = sum((i.end - i.start) for i in self.get_features_of_contig(
                    contig_id=polbs_contigs[1].contig_id,
                    feature_type='polbs'))
                # TODO: comparing just by length is an unreliable way! REWRITE if possible!
                if polb0_length < polb1_length:
                    polbs_fragments = [polb0_fragment, polb1_fragment]
                else:
                    polbs_fragments = [polb1_fragment, polb0_fragment]
            else:
                raise NotImplementedError

        elif len(polbs_contigs) == 1:
            polbs_fragments = [PipolinFragment(contig=self.get_contig_by_id(self.polbs[0].contig.contig_id),
                                               start=0, end=self.polbs[0].contig.contig_length)]
        else:
            raise NotImplementedError

        return polbs_fragments

    @staticmethod
    def _get_direction_of_unchangeable(polbs_sorted, atts_sorted):
        if polbs_sorted[0].start > atts_sorted[-1].end:
            return 'right'
        elif polbs_sorted[-1].end < atts_sorted[0].start:
            return 'left'
        else:
            return 'none'

    def _get_unchangeable_contigs(self) -> MutableSequence[Contig]:
        polbs_contigs = set([i.contig for i in self.polbs])
        if len(polbs_contigs) > 2:
            raise NotImplementedError

        contigs_to_return = []
        for contig in polbs_contigs:
            atts_next_polbs = self.get_features_of_contig(contig_id=contig.contig_id, feature_type='atts')
            if len(atts_next_polbs) != 0:
                contigs_to_return.append(contig)
        return contigs_to_return

    def _order_att_only_contigs(self) -> MutableSequence[Contig]:
        att_only_contigs = self._get_att_only_contigs()

        if len(att_only_contigs) == 1:
            print('The single record can be created!!!\n')
            return [att_only_contigs.pop()]
        elif len(att_only_contigs) == 2:
            contigs = list(att_only_contigs)

            if len(self.get_features_of_contig(contig_id=contigs[0].contig_id, feature_type='target_trnas')) != 0:
                print('The single record can be created!!!\n')
                return [contigs[1], contigs[0]]
            elif len(self.get_features_of_contig(contig_id=contigs[1].contig_id, feature_type='target_trnas')) != 0:
                print('The single record can be created!!!\n')
                return [contigs[0], contigs[1]]
            else:
                raise NotImplementedError

        else:
            raise NotImplementedError

    def _get_att_only_contigs(self) -> Set[Contig]:
        att_only_contigs = set()
        for att in self.atts:
            polbs_next_att = self.get_features_of_contig(contig_id=att.contig.contig_id, feature_type='polbs')
            if len(polbs_next_att) == 0:
                att_only_contigs.add(att.contig)

        return att_only_contigs

    @staticmethod
    def _is_overlapping(range1, range2):
        max_start = max(range1[0], range2[0])
        min_end = min(range1[1], range2[1])
        return max_start <= min_end

    @staticmethod
    def _is_polymerase_inside(atts, polymerases):
        return atts[0].start < polymerases[0].start and polymerases[-1].end < atts[-1].end

    @staticmethod
    def _dict_by_contig_normalized(features):
        return {contig: sorted(list(ps), key=lambda p: p.start) for contig, ps
                in groupby((i for i in features), key=lambda x: x.contig.contig_id)}

    def create_att_feature(self, start: int, end: int, frame: Orientation, records_format: str) -> SeqFeature:
        random_number = randrange(10000, 99999)
        gb_qualifiers = {'inference': ['HMM:custom'], 'locus_tag': [f'{self.gquery_id}_{random_number}'],
                         'rpt_family': ['Att'], 'rpt_type': ['direct']}
        gff_qualifiers = {'phase': ['.'], 'source': ['HMM:custom'],
                          'ID': [f'{self.gquery_id}_{random_number}'], 'inference': ['HMM:custom'],
                          'locus_tag': [f'{self.gquery_id}_{random_number}'],
                          'rpt_family': ['Att'], 'rpt_type': ['direct']}
        att_feature = SeqFeature(type='repeat_region',
                                 location=FeatureLocation(start=start, end=end, strand=frame.to_pm_one_encoding()),
                                 qualifiers=gb_qualifiers if records_format == 'gb' else gff_qualifiers)
        return att_feature


def read_blastxml(blast_xml):
    return SearchIO.read(blast_xml, 'blast-xml')


# CHECK
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


def read_gff_records(file) -> SeqIORecords:
    gff_records = {}
    with open(file) as inf:
        for entry in GFF.parse(inf):
            gff_records[entry.id] = entry

    return gff_records


def run_aragorn(genome, aragorn_results):
    os.makedirs(aragorn_results, exist_ok=True)
    with open(os.path.join(aragorn_results, define_gquery_id(genome=genome) + '.batch'), 'w') as ouf:
        subprocess.run(['aragorn', '-l', '-w', genome], stdout=ouf)


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


# CHECK
def write_fna_records(gb_records, out_dir):
    for key, value in gb_records.items():
        if len(value) > 1:
            raise AssertionError('The only record is expected!')
        record = list(value.values()).pop()
        record.id = key
        with open(os.path.join(out_dir, f'{key}.fa'), 'w') as ouf:
            SeqIO.write(record, ouf, format='fasta')


# CHECK
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
                        '-outfmt', '5', '-perc_identity', '100', '-word_size', '6',
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
