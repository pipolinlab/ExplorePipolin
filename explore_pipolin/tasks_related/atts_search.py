import os
from collections import defaultdict
from itertools import chain
from random import random
from typing import List, Mapping, Union

from Bio import SeqIO
from prefect import task

from explore_pipolin.common import Genome, FeatureType, Range, PairedLocation, Strand, AttFeature, ContigID, \
    MultiLocation
from explore_pipolin.tasks_related.misc import get_ranges_around_pipolbs
from explore_pipolin.utilities.external_tools import blastn_against_ref_att, blast_for_repeats
from explore_pipolin.utilities.io import read_blastxml
from explore_pipolin.utilities.logging import genome_specific_logging


@task()
@genome_specific_logging
def find_atts(genome: Genome, out_dir: str) -> Genome:
    atts_dir = os.path.join(out_dir, 'atts_search')
    os.makedirs(atts_dir, exist_ok=True)

    finder = AttFinder(genome=genome, output_dir=atts_dir)

    entries = finder.blast_for_atts()
    finder.add_att_features_from_blast_entries(entries)

    return genome


@task()
@genome_specific_logging
def find_atts_denovo(genome: Genome, out_dir: str) -> Genome:
    atts_denovo_dir = os.path.join(out_dir, 'atts_denovo_search')
    os.makedirs(atts_denovo_dir, exist_ok=True)

    finder = AttDenovoFinder(genome=genome, output_dir=atts_denovo_dir)

    repeats: List[MultiLocation] = finder.find_repeats()
    finder.write_repeats(repeats)

    atts_denovo: List[MultiLocation] = [rep for rep in repeats if finder.is_att_denovo(rep)]
    finder.write_atts_denovo(atts_denovo)

    finder.extend_att_features(atts_denovo)
    finder.extend_target_trna_features()

    return genome


class AttFinder:
    def __init__(self, genome: Genome, output_dir: str):
        self.genome = genome
        self.output_dir = output_dir

    def blast_for_atts(self):
        output_file = os.path.join(self.output_dir, self.genome.id + '.fmt5')
        blastn_against_ref_att(genome_file=self.genome.file, output_file=output_file)
        return read_blastxml(blast_xml=output_file)

    def add_att_features_from_blast_entries(self, entries):
        repeat_id = random()
        for entry in entries:
            for hit in entry:
                feature = self._create_att_feature(hit=hit, contig_id=entry.id, repeat_id=repeat_id)
                self.genome.features.add_features(feature, feature_type=FeatureType.ATT)

    def _create_att_feature(self, hit, contig_id: str, repeat_id) -> AttFeature:
        return AttFeature(location=Range(start=hit.hit_start, end=hit.hit_end),
                          strand=Strand.from_pm_one_encoding(hit.hit_strand),
                          contig_id=contig_id, genome=self.genome, repeat_id=repeat_id)


class AttDenovoFinder:
    def __init__(self, genome: Genome, output_dir: str):
        self.genome = genome
        self.output_dir = output_dir

    def find_repeats(self) -> List[MultiLocation]:
        ranges_around_pipolbs = get_ranges_around_pipolbs(self.genome)
        self._save_seqs_around_pipolbs(ranges_around_pipolbs)
        blast_for_repeats(ranges_around_pipolbs, genome_id=self.genome.id, repeats_dir=self.output_dir)
        paired_repeats: List[PairedLocation] = self._extract_repeats(ranges_around_pipolbs)
        repeats = self._regroup_paired_repeats(paired_repeats)
        return repeats

    def write_repeats(self, repeats: List[MultiLocation]):
        with open(os.path.join(self.output_dir, self.genome.id + '.repeats'), 'w') as ouf:
            for loc in repeats:
                ranges = [f'({i.start},{i.end})' for i in loc.ranges]
                print(loc.contig_id, ' '.join(ranges), sep='\t', file=ouf)

    def is_att_denovo(self, repeat: Union[MultiLocation, PairedLocation]) -> bool:
        if isinstance(repeat, PairedLocation):
            return self.is_att_denovo(MultiLocation(
                ranges=[repeat.left, repeat.right], contig_id=repeat.contig_id))
        atts_of_contig = self.genome.features.atts_dict()[repeat.contig_id]
        for att in atts_of_contig:
            if repeat.ranges[0].is_overlapping(att):
                return False

        trnas_of_contig = self.genome.features.trnas_dict()[repeat.contig_id]
        for trna in trnas_of_contig:
            if trna.location.is_overlapping_any(repeat.ranges):
                return True
        return False

    def extend_att_features(self, atts_denovo: List[MultiLocation]):
        for att in atts_denovo:
            repeat_id = random()
            for r in att.ranges:
                new_att = AttFeature(r, Strand.FORWARD, att.contig_id, self.genome, repeat_id)
                self.genome.features.add_features(new_att, feature_type=FeatureType.ATT)

    def extend_target_trna_features(self):
        target_trnas_dict = self.genome.features.target_trnas_dict()
        for att in self.genome.features.get_features(FeatureType.ATT):
            target_trna = self.genome.features.get_features(FeatureType.TRNA).get_overlapping(att)
            if target_trna is not None:
                if target_trna not in target_trnas_dict[target_trna.contig_id]:
                    self.genome.features.add_features(target_trna, feature_type=FeatureType.TARGET_TRNA)

    def write_atts_denovo(self, atts_denovo: List[MultiLocation]):
        with open(os.path.join(self.output_dir, self.genome.id + '.atts_denovo'), 'w') as ouf:
            print('att_start', 'att_end', sep='\t', file=ouf)
            for att in atts_denovo:
                for r in att.ranges:
                    print(r.start, r.end, sep='\t', file=ouf)

    def _save_seqs_around_pipolbs(self, ranges_around_pipolbs: List[PairedLocation]):
        genome_seq = SeqIO.read(handle=self.genome.file, format='fasta')

        for i, range_pair in enumerate(ranges_around_pipolbs):
            left_seq = genome_seq[range_pair.left.start:range_pair.left.end]
            right_seq = genome_seq[range_pair.right.start:range_pair.right.end]
            SeqIO.write(sequences=left_seq, format='fasta',
                        handle=os.path.join(self.output_dir, self.genome.id + f'_{i}.left'))
            SeqIO.write(sequences=right_seq, format='fasta',
                        handle=os.path.join(self.output_dir, self.genome.id + f'_{i}.right'))

    def _extract_repeats(self, ranges_around_pipolbs: List[PairedLocation]) -> List[PairedLocation]:
        repeats = []
        for i, range_pair in enumerate(ranges_around_pipolbs):
            repeats_xml = read_blastxml(os.path.join(self.output_dir, self.genome.id + f'_{i}.fmt5'))
            for entry in repeats_xml:
                for hit in entry:
                    left_repeat = Range(start=hit.query_start, end=hit.query_end).shift(range_pair.left.start)
                    right_repeat = Range(start=hit.hit_start, end=hit.hit_end).shift(range_pair.right.start)
                    repeat_pair = PairedLocation(left_repeat, right_repeat, range_pair.contig_id)
                    repeats.append(repeat_pair)
        return repeats

    @staticmethod
    def _regroup_paired_repeats(paired_repeats: List[PairedLocation]) -> List[MultiLocation]:
        result: Mapping[ContigID, List[MultiLocation]] = defaultdict(list)
        for repeat in paired_repeats:
            for rs in result[repeat.contig_id]:
                if repeat.left.is_overlapping_any(rs) and not repeat.right.is_overlapping_any(rs.ranges):
                    rs.ranges.append(repeat.right)
                    break
                elif repeat.right.is_overlapping_any(rs) and not repeat.left.is_overlapping_any(rs.ranges):
                    rs.ranges.append(repeat.left)
                    break
            else:
                result[repeat.contig_id].append(MultiLocation(
                    ranges=[repeat.right, repeat.left], contig_id=repeat.contig_id))
        return list(chain(*result.values()))
