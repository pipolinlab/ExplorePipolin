import os
from collections import defaultdict
from itertools import chain
from typing import List, Mapping, Sequence

from Bio import SeqIO
from prefect import task, context

from explore_pipolin.common import Genome, FeatureType, Range, PairedLocation, Strand, AttFeature, ContigID, \
    MultiLocation
from explore_pipolin.tasks_related.misc import get_ranges_around_pipolbs
from explore_pipolin.utilities.external_tools import blastn_against_ref_att, blast_for_repeats
from explore_pipolin.utilities.io import read_blastxml, create_seqio_records_dict
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
        att_id = self.genome.features.get_features(FeatureType.ATT).get_next_att_id()
        for entry in entries:
            att_features = self._create_att_features(entry, att_id)
            self.genome.features.add_features(*att_features, feature_type=FeatureType.ATT)

    def _create_att_features(self, entry, att_id) -> Sequence[AttFeature]:
        new_att_features = []
        for hit in entry:
            new_att_features.append(AttFeature(location=Range(start=hit.hit_start, end=hit.hit_end),
                                               strand=Strand.from_pm_one_encoding(hit.hit_strand),
                                               contig_id=entry.id, genome=self.genome, att_id=att_id))
        return new_att_features


class AttDenovoFinder:
    def __init__(self, genome: Genome, output_dir: str):
        self.genome = genome
        self.output_dir = output_dir

    def find_repeats(self) -> List[MultiLocation]:
        ranges_around_pipolbs = get_ranges_around_pipolbs(self.genome)
        self._save_seqs_around_pipolbs(ranges_around_pipolbs)
        blast_for_repeats(ranges_around_pipolbs, genome_id=self.genome.id, repeats_dir=self.output_dir)
        paired_repeats: List[PairedLocation] = self._extract_repeats(ranges_around_pipolbs)
        return self._regroup_paired_repeats(paired_repeats)

    def write_repeats(self, repeats: List[MultiLocation]):
        with open(os.path.join(self.output_dir, self.genome.id + '.repeats'), 'w') as ouf:
            for loc in repeats:
                ranges = [f'({i.start},{i.end})' for i in loc.ranges]
                print(loc.contig_id, ' '.join(ranges), sep='\t', file=ouf)

    def is_att_denovo(self, repeat: MultiLocation) -> bool:
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
            att_id = self.genome.features.get_features(FeatureType.ATT).get_next_att_id()
            for r in att.ranges:
                new_att = AttFeature(r, Strand.FORWARD, att.contig_id, self.genome, att_id)
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
            for att in atts_denovo:
                ranges = [f'({i.start},{i.end}' for i in att.ranges]
                print(att.contig_id, ' '.join(ranges), sep='\t', file=ouf)

    def _save_seqs_around_pipolbs(self, ranges_around_pipolbs: List[PairedLocation]):
        genome_dict = create_seqio_records_dict(file=self.genome.file, file_format='fasta')

        for i, range_pair in enumerate(ranges_around_pipolbs):
            left_seq = genome_dict[range_pair.contig_id][range_pair.left_range.start:range_pair.left_range.end]
            right_seq = genome_dict[range_pair.contig_id][range_pair.right_range.start:range_pair.right_range.end]
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
                    left_repeat = Range(start=hit.query_start, end=hit.query_end).shift(range_pair.left_range.start)
                    right_repeat = Range(start=hit.hit_start, end=hit.hit_end).shift(range_pair.right_range.start)
                    repeat_pair = PairedLocation(left_repeat, right_repeat, range_pair.contig_id)
                    repeats.append(repeat_pair)
        return repeats

    @staticmethod
    def _regroup_paired_repeats(paired_repeats: List[PairedLocation]) -> List[MultiLocation]:
        result: Mapping[ContigID, List[MultiLocation]] = defaultdict(list)
        for repeat in paired_repeats:
            for rs in result[repeat.contig_id]:
                if repeat.left_range.is_overlapping_any(rs.ranges) and \
                        not repeat.right_range.is_overlapping_any(rs.ranges):
                    rs.ranges.append(repeat.right_range)
                    break
                elif repeat.right_range.is_overlapping_any(rs.ranges) and \
                        not repeat.left_range.is_overlapping_any(rs.ranges):
                    rs.ranges.append(repeat.left_range)
                    break
                elif repeat.left_range.is_overlapping_any(rs.ranges) and \
                        repeat.right_range.is_overlapping_any(rs.ranges):
                    break
            else:
                result[repeat.contig_id].append(MultiLocation(
                    ranges=[repeat.right_range, repeat.left_range], contig_id=repeat.contig_id))
        return list(chain(*result.values()))


@task()
@genome_specific_logging
def are_atts_present(genome: Genome) -> Genome:
    logger = context.get('logger')

    num_atts = len(genome.features.get_features(FeatureType.ATT))
    if num_atts == 0:
        logger.warning('\n\n>>>No atts were found! Not able to define pipolin bounds!\n')
    return genome
