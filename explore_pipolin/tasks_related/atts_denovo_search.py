import os

from typing import List, Sequence

from Bio import SeqIO
from prefect import task

from explore_pipolin.common import Genome, FeatureType, Range, RangePair, Strand, AttFeature, Feature
from explore_pipolin.utilities.external_tools import blast_for_repeats
from explore_pipolin.utilities.io import read_blastxml
from explore_pipolin.tasks_related.misc import get_ranges_around_pipolbs
from explore_pipolin.utilities.logging import genome_specific_logging


@task()
@genome_specific_logging
def find_atts_denovo(genome: Genome, out_dir) -> Genome:
    atts_denovo_dir = os.path.join(out_dir, 'atts_denovo_search')
    os.makedirs(atts_denovo_dir, exist_ok=True)

    finder = AttsDenovoFinder(genome=genome, output_dir=atts_denovo_dir)

    repeats: List[RangePair] = finder.find_repeats()
    finder.write_repeats(repeats)

    atts_denovo: List[RangePair] = [rep for rep in repeats if finder.is_att_denovo(rep)]

    finder.add_atts_denovo_features(atts_denovo)
    finder.add_target_trnas_denovo_features()
    finder.write_atts_denovo()

    return genome


class AttsDenovoFinder:
    def __init__(self, genome: Genome, output_dir: str):
        self.genome = genome
        self.output_dir = output_dir

    def find_repeats(self) -> List[RangePair]:
        ranges_around_pipolbs = get_ranges_around_pipolbs(self.genome)
        self._save_seqs_around_pipolbs(ranges_around_pipolbs)
        blast_for_repeats(ranges_around_pipolbs, genome_id=self.genome.id, repeats_dir=self.output_dir)
        repeats = self._extract_repeats(ranges_around_pipolbs)
        return repeats

    def write_repeats(self, repeats: List[RangePair]):
        with open(os.path.join(self.output_dir, self.genome.id + '.repeats'), 'w') as ouf:
            print('left_rep_range', 'right_rep_range', 'length', 'contig_id', sep='\t', file=ouf)
            for repeat in repeats:
                repeat_length = max(repeat.left.end - repeat.left.start, repeat.right.end - repeat.right.start)
                print([repeat.left.start, repeat.right.end], (repeat.right.start, repeat.right.end),
                      repeat_length, repeat.contig_id, sep='\t', file=ouf)

    def is_att_denovo(self, repeat_pair: RangePair) -> bool:
        atts_of_contig = self.genome.features.atts_dict()[repeat_pair.contig_id]
        for att in atts_of_contig:
            if repeat_pair.left.is_overlapping(att):
                return False

        trnas_of_contig = self.genome.features.trnas_dict()[repeat_pair.contig_id]
        for trna in trnas_of_contig:
            if repeat_pair.left.is_overlapping(trna) or repeat_pair.right.is_overlapping(trna):
                return True
        return False

    def add_atts_denovo_features(self, atts_denovo: List[RangePair]):
        for att in atts_denovo:
            self.genome.features.add_features(AttFeature(location=att.left, strand=Strand.FORWARD, repeat_id='',
                                                         contig_id=self.genome.contigs[0].id, genome=self.genome),
                                              feature_type=FeatureType.ATT_DENOVO)
            self.genome.features.add_features(AttFeature(location=att.right, strand=Strand.FORWARD, repeat_id='',
                                                         contig_id=self.genome.contigs[0].id, genome=self.genome),
                                              feature_type=FeatureType.ATT_DENOVO)

    def add_target_trnas_denovo_features(self):
        for att in self.genome.features.get_features(FeatureType.ATT_DENOVO):
            target_trna = self.genome.features.get_features(FeatureType.TRNA).get_overlapping(att)
            if target_trna is not None:
                self.genome.features.add_features(target_trna, feature_type=FeatureType.TARGET_TRNA_DENOVO)

    def write_atts_denovo(self):
        atts_denovo = self.genome.features.get_features(FeatureType.ATT_DENOVO)
        with open(os.path.join(self.output_dir, self.genome.id + '.atts_denovo'), 'w') as ouf:
            print('att_start', 'att_end', sep='\t', file=ouf)
            for att in atts_denovo:
                print(att.start, att.end, sep='\t', file=ouf)

    def _save_seqs_around_pipolbs(self, ranges_around_pipolbs: List[RangePair]):
        genome_seq = SeqIO.read(handle=self.genome.file, format='fasta')

        for i, range_pair in enumerate(ranges_around_pipolbs):
            left_seq = genome_seq[range_pair.left.start:range_pair.left.end]
            right_seq = genome_seq[range_pair.right.start:range_pair.right.end]
            SeqIO.write(sequences=left_seq, format='fasta',
                        handle=os.path.join(self.output_dir, self.genome.id + f'_{i}.left'))
            SeqIO.write(sequences=right_seq, format='fasta',
                        handle=os.path.join(self.output_dir, self.genome.id + f'_{i}.right'))

    def _extract_repeats(self, ranges_around_pipolbs: List[RangePair]) -> List[RangePair]:
        repeats = []
        for i, range_pair in enumerate(ranges_around_pipolbs):
            repeats_xml = read_blastxml(os.path.join(self.output_dir, self.genome.id + f'_{i}.fmt5'))
            for entry in repeats_xml:
                for hit in entry:
                    left_repeat = Range(start=hit.query_start, end=hit.query_end).shift(range_pair.left.start)
                    right_repeat = Range(start=hit.hit_start, end=hit.hit_end).shift(range_pair.right.start)
                    repeat_pair = RangePair(left_repeat, right_repeat, range_pair.contig_id)
                    repeats.append(repeat_pair)
        return repeats
