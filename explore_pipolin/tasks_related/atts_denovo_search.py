import os

from typing import List, Iterable

from Bio import SeqIO

from explore_pipolin.common import Genome, FeatureType, Range, RangePair, Feature, Strand, AttFeature
from explore_pipolin.utilities.external_tools import blast_for_repeats
from explore_pipolin.utilities.io import read_blastxml
from explore_pipolin.tasks_related.misc import get_windows_around_pipolbs


class AttsDenovoFinder:
    def __init__(self, genome: Genome, output_dir: str):
        self.genome = genome
        self.output_dir = output_dir

    def find_atts_denovo(self):
        repeats: List[RangePair] = self._find_repeats()
        self._write_repeats(repeats=repeats)

        atts_denovo: List[RangePair] = [rep for rep in repeats if self._is_att_denovo(rep)]

        self._add_atts_denovo_features(atts_denovo)
        self._add_target_trnas_denovo_features()
        self._write_atts_denovo(self.genome.features.get_features(FeatureType.ATT_DENOVO))

    def _find_repeats(self) -> List[RangePair]:
        windows = get_windows_around_pipolbs(self.genome)
        self._save_left_right_subsequences(windows)
        blast_for_repeats(windows, genome_id=self.genome.id, repeats_dir=self.output_dir)
        repeats = self._extract_repeats(windows)
        return repeats

    def _write_repeats(self, repeats: List[RangePair]):
        with open(os.path.join(self.output_dir, self.genome.id + '.repeats'), 'w') as ouf:
            print('left_rep_range', 'right_rep_range', 'length', 'contig_id', sep='\t', file=ouf)
            for repeat in repeats:
                print([repeat.left.start, repeat.right.end], (repeat.right.start, repeat.right.end),
                      repeat.left.end - repeat.left.start, repeat.contig_id, sep='\t', file=ouf)

    def _is_att_denovo(self, repeat_pair: RangePair) -> bool:
        if self.genome.features.is_overlapping_with(repeat_pair.left, FeatureType.ATT):
            return False
        left_overlaps = self.genome.features.is_overlapping_with(repeat_pair.left, FeatureType.TRNA)
        right_overlaps = self.genome.features.is_overlapping_with(repeat_pair.right, FeatureType.TRNA)
        return left_overlaps or right_overlaps

    def _add_atts_denovo_features(self, atts_denovo: List[RangePair]):
        for att in atts_denovo:
            self.genome.features.add_feature(feature=AttFeature(
                location=att[0],
                strand=Strand.FORWARD,
                repeat_id='',
                contig_id=self.genome.contigs[0].id,
                genome=self.genome), feature_type=FeatureType.ATT_DENOVO)
            self.genome.features.add_feature(feature=AttFeature(
                location=att[1],
                strand=Strand.FORWARD,
                repeat_id='',
                contig_id=self.genome.contigs[0].id,
                genome=self.genome), feature_type=FeatureType.ATT_DENOVO)

    def _add_target_trnas_denovo_features(self):
        for att in self.genome.features.get_features(FeatureType.ATT_DENOVO):
            target_trna = self.genome.features.get_features(FeatureType.TRNA).get_overlapping(att)
            if target_trna is not None:
                self.genome.features.add_feature(feature=target_trna, feature_type=FeatureType.TARGET_TRNA_DENOVO)

    def _write_atts_denovo(self, atts_denovo: Iterable[Feature]):
        with open(os.path.join(self.output_dir, self.genome.id + '.atts_denovo'), 'w') as ouf:
            print('att_start', 'att_end', sep='\t', file=ouf)
            for att in atts_denovo:
                print(att.start, att.end, sep='\t', file=ouf)

    def _save_left_right_subsequences(self, windows: List[RangePair]):
        genome_seq = SeqIO.read(handle=self.genome.file, format='fasta')

        for i, window in enumerate(windows):
            left_seq = genome_seq[window.left.start:window.left.end]
            right_seq = genome_seq[window.right.start:window.right.end]
            SeqIO.write(sequences=left_seq, format='fasta',
                        handle=os.path.join(self.output_dir, self.genome.id + f'_{i}.left'))
            SeqIO.write(sequences=right_seq, format='fasta',
                        handle=os.path.join(self.output_dir, self.genome.id + f'_{i}.right'))

    def _extract_repeats(self, windows: List[RangePair]) -> List[RangePair]:
        repeats = []
        for i, window in enumerate(windows):
            repeats_xml = read_blastxml(os.path.join(self.output_dir, self.genome.id + f'_{i}.fmt5'))
            for entry in repeats_xml:
                for hit in entry:
                    left_repeat = Range(start=hit.query_start, end=hit.query_end).shift(window.left.start)
                    right_repeat = Range(start=hit.hit_start, end=hit.hit_end).shift(window.right.start)
                    repeat_pair = RangePair(left_repeat, right_repeat, window.contig_id)
                    repeats.append(repeat_pair)
        return repeats
