from __future__ import annotations
import os

from typing import Sequence

from explore_pipolin.common import RepeatPair, Feature, Orientation, Genome, FeatureType
from explore_pipolin.utilities.external_tools import blast_for_repeats
from explore_pipolin.utilities.io import read_blastxml
from explore_pipolin.utilities.io import save_left_right_subsequences
from explore_pipolin.tasks_related.misc import get_left_right_windows


def is_att_denovo(genome: Genome, repeat_pair: RepeatPair) -> bool:
    if genome.features.is_overlapping_with(feature=repeat_pair.left, feature_type=FeatureType.ATT):
        return False
    left_overlaps = genome.features.is_overlapping_with(repeat_pair.left, FeatureType.TRNA)
    right_overlaps = genome.features.is_overlapping_with(repeat_pair.right, FeatureType.TRNA)
    return left_overlaps or right_overlaps


def find_repeats(genome: Genome, repeats_dir: str) -> Sequence[RepeatPair]:
    left_window, right_window = get_left_right_windows(genome, FeatureType.PIPOLB)
    save_left_right_subsequences(left_window, right_window, repeats_dir)
    blast_for_repeats(genome_id=genome.genome_id, repeats_dir=repeats_dir)
    repeats = _extract_repeats(file=os.path.join(repeats_dir, genome.genome_id + '.fmt5'), genome=genome)
    repeats = [rep.shift(left_window.start, right_window.start) for rep in repeats]
    return repeats


def _extract_repeats(file: str, genome: Genome) -> Sequence[RepeatPair]:
    repeats_xml = read_blastxml(file)
    repeats = []
    for entry in repeats_xml:
        for hit in entry:
            left_repeat = Feature(start=hit.query_start, end=hit.query_end, strand=Orientation.FORWARD,
                                  contig_id=genome.get_complete_genome_contig_id(), genome=genome)
            right_repeat = Feature(start=hit.hit_start, end=hit.hit_end, strand=Orientation.FORWARD,
                                   contig_id=genome.get_complete_genome_contig_id(), genome=genome)
            repeats.append(RepeatPair(left=left_repeat, right=right_repeat, seq=str(hit.hit.seq)))
    return repeats