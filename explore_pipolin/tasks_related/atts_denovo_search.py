import os

from typing import Sequence, List

from explore_pipolin.common import RepeatPair, Genome, FeatureType, Range
from explore_pipolin.utilities.external_tools import blast_for_repeats
from explore_pipolin.utilities.io import read_blastxml
from explore_pipolin.utilities.io import save_left_right_subsequences
from explore_pipolin.tasks_related.misc import get_windows, Window


def is_att_denovo(genome: Genome, repeat_pair: RepeatPair) -> bool:
    if genome.features.is_overlapping_with(repeat_pair.left, FeatureType.ATT):
        return False
    left_overlaps = genome.features.is_overlapping_with(repeat_pair.left, FeatureType.TRNA)
    right_overlaps = genome.features.is_overlapping_with(repeat_pair.right, FeatureType.TRNA)
    return left_overlaps or right_overlaps


def find_repeats(genome: Genome, repeats_dir: str) -> Sequence[RepeatPair]:
    windows = get_windows(genome)
    save_left_right_subsequences(windows, repeats_dir)
    blast_for_repeats(windows, repeats_dir)
    repeats = _extract_repeats(windows, repeats_dir)
    return repeats


def _extract_repeats(windows: List[Window], repeats_dir: str) -> Sequence[RepeatPair]:
    genome = windows[0].pipolbs[0].genome
    repeats = []
    for i, window in enumerate(windows):
        repeats_xml = read_blastxml(os.path.join(repeats_dir, genome.genome_id + f'_{i}.fmt5'))
        for entry in repeats_xml:
            for hit in entry:
                left_repeat = Range(start=hit.query_start, end=hit.query_end)
                right_repeat = Range(start=hit.hit_start, end=hit.hit_end)
                repeat_pair = RepeatPair(left=left_repeat, right=right_repeat, left_seq=str(hit.query.seq.ungap("-")),
                                         right_seq=str(hit.hit.seq.ungap("-")), pipolbs=window.pipolbs)
                repeat_pair.shift(left_shift=window.left.start, right_shift=window.right.start)
                repeats.append(repeat_pair)
    return repeats
