import os

from explore_pipolin.utilities.external_tools_run import blast_for_repeats
from explore_pipolin.utilities.io import read_blastxml
from explore_pipolin.utilities.io import save_left_right_subsequences


def find_repeats(genome, gquery, repeats_dir):
    left_window, right_window = gquery.get_left_right_windows()
    save_left_right_subsequences(genome=genome, left_window=left_window, right_window=right_window,
                                 repeats_dir=repeats_dir)
    blast_for_repeats(gquery_id=gquery.gquery_id, repeats_dir=repeats_dir)
    left_repeats, right_repeats, sequences = _extract_repeats(file=os.path.join(repeats_dir, gquery.gquery_id + '.fmt5'))
    left_repeats = [_set_proper_location(seq_range=i, shift=left_window[0]) for i in left_repeats]
    right_repeats = [_set_proper_location(seq_range=i, shift=right_window[0]) for i in right_repeats]
    return left_repeats, right_repeats, sequences


def _extract_repeats(file):
    repeats = read_blastxml(file)
    left_repeats = []
    right_repeats = []
    sequences = []
    for entry in repeats:
        for hit in entry:
            left_repeats.append((hit.query_start, hit.query_end))
            right_repeats.append((hit.hit_start, hit.hit_end))
            sequences.append(str(hit.hit.seq))
    return left_repeats, right_repeats, sequences


def _set_proper_location(seq_range, shift):
    return seq_range[0] + shift, seq_range[1] + shift
