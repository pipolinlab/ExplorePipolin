from __future__ import annotations
import os

from typing import Sequence

from explore_pipolin.utilities.common import Repeat
from explore_pipolin.utilities.external_tools import blast_for_repeats
from explore_pipolin.utilities.io import read_blastxml
from explore_pipolin.utilities.io import save_left_right_subsequences


def find_repeats(genome, gquery, repeats_dir) -> Sequence[Repeat]:
    left_window, right_window = gquery.get_left_right_windows()
    save_left_right_subsequences(genome=genome, left_window=left_window, right_window=right_window,
                                 repeats_dir=repeats_dir)
    blast_for_repeats(gquery_id=gquery.gquery_id, repeats_dir=repeats_dir)
    repeats = _extract_repeats(file=os.path.join(repeats_dir, gquery.gquery_id + '.fmt5'))
    repeats = [rep.shift(left_window[0], right_window[0]) for rep in repeats]
    return repeats


def _extract_repeats(file) -> Sequence[Repeat]:
    repeats_xml = read_blastxml(file)
    repeats = []
    for entry in repeats_xml:
        for hit in entry:
            repeats.append(Repeat((hit.query_start, hit.query_end), (hit.hit_start, hit.hit_end), str(hit.hit.seq)))
    return repeats
