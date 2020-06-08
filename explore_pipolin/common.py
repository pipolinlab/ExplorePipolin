from __future__ import annotations

import os
from enum import Enum, auto
from typing import Tuple, MutableSequence, Optional


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

    def to_string(self):
        return 'forward' if self.FORWARD else 'reverse'

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


class Genome:
    def __init__(self, genome_id: str):
        self.genome_id = genome_id
        self.contigs: MutableSequence[Contig] = []

    def get_contig_by_id(self, contig_id: str) -> Optional[Contig]:
        for contig in self.contigs:
            if contig.contig_id == contig_id:
                return contig
        return None

    def is_single_contig(self):
        return len(self.contigs) == 1

    def get_complete_genome_contig_id(self):
        if not self.is_single_contig():
            raise AssertionError('Unsupported! Not a complete genome!')
        return self.contigs[0].contig_id

    def get_complete_genome_length(self) -> int:
        if not self.is_single_contig():
            raise AssertionError('Unsupported! Not a complete genome!')
        return self.contigs[0].contig_length


class Feature:
    def __init__(self, start: int, end: int, frame: Orientation, contig_id: str, genome: Genome):
        self.start = start
        self.end = end
        self.frame = frame
        self.contig_id = contig_id
        self.genome = genome

    @property
    def contig(self):
        return self.genome.get_contig_by_id(self.contig_id)


class FeatureType(Enum):
    PIPOLB = auto()
    ATT = auto()
    TARGET_TRNA = auto()
    TRNA = auto()


class Repeat:
    def __init__(self, left: Tuple[int, int], right: Tuple[int, int], seq: str):
        self.left = left
        self.right = right
        self.seq = seq

    def shift(self, left_shift: int, right_shift: int) -> Repeat:
        return Repeat(self._shift_range(self.left, left_shift), self._shift_range(self.right, right_shift), self.seq)

    @staticmethod
    def _shift_range(seq_range: Tuple[int, int], shift):
        return seq_range[0] + shift, seq_range[1] + shift


class PipolinFragment:
    def __init__(self, contig_id: str, genome: Genome, start: int, end: int):
        self.contig_id = contig_id
        self.start = start
        self.end = end
        self.atts: MutableSequence[Feature] = []
        self.genome = genome

    @property
    def contig(self):
        return self.genome.get_contig_by_id(self.contig_id)


def define_genome_id(genome):
    return os.path.splitext(os.path.basename(genome))[0]
