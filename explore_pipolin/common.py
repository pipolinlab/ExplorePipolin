from __future__ import annotations

import os
from enum import Enum, auto
from typing import Tuple, MutableSequence, Optional


class Orientation(Enum):
    FORWARD = auto()
    REVERSE = auto()

    @staticmethod
    def orientation_from_blast(hit_strand):
        if hit_strand == 1:
            return Orientation.FORWARD
        elif hit_strand == -1:
            return Orientation.REVERSE
        else:
            raise AssertionError(f'Unknown hit_frame: {hit_strand}! 0 is also unexpected in our case!')

    def to_pm_one_encoding(self):
        return +1 if self is self.FORWARD else -1

    def __neg__(self):
        if self is self.FORWARD:
            return self.REVERSE
        else:
            return self.FORWARD


class Contig:
    def __init__(self, contig_id: str, contig_length: int, orientation=Orientation.FORWARD):
        self.contig_id = contig_id
        self.contig_length = contig_length
        self.contig_orientation: Orientation = orientation


class Genome:
    def __init__(self, genome_id: str, genome_file: str):
        self.genome_id = genome_id
        self.genome_file = genome_file
        self.contigs: MutableSequence[Contig] = []

    def get_contig_by_id(self, contig_id: str) -> Optional[Contig]:
        for contig in self.contigs:
            if contig.contig_id == contig_id:
                return contig
        raise AssertionError('It is unexpected that you ask for a non-existent contig!')

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
    def __init__(self, start: int, end: int, strand: Orientation, contig_id: str, genome: Genome):
        self.start = start
        self.end = end
        self.strand = strand
        self.contig_id = contig_id
        self.genome = genome

        if self.start > self.end:
            raise AssertionError('Feature start cannot be greater than feature end!')

        if self.end > self.contig.contig_length:
            raise AssertionError('Feature end cannot be greater than contig length!')

    @property
    def contig(self):
        return self.genome.get_contig_by_id(self.contig_id)


class FeatureType(Enum):
    PIPOLB = auto()
    ATT = auto()
    TRNA = auto()
    TARGET_TRNA = auto()


class Repeat:
    def __init__(self, left: Tuple[int, int], right: Tuple[int, int], seq: str):
        self.left = left
        self.right = right
        self.seq = seq

        if self.left[0] > self.left[1] or self.right[0] > self.right[1]:
            raise AssertionError('Repeat start cannot be greater than repeat end!')

        left_repeat_length = self.left[1] - self.left[0]
        right_repeat_length = self.right[1] - self.right[0]
        seq_length = len(self.seq)
        if seq_length > left_repeat_length or seq_length > right_repeat_length:
            raise AssertionError('Repeat sequence length cannot be greater than repeat ranges!')

    def shift(self, left_shift: int, right_shift: int) -> Repeat:
        if left_shift > right_shift:
            raise AssertionError('Left shift cannot be greater than right shift!')

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

        if self.start > self.end:
            raise AssertionError('Fragment start cannot be greater than fragment end!')

        if self.end > self.contig.contig_length:
            raise AssertionError('Fragment end cannot be greater than contig length!')

    @property
    def contig(self):
        return self.genome.get_contig_by_id(self.contig_id)


def define_genome_id(genome_path: str):
    genome_id = os.path.splitext(os.path.basename(genome_path))[0]
    check_genome_id_length(genome_id)
    return genome_id


def check_genome_id_length(genome_id):
    if len(genome_id) > 16:
        raise AssertionError('Genome file basename is going to be used as a genome identifier. '
                             'Due to Biopython restrictions, it cannot be longer than 16 characters. '
                             'Please, rename the file, so that its basename does not exceed the limit!')
