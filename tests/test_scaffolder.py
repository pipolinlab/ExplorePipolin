import unittest
from typing import Sequence

from explore_pipolin.common import Genome, Contig, Feature, Range, Strand, FeatureType
from explore_pipolin.tasks_related.scaffolding import Scaffolder


_GENOME_ID = 'GENOME'
_GENOME_FILE = 'GENOME.fa'

_ATT = 'att'
_PIPOLB = 'pol'
_TRNA = '(t)'   # goes always after att!
# intragenic = '---'
_GAP = '...'


def _create_genome_from_scheme(scheme: str) -> Genome:
    if len(scheme) % 3 != 0:
        raise AssertionError('Wrong scheme! Its length should be multiple of 3!')

    contigs_schemes = scheme.split(_GAP)
    genome = _create_genome_with_contigs(contigs_schemes)

    _add_features_to_genome(contigs_schemes, genome)

    return genome


def _create_genome_with_contigs(contigs_schemes: Sequence[str]) -> Genome:
    contigs = []
    for i in range(len(contigs_schemes)):
        contigs.append(Contig(f'CONTIG_{i}', len(contigs_schemes[i]) // 3 * 10))
    genome = Genome(_GENOME_ID, _GENOME_FILE, contigs)
    return genome


def _add_features_to_genome(contigs_schemes, genome):
    for i, sch in enumerate(contigs_schemes):
        for triplet_start in range(0, len(sch), 3):
            triplet = sch[triplet_start: triplet_start + 3]
            feature_start = (triplet_start // 3 * 10 - 5) if triplet == _TRNA else (triplet_start // 3 * 10)
            feature = Feature(Range(feature_start, feature_start + 10), Strand.FORWARD, f'CONTIG_{i}', genome)
            if triplet == _ATT:
                genome.features.add_feature(feature, FeatureType.ATT)
            if triplet == _PIPOLB:
                genome.features.add_feature(feature, FeatureType.PIPOLB)
            if triplet == _TRNA:
                genome.features.add_feature(feature, FeatureType.TARGET_TRNA)


class SetUpPipolinsToScaffold(unittest.TestCase):
    def setUp(self) -> None:
        self.genome1 = _create_genome_from_scheme('---pol---')
        self.genome2 = _create_genome_from_scheme('---att---pol---')   # ---att(t)---pol---
        self.genome3 = _create_genome_from_scheme('---att---pol---att---')   # ---att(t)---polb---att---
        self.genome4 = _create_genome_from_scheme('---att---pol---...---att---')
        self.genome5 = _create_genome_from_scheme('---att---...---pol---...---att---')
        self.genome5_trna = _create_genome_from_scheme('---att---...---pol---...---att(t)---')

        self.genome6 = _create_genome_from_scheme('---att---pol---att---att(t)---')
        self.genome7 = _create_genome_from_scheme('---att---...---pol---att---att(t)---')
        # ---att---pol---...---att---att(t)---
        # ---att---pol---att---...---att(t)---
        # ---att---...---pol---...---att---att(t)---
        # ---att---pol---...---att---...---att(t)---
        # ---att---...---pol---att---...---att(t)---


class TestScaffolder(SetUpPipolinsToScaffold):
    def test_genome1(self):
        with self.assertRaises(NotImplementedError):
            Scaffolder(self.genome1).scaffold()

    def test_genome2(self):
        with self.assertRaises(NotImplementedError):
            Scaffolder(self.genome2).scaffold()

    def test_genome3(self):
        with self.assertRaises(NotImplementedError):
            Scaffolder(self.genome3).scaffold()

    def test_genome4(self):
        Scaffolder(self.genome4).scaffold()

    def test_genome5(self):
        with self.assertRaises(NotImplementedError):
            Scaffolder(self.genome5).scaffold()
        Scaffolder(self.genome5_trna).scaffold()

    def test_genome6(self):
        with self.assertRaises(NotImplementedError):
            Scaffolder(self.genome6).scaffold()

    def test_genome7(self):
        Scaffolder(self.genome7).scaffold()
