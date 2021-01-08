import unittest
from typing import Sequence

from explore_pipolin.common import Genome, Contig, Feature, Range, Strand, FeatureType, Pipolin, PipolinFragment, \
    AttFeature
from explore_pipolin.tasks_related.scaffolding import scaffold


_GENOME_ID = 'GENOME'
_GENOME_FILE = 'GENOME.fa'

_ATT = 'att'
_PIPOLB = 'pol'
_TRNA = '(t)'   # goes always after att!
# intragenic = '---'
_GAP = '...'


def create_genome_from_scheme(scheme: str) -> Genome:
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
                feature = AttFeature(Range(feature_start, feature_start + 10), Strand.FORWARD,
                                     f'CONTIG_{i}', genome, repeat_id=0)
                genome.features.add_features(feature, feature_type=FeatureType.ATT)
            if triplet == _PIPOLB:
                genome.features.add_features(feature, feature_type=FeatureType.PIPOLB)
            if triplet == _TRNA:
                genome.features.add_features(feature, feature_type=FeatureType.TARGET_TRNA)


class TestScaffolder(unittest.TestCase):
    def _check_pipolins(self, expected: Sequence[Pipolin], obtained: Sequence[Pipolin]):
        self.assertEqual(len(expected), len(obtained))
        # TODO: how to check what's inside each Pipolin?

    def test_genome1(self):
        genome = create_genome_from_scheme('---pol---')
        exp = [Pipolin(PipolinFragment(Range(0, 30), 'CONTIG_0', genome))]

        self._check_pipolins(exp, scaffold(genome))

    def test_genome2(self):
        genome = create_genome_from_scheme('---att---pol---')
        exp = [Pipolin(PipolinFragment(Range(0, 50), 'CONTIG_0', genome))]

        self._check_pipolins(exp, scaffold(genome))

    def test_genome3(self):
        genome = create_genome_from_scheme('---att---pol---att---')
        exp = [Pipolin(PipolinFragment(Range(0, 70), 'CONTIG_0', genome))]

        self._check_pipolins(exp, scaffold(genome))

    def test_genome4(self):
        genome = create_genome_from_scheme('---att---pol---...---att---')
        f1 = PipolinFragment(Range(0, 50), 'CONTIG_0', genome)
        f2 = PipolinFragment(Range(0, 30), 'CONTIG_1', genome)
        exp = [Pipolin(f1, f2)]

        self._check_pipolins(exp, scaffold(genome))

    def test_genome5(self):
        genome = create_genome_from_scheme('---att---...---pol---...---att---')
        f1 = PipolinFragment(Range(0, 30), 'CONTIG_0', genome)
        f2 = PipolinFragment(Range(0, 30), 'CONTIG_1', genome)
        f3 = PipolinFragment(Range(0, 30), 'CONTIG_2', genome)
        exp = [Pipolin(f1, f2, f3)]

        self._check_pipolins(exp, scaffold(genome))

    def test_genome5_trna(self):
        genome = create_genome_from_scheme('---att---...---pol---...---att(t)---')
        f1 = PipolinFragment(Range(0, 30), 'CONTIG_0', genome)
        f2 = PipolinFragment(Range(0, 30), 'CONTIG_1', genome)
        f3 = PipolinFragment(Range(0, 30), 'CONTIG_2', genome)
        exp = [Pipolin(f1, f2, f3)]

        self._check_pipolins(exp, scaffold(genome))

    def test_genome6(self):
        genome = create_genome_from_scheme('---att---pol---att---att(t)---')
        exp = [Pipolin(PipolinFragment(Range(0, 90), 'CONTIG_0', genome))]
        self._check_pipolins(exp, scaffold(genome))

    def test_genome7(self):
        genome = create_genome_from_scheme('---att---...---pol---att---att(t)---')
        f1 = PipolinFragment(Range(0, 30), 'CONTIG_0', genome)
        f2 = PipolinFragment(Range(0, 70), 'CONTIG_1', genome)
        exp = [Pipolin(f1, f2)]

        self._check_pipolins(exp, scaffold(genome))

    def test_genome8(self):
        genome = create_genome_from_scheme('---att---pol---...---att---att(t)---')
        f1 = PipolinFragment(Range(0, 50), 'CONTIG_0', genome)
        f2 = PipolinFragment(Range(0, 50), 'CONTIG_1', genome)
        exp = [Pipolin(f1, f2)]

        self._check_pipolins(exp, scaffold(genome))

    def test_genome9(self):
        genome = create_genome_from_scheme('---att---pol---att---...---att(t)---')
        f1 = PipolinFragment(Range(0, 70), 'CONTIG_0', genome)
        f2 = PipolinFragment(Range(0, 30), 'CONTIG_1', genome)
        exp = [Pipolin(f1, f2)]

        self._check_pipolins(exp, scaffold(genome))

    def test_genome10(self):
        genome = create_genome_from_scheme('---att---...---pol---...---att---att(t)---')
        f1 = PipolinFragment(Range(0, 30), 'CONTIG_0', genome)
        f2 = PipolinFragment(Range(0, 30), 'CONTIG_1', genome)
        f3 = PipolinFragment(Range(0, 50), 'CONTIG_2', genome)
        exp = [Pipolin(f1, f2, f3)]

        self._check_pipolins(exp, scaffold(genome))

    def test_genome11(self):
        genome = create_genome_from_scheme('---att---pol---...---att---...---att(t)---')
        f1 = PipolinFragment(Range(0, 50), 'CONTIG_0', genome)
        f2 = PipolinFragment(Range(0, 30), 'CONTIG_1', genome)
        f3 = PipolinFragment(Range(0, 30), 'CONTIG_2', genome)
        exp = [Pipolin(f1, f2, f3)]

        self._check_pipolins(exp, scaffold(genome))

    def test_genome12(self):
        genome = create_genome_from_scheme('---att---...---pol---att---...---att(t)---')
        f1 = PipolinFragment(Range(0, 30), 'CONTIG_0', genome)
        f2 = PipolinFragment(Range(0, 50), 'CONTIG_1', genome)
        f3 = PipolinFragment(Range(0, 30), 'CONTIG_2', genome)
        exp = [Pipolin(f1, f2, f3)]

        self._check_pipolins(exp, scaffold(genome))

    def test_genome13(self):   # consider as two pipolins
        genome = create_genome_from_scheme('---att---pol---att---pol---att(t)---')
        p1 = Pipolin(PipolinFragment(Range(0, 70), 'CONTIG_0', genome))
        p2 = Pipolin(PipolinFragment(Range(40, 110), 'CONTIG_0', genome))

        self._check_pipolins([p1, p2], scaffold(genome))

    def test_genome14(self):   # consider as two pipolins
        genome = create_genome_from_scheme('---att---...---pol---att---pol---att(t)---')
        p1f1 = PipolinFragment(Range(0, 30), 'CONTIG_0', genome)
        p1f2 = PipolinFragment(Range(0, 50), 'CONTIG_1', genome)
        p2 = Pipolin(PipolinFragment(Range(20, 90), 'CONTIG_1', genome))

        self._check_pipolins([Pipolin(p1f1, p1f2), p2], scaffold(genome))

    def test_genome15(self):   # consider as two pipolins
        genome = create_genome_from_scheme('---att---pol---...---att---pol---att(t)---')
        p1f1 = PipolinFragment(Range(0, 50), 'CONTIG_0', genome)
        p1f2 = PipolinFragment(Range(0, 30), 'CONTIG_1', genome)
        p2 = Pipolin(PipolinFragment(Range(0, 70), 'CONTIG_1', genome))

        self._check_pipolins([Pipolin(p1f1, p1f2), p2], scaffold(genome))

    def test_genome16(self):   # consider as two pipolins
        genome = create_genome_from_scheme('---att---pol---att---...---pol---att(t)---')
        p1 = Pipolin(PipolinFragment(Range(0, 70), 'CONTIG_0', genome))
        p2f1 = PipolinFragment(Range(40, 70), 'CONTIG_0', genome)
        p2f2 = PipolinFragment(Range(0, 50), 'CONTIG_1', genome)

        self._check_pipolins([p1, Pipolin(p2f1, p2f2)], scaffold(genome))

    def test_genome18(self):   # two pipolins
        genome = create_genome_from_scheme('---att---pol---att(t)---att---pol---att(t)---')
        p1 = Pipolin(PipolinFragment(Range(0, 70), 'CONTIG_0', genome))
        p2 = Pipolin(PipolinFragment(Range(60, 130), 'CONTIG_0', genome))

        self._check_pipolins([p1, p2], scaffold(genome))

    def test_genome19(self):   # two pipolins
        genome = create_genome_from_scheme('---att---pol---att(t)---...---att---pol---att(t)---')
        p1 = Pipolin(PipolinFragment(Range(0, 70), 'CONTIG_0', genome))
        p2 = Pipolin(PipolinFragment(Range(0, 70), 'CONTIG_1', genome))

        self._check_pipolins([p1, p2], scaffold(genome))

    def test_genome20(self):   # consider as two pipolins
        genome = create_genome_from_scheme('---att---pol---att---pol---')
        p1 = Pipolin(PipolinFragment(Range(0, 70), 'CONTIG_0', genome))
        p2 = Pipolin(PipolinFragment(Range(60, 90), 'CONTIG_0', genome))

        self._check_pipolins([p1, p2], scaffold(genome))

    def test_genome21(self):
        genome = create_genome_from_scheme('---att---pol---pol---att---')
        exp = [Pipolin(PipolinFragment(Range(0, 90), 'CONTIG_0', genome))]

        self._check_pipolins(exp, scaffold(genome))

    def test_genome22(self):
        genome = create_genome_from_scheme('---att---att---pol---pol---')
        exp = [Pipolin(PipolinFragment(Range(0, 90), 'CONTIG_0', genome))]

        self._check_pipolins(exp, scaffold(genome))

    def test_genome23(self):
        genome = create_genome_from_scheme('---att---pol---pol---att---att---')
        exp = [Pipolin(PipolinFragment(Range(0, 110), 'CONTIG_0', genome))]

        self._check_pipolins(exp, scaffold(genome))
