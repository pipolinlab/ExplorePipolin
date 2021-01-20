import unittest
from typing import Sequence, Tuple, List

from explore_pipolin.common import Genome, Contig, Feature, Range, Strand, FeatureType, \
    Pipolin, PipolinFragment, AttFeature, ContigID
from explore_pipolin.tasks_related.find_pipolins import find_pipolins

_GENOME_ID = 'GENOME'
_GENOME_FILE = 'GENOME.fa'

_ATT = 'at'
_PIPOLB = 'pol'
_TRNA = '(t)'   # goes always after att!
# intragenic = '---'
_GAP = '...'


def create_genome_from_scheme(scheme: str) -> Genome:
    contigs_schemes = scheme.split(_GAP)
    genome = _create_genome_with_contigs(contigs_schemes)

    features = _create_features_for_genome(genome, contigs_schemes)
    _add_features_to_genome(features, genome)

    return genome


def _create_genome_with_contigs(contigs_schemes: Sequence[str]) -> Genome:
    contigs = []
    for i in range(len(contigs_schemes)):
        contigs.append(Contig(ContigID(f'CONTIG_{i}'), len(contigs_schemes[i]) // 3 * 100))
    genome = Genome(_GENOME_ID, _GENOME_FILE, contigs)
    return genome


def _create_features_for_genome(
        genome: Genome, contigs_schemes: List[str]) -> Sequence[Tuple[Feature, FeatureType]]:

    features = []
    for i, sch in enumerate(contigs_schemes):
        for triplet_start in range(0, len(sch), 3):
            triplet = sch[triplet_start: triplet_start + 3]
            feature_start = (triplet_start // 3 * 100 - 5) if triplet == _TRNA else (triplet_start // 3 * 100)
            if triplet[:2] == _ATT:
                feature = AttFeature(Range(feature_start, feature_start + 100), Strand.FORWARD,
                                     ContigID(f'CONTIG_{i}'), genome, att_id=int(triplet[2]))
                features.append((feature, FeatureType.ATT))
            else:
                feature = Feature(Range(feature_start, feature_start + 100),
                                  Strand.FORWARD, ContigID(f'CONTIG_{i}'), genome)
                if triplet == _PIPOLB:
                    features.append((feature, FeatureType.PIPOLB))
                if triplet == _TRNA:
                    features.append((feature, FeatureType.TARGET_TRNA))

    return features


def _add_features_to_genome(
        features: Sequence[Tuple[Feature, FeatureType]], genome: Genome) -> None:
    for feature in features:
        genome.features.add_features(feature[0], feature_type=feature[1])


def create_pipolin(genome: Genome, scheme: str, *fragments: PipolinFragment) -> Pipolin:
    new_fragments = []
    for fragment in fragments:
        features = _create_features_for_pipolin_fragment(genome, scheme, fragment)
        new_fragments.append(PipolinFragment(fragment.location, fragment.contig_id, features))
    return Pipolin.from_fragments(*new_fragments)


def _create_features_for_pipolin_fragment(
        genome: Genome, scheme: str, fragment: PipolinFragment) -> Sequence[Tuple[Feature, FeatureType]]:

    contigs_schemes = scheme.split(_GAP)
    features = _create_features_for_genome(genome, contigs_schemes)

    fragment_features = []
    for feature in features:
        if feature[0].contig_id == fragment.contig_id:
            if feature[0].start >= fragment.start and feature[0].end <= fragment.end:
                fragment_features.append(feature)

    return tuple(fragment_features)


class TestPipolinFinder(unittest.TestCase):
    def _check_pipolins(self, expected: Sequence[Pipolin], obtained: Sequence[Pipolin]):
        self.assertEqual(len(expected), len(obtained))
        expected = self._order_pipolin_fragments_for_each_pipolin(expected)
        obtained = self._order_pipolin_fragments_for_each_pipolin(obtained)
        self.assertEqual(set(expected), set(obtained))

    @staticmethod
    def _order_pipolin_fragments_for_each_pipolin(pipolins: Sequence[Pipolin]) -> Sequence[Pipolin]:
        new_pipolins = []
        for pipolin in pipolins:
            sorted_pipolin = Pipolin.from_fragments(*sorted(pipolin.fragments, key=lambda x: x.contig_id))
            new_pipolins.append(sorted_pipolin)
        return new_pipolins

    def test_genome1(self):
        scheme = '---pol---'
        genome = create_genome_from_scheme(scheme)
        f1 = PipolinFragment(Range(100, 200), ContigID('CONTIG_0'))
        exp = [create_pipolin(genome, scheme, f1)]

        self._check_pipolins(exp, find_pipolins.run(genome))

    def test_genome1_1(self):
        scheme = '---pol---pol---'
        genome = create_genome_from_scheme(scheme)
        p1f1 = PipolinFragment(Range(100, 200), ContigID('CONTIG_0'))
        p2f1 = PipolinFragment(Range(300, 400), ContigID('CONTIG_0'))
        exp = [create_pipolin(genome, scheme, p1f1), create_pipolin(genome, scheme, p2f1)]

        self._check_pipolins(exp, find_pipolins.run(genome))

    def test_genome1_2(self):
        pass   # genome = create_genome_from_scheme('---pol===pol---')

    def test_genome2(self):
        scheme = '---at1---pol---'
        genome = create_genome_from_scheme(scheme)
        f1 = PipolinFragment(Range(100, 400), ContigID('CONTIG_0'))
        exp = [create_pipolin(genome, scheme, f1)]

        self._check_pipolins(exp, find_pipolins.run(genome))

    def test_genome3(self):
        scheme = '---at1---pol---at1---'
        genome = create_genome_from_scheme(scheme)
        f1 = PipolinFragment(Range(100, 600), ContigID('CONTIG_0'))
        exp = [create_pipolin(genome, scheme, f1)]

        self._check_pipolins(exp, find_pipolins.run(genome))

    def test_genome4(self):
        scheme = '---at1---pol---...---at1---'
        genome = create_genome_from_scheme(scheme)
        f1 = PipolinFragment(Range(100, 400), ContigID('CONTIG_0'))
        f2 = PipolinFragment(Range(100, 200), ContigID('CONTIG_1'))
        exp = [create_pipolin(genome, scheme, f1, f2)]

        self._check_pipolins(exp, find_pipolins.run(genome))

    def test_genome5(self):
        scheme = '---at1---...---pol---...---at1---'
        genome = create_genome_from_scheme(scheme)
        f1 = PipolinFragment(Range(100, 200), ContigID('CONTIG_0'))
        f2 = PipolinFragment(Range(100, 200), ContigID('CONTIG_1'))
        f3 = PipolinFragment(Range(100, 200), ContigID('CONTIG_2'))
        exp = [create_pipolin(genome, scheme, f1, f2, f3)]

        self._check_pipolins(exp, find_pipolins.run(genome))

    def test_genome5_trna(self):
        scheme = '---at1---...---pol---...---at1(t)---'
        genome = create_genome_from_scheme(scheme)
        f1 = PipolinFragment(Range(100, 200), ContigID('CONTIG_0'))
        f2 = PipolinFragment(Range(100, 200), ContigID('CONTIG_1'))
        f3 = PipolinFragment(Range(100, 200), ContigID('CONTIG_2'))
        exp = [create_pipolin(genome, scheme, f1, f2, f3)]

        self._check_pipolins(exp, find_pipolins.run(genome))

    def test_genome6(self):
        scheme = '---at1---pol---at1---at1(t)---'
        genome = create_genome_from_scheme(scheme)
        f1 = PipolinFragment(Range(100, 800), ContigID('CONTIG_0'))
        exp = [create_pipolin(genome, scheme, f1)]

        self._check_pipolins(exp, find_pipolins.run(genome))

    def test_genome7(self):
        scheme = '---at1---...---pol---at1---at1(t)---'
        genome = create_genome_from_scheme(scheme)
        f1 = PipolinFragment(Range(100, 200), ContigID('CONTIG_0'))
        f2 = PipolinFragment(Range(100, 600), ContigID('CONTIG_1'))
        exp = [create_pipolin(genome, scheme, f1, f2)]

        self._check_pipolins(exp, find_pipolins.run(genome))

    def test_genome8(self):
        scheme = '---at1---pol---...---at1---at1(t)---'
        genome = create_genome_from_scheme(scheme)
        f1 = PipolinFragment(Range(100, 400), ContigID('CONTIG_0'))
        f2 = PipolinFragment(Range(100, 400), ContigID('CONTIG_1'))
        exp = [create_pipolin(genome, scheme, f1, f2)]

        self._check_pipolins(exp, find_pipolins.run(genome))

    def test_genome9(self):
        scheme = '---at1---pol---at1---...---at1(t)---'
        genome = create_genome_from_scheme(scheme)
        f1 = PipolinFragment(Range(100, 600), ContigID('CONTIG_0'))
        f2 = PipolinFragment(Range(100, 200), ContigID('CONTIG_1'))
        exp = [create_pipolin(genome, scheme, f1, f2)]

        self._check_pipolins(exp, find_pipolins.run(genome))

    def test_genome10(self):
        scheme = '---at1---...---pol---...---at1---at1(t)---'
        genome = create_genome_from_scheme(scheme)
        f1 = PipolinFragment(Range(100, 200), ContigID('CONTIG_0'))
        f2 = PipolinFragment(Range(100, 200), ContigID('CONTIG_1'))
        f3 = PipolinFragment(Range(100, 400), ContigID('CONTIG_2'))
        exp = [create_pipolin(genome, scheme, f1, f2, f3)]

        self._check_pipolins(exp, find_pipolins.run(genome))

    def test_genome11(self):
        scheme = '---at1---pol---...---at1---...---at1(t)---'
        genome = create_genome_from_scheme(scheme)
        f1 = PipolinFragment(Range(100, 400), ContigID('CONTIG_0'))
        f2 = PipolinFragment(Range(100, 200), ContigID('CONTIG_1'))
        f3 = PipolinFragment(Range(100, 200), ContigID('CONTIG_2'))
        exp = [create_pipolin(genome, scheme, f1, f2, f3)]

        self._check_pipolins(exp, find_pipolins.run(genome))

    def test_genome12(self):
        scheme = '---at1---...---pol---at1---...---at1(t)---'
        genome = create_genome_from_scheme(scheme)
        f1 = PipolinFragment(Range(100, 200), ContigID('CONTIG_0'))
        f2 = PipolinFragment(Range(100, 400), ContigID('CONTIG_1'))
        f3 = PipolinFragment(Range(100, 200), ContigID('CONTIG_2'))
        exp = [create_pipolin(genome, scheme, f1, f2, f3)]

        self._check_pipolins(exp, find_pipolins.run(genome))

    def test_genome13(self):
        scheme = '---at1---pol---...---pol---at1(t)---'
        genome = create_genome_from_scheme(scheme)
        p1f1 = PipolinFragment(Range(100, 400), ContigID('CONTIG_0'))
        p1f2 = PipolinFragment(Range(300, 400), ContigID('CONTIG_1'))
        p2f1 = PipolinFragment(Range(100, 200), ContigID('CONTIG_1'))
        exp = [create_pipolin(genome, scheme, p1f1, p1f2), create_pipolin(genome, scheme, p2f1)]

        self._check_pipolins(exp, find_pipolins.run(genome))

    def test_genome18(self):
        scheme = '---at1---pol---at1(t)---at2---pol---at2(t)---'
        genome = create_genome_from_scheme(scheme)
        p1f1 = PipolinFragment(Range(100, 600), ContigID('CONTIG_0'))
        p2f1 = PipolinFragment(Range(800, 1300), ContigID('CONTIG_0'))
        exp = [create_pipolin(genome, scheme, p1f1), create_pipolin(genome, scheme, p2f1)]

        self._check_pipolins(exp, find_pipolins.run(genome))

    def test_genome19(self):
        scheme = '---at1---pol---at1(t)---...---at2---pol---at2(t)---'
        genome = create_genome_from_scheme(scheme)
        p1f1 = PipolinFragment(Range(100, 600), ContigID('CONTIG_0'))
        p2f1 = PipolinFragment(Range(100, 600), ContigID('CONTIG_1'))
        exp = [create_pipolin(genome, scheme, p1f1), create_pipolin(genome, scheme, p2f1)]

        self._check_pipolins(exp, find_pipolins.run(genome))

    def test_genome21(self):
        scheme = '---at1---pol---pol---at1---'
        genome = create_genome_from_scheme(scheme)
        f1 = PipolinFragment(Range(100, 800), ContigID('CONTIG_0'))
        exp = [create_pipolin(genome, scheme, f1)]

        self._check_pipolins(exp, find_pipolins.run(genome))

    def test_genome22(self):
        scheme = '---at1---at1---pol---pol---'
        genome = create_genome_from_scheme(scheme)
        f1 = PipolinFragment(Range(100, 800), ContigID('CONTIG_0'))
        exp = [create_pipolin(genome, scheme, f1)]

        self._check_pipolins(exp, find_pipolins.run(genome))

    def test_genome23(self):
        scheme = '---at1---pol---pol---at1---at1---'
        genome = create_genome_from_scheme(scheme)
        f1 = PipolinFragment(Range(100, 1000), ContigID('CONTIG_0'))
        exp = [create_pipolin(genome, scheme, f1)]

        self._check_pipolins(exp, find_pipolins.run(genome))

    def test_genome_strange1(self):
        scheme = '---at1---pol---at1---pol---at1(t)---'
        genome = create_genome_from_scheme(scheme)
        f1 = PipolinFragment(Range(100, 1000), ContigID('CONTIG_0'))
        exp = [create_pipolin(genome, scheme, f1)]

        self._check_pipolins(exp, find_pipolins.run(genome))

    def test_genome_strange2(self):
        scheme = '---at1---...---pol---at1---pol---at1(t)---'
        genome = create_genome_from_scheme(scheme)
        p1f1 = PipolinFragment(Range(100, 200), ContigID('CONTIG_0'))
        p1f2 = PipolinFragment(Range(300, 800), ContigID('CONTIG_1'))
        p2f1 = PipolinFragment(Range(100, 200), ContigID('CONTIG_1'))
        exp = [create_pipolin(genome, scheme, p1f1, p1f2), create_pipolin(genome, scheme, p2f1)]

        self._check_pipolins(exp, find_pipolins.run(genome))

    def test_genome_strange3(self):
        scheme = '---at1---pol---...---at1---pol---at1(t)---'
        genome = create_genome_from_scheme(scheme)
        p1f1 = PipolinFragment(Range(100, 200), ContigID('CONTIG_0'))
        p1f2 = PipolinFragment(Range(100, 600), ContigID('CONTIG_1'))
        p2f1 = PipolinFragment(Range(300, 400), ContigID('CONTIG_0'))

        exp = [create_pipolin(genome, scheme, p1f1, p1f2), create_pipolin(genome, scheme, p2f1)]

        self._check_pipolins(exp, find_pipolins.run(genome))

    def test_genome_strange4(self):
        scheme = '---at1---pol---at1---...---pol---at1(t)---'
        genome = create_genome_from_scheme(scheme)
        p1f1 = PipolinFragment(Range(100, 600), ContigID('CONTIG_0'))
        p1f2 = PipolinFragment(Range(300, 400), ContigID('CONTIG_1'))
        p2f1 = PipolinFragment(Range(100, 200), ContigID('CONTIG_1'))
        exp = [create_pipolin(genome, scheme, p1f1, p1f2), create_pipolin(genome, scheme, p2f1)]

        self._check_pipolins(exp, find_pipolins.run(genome))

    def test_genome_strange5(self):
        scheme = '---at1---pol---at1---pol---'
        genome = create_genome_from_scheme(scheme)
        p1f1 = PipolinFragment(Range(100, 600), ContigID('CONTIG_0'))
        p2f1 = PipolinFragment(Range(700, 800), ContigID('CONTIG_0'))
        exp = [create_pipolin(genome, scheme, p1f1), create_pipolin(genome, scheme, p2f1)]

        self._check_pipolins(exp, find_pipolins.run(genome))

    def test_genome_strange6(self):
        scheme = '---pol---at1---pol---at1---'
        genome = create_genome_from_scheme(scheme)
        p1f1 = PipolinFragment(Range(100, 200), ContigID('CONTIG_0'))
        p2f1 = PipolinFragment(Range(300, 800), ContigID('CONTIG_0'))
        exp = [create_pipolin(genome, scheme, p1f1), create_pipolin(genome, scheme, p2f1)]

        self._check_pipolins(exp, find_pipolins.run(genome))

    def test_test(self):
        scheme = '---at0---pol---at1---pol---at0---at1---pol---at2---pol---at2---pol---'
        genome = create_genome_from_scheme(scheme)
        p1f1 = PipolinFragment(Range(100, 1200), ContigID('CONTIG_0'))
        p2f1 = PipolinFragment(Range(1300, 1400), ContigID('CONTIG_0'))
        p3f1 = PipolinFragment(Range(1500, 2000), ContigID('CONTIG_0'))
        p4f1 = PipolinFragment(Range(2100, 2200), ContigID('CONTIG_0'))

        exp = [create_pipolin(genome, scheme, p1f1), create_pipolin(genome, scheme, p2f1),
               create_pipolin(genome, scheme, p3f1), create_pipolin(genome, scheme, p4f1)]
        self._check_pipolins(exp, find_pipolins.run(genome))
