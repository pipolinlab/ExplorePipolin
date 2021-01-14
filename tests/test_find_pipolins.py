import unittest
from typing import Sequence

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

    _add_features_to_genome(contigs_schemes, genome)

    return genome


def _create_genome_with_contigs(contigs_schemes: Sequence[str]) -> Genome:
    contigs = []
    for i in range(len(contigs_schemes)):
        contigs.append(Contig(ContigID(f'CONTIG_{i}'), len(contigs_schemes[i]) // 3 * 100))
    genome = Genome(_GENOME_ID, _GENOME_FILE, contigs)
    return genome


def _add_features_to_genome(contigs_schemes, genome):
    for i, sch in enumerate(contigs_schemes):
        for triplet_start in range(0, len(sch), 3):
            triplet = sch[triplet_start: triplet_start + 3]
            feature_start = (triplet_start // 3 * 100 - 5) if triplet == _TRNA else (triplet_start // 3 * 100)
            feature = Feature(Range(feature_start, feature_start + 100),
                              Strand.FORWARD, ContigID(f'CONTIG_{i}'), genome)
            if triplet[:2] == _ATT:
                feature = AttFeature(Range(feature_start, feature_start + 100), Strand.FORWARD,
                                     ContigID(f'CONTIG_{i}'), genome, att_id=int(triplet[2]))
                genome.features.add_features(feature, feature_type=FeatureType.ATT)
            if triplet == _PIPOLB:
                genome.features.add_features(feature, feature_type=FeatureType.PIPOLB)
            if triplet == _TRNA:
                genome.features.add_features(feature, feature_type=FeatureType.TARGET_TRNA)


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
        genome = create_genome_from_scheme('---pol---')
        exp = [Pipolin.from_fragments(PipolinFragment(Range(0, 300), ContigID('CONTIG_0')))]

        self._check_pipolins(exp, find_pipolins.run(genome))

    def test_genome1_1(self):
        genome = create_genome_from_scheme('---pol---pol---')
        exp = [Pipolin.from_fragments(PipolinFragment(Range(0, 500), ContigID('CONTIG_0')))]

        self._check_pipolins(exp, find_pipolins.run(genome))

    def test_genome1_2(self):
        pass   # genome = create_genome_from_scheme('---pol===pol---')

    def test_genome2(self):
        genome = create_genome_from_scheme('---at1---pol---')
        exp = [Pipolin.from_fragments(PipolinFragment(Range(0, 500), ContigID('CONTIG_0')))]

        self._check_pipolins(exp, find_pipolins.run(genome))

    def test_genome3(self):
        genome = create_genome_from_scheme('---at1---pol---at1---')
        exp = [Pipolin.from_fragments(PipolinFragment(Range(50, 650), ContigID('CONTIG_0')))]

        self._check_pipolins(exp, find_pipolins.run(genome))

    def test_genome4(self):
        genome = create_genome_from_scheme('---at1---pol---...---at1---')
        f1 = PipolinFragment(Range(0, 500), ContigID('CONTIG_0'))
        f2 = PipolinFragment(Range(0, 300), ContigID('CONTIG_1'))
        exp = [Pipolin.from_fragments(f1, f2)]

        self._check_pipolins(exp, find_pipolins.run(genome))

    def test_genome5(self):
        genome = create_genome_from_scheme('---at1---...---pol---...---at1---')
        f1 = PipolinFragment(Range(0, 300), ContigID('CONTIG_0'))
        f2 = PipolinFragment(Range(0, 300), ContigID('CONTIG_1'))
        f3 = PipolinFragment(Range(0, 300), ContigID('CONTIG_2'))
        exp = [Pipolin.from_fragments(f1, f2, f3)]

        self._check_pipolins(exp, find_pipolins.run(genome))

    def test_genome5_trna(self):
        genome = create_genome_from_scheme('---at1---...---pol---...---at1(t)---')
        f1 = PipolinFragment(Range(0, 300), ContigID('CONTIG_0'))
        f2 = PipolinFragment(Range(0, 300), ContigID('CONTIG_1'))
        f3 = PipolinFragment(Range(0, 250), ContigID('CONTIG_2'))   # (0, 400)
        exp = [Pipolin.from_fragments(f1, f2, f3)]

        self._check_pipolins(exp, find_pipolins.run(genome))

    def test_genome6(self):
        genome = create_genome_from_scheme('---at1---pol---at1---at1(t)---')
        exp = [Pipolin.from_fragments(PipolinFragment(Range(50, 850), ContigID('CONTIG_0')))]
        self._check_pipolins(exp, find_pipolins.run(genome))

    def test_genome7(self):
        genome = create_genome_from_scheme('---at1---...---pol---at1---at1(t)---')
        f1 = PipolinFragment(Range(0, 300), ContigID('CONTIG_0'))
        f2 = PipolinFragment(Range(0, 650), ContigID('CONTIG_1'))   # (0, 800)
        exp = [Pipolin.from_fragments(f1, f2)]

        self._check_pipolins(exp, find_pipolins.run(genome))

    def test_genome8(self):
        genome = create_genome_from_scheme('---at1---pol---...---at1---at1(t)---')
        f1 = PipolinFragment(Range(0, 500), ContigID('CONTIG_0'))
        f2 = PipolinFragment(Range(0, 450), ContigID('CONTIG_1'))   # (0, 600)
        exp = [Pipolin.from_fragments(f1, f2)]

        self._check_pipolins(exp, find_pipolins.run(genome))

    def test_genome9(self):
        genome = create_genome_from_scheme('---at1---pol---at1---...---at1(t)---')
        f1 = PipolinFragment(Range(0, 700), ContigID('CONTIG_0'))
        f2 = PipolinFragment(Range(0, 250), ContigID('CONTIG_1'))   # (0, 400)
        exp = [Pipolin.from_fragments(f1, f2)]

        self._check_pipolins(exp, find_pipolins.run(genome))

    def test_genome10(self):
        genome = create_genome_from_scheme('---at1---...---pol---...---at1---at1(t)---')
        f1 = PipolinFragment(Range(0, 300), ContigID('CONTIG_0'))
        f2 = PipolinFragment(Range(0, 300), ContigID('CONTIG_1'))
        f3 = PipolinFragment(Range(0, 450), ContigID('CONTIG_2'))   # (0, 600)
        exp = [Pipolin.from_fragments(f1, f2, f3)]

        self._check_pipolins(exp, find_pipolins.run(genome))

    def test_genome11(self):
        genome = create_genome_from_scheme('---at1---pol---...---at1---...---at1(t)---')
        f1 = PipolinFragment(Range(0, 500), ContigID('CONTIG_0'))
        f2 = PipolinFragment(Range(0, 300), ContigID('CONTIG_1'))
        f3 = PipolinFragment(Range(0, 250), ContigID('CONTIG_2'))   # (0, 400)
        exp = [Pipolin.from_fragments(f1, f2, f3)]

        self._check_pipolins(exp, find_pipolins.run(genome))

    def test_genome12(self):
        genome = create_genome_from_scheme('---at1---...---pol---at1---...---at1(t)---')
        f1 = PipolinFragment(Range(0, 300), ContigID('CONTIG_0'))
        f2 = PipolinFragment(Range(0, 500), ContigID('CONTIG_1'))
        f3 = PipolinFragment(Range(0, 250), ContigID('CONTIG_2'))   # (0, 400)
        exp = [Pipolin.from_fragments(f1, f2, f3)]

        self._check_pipolins(exp, find_pipolins.run(genome))

    def test_genome13(self):
        genome = create_genome_from_scheme('---at1---pol---...---pol---at1(t)---')
        f1 = PipolinFragment(Range(0, 500), ContigID('CONTIG_0'))
        f2 = PipolinFragment(Range(0, 450), ContigID('CONTIG_1'))   # (0, 600)

        self._check_pipolins([Pipolin.from_fragments(f1, f2)], find_pipolins.run(genome))

    def test_genome18(self):
        genome = create_genome_from_scheme('---at1---pol---at1(t)---at2---pol---at2(t)---')
        p1 = Pipolin.from_fragments(PipolinFragment(Range(50, 650), ContigID('CONTIG_0')))
        p2 = Pipolin.from_fragments(PipolinFragment(Range(750, 1350), ContigID('CONTIG_0')))

        self._check_pipolins([p1, p2], find_pipolins.run(genome))

    def test_genome19(self):
        genome = create_genome_from_scheme('---at1---pol---at1(t)---...---at2---pol---at2(t)---')
        p1 = Pipolin.from_fragments(PipolinFragment(Range(50, 650), ContigID('CONTIG_0')))
        p2 = Pipolin.from_fragments(PipolinFragment(Range(50, 650), ContigID('CONTIG_1')))

        self._check_pipolins([p1, p2], find_pipolins.run(genome))

    def test_genome21(self):
        genome = create_genome_from_scheme('---at1---pol---pol---at1---')
        exp = [Pipolin.from_fragments(PipolinFragment(Range(50, 850), ContigID('CONTIG_0')))]

        self._check_pipolins(exp, find_pipolins.run(genome))

    def test_genome22(self):
        genome = create_genome_from_scheme('---at1---at1---pol---pol---')
        exp = [Pipolin.from_fragments(PipolinFragment(Range(0, 900), ContigID('CONTIG_0')))]

        self._check_pipolins(exp, find_pipolins.run(genome))

    def test_genome23(self):
        genome = create_genome_from_scheme('---at1---pol---pol---at1---at1---')
        exp = [Pipolin.from_fragments(PipolinFragment(Range(50, 1050), ContigID('CONTIG_0')))]

        self._check_pipolins(exp, find_pipolins.run(genome))

    def test_genome_strange1(self):
        genome = create_genome_from_scheme('---at1---pol---at1---pol---at1(t)---')
        exp = [Pipolin.from_fragments(PipolinFragment(Range(50, 1050), ContigID('CONTIG_0')))]

        self._check_pipolins(exp, find_pipolins.run(genome))

    def test_genome_strange2(self):
        genome = create_genome_from_scheme('---at1---...---pol---at1---pol---at1(t)---')
        f1 = PipolinFragment(Range(0, 300), ContigID('CONTIG_0'))
        f2 = PipolinFragment(Range(0, 850), ContigID('CONTIG_1'))   # (0, 1000)

        self._check_pipolins([Pipolin.from_fragments(f1, f2)], find_pipolins.run(genome))

    def test_genome_strange3(self):
        genome = create_genome_from_scheme('---at1---pol---...---at1---pol---at1(t)---')
        f1 = PipolinFragment(Range(0, 500), ContigID('CONTIG_0'))
        f2 = PipolinFragment(Range(0, 800), ContigID('CONTIG_1'))   # or 650? ERROR HAPPENS!!!
        self._check_pipolins([Pipolin.from_fragments(f1, f2)], find_pipolins.run(genome))

    def test_genome_strange4(self):
        genome = create_genome_from_scheme('---at1---pol---at1---...---pol---at1(t)---')
        f1 = PipolinFragment(Range(0, 700), ContigID('CONTIG_0'))
        f2 = PipolinFragment(Range(0, 450), ContigID('CONTIG_1'))   # (0, 600)

        self._check_pipolins([Pipolin.from_fragments(f1, f2)], find_pipolins.run(genome))

    def test_genome_strange5(self):
        genome = create_genome_from_scheme('---at1---pol---at1---pol---')
        p1 = Pipolin.from_fragments(PipolinFragment(Range(50, 650), ContigID('CONTIG_0')))
        # It's weird here as the default inflate size is 100000
        p2 = Pipolin.from_fragments(PipolinFragment(Range(0, 900), ContigID('CONTIG_0')))

        self._check_pipolins([p1, p2], find_pipolins.run(genome))
