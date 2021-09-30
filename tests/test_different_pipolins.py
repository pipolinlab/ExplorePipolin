import shutil
import tempfile
import unittest
from typing import Sequence, List

from explore_pipolin.common import Genome, Contig, Feature, Range, Strand, FeatureType, \
    Pipolin, PipolinFragment, AttFeature, ContigID, AttType
import explore_pipolin.settings as settings
from explore_pipolin.settings import GlobalSettings
from explore_pipolin.tasks.find_pipolins import PipolinFinder, find_pipolins

from explore_pipolin.tasks.reconstruct_pipolins import reconstruct_pipolins

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
    genome.features.add_features(*features)

    return genome


def _create_genome_with_contigs(contigs_schemes: Sequence[str]) -> Genome:
    contigs = []
    for i in range(len(contigs_schemes)):
        contigs.append(Contig(ContigID(f'CONTIG_{i}'), len(contigs_schemes[i]) // 3 * 100))
    genome = Genome(_GENOME_ID, _GENOME_FILE, contigs)
    return genome


def _create_features_for_genome(
        genome: Genome, contigs_schemes: List[str]) -> Sequence[Feature]:

    features = []
    for i, sch in enumerate(contigs_schemes):
        for triplet_start in range(0, len(sch), 3):
            triplet = sch[triplet_start: triplet_start + 3]
            feature_start = (triplet_start // 3 * 100 - 60) if triplet == _TRNA else (triplet_start // 3 * 100)
            if triplet[:2] == _ATT:
                feature = AttFeature(Range(feature_start, feature_start + 100), Strand.FORWARD, FeatureType.ATT,
                                     ContigID(f'CONTIG_{i}'), genome,
                                     att_id=int(triplet[2]), att_type=AttType.CONSERVED)
                features.append(feature)
            else:
                if triplet == _PIPOLB:
                    feature = Feature(Range(feature_start, feature_start + 100), Strand.FORWARD,
                                      FeatureType.PIPOLB, ContigID(f'CONTIG_{i}'), genome)
                    features.append(feature)
                if triplet == _TRNA:
                    feature = Feature(Range(feature_start, feature_start + 100), Strand.REVERSE,
                                      FeatureType.TARGET_TRNA, ContigID(f'CONTIG_{i}'), genome)
                    features.append(feature)

    return tuple(features)


def create_pipolin(scheme: str, *fragments: PipolinFragment) -> Pipolin:
    new_fragments = []
    for fragment in fragments:
        features = _create_features_for_pipolin_fragment(scheme, fragment)
        new_fragments.append(PipolinFragment(fragment.location, fragment.contig_id, fragment.genome, tuple(features)))
    return Pipolin.from_fragments(*new_fragments)


def _create_features_for_pipolin_fragment(
        scheme: str, fragment: PipolinFragment) -> Sequence[Feature]:

    contigs_schemes = scheme.split(_GAP)
    features = _create_features_for_genome(fragment.genome, contigs_schemes)

    fragment_features = []
    for feature in features:
        if feature.contig_id == fragment.contig_id:
            if feature.start >= fragment.start and feature.end <= fragment.end:
                fragment_features.append(feature)

    return fragment_features


class TestPipolinFinder(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp_dir = tempfile.TemporaryDirectory()
        self.settings = GlobalSettings.create_instance(
            out_dir_prefix=None,
            user_defined_out_dir=self.tmp_dir.name,
            user_defined_profile=None,
            user_defined_att=None,
            percent_identity=85,
            max_inflate=30_000,
            user_defined_proteins=None,
            prokka_cpus=0,
            skip_colours=False,
        )
        settings.set_instance(self.settings)

    def tearDown(self) -> None:
        self.tmp_dir.cleanup()

    def _check_pipolins(self, expected: Sequence[Pipolin], obtained: Sequence[Pipolin], ordered=False):
        self.assertEqual(len(expected), len(obtained))
        if not ordered:
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

    def check_found_pipolins(self, genome, scheme, *pipolin: Sequence[PipolinFragment]):
        exp_found = [create_pipolin(scheme, *i) for i in pipolin]
        obt_found = PipolinFinder(genome).find_pipolins()
        self._check_pipolins(exp_found, obt_found)

        return obt_found

    def test_genome1(self):
        scheme = '---pol---'
        genome = create_genome_from_scheme(scheme)
        f1 = PipolinFragment(Range(100, 200), ContigID('CONTIG_0'), genome)
        self.check_found_pipolins(genome, scheme, [f1])

        reconstruct_pipolins.run(genome=genome, pipolins=[create_pipolin(scheme, f1)])

    def test_genome1_1(self):
        scheme = '---pol---pol---'
        genome = create_genome_from_scheme(scheme)
        p1f1 = PipolinFragment(Range(100, 200), ContigID('CONTIG_0'), genome)
        p2f1 = PipolinFragment(Range(300, 400), ContigID('CONTIG_0'), genome)
        self.check_found_pipolins(genome, scheme, [p1f1], [p2f1])

        reconstruct_pipolins.run(
            genome=genome, pipolins=[create_pipolin(scheme, p1f1), create_pipolin(scheme, p2f1)]
        )

    def test_genome2(self):
        scheme = '---at1---pol---'
        genome = create_genome_from_scheme(scheme)
        f1 = PipolinFragment(Range(100, 400), ContigID('CONTIG_0'), genome)
        self.check_found_pipolins(genome, scheme, [f1])

        reconstruct_pipolins.run(genome=genome, pipolins=[create_pipolin(scheme, f1)])

    def test_genome3(self):
        scheme = '---at1---pol---at1---'
        genome = create_genome_from_scheme(scheme)
        f1 = PipolinFragment(Range(100, 600), ContigID('CONTIG_0'), genome)
        self.check_found_pipolins(genome, scheme, [f1])

        reconstruct_pipolins.run(genome=genome, pipolins=[create_pipolin(scheme, f1)])

    def test_genome3_strange(self):
        scheme = '---at1(t)---pol---at1(t)---'
        genome = create_genome_from_scheme(scheme)
        f1 = PipolinFragment(Range(100, 740), ContigID('CONTIG_0'), genome)
        self.check_found_pipolins(genome, scheme, [f1])

        reconstruct_pipolins.run(genome=genome, pipolins=[create_pipolin(scheme, f1)])

    def test_genome4(self):
        scheme = '---at1---pol---...---at1---'
        genome = create_genome_from_scheme(scheme)
        f1 = PipolinFragment(Range(100, 400), ContigID('CONTIG_0'), genome)
        f2 = PipolinFragment(Range(100, 200), ContigID('CONTIG_1'), genome)
        self.check_found_pipolins(genome, scheme, [f1, f2])

        reconstruct_pipolins.run(genome=genome, pipolins=[create_pipolin(scheme, f1, f2)])

    def test_genome5(self):
        scheme = '---at1---...---pol---...---at1---'
        genome = create_genome_from_scheme(scheme)
        f1 = PipolinFragment(Range(100, 200), ContigID('CONTIG_0'), genome)
        f2 = PipolinFragment(Range(100, 200), ContigID('CONTIG_1'), genome)
        f3 = PipolinFragment(Range(100, 200), ContigID('CONTIG_2'), genome)
        self.check_found_pipolins(genome, scheme, [f1, f2, f3])

        reconstruct_pipolins.run(genome=genome, pipolins=[create_pipolin(scheme, f1, f2, f3)])

    def test_genome5_trna(self):
        scheme = '---at1---...---pol---...---at1(t)---'
        genome = create_genome_from_scheme(scheme)
        f1 = PipolinFragment(Range(100, 200), ContigID('CONTIG_0'), genome)
        f2 = PipolinFragment(Range(100, 200), ContigID('CONTIG_1'), genome)
        f3 = PipolinFragment(Range(100, 240), ContigID('CONTIG_2'), genome)
        self.check_found_pipolins(genome, scheme, [f1, f2, f3])

        reconstruct_pipolins.run(genome=genome, pipolins=[create_pipolin(scheme, f1, f2, f3)])

    def test_genome6(self):
        scheme = '---at1---pol---at1---at1(t)---'
        genome = create_genome_from_scheme(scheme)
        f1 = PipolinFragment(Range(100, 840), ContigID('CONTIG_0'), genome)
        self.check_found_pipolins(genome, scheme, [f1])

        reconstruct_pipolins.run(genome=genome, pipolins=[create_pipolin(scheme, f1)])

    def test_genome7(self):
        scheme = '---at1---...---pol---at1---at1(t)---'
        genome = create_genome_from_scheme(scheme)
        f1 = PipolinFragment(Range(100, 200), ContigID('CONTIG_0'), genome)
        f2 = PipolinFragment(Range(100, 640), ContigID('CONTIG_1'), genome)
        self.check_found_pipolins(genome, scheme, [f1, f2])

        reconstruct_pipolins.run(genome=genome, pipolins=[create_pipolin(scheme, f1, f2)])

    def test_genome8(self):
        scheme = '---at1---pol---...---at1---at1(t)---'
        genome = create_genome_from_scheme(scheme)
        f1 = PipolinFragment(Range(100, 400), ContigID('CONTIG_0'), genome)
        f2 = PipolinFragment(Range(100, 440), ContigID('CONTIG_1'), genome)
        self.check_found_pipolins(genome, scheme, [f1, f2])

        reconstruct_pipolins.run(genome=genome, pipolins=[create_pipolin(scheme, f1, f2)])

    def test_genome9(self):
        scheme = '---at1---pol---at1---...---at1(t)---'
        genome = create_genome_from_scheme(scheme)
        f1 = PipolinFragment(Range(100, 600), ContigID('CONTIG_0'), genome)
        f2 = PipolinFragment(Range(100, 240), ContigID('CONTIG_1'), genome)
        self.check_found_pipolins(genome, scheme, [f1, f2])

        reconstruct_pipolins.run(genome=genome, pipolins=[create_pipolin(scheme, f1, f2)])

    def test_genome10(self):
        scheme = '---at1---...---pol---...---at1---at1(t)---'
        genome = create_genome_from_scheme(scheme)
        f1 = PipolinFragment(Range(100, 200), ContigID('CONTIG_0'), genome)
        f2 = PipolinFragment(Range(100, 200), ContigID('CONTIG_1'), genome)
        f3 = PipolinFragment(Range(100, 440), ContigID('CONTIG_2'), genome)
        self.check_found_pipolins(genome, scheme, [f1, f2, f3])

        reconstruct_pipolins.run(genome=genome, pipolins=[create_pipolin(scheme, f1, f2, f3)])

    def test_genome11(self):
        scheme = '---at1---pol---...---at1---...---at1(t)---'
        genome = create_genome_from_scheme(scheme)
        f1 = PipolinFragment(Range(100, 400), ContigID('CONTIG_0'), genome)
        f2 = PipolinFragment(Range(100, 200), ContigID('CONTIG_1'), genome)
        f3 = PipolinFragment(Range(100, 240), ContigID('CONTIG_2'), genome)
        self.check_found_pipolins(genome, scheme, [f1, f2, f3])

        reconstruct_pipolins.run(genome=genome, pipolins=[create_pipolin(scheme, f1, f2, f3)])

    def test_genome12(self):
        scheme = '---at1---...---pol---at1---...---at1(t)---'
        genome = create_genome_from_scheme(scheme)
        f1 = PipolinFragment(Range(100, 200), ContigID('CONTIG_0'), genome)
        f2 = PipolinFragment(Range(100, 400), ContigID('CONTIG_1'), genome)
        f3 = PipolinFragment(Range(100, 240), ContigID('CONTIG_2'), genome)
        self.check_found_pipolins(genome, scheme, [f1, f2, f3])

        reconstruct_pipolins.run(genome=genome, pipolins=[create_pipolin(scheme, f1, f2, f3)])

    # fails from time to time (it's expected)
    def test_genome13(self):
        scheme = '---at1---pol---...---pol---at1(t)---'
        genome = create_genome_from_scheme(scheme)
        p1f1 = PipolinFragment(Range(100, 400), ContigID('CONTIG_0'), genome)
        p1f2 = PipolinFragment(Range(300, 440), ContigID('CONTIG_1'), genome)
        p2f1 = PipolinFragment(Range(100, 200), ContigID('CONTIG_1'), genome)
        self.check_found_pipolins(genome, scheme, [p1f1, p1f2], [p2f1])

        reconstruct_pipolins.run(genome=genome, pipolins=[create_pipolin(scheme, p1f1, p1f2),
                                                          create_pipolin(scheme, p2f1)])

    def test_genome18(self):
        scheme = '---at1---pol---at1(t)---at2---pol---at2(t)---'
        genome = create_genome_from_scheme(scheme)
        p1f1 = PipolinFragment(Range(100, 640), ContigID('CONTIG_0'), genome)
        p2f1 = PipolinFragment(Range(800, 1340), ContigID('CONTIG_0'), genome)
        self.check_found_pipolins(genome, scheme, [p1f1], [p2f1])

        reconstruct_pipolins.run(genome=genome, pipolins=[create_pipolin(scheme, p1f1), create_pipolin(scheme, p2f1)])

    def test_genome19(self):
        scheme = '---at1---pol---at1(t)---...---at2---pol---at2(t)---'
        genome = create_genome_from_scheme(scheme)
        p1f1 = PipolinFragment(Range(100, 640), ContigID('CONTIG_0'), genome)
        p2f1 = PipolinFragment(Range(100, 640), ContigID('CONTIG_1'), genome)
        self.check_found_pipolins(genome, scheme, [p1f1], [p2f1])

        reconstruct_pipolins.run(genome=genome, pipolins=[create_pipolin(scheme, p1f1), create_pipolin(scheme, p2f1)])

    def test_genome21(self):
        scheme = '---at1---pol---pol---at1---'
        genome = create_genome_from_scheme(scheme)
        f1 = PipolinFragment(Range(100, 800), ContigID('CONTIG_0'), genome)
        self.check_found_pipolins(genome, scheme, [f1])

        reconstruct_pipolins.run(genome=genome, pipolins=[create_pipolin(scheme, f1)])

    def test_genome22(self):
        scheme = '---at1---at1---pol---pol---'
        genome = create_genome_from_scheme(scheme)
        f1 = PipolinFragment(Range(100, 800), ContigID('CONTIG_0'), genome)
        self.check_found_pipolins(genome, scheme, [f1])

        reconstruct_pipolins.run(genome=genome, pipolins=[create_pipolin(scheme, f1)])

    def test_genome23(self):
        scheme = '---at1---pol---pol---at1---at1---'
        genome = create_genome_from_scheme(scheme)
        f1 = PipolinFragment(Range(100, 1000), ContigID('CONTIG_0'), genome)
        self.check_found_pipolins(genome, scheme, [f1])

        reconstruct_pipolins.run(genome=genome, pipolins=[create_pipolin(scheme, f1)])

    def test_genome_strange1(self):
        scheme = '---at1---pol---at1---pol---at1(t)---'
        genome = create_genome_from_scheme(scheme)
        f1 = PipolinFragment(Range(100, 1040), ContigID('CONTIG_0'), genome)
        self.check_found_pipolins(genome, scheme, [f1])

        reconstruct_pipolins.run(genome=genome, pipolins=[create_pipolin(scheme, f1)])

    def test_genome_strange2(self):
        scheme = '---at1---...---pol---at1---pol---at1(t)---'
        genome = create_genome_from_scheme(scheme)
        p1f1 = PipolinFragment(Range(100, 200), ContigID('CONTIG_0'), genome)
        p1f2 = PipolinFragment(Range(300, 840), ContigID('CONTIG_1'), genome)
        p2f1 = PipolinFragment(Range(100, 200), ContigID('CONTIG_1'), genome)
        self.check_found_pipolins(genome, scheme, [p1f1, p1f2], [p2f1])

        reconstruct_pipolins.run(
            genome=genome,
            pipolins=[create_pipolin(scheme, p1f1, p1f2), create_pipolin(scheme, p2f1)]
        )

    def test_genome_strange3(self):
        scheme = '---at1---pol---...---at1---pol---at1(t)---'
        genome = create_genome_from_scheme(scheme)
        p1f1 = PipolinFragment(Range(100, 200), ContigID('CONTIG_0'), genome)
        p1f2 = PipolinFragment(Range(100, 640), ContigID('CONTIG_1'), genome)
        p2f1 = PipolinFragment(Range(300, 400), ContigID('CONTIG_0'), genome)
        self.check_found_pipolins(genome, scheme, [p1f1, p1f2], [p2f1])

        reconstruct_pipolins.run(genome=genome, pipolins=[create_pipolin(scheme, p1f1, p1f2),
                                                          create_pipolin(scheme, p2f1)])

    def test_genome_strange4(self):
        scheme = '---at1---pol---at1---...---pol---at1(t)---'
        genome = create_genome_from_scheme(scheme)
        p1f1 = PipolinFragment(Range(100, 600), ContigID('CONTIG_0'), genome)
        p1f2 = PipolinFragment(Range(300, 440), ContigID('CONTIG_1'), genome)
        p2f1 = PipolinFragment(Range(100, 200), ContigID('CONTIG_1'), genome)
        self.check_found_pipolins(genome, scheme, [p1f1, p1f2], [p2f1])

        reconstruct_pipolins.run(genome=genome, pipolins=[create_pipolin(scheme, p1f1, p1f2),
                                                          create_pipolin(scheme, p2f1)])

    def test_genome_strange5(self):
        scheme = '---at1---pol---at1---pol---'
        genome = create_genome_from_scheme(scheme)
        p1f1 = PipolinFragment(Range(100, 600), ContigID('CONTIG_0'), genome)
        p2f1 = PipolinFragment(Range(700, 800), ContigID('CONTIG_0'), genome)
        self.check_found_pipolins(genome, scheme, [p1f1], [p2f1])

        reconstruct_pipolins.run(genome=genome, pipolins=[create_pipolin(scheme, p1f1), create_pipolin(scheme, p2f1)])

    def test_genome_strange6(self):
        scheme = '---pol---at1---pol---at1---'
        genome = create_genome_from_scheme(scheme)
        # p1f1 = PipolinFragment(Range(100, 200), ContigID('CONTIG_0'), genome)
        # p2f1 = PipolinFragment(Range(300, 800), ContigID('CONTIG_0'), genome)
        # self.check_found_pipolins(genome, scheme, [p1f1], [p2f1])
        found = find_pipolins.run(genome=genome)
        reconstruct_pipolins.run(genome=genome, pipolins=found)

    def test_test(self):
        scheme = '---at0---pol---at1---pol---at0---at1---pol---at2---pol---at2---pol---'
        genome = create_genome_from_scheme(scheme)
        p1f1 = PipolinFragment(Range(100, 1200), ContigID('CONTIG_0'), genome)
        p2f1 = PipolinFragment(Range(1300, 1400), ContigID('CONTIG_0'), genome)
        p3f1 = PipolinFragment(Range(1500, 2000), ContigID('CONTIG_0'), genome)
        p4f1 = PipolinFragment(Range(2100, 2200), ContigID('CONTIG_0'), genome)
        self.check_found_pipolins(genome, scheme, [p1f1], [p2f1], [p3f1], [p4f1])

        reconstruct_pipolins.run(genome=genome, pipolins=[create_pipolin(scheme, p1f1),
                                                          create_pipolin(scheme, p2f1),
                                                          create_pipolin(scheme, p3f1),
                                                          create_pipolin(scheme, p4f1)])

    def test_two_contigs_with_ttrnas(self):
        scheme = '---at2(t)---...---pol---at2(t)'
        genome = create_genome_from_scheme(scheme)
        f1 = PipolinFragment(Range(100, 240), ContigID('CONTIG_0'), genome)
        f2 = PipolinFragment(Range(100, 440), ContigID('CONTIG_1'), genome)
        self.check_found_pipolins(genome, scheme, [f1, f2])

        reconstruct_pipolins.run(genome=genome, pipolins=[create_pipolin(scheme, f1, f2)])

    def test_contig_with_two_ttrnas(self):
        scheme = '---at2---at1---...---pol---at2(t)---at1(t)---'   # ???
        genome = create_genome_from_scheme(scheme)
        f1 = PipolinFragment(Range(100, 200), ContigID('CONTIG_0'), genome)
        f2 = PipolinFragment(Range(100, 440), ContigID('CONTIG_1'), genome)
        self.check_found_pipolins(genome, scheme, [f1, f2])

        reconstruct_pipolins.run(genome=genome, pipolins=[create_pipolin(scheme, f1, f2)])
