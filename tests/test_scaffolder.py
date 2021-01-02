import unittest
from typing import Sequence

from explore_pipolin.common import Genome, Contig, Feature, Range, Strand, FeatureType
from explore_pipolin.tasks_related.scaffolding import Scaffolder


_ATT_LEN = 10
_PIPOLB_LEN = 10
_INTRAGENIC_LEN = 10

_GENOME_ID = 'GENOME'
_GENOME_FILE = 'GENOME.fa'
_CONTIG_ID_1 = 'CONTIG_1'
_CONTIG_ID_2 = 'CONTIG_2'
_CONTIG_ID_3 = 'CONTIG_3'


def add_features(genome: Genome, features: Sequence[Feature], feature_types: Sequence[FeatureType]):
    for f, t in zip(features, feature_types):
        genome.features.add_feature(f, t)


class SetUpPipolinsToScaffold(unittest.TestCase):
    def setUp(self) -> None:
        self.create_genome1()   # ---pipolb---
        self.create_genome2()   # ---att---pipolb---
        self.create_genome3()   # ---att---pipolb---att---
        self.create_genome4()   # ---att---pipolb---...---att---
        self.create_genome5()   # ---att---...---pipolb---...---att---

        self.create_genome6()   # ---att---pipolb---att---att---
        self.create_genome7()   # ---att---...---pipolb---att---att---
        # ---att---pipolb---...---att---att---
        # ---att---pipolb---att---...---att---
        # ---att---...---pipolb---...---att---att---
        # ---att---pipolb---...---att---...---att---
        # ---att---...---pipolb---att---...---att---

    @staticmethod
    def create_genome1() -> Genome:
        # ---pipolb---
        contig_length = _INTRAGENIC_LEN * 2 + _PIPOLB_LEN

        genome = Genome(_GENOME_ID, _GENOME_FILE, contigs=[Contig(_CONTIG_ID_1, contig_length)])

        pipolb = Feature(Range(11, 21), Strand.FORWARD, _CONTIG_ID_1, genome)
        genome.features.add_feature(pipolb, FeatureType.PIPOLB)
        return genome

    @staticmethod
    def create_genome2() -> Genome:
        # ---att---pipolb---
        contig_length = _INTRAGENIC_LEN * 3 + _ATT_LEN + _PIPOLB_LEN

        genome = Genome(_GENOME_ID, _GENOME_FILE, contigs=[Contig(_CONTIG_ID_1, contig_length)])

        att = Feature(Range(11, 21), Strand.FORWARD, _CONTIG_ID_1, genome)
        pipolb = Feature(Range(31, 41), Strand.FORWARD, _CONTIG_ID_1, genome)
        add_features(genome, [att, pipolb], [FeatureType.ATT, FeatureType.PIPOLB])
        return genome

    @staticmethod
    def create_genome3() -> Genome:
        # genome1: ---att---pipolb---att---
        contig_length = _INTRAGENIC_LEN * 4 + _ATT_LEN * 2 + _PIPOLB_LEN

        genome = Genome(_GENOME_ID, _GENOME_FILE, contigs=[Contig(_CONTIG_ID_1, contig_length)])

        att1 = Feature(Range(11, 21), Strand.FORWARD, _CONTIG_ID_1, genome)
        pipolb = Feature(Range(31, 41), Strand.FORWARD, _CONTIG_ID_1, genome)
        att2 = Feature(Range(51, 61), Strand.FORWARD, _CONTIG_ID_1, genome)
        add_features(genome, [att1, att2, pipolb], [FeatureType.ATT, FeatureType.ATT, FeatureType.PIPOLB])
        return genome

    @staticmethod
    def create_genome4() -> Genome:
        # genome2: ---att---pipolb---...---att---
        contig1_length = _INTRAGENIC_LEN * 3 + _ATT_LEN + _PIPOLB_LEN
        contig2_length = _INTRAGENIC_LEN * 2 + _ATT_LEN

        genome = Genome(_GENOME_ID, _GENOME_FILE, contigs=[Contig(_CONTIG_ID_1, contig1_length),
                                                           Contig(_CONTIG_ID_2, contig2_length)])

        att1 = Feature(Range(11, 21), Strand.FORWARD, _CONTIG_ID_1, genome)
        pipolb = Feature(Range(31, 41), Strand.FORWARD, _CONTIG_ID_1, genome)
        att2 = Feature(Range(11, 21), Strand.FORWARD, _CONTIG_ID_2, genome)
        add_features(genome, [att1, att2, pipolb], [FeatureType.ATT, FeatureType.ATT, FeatureType.PIPOLB])
        return genome

    @staticmethod
    def create_genome5() -> Genome:
        # genome3: ---att---...---pipolb---...---att---
        contig13_length = _INTRAGENIC_LEN * 2 + _ATT_LEN
        contig2_length = _INTRAGENIC_LEN * 2 + _PIPOLB_LEN

        genome = Genome(_GENOME_ID, _GENOME_FILE, contigs=[Contig(_CONTIG_ID_1, contig13_length),
                                                           Contig(_CONTIG_ID_2, contig2_length),
                                                           Contig(_CONTIG_ID_3, contig13_length)])
        att1 = Feature(Range(11, 21), Strand.FORWARD, _CONTIG_ID_1, genome)
        pipolb = Feature(Range(11, 21), Strand.FORWARD, _CONTIG_ID_2, genome)
        att2 = Feature(Range(11, 21), Strand.FORWARD, _CONTIG_ID_3, genome)
        add_features(genome, [att1, att2, pipolb], [FeatureType.ATT, FeatureType.ATT, FeatureType.PIPOLB])
        return genome

    @staticmethod
    def create_genome6() -> Genome:
        # ---att---pipolb---att---att---
        contig_length = _INTRAGENIC_LEN * 5 + _ATT_LEN * 3 + _PIPOLB_LEN

        genome = Genome(_GENOME_ID, _GENOME_FILE, contigs=[Contig(_CONTIG_ID_1, contig_length)])

        att1 = Feature(Range(11, 21), Strand.FORWARD, _CONTIG_ID_1, genome)
        att2 = Feature(Range(51, 61), Strand.FORWARD, _CONTIG_ID_1, genome)
        att3 = Feature(Range(71, 81), Strand.FORWARD, _CONTIG_ID_1, genome)
        pipolb = Feature(Range(31, 41), Strand.FORWARD, _CONTIG_ID_1, genome)
        add_features(genome, [att1, att2, att3, pipolb],
                     [FeatureType.ATT, FeatureType.ATT, FeatureType.ATT, FeatureType.PIPOLB])
        return genome

    @staticmethod
    def create_genome7() -> Genome:
        # ---att---...---pipolb---att---att---
        contig1_length = _INTRAGENIC_LEN * 2 + _ATT_LEN
        contig2_length = _INTRAGENIC_LEN * 4 + _ATT_LEN * 2 + _PIPOLB_LEN

        genome = Genome(_GENOME_ID, _GENOME_FILE, contigs=[Contig(_CONTIG_ID_1, contig1_length),
                                                           Contig(_CONTIG_ID_2, contig2_length)])

        att1 = Feature(Range(11, 21), Strand.FORWARD, _CONTIG_ID_1, genome)
        att2 = Feature(Range(31, 41), Strand.FORWARD, _CONTIG_ID_2, genome)
        att3 = Feature(Range(51, 61), Strand.FORWARD, _CONTIG_ID_2, genome)
        pipolb = Feature(Range(11, 21), Strand.FORWARD, _CONTIG_ID_2, genome)
        add_features(genome, [att1, att2, att3, pipolb],
                     [FeatureType.ATT, FeatureType.ATT, FeatureType.ATT, FeatureType.PIPOLB])
        return genome


class TestScaffolder(SetUpPipolinsToScaffold):
    def test_genome1(self):
        with self.assertRaises(NotImplementedError):
            Scaffolder(self.create_genome1()).scaffold()

    def test_genome2(self):
        with self.assertRaises(NotImplementedError):
            Scaffolder(self.create_genome2()).scaffold()

    def test_genome3(self):
        with self.assertRaises(NotImplementedError):
            Scaffolder(self.create_genome3()).scaffold()

    def test_genome4(self):
        Scaffolder(self.create_genome4()).scaffold()

    def test_genome5(self):
        genome5 = self.create_genome5()
        with self.assertRaises(NotImplementedError):
            Scaffolder(genome5).scaffold()

        target_trna = Feature(Range(16, 26), Strand.REVERSE, _CONTIG_ID_3, genome5)
        genome5.features.add_feature(target_trna, FeatureType.TARGET_TRNA)
        Scaffolder(genome5).scaffold()

    def test_genome6(self):
        with self.assertRaises(NotImplementedError):
            Scaffolder(self.create_genome6()).scaffold()

    def test_genome7(self):
        Scaffolder(self.create_genome7()).scaffold()
