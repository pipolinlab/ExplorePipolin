import unittest

from explore_pipolin.common import Contig, Genome, Feature, Range, Strand, FeatureType, AttFeature, \
    ContigID, MultiLocation
from explore_pipolin.tasks.find_atts import AttDenovoFinder


class TestAttsDenovoSearch(unittest.TestCase):
    def setUp(self) -> None:
        self.contig1_id = ContigID('CONTIG_1')
        self.contig2_id = ContigID('CONTIG_2')
        contig1 = Contig(self.contig1_id, 1000000)
        contig2 = Contig(self.contig2_id, 10000)
        self.genome = Genome('GENOME_ID', 'genome.fa', contigs=[contig1, contig2])

        pipolb1 = Feature(Range(100, 1100), Strand.FORWARD, FeatureType.PIPOLB, self.contig1_id, self.genome)
        pipolb2 = Feature(Range(3000, 4000), Strand.FORWARD, FeatureType.PIPOLB, self.contig1_id, self.genome)
        pipolb3 = Feature(Range(500, 1500), Strand.REVERSE, FeatureType.PIPOLB, self.contig2_id, self.genome)

        self.genome.features.add_features(pipolb1, pipolb2, pipolb3)

    def test_is_att_denovo(self):
        att1 = AttFeature(Range(30, 70), Strand.REVERSE, FeatureType.ATT, self.contig1_id, self.genome, 0)
        att2 = AttFeature(Range(2000, 2040), Strand.REVERSE, FeatureType.ATT, self.contig1_id, self.genome, att1.att_id)
        trna1 = Feature(Range(2030, 2060), Strand.FORWARD, FeatureType.TRNA, self.contig1_id, self.genome)
        trna2 = Feature(Range(5030, 5060), Strand.FORWARD, FeatureType.TRNA, self.contig1_id, self.genome)
        trna3 = Feature(Range(2000, 2040), Strand.REVERSE, FeatureType.TRNA, self.contig2_id, self.genome)
        self.genome.features.add_features(att1, att2, trna1, trna2, trna3)

        finder = AttDenovoFinder(self.genome, 'output_dir')
        # overlaps with att1
        self.assertFalse(finder._is_att_denovo(MultiLocation([Range(35, 65), Range(1995, 2035)], self.contig1_id)))
        self.assertTrue(finder._is_att_denovo(MultiLocation([Range(2800, 2850), Range(5000, 5050)], self.contig1_id)))
        self.assertTrue(finder._is_att_denovo(MultiLocation([Range(35, 65), Range(1995, 2035)], self.contig2_id)))
        # does not overlap with a tRNA gene
        self.assertFalse(finder._is_att_denovo(MultiLocation([Range(300, 350), Range(3000, 3050)], self.contig2_id)))
