import unittest
from typing import List

from explore_pipolin.common import Contig, Genome, Feature, Range, Strand, FeatureType, RangePair, AttFeature
from explore_pipolin.tasks_related.misc import get_ranges_around_pipolbs
from explore_pipolin.tasks_related.atts_denovo_search import AttsDenovoFinder


class TestAttsDenovoSearch(unittest.TestCase):
    def setUp(self) -> None:
        self.contig1_id = 'CONTIG_1'
        self.contig2_id = 'CONTIG_2'
        contig1 = Contig(self.contig1_id, 1000000)
        contig2 = Contig(self.contig2_id, 10000)
        self.genome = Genome('GENOME_ID', 'genome.fa', contigs=[contig1, contig2])

        pipolb1 = Feature(Range(100, 1100), Strand.FORWARD, self.contig1_id, self.genome)
        pipolb2 = Feature(Range(3000, 4000), Strand.FORWARD, self.contig1_id, self.genome)
        pipolb3 = Feature(Range(500, 1500), Strand.REVERSE, self.contig2_id, self.genome)

        self.genome.features.add_features(pipolb1, pipolb2, pipolb3, feature_type=FeatureType.PIPOLB)

    def test_get_ranges_around_pipolbs(self):
        range_pair1 = RangePair(left=Range(0, 100), right=Range(1100, 1100 + 100000), contig_id=self.contig1_id)
        range_pair2 = RangePair(left=Range(0, 3000), right=Range(4000, 4000 + 100000), contig_id=self.contig1_id)
        range_pair3 = RangePair(left=Range(0, 100), right=Range(4000, 4000 + 100000), contig_id=self.contig1_id)
        range_pair4 = RangePair(left=Range(0, 500), right=Range(1500, 10000), contig_id=self.contig2_id)

        obt: List[RangePair] = get_ranges_around_pipolbs(self.genome)
        self.assertCountEqual([range_pair1, range_pair2, range_pair3, range_pair4], obt)

    def test_is_att_denovo(self):
        att1 = AttFeature(Range(30, 70), Strand.REVERSE, 'repeat_01', self.contig1_id, self.genome)
        att2 = AttFeature(Range(2000, 2040), Strand.REVERSE, 'repeat_01', self.contig1_id, self.genome)
        trna1 = Feature(Range(2030, 2060), Strand.FORWARD, self.contig1_id, self.genome)
        trna2 = Feature(Range(5030, 5060), Strand.FORWARD, self.contig1_id, self.genome)
        trna3 = Feature(Range(2000, 2040), Strand.REVERSE, self.contig2_id, self.genome)
        self.genome.features.add_features(att1, att2, feature_type=FeatureType.ATT)
        self.genome.features.add_features(trna1, trna2, trna3, feature_type=FeatureType.TRNA)

        finder = AttsDenovoFinder(self.genome, 'output_dir')
        self.assertFalse(finder.is_att_denovo(RangePair(Range(35, 65), Range(2005, 2035), self.contig1_id)))
        self.assertTrue(finder.is_att_denovo(RangePair(Range(2800, 2850), Range(5000, 5050), self.contig1_id)))
        self.assertTrue(finder.is_att_denovo(RangePair(Range(35, 65), Range(2005, 2035), self.contig2_id)))
        self.assertFalse(finder.is_att_denovo(RangePair(Range(300, 350), Range(3000, 3050), self.contig2_id)))
