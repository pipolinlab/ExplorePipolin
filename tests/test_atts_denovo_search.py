import unittest
from typing import List

from explore_pipolin.common import Contig, Genome, Feature, Range, Strand, FeatureType, RangePair
from explore_pipolin.tasks_related.misc import get_ranges_around_pipolbs


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

        self.genome.features.add_feature(pipolb1, FeatureType.PIPOLB)
        self.genome.features.add_feature(pipolb2, FeatureType.PIPOLB)
        self.genome.features.add_feature(pipolb3, FeatureType.PIPOLB)

    def test_get_ranges_around_pipolbs(self):
        range_pair1 = RangePair(left=Range(0, 100), right=Range(1100, 1100 + 100000), contig_id=self.contig1_id)
        range_pair2 = RangePair(left=Range(0, 3000), right=Range(4000, 4000 + 100000), contig_id=self.contig1_id)
        range_pair3 = RangePair(left=Range(0, 100), right=Range(4000, 4000 + 100000), contig_id=self.contig1_id)
        range_pair4 = RangePair(left=Range(0, 500), right=Range(1500, 10000), contig_id=self.contig2_id)

        obt: List[RangePair] = get_ranges_around_pipolbs(self.genome)
        self.assertCountEqual([range_pair1, range_pair2, range_pair3, range_pair4], obt)
