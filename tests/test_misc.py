import unittest

from explore_pipolin.common import Genome, Feature, Orientation, FeatureType, Contig
from explore_pipolin.utilities.misc import GQuery


class MiscTestCase(unittest.TestCase):
    def test_get_overlapping_feature(self):
        genome = Genome('gen', 'gen.fa')
        genome.contigs.append(Contig('cid', 10))
        gquery: GQuery = GQuery(genome)
        trna: Feature = Feature(0, 1, Orientation.FORWARD, 'cid', genome)
        gquery.trnas.append(trna)
        att: Feature = Feature(0, 1, Orientation.FORWARD, 'cid', genome)
        self.assertIs(trna, gquery.find_overlapping_feature(att, FeatureType.TRNA))

