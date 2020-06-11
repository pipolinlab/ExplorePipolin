import unittest

from explore_pipolin.common import Genome, Feature, Orientation, FeatureType, Contig, GenomeFeatures


class MiscTestCase(unittest.TestCase):
    def test_get_overlapping_feature(self):
        genome = Genome('gen', 'gen.fa')
        genome.contigs.append(Contig('cid', 10))
        features_container: GenomeFeatures = GenomeFeatures(genome)
        trna: Feature = Feature(0, 1, Orientation.FORWARD, 'cid', genome)
        features_container.get_features(FeatureType.TRNA).append(trna)
        att: Feature = Feature(0, 1, Orientation.FORWARD, 'cid', genome)
        self.assertIs(trna, features_container.find_overlapping_feature(att, FeatureType.TRNA))

