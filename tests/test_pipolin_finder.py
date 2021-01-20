import unittest

from explore_pipolin.common import ContigID, Genome, Contig, Feature, AttFeature, Range, Strand, FeatureType, \
    PipolinFragment, Pipolin
from explore_pipolin.tasks_related.find_pipolins import PipolinFinder


class TestPipolinFinder(unittest.TestCase):
    def setUp(self) -> None:
        self.contig1_id = ContigID('contig1')
        self.contig2_id = ContigID('contig2')
        self.genome = Genome('genome', 'genome.fa', contigs=[Contig(self.contig1_id, 10000),
                                                             Contig(self.contig2_id, 5000)])
        self.att1 = AttFeature(Range(50, 100), Strand.FORWARD, self.contig1_id, self.genome, att_id=1)
        self.pipolb = Feature(Range(3000, 5000), Strand.FORWARD, self.contig1_id, self.genome)
        self.att2 = AttFeature(Range(9000, 9050), Strand.FORWARD, self.contig1_id, self.genome, att_id=1)

        self.genome.features.add_features(self.att1, self.att2, feature_type=FeatureType.ATT)
        self.genome.features.add_features(self.pipolb, feature_type=FeatureType.PIPOLB)

        # analysis
        self.finder = PipolinFinder(self.genome)
        self.fragments = self.finder.find_pipolin_fragment_candidates()
        self.pipolins = self.finder.find_pipolin_candidates(self.fragments)

        # expected
        f1 = PipolinFragment(self.pipolb.location, self.contig1_id,
                             features=((self.pipolb, FeatureType.PIPOLB),))
        f2 = PipolinFragment(Range(self.att1.start, self.pipolb.end), self.contig1_id,
                             features=((self.att1, FeatureType.ATT), (self.pipolb, FeatureType.PIPOLB)))
        f3 = PipolinFragment(Range(self.pipolb.start, self.att2.end), self.contig1_id,
                             features=((self.pipolb, FeatureType.PIPOLB), (self.att2, FeatureType.ATT)))
        f4 = PipolinFragment(Range(self.att1.start, self.att2.end), self.contig1_id,
                             features=((self.att1, FeatureType.ATT),
                                       (self.pipolb, FeatureType.PIPOLB),
                                       (self.att2, FeatureType.ATT)))
        self.exp_fragments = [f1, f2, f3, f4]
        self.exp_pipolins = [Pipolin.from_fragments(f1), Pipolin.from_fragments(f2),
                             Pipolin.from_fragments(f3), Pipolin.from_fragments(f4)]

    def test_find_pipolin_fragment_candidates(self):
        self.assertEqual(set(self.exp_fragments), set(self.fragments))

        # additional pipolb on the same contig

        # additional pipolb on a separate contig
        pipolb2 = Feature(Range(500, 1500), Strand.FORWARD, self.contig2_id, self.genome)
        self.genome.features.add_features(pipolb2, feature_type=FeatureType.PIPOLB)

        obs = PipolinFinder(self.genome).find_pipolin_fragment_candidates()

        additional_fragment = PipolinFragment(pipolb2.location, self.contig2_id,
                                              features=((pipolb2, FeatureType.PIPOLB),))

        self.assertEqual(set(self.exp_fragments + [additional_fragment]), set(obs))

    def test_find_pipolin_candidates(self):
        self.assertEqual(set(self.exp_pipolins), set(self.pipolins))
