from tests.test_common import SetUpGenome


class TestReconstructPipolins(SetUpGenome):
    def test_filter_redundant(self):
        p1_features = {self.long_contig_feature2_pipolb, self.long_contig_feature1_pipolb}
        p2_features = {self.long_contig_feature1_pipolb}
        self.assertTrue(p2_features.issubset(p1_features))
