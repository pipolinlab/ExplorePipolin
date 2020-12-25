import unittest

from explore_pipolin.common import Orientation, Contig, Genome, Feature, RepeatPair, PipolinFragment, Range
from explore_pipolin.common import define_genome_id


class TestOrientation(unittest.TestCase):
    def test_default_contig_orientation(self):
        self.assertEqual(Contig('foo', 100).orientation, Orientation.FORWARD)

    def test_from_pm_one_encoding(self):
        self.assertEqual(Orientation.from_pm_one_encoding(1), Orientation.FORWARD)
        self.assertEqual(Orientation.from_pm_one_encoding(-1), Orientation.REVERSE)
        with self.assertRaises(AssertionError):
            Orientation.from_pm_one_encoding(0)

    def test_to_pm_one_encoding(self):
        self.assertEqual(Orientation.FORWARD.to_pm_one_encoding(), 1)
        self.assertEqual(Orientation.REVERSE.to_pm_one_encoding(), -1)

    def test_negation(self):
        self.assertEqual(-Orientation.FORWARD, Orientation.REVERSE)
        self.assertEqual(-Orientation.REVERSE, Orientation.FORWARD)


class SetUpGenome(unittest.TestCase):
    def setUp(self) -> None:
        self.short_contig_id = 'foo'
        self.short_contig = Contig(contig_id=self.short_contig_id, contig_length=100)

        self.long_contig_id = 'boo'
        self.long_contig = Contig(contig_id=self.long_contig_id, contig_length=500)

        self.single_contig_genome = Genome(genome_id='bar', genome_file='dir/bar.fa',
                                           contigs=[self.short_contig])

        self.multi_contig_genome = Genome(genome_id='car', genome_file='dir/car.fa',
                                          contigs=[self.short_contig, self.long_contig])

        self.long_contig_feature = Feature(Range(start=123, end=321), strand=Orientation.REVERSE,
                                           contig_id=self.long_contig_id, genome=self.multi_contig_genome)

        self.pipolin = PipolinFragment(contig_id='boo', genome=self.multi_contig_genome,
                                       coords=Range(start=300, end=400))

        self.repeat_f1 = Range(start=10, end=20)
        self.repeat_f2 = Range(start=60, end=70)
        self.repeat_f3 = Range(start=65, end=85)
        self.repeat = RepeatPair(self.repeat_f1, self.repeat_f2,
                                 right_seq='GATTACAATC', left_seq='GATTACAATC',
                                 pipolbs=[self.long_contig_feature])


class TestGenome(SetUpGenome):
    def test_is_single_contig(self):
        self.assertTrue(self.single_contig_genome.is_single_contig())
        self.assertFalse(self.multi_contig_genome.is_single_contig())
        # TODO: implement repeats search for incomplete genomes and delete this!
        self.assertEqual(self.single_contig_genome.get_complete_genome_contig_id(), self.short_contig_id)
        with self.assertRaises(AssertionError):
            self.multi_contig_genome.get_complete_genome_contig_id()

    def test_get_contig_by_id(self):
        self.assertEqual(self.single_contig_genome.get_contig_by_id(contig_id=self.short_contig_id), self.short_contig)
        with self.assertRaises(AssertionError):
            self.single_contig_genome.get_contig_by_id(contig_id=self.long_contig_id)


class TestRange(SetUpGenome):
    def test_start_greater_than_end(self):
        with self.assertRaises(AssertionError):
            Range(start=20, end=10)

    def test_start_less_than_zero(self):
        with self.assertRaises(AssertionError):
            Range(start=-1, end=10)

    def test_shift(self):
        size = 5
        self.assertEqual(self.repeat_f1.shift(size),
                         Range(self.repeat_f1.start + size, self.repeat_f1.end + size))

    def test_inflate(self):
        size = 5
        self.assertEqual(self.repeat_f1.inflate(size),
                         Range(self.repeat_f1.start - size, self.repeat_f1.end + size))

    def test_inflate_beyond_zero(self):
        size = 20
        self.assertEqual(self.repeat_f1.inflate(size).start, 0)

    def test_clamp(self):
        self.assertEqual(self.repeat_f1.clamp(0, self.short_contig.length), self.repeat_f1)

    def test_is_overlapping(self):
        self.assertFalse(self.repeat_f1.is_overlapping(self.repeat_f2))
        self.assertTrue(self.repeat_f2.is_overlapping(self.repeat_f3))


class TestFeature(SetUpGenome):
    def test_feature_contig_property(self):
        self.assertEqual(self.long_contig_feature.contig, self.long_contig)

    def test_feature_end_not_greater_contig_length(self):
        with self.assertRaises(AssertionError):
            Feature(Range(start=123, end=321), strand=Orientation.FORWARD,
                    contig_id=self.short_contig.id, genome=self.multi_contig_genome)


class UtilitiesTestCase(SetUpGenome):
    def test_repeat_seq_not_greater_range(self):
        with self.assertRaises(AssertionError):
            RepeatPair(self.repeat_f1, self.repeat_f2, 'ATCGTAGCTTGACTTCG', 'ATCGTAGCTTGACTTCG', [self.long_contig_feature])

    def test_repeat_left_shift_greater_right_shift(self):
        with self.assertRaises(AssertionError):
            self.repeat.shift(left_shift=250, right_shift=150)

    def test_pipolin_contig_property(self):
        self.assertEqual(self.pipolin.contig, self.long_contig)

    def test_pipolin_start_not_greater_end(self):
        with self.assertRaises(AssertionError):
            PipolinFragment(contig_id='boo', genome=self.multi_contig_genome,
                            coords=Range(start=400, end=300))

    def test_pipolin_end_not_greater_contig_length(self):
        with self.assertRaises(AssertionError):
            PipolinFragment(contig_id='foo', genome=self.multi_contig_genome,
                            coords=Range(start=300, end=400))

    def test_define_genome_id(self):
        self.assertEqual(define_genome_id('my_genome.fa'), 'my_genome')
        with self.assertRaises(AssertionError):
            define_genome_id('thisisverylongfilebasename.fa')


if __name__ == '__main__':
    unittest.main()
