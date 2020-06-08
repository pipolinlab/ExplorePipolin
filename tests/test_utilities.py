import unittest

from explore_pipolin.common import Orientation, Contig


class UtilitiesTestCase(unittest.TestCase):
    def test_orientation(self):
        self.assertEqual(-Orientation.FORWARD, Orientation.REVERSE)
        self.assertEqual(-Orientation.REVERSE, Orientation.FORWARD)

    def test_default_contig_orientation(self):
        self.assertEqual(Contig('foo', 10).contig_orientation, Orientation.FORWARD)

    def test_blast_orientation(self):
        self.assertEqual(Orientation.orientation_from_blast(1), Orientation.FORWARD)
        self.assertEqual(Orientation.orientation_from_blast(-1), Orientation.REVERSE)


if __name__ == '__main__':
    unittest.main()
