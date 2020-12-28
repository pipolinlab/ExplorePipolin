import unittest
import tempfile

from explore_pipolin.utilities.io import create_seqio_records_dict


class TestIO(unittest.TestCase):
    def setUp(self):
        self.contig1 = 'CONTIG1'
        self.contig2 = 'CONTIG2'

    def test_create_seqio_records_genbank(self):
        with tempfile.NamedTemporaryFile('w') as inf:
            gb = (
                'LOCUS       {id}               100 bp    DNA     linear       18-DEC-2020\n'
                'FEATURES             Location/Qualifiers\n'
                '     source          1..100\n'
                'ORIGIN\n'
                '        1 tttgagtggc tttgccggtg attaaaaatt aaggagggtg taacgacaag ttgcaggcac\n'
                '       61 aaaaaaacca cccgaaggtg gtttcacgac actgcttatt\n'
                '//')
            print(gb.format(id=self.contig1), file=inf)
            print(gb.format(id=self.contig2), file=inf)
            inf.flush()
            records = create_seqio_records_dict(file=inf.name, file_format='genbank')
        self.assertEqual(records[self.contig1].id, self.contig1)
        self.assertEqual(records[self.contig2].id, self.contig2)
        self.assertEqual(len(records[self.contig1].seq), 100)

    def test_create_seqio_records_fasta(self):
        with tempfile.NamedTemporaryFile('w') as inf:
            fa = ('>{id} Some information here'
                  'CAGGTACCGGATTCGGAGCTGGATTCGTGGATGGAAAGCCGGATTTATCCGGCGATGACT'
                  'GCGATCCCGGCACTGGCAGACCTGATTACCACGATGGTTACGCAGGGCTATGA')
            print(fa.format(id=self.contig1), file=inf)
            print(fa.format(id=self.contig2), file=inf)
            inf.flush()
            records = create_seqio_records_dict(file=inf.name, file_format='fasta')
        self.assertEqual(records[self.contig1].id, self.contig1)
        self.assertEqual(records[self.contig2].id, self.contig2)
