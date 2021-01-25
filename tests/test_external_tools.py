import unittest
from contextlib import contextmanager
from tempfile import NamedTemporaryFile, TemporaryDirectory

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from prefect import Flow

from explore_pipolin.tasks.find_pipolbs import find_pipolbs
from explore_pipolin.common import Genome, Contig, ContigID
from explore_pipolin.utilities.external_tools import subprocess_with_retries, ExternalTools
from explore_pipolin.utilities.logging import set_logging_dir


@unittest.skip('Will be probably replaced with the EMBL-EBI API!')
class ExternalToolsTestCase(unittest.TestCase):
    def test_retries_success(self):
        proc = subprocess_with_retries(['echo', 'HI', '>', '/dev/null'])
        self.assertIsNone(proc.check_returncode())

    def test_retries_fail(self):
        proc = subprocess_with_retries(['cat', 'asdf.asdf'])
        self.assertIsNone(proc)


@contextmanager
def temp_genome(genome_id: str, seq: Seq):
    with NamedTemporaryFile() as genome_file:
        SeqIO.write([SeqRecord(seq, id=genome_id)], genome_file.name, format='fasta')
        yield Genome(genome_id, genome_file.name, 'output/' + genome_id,
                     [Contig(contig_id=ContigID('pipolb'), contig_length=len(seq))])


class MockExternalTools(ExternalTools):
    def find_cdss(self, genome_file: str, output_file: str):
        with open(output_file, 'w'):
            SeqIO.write([SeqRecord(Seq('MI*'), description='# 8 # 16 # 1 # ID=1', id='pipolb_1')],
                        output_file, format='fasta')

    def find_pipolbs(self, cds_file: str, output_file: str):
        with open(output_file, 'w') as ouf:
            print('pipolb_1 - pipolb_expanded_definitive - 0 1284.5 0.0 0 1284.4 0.0 1.0 1 0 0 1 1 1 1 '
                  '# 8 # 16 # 1 # ID=1', file=ouf)

    def blastn_against_ref_att(self, genome_file: str, output_file: str):
        pass


class TestingGenome(unittest.TestCase):
    @unittest.skip
    def test_genome(self):
        with temp_genome('G1', Seq('GCATCGGATGATCCGAGGCATCGGGGCATG')) as genome, TemporaryDirectory() as results_dir:
            set_logging_dir(results_dir)
            with Flow('MAIN') as flow:
                find_pipolbs(genome=genome)
            state = flow.run()
            assert state.is_successful()
