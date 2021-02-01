import os
import unittest
from contextlib import contextmanager
from subprocess import CalledProcessError
from tempfile import NamedTemporaryFile, TemporaryDirectory

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from explore_pipolin.utilities.external_tools import run_prodigal, run_hmmsearch, \
    blastn_against_ref_att, blast_for_repeats, run_aragorn, run_prokka
from explore_pipolin.utilities.external_tools import subprocess_with_retries


class TestExternalRuns(unittest.TestCase):
    def setUp(self) -> None:
        self.nucl_seq = Seq('AGCTAGCAGCGAGGTATGCGATGCGAC')
        self.aa_seq = Seq('MCNIKHVDGGSAEYHVDGGSAEYWD*')

    def test_prodigal_nofile(self):
        with TemporaryDirectory() as tmp_dir:
            with self.assertRaises(FileNotFoundError):
                run_prodigal(os.path.join(tmp_dir, 'genome.fa'), os.path.join(tmp_dir, 'genome.faa'))

    def test_prodigal(self):
        with temp_fasta_file('genome', self.nucl_seq) as genome_file:
            run_prodigal(genome_file.name, NamedTemporaryFile().name)

    def test_hmmsearch_nofile(self):
        with TemporaryDirectory() as tmp_dir:
            with self.assertRaises(CalledProcessError):
                run_hmmsearch(os.path.join(tmp_dir, 'genome.faa'), os.path.join(tmp_dir, 'genome.tbl'))

    def test_hmmsearch(self):
        with temp_fasta_file('protein', self.aa_seq) as proteins_file:
            run_hmmsearch(proteins_file.name, NamedTemporaryFile().name)

    def test_blastn_nofile(self):
        with TemporaryDirectory() as tmp_dir:
            with self.assertRaises(CalledProcessError):
                blastn_against_ref_att(os.path.join(tmp_dir, 'genome.fa'), os.path.join(tmp_dir, 'genome.fmt5'))

    def test_blastn(self):
        with temp_fasta_file('genome', self.nucl_seq) as genome_file:
            blastn_against_ref_att(genome_file.name, NamedTemporaryFile().name)

    def test_blast_for_repeats(self):
        with TemporaryDirectory() as tmp_dir:
            with NamedTemporaryFile(suffix='_0.left', prefix='genome', dir=tmp_dir) as left, \
                    NamedTemporaryFile(suffix='_0.right', prefix='genome', dir=tmp_dir) as right:
                SeqIO.write(SeqRecord(self.nucl_seq, id='genome'), left.name, format='fasta')
                SeqIO.write(SeqRecord(self.nucl_seq, id='genome'), right.name, format='fasta')
                blast_for_repeats('genome', tmp_dir)
            self.assertTrue('genome_0.fmt5' in os.listdir(tmp_dir))

    def test_blast_for_repeats_no_right(self):
        with TemporaryDirectory() as tmp_dir:
            with NamedTemporaryFile(suffix='_0.left', prefix='genome', dir=tmp_dir) as left:
                SeqIO.write(SeqRecord(self.nucl_seq, id='genome'), left.name, format='fasta')
                with self.assertRaises(AssertionError) as context:
                    blast_for_repeats('genome', tmp_dir)
                self.assertTrue('Number of .left files is not equal' in str(context.exception))

    def test_blast_for_repeats_wrong_pair(self):
        with TemporaryDirectory() as tmp_dir:
            with NamedTemporaryFile(suffix='_0.left', prefix='genome', dir=tmp_dir) as left,\
                    NamedTemporaryFile(suffix='_1.right', prefix='genome', dir=tmp_dir) as right:
                SeqIO.write(SeqRecord(self.nucl_seq, id='genome'), left.name, format='fasta')
                SeqIO.write(SeqRecord(self.nucl_seq, id='genome'), right.name, format='fasta')
                with self.assertRaises(AssertionError) as context:
                    blast_for_repeats('genome', tmp_dir)
                self.assertTrue(f'Wrong pair for file ' in str(context.exception))

    def test_aragorn_nofile(self):
        with TemporaryDirectory() as tmp_dir:
            with self.assertRaises(FileNotFoundError):
                run_aragorn(os.path.join(tmp_dir, 'genome.fa'), os.path.join(tmp_dir, 'genome.batch'))

    def test_aragorn(self):
        with temp_fasta_file('genome', self.nucl_seq) as genome_file:
            run_aragorn(genome_file.name, NamedTemporaryFile().name)

    def test_run_prokka_no_file(self):
        with TemporaryDirectory() as tmp_dir:
            input_file = os.path.join(tmp_dir, 'genome.fa')
            with self.assertRaises(CalledProcessError):
                run_prokka(input_file, tmp_dir)

    def test_run_prokka(self):
        with TemporaryDirectory() as tmp_dir:
            with temp_fasta_file('genome', self.nucl_seq) as pipolin_file:
                run_prokka(pipolin_file.name, tmp_dir)


@contextmanager
def temp_fasta_file(seq_id: str, seq: Seq):
    with NamedTemporaryFile() as fasta_file:
        SeqIO.write([SeqRecord(seq, id=seq_id)], fasta_file.name, format='fasta')
        yield fasta_file


@unittest.skip
class ExternalToolsTestCase(unittest.TestCase):
    def test_retries_success(self):
        proc = subprocess_with_retries(['echo', 'HI', '>', '/dev/null'])
        self.assertIsNone(proc.check_returncode())

    def test_retries_fail(self):
        proc = subprocess_with_retries(['cat', 'asdf.asdf'])
        self.assertIsNone(proc)
