from click.testing import CliRunner, Result
import tempfile
import unittest
from contextlib import contextmanager

from explore_pipolin.main import main


class TestWrongInvocation(unittest.TestCase):
    def setUp(self) -> None:
        self.runner = CliRunner()

    def test_identical_genome_file_names(self):
        with temp_genome_file() as genome_file:
            with self.assertLogs() as logs:
                result: Result = self.runner.invoke(main, [genome_file, genome_file])
                assert result.exit_code == 1
            self.assertEqual(len(logs.output), 1)
            self.assertIn('GENOMES should have different names!', logs.output[0])

    def test_out_dirs_conflict(self):
        with temp_genome_file() as genome_file, tempfile.TemporaryDirectory() as tmp:
            with self.assertLogs() as logs:
                result = self.runner.invoke(
                    main, ['--out-dir-prefix', 'my_output', '--out-dir', tmp, genome_file]
                )
                assert result.exit_code == 1
            self.assertEqual(len(logs.output), 1)
            self.assertIn('Options --out-dir-prefix and --out-dir are mutually exclusive!', logs.output[0])


@contextmanager
def temp_genome_file():
    with tempfile.NamedTemporaryFile() as genome_file:
        yield genome_file.name
