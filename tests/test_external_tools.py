import unittest

from explore_pipolin.utilities.external_tools import subprocess_with_retries


class ExternalToolsTestCase(unittest.TestCase):
    def test_retries_success(self):
        proc = subprocess_with_retries(['echo', 'HI', '>', '/dev/null'])
        self.assertIsNone(proc.check_returncode())

    def test_retries_fail(self):
        proc = subprocess_with_retries(['cat', 'asdf.asdf'])
        self.assertIsNone(proc)
