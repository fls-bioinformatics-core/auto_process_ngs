#######################################################################
# Tests for config.py module
#######################################################################

import unittest
from config import fetch_runner

class TestFetchRunnerFunction(unittest.TestCase):
    """Tests for the fetch_runner function
    """

    def test_fetch_simple_job_runner(self):
        """fetch_runner returns a SimpleJobRunner
        """
        runner = fetch_runner("SimpleJobRunner")
        self.assertTrue(isinstance(runner,SimpleJobRunner))

    def test_fetch_ge_job_runner(self):
        """fetch_runner returns a GEJobRunner
        """
        runner = fetch_runner("GEJobRunner")
        self.assertTrue(isinstance(runner,GEJobRunner))

    def test_fetch_ge_job_runner_with_extra_args(self):
        """fetch_runner returns a GEJobRunner with additional arguments
        """
        runner = fetch_runner("GEJobRunner(-j y)")
        self.assertTrue(isinstance(runner,GEJobRunner))
        self.assertEqual(runner.ge_extra_args,['-j','y'])

    def test_fetch_bad_runner_raises_exception(self):
        """fetch_runner raises exception for unknown runner
        """
        self.assertRaises(Exception,fetch_runner,"SimpleRunner")
