#######################################################################
# Tests for conda.py module
#######################################################################

import unittest
import tempfile
import shutil
import os
from auto_process_ngs.command import Command
from auto_process_ngs.conda import CondaWrapper
from auto_process_ngs.conda import CondaCreateEnvError
from auto_process_ngs.conda import make_conda_env_name
from auto_process_ngs.mock import MockConda

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

# Tests

class TestCondaWrapper(unittest.TestCase):

    def setUp(self):
        # Make a temporary working dir
        self.working_dir = tempfile.mkdtemp(
            suffix='TestCondaWrapper')
        # Save PATH
        self.save_path = os.environ['PATH']

    def tearDown(self):
        # Remove temp dir
        if os.path.exists(self.working_dir):
            shutil.rmtree(self.working_dir)
        # Restore PATH
        os.environ['PATH'] = self.save_path

    def _make_mock_conda(self,create_fails=False,
                         activate_fails=False):
        # Internal: make a mock conda installation
        self.conda_dir = os.path.join(self.working_dir,
                                      "conda")
        self.conda_bin_dir = os.path.join(self.conda_dir,
                                          "bin")
        self.conda_env_dir = os.path.join(self.conda_dir,
                                          "envs")
        # Create mock conda using supplied options
        MockConda.create(self.conda_dir,
                         create_fails=create_fails,
                         activate_fails=activate_fails)
        self.conda = os.path.join(self.conda_bin_dir,"conda")
        # Update PATH (put mock conda first)
        os.environ['PATH'] = self.conda_bin_dir + \
                             os.pathsep + \
                             os.environ['PATH']

    def test_conda_wrapper_conda_version(self):
        """
        CondaWrapper: get conda version
        """
        self._make_mock_conda()
        conda = CondaWrapper(conda=self.conda)
        self.assertEqual(conda.version,"4.10.3")

    def test_conda_wrapper_conda_not_specified(self):
        """
        CondaWrapper: check properties when conda exe not specified
        """
        self._make_mock_conda()
        print(os.environ['PATH'])
        conda = CondaWrapper()
        self.assertEqual(conda.conda,self.conda)
        self.assertTrue(conda.is_installed)
        self.assertEqual(conda.env_dir,self.conda_env_dir)
        self.assertEqual(conda.list_envs,[])

    def test_conda_wrapper_conda_defaults(self):
        """
        CondaWrapper: check properties for default env dir
        """
        self._make_mock_conda()
        conda = CondaWrapper(conda=self.conda)
        self.assertEqual(conda.conda,self.conda)
        self.assertTrue(conda.is_installed)
        self.assertEqual(conda.env_dir,self.conda_env_dir)
        self.assertEqual(conda.list_envs,[])

    def test_conda_wrapper_conda_custom_env_dir(self):
        """
        CondaWrapper: check properties for custom env dir
        """
        self._make_mock_conda()
        custom_env_dir = os.path.join(self.working_dir,
                                      "local_conda_envs")
        os.makedirs(custom_env_dir)
        conda = CondaWrapper(conda=self.conda,
                             env_dir=custom_env_dir)
        self.assertEqual(conda.conda,self.conda)
        self.assertTrue(conda.is_installed)
        self.assertEqual(conda.env_dir,custom_env_dir)
        self.assertEqual(conda.list_envs,[])

    def test_conda_wrapper_null_wrapper(self):
        """
        CondaWrapper: check properties for 'null' wrapper
        """
        save_path = os.environ['PATH']
        os.environ['PATH'] = ''
        conda = CondaWrapper()
        self.assertEqual(conda.conda,None)
        self.assertFalse(conda.is_installed)
        self.assertEqual(conda.env_dir,None)
        self.assertEqual(conda.list_envs,[])
        self.assertEqual(conda.version,None)
        os.environ['PATH'] = save_path

    def test_conda_wrapper_missing_conda_executable(self):
        """
        CondaWrapper: check properties when conda is missing
        """
        conda = CondaWrapper(
            conda="/usr/local/missing/conda/bin/conda")
        self.assertEqual(conda.conda,
                         "/usr/local/missing/conda/bin/conda")
        self.assertFalse(conda.is_installed)
        self.assertEqual(conda.env_dir,
                         "/usr/local/missing/conda/envs")
        self.assertEqual(conda.list_envs,[])
        self.assertEqual(conda.version,None)

    def test_conda_wrapper_list_envs(self):
        """
        CondaWrapper: check listing environments
        """
        self._make_mock_conda()
        # Default env dir
        conda = CondaWrapper(conda=self.conda)
        self.assertEqual(conda.list_envs,[])
        for env in ('fastqc','picard','star'):
            os.makedirs(os.path.join(self.conda_env_dir,env))
        self.assertEqual(conda.list_envs,
                         ['fastqc','picard','star'])
        # Custom env dir
        custom_env_dir = os.path.join(self.working_dir,
                                      "local_conda_envs")
        os.makedirs(custom_env_dir)
        conda = CondaWrapper(conda=self.conda,
                             env_dir=custom_env_dir)
        self.assertEqual(conda.list_envs,[])
        for env in ('fastqc','picard','star'):
            os.makedirs(os.path.join(custom_env_dir,env))
        self.assertEqual(conda.list_envs,
                         ['fastqc','picard','star'])

    def test_conda_wrapper_create_env(self):
        """
        CondaWrapper: check create new environment
        """
        self._make_mock_conda()
        conda = CondaWrapper(conda=self.conda)
        conda.create_env("qc",
                         "fastqc=0.11.3",
                         "fastq-screen=0.14.0",
                         "bowtie=1.2.3")
        env_dir = os.path.join(self.conda_env_dir,"qc")
        self.assertTrue(os.path.exists(env_dir))
        packages = os.path.join(env_dir,"packages.txt")
        self.assertTrue(os.path.exists(packages))
        with open(packages,'rt') as fp:
            contents = fp.read().strip()
            self.assertEqual(contents,
                             "fastqc=0.11.3 fastq-screen=0.14.0 bowtie=1.2.3")

    def test_conda_wrapper_activate_command(self):
        """
        CondaWrapper: check environment activation command
        """
        self._make_mock_conda()
        conda = CondaWrapper(conda=self.conda)
        self.assertEqual(conda.activate_env_cmd("fastqc@0.11.3").command_line,
                         Command(
                             'source',
                             os.path.join(self.conda_bin_dir,'activate'),
                             'fastqc@0.11.3').command_line)

    def test_conda_wrapper_verify_environment(self):
        """
        CondaWrapper: verify environment can be activated
        """
        self._make_mock_conda()
        os.makedirs(os.path.join(self.conda_env_dir,'fastqc@0.11.3'))
        conda = CondaWrapper(conda=self.conda)
        self.assertTrue(conda.verify_env("fastqc@0.11.3"))

    def test_conda_wrapper_verify_fails_for_bad_environment(self):
        """
        CondaWrapper: verify fails for 'bad' environment
        """
        self._make_mock_conda(activate_fails=True)
        os.makedirs(os.path.join(self.conda_env_dir,'fastqc@0.11.3'))
        conda = CondaWrapper(conda=self.conda)
        self.assertFalse(conda.verify_env("fastqc@0.11.3"))

    def test_conda_wrapper_create_env_raise_exception_on_error(self):
        """
        CondaWrapper: raise exception on error when creating environment
        """
        self._make_mock_conda(create_fails=True)
        conda = CondaWrapper(conda=self.conda)
        self.assertRaises(CondaCreateEnvError,
                          conda.create_env,
                          "qc",
                          "fastqc=0.11.3",
                          "fastq-screen=0.14.0",
                          "bowtie=1.2.3")

class TestMakeCondaEnvName(unittest.TestCase):

    def test_make_conda_env_name(self):
        """
        make_conda_env_name: returns consistent environment name
        """
        self.assertEqual(make_conda_env_name("fastq-screen=0.14.0",
                                             "bowtie=1.2.3"),
                         "bowtie@1.2.3+fastq-screen@0.14.0")
        self.assertEqual(make_conda_env_name("bowtie=1.2.3",
                                             "fastq-screen=0.14.0"),
                         "bowtie@1.2.3+fastq-screen@0.14.0")
