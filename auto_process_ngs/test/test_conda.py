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

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

# Helpers

class _Mock(object):
    """
    Helper class for mocking executables for testing
    """

    @staticmethod
    def conda(bin_dir):
        """
        Make a mock conda executable for testing

        Arguments:
          bin_dir (str): path to directory to put mock
            'conda' executable into (must already exist)

        Returns:
           String: path to mock 'conda' executable.
        """
        conda_bin_dir = os.path.abspath(bin_dir)
        conda_ = os.path.join(bin_dir,"conda")
        with open(conda_,'wt') as fp:
            fp.write("""#!/bin/bash
if [ "$1" == "--version" ] ; then
   echo "conda 4.10.3"
   exit 0
elif [ "$1" == "create" ] ; then
   YES=
   PREFIX=
   PACKAGES=
   while [ ! -z "$2" ] ; do
     case "$2" in
       -n)
         shift
         PREFIX=$(dirname $(dirname $0))/envs/${2}
         ;;
       --prefix)
         shift
         PREFIX=$2
         ;;
       -y)
         YES=yes
         ;;
       -c)
         shift
         ;;
       --override-channels)
         ;;
       *)
         PACKAGES="$PACKAGES $2"
         ;;
     esac
     shift
   done
fi
if [ -z "$YES" ] ; then
   echo "Need to supply -y option"
   exit 1
fi
if [ -z "$PREFIX" ] ; then
   echo "Need to supply either -n or --prefix"
   exit 1
fi
# Make directory for new environment
mkdir -p $PREFIX
# Write package list to a 'packages.txt' file
echo $PACKAGES >${PREFIX}/packages.txt
# Make an executable script for each package name
for pkg in $PACKAGES ; do
   name=$(echo $pkg | cut -f1 -d=)
   cat >${PREFIX}/${name} <<EOF
#!/bin/bash
echo \$1
exit 0
EOF
   chmod +x ${PREFIX}/${name}
done
""")
            os.chmod(conda_,0o755)
        activate_ = os.path.join(bin_dir,"activate")
        with open(activate_,'wt') as fp:
            fp.write("""#!/bin/bash
export PATH=$PATH:${1}
""")
            os.chmod(activate_,0o755)
        # Return path to mock conda
        return conda_

    @staticmethod
    def conda_with_failing_create(bin_dir):
        """
        Make a mock conda executable with failing 'create' command

        Arguments:
          bin_dir (str): path to directory to put mock
            'conda' executable into (must already exist)

        Returns:
           String: path to mock 'conda' executable.
        """
        conda_ = os.path.join(bin_dir,"conda")
        with open(conda_,'wt') as fp:
            fp.write("""#!/bin/bash
if [ "$1" == "--version" ] ; then
   echo "conda 4.10.3"
   exit 0
elif [ "$1" == "create" ] ; then
   echo "!!!! Failed to create environment !!!!"
   exit 1
fi
""")
            os.chmod(conda_,0o755)
        activate_ = os.path.join(bin_dir,"activate")
        with open(activate_,'wt') as fp:
            fp.write("""#!/bin/bash
export PATH=$PATH:${1}
""")
            os.chmod(activate_,0o755)
        # Return path to mock conda
        return conda_

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

    def _make_mock_conda(self,mock_conda_func):
        # Internal: make a mock conda executable
        # for use in tests
        self.conda_dir = os.path.join(self.working_dir,
                                      "conda")
        self.conda_bin_dir = os.path.join(self.conda_dir,
                                          "bin")
        self.conda_env_dir = os.path.join(self.conda_dir,
                                          "envs")
        for d in (self.conda_dir,
                  self.conda_bin_dir,
                  self.conda_env_dir):
            os.makedirs(d)
        # Create mock conda using supplied function
        self.conda = mock_conda_func(self.conda_bin_dir)
        # Update PATH (put mock conda first)
        os.environ['PATH'] = self.conda_bin_dir + \
                             os.pathsep + \
                             os.environ['PATH']

    def test_conda_wrapper_conda_version(self):
        """
        CondaWrapper: get conda version
        """
        self._make_mock_conda(_Mock.conda)
        conda = CondaWrapper(conda=self.conda)
        self.assertEqual(conda.version,"4.10.3")

    def test_conda_wrapper_conda_not_specified(self):
        """
        CondaWrapper: check properties when conda exe not specified
        """
        self._make_mock_conda(_Mock.conda)
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
        self._make_mock_conda(_Mock.conda)
        conda = CondaWrapper(conda=self.conda)
        self.assertEqual(conda.conda,self.conda)
        self.assertTrue(conda.is_installed)
        self.assertEqual(conda.env_dir,self.conda_env_dir)
        self.assertEqual(conda.list_envs,[])

    def test_conda_wrapper_conda_custom_env_dir(self):
        """
        CondaWrapper: check properties for custom env dir
        """
        self._make_mock_conda(_Mock.conda)
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
        self._make_mock_conda(_Mock.conda)
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
        self._make_mock_conda(_Mock.conda)
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
        self._make_mock_conda(_Mock.conda)
        conda = CondaWrapper(conda=self.conda)
        self.assertEqual(conda.activate_env_cmd("fastqc@0.11.3").command_line,
                         Command(
                             'source',
                             os.path.join(self.conda_bin_dir,'activate'),
                             'fastqc@0.11.3').command_line)

    def test_conda_wrapper_create_env_raise_exception_on_error(self):
        """
        CondaWrapper: raise exception on error when creating environment
        """
        self._make_mock_conda(_Mock.conda_with_failing_create)
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
