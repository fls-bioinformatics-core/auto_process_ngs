#######################################################################
# Tests for bcl2fastq_utils.py module
#######################################################################

import unittest
import tempfile
import os
import shutil
from auto_process_ngs.bcl2fastq_utils import available_bcl2fastq_versions
from auto_process_ngs.bcl2fastq_utils import get_nmismatches

class MockBcl2fastq:
    """
    Class for setting up fake bcl2fastq conversion installs

    Usage:

    Create a new MockBcl2fastq instance:

    >>> m = MockBcl2fastq()

    Set up mock bcl2fastq 1.8.3 and 1.8.4 installations:

    >>> m.bcl2fastq_183()
    >>> m.bcl2fastq_184()

    Get a list of the executables that are defined:

    >>> exes = m.exes

    Append the paths for these installations to the PATH
    environment variable:

    >>> m.set_path(append=True)

    Remove the base directory and contents:

    >>> del(m)

    """
    def __init__(self):
        self.dirn = tempfile.mkdtemp(suffix='mockbcl2fastq')
        self._exes = []

    def __del__(self):
        print "Deleting %s" % self.dirn
        shutil.rmtree(self.dirn)

    def _makedirs(self,*args):
        """
        Internal: create a subdirectory under the base dir

        Creates a new directory using the supplied
        arguments to specify the path relative to the
        base directory, e.g.

        >>> self._makedirs('package','lib')

        will create directories 'package/lib/' under the
        base directory.

        """
        os.makedirs(os.path.join(self.dirn,*args))

    def _make_exe(self,*args,**kws):
        """
        Internal: create an executable file

        Creates a new executable file using the supplied
        arguments to specify the path relative to the
        base directory, e.g.

        >>> self._make_exe('bin','test.sh')

        will create an executable file called 'test.sh'
        in the 'bin' subdir of the base directory.

        The file itself will be empty by default; content
        can be supplied via the 'content' keyword, e.g.

        >>> self._make_exe('bin','test.sh',
        ...                content="#!/bin/bash\necho hello")

        """
        if 'content' in kws:
            content = kws['content']
        else:
            content = ""
        print "Content = %s" % content
        with open(os.path.join(self.dirn,*args),'w') as fp:
            fp.write(content)
        os.chmod(os.path.join(self.dirn,*args),0775)
        self._exes.append(os.path.join(self.dirn,*args))

    def set_path(self,append=False):
        """
        Set the PATH env variable to the paths of the exes

        Arguments:
          append (boolean): if True then append to the
            existing PATH, otherwise reset the PATH to
            only contain the defined paths.

        """
        new_path = os.pathsep.join(self.paths)
        if append:
            new_path = os.environ['PATH'] + os.pathsep + new_path
        os.environ['PATH'] = new_path
        print "PATH = %s" % os.environ['PATH']

    @property
    def paths(self):
        """
        Return paths defined by executables

        """
        return [os.path.dirname(x) for x in self._exes]

    @property
    def exes(self):
        """
        Return executables (full paths)

        """
        return [x for x in self._exes]

    def casava_182(self):
        """
        Add a mock CASAVA 1.8.2 installation

        """
        # CASAVA 1.8.2
        self._makedirs('casava','1.8.2','bin')
        self._makedirs('casava','1.8.2','etc','CASAVA-1.8.2')
        self._make_exe('casava','1.8.2','bin','configureBclToFastq.pl')

    def bcl2fastq_183(self):
        """
        Add a mock bcl2fastq 1.8.3 installation

        """
        # bcl2fastq 1.8.3
        self._makedirs('bcl2fastq','1.8.3','bin')
        self._makedirs('bcl2fastq','1.8.3','etc','bcl2fastq-1.8.3')
        self._make_exe('bcl2fastq','1.8.3','bin','configureBclToFastq.pl')

    def bcl2fastq_184(self):
        """
        Add a mock bcl2fastq 1.8.4 installation

        """
        # bcl2fastq 1.8.4
        self._makedirs('bcl2fastq','1.8.4','bin')
        self._makedirs('bcl2fastq','1.8.4','etc','bcl2fastq-1.8.4')
        self._make_exe('bcl2fastq','1.8.4','bin','configureBclToFastq.pl')

    def bcl2fastq_217(self):
        """
        Add a mock bcl2fastq 2.17 installation

        """
        # bcl2fastq 2.17.1.14
        self._makedirs('bcl2fastq','2.17.1.14','bin')
        self._makedirs('bcl2fastq','2.17.1.14','etc','bcl2fastq-2.17.1.14')
        self._make_exe('bcl2fastq','2.17.1.14','bin','bcl2fastq',
                       content="#!/bin/bash\nif [ \"$1\" == \"--version\" ] ; then cat >&2 <<EOF\nBCL to FASTQ file converter\nbcl2fastq v2.17.1.14\nCopyright (c) 2007-2015 Illumina, Inc.\n\n2015-12-17 14:08:00 [7fa113f3f780] Command-line invocation: bcl2fastq --version\nEOF\nfi")

class TestAvailableBcl2fastqVersions(unittest.TestCase):
    """
    Tests for the available_bcl2fastq_versions function

    """
    def setUp(self):
        self.mockbcl2fastq = MockBcl2fastq()
        self.original_path = os.environ['PATH']

    def tearDown(self):
        os.environ['PATH'] = self.original_path
        del(self.mockbcl2fastq)

    def test_no_versions(self):
        """
        available_bcl2fastq_versions: no conversion software present

        """
        os.environ['PATH'] = ''
        self.assertEqual(available_bcl2fastq_versions(),[])

    def test_casava_182_only(self):
        """
        available_bcl2fastq_versions: only CASAVA 1.8.2 present

        """
        self.mockbcl2fastq.casava_182()
        self.mockbcl2fastq.set_path()
        self.assertEqual(available_bcl2fastq_versions(),
                         self.mockbcl2fastq.exes)

    def test_bcl2fastq_183_only(self):
        """
        available_bcl2fastq_versions: only bcl2fastq 1.8.3 present
        
        """
        self.mockbcl2fastq.bcl2fastq_183()
        self.mockbcl2fastq.set_path()
        self.assertEqual(available_bcl2fastq_versions(),
                         self.mockbcl2fastq.exes)

    def test_bcl2fastq_184_only(self):
        """
        available_bcl2fastq_versions: only bcl2fastq 1.8.4 present
        
        """
        self.mockbcl2fastq.bcl2fastq_184()
        self.mockbcl2fastq.set_path()
        self.assertEqual(available_bcl2fastq_versions(),
                         self.mockbcl2fastq.exes)

    def test_bcl2fastq_217_only(self):
        """
        available_bcl2fastq_versions: only bcl2fastq 2.17 available
        
        """
        self.mockbcl2fastq.bcl2fastq_217()
        self.mockbcl2fastq.set_path()
        self.assertEqual(available_bcl2fastq_versions(),
                         self.mockbcl2fastq.exes)

    def test_bcl2fastq_184_and_bcl2fastq_217(self):
        """
        available_bcl2fastq_versions: both bcl2fastq 1.8.4 and 2.17 present
        
        """
        self.mockbcl2fastq.bcl2fastq_184()
        self.mockbcl2fastq.bcl2fastq_217()
        self.mockbcl2fastq.set_path()
        self.assertEqual(available_bcl2fastq_versions(),
                         self.mockbcl2fastq.exes)

    def test_bcl2fastq_217_and_bcl2fastq_184(self):
        """
        available_bcl2fastq_versions: both bcl2fastq 2.17 and 1.8.4 present
        
        """
        self.mockbcl2fastq.bcl2fastq_217()
        self.mockbcl2fastq.bcl2fastq_184()
        self.mockbcl2fastq.set_path()
        self.assertEqual(available_bcl2fastq_versions(),
                         self.mockbcl2fastq.exes)

class TestGetNmismatches(unittest.TestCase):
    """Tests for the get_nmismatches function

    """
    def test_n_mismatches(self):
        self.assertEqual(get_nmismatches('y50'),0)
        self.assertEqual(get_nmismatches('y50,I4'),0)
        self.assertEqual(get_nmismatches('y50,I6'),1)
        self.assertEqual(get_nmismatches('y101,I6,y101'),1)
        self.assertEqual(get_nmismatches('y250,I8,I8,y250'),1)
        self.assertEqual(get_nmismatches('y250,I6nn,I6nn,y250'),1)
        self.assertEqual(get_nmismatches('y250,I6n2,I6n2,y250'),1)
        self.assertEqual(get_nmismatches('y250,I16'),1)

