#######################################################################
# Tests for bcl2fastq_utils.py module
#######################################################################

import unittest
import tempfile
import os
import shutil
from auto_process_ngs.bcl2fastq_utils import *

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

    Prepend the paths for these installations to the PATH
    environment variable:

    >>> m.set_path(prepend=True)

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

    def set_path(self,prepend=True):
        """
        Set the PATH env variable to the paths of the exes

        Arguments:
          prepend (boolean): if True then prepend to the
            existing PATH, otherwise reset the PATH to
            only contain the defined paths.

        """
        new_path = os.pathsep.join(self.paths)
        if prepend:
            new_path = new_path + os.pathsep + os.environ['PATH']
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

    def casava_no_version(self):
        """
        Add a mock CASAVA installation with no discernible version

        """
        self._makedirs('casava','bin')
        self._make_exe('casava','bin','configureBclToFastq.pl')

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

    def test_require_bcl2fastq_version_eq_184(self):
        """
        available_bcl2fastq_versions: require version == '1.8.4'
        """
        self.mockbcl2fastq.bcl2fastq_217()
        self.mockbcl2fastq.bcl2fastq_184()
        self.mockbcl2fastq.set_path()
        # Only second executable on PATH will be returned (1.8.4)
        self.assertEqual(available_bcl2fastq_versions('==1.8.4'),
                         self.mockbcl2fastq.exes[1:])

    def test_require_bcl2fastq_version_ge_184(self):
        """
        available_bcl2fastq_versions: require version >= '1.8.4'
        """
        self.mockbcl2fastq.bcl2fastq_184()
        self.mockbcl2fastq.bcl2fastq_217()
        self.mockbcl2fastq.set_path()
        # Executables returned order by version (2.17 then 1.8.4)
        # (i.e. reverse of order on PATH)
        self.assertEqual(available_bcl2fastq_versions('>=1.8.4'),
                         self.mockbcl2fastq.exes[::-1])

    def test_require_bcl2fastq_version_gt_184(self):
        """
        available_bcl2fastq_versions: require version > '1.8.4'
        """
        self.mockbcl2fastq.bcl2fastq_217()
        self.mockbcl2fastq.bcl2fastq_184()
        self.mockbcl2fastq.set_path()
        # Only first executable will be returned (2.17)
        self.assertEqual(available_bcl2fastq_versions('>1.8.4'),
                         self.mockbcl2fastq.exes[0:1])

    def test_require_bcl2fastq_version_lt_217(self):
        """
        available_bcl2fastq_versions: require version < '2.17'
        """
        self.mockbcl2fastq.bcl2fastq_217()
        self.mockbcl2fastq.bcl2fastq_184()
        self.mockbcl2fastq.set_path()
        # Only second executable on PATH will be returned (1.8.4)
        self.assertEqual(available_bcl2fastq_versions('<2.17'),
                         self.mockbcl2fastq.exes[1:])

    def test_require_bcl2fastq_version_le_217(self):
        """
        available_bcl2fastq_versions: require version <= '2.17.1.14'
        """
        self.mockbcl2fastq.bcl2fastq_217()
        self.mockbcl2fastq.bcl2fastq_184()
        self.mockbcl2fastq.set_path()
        # Executables returned order by version (2.17 then 1.8.4)
        # i.e. same as order on PATH
        self.assertEqual(available_bcl2fastq_versions('<=2.17.1.14'),
                         self.mockbcl2fastq.exes)

    def test_require_bcl2fastq_version_gt_182_lt_184(self):
        """
        available_bcl2fastq_versions: require version > '1.8.2, < '1.8.4''
        """
        self.mockbcl2fastq.bcl2fastq_217()
        self.mockbcl2fastq.bcl2fastq_184()
        self.mockbcl2fastq.bcl2fastq_183()
        self.mockbcl2fastq.casava_182()
        self.mockbcl2fastq.set_path()
        # Only 3rd executable will be returned (2.17)
        self.assertEqual(available_bcl2fastq_versions('>1.8.2,<1.8.4'),
                         self.mockbcl2fastq.exes[2:3])

class TestBclToFastqInfo(unittest.TestCase):
    """
    Tests for the bcl_to_fastq_info function

    """
    def setUp(self):
        # Make some fake directories for different
        # software versions
        self.mockbcl2fastq = MockBcl2fastq()
        self.original_path = os.environ['PATH']

    def tearDown(self):
        # Remove the temporary test directory
        os.environ['PATH'] = self.original_path
        del(self.mockbcl2fastq)

    def test_casava_1_8_2(self):
        """
        Collect info for CASAVA 1.8.2

        """
        # CASAVA 1.8.2
        self.mockbcl2fastq.casava_182()
        self.mockbcl2fastq.set_path()
        exe = self.mockbcl2fastq.exes[0]
        path,name,version = bcl_to_fastq_info()
        self.assertEqual(path,exe)
        self.assertEqual(name,'CASAVA')
        self.assertEqual(version,'1.8.2')

    def test_bcl2fastq_1_8_3(self):
        """
        Collect info for bcl2fastq 1.8.3

        """
        # bcl2fastq 1.8.3
        self.mockbcl2fastq.bcl2fastq_183()
        self.mockbcl2fastq.set_path()
        exe = self.mockbcl2fastq.exes[0]
        path,name,version = bcl_to_fastq_info()
        self.assertEqual(path,exe)
        self.assertEqual(name,'bcl2fastq')
        self.assertEqual(version,'1.8.3')

    def test_bcl2fastq_1_8_4(self):
        """
        Collect info for bcl2fastq 1.8.4

        """
        # bcl2fastq 1.8.4
        self.mockbcl2fastq.bcl2fastq_184()
        self.mockbcl2fastq.set_path()
        exe = self.mockbcl2fastq.exes[0]
        path,name,version = bcl_to_fastq_info()
        self.assertEqual(path,exe)
        self.assertEqual(name,'bcl2fastq')
        self.assertEqual(version,'1.8.4')

    def test_bcl2fastq_2_17_1_14(self):
        """
        Collect info for bcl2fastq 2.17.1.14

        """
        # bcl2fastq 2.17.1.14
        self.mockbcl2fastq.bcl2fastq_217()
        self.mockbcl2fastq.set_path()
        exe = self.mockbcl2fastq.exes[0]
        path,name,version = bcl_to_fastq_info()
        self.assertEqual(path,exe)
        self.assertEqual(name,'bcl2fastq')
        self.assertEqual(version,'2.17.1.14')

    def test_bcl2fastq_1_8_4_supplied_exe(self):
        """
        Collect info for bcl2fastq 1.8.4 when executable is supplied

        """
        # bcl2fastq 1.8.4
        self.mockbcl2fastq.bcl2fastq_183()
        self.mockbcl2fastq.bcl2fastq_184()
        self.mockbcl2fastq.set_path()
        bcl2fastq_184 = self.mockbcl2fastq.exes[1]
        path,name,version = bcl_to_fastq_info(bcl2fastq_184)
        self.assertEqual(path,bcl2fastq_184)
        self.assertEqual(name,'bcl2fastq')
        self.assertEqual(version,'1.8.4')

    def test_bcl2fastq_2_17_1_14_supplied_exe(self):
        """
        Collect info for bcl2fastq 2.17.1.14 when executable is supplied

        """
        # bcl2fastq 2.17.1.14
        self.mockbcl2fastq.bcl2fastq_217()
        self.mockbcl2fastq.set_path()
        bcl2fastq_217 = self.mockbcl2fastq.exes[0]
        path,name,version = bcl_to_fastq_info(bcl2fastq_217)
        self.assertEqual(path,bcl2fastq_217)
        self.assertEqual(name,'bcl2fastq')
        self.assertEqual(version,'2.17.1.14')

    def test_configurebcltofastq_not_found(self):
        """
        Collect info when configureBclToFastq.pl is not found

        """
        path,name,version = bcl_to_fastq_info()
        self.assertEqual(path,None)
        self.assertEqual(name,None)
        self.assertEqual(version,None)

    def test_configurebcltofastq_no_package(self):
        """
        Collect info when there is no package info for configureBclToFastq.pl

        """
        # No version
        self.mockbcl2fastq.casava_no_version()
        self.mockbcl2fastq.set_path()
        exe = self.mockbcl2fastq.exes[0]
        path,name,version = bcl_to_fastq_info()
        self.assertEqual(path,exe)
        self.assertEqual(name,None)
        self.assertEqual(version,None)

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

