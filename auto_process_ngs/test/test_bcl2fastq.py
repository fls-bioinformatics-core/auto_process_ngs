#######################################################################
# Tests for bcl2fastq_utils.py module
#######################################################################

import unittest
import tempfile
import os
import shutil
from bcftbx.mock import RunInfoXml
from auto_process_ngs.settings import Settings
from auto_process_ngs.bcl2fastq_utils import *

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

class MockBcl2fastq(object):
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
        print("Deleting %s" % self.dirn)
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
        print("Content = %s" % content)
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
        print("PATH = %s" % os.environ['PATH'])

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

    def bcl2fastq_220(self):
        """
        Add a mock bcl2fastq 2.20 installation

        """
        # bcl2fastq 2.20.0.422
        self._makedirs('bcl2fastq','2.20.0.422','bin')
        self._makedirs('bcl2fastq','2.20.0.422','etc','bcl2fastq-2.20.0.422')
        self._make_exe('bcl2fastq','2.20.0.422','bin','bcl2fastq',
                       content="#!/bin/bash\nif [ \"$1\" == \"--version\" ] ; then cat >&2 <<EOF\nBCL to FASTQ file converter\nbcl2fastq v2.20.0.422\nCopyright (c) 2007-2017 Illumina, Inc.\n\nEOF\nfi")

class TestGetSequencerPlatform(unittest.TestCase):
    """
    Tests for the get_sequencer_platform function
    """
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestAutoProcessRunQc')

    def tearDown(self):
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.dirn)

    def test_get_sequencer_platform(self):
        """
        get_sequencer_platform: use run directory with known sequencers
        """
        self.assertEqual(get_sequencer_platform(
            "/mnt/data/120919_SN7001250_0035_BC133VACXX"),"hiseq")
        self.assertEqual(get_sequencer_platform(
            "/mnt/data/121210_M00879_0001_000000000-A2Y1L"),"miseq")
        self.assertEqual(get_sequencer_platform(
            "/mnt/data/120518_ILLUMINA-73D9FA_00002_FC"),"illumina-ga2x")
        self.assertEqual(get_sequencer_platform(
            "/mnt/data/151216_NB500968_0008_AH5CFGAFXX"),"nextseq")

    def test_get_sequencer_platform_from_settings(self):
        """
        get_sequencer_platform: use run directory and settings
        """
        settings_ini = os.path.join(self.dirn,"settings.ini")
        with open(settings_ini,'w') as s:
            s.write("""[sequencers]
SN7001251 = hiseq
M00880 = miseq
FS10000171 = iseq
""")
        settings = Settings(settings_ini)
        print(settings.sequencers)
        self.assertEqual(get_sequencer_platform(
            "/mnt/data/120919_SN7001251_0035_BC133VACXX",
            settings=settings),"hiseq")
        self.assertEqual(get_sequencer_platform(
            "/mnt/data/121210_M00880_0001_000000000-A2Y1L",
            settings=settings),"miseq")
        self.assertEqual(get_sequencer_platform(
            "/mnt/data/20180829_FS10000171_3_BNT40323-1530",
            settings=settings),"iseq")

    def test_get_sequencer_platform_from_instrument_name_and_settings(self):
        """
        get_sequencer_platform: use instrument name and settings
        """
        settings_ini = os.path.join(self.dirn,"settings.ini")
        with open(settings_ini,'w') as s:
            s.write("""[sequencers]
SN7001251 = hiseq
M00880 = miseq
FS10000171 = iseq
""")
        settings = Settings(settings_ini)
        print(settings.sequencers)
        self.assertEqual(get_sequencer_platform(
            "/mnt/data/120919_BLAH_0035_BC133VACXX",
            instrument="SN7001251",
            settings=settings),"hiseq")
        self.assertEqual(get_sequencer_platform(
            "/mnt/data/121210_BLAH_0001_000000000-A2Y1L",
            instrument="M00880",
            settings=settings),"miseq")
        self.assertEqual(get_sequencer_platform(
            "/mnt/data/20180829_BLAH_3_BNT40323-1530",
            instrument="FS10000171",
            settings=settings),"iseq")

    def test_get_sequencer_platform_unknown_instrument(self):
        """
        get_sequencer_platform: handle unknown instrument name
        """
        settings_ini = os.path.join(self.dirn,"settings.ini")
        with open(settings_ini,'w') as s:
            s.write("""[sequencers]
SN7001251 = hiseq
M00880 = miseq
FS10000171 = iseq
""")
        settings = Settings(settings_ini)
        print(settings.sequencers)
        self.assertEqual(get_sequencer_platform(
            "/mnt/data/120919_BLAH_0035_BC133VACXX"),None)
        self.assertEqual(get_sequencer_platform(
            "/mnt/data/120919_BLAH_0035_BC133VACXX",
            settings=settings),None)
        self.assertEqual(get_sequencer_platform(
            "/mnt/data/120919_BLAH_0035_BC133VACXX",
            instrument="BLEURGH",
            settings=settings),None)

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

    def test_require_bcl2fastq_version_eq_217_no_operator(self):
        """
        available_bcl2fastq_versions: require version '2.17.1.14' (no '==')
        """
        self.mockbcl2fastq.bcl2fastq_217()
        self.mockbcl2fastq.set_path()
        self.assertEqual(available_bcl2fastq_versions('2.17.1.14'),
                         self.mockbcl2fastq.exes)

    def test_require_bcl2fastq_version_with_versionless_package_on_path(self):
        """
        available_bcl2fastq_versions: require specific version when 'versionless' package is also present
        """
        self.mockbcl2fastq.casava_no_version()
        self.mockbcl2fastq.bcl2fastq_184()
        self.mockbcl2fastq.bcl2fastq_217()
        self.mockbcl2fastq.set_path()
        # Only last two executables will be returned, ordered
        # by version (2.17 then 1.8.4, i.e. reverse of order on
        # PATH)
        self.assertEqual(available_bcl2fastq_versions('>=1.8.4'),
                         self.mockbcl2fastq.exes[1:][::-1])

    def test_specify_path_for_bcl2fastq_184(self):
        """
        available_bcl2fastq_versions: specify path for bcl2fastq 1.8.4

        """
        self.mockbcl2fastq.bcl2fastq_217()
        self.mockbcl2fastq.bcl2fastq_184()
        self.mockbcl2fastq.set_path()
        # Get the full path for the bcl2fastq 1.8.4 exe
        bcl2fastq_184 = filter(lambda x: x.count('1.8.4') == 1,
                               self.mockbcl2fastq.exes)[0]
        self.assertEqual(available_bcl2fastq_versions(
            paths=(os.path.dirname(bcl2fastq_184),)),
                         [bcl2fastq_184,])

    def test_specify_path_for_bcl2fastq_184_exe(self):
        """
        available_bcl2fastq_versions: specify path for bcl2fastq 1.8.4 exe

        """
        self.mockbcl2fastq.bcl2fastq_217()
        self.mockbcl2fastq.bcl2fastq_184()
        self.mockbcl2fastq.set_path()
        # Get the full path for the bcl2fastq 1.8.4 exe
        bcl2fastq_184 = filter(lambda x: x.count('1.8.4') == 1,
                               self.mockbcl2fastq.exes)[0]
        self.assertEqual(available_bcl2fastq_versions(
            paths=(bcl2fastq_184,)),
                         [bcl2fastq_184,])

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

    def test_bcl2fastq_2_20_0_422(self):
        """
        Collect info for bcl2fastq 2.20.0.422

        """
        # bcl2fastq 2.20.0.422
        self.mockbcl2fastq.bcl2fastq_220()
        self.mockbcl2fastq.set_path()
        exe = self.mockbcl2fastq.exes[0]
        path,name,version = bcl_to_fastq_info()
        self.assertEqual(path,exe)
        self.assertEqual(name,'bcl2fastq')
        self.assertEqual(version,'2.20.0.422')

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
        self.assertEqual(path,'')
        self.assertEqual(name,'')
        self.assertEqual(version,'')

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
        self.assertEqual(name,'')
        self.assertEqual(version,'')

class TestMakeCustomSampleSheet(unittest.TestCase):
    """Tests for the make_custom_sample_sheet function
    """
    def setUp(self):
        # Create a temporary working dir
        self.wd = tempfile.mkdtemp()
        self.iem_content = """[Header]
IEMFileVersion,4

[Reads]
101
101

[Settings]
ReverseComplement,0
Adapter,CTGTCTCTTATACACATCT

[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,PJB1-1579,PJB1-1579,,,N701,CGATGTAT,N501,TCTTTCCC,PeterBriggs,
1,PJB2-1580,PJB2-1580,,,N702,TGACCAAT,N502,TCTTTCCC,PeterBriggs,
2,PJB1-1579,PJB1-1579,,,N701,CGATGTAT,N501,TCTTTCCC,PeterBriggs,
2,PJB2-1580,PJB2-1580,,,N702,TGACCAAT,N502,TCTTTCCC,PeterBriggs,
3,PJB1-1579,PJB1-1579,,,N701,CGATGTAT,N501,TCTTTCCC,PeterBriggs,
3,PJB2-1580,PJB2-1580,,,N702,TGACCAAT,N502,TCTTTCCC,PeterBriggs,
4,PJB1-1579,PJB1-1579,,,N701,CGATGTAT,N501,TCTTTCCC,PeterBriggs,
4,PJB2-1580,PJB2-1580,,,N702,TGACCAAT,N502,TCTTTCCC,PeterBriggs,
"""
    def tearDown(self):
        # Remove working dir
        if self.wd is not None:
            shutil.rmtree(self.wd)
    def test_make_custom_sample_sheet_lanes(self):
        """
        make_custom_sample_sheet: output subset of lanes
        """
        sample_sheet_in = os.path.join(self.wd,
                                       "SampleSheet.csv")
        with open(sample_sheet_in,'w') as fp:
            fp.write(self.iem_content)
        sample_sheet_out = os.path.join(self.wd,
                                       "custom_SampleSheet.csv")
        make_custom_sample_sheet(sample_sheet_in,
                                 output_sample_sheet=sample_sheet_out,
                                 lanes=[3,4])
        expected = """[Header]
IEMFileVersion,4

[Reads]
101
101

[Settings]
ReverseComplement,0
Adapter,CTGTCTCTTATACACATCT

[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
3,PJB1-1579,PJB1-1579,,,N701,CGATGTAT,N501,TCTTTCCC,PeterBriggs,
3,PJB2-1580,PJB2-1580,,,N702,TGACCAAT,N502,TCTTTCCC,PeterBriggs,
4,PJB1-1579,PJB1-1579,,,N701,CGATGTAT,N501,TCTTTCCC,PeterBriggs,
4,PJB2-1580,PJB2-1580,,,N702,TGACCAAT,N502,TCTTTCCC,PeterBriggs,
"""
        self.assertTrue(os.path.exists(sample_sheet_out))
        self.assertEqual(open(sample_sheet_out,'r').read(),
                         expected)

class TestBasesMaskIsValid(unittest.TestCase):
    """Tests for the bases_mask_is_valid function
    """
    def test_bases_mask_is_valid(self):
        self.assertTrue(bases_mask_is_valid('y50'))
        self.assertTrue(bases_mask_is_valid('y50,I4'))
        self.assertTrue(bases_mask_is_valid('y50,I6'))
        self.assertTrue(bases_mask_is_valid('y101,I6,y101'))
        self.assertTrue(bases_mask_is_valid('y250,I8,I8,y250'))
        self.assertTrue(bases_mask_is_valid('y250,I6nn,I6nn,y250'))
        self.assertTrue(bases_mask_is_valid('y250,I6n2,I6n2,y250'))
        self.assertTrue(bases_mask_is_valid('y25n76,I8,I8,y101'))
        self.assertTrue(bases_mask_is_valid('y250,I16'))
        self.assertTrue(bases_mask_is_valid('yyyyyyyyy,IIIIII'))
        self.assertTrue(bases_mask_is_valid('yyyyyyyyy,IIIInn'))
        self.assertFalse(bases_mask_is_valid(123))

class TestGetBasesMask(unittest.TestCase):
    """Tests for the get_bases_mask function
    """
    def setUp(self):
        # Create a temporary working dir
        self.wd = tempfile.mkdtemp()

    def tearDown(self):
        # Remove working dir
        if self.wd is not None:
            shutil.rmtree(self.wd)

    def test_get_bases_mask_single_index(self):
        """get_bases_mask: handle single index
        """
        # Make a single index RunInfo.xml file
        run_info_xml = os.path.join(self.wd,"RunInfo.xml")
        with open(run_info_xml,'w') as fp:
            fp.write(RunInfoXml.nextseq("171020_NB500968_00002_AHGXXXX"))
        # Make a matching sample sheet
        sample_sheet_content = """[Header]
IEMFileVersion,4
Date,11/23/2015
Workflow,GenerateFASTQ
Application,FASTQ Only
Assay,TruSeq HT
Description,
Chemistry,Amplicon

[Reads]
76
76

[Settings]
ReverseComplement,0
Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description
AB1,AB1,,,D701,CGTGTA,AB,
AB2,AB2,,,D702,ATTCAG,AB,
"""
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'w') as fp:
            fp.write(sample_sheet_content)
        # Check the bases mask
        self.assertEqual(get_bases_mask(run_info_xml,sample_sheet),
                         "y76,I6,y76")

    def test_get_bases_mask_dual_index(self):
        """get_bases_mask: handle dual index
        """
        # Make a RunInfo.xml file
        run_info_xml = os.path.join(self.wd,"RunInfo.xml")
        with open(run_info_xml,'w') as fp:
            fp.write(RunInfoXml.hiseq("171020_SN7001250_00002_AHGXXXX"))
        # Make a matching sample sheet
        sample_sheet_content = """[Header]
IEMFileVersion,4
Date,11/23/2015
Workflow,GenerateFASTQ
Application,FASTQ Only
Assay,TruSeq HT
Description,
Chemistry,Amplicon

[Reads]
76
76

[Settings]
ReverseComplement,0
Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index,index2,Sample_Project,Description
1,AB1,AB1,,,D701,CGTGTAGG,D501,GACCTGTA,AB,
2,AB2,AB2,,,D701,CGTGTAGG,D501,GACCTGTA,AB,
3,AB1,AB1,,,D701,CGTGTAGG,D501,GACCTGTA,AB,
4,AB2,AB2,,,D701,CGTGTAGG,D501,GACCTGTA,AB,
5,CD1,CD1,,,D701,CGTGTAGG,D501,GACCTGTA,CD,
6,CD2,CD2,,,D701,CGTGTAGG,D501,GACCTGTA,CD,
7,CD1,CD1,,,D701,CGTGTAGG,D501,GACCTGTA,CD,
8,CD2,CD2,,,D701,CGTGTAGG,D501,GACCTGTA,CD,
"""
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'w') as fp:
            fp.write(sample_sheet_content)
        # Check the bases mask
        self.assertEqual(get_bases_mask(run_info_xml,sample_sheet),
                         "y101,I8,I8,y101")

    def test_get_bases_mask_single_index_truncated(self):
        """get_bases_mask: handle truncated single index
        """
        # Make a single index RunInfo.xml file
        run_info_xml = os.path.join(self.wd,"RunInfo.xml")
        with open(run_info_xml,'w') as fp:
            fp.write(RunInfoXml.nextseq("171020_NB500968_00002_AHGXXXX"))
        # Make a matching sample sheet
        sample_sheet_content = """[Header]
IEMFileVersion,4
Date,11/23/2015
Workflow,GenerateFASTQ
Application,FASTQ Only
Assay,TruSeq HT
Description,
Chemistry,Amplicon

[Reads]
76
76

[Settings]
ReverseComplement,0
Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description
AB1,AB1,,,D701,CGTGT,AB,
AB2,AB2,,,D702,ATTCA,AB,
"""
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'w') as fp:
            fp.write(sample_sheet_content)
        # Check the bases mask
        self.assertEqual(get_bases_mask(run_info_xml,sample_sheet),
                         "y76,I5n,y76")

    def test_get_bases_mask_dual_index_truncated(self):
        """get_bases_mask: handle truncated dual index
        """
        # Make a RunInfo.xml file
        run_info_xml = os.path.join(self.wd,"RunInfo.xml")
        with open(run_info_xml,'w') as fp:
            fp.write(RunInfoXml.hiseq("171020_SN7001250_00002_AHGXXXX"))
        # Make a matching sample sheet
        sample_sheet_content = """[Header]
IEMFileVersion,4
Date,11/23/2015
Workflow,GenerateFASTQ
Application,FASTQ Only
Assay,TruSeq HT
Description,
Chemistry,Amplicon

[Reads]
76
76

[Settings]
ReverseComplement,0
Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index,index2,Sample_Project,Description
1,AB1,AB1,,,D701,CGTGTA,,,AB,
2,AB2,AB2,,,D701,CGTGTA,,,AB,
3,AB1,AB1,,,D701,CGTGTA,,,AB,
4,AB2,AB2,,,D701,CGTGTA,,,AB,
5,CD1,CD1,,,D701,CGTGTA,,,CD,
6,CD2,CD2,,,D701,CGTGTA,,,CD,
7,CD1,CD1,,,D701,CGTGTA,,,CD,
8,CD2,CD2,,,D701,CGTGTA,,,CD,
"""
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'w') as fp:
            fp.write(sample_sheet_content)
        # Check the bases mask
        self.assertEqual(get_bases_mask(run_info_xml,sample_sheet),
                         "y101,I6nn,nnnnnnnn,y101")

    def test_get_bases_mask_single_index_no_barcode(self):
        """get_bases_mask: handle single index with no barcode
        """
        # Make a single index RunInfo.xml file
        run_info_xml = os.path.join(self.wd,"RunInfo.xml")
        with open(run_info_xml,'w') as fp:
            fp.write(RunInfoXml.nextseq("171020_NB500968_00002_AHGXXXX"))
        # Make a matching sample sheet
        sample_sheet_content = """[Header]
IEMFileVersion,4
Date,11/23/2015
Workflow,GenerateFASTQ
Application,FASTQ Only
Assay,TruSeq HT
Description,
Chemistry,Amplicon

[Reads]
76
76

[Settings]
ReverseComplement,0
Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description
AB1,AB1,,,,,AB,
"""
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'w') as fp:
            fp.write(sample_sheet_content)
        # Check the bases mask
        self.assertEqual(get_bases_mask(run_info_xml,sample_sheet),
                         "y76,nnnnnn,y76")

    def test_get_bases_mask_dual_index_no_barcode(self):
        """get_bases_mask: handle dual index with no barcode
        """
        # Make a RunInfo.xml file
        run_info_xml = os.path.join(self.wd,"RunInfo.xml")
        with open(run_info_xml,'w') as fp:
            fp.write(RunInfoXml.hiseq("171020_SN7001250_00002_AHGXXXX"))
        # Make a matching sample sheet
        sample_sheet_content = """[Header]
IEMFileVersion,4
Date,11/23/2015
Workflow,GenerateFASTQ
Application,FASTQ Only
Assay,TruSeq HT
Description,
Chemistry,Amplicon

[Reads]
76
76

[Settings]
ReverseComplement,0
Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index,index2,Sample_Project,Description
1,AB1,AB1,,,,,,,AB,
2,AB2,AB2,,,,,,,AB,
3,AB1,AB1,,,,,,,AB,
4,AB2,AB2,,,,,,,AB,
5,CD1,CD1,,,,,,,CD,
6,CD2,CD2,,,,,,,CD,
7,CD1,CD1,,,,,,,CD,
8,CD2,CD2,,,,,,,CD,
"""
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'w') as fp:
            fp.write(sample_sheet_content)
        # Check the bases mask
        self.assertEqual(get_bases_mask(run_info_xml,sample_sheet),
                         "y101,nnnnnnnn,nnnnnnnn,y101")

    def test_get_bases_mask_10xgenomics_sample_set(self):
        """get_bases_mask: handle 10xGenomics sample set IDs
        """
        # Make a RunInfo.xml file
        run_info_xml = os.path.join(self.wd,"RunInfo.xml")
        with open(run_info_xml,'w') as fp:
            fp.write(RunInfoXml.nextseq("171020_NB500968_00002_AHGXXXX"))
        # Make a matching sample sheet
        sample_sheet_content = """[Header]
IEMFileVersion,4
Date,4/24/2018
Workflow,GenerateFASTQ
Application,NextSeq FASTQ Only
Assay,Nextera XT v2 Set A
Description,
Chemistry,Default

[Reads]
76
76

[Settings]
Adapter,CTGTCTCTTATACACATCT

[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description
SL1,SL1,,,N701,SI-GA-A2,SL,
SL2,SL2,,,N702,SI-GA-B2,SL,
SL3,SL3,,,N703,SI-GA-C2,SL,
SL4,SL4,,,N704,SI-GA-D2,SL,
"""
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'w') as fp:
            fp.write(sample_sheet_content)
        # Check the bases mask
        self.assertEqual(get_bases_mask(run_info_xml,sample_sheet),
                         "y76,I6,y76")

    def test_get_bases_mask_single_index_no_sample_sheet(self):
        """get_bases_mask: handle single index (no samplesheet)
        """
        # Make a single index RunInfo.xml file
        run_info_xml = os.path.join(self.wd,"RunInfo.xml")
        with open(run_info_xml,'w') as fp:
            fp.write(RunInfoXml.nextseq("171020_NB500968_00002_AHGXXXX"))
        # Check the bases mask
        self.assertEqual(get_bases_mask(run_info_xml),"y76,I6,y76")

    def test_get_bases_mask_dual_index_no_sample_sheet(self):
        """get_bases_mask: handle dual index (no samplesheet)
        """
        # Make a RunInfo.xml file
        run_info_xml = os.path.join(self.wd,"RunInfo.xml")
        with open(run_info_xml,'w') as fp:
            fp.write(RunInfoXml.hiseq("171020_SN7001250_00002_AHGXXXX"))
        # Check the bases mask
        self.assertEqual(get_bases_mask(run_info_xml),"y101,I8,I8,y101")

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
        self.assertEqual(get_nmismatches('yyyyyyyyy,IIIIII'),1)
        self.assertEqual(get_nmismatches('yyyyyyyyy,IIIINN'),0)
        self.assertEqual(get_nmismatches('yyyyyyyyy,IIII'),0)
        self.assertEqual(get_nmismatches('yyyyyyyyy,IIIIN4'),0)
        self.assertEqual(get_nmismatches('y250,I4,I4,y250'),1)

    def test_n_mismatches_invalid_input(self):
        self.assertRaises(Exception,get_nmismatches,'auto')
        self.assertRaises(Exception,get_nmismatches,123)
