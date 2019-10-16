#######################################################################
# Tests for archive_cmd.py module
#######################################################################

import unittest
import tempfile
import shutil
import os
import pwd
import grp
from auto_process_ngs.settings import Settings
from auto_process_ngs.auto_processor import AutoProcess
from auto_process_ngs.mock import MockAnalysisDirFactory
from auto_process_ngs.mock import UpdateAnalysisProject
from auto_process_ngs.commands.archive_cmd import archive

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

# Unit tests

class TestArchiveCommand(unittest.TestCase):
    """
    Tests for the 'archive' command
    """
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestArchiveCommand')
        # Create settings instance
        # This allows us to set the polling interval for the
        # unit tests
        settings_ini = os.path.join(self.dirn,"auto_process.ini")
        with open(settings_ini,'w') as s:
            s.write("""[general]
poll_interval = 0.5
""")
        self.settings = Settings(settings_ini)
        # Store original location so we can get back at the end
        self.pwd = os.getcwd()
        # Move to working dir
        os.chdir(self.dirn)

    def tearDown(self):
        # Return to original dir
        os.chdir(self.pwd)
        # Remove the temporary test directory
        def del_rw(action,name,excinfo):
            # Explicitly remove read only files/
            # dirs
            if os.path.isfile(name):
                os.chmod(os.path.dirname(name),0755)
                os.chmod(name,0655)
                os.remove(name)
            elif os.path.isdir(name):
                os.chmod(os.path.dirname(name),0755)
                os.chmod(name,0755)
                os.rmdir(name)
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.dirn,onerror=del_rw)

    def test_archive_to_staging(self):
        """archive: test copying to staging archive dir
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            top_dir=self.dirn)
        mockdir.create()
        # Make a mock archive directory
        archive_dir = os.path.join(self.dirn,"archive")
        final_dir = os.path.join(archive_dir,
                                 "2017",
                                 "miseq")
        os.makedirs(final_dir)
        self.assertTrue(os.path.isdir(final_dir))
        self.assertEqual(len(os.listdir(final_dir)),0)
        # Make autoprocess instance and set required metadata
        ap = AutoProcess(analysis_dir=mockdir.dirn,
                         settings=self.settings)
        ap.set_metadata("source","testing")
        ap.set_metadata("run_number","87")
        # Do archiving op
        status = archive(ap,
                         archive_dir=archive_dir,
                         year='2017',platform='miseq',
                         read_only_fastqs=False,
                         final=False)
        self.assertEqual(status,0)
        # Check that staging dir exists
        staging_dir = os.path.join(
            final_dir,
            "__170901_M00879_0087_000000000-AGEW9_analysis.pending")
        self.assertTrue(os.path.exists(staging_dir))
        self.assertEqual(len(os.listdir(final_dir)),1)
        # Check contents
        dirs = ("AB","CDE","logs","undetermined")
        for d in dirs:
            d = os.path.join(staging_dir,d)
            self.assertTrue(os.path.exists(d))
        files = ("auto_process.info",
                 "custom_SampleSheet.csv",
                 "metadata.info",
                 "projects.info",
                 "SampleSheet.orig.csv")
        for f in files:
            f = os.path.join(staging_dir,f)
            self.assertTrue(os.path.exists(f))

    def test_archive_to_staging_set_group(self):
        """archive: test copying to staging archive dir and set group
        """
        # Find groups for current user
        current_user = pwd.getpwuid(os.getuid()).pw_name
        groups = [g.gr_gid
                  for g in grp.getgrall()
                  if current_user in g.gr_mem]
        if len(groups) < 2:
            raise unittest.SkipTest("user '%s' must be in at least "
                                    "two groups for this test" %
                                    current_user)
        # Find a group to set archived files to
        current_gid = os.stat(self.dirn).st_gid
        new_group = None
        for gid in groups:
            if gid != current_gid:
                new_group = gid
                break
        self.assertTrue(new_group is not None)
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            top_dir=self.dirn)
        mockdir.create()
        # Make a mock archive directory
        archive_dir = os.path.join(self.dirn,"archive")
        final_dir = os.path.join(archive_dir,
                                 "2017",
                                 "miseq")
        os.makedirs(final_dir)
        self.assertTrue(os.path.isdir(final_dir))
        self.assertEqual(len(os.listdir(final_dir)),0)
        # Make autoprocess instance and set required metadata
        ap = AutoProcess(analysis_dir=mockdir.dirn,
                         settings=self.settings)
        ap.set_metadata("source","testing")
        ap.set_metadata("run_number","87")
        # Do archiving op
        status = archive(ap,
                         archive_dir=archive_dir,
                         year='2017',platform='miseq',
                         read_only_fastqs=False,
                         group=new_group,
                         final=False)
        self.assertEqual(status,0)
        # Check that staging dir exists
        staging_dir = os.path.join(
            final_dir,
            "__170901_M00879_0087_000000000-AGEW9_analysis.pending")
        self.assertTrue(os.path.exists(staging_dir))
        self.assertEqual(len(os.listdir(final_dir)),1)
        # Check group of staging dir
        self.assertEqual(os.stat(staging_dir).st_gid,new_group)
        # Check contents
        dirs = ("AB","CDE","logs","undetermined")
        for d in dirs:
            d = os.path.join(staging_dir,d)
            self.assertTrue(os.path.exists(d))
            self.assertEqual(os.stat(d).st_gid,new_group)
        files = ("auto_process.info",
                 "custom_SampleSheet.csv",
                 "metadata.info",
                 "projects.info",
                 "SampleSheet.orig.csv")
        for f in files:
            f = os.path.join(staging_dir,f)
            self.assertTrue(os.path.exists(f))
            self.assertEqual(os.stat(f).st_gid,new_group)

    def test_archive_to_final(self):
        """archive: test copying to final archive dir
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            top_dir=self.dirn)
        mockdir.create()
        # Make a mock archive directory
        archive_dir = os.path.join(self.dirn,"archive")
        final_dir = os.path.join(archive_dir,
                                 "2017",
                                 "miseq")
        os.makedirs(final_dir)
        self.assertTrue(os.path.isdir(final_dir))
        self.assertEqual(len(os.listdir(final_dir)),0)
        # Make autoprocess instance and set required metadata
        ap = AutoProcess(analysis_dir=mockdir.dirn,
                         settings=self.settings)
        ap.set_metadata("source","testing")
        ap.set_metadata("run_number","87")
        # Do archiving op
        status = archive(ap,
                         archive_dir=archive_dir,
                         year='2017',platform='miseq',
                         read_only_fastqs=False,
                         final=True)
        self.assertEqual(status,0)
        # Check that final dir exists
        final_archive_dir = os.path.join(
            final_dir,
            "170901_M00879_0087_000000000-AGEW9_analysis")
        self.assertTrue(os.path.exists(final_archive_dir))
        self.assertEqual(len(os.listdir(final_dir)),1)
        # Check contents
        dirs = ("AB","CDE","logs","undetermined")
        for d in dirs:
            d = os.path.join(final_archive_dir,d)
            self.assertTrue(os.path.exists(d))
        files = ("auto_process.info",
                 "custom_SampleSheet.csv",
                 "metadata.info",
                 "projects.info",
                 "SampleSheet.orig.csv")
        for f in files:
            f = os.path.join(final_archive_dir,f)
            self.assertTrue(os.path.exists(f))
        # Check that Fastqs are not writable
        for project in ("AB","CDE","undetermined"):
            fq_dir = os.path.join(final_archive_dir,
                                  project,
                                  "fastqs")
            self.assertTrue(os.path.exists(fq_dir))
            fqs = os.listdir(fq_dir)
            self.assertTrue(len(fqs) > 0)
            for fq in fqs:
                fq = os.path.join(fq_dir,fq)
                self.assertTrue(os.access(fq,os.R_OK))
                self.assertTrue(os.access(fq,os.W_OK))

    def test_archive_to_final_read_only_fastqs(self):
        """archive: test copying to final archive dir (read-only Fastqs)
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            top_dir=self.dirn)
        mockdir.create()
        # Make a mock archive directory
        archive_dir = os.path.join(self.dirn,"archive")
        final_dir = os.path.join(archive_dir,
                                 "2017",
                                 "miseq")
        os.makedirs(final_dir)
        self.assertTrue(os.path.isdir(final_dir))
        self.assertEqual(len(os.listdir(final_dir)),0)
        # Make autoprocess instance and set required metadata
        ap = AutoProcess(analysis_dir=mockdir.dirn,
                         settings=self.settings)
        ap.set_metadata("source","testing")
        ap.set_metadata("run_number","87")
        # Do archiving op
        status = archive(ap,
                         archive_dir=archive_dir,
                         year='2017',platform='miseq',
                         read_only_fastqs=True,
                         final=True)
        self.assertEqual(status,0)
        # Check that final dir exists
        final_archive_dir = os.path.join(
            final_dir,
            "170901_M00879_0087_000000000-AGEW9_analysis")
        self.assertTrue(os.path.exists(final_archive_dir))
        self.assertEqual(len(os.listdir(final_dir)),1)
        # Check contents
        dirs = ("AB","CDE","logs","undetermined")
        for d in dirs:
            d = os.path.join(final_archive_dir,d)
            self.assertTrue(os.path.exists(d))
        files = ("auto_process.info",
                 "custom_SampleSheet.csv",
                 "metadata.info",
                 "projects.info",
                 "SampleSheet.orig.csv")
        for f in files:
            f = os.path.join(final_archive_dir,f)
            self.assertTrue(os.path.exists(f))
        # Check that Fastqs are not writable
        for project in ("AB","CDE","undetermined"):
            fq_dir = os.path.join(final_archive_dir,
                                  project,
                                  "fastqs")
            self.assertTrue(os.path.exists(fq_dir))
            fqs = os.listdir(fq_dir)
            self.assertTrue(len(fqs) > 0)
            for fq in fqs:
                fq = os.path.join(fq_dir,fq)
                self.assertTrue(os.access(fq,os.R_OK))
                self.assertFalse(os.access(fq,os.W_OK))

    def test_archive_to_final_multiple_fastq_sets_read_only_fastqs(self):
        """archive: test copying multiple fastq sets to final archive dir (read-only Fastqs)
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            top_dir=self.dirn)
        mockdir.create()
        # Make a mock archive directory
        archive_dir = os.path.join(self.dirn,"archive")
        final_dir = os.path.join(archive_dir,
                                 "2017",
                                 "miseq")
        os.makedirs(final_dir)
        self.assertTrue(os.path.isdir(final_dir))
        self.assertEqual(len(os.listdir(final_dir)),0)
        # Make autoprocess instance and set required metadata
        ap = AutoProcess(analysis_dir=mockdir.dirn,
                         settings=self.settings)
        ap.set_metadata("source","testing")
        ap.set_metadata("run_number","87")
        # Add additional fastq set for first project
        multi_fastqs_project = ap.get_analysis_projects()[0]
        UpdateAnalysisProject(multi_fastqs_project).add_fastq_set(
            "fastqs.extra",
            ("Alt1.r1.fastq.gz","Alt2.r1.fastq.gz"))
        # Do archiving op
        status = archive(ap,
                         archive_dir=archive_dir,
                         year='2017',platform='miseq',
                         read_only_fastqs=True,
                         final=True)
        self.assertEqual(status,0)
        # Check that final dir exists
        final_archive_dir = os.path.join(
            final_dir,
            "170901_M00879_0087_000000000-AGEW9_analysis")
        self.assertTrue(os.path.exists(final_archive_dir))
        self.assertEqual(len(os.listdir(final_dir)),1)
        # Check contents
        dirs = ("AB","CDE","logs","undetermined")
        for d in dirs:
            d = os.path.join(final_archive_dir,d)
            self.assertTrue(os.path.exists(d))
        files = ("auto_process.info",
                 "custom_SampleSheet.csv",
                 "metadata.info",
                 "projects.info",
                 "SampleSheet.orig.csv")
        for f in files:
            f = os.path.join(final_archive_dir,f)
            self.assertTrue(os.path.exists(f))
        # Check that Fastqs are not writable
        for project in ("AB","CDE","undetermined"):
            fq_dir = os.path.join(final_archive_dir,
                                  project,
                                  "fastqs")
            self.assertTrue(os.path.exists(fq_dir))
            fqs = os.listdir(fq_dir)
            self.assertTrue(len(fqs) > 0)
            for fq in fqs:
                fq = os.path.join(fq_dir,fq)
                self.assertTrue(os.access(fq,os.R_OK))
                self.assertFalse(os.access(fq,os.W_OK))
        # Check additional Fastqs are not writable
        fq_dir = os.path.join(final_archive_dir,
                              multi_fastqs_project.name,
                              "fastqs.extra")
        self.assertTrue(os.path.exists(fq_dir))
        fqs = os.listdir(fq_dir)
        self.assertTrue(len(fqs) > 0)
        for fq in fqs:
            fq = os.path.join(fq_dir,fq)
            self.assertTrue(os.access(fq,os.R_OK))
            self.assertFalse(os.access(fq,os.W_OK))

    def test_archive_to_final_via_staging(self):
        """archive: test copying to staging then final archive dir
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            top_dir=self.dirn)
        mockdir.create()
        # Make a mock archive directory
        archive_dir = os.path.join(self.dirn,"archive")
        final_dir = os.path.join(archive_dir,
                                 "2017",
                                 "miseq")
        os.makedirs(final_dir)
        self.assertTrue(os.path.isdir(final_dir))
        self.assertEqual(len(os.listdir(final_dir)),0)
        # Make autoprocess instance and set required metadata
        ap = AutoProcess(analysis_dir=mockdir.dirn,
                         settings=self.settings)
        ap.set_metadata("source","testing")
        ap.set_metadata("run_number","87")
        # Do staging archiving op
        status = archive(ap,
                         archive_dir=archive_dir,
                         year='2017',platform='miseq',
                         read_only_fastqs=False,
                         final=False)
        self.assertEqual(status,0)
        # Check that staging dir exists
        staging_dir = os.path.join(
            final_dir,
            "__170901_M00879_0087_000000000-AGEW9_analysis.pending")
        final_archive_dir = os.path.join(
            final_dir,
            "170901_M00879_0087_000000000-AGEW9_analysis")
        self.assertTrue(os.path.exists(staging_dir))
        self.assertFalse(os.path.exists(final_archive_dir))
        self.assertEqual(len(os.listdir(final_dir)),1)
        # Do final archiving op
        status = archive(ap,
                         archive_dir=archive_dir,
                         year='2017',platform='miseq',
                         read_only_fastqs=False,
                         final=True)
        self.assertEqual(status,0)
        self.assertFalse(os.path.exists(staging_dir))
        self.assertTrue(os.path.exists(final_archive_dir))
        self.assertEqual(len(os.listdir(final_dir)),1)
        # Check contents
        dirs = ("AB","CDE","logs","undetermined")
        for d in dirs:
            d = os.path.join(final_archive_dir,d)
            self.assertTrue(os.path.exists(d))
        files = ("auto_process.info",
                 "custom_SampleSheet.csv",
                 "metadata.info",
                 "projects.info",
                 "SampleSheet.orig.csv")
        for f in files:
            f = os.path.join(final_archive_dir,f)
            self.assertTrue(os.path.exists(f))

    def test_archive_staging_to_final(self):
        """archive: test archiving directly from staging dir
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            top_dir=self.dirn)
        mockdir.create()
        # Make a mock archive directory
        archive_dir = os.path.join(self.dirn,"archive")
        final_dir = os.path.join(archive_dir,
                                 "2017",
                                 "miseq")
        os.makedirs(final_dir)
        self.assertTrue(os.path.isdir(final_dir))
        self.assertEqual(len(os.listdir(final_dir)),0)
        # Make autoprocess instance and set required metadata
        ap = AutoProcess(analysis_dir=mockdir.dirn,
                         settings=self.settings)
        ap.set_metadata("source","testing")
        ap.set_metadata("run_number","87")
        ap.save_metadata()
        # Move to the archive area as a "pending" directory
        os.rename(mockdir.dirn,
                  os.path.join(
                    archive_dir,
                    "2017",
                    "miseq",
                    "__170901_M00879_0087_000000000-AGEW9_analysis.pending"))
        # Load pending dir into a new autoprocess instance
        ap = AutoProcess(
            analysis_dir=os.path.join(
                archive_dir,
                "2017",
                "miseq",
                "__170901_M00879_0087_000000000-AGEW9_analysis.pending"))
        # Staging archiving attempt should fail
        self.assertRaises(Exception,
                          archive,
                          ap,
                          archive_dir=archive_dir,
                          year='2017',platform='miseq',
                          read_only_fastqs=False,
                          final=False)
        staging_dir = os.path.join(
            final_dir,
            "__170901_M00879_0087_000000000-AGEW9_analysis.pending")
        final_archive_dir = os.path.join(
            final_dir,
            "170901_M00879_0087_000000000-AGEW9_analysis")
        self.assertFalse(os.path.exists(final_archive_dir))
        self.assertEqual(len(os.listdir(final_dir)),1)
        # Copy to final should work
        status = archive(ap,
                         archive_dir=archive_dir,
                         year='2017',platform='miseq',
                         read_only_fastqs=False,
                         final=True)
        self.assertEqual(status,0)
        self.assertTrue(os.path.exists(final_archive_dir))
        self.assertFalse(os.path.exists(staging_dir))
        self.assertEqual(len(os.listdir(final_dir)),1)
        # Check contents
        dirs = ("AB","CDE","logs","undetermined")
        for d in dirs:
            d = os.path.join(final_archive_dir,d)
            self.assertTrue(os.path.exists(d))
        files = ("auto_process.info",
                 "custom_SampleSheet.csv",
                 "metadata.info",
                 "projects.info",
                 "SampleSheet.orig.csv")
        for f in files:
            f = os.path.join(final_archive_dir,f)
            self.assertTrue(os.path.exists(f))

    def test_archive_automatically_sets_correct_year(self):
        """archive: test archiving sets the year correctly if not specified
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            top_dir=self.dirn)
        mockdir.create()
        # Make a mock archive directory
        archive_dir = os.path.join(self.dirn,"archive")
        final_dir = os.path.join(archive_dir,
                                 "2017",
                                 "miseq")
        os.makedirs(final_dir)
        self.assertTrue(os.path.isdir(final_dir))
        self.assertEqual(len(os.listdir(final_dir)),0)
        # Make autoprocess instance and set required metadata
        ap = AutoProcess(analysis_dir=mockdir.dirn,
                         settings=self.settings)
        ap.set_metadata("source","testing")
        ap.set_metadata("run_number","87")
        # Do staging archiving op with no year
        status = archive(ap,
                         archive_dir=archive_dir,
                         platform='miseq',
                         read_only_fastqs=False,
                         final=False)
        self.assertEqual(status,0)
        # Check that staging dir exists
        staging_dir = os.path.join(
            final_dir,
            "__170901_M00879_0087_000000000-AGEW9_analysis.pending")
        final_archive_dir = os.path.join(
            final_dir,
            "170901_M00879_0087_000000000-AGEW9_analysis")
        self.assertTrue(os.path.exists(staging_dir))
        self.assertFalse(os.path.exists(final_archive_dir))
        self.assertEqual(len(os.listdir(final_dir)),1)
        # Do final archiving op with no year
        status = archive(ap,
                         archive_dir=archive_dir,
                         platform='miseq',
                         read_only_fastqs=False,
                         final=True)
        self.assertEqual(status,0)
        self.assertFalse(os.path.exists(staging_dir))
        self.assertTrue(os.path.exists(final_archive_dir))
        self.assertEqual(len(os.listdir(final_dir)),1)
        # Check contents
        dirs = ("AB","CDE","logs","undetermined")
        for d in dirs:
            d = os.path.join(final_archive_dir,d)
            self.assertTrue(os.path.exists(d))
        files = ("auto_process.info",
                 "custom_SampleSheet.csv",
                 "metadata.info",
                 "projects.info",
                 "SampleSheet.orig.csv")
        for f in files:
            f = os.path.join(final_archive_dir,f)
            self.assertTrue(os.path.exists(f))

    def test_archive_handles_four_digit_year_in_datestamp(self):
        """archive: test archiving handles 4-digit year in datestamp
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '20170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "20170901" },
            top_dir=self.dirn)
        mockdir.create()
        # Make a mock archive directory
        archive_dir = os.path.join(self.dirn,"archive")
        final_dir = os.path.join(archive_dir,
                                 "2017",
                                 "miseq")
        os.makedirs(final_dir)
        self.assertTrue(os.path.isdir(final_dir))
        self.assertEqual(len(os.listdir(final_dir)),0)
        # Make autoprocess instance and set required metadata
        ap = AutoProcess(analysis_dir=mockdir.dirn,
                         settings=self.settings)
        ap.set_metadata("source","testing")
        ap.set_metadata("run_number","87")
        # Do staging archiving op with no year
        status = archive(ap,
                         archive_dir=archive_dir,
                         platform='miseq',
                         read_only_fastqs=False,
                         final=False)
        self.assertEqual(status,0)
        # Check that staging dir exists
        staging_dir = os.path.join(
            final_dir,
            "__20170901_M00879_0087_000000000-AGEW9_analysis.pending")
        final_archive_dir = os.path.join(
            final_dir,
            "20170901_M00879_0087_000000000-AGEW9_analysis")
        self.assertTrue(os.path.exists(staging_dir))
        self.assertFalse(os.path.exists(final_archive_dir))
        self.assertEqual(len(os.listdir(final_dir)),1)
        # Do final archiving op with no year
        status = archive(ap,
                         archive_dir=archive_dir,
                         platform='miseq',
                         read_only_fastqs=False,
                         final=True)
        self.assertEqual(status,0)
        self.assertFalse(os.path.exists(staging_dir))
        self.assertTrue(os.path.exists(final_archive_dir))
        self.assertEqual(len(os.listdir(final_dir)),1)
        # Check contents
        dirs = ("AB","CDE","logs","undetermined")
        for d in dirs:
            d = os.path.join(final_archive_dir,d)
            self.assertTrue(os.path.exists(d))
        files = ("auto_process.info",
                 "custom_SampleSheet.csv",
                 "metadata.info",
                 "projects.info",
                 "SampleSheet.orig.csv")
        for f in files:
            f = os.path.join(final_archive_dir,f)
            self.assertTrue(os.path.exists(f))

    def test_archive_oserror_if_destination_doesnt_exist(self):
        """archive: test archiving raises OSError if destination doesn't exist
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            top_dir=self.dirn)
        mockdir.create()
        # Make a mock archive directory
        archive_dir = os.path.join(self.dirn,"archive")
        os.makedirs(archive_dir)
        final_dir = os.path.join(archive_dir,
                                 "2017",
                                 "miseq")
        self.assertFalse(os.path.isdir(final_dir))
        staging_dir = os.path.join(
            final_dir,
            "__170901_M00879_0087_000000000-AGEW9_analysis.pending")
        final_archive_dir = os.path.join(
            final_dir,
            "170901_M00879_0087_000000000-AGEW9_analysis")
        # Make autoprocess instance and set required metadata
        ap = AutoProcess(analysis_dir=mockdir.dirn,
                         settings=self.settings)
        ap.set_metadata("source","testing")
        ap.set_metadata("run_number","87")
        # Staging attempt should fail
        self.assertRaises(OSError,
                          archive,
                          ap,
                          archive_dir=archive_dir,
                          year='2017',platform='miseq',
                          read_only_fastqs=False,
                          final=False)
        self.assertFalse(os.path.exists(staging_dir))
        # Final archiving attempt should fail
        self.assertRaises(OSError,
                          archive,
                          ap,
                          archive_dir=archive_dir,
                          year='2017',platform='miseq',
                          read_only_fastqs=False,
                          final=True)
        self.assertFalse(os.path.exists(final_archive_dir))
        # Make "YEAR" level in archive dir
        os.makedirs(os.path.join(archive_dir,"2017"))
        # Staging attempt should fail
        self.assertRaises(OSError,
                          archive,
                          ap,
                          archive_dir=archive_dir,
                          year='2017',platform='miseq',
                          read_only_fastqs=False,
                          final=False)
        self.assertFalse(os.path.exists(staging_dir))
        # Final archiving attempt should fail
        self.assertRaises(OSError,
                          archive,
                          ap,
                          archive_dir=archive_dir,
                          year='2017',platform='miseq',
                          read_only_fastqs=False,
                          final=True)
        self.assertFalse(os.path.exists(final_archive_dir))

    def test_archive_to_staging_ignores_bak_projects(self):
        """archive: check staging ignores .bak directories
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            top_dir=self.dirn)
        mockdir.create()
        # Make a mock archive directory
        archive_dir = os.path.join(self.dirn,"archive")
        final_dir = os.path.join(archive_dir,
                                 "2017",
                                 "miseq")
        os.makedirs(final_dir)
        self.assertTrue(os.path.isdir(final_dir))
        self.assertEqual(len(os.listdir(final_dir)),0)
        # Make autoprocess instance and set required metadata
        ap = AutoProcess(analysis_dir=mockdir.dirn,
                         settings=self.settings)
        ap.set_metadata("source","testing")
        ap.set_metadata("run_number","87")
        # Add a .bak project directory
        shutil.copytree(os.path.join(mockdir.dirn,"AB"),
                        os.path.join(mockdir.dirn,"AB.bak"))
        # Do archiving op
        status = archive(ap,
                         archive_dir=archive_dir,
                         year='2017',platform='miseq',
                         read_only_fastqs=False,
                         final=False)
        self.assertEqual(status,0)
        # Check that staging dir exists
        staging_dir = os.path.join(
            final_dir,
            "__170901_M00879_0087_000000000-AGEW9_analysis.pending")
        self.assertTrue(os.path.exists(staging_dir))
        self.assertEqual(len(os.listdir(final_dir)),1)
        # Check contents
        dirs = ("AB","CDE","logs","undetermined")
        for d in dirs:
            d = os.path.join(staging_dir,d)
            self.assertTrue(os.path.exists(d))
        files = ("auto_process.info",
                 "custom_SampleSheet.csv",
                 "metadata.info",
                 "projects.info",
                 "SampleSheet.orig.csv")
        for f in files:
            f = os.path.join(staging_dir,f)
            self.assertTrue(os.path.exists(f))
        # Check .bak directory wasn't copied
        dirs = ("AB.bak",)
        for d in dirs:
            d = os.path.join(staging_dir,d)
            self.assertFalse(os.path.exists(d))

    def test_archive_to_final_ignores_bak_projects(self):
        """archive: check final archiving ignores .bak directories
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            top_dir=self.dirn)
        mockdir.create()
        # Make a mock archive directory
        archive_dir = os.path.join(self.dirn,"archive")
        final_dir = os.path.join(archive_dir,
                                 "2017",
                                 "miseq")
        os.makedirs(final_dir)
        self.assertTrue(os.path.isdir(final_dir))
        self.assertEqual(len(os.listdir(final_dir)),0)
        # Make autoprocess instance and set required metadata
        ap = AutoProcess(analysis_dir=mockdir.dirn,
                         settings=self.settings)
        ap.set_metadata("source","testing")
        ap.set_metadata("run_number","87")
        # Add a .bak project directory
        shutil.copytree(os.path.join(mockdir.dirn,"AB"),
                        os.path.join(mockdir.dirn,"AB.bak"))
        # Do archiving op
        status = archive(ap,
                         archive_dir=archive_dir,
                         year='2017',platform='miseq',
                         read_only_fastqs=True,
                         final=True)
        self.assertEqual(status,0)
        # Check that final dir exists
        final_archive_dir = os.path.join(
            final_dir,
            "170901_M00879_0087_000000000-AGEW9_analysis")
        self.assertTrue(os.path.exists(final_archive_dir))
        self.assertEqual(len(os.listdir(final_dir)),1)
        # Check contents
        dirs = ("AB","CDE","logs","undetermined")
        for d in dirs:
            d = os.path.join(final_archive_dir,d)
            self.assertTrue(os.path.exists(d))
        files = ("auto_process.info",
                 "custom_SampleSheet.csv",
                 "metadata.info",
                 "projects.info",
                 "SampleSheet.orig.csv")
        for f in files:
            f = os.path.join(final_archive_dir,f)
            self.assertTrue(os.path.exists(f))
        # Check that Fastqs are not writable
        for project in ("AB","CDE","undetermined"):
            fq_dir = os.path.join(final_archive_dir,
                                  project,
                                  "fastqs")
            self.assertTrue(os.path.exists(fq_dir))
            fqs = os.listdir(fq_dir)
            self.assertTrue(len(fqs) > 0)
            for fq in fqs:
                fq = os.path.join(fq_dir,fq)
                self.assertTrue(os.access(fq,os.R_OK))
                self.assertFalse(os.access(fq,os.W_OK))
        # Check .bak directory wasn't copied
        dirs = ("AB.bak",)
        for d in dirs:
            d = os.path.join(final_archive_dir,d)
            self.assertFalse(os.path.exists(d))

    def test_archive_force_stage_no_projects(self):
        """archive: force staging of run with no projects
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Make a mock archive directory
        archive_dir = os.path.join(self.dirn,"archive")
        final_dir = os.path.join(archive_dir,
                                 "2017",
                                 "miseq")
        os.makedirs(final_dir)
        self.assertTrue(os.path.isdir(final_dir))
        self.assertEqual(len(os.listdir(final_dir)),0)
        # Make autoprocess instance and set required metadata
        ap = AutoProcess(analysis_dir=mockdir.dirn,
                         settings=self.settings)
        ap.set_metadata("source","testing")
        ap.set_metadata("run_number","87")
        # Do archiving op to staging
        status = archive(ap,
                         archive_dir=archive_dir,
                         year='2017',platform='miseq',
                         read_only_fastqs=True,
                         final=False,
                         force=True)
        self.assertEqual(status,0)
        # Check that staging dir exists
        staging_dir = os.path.join(
            final_dir,
            "__170901_M00879_0087_000000000-AGEW9_analysis.pending")
        self.assertTrue(os.path.exists(staging_dir))
        self.assertEqual(len(os.listdir(final_dir)),1)
        # Check contents
        dirs = ("bcl2fastq","logs",)
        for d in dirs:
            d = os.path.join(staging_dir,d)
            self.assertTrue(os.path.exists(d))
        files = ("auto_process.info",
                 "custom_SampleSheet.csv",
                 "metadata.info",
                 "projects.info",
                 "SampleSheet.orig.csv")
        for f in files:
            f = os.path.join(staging_dir,f)
            self.assertTrue(os.path.exists(f))

    def test_archive_force_no_projects(self):
        """archive: force final archiving of run with no projects
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Make a mock archive directory
        archive_dir = os.path.join(self.dirn,"archive")
        final_dir = os.path.join(archive_dir,
                                 "2017",
                                 "miseq")
        os.makedirs(final_dir)
        self.assertTrue(os.path.isdir(final_dir))
        self.assertEqual(len(os.listdir(final_dir)),0)
        # Make autoprocess instance and set required metadata
        ap = AutoProcess(analysis_dir=mockdir.dirn,
                         settings=self.settings)
        ap.set_metadata("source","testing")
        ap.set_metadata("run_number","87")
        # Do archiving op to staging
        status = archive(ap,
                         archive_dir=archive_dir,
                         year='2017',platform='miseq',
                         read_only_fastqs=True,
                         final=True,
                         force=True)
        self.assertEqual(status,0)
        # Check that final dir exists
        final_archive_dir = os.path.join(
            final_dir,
            "170901_M00879_0087_000000000-AGEW9_analysis")
        self.assertTrue(os.path.exists(final_archive_dir))
        self.assertEqual(len(os.listdir(final_dir)),1)
        # Check contents
        dirs = ("bcl2fastq","logs",)
        for d in dirs:
            d = os.path.join(final_archive_dir,d)
            self.assertTrue(os.path.exists(d))
        files = ("auto_process.info",
                 "custom_SampleSheet.csv",
                 "metadata.info",
                 "projects.info",
                 "SampleSheet.orig.csv")
        for f in files:
            f = os.path.join(final_archive_dir,f)
            self.assertTrue(os.path.exists(f))

    def test_archive_force_fails_for_no_data(self):
        """archive: force archiving fails if run has no data
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Remove the 'bcl2fastq' subdirectory
        shutil.rmtree(os.path.join(
            self.dirn,
            '170901_M00879_0087_000000000-AGEW9_analysis',
            'bcl2fastq'))
        # Make a mock archive directory
        archive_dir = os.path.join(self.dirn,"archive")
        final_dir = os.path.join(archive_dir,
                                 "2017",
                                 "miseq")
        os.makedirs(final_dir)
        self.assertTrue(os.path.isdir(final_dir))
        self.assertEqual(len(os.listdir(final_dir)),0)
        # Make autoprocess instance and set required metadata
        ap = AutoProcess(analysis_dir=mockdir.dirn,
                         settings=self.settings)
        ap.set_metadata("source","testing")
        ap.set_metadata("run_number","87")
        # Do archiving op to staging
        self.assertRaises(Exception,
                          archive,
                          ap,
                          archive_dir=archive_dir,
                          year='2017',platform='miseq',
                          read_only_fastqs=True,
                          final=False,
                          force=True)
        # Check that staging dir exists
        staging_dir = os.path.join(
            final_dir,
            "__170901_M00879_0087_000000000-AGEW9_analysis.pending")
        self.assertFalse(os.path.exists(staging_dir))

    def test_archive_removes_redundant_fastqs(self):
        """archive: remove redundant Fastqs from final archive (read-only Fastqs)
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            top_dir=self.dirn)
        mockdir.create()
        # Make a mock archive directory
        archive_dir = os.path.join(self.dirn,"archive")
        final_dir = os.path.join(archive_dir,
                                 "2017",
                                 "miseq")
        os.makedirs(final_dir)
        self.assertTrue(os.path.isdir(final_dir))
        self.assertEqual(len(os.listdir(final_dir)),0)
        # Make autoprocess instance and set required metadata
        ap = AutoProcess(analysis_dir=mockdir.dirn,
                         settings=self.settings)
        ap.set_metadata("source","testing")
        ap.set_metadata("run_number","87")
        # Do archiving to staging
        status = archive(ap,
                         archive_dir=archive_dir,
                         year='2017',platform='miseq',
                         read_only_fastqs=True,
                         final=False)
        self.assertEqual(status,0)
        # Check staging dir
        staging_dir = os.path.join(
            final_dir,
            "__170901_M00879_0087_000000000-AGEW9_analysis.pending")
        self.assertTrue(os.path.exists(staging_dir))
        dirs = ("AB","CDE","logs","undetermined")
        for d in dirs:
            d = os.path.join(staging_dir,d)
            self.assertTrue(os.path.exists(d))
        # Check Fastqs for project 'AB'
        initial_fastqs = ("AB1_S1_R1_001.fastq.gz",
                          "AB1_S1_R2_001.fastq.gz",
                          "AB2_S2_R1_001.fastq.gz",
                          "AB2_S2_R2_001.fastq.gz")
        for fq in initial_fastqs:
            fastq = os.path.join(staging_dir,
                                 "AB",
                                 "fastqs",
                                 fq)
            self.assertTrue(os.path.exists(fastq))
        # Update Fastqs in project 'AB': remove existing ones
        # and make new ones with different names
        new_fastqs = ("AB3_S3_R1_001.fastq.gz",
                      "AB3_S3_R2_001.fastq.gz",
                      "AB4_S4_R1_001.fastq.gz",
                      "AB4_S4_R2_001.fastq.gz")
        for fq in initial_fastqs:
            fastq = os.path.join(mockdir.dirn,
                                 "AB",
                                 "fastqs",
                                 fq)
            os.remove(fastq)
        for fq in new_fastqs:
            fastq = os.path.join(mockdir.dirn,
                                 "AB",
                                 "fastqs",
                                 fq)
            with open(fastq,'w') as fp:
                fp.write("")
        # Do second archive operation to final
        status = archive(ap,
                         archive_dir=archive_dir,
                         year='2017',platform='miseq',
                         read_only_fastqs=True,
                         final=True)
        self.assertEqual(status,0)
        # Check final dir
        final_archive_dir = os.path.join(
            final_dir,
            "170901_M00879_0087_000000000-AGEW9_analysis")
        self.assertTrue(os.path.exists(final_archive_dir))
        dirs = ("AB","CDE","logs","undetermined")
        for d in dirs:
            d = os.path.join(final_archive_dir,d)
            self.assertTrue(os.path.exists(d))
        # Check Fastqs have been updated for project 'AB'
        for fq in initial_fastqs:
            fastq = os.path.join(final_archive_dir,
                                 "AB",
                                 "fastqs",
                                 fq)
            self.assertFalse(os.path.exists(fastq))
        for fq in new_fastqs:
            fastq = os.path.join(final_archive_dir,
                                 "AB",
                                 "fastqs",
                                 fq)
            self.assertTrue(os.path.exists(fastq))

    def test_archive_excludes_extra_bcl2fastq_dirs(self):
        """archive: excludes extra bcl2fastq directories
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            top_dir=self.dirn)
        mockdir.create()
        # Add extra bcl2fastq directories
        for i,bcl2fastq_dir in enumerate(('bcl2fastq.1',
                                          'bcl2fastq.2')):
            mockdir1 = MockAnalysisDirFactory.bcl2fastq2(
                '170901_M00879_0087_000000000-AGEW9_%d' % i,
                'miseq',
                unaligned_dir=bcl2fastq_dir,
                metadata={ "instrument_datestamp": "170901" },
                top_dir=self.dirn)
            mockdir1.create()
            os.rename(os.path.join(mockdir1.dirn,bcl2fastq_dir),
                      os.path.join(mockdir.dirn,bcl2fastq_dir))
        # Make a mock archive directory
        archive_dir = os.path.join(self.dirn,"archive")
        final_dir = os.path.join(archive_dir,
                                 "2017",
                                 "miseq")
        os.makedirs(final_dir)
        self.assertTrue(os.path.isdir(final_dir))
        self.assertEqual(len(os.listdir(final_dir)),0)
        # Make autoprocess instance and set required metadata
        ap = AutoProcess(analysis_dir=mockdir.dirn,
                         settings=self.settings)
        ap.set_metadata("source","testing")
        ap.set_metadata("run_number","87")
        # Do archiving to staging
        status = archive(ap,
                         archive_dir=archive_dir,
                         year='2017',platform='miseq',
                         read_only_fastqs=True,
                         final=False)
        self.assertEqual(status,0)
        # Check staging dir
        staging_dir = os.path.join(
            final_dir,
            "__170901_M00879_0087_000000000-AGEW9_analysis.pending")
        self.assertTrue(os.path.exists(staging_dir))
        # Check bcl2fastq dirs not copied across
        for d in ('bcl2fastq',
                  'bcl2fastq.1',
                  'bcl2fastq.2',):
            self.assertFalse(os.path.exists(os.path.join(staging_dir,d)))
        # Do second archive operation to final
        status = archive(ap,
                         archive_dir=archive_dir,
                         year='2017',platform='miseq',
                         read_only_fastqs=True,
                         final=True)
        self.assertEqual(status,0)
        # Check final dir
        final_archive_dir = os.path.join(
            final_dir,
            "170901_M00879_0087_000000000-AGEW9_analysis")
        self.assertTrue(os.path.exists(final_archive_dir))
        # Check bcl2fastq dirs not copied across
        for d in ('bcl2fastq',
                  'bcl2fastq.1',
                  'bcl2fastq.2',):
            self.assertFalse(os.path.exists(os.path.join(final_archive_dir,d)))

    def test_archive_excludes_unwanted_content(self):
        """archive: excludes unwanted content
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            top_dir=self.dirn)
        mockdir.create()
        # Add unwanted content
        exclude_dirs = ("primary_data",
                        "save.bcl2fastqL1234",
                        "save.bcl2fastqL5678",
                        "tmp.illumina_qc.XYZ",
                        "__qc.XYZ123.tmp",)
        exclude_files = ("custom_SampleSheet.csv.bak",)
        for d in exclude_dirs:
            os.makedirs(os.path.join(mockdir.dirn,d))
            with open(os.path.join(mockdir.dirn,d,"placeholder"),'wt')  as fp:
                fp.write("")
        for f in exclude_files:
            with open(os.path.join(mockdir.dirn,f),'wt') as fp:
                fp.write("")
        # Make a mock archive directory
        archive_dir = os.path.join(self.dirn,"archive")
        final_dir = os.path.join(archive_dir,
                                 "2017",
                                 "miseq")
        os.makedirs(final_dir)
        self.assertTrue(os.path.isdir(final_dir))
        self.assertEqual(len(os.listdir(final_dir)),0)
        # Make autoprocess instance and set required metadata
        ap = AutoProcess(analysis_dir=mockdir.dirn,
                         settings=self.settings)
        ap.set_metadata("source","testing")
        ap.set_metadata("run_number","87")
        # Do archiving to staging
        status = archive(ap,
                         archive_dir=archive_dir,
                         year='2017',platform='miseq',
                         read_only_fastqs=True,
                         final=False)
        self.assertEqual(status,0)
        # Check staging dir
        staging_dir = os.path.join(
            final_dir,
            "__170901_M00879_0087_000000000-AGEW9_analysis.pending")
        self.assertTrue(os.path.exists(staging_dir))
        # Check unwanted content not copied across
        for d in exclude_dirs:
            self.assertFalse(os.path.exists(os.path.join(staging_dir,d)))
        for f in exclude_files:
            self.assertFalse(os.path.exists(os.path.join(staging_dir,f)))
        # Do second archive operation to final
        status = archive(ap,
                         archive_dir=archive_dir,
                         year='2017',platform='miseq',
                         read_only_fastqs=True,
                         final=True)
        self.assertEqual(status,0)
        # Check final dir
        final_archive_dir = os.path.join(
            final_dir,
            "170901_M00879_0087_000000000-AGEW9_analysis")
        self.assertTrue(os.path.exists(final_archive_dir))
        # Check unwanted content not copied across
        for d in exclude_dirs:
            self.assertFalse(os.path.exists(os.path.join(final_archive_dir,d)))
        for f in exclude_files:
            self.assertFalse(os.path.exists(os.path.join(final_archive_dir,f)))
