#######################################################################
# Tests for commands.py module
#######################################################################

import unittest
import tempfile
import shutil
import os
import pwd
import grp
from auto_process_ngs.auto_processor import AutoProcess
from auto_process_ngs.mock import MockAnalysisDirFactory
from auto_process_ngs.commands import archive

# Unit tests

class TestArchiveCommand(unittest.TestCase):
    """
    Tests for the 'archive' command
    """
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestArchiveCommand')
        # Store original location so we can get back at the end
        self.pwd = os.getcwd()
        # Move to working dir
        os.chdir(self.dirn)

    def tearDown(self):
        # Return to original dir
        os.chdir(self.pwd)
        # Remove the temporary test directory
        shutil.rmtree(self.dirn)

    def test_archive_to_staging(self):
        """archive: test copying to staging archive dir
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
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
        ap = AutoProcess(analysis_dir=mockdir.dirn)
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
        ap = AutoProcess(analysis_dir=mockdir.dirn)
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

    def test_archive_to_final(self):
        """archive: test copying to final archive dir
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
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
        ap = AutoProcess(analysis_dir=mockdir.dirn)
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

    def test_archive_to_final_via_staging(self):
        """archive: test copying to staging then final archive dir
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
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
        ap = AutoProcess(analysis_dir=mockdir.dirn)
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

    def test_archive_staging_to_final(self):
        """archive: test archiving directly from staging dir
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
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
        ap = AutoProcess(analysis_dir=mockdir.dirn)
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
