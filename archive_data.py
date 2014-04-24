#!/bin/env python
#
#     archive_data.py: archiving utility for sequencing data
#     Copyright (C) University of Manchester 2014 Peter Briggs
#
########################################################################
#
# archive_data.py
#
#########################################################################


"""Perform archiving operations for sequencing data directories

"""

#######################################################################
# Module metadata
#######################################################################

__version__ = '0.1.0'

#######################################################################
# Import modules that this module depends on
#######################################################################

import platforms
import applications
import simple_scheduler
import JobRunner
import bcf_utils
import logging
import optparse
import sys
import os
import time

#######################################################################
# Classes
#######################################################################

class DataArchiver:
    # Class to handle archiving and verification

    def __init__(self,archive_dir,new_group=None,log_dir=None):
        # Top-level archive dir
        self._archive_dir = archive_dir
        # List of data directories
        self._data_dirs = []
        # Other options
        self._new_group = new_group
        self._log_dir = log_dir
        # Set up runner(s)
        self._runner = JobRunner.SimpleJobRunner(join_logs=True,log_dir=log_dir)
        # Set up scheduler
        self._sched = simple_scheduler.SimpleScheduler(runner=self._runner)

    def add_data_dir(self,data_dir):
        if data_dir not in self._data_dirs:
            self._data_dirs.append(data_dir)
        else:
            sys.stderr.write("%s already set up to be archived" % data_dir)

    def archive_dirs(self):
        # Check for data_manger
        if bcf_utils.find_program('data_manager.py') is None:
            raise Exception,"data_manager.py not found"
        # Start scheduler
        self._sched.start()
        # Set up archive jobs
        for data_dir in self._data_dirs:
            # Schedule copy
            archive_to = get_archive_dir(self._archive_dir,data_dir)
            copy_name="copy.%s" % os.path.basename(data_dir)
            print "Starting copy of %s to %s" % (data_dir,self._archive_dir)
            job = self._sched.submit(['data_manager.py','--copy-to=%s' % archive_to,data_dir],
                                     name=copy_name)
            # Schedule group name reset
            if self._new_group is not None:
                set_group_name="set_group.%s" % os.path.basename(data_dir)
                job = self._sched.submit(['data_manager.py','--set-group=%s' % self._new_group,
                                          os.path.join(archive_to,os.path.basename(data_dir))],
                                         name=set_group_name,wait_for=(copy_name,))

            # Run verification
            verify_name="verify.%s" % os.path.basename(data_dir)
            job = self._sched.submit(['data_manager.py','--verify=%s' % data_dir,
                                      os.path.join(archive_to,os.path.basename(data_dir))],
                                     name=verify_name,
                                     wait_for=(copy_name,),
                                     callbacks=(self.verification_complete,))
            # Final completion
            self._sched.callback("Finished",self.report_job_complete,
                                 wait_for=(verify_name,))
            # Wait for all jobs to complete
            self._sched.wait()

    def report_job_complete(self,name,jobs,sched):
        # Generic report completion of scheduled job(s)
        print "Job(s) completed:"
        for job in jobs:
            print "\t%s" % job.job_name

    def verification_complete(self,name,jobs,sched):
        # Callback handler for verification
        assert(len(jobs) == 1)
        job = jobs[0]
        print "Processing verification results from job %s" % job.job_name
        fp = open(job.log,'rU')
        for line in fp:
            if line.startswith('\t'):
                print line.strip()
        fp.close()

#######################################################################
# Functions
#######################################################################

def extract_year_and_platform(dirname):
    # Given directory path return the year and platform
    # if found and return as a tuple (year,platform)
    year = None
    platform = None
    # Split full path on OS separator
    dirs = os.path.abspath(dirname).split(os.sep)
    # Test for platform
    platform = dirs[-2]
    if platform not in platforms.list_platforms():
        platform = None
        year = dirs[-2]
    else:
        year = dirs[-3]
    # Test for year
    if not year.isdigit or len(year) != 4:
        year = None
    return (year,platform)

def get_archive_dir(archive_base,data_dirname):
    year,platform = extract_year_and_platform(data_dirname)
    archive_dir = archive_base
    if year is not None:
        archive_dir = os.path.join(archive_dir,year)
    if platform is not None:
        archive_dir = os.path.join(archive_dir,platform)
    return archive_dir

#######################################################################
# Tests
#######################################################################
import unittest
from mock_data import TestUtils,ExampleDirLanguages

class TestDataArchiver(unittest.TestCase):
    """Tests for DataArchiver class
    """
    def setUp(self):
        # Make a test data directory structure
        self.example_dir = ExampleDirLanguages()
        self.data_dir = self.example_dir.create_directory()
        self.archive_dir = TestUtils.make_dir()

    def tearDown(self):
        # Remove the test data directory (and copy)
        self.example_dir.delete_directory()
        TestUtils.remove_dir(self.archive_dir)

    def test_archive_data(self):
        """DataArchiver copies and verifies data directory

        """
        # Initial checks
        if bcf_utils.find_program('data_manager.py') is None:
            raise unittest.SkipTest("'data_manager.py' not found, test cannot run")
        dest_dir = os.path.join(self.archive_dir,os.path.basename(self.data_dir))
        self.assertFalse(os.path.exists(dest_dir))
        # Run the archiver
        archiver = DataArchiver(self.archive_dir)
        archiver.add_data_dir(self.data_dir)
        archiver.archive_dirs()
        # Check the copy
        self.assertTrue(os.path.exists(dest_dir))
        for dirpath,dirnames,filenames in os.walk(self.data_dir):
            for d in dirnames:
                d2 = os.path.join(dest_dir,
                                  os.path.relpath(os.path.join(dirpath,d),self.data_dir))
                self.assertTrue(os.path.isdir(d2))
            for f in filenames:
                f2 = os.path.join(dest_dir,
                                  os.path.relpath(os.path.join(dirpath,f),self.data_dir))
                self.assertTrue(os.path.isfile(f2))

class TestExtractYearAndPlatformFunction(unittest.TestCase):
    """Tests for test_extract_year_and_platform() function
    """
    def test_extract_year_and_platform(self):
        """extract_year_and_platform extracts year and platform data correctly
        """
        self.assertEqual(extract_year_and_platform('/no/year/or/platform/data'),(None,None))
        self.assertEqual(extract_year_and_platform('/year/2012/data'),('2012',None))
        self.assertEqual(extract_year_and_platform('/year/2012/solid4/data'),('2012','solid4'))
        self.assertEqual(extract_year_and_platform('/year/20121/solid4/data'),(None,'solid4'))
        self.assertEqual(extract_year_and_platform('/year/2012/nonesense/data'),(None,None))

#######################################################################
# Main program
#######################################################################

if __name__ == "__main__":

    # Process command line
    p = optparse.OptionParser(usage="%prog DATA_DIR [ DATA_DIR ... ] ARCHIVE_DIR",
                              version="%prog "+__version__,
                              description="Copy multiple sequence data directories "
                              "to ARCHIVE_DIR, then run MD5 checksums to verify.")
    p.add_option("--dry-run",action='store_true',dest='dry_run',
                 help="show where data would be copied but don't perform any operations")
    p.add_option("--log-dir",action='store',dest='log_dir',default=os.getcwd(),
                 help="write log files from archiving operations to LOG_DIR (defaults "
                 "to cwd)")
    p.add_option("--set-group",action='store',dest='new_group',default=None,
                 help="set group to NEW_GROUP on all copied files")
    p.add_option("--use-grid-engine",action='store_true',dest='use_grid_engine',default=False,
                 help="submit computationally-intensive jobs (e.g. MD5 checksumming) "
                 "to Grid Engine")
    p.add_option("--no-checks",action='store_true',dest="no_checks",default=False,
                 help="only copy, don't run MD5 checksums or link checks")
    p.add_option("--debug",action='store_true',dest='debug',
                 help="turn on debugging output (nb very verbose!)")
    options,args = p.parse_args()
    if len(args) < 2:
        p.error("Need to supply at least one DATA_DIR and an ARCHIVE_DIR")

    # Debugging output
    if options.debug:
        logging.getLogger().setLevel(logging.DEBUG)

    # Sort out directories and options
    data_dirs = args[:-1]
    archive_dir = os.path.abspath(args[-1])
    log_dir = os.path.abspath(options.log_dir)
    new_group = options.new_group
    print "Data dirs  : %s" % '\n\t'.join(data_dirs)
    print "Archive dir: %s" % archive_dir
    print "Log dir    : %s" % log_dir
    print "New group  : %s" % new_group
    print "Use GE     : %s" % options.use_grid_engine

    if options.dry_run:
        # Just print what would be archived where
        for data_dir in data_dirs:
            print "%s -> %s" % (data_dir,get_archive_dir(archive_dir,data_dir))
        sys.exit()

    # Make a job runner for computationally-intensive (CI) jobs
    ##if options.use_grid_engine:
    ##    ci_runner = JobRunner.GEJobRunner(log_dir=log_dir,
    ##                                      ge_extra_args=['-j','y'])
    ##else:
    ##    ci_runner = s.default_runner

    # Construct, populate and run archiver
    archiver = DataArchiver(archive_dir,new_group=new_group)
    for data_dir in data_dirs:
        archiver.add_data_dir(data_dir)
    archiver.archive_dirs()
    print "Archiving complete"
