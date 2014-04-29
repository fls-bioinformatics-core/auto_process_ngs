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

__version__ = '0.1.5'

#######################################################################
# Import modules that this module depends on
#######################################################################

import platforms
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
        self._sched = simple_scheduler.SimpleScheduler(runner=self._runner,
                                                       max_concurrent=4)

    def add_data_dir(self,dirn):
        data_dir = bcf_utils.AttributeDictionary()
        data_dir['name'] = os.path.basename(dirn)
        data_dir['dirn'] = os.path.abspath(dirn)
        data_dir['result'] = bcf_utils.AttributeDictionary(completed=1,
                                                           copy_status=None,
                                                           group_status=None,
                                                           verify_status=None)
        self._data_dirs.append(data_dir)

    def archive_dirs(self):
        # Check for data_manger
        if bcf_utils.find_program('data_manager.py') is None:
            raise Exception,"data_manager.py not found"
        # Start scheduler
        self._sched.start()
        # Set up archive jobs
        previous_copy = None
        for data_dir in self._data_dirs:
            # Schedule copy
            archive_to = get_archive_dir(self._archive_dir,data_dir.dirn)
            print "Setting up archiving of  %s to %s" % (data_dir.dirn,self._archive_dir)
            group = self._sched.group(data_dir.name,
                                      callbacks=(self.archiving_complete,))
            # Set up copy operation
            copy_name="copy.%s" % data_dir.name
            wait_for = []
            if previous_copy is not None:
                # Copy jobs should run sequentially
                wait_for.append(previous_copy)
                previous_copy = copy_name
            job = group.add(['data_manager.py','--copy-to=%s' % archive_to,data_dir.dirn],
                            name=copy_name,wait_for=wait_for)
            # Schedule group name reset
            if self._new_group is not None:
                set_group_name="set_group.%s" % data_dir.name
                job = group.add(['data_manager.py','--set-group=%s' % self._new_group,
                                 os.path.join(archive_to,data_dir.name)],
                                name=set_group_name,wait_for=(copy_name,))
                # Verify group name reset
                check_group_name="check_group.%s" % data_dir.name
                job = group.add(['data_manager.py','--check-group=%s' % self._new_group,
                                 os.path.join(archive_to,data_dir.name)],
                                name=check_group_name,wait_for=(set_group_name,))

            # Run verification
            verify_name="verify.%s" % data_dir.name
            job = group.add(['data_manager.py','--verify=%s' % data_dir.dirn,
                             os.path.join(archive_to,data_dir.name)],
                            name=verify_name,
                            wait_for=(copy_name,))
            group.close()
            # Final completion
            ##self._sched.callback("Finished",self.report_job_complete,
            ##                     wait_for=(verify_name,))
        # Wait for all jobs to complete
        self._sched.wait()

    def report_job_complete(self,name,jobs,sched):
        # Generic report completion of scheduled job(s)
        print "Job(s) completed:"
        for job in jobs:
            print "\t%s" % job.job_name

    def archiving_complete(self,name,groups,sched):
        # Callback handler for completion of archiving operations
        assert(len(groups) == 1)
        group = groups[0]
        # Locate the matching data dir
        name = group.group_name
        data_dir = None
        for d in self._data_dirs:
            if name == d.name:
                data_dir = d
        assert(data_dir is not None)
        data_dir.result['completed'] = 0
        # Check the operations
        print "Checking results from archiving operations for '%s':" % data_dir.dirn
        for job in group.jobs:
            print "\t%s" % job.job_name
            if job.job_name.startswith("copy."):
                # Check copy operation
                fp = open(job.log,'rU')
                for line in fp:
                    if line.startswith("Copy completed with status "):
                        data_dir.result['copy_status'] = int(line.split()[-1])
                        break
                    ##print line.strip()
                fp.close()
            elif job.job_name.startswith("check_group."):
                # Check group setting operation
                fp = open(job.log,'rU')
                for line in fp:
                    if line.startswith("Group/permissions check completed with status "):
                        data_dir.result['group_status'] = int(line.split()[-1])
                        break
                    ##print line.strip()
                fp.close()
            elif job.job_name.startswith("verify."):
                # Check verification operation
                fp = open(job.log,'rU')
                for line in fp:
                    if line.startswith("Verification completed with status "):
                        data_dir.result['verify_status'] = int(line.split()[-1])
                        break
                    ##print line.strip()
                fp.close()

    def result(self,data_dir):
        # Fetch the results for the specified data dir
        result = None
        for d in self._data_dirs:
            if d.dirn == data_dir or d.name == data_dir:
                return d.result
        # Didn't find the data dir
        raise KeyError,"Couldn't find '%s'" % data_dir

    def report(self,report_file=None):
        # Report the results
        if report_file is None:
            fp = sys.stdout
        else:
            print "Writing report to %s" % report_file
            fp = open(report_file,'w')
        for data_dir in self._data_dirs:
            fp.write("Data_dir: %s\n" % data_dir.name)
            status = 0
            for op_status in ('completed','copy_status','group_status','verify_status'):
                if data_dir.result[op_status] is not None:
                    fp.write("\t%s\t%s\n" % (op_status,
                                             'FAILED' if data_dir.result[op_status] else 'ok'))
                    if data_dir.result[op_status]: status = 1
            fp.write("Archiving status\t%s\n" % ('FAILED' if status else 'ok'))

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
import pwd
import grp
from mock_data import TestUtils,ExampleDirLanguages

class TestDataArchiver(unittest.TestCase):
    """Tests for DataArchiver class
    """
    def setUp(self):
        # Make a test data directory structure
        self.example_dir = ExampleDirLanguages()
        self.data_dir = self.example_dir.create_directory()
        self.archive_dir = TestUtils.make_dir()
        self.log_dir = TestUtils.make_dir()
        # Placeholders for additional data structures
        # These should be set in the tests where multiple
        # directories are required
        self.example_dir2 = None
        self.data_dir2 = None

    def tearDown(self):
        # Remove the test data directory (and copy)
        self.example_dir.delete_directory()
        if self.example_dir2 is not None:
            self.example_dir2.delete_directory()
        TestUtils.remove_dir(self.archive_dir)
        TestUtils.remove_dir(self.log_dir)

    def check_directory_contents(self,dir1,dir2):
        # Check that contents of dir1 are also in dir2
        # Checks that file and directory names match
        for dirpath,dirnames,filenames in os.walk(dir1):
            for d in dirnames:
                d2 = os.path.join(dir2,os.path.relpath(os.path.join(dirpath,d),dir1))
                self.assertTrue(os.path.isdir(d2))
            for f in filenames:
                f2 = os.path.join(dir2,os.path.relpath(os.path.join(dirpath,f),dir1))
                self.assertTrue(os.path.isfile(f2))

    def check_directory_group(self,dirn,group):
        # Check that contents of dirn belong to specified group
        gid = bcf_utils.get_gid_from_group(group)
        ##print "Checking: group %s GID %s" % (group,gid)
        for dirpath,dirnames,filenames in os.walk(dirn):
            for d in dirnames:
                d2 = os.path.join(dirpath,d)
                if os.path.islink(d2):
                    # Ignore links
                    continue
                ##print "%s/ %s" % (d2,bcf_utils.PathInfo(d2).gid)
                self.assertEqual(bcf_utils.PathInfo(d2).gid,gid)
            for f in filenames:
                f2 = os.path.join(dirpath,f)
                if os.path.islink(f2):
                    # Ignore links
                    continue
                ##print "%s %s" % (f2,bcf_utils.PathInfo(f2).gid)
                self.assertEqual(bcf_utils.PathInfo(f2).gid,gid)

    def test_archive_single_data_dir(self):
        """DataArchiver copies and verifies single data directory

        """
        # Initial checks
        dest_dir = os.path.join(self.archive_dir,os.path.basename(self.data_dir))
        self.assertFalse(os.path.exists(dest_dir))
        # Run the archiver
        archiver = DataArchiver(self.archive_dir,log_dir=self.log_dir)
        archiver.add_data_dir(self.data_dir)
        archiver.archive_dirs()
        # Check the copy
        self.assertTrue(os.path.exists(dest_dir))
        self.check_directory_contents(self.data_dir,dest_dir)
        # Check the status of each operation
        self.assertEqual(archiver.result(self.data_dir).completed,0)
        self.assertEqual(archiver.result(self.data_dir).copy_status,0)
        self.assertEqual(archiver.result(self.data_dir).group_status,None)
        self.assertEqual(archiver.result(self.data_dir).verify_status,0)

    def test_archive_multiple_data_dirs(self):
        """DataArchiver copies and verifies multiple data directories

        """
        # Make additional data dir
        self.example_dir2 = ExampleDirLanguages()
        self.data_dir2 = self.example_dir2.create_directory()
        # Initial checks
        dest_dir = os.path.join(self.archive_dir,os.path.basename(self.data_dir))
        self.assertFalse(os.path.exists(dest_dir))
        dest_dir2 = os.path.join(self.archive_dir,os.path.basename(self.data_dir2))
        self.assertFalse(os.path.exists(dest_dir2))
        # Run the archiver
        archiver = DataArchiver(self.archive_dir,log_dir=self.log_dir)
        archiver.add_data_dir(self.data_dir)
        archiver.add_data_dir(self.data_dir2)
        archiver.archive_dirs()
        # Check the copies
        self.assertTrue(os.path.exists(dest_dir))
        self.check_directory_contents(self.data_dir,dest_dir)
        self.assertTrue(os.path.exists(dest_dir2))
        self.check_directory_contents(self.data_dir,dest_dir2)
        # Check the status of each operation
        self.assertEqual(archiver.result(self.data_dir).completed,0)
        self.assertEqual(archiver.result(self.data_dir).copy_status,0)
        self.assertEqual(archiver.result(self.data_dir).group_status,None)
        self.assertEqual(archiver.result(self.data_dir).verify_status,0)
        self.assertEqual(archiver.result(self.data_dir2).completed,0)
        self.assertEqual(archiver.result(self.data_dir2).copy_status,0)
        self.assertEqual(archiver.result(self.data_dir2).group_status,None)
        self.assertEqual(archiver.result(self.data_dir2).verify_status,0)

    def test_archive_data_dir_set_group(self):
        """DataArchiver sets group when copying data directory

        """
        # Get a list of groups
        current_user = pwd.getpwuid(os.getuid()).pw_name
        groups = [g.gr_gid for g in grp.getgrall() if current_user in g.gr_mem]
        if len(groups) < 2:
            raise unittest.SkipTest("user '%s' must be in at least two groups" % current_user)
        # Get a second group
        gid = bcf_utils.PathInfo(self.archive_dir).gid
        new_gid = None
        for group in groups:
            if group != gid:
                new_gid = group
                break
        self.assertNotEqual(new_gid,gid)
        new_group = bcf_utils.get_group_from_gid(new_gid)
        ##print "Group %s GID %s" % (new_group,new_gid)
        # Initial checks
        dest_dir = os.path.join(self.archive_dir,os.path.basename(self.data_dir))
        self.assertFalse(os.path.exists(dest_dir))
        # Run the archiver
        archiver = DataArchiver(self.archive_dir,new_group=new_group,log_dir=self.log_dir)
        archiver.add_data_dir(self.data_dir)
        archiver.archive_dirs()
        # Check the copy
        self.assertTrue(os.path.exists(dest_dir))
        self.check_directory_contents(self.data_dir,dest_dir)
        self.check_directory_group(dest_dir,new_group)
        # Check the status of each operation
        self.assertEqual(archiver.result(self.data_dir).completed,0)
        self.assertEqual(archiver.result(self.data_dir).copy_status,0)
        self.assertEqual(archiver.result(self.data_dir).group_status,0)
        self.assertEqual(archiver.result(self.data_dir).verify_status,0)

    def test_archive_single_data_dir_with_unreadable_file(self):
        """DataArchiver reports failure when copying directory with unreadable file

        """
        # Remove read permission from a file in source dir
        os.chmod(self.example_dir.path('hello'),0244)
        # Initial checks
        dest_dir = os.path.join(self.archive_dir,os.path.basename(self.data_dir))
        self.assertFalse(os.path.exists(dest_dir))
        # Run the archiver
        archiver = DataArchiver(self.archive_dir,log_dir=self.log_dir)
        archiver.add_data_dir(self.data_dir)
        archiver.archive_dirs()
        # Check the copy
        self.assertTrue(os.path.exists(dest_dir))
        ##self.check_directory_contents(self.data_dir,dest_dir)
        # Check the status of each operation
        self.assertEqual(archiver.result(self.data_dir).completed,0)
        self.assertEqual(archiver.result(self.data_dir).copy_status,23)
        self.assertEqual(archiver.result(self.data_dir).group_status,None)
        self.assertEqual(archiver.result(self.data_dir).verify_status,1)

class TestExtractYearAndPlatformFunction(unittest.TestCase):
    """Tests for extract_year_and_platform() function
    """
    def test_extract_year_and_platform(self):
        """extract_year_and_platform extracts year and platform data correctly
        """
        self.assertEqual(extract_year_and_platform('/no/year/or/platform/data'),
                         (None,None))
        self.assertEqual(extract_year_and_platform('/year/2012/data'),
                         ('2012',None))
        self.assertEqual(extract_year_and_platform('/year/2012/solid4/data'),
                         ('2012','solid4'))
        self.assertEqual(extract_year_and_platform('/year/20121/solid4/data'),
                         (None,'solid4'))
        self.assertEqual(extract_year_and_platform('/year/2012/nonesense/data'),
                         (None,None))

class TestGetArchiveDirFunction(unittest.TestCase):
    """Tests for get_archive_dir() function
    """
    def test_get_archive_dir(self):
        """get_archive_dir returns correct archive directory
        """
        self.assertEqual(get_archive_dir('/mnt/archive','/no/year/or/platform/data_dir'),
                         '/mnt/archive')
        self.assertEqual(get_archive_dir('/mnt/archive','/year/2012/data_dir'),
                         '/mnt/archive/2012')
        self.assertEqual(get_archive_dir('/mnt/archive','/year/2012/solid4/data'),
                         '/mnt/archive/2012/solid4')
        self.assertEqual(get_archive_dir('/mnt/archive','/year/20121/solid4/data'),
                         '/mnt/archive/solid4')
        self.assertEqual(get_archive_dir('/mnt/archive','/year/2012/nonesense/data'),
                         '/mnt/archive')

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
    #p.add_option("--use-grid-engine",action='store_true',dest='use_grid_engine',default=False,
    #             help="submit computationally-intensive jobs (e.g. MD5 checksumming) "
    #             "to Grid Engine")
    p.add_option("--no-checks",action='store_true',dest="no_checks",default=False,
                 help="only copy, don't run MD5 checksums or link checks")
    p.add_option("--report",action='store',dest="report_file",default=None,
                 help="write final report to REPORT_FILE (otherwise write to stdout)")
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
    if options.report_file is not None:
        report_file = os.path.abspath(options.report_file)
    else:
        report_file = None
    print "Data dirs  : %s" % '\n\t'.join(data_dirs)
    print "Archive dir: %s" % archive_dir
    print "Log dir    : %s" % log_dir
    print "New group  : %s" % new_group
    print "Report file: %s" % report_file
    ##print "Use GE     : %s" % options.use_grid_engine

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
    archiver = DataArchiver(archive_dir,new_group=new_group,log_dir=log_dir)
    for data_dir in data_dirs:
        archiver.add_data_dir(data_dir)
    archiver.archive_dirs()
    archiver.report(report_file=report_file)
    print "Archiving complete"
