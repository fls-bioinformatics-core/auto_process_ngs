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

__version__ = '0.0.1'

#######################################################################
# Import modules that this module depends on
#######################################################################

import applications
import simple_scheduler
import JobRunner
import logging
import optparse
import sys
import os
import time

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
    if platform not in ('solid4','solid5500','hiseq','miseq','illumina-ga2x','other'):
        platform = None
        year = dirs[-2]
    else:
        year = dirs[-3]
    # Test for year
    if not year.isdigit or len(year) != 4:
        year = None
    return (year,platform)

def test_extract_year_and_platform():
    assert(extract_year_and_platform('/no/year/or/platform/data')==(None,None))
    assert(extract_year_and_platform('/year/2012/data')==('2012',None))
    assert(extract_year_and_platform('/year/2012/solid4/data')==('2012','solid4'))
    assert(extract_year_and_platform('/year/20121/solid4/data')==(None,'solid4'))
    assert(extract_year_and_platform('/year/2012/nonesense/data')==(None,None))

def get_archive_dir(archive_base,data_dirname):
    year,platform = extract_year_and_platform(data_dirname)
    archive_dir = archive_base
    if year is not None:
        archive_dir = os.path.join(archive_dir,year)
    if platform is not None:
        archive_dir = os.path.join(archive_dir,platform)
    return archive_dir

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

    # Dry run: just display where data dirs would be copied to
    if options.dry_run:
        for data_dir in data_dirs:
            print "%s -> %s" % (data_dir,get_archive_dir(archive_dir,data_dir))
        sys.exit()

    # Start scheduler
    s = simple_scheduler.SimpleScheduler()
    s.default_runner.log_dir(log_dir)
    s.start()
    # Make a job runner for computationally-intensive jobs
    if options.use_grid_engine:
        ci_runner = JobRunner.GEJobRunner(log_dir=log_dir)
    else:
        ci_runner = s.default_runner
    # Dictionary to store job references
    jobs = {}
    # Do copying for each directory
    for data_dir in data_dirs:
        # Schedule rsync job
        archive_to = get_archive_dir(archive_dir,data_dir)
        copy_name="copy.%s" % os.path.basename(data_dir)
        print "Starting copy of %s to %s" % (data_dir,archive_dir)
        job = s.submit(['data_manager.py','--copy-to=%s' % archive_to,data_dir],
                       name=copy_name)
        jobs[copy_name] = job
        # Schedule group name reset
        if new_group is not None:
            set_group_name="set_group.%s" % os.path.basename(data_dir)
            job = s.submit(['data_manager.py','--set-group=%s' % new_group,
                            os.path.join(archive_to,os.path.basename(data_dir))],
                           name=set_group_name,wait_for=(copy_name,))
            jobs[set_group_name] = job
        # Schedule MD5 check
        md5check_name="md5check.%s" % os.path.basename(data_dir)
        job = s.submit(['md5checker.py','-d',
                        data_dir,
                        os.path.join(archive_to,os.path.basename(data_dir))],
                       name=md5check_name,
                       wait_for=(copy_name,),
                       runner=ci_runner)
        jobs[md5check_name] = job
    while not s.is_empty():
        time.sleep(5)
    # Check MD5sums
    for data_dir in data_dirs:
        md5check_name="md5check.%s" % os.path.basename(data_dir)
        job = jobs[md5check_name]
        fp = open(job.log,'rU')
        skip_first_line = True
        for line in fp:
            if skip_first_line:
                skip_first_line = False
                continue
            elif line.strip('\n') == "Summary:":
                break
            elif not line.strip('\n').endswith('OK'):
                logging.error("Checksums failed for %s" % data_dir)
                break
        fp.close()
    print "Finished"
