#!/usr/bin/env python
#
#     commands.py: implement auto process commands
#     Copyright (C) University of Manchester 2017 Peter Briggs
#
#########################################################################

#######################################################################
# Imports
#######################################################################

import os
import time
import applications
import simple_scheduler
import tenx_genomics_utils
import fileops
import logging
from bcftbx.JobRunner import fetch_runner

# Fetch configuration settings
import auto_process_ngs.settings
__settings = auto_process_ngs.settings.Settings()

# Module specific logger
logger = logging.getLogger(__name__)

import logging

#######################################################################
# Functions
#######################################################################

def archive(ap,archive_dir=None,platform=None,year=None,
            perms=None,group=None,include_bcl2fastq=False,
            read_only_fastqs=True,runner=None,
            final=False,force=False,dry_run=False):
    """
    Copy an analysis directory and contents to an archive area

    Copies the contents of the analysis directory to an archive
    area, which can be on a local or remote system.

    The archive directory is constructed in the form

    <TOP_DIR>/<YEAR>/<PLATFORM>/<DIR>/...

    The YEAR and PLATFORM can be overriden using the appropriate
    arguments.

    By default the data is copied to a 'staging' directory
    called '__ANALYSIS_DIR.pending' in the archive directory.
    The archiving can be finalised by setting the 'final'
    argumente to 'True', which performs a last update of the
    staging area before moving the data to its final location.

    Once the archive has been finalised any further archiving
    attempts will be refused.

    Copying of the data is performed using 'rsync'; multiple
    archive operations mirror the contents of the analysis
    directory (so any data removed from the source will also
    be removed from the archive).

    By default the 'bcl2fastq' directory is omitted from the
    archive, unless the fastq files in any projects are links to
    the data. Inclusion of this directory can be forced by
    setting the appropriate argument.

    The fastqs will be switched to be read-only in the archive
    by default.

    Arguments:
      ap (AutoProcessor): autoprocessor pointing to the
        analysis directory to be archived
      archive_dir (str): top level archive directory, of the
        form '[[user@]host:]dir' (if not set then use the value
        from the settings.ini file).
      platform (str): set the value of the <PLATFORM> level in
        the archive (if not set then taken from the supplied
        autoprocessor instance).
      year (str): set the value of the <YEAR> level in the
        archive (if not set then defaults to the current year)
        (4 digits)
      perms (str): change the permissions of the destination
        files and directories according to the supplied
        argument (e.g. 'g+w') (if not set then use the value
         from the settings.ini file).
      group (str): set the group of the destination files to
        the supplied argument (if not set then use the value
        from the settings.ini file).
      include_bcl2fastq (bool): if True then force inclusion
        of the 'bcl2fastq' subdirectory; otherwise only include
        it if fastq files in project subdirectories are symlinks.
      read_only_fastqs (bool): if True then make the fastqs
        read-only in the destination directory; otherwise keep
        the original permissions.
      runner: (optional) specify a non-default job runner to use
        for primary data rsync
      final (bool): if True then finalize the archive by
        moving the '.pending' temporary archive to the final
        location
      force (bool): if True then do archiving even if key
        metadata items are not set; otherwise abort archiving
        operation.
      dry_run (bool): report what would be done but don't
        perform any operations.

    Returns:
      UNIX-style integer returncode: 0 = successful termination,
        non-zero indicates an error occurred.
    """
    # Return value
    retval = 0
    # Fetch archive location
    if archive_dir is None:
        archive_dir = ap.settings.archive.dirn
    if archive_dir is None:
        raise Exception("No archive directory specified (use "
                        "--archive_dir option?)")
    # Construct subdirectory structure i.e. platform and year
    if platform is None:
        platform = ap.metadata.platform
    if platform is None:
        raise Exception("No platform specified (use --platform "
                        "option?)")
    if year is None:
        year = time.strftime("%Y")
    archive_dir = os.path.join(archive_dir,year,platform)
    # Determine target directory
    final_dest = os.path.basename(ap.analysis_dir)
    staging = "__%s.pending" % final_dest
    if final:
        dest = final_dest
    else:
        dest = staging
    print "Copying to archive directory: %s" % archive_dir
    print "Platform   : %s" % platform
    print "Year       : %s" % year
    print "Destination: %s%s" % (dest,
                                 " (final)" if final else "")
    # Check if final archive already exists
    if fileops_exists(os.path.join(archive_dir,final_dest)):
        logging.fatal("Final archive already exists, stopping")
        return 1
    # Are there any projects to?
    projects = ap.get_analysis_projects()
    if not projects:
        raise Exception("No project directories found, nothing "
                        "to archive")
    # Check metadata
    check_metadata = ap.check_metadata(('source','run_number'))
    if not check_metadata:
        if not force:
            logging.fatal("Some metadata items not set, stopping")
            return 1
        logging.warning("Some metadata items not set, proceeding")
    # Determine which directories to exclude
    excludes = ['--exclude=primary_data',
                '--exclude=save.*',
                '--exclude=*.bak',
                '--exclude=tmp.*']
    if not include_bcl2fastq:
        # Determine whether bcl2fastq dir should be included implicitly
        # because there are links from the analysis directories
        for project in projects:
            if project.fastqs_are_symlinks:
                print "Found at least one project with fastq " \
                    "symlinks (%s)" % project.name
                include_bcl2fastq = True
                break
    if not include_bcl2fastq:
        print "Excluding '%s' directory from archive" % \
            ap.params.unaligned_dir
        excludes.append('--exclude=%s' % ap.params.unaligned_dir)
    # 10xgenomics products to exclude
    excludes.append('--exclude=*.mro')
    excludes.append('--exclude="%s*"' %
                    tenx_genomics_utils.flow_cell_id(self.run_name))
    # Log dir
    ap.set_log_dir(ap.get_log_subdir('archive'))
    # Set up runner
    if runner is not None:
        runner = fetch_runner(runner)
    else:
        runner = ap.settings.runners.rsync
    runner.set_log_dir(ap.log_dir)
    # Setup a scheduler for multiple rsync jobs
    sched = simple_scheduler.SimpleScheduler(
        runner=runner,
        max_concurrent=ap.settings.general.max_concurrent_jobs)
    sched.start()
    # If making fastqs read-only then transfer them separately
    if read_only_fastqs and final:
        rsync_fastqs = applications.general.rsync(
            "%s/" % ap.analysis_dir,
            os.path.join(archive_dir,staging),
            prune_empty_dirs=True,
            dry_run=dry_run,
            chmod='ugo-w',
            extra_options=(
                '--include=*/',
                '--include=fastqs/**',
                '--exclude=*',))
        print "Running %s" % rsync_fastqs
        rsync_fastqs_job = sched.submit(rsync_fastqs,
                                        name="rsync.archive_fastqs")
        # Exclude fastqs from main rsync
        excludes.append('--exclude=fastqs')
        wait_for = [rsync_fastqs_job.job_name]
    else:
        rsync_fastqs_job = None
        wait_for = ()
    # Main rsync command
    rsync = applications.general.rsync(
        "%s/" % ap.analysis_dir,
        os.path.join(archive_dir,staging),
        prune_empty_dirs=True,
        mirror=True,
        dry_run=dry_run,
        chmod=perms,
        extra_options=excludes)
    print "Running %s" % rsync
    rsync_job = sched.submit(rsync,name="rsync.archive",
                             wait_for=wait_for)
    # Wait for scheduler to complete
    sched.wait()
    sched.stop()
    # Check exit status on job(s)
    if rsync_fastqs_job:
        exit_code = rsync_fastqs_job.exit_code
        print "rsync of FASTQs completed: exit code %s" % exit_code
        if exit_code != 0:
            logging.error("Failed to copy FASTQs to archive location "
                          "(non-zero exit code returned)")
    else:
        exit_code = 0
    print "rsync of data completed: exit code %s" % rsync_job.exit_code
    if rsync_job.exit_code != 0:
        logging.error("Failed to copy data to archive location "
                      "(non-zero exit code returned)")
    # Set returncode
    retval = exit_code if exit_code != 0 else rsync_job.exit_code
    # Set the group
    if group is not None:
        print "Setting group of archived files to '%s'" % group
        fileops.set_group(
            group,
            os.path.join(archive_dir,staging))
    # Move to final location
    if final:
        print "Moving to final location: %s" % final_dest
        fileops_rename(os.path.join(archive_dir,staging),
                       os.path.join(archive_dir,final_dest))
    # Finish
    return retval

def fileops_rename(src,dest):
    # Placeholder for renaming operation
    os.rename(src,dest)

def fileops_exists(path):
    # Placeholder for path existence checker
    return os.path.exists(path)
