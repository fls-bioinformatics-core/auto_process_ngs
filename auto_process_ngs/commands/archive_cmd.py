
#!/usr/bin/env python
#
#     archive_cmd.py: implement auto process archive command
#     Copyright (C) University of Manchester 2017-2018 Peter Briggs
#
#########################################################################

#######################################################################
# Imports
#######################################################################

import os
import time
import logging
import auto_process_ngs.applications as applications
import auto_process_ngs.fileops as fileops
import auto_process_ngs.simple_scheduler as simple_scheduler
import auto_process_ngs.tenx_genomics_utils as tenx_genomics_utils
from bcftbx.JobRunner import fetch_runner
from bcftbx.utils import format_file_size

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Command functions
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
      force (bool): if True then do archiving even if there are
        errors (e.g. key metadata items not set, permission error
        when setting group etc); otherwise abort archiving
        operation.
      dry_run (bool): report what would be done but don't
        perform any operations.

    Returns:
      UNIX-style integer returncode: 0 = successful termination,
        non-zero indicates an error occurred.
    """
    # Return value
    retval = 0
    # Check if analysis dir is actually staging directory
    analysis_dir = os.path.basename(ap.analysis_dir)
    is_staging = False
    if analysis_dir.startswith("__") and analysis_dir.endswith(".pending"):
        logger.warning("Operating directly on staged directory")
        if not final:
            raise Exception("Cannot re-stage already staged "
                            "analysis directory")
        else:
            is_staging = True
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
        datestamp = str(ap.metadata.instrument_datestamp)
        if len(datestamp) == 6:
            # Assume YYMMDD datestamp format
            year = "20%s" % datestamp[0:2]
        elif len(datestamp) == 8:
            # Assume YYYYMMDD datestamp format
            year = datestamp[0:4]
        else:
            raise Exception("Invalid datestamp '%s' (use "
                            "--year option)" % datestamp)
    archive_dir = os.path.join(archive_dir,year,platform)
    if not fileops.exists(archive_dir):
        raise OSError("Archive directory '%s' doesn't exist" %
                      archive_dir)
    # Determine target directory
    if not is_staging:
        final_dest = analysis_dir
        staging = "__%s.pending" % analysis_dir
    else:
        final_dest = analysis_dir[len("__"):-len(".pending")]
        staging = analysis_dir
    if final:
        dest = final_dest
    else:
        dest = staging
    print "Copying to archive directory: %s" % archive_dir
    print "Platform   : %s" % platform
    print "Year       : %s" % year
    print "Destination: %s %s" % (dest,
                                  "(final)" if final else
                                  "(staging)")
    # Check if final archive already exists
    if fileops.exists(os.path.join(archive_dir,final_dest)):
        raise Exception("Final archive already exists, stopping")
    # Report available space on target filesystem
    usage = fileops.disk_usage(archive_dir)
    print "Available  : %s/%s (%s%%)" % (format_file_size(usage.free),
                                         format_file_size(usage.total),
                                         usage.percent)
    # Check metadata
    check_metadata = ap.check_metadata(('source','run_number'))
    if not check_metadata:
        if not force or not is_staging:
            raise Exception("Some metadata items not set, stopping")
        logger.warning("Some metadata items not set, proceeding")
    if not is_staging:
        # Are there any projects to archive?
        projects = ap.get_analysis_projects()
        if not projects:
            raise Exception("No project directories found, nothing "
                            "to archive")
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
        excludes.append('--exclude=%s*' %
                        tenx_genomics_utils.flow_cell_id(ap.run_name))
        # Log dir
        log_dir = 'archive%s' % ('_final' if final else '_staging')
        if dry_run:
            log_dir += '_dry_run'
        ap.set_log_dir(ap.get_log_subdir(log_dir))
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
        # Keep track of jobs
        archiving_jobs = []
        # If making fastqs read-only then transfer them separately
        if read_only_fastqs and final:
            # Make sure excluded directories are excluded
            extra_options =  [ex for ex in excludes]
            # Set up to include only the fastq directories in
            # projects
            fastq_dirs = []
            for project in ap.get_analysis_projects():
                for fastq_dir in project.fastq_dirs:
                    fastq_dirs.append(os.path.join(
                        os.path.basename(project.dirn),
                        fastq_dir))
            # Update the extra options with includes/excludes
            extra_options.append('--include=*/')
            for fastq_dir in fastq_dirs:
                extra_options.append('--include=%s/**' % fastq_dir)
            extra_options.append('--exclude=*')
            # Execute the rsync
            rsync_fastqs = applications.general.rsync(
                "%s/" % ap.analysis_dir,
                os.path.join(archive_dir,staging),
                prune_empty_dirs=True,
                dry_run=dry_run,
                chmod='ugo-w',
                extra_options=extra_options)
            print "Running %s" % rsync_fastqs
            rsync_fastqs_job = sched.submit(rsync_fastqs,
                                            name="rsync.archive_fastqs")
            # Exclude fastqs from main rsync
            for fastq_dir in fastq_dirs:
                excludes.append('--exclude=%s' % fastq_dir)
            wait_for = [rsync_fastqs_job.job_name]
            # Add to list of jobs
            archiving_jobs.append(rsync_fastqs_job)
        else:
            # No separate Fastq rsync
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
        archiving_jobs.append(rsync_job)
        # Wait for jobs to complete
        rsync_job.wait()
        # Check exit status on jobs
        for job in archiving_jobs:
            print "%s completed: exit code %s" % (job.name,
                                                  job.exit_code)
        retval = sum([j.exit_code for j in archiving_jobs])
        if retval != 0:
            logger.warning("One or more archiving jobs failed "
                           "(non-zero exit code returned)")
        else:
            # Set the group
            if group is not None:
                print "Setting group of archived files to '%s'" % group
                if not dry_run:
                    set_group = fileops.set_group_command(
                        group,
                        os.path.join(archive_dir,staging),
                        safe=force,
                        verbose=True)
                    print "Running %s" % set_group
                    set_group_job = sched.submit(
                        set_group,
                        name="set_group.archive")
                    set_group_job.wait()
                    # Check exit status
                    exit_code = set_group_job.exit_code
                    print "%s completed: exit code %s" % (
                        set_group_job.name,
                        exit_code)
                    if exit_code != 0:
                        logger.warning("Setting group failed (non-zero "
                                       "exit status code returned)")
                    retval = retval + exit_code
        # Finish with scheduler
        sched.wait()
        sched.stop()
        # Bail out if there was a problem
        if retval != 0:
            if not force:
                raise Exception("Staging to archive failed")
            else:
                logger.warning("Staging to archive failed (ignored)")
    # Move to final location
    if final:
        print "Moving to final location: %s" % final_dest
        if not dry_run:
            fileops.rename(os.path.join(archive_dir,staging),
                           os.path.join(archive_dir,final_dest))
    # Report usage of target filesystem
    usage = fileops.disk_usage(archive_dir)
    print "Usage of archive: %s available (of %s) (%s%%)" % \
        (format_file_size(usage.free),
         format_file_size(usage.total),
         usage.percent)
    # Finish
    return retval
