#!/usr/bin/env python
#
#     archive_cmd.py: implement auto process archive command
#     Copyright (C) University of Manchester 2017-2024 Peter Briggs
#
#########################################################################

#######################################################################
# Imports
#######################################################################

import os
import time
import logging
from ..analysis import AnalysisDir
from ..metadata import AnalysisDirParameters
from .. import applications
from .. import fileops
from .. import simple_scheduler
from ..tenx.utils import flow_cell_id
from bcftbx.IlluminaData import IlluminaData
from bcftbx.utils import format_file_size
from bcftbx.utils import list_dirs

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
        from the auto_process.ini file).
      platform (str): set the value of the <PLATFORM> level in
        the archive (if not set then taken from the supplied
        autoprocessor instance).
      year (str): set the value of the <YEAR> level in the
        archive (if not set then defaults to the current year)
        (4 digits)
      perms (str): change the permissions of the destination
        files and directories according to the supplied
        argument (e.g. 'g+w') (if not set then use the value
         from the auto_process.ini file).
      group (str): set the group of the destination files to
        the supplied argument (if not set then use the value
        from the auto_process.ini file).
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
    print("Copying to archive directory: %s" % archive_dir)
    print("Platform   : %s" % platform)
    print("Year       : %s" % year)
    print("Destination: %s %s" % (dest,
                                  "(final)" if final else
                                  "(staging)"))
    # Check if final archive already exists
    if fileops.exists(os.path.join(archive_dir,final_dest)):
        raise Exception("Final archive already exists, stopping")
    # Report available space on target filesystem
    usage = fileops.disk_usage(archive_dir)
    print("Available  : %s/%s (%s%% in use)" %
          (format_file_size(usage.free),
           format_file_size(usage.total),
           usage.percent))
    # Check metadata
    check_metadata = ap.check_metadata(('source','run_number'))
    if not check_metadata:
        if not force or not is_staging:
            raise Exception("Some metadata items not set, stopping")
        logger.warning("Some metadata items not set, proceeding")
    # Locate extra bcl2fastq directories
    extra_bcl2fastq_dirs = list()
    for dirn in list_dirs(ap.analysis_dir):
        if dirn.endswith(".bak") or dirn.startswith("save."):
            # Ignore
            continue
        elif dirn == os.path.basename(ap.params.unaligned_dir):
            continue
        # Try to load data from the directory
        try:
            illumina_data = IlluminaData(ap.analysis_dir,
                                         unaligned_dir=dirn)
            extra_bcl2fastq_dirs.append(dirn)
        except Exception:
            pass
    if not is_staging:
        # Are there any projects to archive?
        try:
            projects = ap.get_analysis_projects()
        except Exception as ex:
            logging.warning("Error trying to fetch analysis projects: "
                            "%s" % ex)
            projects = []
        # Check projects
        empty_visium_dirs = False
        for project in projects:
            # Check for empty 'Visium_images' subdirs
            visium_images = os.path.join(project.dirn,"Visium_images")
            if os.path.isdir(visium_images):
                if not os.listdir(visium_images):
                    logger.warning("'%s': project contains an empty "
                                   "'Visium_images' subdirectory" %
                                   project.name)
                    empty_visium_dirs = True
        # Check for problems before final archiving
        if final and empty_visium_dirs:
            if not force:
                raise Exception("Empty 'Visium_images' subdirectories "
                                "detected in one or more projects; "
                                "either populate or remove (or use "
                                "--force)")
        if not projects:
            if not force:
                raise Exception("No project directories found, nothing "
                                "to archive")
            # Check if there is a bcl2fastq directory instead
            unaligned_dir = ap.params.unaligned_dir
            if not os.path.isabs(unaligned_dir):
                unaligned_dir = os.path.join(ap.analysis_dir,
                                             unaligned_dir)
            if os.path.exists(unaligned_dir):
                logging.warning("No project directories found, forcing "
                                "archiving of bcl2fastq output directory "
                                "'%s' instead" % ap.params.unaligned_dir)
                include_bcl2fastq = True
            else:
                raise Exception("No project directories or bcl2fastq "
                                "directory output found, nothing to "
                                "archive (even with --force)")
        # Determine which directories to exclude
        excludes = ['--exclude=primary_data',
                    '--exclude=save.*',
                    '--exclude=*.bak',
                    '--exclude=*.tmp',
                    '--exclude=tmp.*',
                    '--exclude=__*',]
        if not include_bcl2fastq:
            # Determine whether bcl2fastq dir should be included implicitly
            # because there are links from the analysis directories
            for project in projects:
                if project.fastqs_are_symlinks:
                    print("Found at least one project with fastq "
                          "symlinks (%s)" % project.name)
                    include_bcl2fastq = True
                    break
        if not include_bcl2fastq:
            print("Excluding '%s' directory from archive" %
                  ap.params.unaligned_dir)
            excludes.append('--exclude=%s' % ap.params.unaligned_dir)
        # Exclude extra bcl2fastq dirs
        for dirn in extra_bcl2fastq_dirs:
            print("Excluding '%s' directory from archive" % dirn)
            excludes.append('--exclude=%s' % dirn)
        # 10xgenomics products to exclude
        excludes.append('--exclude=*.mro')
        excludes.append('--exclude=%s*' % flow_cell_id(ap.run_name))
        # Log dir
        log_dir = 'archive%s' % ('_final' if final else '_staging')
        if dry_run:
            log_dir += '_dry_run'
        ap.set_log_dir(ap.get_log_subdir(log_dir))
        # Set up runners
        if runner is None:
            rsync_runner = ap.settings.runners.rsync
            default_runner = ap.settings.general.default_runner
        else:
            rsync_runner = runner
            default_runner = runner
        # Set log directory
        for r in (rsync_runner,
                  default_runner,):
            r.set_log_dir(ap.log_dir)
        # Setup a scheduler for multiple rsync jobs
        sched = simple_scheduler.SimpleScheduler(
            runner=default_runner,
            max_concurrent=ap.settings.general.max_concurrent_jobs,
            poll_interval=ap.settings.general.poll_interval)
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
            for project in projects:
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
                prune_empty_dirs=False,
                mirror=True,
                dry_run=dry_run,
                chmod='ugo-w',
                extra_options=extra_options)
            print("Running %s" % rsync_fastqs)
            rsync_fastqs_job = sched.submit(rsync_fastqs,
                                            name="rsync.archive_fastqs",
                                            runner=rsync_runner)
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
        print("Running %s" % rsync)
        rsync_job = sched.submit(rsync,name="rsync.archive",
                                 runner=rsync_runner,
                                 wait_for=wait_for)
        archiving_jobs.append(rsync_job)
        # Wait for jobs to complete
        rsync_job.wait()
        # Check exit status on jobs
        for job in archiving_jobs:
            print("%s completed: exit code %s" % (job.name,
                                                  job.exit_code))
        retval = sum([j.exit_code for j in archiving_jobs])
        if retval != 0:
            logger.warning("One or more archiving jobs failed "
                           "(non-zero exit code returned)")
        # Finish with scheduler
        sched.wait()
        sched.stop()
        # Bail out if there was a problem
        if retval != 0:
            if not force:
                raise Exception("Staging to archive area failed")
            else:
                logger.warning("Staging to archive area failed (ignored)")
    # Set the group
    if group is not None:
        print("Setting group of archived files to '%s'" % group)
        if not dry_run:
            # Setup a scheduler for file operations
            sched = simple_scheduler.SimpleScheduler(
                runner=default_runner,
                max_concurrent=ap.settings.general.max_concurrent_jobs,
                poll_interval=ap.settings.general.poll_interval)
            sched.start()
            # Set the group
            set_group = fileops.set_group_command(
                group,
                os.path.join(archive_dir,staging),
                safe=force,
                verbose=True)
            print("Running %s" % set_group)
            set_group_job = sched.submit(
                set_group,
                name="set_group.archive")
            set_group_job.wait()
            # Check exit status
            exit_code = set_group_job.exit_code
            print("%s completed: exit code %s" % (
                set_group_job.name,
                exit_code))
            if exit_code != 0:
                logger.warning("Setting group failed (non-zero "
                               "exit status code returned)")
                retval = retval + exit_code
            # Finish with scheduler
            sched.wait()
            sched.stop()
    # Perform final archiving operations
    if final:
        # Update the final stored Fastq paths and metadata
        # FIXME this is essentially duplicating functionality
        # FIXME in the 'update' command
        # FIXME (Also probably shouldn't update metadata for
        # FIXME 'dry_run' mode?)
        print("Updating stored paths and metadata")
        staged_analysis_dir = os.path.join(archive_dir,staging)
        archived_analysis_dir = os.path.join(archive_dir,final_dest)
        parameter_file = os.path.join(staged_analysis_dir,
                                      "auto_process.info")
        if os.path.exists(parameter_file):
            params = AnalysisDirParameters()
            params.load(parameter_file,strict=False)
            base_path = params.analysis_dir
            print("Stored base path: %s" % base_path)
            for p in ('analysis_dir',
                      'primary_data_dir',
                      'sample_sheet'):
                if not params[p]:
                    continue
                params[p] = os.path.normpath(
                    os.path.join(archived_analysis_dir,
                                 os.path.relpath(params[p],
                                                 base_path)))
                print("...updated '%s' (set to '%s')" % (p,params[p]))
            params.save()
        else:
            base_path = ap.analysis_dir
            logger.warning("Unable to get old base path from parameters")
            logger.warning("Using base path: %s (may be incorrect)" %
                           base_path)
        # Paths in QC info
        analysis_dir =  AnalysisDir(staged_analysis_dir)
        # FIXME AnalysisDir.get_projects method might not get all
        # FIXME the projects?
        for project in analysis_dir.get_projects():
            # FIXME should do all QC dirs (not just the primary one)
            qc_info = project.qc_info(project.qc_dir)
            if qc_info.fastq_dir:
                print("Project '%s': updating stored Fastq directory for QC" %
                      project.name)
                # FIXME could we just set it to the current Fastq path?
                new_fastq_dir = os.path.normpath(
                    os.path.join(archived_analysis_dir,
                                 os.path.relpath(qc_info.fastq_dir,
                                                 base_path)))
                print("...updated Fastq directory: %s" % new_fastq_dir)
                qc_info['fastq_dir'] = new_fastq_dir
                if not dry_run:
                    qc_info.save()
                # Bail out if there was a problem
                if retval != 0:
                    if not force:
                        raise Exception("Finalising archive failed")
                    else:
                        logger.warning("Finalising archive failed (ignored)")
        # Complete archiving
        print("Moving to final location: %s" % final_dest)
        if not dry_run:
            fileops.rename(os.path.join(archive_dir,staging),
                           os.path.join(archive_dir,final_dest))
    # Report usage of target filesystem
    usage = fileops.disk_usage(archive_dir)
    print("Usage of archive: %s available (of %s) (%s%% in use)" %
          (format_file_size(usage.free),
           format_file_size(usage.total),
           usage.percent))
    # Finish
    return retval
