#!/usr/bin/env python
#
#     cli/transfer_data.py: utility for copying data for sharing
#     Copyright (C) University of Manchester 2019-2024 Peter Briggs
#
#########################################################################
#
# transfer_data.py
#
#########################################################################

import sys
import os
import re
import argparse
import time
import shutil
from random import shuffle
from datetime import date
from fnmatch import fnmatch
from bcftbx.JobRunner import fetch_runner
from bcftbx.JobRunner import SimpleJobRunner
from bcftbx.utils import find_program
from bcftbx.utils import format_file_size
from bcftbx.utils import walk
from ..analysis import AnalysisDir
from ..analysis import AnalysisProject
from ..analysis import AnalysisFastq
from ..command import Command
from ..simple_scheduler import SimpleScheduler
from ..simple_scheduler import SchedulerReporter
from ..fileops import exists
from ..fileops import mkdir
from ..fileops import listdir
from ..fileops import copy_command
from ..fileops import set_permissions_command
from ..settings import Settings
from ..settings import get_config_dir
from ..utils import Location
from .. import get_version

# Logging
import logging
logging.basicConfig()
logger = logging.getLogger(__name__)

CELLRANGER_DIRS = ("cellranger_count",
                   "cellranger_multi")

def get_templates_dir():
    """
    Return location of 'templates' directory

    Returns the path to the 'templates' directory, or None if
    it doesn't exist.

    """
    path = get_config_dir()
    path = os.path.join(os.path.dirname(path),"templates")
    logging.debug("Putative templates dir: %s" % path)
    if os.path.isdir(path):
        return path
    else:
        return None


class TransferData:
    """
    Class for handling workflow for data transfer

    Wraps the "SimpleScheduler" class and provides the
    following methods:

    - run_job: runs a "Command" instance
    - wait: blocks until specified jobs (or all jobs)
      have completed
    - finish: checks the exit status of all jobs

    Arguments:
      runner (JobRunner): the default job runner to use
        when executing jobs
      working_dir (str): path to the working directory
      poll_interval (int): polling interval for
        scheduler
    """

    def __init__(self, runner, working_dir, poll_interval):
        print("Initialising job handler...")
        self._working_dir = working_dir
        self._runner = runner
        self._sched = SimpleScheduler(
            runner=self._runner,
            reporter=TransferDataSchedulerReporter(),
            poll_interval=poll_interval)
        self._sched.start()
        self._check_jobs = {}

    def run_job(self, name, cmd, runner=None, wait_for=None):
        """
        Submits the specified command to the scheduler

        Arguments:
          name (str): name for the command
          cmd (Command): the command to run
          runner (JobRunner): optional, job runner to use
          wait_for (list): optional, list of jobs to wait for
            before scheduling this command

        Returns:
          SchedulerJob: a SchedulerJob instance for the submitted
            job.
        """
        print(f"Running '{cmd}'")
        if runner is None:
            runner = self._runner
        if wait_for:
            wait_for_ = []
            for job in wait_for:
                wait_for_.append(job.job_name)
            wait_for = wait_for_
        job = self._sched.submit(cmd.command_line,
                                 name=name,
                                 runner=runner,
                                 wd=self._working_dir,
                                 wait_for=wait_for)
        self._check_jobs[job.name] = job
        return job

    def wait(self, *jobs):
        """
        Block until specified jobs have completed.

        If one or more jobs are supplied then wait for those
        jobs to complete before returning; otherwise wait for
        all jobs currently scheduled to complete.

        Arguments:
          jobs (list): optional list of SchedulerJobs
            to wait for (otherwise wait for all jobs)
        """
        if not jobs:
            print("Waiting for all scheduled jobs to complete...")
            self._sched.wait()
        else:
            print(f"Waiting for jobs:")
            for job in jobs:
                print(f"...{job.job_name}")
                job.wait()
        print("All jobs completed")

    def finish(self):
        """
        Wait for all jobs then check and return status

        Returns:
          Integer: status code (0 for success, 1 if any
            of the jobs failed).
        """
        print("Finishing...")
        self._sched.wait()
        status = 0
        for name in self._check_jobs:
            if self._check_jobs[name].exit_code != 0:
                logger.warning(f"'{name}': job failed")
                status = 1
        if status != 0:
            logger.error("some transfer operations did not complete "
                         "successfully (see warnings above)")
        return status


class TransferDataSchedulerReporter(SchedulerReporter):
    """
    Custom reporter for scheduler
    """
    def __init__(self):
        SchedulerReporter.__init__(
            self,
            scheduler_status="[%(time_stamp)s] %(n_running)d running, "
            "%(n_waiting)d waiting, %(n_finished)d finished",
            job_scheduled="[%(time_stamp)s] scheduled job '%(job_name)s' "
            "#%(job_number)d",
            job_start="[%(time_stamp)s] started job '%(job_name)s' "
            "#%(job_number)d (%(job_id)s)",
            job_end="[%(time_stamp)s] completed job '%(job_name)s' "
            "#%(job_number)d (%(job_id)s)"
        )

# Main function

def main(argv=None):
    """
    Run the 'transfer_data' utility

    Arguments:
      argv (list): optional, command line arguments to
        process (otherwise take arguments from
        'sys.argv')

    Returns:
      Integer: 0 on success, 1 on failure.
    """
    # Command line arguments
    if argv is None:
        argv = sys.argv[1:]

    # Load configuration
    settings = Settings()

    # Collect defaults
    default_runner = settings.runners.rsync

    # Get pre-defined destinations
    destinations = [name for name in settings.destination]

    # Command line
    p = argparse.ArgumentParser(
        description="Transfer copies of Fastq data from an analysis "
        "project to an arbitrary destination for sharing with other "
        "people")
    p.add_argument('--version', action='version',
                   version=("%%(prog)s %s" % get_version()))
    p.add_argument('--subdir',action='store',
                   choices=('random_bin','run_id'),
                   default=None,
                   help="subdirectory naming scheme: 'random_bin' "
                   "locates a random pre-existing empty subdirectory "
                   "under the target directory; 'run_id' creates a "
                   "new subdirectory "
                   "'PLATFORM_DATESTAMP.RUN_ID-PROJECT'. If this "
                   "option is not set then no subdirectory will be "
                   "used")
    sp = p.add_argument_group("Fastq selection")
    sp.add_argument('--samples',action='store',dest='sample_list',
                    default=None,
                    help="list of names of samples to transfer")
    sp.add_argument('--filter',action='store',dest='filter_pattern',
                    default=None,
                    help="filter Fastq file names based on PATTERN")
    sp.add_argument('--no_fastqs',action='store_true',
                    help="don't copy Fastqs (other artefacts will be "
                    "copied, if specified)")
    sp = p.add_argument_group("ZIP file archives")
    sp.add_argument('--zip_fastqs',action='store_true',
                    help="put Fastqs into a ZIP file")
    sp.add_argument('--max_zip_size',action='store',dest='max_zip_size',
                    default=None,
                    help="when using '--zip_fastqs' option, defines the "
                    "maximum size for the output zip file; multiple zip "
                    "files will be created if the data exceeds this "
                    "limit (default is create a single zip file with no "
                    "size limit)")
    sp = p.add_argument_group("README generation")
    sp.add_argument('--readme',action='store',
                    metavar='README_TEMPLATE',dest='readme_template',
                    help="template file to generate README file from; "
                    "can be full path to a template file, or the name "
                    "of a file in the 'templates' directory")
    sp.add_argument('--weburl',action='store',
                    help="base URL for webserver (sets the value of "
                    "the WEBURL variable in the template README)")
    sp = p.add_argument_group("Additional artefacts")
    sp.add_argument('--include_qc_report',action='store_true',
                    help="copy the zipped QC reports to the final "
                    "location")
    sp.add_argument('--include_10x_outputs',action='store_true',
                    help="copy outputs from 10xGenomics pipelines (e.g. "
                    "'cellranger count') to the final location")
    sp.add_argument('--include_cloupe_files',action='store_true',
                    help="copy .cloupe files output from 10xGenomics "
                    "pipelines to the final location")
    sp.add_argument('--include_visium_images', action='store_true',
                    help="copy images for 10x Genomics Visium projects "
                    "to the final location")
    sp.add_argument('--include_downloader',action='store_true',
                   help="copy the 'download_fastqs.py' utility to the "
                    "final location")
    sp = p.add_argument_group("Advanced options")
    sp.add_argument('--link',action='store_true',
                    help="hard link files instead of copying")
    sp.add_argument('--runner',action='store',
                    help="specify the job runner to use for executing "
                    "the checksumming, Fastq copy and tar gzipping "
                    "operations (defaults to job runner defined for "
                    "copying in config file [%s])" % default_runner)
    p.add_argument('dest',action='store',metavar="DEST",
                   help="destination to copy Fastqs to; can be the "
                   "name of a destination defined in the configuration "
                   "file, or an arbitrary location of the form "
                   "'[[USER@]HOST:]DIR' (%s)" %
                   (("available destinations: %s" %
                     (','.join("'%s'" % d for d in sorted(destinations))))
                    if destinations else
                    "no destinations currently defined"))
    p.add_argument('project',action='store',
                   metavar="PROJECT",
                   help="path to project directory (or to a Fastqs "
                   "subdirectory in a project) to copy Fastqs from")

    # Process command line
    args = p.parse_args(argv)

    # Flag for Fastq transfer
    include_fastqs = not args.no_fastqs

    # List of samples to include
    if args.sample_list is not None:
        include_samples = str(args.sample_list).split(',')
    else:
        include_samples = None

    # Check if target is pre-defined destination
    if args.dest in destinations:
        print("Loading settings for destination '%s'" % args.dest)
        dest = settings.destination[args.dest]
        target_dir = dest.directory
        readme_template = dest.readme_template
        subdir = dest.subdir
        zip_fastqs = dest.zip_fastqs
        max_zip_size = dest.max_zip_size
        include_downloader = dest.include_downloader
        include_qc_report = dest.include_qc_report
        hard_links = dest.hard_links
        weburl = dest.url
    else:
        target_dir = args.dest
        readme_template = None
        subdir = None
        zip_fastqs = False
        max_zip_size = None
        include_downloader = False
        include_qc_report = False
        hard_links = False
        weburl = None

    # Update defaults with command line values
    if args.readme_template:
        readme_template = args.readme_template
    if args.subdir:
        subdir = args.subdir
    if args.zip_fastqs:
        zip_fastqs = True
    if args.max_zip_size:
        max_zip_size = args.max_zip_size
    if args.include_downloader:
        include_downloader = True
    if args.include_qc_report:
        include_qc_report = True
    if args.weburl:
        weburl = args.weburl
    if args.link:
        hard_links = args.link

    # Additional artefacts
    include_10x_outputs = bool(args.include_10x_outputs)
    include_cloupe_files = bool(args.include_cloupe_files)
    include_visium_images = bool(args.include_visium_images)

    # Report settings
    print("================= Transfer settings =================")
    print(f"Source project       : {args.project}")
    print(f"Target dir           : {target_dir}")
    print(f"Subdir               : {subdir}")
    print(f"Web URL              : {weburl}")
    print(f"README template      : {readme_template}")
    print(f"Include Fastqs       : {include_fastqs}")
    print(f"Include QC report    : {include_qc_report}")
    print(f"Include 10x outputs  : {include_10x_outputs}")
    print(f"Include .cloupe files: {include_visium_images}")
    print(f"Include Visium images: {include_visium_images}")
    print(f"Include downloader   : {include_downloader}")
    print(f"Hard link Fastqs     : {hard_links}")
    print(f"Zip Fastqs           : {zip_fastqs}")
    print(f"Max ZIP size         : {max_zip_size}")

    # Check at least one artefact is being transferred
    if not (include_fastqs or
            include_downloader or
            include_qc_report or
            include_10x_outputs or
            include_visium_images or
            readme_template):
        logger.error("No artefacts specified for transfer")
        return 1

    # Sort out project directory
    project = AnalysisProject(args.project)
    if not project.is_analysis_dir:
        # Assume it's the Fastq dir
        fastq_dir = os.path.basename(args.project)
        project = AnalysisProject(os.path.dirname(args.project))
    else:
        fastq_dir = None
    if not project.is_analysis_dir:
        logger.error("'%s': project not found" % args.project)
        return 1
    project_name = project.name

    # Parent analysis directory
    analysis_dir = AnalysisDir(os.path.dirname(project.dirn))

    # Fastqs directory
    try:
        project.use_fastq_dir(fastq_dir)
    except Exception as ex:
        logger.error("'%s': failed to load Fastq set '%s': %s" %
                     (project.name,fastq_dir,ex))
        return 1

    # Examine source project
    print("================= Source project =================")
    print("Transferring data from '%s' (%s)" % (project.name,
                                                project.dirn))

    # Locate Fastqs
    print("Looking for Fastqs")
    print(f"...examining {project.fastq_dir}")
    samples = set()
    nfastqs = 0
    nindex_fastqs = 0
    lanes = set()
    fsize = 0
    for sample in project.samples:
        fqs = []
        for fq in sample.fastq:
            if include_samples and not sample.name in include_samples:
                # Sample not in list to include
                continue
            elif args.filter_pattern and \
                 not fnmatch(os.path.basename(fq),args.filter_pattern):
                # Filter pattern specified but Fastq
                # doesn't match so skip
                continue
            fqs.append(fq)
            fq_attrs = AnalysisFastq(fq)
            if fq_attrs.lane_number:
                lanes.add(fq_attrs.lane_number)
            if fq_attrs.is_index_read:
                nindex_fastqs += 1
        if fqs:
            samples.add(sample.name)
            for fq in fqs:
                fsize += os.lstat(fq).st_size
                nfastqs += 1
    nsamples = len(samples)
    if nfastqs:
        print(f"...found {nfastqs} Fastqs ({nsamples} samples, "
              f"{len(lanes)} lanes, {nindex_fastqs} index reads) "
              f"[{format_file_size(fsize)}]")

    # QC reports
    print("Locating for QC reports")
    qc_zips = list()
    # Check QC directories and look for zipped reports
    for qc_dir in project.qc_dirs:
        # Get the associated Fastq set
        # NB only compare the basename of the Fastq dir
        # in case full paths weren't updated
        fq_set = os.path.basename(project.qc_info(qc_dir).fastq_dir)
        if fq_set == os.path.basename(project.fastq_dir):
            for qc_base in ("%s_report.%s.%s" %
                            (qc_dir,project.name,project.info.run),
                            "%s_report.%s.%s" %
                            (qc_dir,project.name,
                             os.path.basename(
                                 analysis_dir.analysis_dir)),):
                qc_zip = os.path.join(project.dirn,
                                      "%s.zip" % qc_base)
                if os.path.exists(qc_zip):
                    qc_zips.append(qc_zip)
                    print(f"...found {qc_zip}")

    # 10xGenomics outputs
    print("Looking for 10xGenomics pipeline outputs")
    cellranger_dirs = list()
    for d in CELLRANGER_DIRS:
        cellranger_dir = os.path.join(project.dirn,d)
        if os.path.isdir(cellranger_dir):
            cellranger_dirs.append(cellranger_dir)
            print(f"...found {cellranger_dir}")
    fsize_10x_outputs = 0
    for cellranger_dir in cellranger_dirs:
        for f in walk(cellranger_dir):
            fsize_10x_outputs += os.lstat(f).st_size
    if cellranger_dirs:
        print(f"...found {len(cellranger_dirs)} 10x pipeline outputs "
              f"[{format_file_size(fsize_10x_outputs)}]")

    # Locate 10xGenomics .cloupe files
    print("Looking for 10xGenomics .cloupe files")
    cloupe_files = list()
    for cellranger_dir in CELLRANGER_DIRS:
        cellranger_dir = os.path.join(project.dirn, cellranger_dir)
        if not os.path.isdir(cellranger_dir):
            continue
        for f in walk(cellranger_dir):
            if f.endswith(".cloupe") and f.split(os.sep)[-2] == "outs":
                cloupe_files.append(f)
                print(f"...found {f}")
    if cloupe_files:
        print("...found {len(cloupe_files)} .cloupe files")

    # Locate Visium images
    print("Looking for Visium images")
    visium_images_dir = os.path.join(project.dirn, "Visium_images")
    fsize_visium_images = 9
    if os.path.isdir(visium_images_dir) and \
       len(os.listdir(visium_images_dir)) > 0:
        print(f"...found {visium_images_dir}")
        for f in walk(visium_images_dir):
            fsize_visium_images += os.lstat(f).st_size
        print(f"...found Visium images "
              f"[{format_file_size(fsize_visium_images)}]")
    else:
        visium_images_dir = None

    # Check target dir
    print("================= Target location =================")

    # Check "base" target dir exists
    if not Location(target_dir).is_remote:
        target_dir = os.path.abspath(target_dir)
    if not exists(target_dir):
        print("'%s': base target directory not found" % target_dir)
        return
    else:
        print("Base target directory %s" % target_dir)

    # Check if hard links are possible
    if hard_links:
        if Location(target_dir).is_remote:
            logger.error("'%s': hard links requested but target directory "
                         "is on a remote filesystem" % target_dir)
            return 1
        else:
            if os.lstat(project.fastq_dir).st_dev != \
               os.lstat(target_dir).st_dev:
                logger.error("'%s': hard links requested but target "
                             "directory is on a different filesystem" %
                             target_dir)
                return 1
        print("Hard linking is ok")

    # Determine subdirectory
    if subdir == "random_bin":
        # Find a random empty directory under the
        # target directory
        print("Locating random empty bin")
        subdirs = [d for d in os.listdir(target_dir)
                   if os.path.isdir(os.path.join(target_dir,d))]
        if not subdirs:
            print("Failed to locate subdirectories")
            return
        shuffle(subdirs)
        subdir = None
        for d in subdirs:
            if not os.listdir(os.path.join(target_dir,d)):
                # Empty bin
                subdir = d
                break
        if subdir is None:
            print("Failed to locate empty subdirectory")
            return
        print("Using subdirectory '%s'" % subdir)
        # Update target dir
        target_dir = os.path.join(target_dir, subdir)
    elif subdir == "run_id":
        # Construct subdirectory name based on the
        # run ID
        run_number = str(analysis_dir.metadata.run_number)
        if analysis_dir.metadata.analysis_number is not None:
            run_number += "_" + str(analysis_dir.metadata.analysis_number)
        subdir = "{platform}_{datestamp}.{run_number}-{project}".format(
            platform=analysis_dir.metadata.platform.upper(),
            datestamp=analysis_dir.metadata.instrument_datestamp,
            run_number=run_number,
            project=project.name)
        print("Using subdirectory '%s'" % subdir)
        # Update target dir
        target_dir = os.path.join(target_dir, subdir)
        # Check it doesn't already exist
        if exists(target_dir):
            logger.error("'%s': directory already exists" % target_dir)
            return
    print(f"Final target directory: {target_dir}")

    # Check artefacts for inclusion
    print("Checking artefacts requested for inclusion")

    # Include Fastqs
    if include_fastqs:
        if nfastqs:
            print("...Fastq files")
        else:
            logger.error("No Fastqs found")
            return 1

    # Include zipped QC reports
    if include_qc_report:
        if qc_zips:
            print("...QC reports")
        else:
            logger.error("No zipped QC reports found")
            return 1
    else:
        qc_zips = None

    # Include 10xGenomics outputs
    if include_10x_outputs:
        if cellranger_dirs:
            print("...10xGenomics pipeline outputs")
        else:
            logger.error("No outputs from 10xGenomics pipelines found")
            return 1
    else:
        cellranger_dirs = None

    # Include 10xGenomics .cloupe files
    if include_cloupe_files:
        if cloupe_files:
            print("...10xGenomics .cloupe files")
        else:
            logger.error("No .cloupe files found")
            return 1
    else:
        cloupe_files = None

    # Include Visium images
    if include_visium_images:
        if visium_images_dir:
            print("...10xGenomics Visium images")
        else:
            logger.error("No Visium images found")
            return 1
    else:
        visium_images_dir = None

    # README template
    if readme_template:
        # Check that template file exists
        print("Locating README template")
        template = None
        for filen in (readme_template,
                      os.path.join(get_templates_dir(),
                                   readme_template),):
            if os.path.exists(filen):
                template = filen
                break
        if template is None:
            logger.error("'%s': template file not found" %
                         readme_template)
            return 1
        else:
            readme_template = template
        print("... found %s" % readme_template)

    # Locate downloader
    if include_downloader:
        print("Locating downloader script")
        downloader = find_program("download_fastqs.py")
        if downloader is None:
            logging.error("Unable to locate download_fastqs.py")
            return 1
        print("... found %s" % downloader)
    else:
        downloader = None

    # Transfer the data
    print("================= Transferring data =================")

    # Get runner for copy job
    if args.runner:
        runner = fetch_runner(args.runner)
    else:
        runner = default_runner

    # Set identifier for jobs
    job_id = "%s%s" % (project_name,
                       (".%s" % fastq_dir
                        if fastq_dir is not None
                        else ''))

    # Set the working directory
    working_dir = os.path.abspath("transfer.%s.%s" % (job_id,
                                                      int(time.time())))
    mkdir(working_dir)
    print("Created working dir %s" % working_dir)

    # Make target directory
    if not exists(target_dir):
        mkdir(target_dir)

    # Construct the README
    if readme_template:
        # Read in template
        with open(readme_template,'rt') as fp:
            readme = fp.read()
        # Substitute template variables
        template_vars = {
            'PLATFORM': analysis_dir.metadata.platform.upper(),
            'RUN_NUMBER': analysis_dir.metadata.run_number,
            'DATESTAMP': analysis_dir.metadata.instrument_datestamp,
            'PROJECT': project_name,
            'WEBURL': weburl,
            'BIN': subdir,
            'DIR': target_dir,
            'TODAY': date.today().strftime("%d/%m/%Y"),
        }
        for var in template_vars:
            value = template_vars[var]
            if value is None:
                value = '?'
            else:
                value = str(value)
            readme = re.sub(r"%{var}%".format(var=var),
                            value,
                            readme)
        # Write out a temporary README file
        readme_file = os.path.join(working_dir,"README")
        with open(readme_file,'wt') as fp:
            fp.write(readme)
    else:
        # No README
        readme_file = None

    # Make a job handler
    td = TransferData(runner, working_dir,
                      settings.general.poll_interval)

    # Transfer Fastqs
    if include_fastqs:
        if not zip_fastqs:
            # Build command to run manage_fastqs.py to copy Fastqs
            copy_cmd = Command("manage_fastqs.py")
            if args.sample_list:
                copy_cmd.add_args("--samples",args.sample_list)
            if args.filter_pattern:
                copy_cmd.add_args("--filter",args.filter_pattern)
            if hard_links:
                copy_cmd.add_args("--link")
            if fastq_dir is not None:
                copy_cmd.add_args("--fastq_dir",
                                  fastq_dir)
            copy_cmd.add_args(analysis_dir.analysis_dir,
                              project_name,
                              "copy",
                              target_dir)
            td.run_job(f"copy_fastqs.{job_id}", copy_cmd)
        else:
            # Build command to zip Fastqs
            zip_cmd = Command("manage_fastqs.py")
            if args.sample_list:
                zip_cmd.add_args("--samples",args.sample_list)
            if args.filter_pattern:
                zip_cmd.add_args("--filter",args.filter_pattern)
            if max_zip_size:
                zip_cmd.add_args("--max_zip_size",max_zip_size)
            if fastq_dir is not None:
                zip_cmd.add_args("--fastq_dir",
                                 fastq_dir)
            zip_cmd.add_args(analysis_dir.analysis_dir,
                             project_name,
                             "zip")
            zip_job = td.run_job(f"zip_fastqs.{job_id}", zip_cmd)
            # Wait for ZIP files to be generated
            td.wait(zip_job)
            # Rename ZIP file(s) and move to final location
            run_number = str(analysis_dir.metadata.run_number)
            if analysis_dir.metadata.analysis_number is not None:
                run_number += "_" + str(analysis_dir.metadata.analysis_number)
            final_zip_basename = \
                "{platform}_{datestamp}.{run_number}-{project}-fastqs".\
                format(
                    platform=analysis_dir.metadata.platform.upper(),
                    datestamp=analysis_dir.metadata.instrument_datestamp,
                    run_number=run_number,
                    project=project.name)
            job_ix = 0
            for f in listdir(working_dir):
                if f == "%s.chksums" % project_name:
                    # Assume it's the checksum file
                    # Copy to final location
                    td.run_job(f"copy_checksums.{job_id}",
                               copy_command(
                                   os.path.join(working_dir,f),
                                   os.path.join(target_dir,
                                                "%s.checksums" %
                                                final_zip_basename)),
                               runner=SimpleJobRunner())
                elif f.endswith(".zip") and \
                   f.startswith("%s." % project_name):
                    # Assume it's ZIP output from packaging process
                    final_zip = "%s%s" % (final_zip_basename,
                                          f[len(project_name):])
                    # Copy to final location
                    job_ix += 1
                    td.run_job(f"copy_zipped_fastqs.{job_ix}.{job_id}",
                               copy_command(
                                   os.path.join(working_dir, f),
                                   os.path.join(target_dir, final_zip)))

    # Copy README
    if readme_file is not None:
        print("Copying README file")
        td.run_job(f"copy_readme.{job_id}",
                   copy_command(readme_file,
                                os.path.join(target_dir,"README")),
                   runner=SimpleJobRunner())

    # Copy download_fastqs.py
    if downloader:
        print("Copying downloader")
        td.run_job(f"copy_downloader.{job_id}",
                   copy_command(downloader,
                                os.path.join(target_dir,
                                             os.path.basename(downloader))),
                   runner=SimpleJobRunner())

    # Copy QC reports
    if qc_zips:
        for qc_zip in qc_zips:
            print("Copying '%s'" % os.path.basename(qc_zip))
            td.run_job(
                f"copy_qc_zip.{job_id}.{os.path.basename(qc_zip)}",
                copy_command(qc_zip,
                             os.path.join(target_dir,
                                          os.path.basename(qc_zip)),
                             link=hard_links),
                runner=SimpleJobRunner())

    # Tar and copy 10xGenomics outputs
    if cellranger_dirs:
        for cellranger_dir in cellranger_dirs:
            print("Tar gzipping and copying '%s'" %
                  os.path.basename(cellranger_dir))
            # Tar & gzip data
            targz = os.path.join(working_dir,
                                 "%s.%s.%s.tgz" % (
                                     os.path.basename(
                                         cellranger_dir),
                                     project_name,
                                     project.info.run))
            targz_job = td.run_job(
                f"targz_10x_output.{job_id}.{os.path.basename(cellranger_dir)}",
                Command("tar",
                        "czvf",
                        targz,
                        "-C",
                        os.path.dirname(cellranger_dir),
                        os.path.basename(cellranger_dir)))
            # Copy the targz file
            td.run_job(
                f"copy_10x_tgz.{job_id}.{os.path.basename(cellranger_dir)}",
                copy_command(targz,
                             os.path.join(target_dir,
                                          os.path.basename(targz))),
                wait_for=(targz_job,))

    # Collect 10xGenomics .cloupe files
    if cloupe_files:
        print("Collecting .cloupe files")
        cloupe_dir = os.path.join(working_dir, "10x_cloupe_files")
        os.mkdir(cloupe_dir)
        for f in cloupe_files:
            # Construct a name based on relative path
            relpath = os.path.relpath(os.path.dirname(f), project.dirn).\
                split(os.sep)
            try:
                relpath.remove("outs")
            except ValueError:
                pass
            # Copy into working dir
            ff = os.path.join(cloupe_dir, f"{'-'.join(relpath)}.cloupe")
            shutil.copy(f, ff)
            print(f"...copied {os.path.basename(ff)}")
        # Make ZIP archive for cloupe files
        zip_file = os.path.join(working_dir, "10x_cloupe_files.zip")
        zip_job = td.run_job(f"zip_10x_cloupes_zip.{job_id}",
                             Command("zip",
                                     "-r",
                                     zip_file,
                                     "10x_cloupe_files"),
                             runner=SimpleJobRunner())
        # Copy ZIP archive
        td.run_job(f"copy_10x_cloupes_zip.{job_id}",
                         copy_command(zip_file,
                                      os.path.join(target_dir,
                                                   os.path.basename(zip_file))),
                         wait_for=(zip_job,))

    # Tar and copy Visium images
    if visium_images_dir:
        print(f"Tar gzipping and copying '{visium_images_dir}'")
        # Tar & gzip data
        targz = os.path.join(working_dir,
                             "%s.%s.%s.tgz" % (
                                 os.path.basename(visium_images_dir),
                                 project_name,
                                 project.info.run))
        targz_job = td.run_job(
            f"targz_visium_images.{job_id}.{os.path.basename(visium_images_dir)}",
            Command("tar",
                    "czvf",
                    targz,
                    "-C",
                    os.path.dirname(visium_images_dir),
                    os.path.basename(visium_images_dir)))
        # Copy the targz file
        td.run_job(
            f"copy_visium_images.{job_id}.{os.path.basename(visium_images_dir)}",
            copy_command(targz, os.path.join(target_dir,
                                             os.path.basename(targz))),
            runner=SimpleJobRunner(),
            wait_for=(targz_job,))

    # Wait for jobs to complete
    td.wait()

    # Update the permissions once everything else has copied
    print("Update permissions to read-write")
    td.run_job(f"set_permissions.{job_id}",
               set_permissions_command("u+rwX,g+rwX,o=rX", target_dir),
               runner=SimpleJobRunner())
    td.wait()

    # Finish and check all jobs completed successfully
    status = td.finish()
    if status != 0:
        logger.error(f"{target_dir}: transfer did not complete "
                     "successfully")
        return 1

    # Summarise transfer
    print("================= Transfer complete =================")
    dataset = "%s%s dataset" % ("%s " % project.info.single_cell_platform
                                if project.info.single_cell_platform else '',
                                project.info.library_type)
    endedness = "paired-end" if project.info.paired_end else "single-end"
    summary = [dataset]
    if nfastqs:
        summary.append("-- %d Fastq%s from %d %s sample%s totalling %s" %
                       (nfastqs,
                        's' if nfastqs != 1 else '',
                        nsamples,
                        endedness,
                        's' if nsamples != 1 else '',
                        format_file_size(fsize)))
    if cellranger_dirs:
        summary.append("-- %d 10x Genomics output director%s totalling %s" %
                       (len(cellranger_dirs),
                        'ies' if len(cellranger_dirs) != 1 else 'y',
                        format_file_size(fsize_10x_outputs)))
    if cloupe_files:
        summary.append("-- %d 10x Genomics '.cloupe' file%s" %
                       (len(cloupe_files),
                        's' if len(cloupe_files) != 1 else ''))
    if visium_images_dir:
        summary.append("-- Visium images totalling %s" %
                       format_file_size(fsize_visium_images))
    print("\n".join(summary))
    print("Files now at %s" % target_dir)
    if weburl:
        url = weburl
        if subdir is not None:
            url = os.path.join(url,subdir)
        print("URL: %s" % url)
    print("Done")
    return 0
