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

import os
import re
import argparse
import time
from random import shuffle
from datetime import date
from fnmatch import fnmatch
from bcftbx.JobRunner import fetch_runner
from bcftbx.JobRunner import SimpleJobRunner
from bcftbx.utils import find_program
from bcftbx.utils import format_file_size
from ..analysis import AnalysisDir
from ..analysis import AnalysisProject
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

def main():
    """
    """
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
    args = p.parse_args()

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

    # Check at least one artefact is being transferred
    if not (include_fastqs or
            include_downloader or
            include_qc_report or
            include_10x_outputs or
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

    # Report
    print("Transferring data from '%s' (%s)" % (project.name,
                                                project.dirn))
    print("Fastqs in %s" % project.fastq_dir)

    # Summarise samples and Fastqs
    samples = set()
    nfastqs = 0
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
        if fqs:
            samples.add(sample.name)
            for fq in fqs:
                fsize += os.lstat(fq).st_size
                nfastqs += 1
    nsamples = len(samples)
    if nsamples == 0:
        # No samples found
        logging.error("No samples found")
        return 1
    dataset = "%s%s dataset" % ("%s " % project.info.single_cell_platform
                                if project.info.single_cell_platform else '',
                                project.info.library_type)
    endedness = "paired-end" if project.info.paired_end else "single-end"
    print("%s with %d Fastqs from %d %s sample%s totalling %s" %
          (dataset,
           nfastqs,
           nsamples,
           endedness,
           's' if nsamples != 1 else '',
           format_file_size(fsize)))

    # Check target dir
    if not Location(target_dir).is_remote:
        target_dir = os.path.abspath(target_dir)
    if not exists(target_dir):
        print("'%s': target directory not found" % target_dir)
        return
    else:
        print("Target directory %s" % target_dir)

    # Check hard links are possible
    if hard_links:
        if Location(target_dir).is_remote:
            print("'%s': hard links requested but target directory "
                  "is on a remote filesystem" % target_dir)
            return
        else:
            if os.lstat(project.fastq_dir).st_dev != \
               os.lstat(target_dir).st_dev:
                print("'%s': hard links requested but target directory "
                  "is on a different filesystem" % target_dir)
                return

    # Locate downloader
    if include_downloader:
        print("Locating downloader for inclusion")
        downloader = find_program("download_fastqs.py")
        if downloader is None:
            logging.error("Unable to locate download_fastqs.py")
            return 1
        print("... found %s" % downloader)
    else:
        downloader = None

    # Locate zipped QC report
    if include_qc_report:
        print("Locating zipped QC reports for inclusion")
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
                        print("... found %s" % qc_zip)
                        qc_zips.append(qc_zip)
        if not qc_zips:
            logger.error("No zipped QC reports found")
            return 1
    else:
        qc_zips = None

    # Locate 10xGenomics outputs
    if include_10x_outputs:
        print("Locating outputs from 10xGenomics pipelines for "
              "inclusion")
        cellranger_dirs = list()
        for d in ('cellranger_count',
                  'cellranger_multi',):
            cellranger_dir = os.path.join(project.dirn,d)
            if os.path.isdir(cellranger_dir):
                print("... found %s" % cellranger_dir)
                cellranger_dirs.append(cellranger_dir)
        if not cellranger_dirs:
            logger.error("No outputs from 10xGenomics pipelines found")
            return 1
    else:
        cellranger_dirs = None

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
        print("... found '%s'" % subdir)
        # Update target dir
        target_dir = os.path.join(target_dir,subdir)
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
        # Check it doesn't already exist
        if exists(os.path.join(target_dir,subdir)):
            logger.error("'%s': subdirectory already exists" % subdir)
            return
        print("Using subdirectory '%s'" % subdir)
        # Update target dir
        target_dir = os.path.join(target_dir,subdir)

    # Make target directory
    if not exists(target_dir):
        mkdir(target_dir)

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

    # Construct the README
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

    # Start a scheduler to run jobs
    sched = SimpleScheduler(runner=runner,
                            reporter=TransferDataSchedulerReporter(),
                            poll_interval=settings.general.poll_interval)
    sched.start()

    # List of jobs to check at the end
    check_jobs = {}

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
            print("Running %s" % copy_cmd)
            copy_job = sched.submit(copy_cmd.command_line,
                                    name="copy_fastqs.%s" % job_id,
                                    wd=working_dir)
            check_jobs[copy_job.name] = copy_job
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
            print("Running %s" % zip_cmd)
            zip_job = sched.submit(zip_cmd.command_line,
                                   name="zip_fastqs.%s" % job_id,
                                   wd=working_dir)
            zip_job.wait()
            check_jobs[zip_job.name] = zip_job
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
                    copy_cmd = copy_command(os.path.join(working_dir,f),
                                            os.path.join(target_dir,
                                                         "%s.checksums" %
                                                         final_zip_basename))
                    copy_job = sched.submit(copy_cmd.command_line,
                                            name="copy_checksums.%s" % job_id,
                                            runner=SimpleJobRunner(),
                                            wd=working_dir)
                    check_jobs[copy_job.name] = copy_job
                elif f.endswith(".zip") and \
                   f.startswith("%s." % project_name):
                    # Assume it's ZIP output from packaging process
                    final_zip = "%s%s" % (final_zip_basename,
                                          f[len(project_name):])
                    # Copy to final location
                    copy_cmd = copy_command(
                        os.path.join(working_dir,f),
                        os.path.join(target_dir,final_zip))
                    job_ix += 1
                    copy_job = sched.submit(
                        copy_cmd.command_line,
                        name="copy_zipped_fastqs.%s.%s" % (job_ix,job_id),
                        wd=working_dir,
                        wait_for=(zip_job.job_name,))
                    check_jobs[copy_job.name] = copy_job

    # Copy README
    if readme_file is not None:
        print("Copying README file")
        copy_cmd = copy_command(readme_file,
                                os.path.join(target_dir,"README"))
        copy_job = sched.submit(copy_cmd.command_line,
                                name="copy_readme.%s" % job_id,
                                runner=SimpleJobRunner(),
                                wd=working_dir)
        check_jobs[copy_job.name] = copy_job

    # Copy download_fastqs.py
    if downloader:
        print("Copying downloader")
        copy_cmd = copy_command(downloader,
                                os.path.join(
                                    target_dir,
                                    os.path.basename(downloader)))
        copy_job = sched.submit(copy_cmd.command_line,
                                name="copy_downloader.%s" % job_id,
                                runner=SimpleJobRunner(),
                                wd=working_dir)
        check_jobs[copy_job.name] = copy_job

    # Copy QC reports
    if qc_zips:
        for qc_zip in qc_zips:
            print("Copying '%s'" % os.path.basename(qc_zip))
            copy_cmd = copy_command(
                qc_zip,
                os.path.join(target_dir,os.path.basename(qc_zip)),
                link=hard_links)
            copy_job = sched.submit(copy_cmd.command_line,
                                    name="copy_qc_zip.%s.%s" %
                                    (job_id,
                                     os.path.basename(qc_zip)),
                                    runner=SimpleJobRunner(),
                                    wd=working_dir)
            check_jobs[copy_job.name] = copy_job

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
            targz_cmd = Command("tar",
                                "czvhf",
                                targz,
                                "-C",
                                os.path.dirname(cellranger_dir),
                                os.path.basename(cellranger_dir))
            print("Running %s" % targz_cmd)
            targz_job = sched.submit(targz_cmd.command_line,
                                     name="targz_10x_output.%s.%s" % (
                                         job_id,
                                         os.path.basename(cellranger_dir)),
                                     wd=working_dir)
            check_jobs[targz_job.name] = targz_job
            # Copy the targz file
            copy_cmd = copy_command(targz,
                                    os.path.join(target_dir,
                                                 os.path.basename(targz)))
            print("Running %s" % copy_cmd)
            copy_job = sched.submit(copy_cmd.command_line,
                                    name="copy_10x_tgz.%s.%s" % (
                                        job_id,
                                        os.path.basename(cellranger_dir)),
                                    runner=SimpleJobRunner(),
                                    wd=working_dir,
                                    wait_for=(targz_job.job_name,))
            check_jobs[copy_job.name] = copy_job

    # Wait for scheduler jobs to complete
    sched.wait()

    # Update the permissions once everything else has copied
    print("Update permissions to read-write")
    permissions_cmd = set_permissions_command("u+rwX,g+rwX,o=rX",
                                              target_dir)
    print("Running %s" % permissions_cmd)
    permissions_job = sched.submit(permissions_cmd.command_line,
                                   name="set_permissions.%s" % job_id,
                                   runner=SimpleJobRunner(),
                                   wd=working_dir)
    permissions_job.wait()
    check_jobs[permissions_job.name] = permissions_job

    # Check all jobs completed successfully
    status = 0
    for name in check_jobs:
        if check_jobs[name].exit_code != 0:
            logger.warning("'%s': operation failed" % name)
            status = 1
    if status != 0:
        logger.error("some transfer operations did not complete "
                     "successfully (see warnings above)")
        return 1
    else:
        print("Files now at %s" % target_dir)
        if weburl:
            url = weburl
            if subdir is not None:
                url = os.path.join(url,subdir)
            print("URL: %s" % url)
        print("Done")
