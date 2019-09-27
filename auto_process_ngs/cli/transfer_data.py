#!/usr/bin/env python
#
#     cli/transfer_data.py: utility for copying data for sharing
#     Copyright (C) University of Manchester 2019 Peter Briggs
#
#########################################################################
#
# transfer_data.py
#
#########################################################################

import os
import re
import argparse
import tempfile
from random import shuffle
from datetime import date
from bcftbx.JobRunner import fetch_runner
from bcftbx.JobRunner import SimpleJobRunner
from bcftbx.utils import find_program
from ..analysis import AnalysisDir
from ..applications import Command
from ..simple_scheduler import SchedulerJob
from ..fileops import exists
from ..fileops import mkdir
from ..fileops import copy
from ..settings import Settings
from ..settings import get_config_dir
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

# Main function

def main():
    """
    """
    # Command line
    p = argparse.ArgumentParser(
        description="Transfer copies of Fastq data from an analysis "
        "project to an arbitrary destination for sharing with other "
        "people",
        version="%(prog)s "+get_version())
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
    p.add_argument('--readme',action='store',
                   metavar='README_TEMPLATE',dest='readme_template',
                   help="template file to generate README file from; "
                   "can be full path to a template file, or the name "
                   "of a file in the 'templates' directory")
    p.add_argument('--weburl',action='store',
                   help="base URL for webserver (sets the value of "
                   "the WEBURL variable in the template README)")
    p.add_argument('--include_downloader',action='store_true',
                   help="copy the 'download_fastqs.py' utility to the "
                   "final location")
    p.add_argument('--include_qc_report',action='store_true',
                   help="copy the zipped QC reports to the final "
                   "location")
    p.add_argument('--runner',action='store',
                   help="specify the job runner to use for executing "
                   "the checksumming and Fastq copy operations "
                   "(defaults to local job runner)")
    p.add_argument('dest',action='store',metavar="DEST",
                   help="destination to copy Fastqs to; can be an "
                   "arbitrary location of the form [[USER@]HOST:]DIR'")
    p.add_argument('analysis_dir',action='store',
                   metavar="ANALYSIS_DIR",
                   help="analysis directory holding the project to "
                   "copy Fastqs from")
    p.add_argument('project',action='store',nargs='?',
                   metavar="PROJECT[:FASTQ_SET]",
                   help="project directory to copy Fastqs from; "
                   "specifying an optional FASTQ_SET copies them "
                   "from the named Fastq subdirectory")

    # Process command line
    args = p.parse_args()

    target_dir = args.dest
    readme_template = args.readme_template

    # Load configuration
    settings = Settings()

    # Load analysis directory and projects
    analysis_dir = AnalysisDir(args.analysis_dir)
    projects = analysis_dir.projects

    # If no project supplied then list projects
    if args.project is None:
        if len(projects) == 0:
            print("No projects found")
        else:
            for project in projects:
                print("%s" % project.name)
                if len(project.fastq_dirs) > 1:
                    # List the fastq sets if there are more than one
                    # and flag the primary set with an asterisk
                    for d in project.fastq_dirs:
                        is_primary = (d == project.info.primary_fastq_dir)
                        print("- %s%s" % (d,
                                          (" *" if is_primary else "")))
        if analysis_dir.undetermined:
            print("_undetermined")
        return

    # Sort out project and Fastq dir
    try:
        project_name,fastq_dir = args.project.split(':')
    except ValueError:
        project_name = args.project
        fastq_dir = None
    project = None
    for p in projects:
        if project_name == p.name:
            project = p
            break
    if project is None:
        logger.error("'%s': project not found" % project_name)
        return 1
    try:
        project.use_fastq_dir(fastq_dir)
    except Exception as ex:
        logger.error("'%s': failed to load Fastq set '%s': %s" %
                     (project.name,fastq_dir,ex))
        return 1
    print("Transferring data from '%s'" % project.name)

    # Check target dir
    if not exists(target_dir):
        print("'%s': target directory not found" % target_dir)
        return
    else:
        print("Target directory %s" % target_dir)

    # Locate downloader
    if args.include_downloader:
        print("Locating downloader for inclusion")
        downloader = find_program("download_fastqs.py")
        if downloader is None:
            logging.error("Unable to locate download_fastqs.py")
            return 1
        print("... found %s" % downloader)
    else:
        downloader = None

    # Locate zipped QC report
    if args.include_qc_report:
        print("Locating zipped QC reports for inclusion")
        qc_zips = list()
        # Check QC directories and look for zipped reports
        for qc_dir in project.qc_dirs:
            fq_dir = project.qc_info(qc_dir).fastq_dir
            if fq_dir == project.fastq_dir:
                qc_zip = os.path.join(project.dirn,
                                      "%s_report.%s.%s.zip" %
                                      (qc_dir,project.name,
                                       os.path.basename(
                                           analysis_dir.analysis_dir)))
                if os.path.exists(qc_zip):
                    print("... found %s" % qc_zip)
                    qc_zips.append(qc_zip)
        if not qc_zips:
            logger.error("No zipped QC reports found")
            return 1
    else:
        qc_zips = None

    # Determine subdirectory
    if args.subdir == "random_bin":
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
    elif args.subdir == "run_id":
        # Construct subdirectory name based on the
        # run ID
        subdir = "{platform}_{datestamp}.{run_number}-{project}".format(
            platform=analysis_dir.metadata.platform.upper(),
            datestamp=analysis_dir.metadata.instrument_datestamp,
            run_number=analysis_dir.metadata.run_number,
            project=project.name)
        # Check it doesn't already exist
        if exists(os.path.join(target_dir,subdir)):
            logger.error("'%s': subdirectory already exists" % subdir)
            return
        print("Using subdirectory '%s'" % subdir)
        # Update target dir
        target_dir = os.path.join(target_dir,subdir)
    else:
        # No subdir
        subdir = None

    # Make target directory
    if not exists(target_dir):
        mkdir(target_dir)

    # Get runner for copy job
    if args.runner:
        runner = fetch_runner(args.runner)
    else:
        runner = settings.general.default_runner

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
            'WEBURL': args.weburl,
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
        tmpdir = tempfile.mkdtemp()
        readme_file = os.path.join(tmpdir,"README")
        with open(readme_file,'wt') as fp:
            fp.write(readme)
    else:
        # No README
        readme_file = None

    # Build command to run manage_fastqs.py
    copy_cmd = Command("manage_fastqs.py",
                       analysis_dir.analysis_dir,
                       project_name)
    if fastq_dir is not None:
        copy_cmd.add_args(fastq_dir)
    copy_cmd.add_args("copy",target_dir)
    print("Running %s" % copy_cmd)
    copy_job = SchedulerJob(runner,
                            copy_cmd.command_line,
                            name="copy.%s%s" % (project_name,
                                                (fastq_dir
                                                 if fastq_dir is not None
                                                 else '')),
                            working_dir=os.getcwd())
    copy_job.start()
    try:
        copy_job.wait(poll_interval=settings.general.poll_interval)
    except KeyboardInterrupt as ex:
        logger.errore("Keyboard interrupt, terminating file copy")
        copy_job.terminate()
        return 1
    exit_code = copy_job.exit_code
    if exit_code != 0:
        logger.error("File copy exited with an error")
        return exit_code

    # Copy README
    if readme_file is not None:
        print("Copying README file")
        copy(readme_file,os.path.join(target_dir,"README"))

    # Copy download_fastqs.py
    if downloader:
        print("Copying downloader")
        copy(downloader,
             os.path.join(target_dir,os.path.basename(downloader)))

    # Copy QC reports
    if qc_zips:
        for qc_zip in qc_zips:
            print("Copying '%s'" % os.path.basename(qc_zip))
            copy(qc_zip,
                 os.path.join(target_dir,os.path.basename(qc_zip)))

    print("Files now at %s" % target_dir)
    print("Done")
