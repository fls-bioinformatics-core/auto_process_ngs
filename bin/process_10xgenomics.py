#!/usr/bin/env python
#
#     process_10xgenomics.py: processing of 10xGenomics Chromium SC data
#     Copyright (C) University of Manchester 2017 Peter Briggs
#
"""
process_10genomics.py

Utility to peform initial processing of data from 10xGenomics Chromium
SC 3'v2 platform.
"""

######################################################################
# Imports
######################################################################

import sys
import os
import argparse
import logging
import shutil
from bcftbx.utils import mkdir
from bcftbx.utils import find_program
from bcftbx.utils import list_dirs
from bcftbx.IlluminaData import IlluminaData
from bcftbx.IlluminaData import IlluminaDataError
from bcftbx.JobRunner import SimpleJobRunner
from auto_process_ngs.applications import Command
from auto_process_ngs.simple_scheduler import SimpleScheduler
from auto_process_ngs.simple_scheduler import SchedulerJob
from auto_process_ngs.simple_scheduler import SchedulerReporter
from auto_process_ngs.docwriter import Document
from auto_process_ngs.docwriter import List
from auto_process_ngs.docwriter import Link
import auto_process_ngs.css_rules as css_rules
from auto_process_ngs.utils import ProjectMetadataFile
from auto_process_ngs.utils import ZipArchive
from auto_process_ngs.tenx_genomics_utils import flow_cell_id
from auto_process_ngs.tenx_genomics_utils import make_qc_summary_html
import auto_process_ngs.envmod as envmod

# Initialise logging
import logging
logger = logging.getLogger(__name__)

# Fetch configuration settings
import auto_process_ngs.settings
__settings = auto_process_ngs.settings.Settings()

######################################################################
# Functions
######################################################################

def mkdir(path):
    """
    """
    if os.path.exists(path):
	return
    parent = os.path.dirname(path)
    if not os.path.exists(parent):
        mkdir(parent)
    print "Making dir: %s" % path
    os.mkdir(path)

def cellranger_mkfastq(samplesheet,
                       primary_data_dir,
                       output_dir,
                       lanes=None,
                       cellranger_jobmode='sge',
                       cellranger_maxjobs=None,
                       cellranger_mempercore=None,
                       cellranger_jobinterval=None,
                       log_dir=None,
                       dry_run=False):
    """
    """
    # Construct the command
    cmd = Command("cellranger","mkfastq",
                  "--samplesheet",samplesheet,
                  "--run",primary_data_dir,
                  "--output-dir",output_dir)
    if lanes is not None:
        cmd.add_args("--lanes=%s" % lanes)
    add_cellranger_args(cmd,
                        jobmode=cellranger_jobmode,
                        mempercore=cellranger_mempercore,
                        maxjobs=cellranger_maxjobs,
                        jobinterval=cellranger_jobinterval)
    # Run the command
    print "Running: %s" % cmd
    if not dry_run:
        # Make a log directory
        if log_dir is None:
            log_dir = os.getcwd()
        log_dir = get_log_subdir(log_dir,"cellranger_mkfastq")
        mkdir(log_dir)
        # Submit the job
        cellranger_mkfastq_job = SchedulerJob(
            SimpleJobRunner(),
            cmd.command_line,
            name='cellranger_mkfastq',
            working_dir=os.getcwd(),
            log_dir=log_dir)
        cellranger_mkfastq_job.start()
        try:
            cellranger_mkfastq_job.wait()
        except KeyboardInterrupt,ex:
            logging.warning("Keyboard interrupt, terminating cellranger")
            cellranger_mkfastq_job.terminate()
            raise ex
        exit_code = cellranger_mkfastq_job.exit_code
        print "cellranger mkfastq completed: exit code %s" % exit_code
        if exit_code != 0:
            logging.error("cellranger mkfastq exited with an error")
            return
        # Deal with the QC summary report
        flow_cell_dir = flow_cell_id(primary_data_dir)
        if lanes is not None:
            lanes_suffix = "_%s" % lanes.replace(',','')
        else:
            lanes_suffix = ""
        flow_cell_dir = "%s%s" % (flow_cell_dir,lanes_suffix)
        if not os.path.isdir(flow_cell_dir):
            logging.error("No output directory '%s'" % flow_cell_dir)
            return
        json_file = os.path.join(flow_cell_dir,
                                 "outs",
                                 "qc_summary.json")
        html_file = "cellranger_qc_summary%s.html" % lanes_suffix
        make_qc_summary_html(json_file,html_file)

def cellranger_count(unaligned_dir,
                     transcriptome,
                     cellranger_jobmode='sge',
                     cellranger_maxjobs=None,
                     cellranger_mempercore=None,
                     cellranger_jobinterval=None,
                     max_jobs=4,
                     log_dir=None,
                     dry_run=False,
                     summary_only=True):
    """
    """
    # Input data
    sample_names = {}
    try:
        illumina_data = IlluminaData(os.getcwd(),
                                     unaligned_dir=unaligned_dir)
        for project in illumina_data.projects:
            sample_names[project.name] = []
            for sample in project.samples:
                sample_names[project.name].append(sample.name)
    except IlluminaDataError:
        logging.critical("Couldn't load data from '%s'" %
                         unaligned_dir)
        sys.exit(1)
    print "Samples: %s" % sample_names
    projects = sample_names.keys()

    # Set up a scheduler
    sched_reporter = SchedulerReporter(
        job_start="SCHEDULER: Started  #%(job_number)d: %(job_name)s:\n-- %(command)s",
        job_end=  "SCHEDULER: Finished #%(job_number)d: %(job_name)s"
    )
    sched_reporter = SchedulerReporter()
    sched = SimpleScheduler(max_concurrent=max_jobs,
                            reporter=sched_reporter)
    sched.start()

    # Make a log directory
    if not dry_run:
        if log_dir is None:
            log_dir = os.getcwd()
        log_dir = get_log_subdir(log_dir,"cellranger_count")
        mkdir(log_dir)

    # Submit the cellranger count jobs
    for project in projects:
        print "Project: %s" % project
        for sample in sample_names[project]:
            print "Sample: %s" % sample
            # Check if outputs already exist
            count_dir = os.path.abspath(
                os.path.join(project,
                             "cellranger_count",
                             sample,
                             "outs"))
            if os.path.isdir(count_dir):
                print "-- %s: outputs exist, nothing to do" % sample
                continue
            else:
                print "-- %s: setting up cellranger count" % sample
            # Set up job for this sample
            work_dir = os.path.abspath("tmp.cellranger_count.%s.%s" %
                                       (project,sample))
            mkdir(work_dir)
            cmd = Command("cellranger","count",
                          "--id",sample,
                          "--fastqs",os.path.abspath(unaligned_dir),
                          "--sample",sample,
                          "--transcriptome",transcriptome)
            add_cellranger_args(cmd,
                                jobmode=cellranger_jobmode,
                                mempercore=cellranger_mempercore,
                                maxjobs=cellranger_maxjobs,
                                jobinterval=cellranger_jobinterval)
            print "Running: %s" % cmd
            if not dry_run:
                sched.submit(cmd,
                             name="cellranger_count.%s.%s" %
                             (project,
                              sample),
                             log_dir=log_dir,
                             wd=work_dir)
    sched.wait()

    # If dry run then stop here
    if dry_run:
        return

    # Finished, handle the outputs
    for project in projects:
        print "Project: %s" % project
        for sample in sample_names[project]:
            print "Sample: %s" % sample
            # Destination for count output
            count_dir = os.path.abspath(
                os.path.join(project,
                             "cellranger_count",
                             sample))
            mkdir(count_dir)
            # Copy the cellranger count outputs
            outs_dir = os.path.join("tmp.cellranger_count.%s.%s"
                                    % (project,sample),
                                    sample,
                                    "outs")
            if not summary_only:
                # Collect all outputs
                print "Copying contents of %s to %s" % (outs_dir,count_dir)
                shutil.copytree(outs_dir,count_dir)
            else:
                # Only collect the web summary
                count_dir = os.path.join(count_dir,"outs")
                mkdir(count_dir)
                print "Copying web_summary.html from %s to %s" % (outs_dir,count_dir)
                shutil.copy(os.path.join(outs_dir,"web_summary.html"),count_dir)

    # Create a report and zip archive for each project
    pwd = os.getcwd()
    analysis_dir = os.path.basename(pwd)
    for project in projects:
        # Descend into project dir
        os.chdir(project)
        # Set up zip file
        report_zip = os.path.join("cellranger_count_report.%s.%s.zip" %
                                  (project,analysis_dir))
        zip_file = ZipArchive(report_zip,
                              prefix="cellranger_count_report.%s.%s" %
                              (project,analysis_dir))
        # Construct index page
        print "Making report for project %s" % project
        count_report = Document("%s: cellranger count" % project)
        count_report.add_css_rule(css_rules.QC_REPORT_CSS_RULES)
        summaries = count_report.add_section()
        summaries.add("Reports from cellranger count for each sample:")
        summary_links = List()
        for sample in sample_names[project]:
            # Link to summary for sample
            web_summary = os.path.join("cellranger_count",
                                       sample,
                                       "outs",
                                       "web_summary.html")
            print "Adding web summary (%s) for %s" % (web_summary,
                                                      sample)
            summary_links.add_item(Link("%s" % sample,
                                        web_summary))
            # Add to the zip file
            zip_file.add_file(web_summary)
        summaries.add(summary_links)
        # Write the report and add to the zip file
        html_file = "cellranger_count_report.html"
        count_report.write(html_file)
        zip_file.add_file(html_file)
        # Finish
        zip_file.close()
        os.chdir(pwd)

def update_project_metadata(unaligned_dir,
                            project_metadata_file):
    """
    """
    analysis_dir = os.path.dirname(os.path.abspath(unaligned_dir))
    unaligned_dir = os.path.basename(unaligned_dir)
    try:
        illumina_data = IlluminaData(analysis_dir,
                                     unaligned_dir=unaligned_dir)
    except IlluminaDataError as ex:
        logging.critical("Failed to load bcl2fastq outputs from "
                         "%s/%s: %s" % (analysis_dir,
                                        unaligned_dir,
                                        ex))
        return
    filen = os.path.abspath(project_metadata_file)
    if os.path.exists(filen):
        # Load data from existing file
        print "Loading project metadata from existing file: %s" % filen
        project_metadata = ProjectMetadataFile(filen)
    else:
        # New (empty) metadata file
        print "Creating new project metadata file: %s" % filen
        project_metadata = ProjectMetadataFile()
    # Populate
    for project in illumina_data.projects:
        project_name = project.name
        sample_names = [s.name for s in project.samples]
        project_metadata.add_project(project_name,sample_names)
    # Save
    project_metadata.save(filen)

def add_cellranger_args(cmd,
                        jobmode='sge',
                        maxjobs=None,
                        mempercore=None,
                        jobinterval=None):
    """
    """
    if jobmode is not None:
        cmd.add_args("--jobmode=%s" % jobmode)
    if mempercore is not None:
        cmd.add_args("--mempercore=%s" % mempercore)
    if maxjobs is not None:
        cmd.add_args("--maxjobs=%s" % maxjobs)
    if jobinterval is not None:
        cmd.add_args("--jobinterval=%s" % jobinterval)
    return cmd

def get_log_subdir(log_dir,name):
    """
    """
    # NB based on 'get_log_subdir' from auto_processor.py
    i = 0
    for d in list_dirs(log_dir):
        try:
            i = max(i,int(d.split('_')[0]))
        except ValueError:
            pass
    # Return the full name
    return os.path.join(log_dir,"%03d_%s" % (i+1,str(name)))

######################################################################
# Main
######################################################################

if __name__ == "__main__":
    # Handle the command line
    parser = argparse.ArgumentParser(
        description="Perform operations on FASTQs from 10xGenomics "
        "Chromium SC 3'v2 using 'cellranger'")
    subparsers = parser.add_subparsers(dest='command',
                                       help='Available commands')

    # Build command-specific subparsers
    # 'mkfastq' command
    mkfastq_parser = subparsers.add_parser("mkfastq",
                                           help="run 'cellranger mkfastq'")
    mkfastq_parser.add_argument("-s","--samplesheet",
                                dest="samplesheet",default=None,
                                help="samplesheet file")
    mkfastq_parser.add_argument("-r","--run",
                                dest="run_dir",default=None,
                                help="run directory")
    mkfastq_parser.add_argument("-o","--output_dir",
                                dest="output_dir",default=None,
                                help="output directory")
    mkfastq_parser.add_argument("-l","--lanes",
                                dest="lanes",default=None,
                                help="comma-separated list of lanes "
                                "(optional)")
    # 'count' parser
    count_parser = subparsers.add_parser("count",
                                         help="run 'cellranger count'")
    count_parser.add_argument("-u","--unaligned",
                              dest="unaligned_dir",default="bcl2fastq",
                              help="'unaligned' dir with output from "
                              "bcl2fastq")
    count_parser.add_argument("-t","--transcriptome",
                              dest="transcriptome",default=None,
                              help="directory with reference data for "
                              "transcriptome of interest")
    count_parser.add_argument("-a","--all-outputs",
                              action="store_true",
                              help="collect all outputs from 'cellranger "
                              "count' (default: only collect the "
                              "'web_summary.html' files)")
    # 'update_projects' parser
    update_projects_parser = subparsers.add_parser(
        "update_projects",
        help="update project metadata file")
    update_projects_parser.add_argument("-u","--unaligned",
                                        dest="unaligned_dir",
                                        default="bcl2fastq",
                                        help="'unaligned' dir with "
                                        "output from bcl2fastq")
    update_projects_parser.add_argument("project_metadata_file",
                                        metavar="PROJECT_METADATA_FILE",
                                        nargs="?",
                                        default="projects.info",
                                        help="name of project "
                                        "metadata file (default: "
                                        "'projects.info')")
    # Add generic options
    for p in (mkfastq_parser,count_parser):
        p.add_argument("--jobmode",
                       dest="job_mode",default="sge",
                       help="job mode to run cellranger in (default: "
                       "'sge'")
        p.add_argument("--mempercore",
                       dest="mem_per_core",default=5,
                       help="memory assumed per core (in Gbs; "
                       "default: 5Gb")
        p.add_argument("--maxjobs",type=int,
                       dest="max_jobs",
                       default= __settings.general.max_concurrent_jobs,
                       help="maxiumum number of concurrent jobs to run "
                       "(default: %d)"
                       % __settings.general.max_concurrent_jobs)
        p.add_argument("--jobinterval",type=int,
                       dest="job_interval",default=100,
                       help="how often jobs are submitted (in ms; "
                       "default: 100ms)")
        p.add_argument('--modulefiles',action='store',
                       dest='modulefiles',default=None,
                       help="comma-separated list of environment "
                       "modules to load before executing commands "
                       "(overrides any modules specified in the global "
                       "settings)")
        p.add_argument('--dry-run',action='store_true',
                       dest='dry_run',
                       help="report commands that would be run "
                       "(don't execute them)")
    args = parser.parse_args()
    print "command: %s" % args.command

    # Deal with module files
    if args.modulefiles is not None:
        modulefiles = args.modulefiles.split(',')
        for modulefile in modulefiles:
            envmod.load(modulefile)

    # Check for underlying programs
    if args.command == "mkfastq":
        required = ["cellranger","bcl2fastq"]
    elif args.command == "count":
        required = ["cellranger"]
    else:
        required = []
    for prog in required:
        if find_program(prog) is None:
            logging.critical("couldn't find '%s'" % prog)
            sys.exit(1)

    # Run the requested command
    if args.command == "mkfastq":
        cellranger_mkfastq(args.samplesheet,
                           args.run_dir,
                           args.output_dir,
                           lanes=args.lanes,
                           cellranger_jobmode=args.job_mode,
                           cellranger_maxjobs=args.max_jobs,
                           cellranger_mempercore=args.mem_per_core,
                           cellranger_jobinterval=args.job_interval,
                           dry_run=args.dry_run,
                           log_dir='logs')
    elif args.command == "count":
        # Run cellranger count over the samples
        cellranger_count(args.unaligned_dir,
                         args.transcriptome,
                         cellranger_jobmode=args.job_mode,
                         cellranger_maxjobs=args.max_jobs,
                         cellranger_mempercore=args.mem_per_core,
                         cellranger_jobinterval=args.job_interval,
                         max_jobs=4,
                         dry_run=args.dry_run,
                         log_dir='logs')
    elif args.command == "update_projects":
        # Generate or update the project metadata file
        update_project_metadata(args.unaligned_dir,
                                args.metadata_file)
