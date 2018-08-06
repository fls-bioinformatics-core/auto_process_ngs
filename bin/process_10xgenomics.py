#!/usr/bin/env python
#
#     process_10xgenomics.py: processing of 10xGenomics Chromium SC data
#     Copyright (C) University of Manchester 2017-2018 Peter Briggs
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
from bcftbx.utils import find_program
from bcftbx.utils import mkdirs
from bcftbx.IlluminaData import IlluminaData
from bcftbx.IlluminaData import IlluminaDataError
from auto_process_ngs.utils import get_numbered_subdir
from auto_process_ngs.analysis import AnalysisProject
from auto_process_ngs.analysis import ProjectMetadataFile
from auto_process_ngs.tenx_genomics_utils import run_cellranger_mkfastq
from auto_process_ngs.tenx_genomics_utils import run_cellranger_count
from auto_process_ngs.tenx_genomics_utils import run_cellranger_count_for_project
import auto_process_ngs.envmod as envmod

# Initialise logging
import logging
logger = logging.getLogger("process_10xgenomics")

# Fetch configuration settings
import auto_process_ngs.settings
__settings = auto_process_ngs.settings.Settings()

######################################################################
# Functions
######################################################################

def cellranger_mkfastq(samplesheet,
                       primary_data_dir,
                       output_dir,
                       lanes=None,
                       cellranger_jobmode='sge',
                       cellranger_maxjobs=None,
                       cellranger_mempercore=None,
                       cellranger_jobinterval=None,
                       log_dir=None,
                       dry_run=False,
                       project_metadata_file='projects.info'):
    """
    Wrapper for running 'cellranger mkfastq'

    Runs the 10xGenomics 'cellranger mkfastq' command to
    generate Fastqs from bcl files for Chromium single-cell
    data.

    Arguments:
      sample_sheet (str): path to input samplesheet with
        10xGenomics barcode indices
      primary_data_dir (str): path to the top-level
        directory holding the sequencing data
      output_dir (str): path to the output directory
      lanes (str): optional, specify the subset of lanes
        to process (default is to process all lanes
        in the run)
      cellranger_jobmode (str): specify the job mode to
        pass to cellranger (default: None)
      cellranger_maxjobs (int): specify the maximum
        number of jobs to pass to cellranger (default:
        None)
      cellranger_mempercore (int): specify the memory
        per core (in Gb) to pass to cellranger (default:
        None)
      cellranger_jobinterval (int): specify the interval
        between launching jobs (in ms) to pass to
        cellranger (default: None)
      log_dir (str): path to a directory to write logs
        (default: current working directory)
      dry_run (bool): if True then only report actions
        that would be performed but don't run anything
      project_metadata_file (str): name of project
        metadata file to create/update with information
        on projects generated by cellranger (default:
        projects.info)

    Returns:
      Integer: exit code from the cellranger command.
    """
    # Make a log directory
    if not dry_run:
        if log_dir is None:
            log_dir = os.getcwd()
        log_dir = get_numbered_subdir("cellranger_mkfastq",
                                      parent_dir=log_dir,
                                      full_path=True)
        mkdirs(log_dir)
    # Run cellranger mkfastq
    retval = run_cellranger_mkfastq(
        samplesheet,
        primary_data_dir,
        output_dir,
        lanes=lanes,
        cellranger_jobmode=cellranger_jobmode,
        cellranger_maxjobs=cellranger_maxjobs,
        cellranger_mempercore=cellranger_mempercore,
        cellranger_jobinterval=cellranger_jobinterval,
        log_dir=log_dir,
        dry_run=dry_run)
    if not dry_run:
        # Update the project metadata file
        update_project_metadata(output_dir,project_metadata_file)
    return retval

def update_project_metadata(unaligned_dir,
                            project_metadata_file):
    """
    Update the entries in a projects.info file

    Adds new entries to a project metadata file (aka
    'projects.info') for projects found in the specified
    bcl2fastq/cellranger mkfastq output directory, or
    updates those which are already present.

    If the metadata file doesn't already exist then a
    new file will be created.

    Arguments:
      unaligned_dir (str): path to bcl2fastq/cellranger
        mkfastq output directory
      project_metadata_file (str): path to a project
        metadata file to update
    """
    analysis_dir = os.path.dirname(os.path.abspath(unaligned_dir))
    unaligned_dir = os.path.basename(unaligned_dir)
    try:
        illumina_data = IlluminaData(analysis_dir,
                                     unaligned_dir=unaligned_dir)
    except IlluminaDataError as ex:
        logger.critical("Failed to load bcl2fastq outputs from "
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
    # Populate/update
    for project in illumina_data.projects:
        project_name = project.name
        sample_names = [s.name for s in project.samples]
        if project_name not in project_metadata:
            project_metadata.add_project(project_name,sample_names)
        else:
            project_metadata.update_project(project_name,
                                            sample_names=sample_names)
    # Save
    project_metadata.save(filen)

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
    count_parser.add_argument("projects",metavar="PROJECT",
                              nargs="*",
                              help="project directory to run "
                              "'cellranger count' on")
    count_parser.add_argument("-u","--unaligned",
                              dest="unaligned_dir",default="bcl2fastq",
                              help="'unaligned' dir with output from "
                              "bcl2fastq (nb ignored if one or more "
                              "projects are supplied)")
    count_parser.add_argument("-t","--transcriptome",
                              dest="transcriptome",default=None,
                              help="directory with reference data for "
                              "transcriptome of interest")
    count_parser.add_argument("-a","--all-outputs",
                              action="store_true",
                              help="collect all outputs from 'cellranger "
                              "count' (default: only collect the "
                              "'web_summary.html' and 'metrics_summary.csv' "
                              "files)")
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
                       dest="job_mode",
                       default=__settings['10xgenomics'].cellranger_jobmode,
                       help="job mode to run cellranger in (default: "
                       "'%s')"
                       % __settings['10xgenomics'].cellranger_jobmode)
        p.add_argument("--mempercore",
                       dest="mem_per_core",
                       default=__settings['10xgenomics'].cellranger_mempercore,
                       help="memory assumed per core (in Gbs; "
                       "default: %d)"
                       % __settings['10xgenomics'].cellranger_mempercore)
        p.add_argument("--maxjobs",type=int,
                       dest="max_jobs",
                       default=__settings.general.max_concurrent_jobs,
                       help="maxiumum number of concurrent jobs to run "
                       "(default: %d)"
                       % __settings.general.max_concurrent_jobs)
        p.add_argument("--jobinterval",type=int,
                       dest="job_interval",
                       default=__settings['10xgenomics'].cellranger_jobinterval,
                       help="how often jobs are submitted (in ms; "
                       "default: %d)"
                       % __settings['10xgenomics'].cellranger_jobinterval)
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

    # Deal with environment modules
    if args.command in ("mkfastq","count"):
        modulefiles = args.modulefiles
        if modulefiles is None:
            try:
                modulefiles = __settings.modulefiles['process_10xgenomics']
            except KeyError:
                # No environment modules specified
                pass
        if modulefiles is not None:
            for modulefile in modulefiles.split(','):
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
            logger.critical("couldn't find '%s'" % prog)
            sys.exit(1)

    # Run the requested command
    if args.command == "mkfastq":
        # Warn that this command is now deprecated
        logger.warning("This command is deprecated, use 'auto_process "
                       "make_fastqs --protocol=10x_chromium_sc' instead")
        cellranger_mkfastq(args.samplesheet,
                           args.run_dir,
                           args.output_dir,
                           lanes=args.lanes,
                           cellranger_jobmode=args.job_mode,
                           cellranger_maxjobs=args.max_jobs,
                           cellranger_mempercore=args.mem_per_core,
                           cellranger_jobinterval=args.job_interval,
                           dry_run=args.dry_run,
                           log_dir='logs',
                           project_metadata_file='projects.info')
    elif args.command == "count":
        # Run cellranger count over the samples
        if args.projects:
            for project in args.projects:
                # Fetch transcriptome for project
                print "Project: %s" % project
                transcriptome = args.transcriptome
                if transcriptome is None:
                    try:
                        organisms = AnalysisProject(
                            os.path.basename(project),
                            project).\
                            info.organism.lower().split(',')
                    except AttributeError:
                        organisms = None
                    if not organisms:
                        raise Exception("%s: can't look up transcriptome "
                                        " (no organism specified); use "
                                        " -t/--transcriptome option" %
                                        project)
                    elif len(organisms) > 1:
                        raise Exception("%s: can't look up transcriptome "
                                        " (multiple organisms); use "
                                        " -t/--transcriptome option" %
                                        project)
                    try:
                        transcriptome = __settings['10xgenomics_transcriptomes']\
                                        [organisms[0]]
                    except KeyError:
                        raise Exception("%s: no transcriptome found for "
                                        "organism '%s'; use -t/--transcriptome"
                                        "option" % (project,organism[0]))
                print "Transcriptome: %s" % transcriptome
                # Run single library analysis
                run_cellranger_count_for_project(
                    project,
                    transcriptome,
                    cellranger_jobmode=args.job_mode,
                    cellranger_maxjobs=args.max_jobs,
                    cellranger_mempercore=args.mem_per_core,
                    cellranger_jobinterval=args.job_interval,
                    max_jobs=args.max_jobs,
                    dry_run=args.dry_run,
                    log_dir='logs')
        else:
            run_cellranger_count(args.unaligned_dir,
                                 args.transcriptome,
                                 cellranger_jobmode=args.job_mode,
                                 cellranger_maxjobs=args.max_jobs,
                                 cellranger_mempercore=args.mem_per_core,
                                 cellranger_jobinterval=args.job_interval,
                                 max_jobs=args.max_jobs,
                                 dry_run=args.dry_run,
                                 log_dir='logs')
    elif args.command == "update_projects":
        # Generate or update the project metadata file
        # Warn that this command is now deprecated
        logger.warning("This command is deprecated")
        update_project_metadata(args.unaligned_dir,
                                args.project_metadata_file)
