#!/usr/bin/env python
#
#     cli/auto_process.py: command line interface for auto_process_ngs
#     Copyright (C) University of Manchester 2013-2020 Peter Briggs
#
#########################################################################
#
# auto_process.py
#
#########################################################################

"""
Automated data processing & QC pipeline for Illumina sequence data

Implements a program for automating stages of a standard protocol for
processing and QC'ing Illumina sequencing data.

The stages are:

    setup
    make_fastqs
    setup_analysis_dirs
    run_qc
    publish_qc
    archive
    report

The 'setup' stage creates an analysis directory and acquires the basic
data about the sequencing run from a source directory. Subsequent stages
should be run in sequence to create fastq files, set up analysis
directories for each project, and run QC scripts for each sample in
each project.

The following commands enable the querying and setting of configuration
settings and project metadata:

    config
    params
    metadata

Additional commands are available:

    clone
    samplesheet
    analyse_barcodes
    merge_fastq_dirs
    update_fastq_stats
    import_project
    readme

but these are not part of the standard workflow - they are used for
special cases and testing.
"""

#######################################################################
# Imports
#######################################################################

import sys
import os
import subprocess
import argparse
import time
import logging
import bcftbx.utils as bcf_utils
import bcftbx.platforms
from bcftbx.cmdparse import CommandParser
from bcftbx.cmdparse import add_debug_option
from bcftbx.cmdparse import add_no_save_option
from bcftbx.cmdparse import add_dry_run_option
from bcftbx.cmdparse import add_nprocessors_option
from bcftbx.cmdparse import add_runner_option
from bcftbx.cmdparse import add_arg
from bcftbx.JobRunner import fetch_runner
from .. import get_version
from ..auto_processor import AutoProcess
from ..bcl2fastq.pipeline import PROTOCOLS
from ..bcl2fastq.pipeline import subset
from ..commands.make_fastqs_cmd import BCL2FASTQ_DEFAULTS
from ..commands.report_cmd import ReportingMode
from ..envmod import load as envmod_load
from ..samplesheet_utils import predict_outputs
from ..settings import Settings
from ..settings import locate_settings_file
from ..tenx_genomics_utils import CELLRANGER_ASSAY_CONFIGS
from ..utils import paginate

# Module specific logger
logger = logging.getLogger(__name__)

# Version and configuration
__version__ = get_version()
__settings = Settings()

#######################################################################
# Functions
#######################################################################

# Command line parsers

def add_setup_command(cmdparser):
    """Create a parser for the 'setup' command
    """
    p = cmdparser.add_command('setup',
                              help="Set up a new analysis directory",
                              description="Set up automatic processing of "
                              "Illumina sequencing data from RUN_DIR.")
    p.add_argument('--sample-sheet',action='store',dest='sample_sheet',
                   default=None,
                   help="Copy sample sheet file from name and location "
                   "SAMPLE_SHEET (default is to look for SampleSheet.csv "
                   "inside DIR). SAMPLE_SHEET can be a local or remote "
                   "file, or a URL")
    p.add_argument('-f','--file',action='append',dest='extra_files',
                   metavar="FILE",default=None,
                   help="Additional file(s) to copy into new analysis "
                   "directory (e.g. ICELL8 well list). FILE can be a "
                   "local or remote file, or a URL")
    p.add_argument('--fastq-dir',action='store',dest='unaligned_dir',
                   default=None,
                   help="Import fastq.gz files from FASTQ_DIR (which should "
                   "be a subdirectory of DIR with the same structure as that "
                   "the 'Unaligned' or 'bcl2fastq2' output directory produced "
                   "by CASAVA/bcl2fastq)")
    p.add_argument('--analysis-dir',action='store',dest='analysis_dir',
                   default=None,
                   help="Make new directory called ANALYSIS_DIR (otherwise "
                   "default is '<RUN_DIR>_analysis')")
    p.add_argument('run_dir',metavar="RUN_DIR",
                   help="directory with the output from an Illumina "
                   "sequencer")
    add_debug_option(p)

def add_config_command(cmdparser):
    """
    Create a parser for the 'config' command
    """
    p = cmdparser.add_command('config',
                              help="Query and change global configuration",
                              description="Query and change global "
                              "configuration.")
    p.add_argument('--init',action='store_true',dest='init',default=False,
                   help="Create a new configuration file from the sample.")
    p.add_argument('--set',action='append',dest='key_value',default=None,
                   help="Set the value of a parameter. KEY_VALUE should be "
                   "of the form '<param>=<value>'. Multiple --set options "
                   "can be specified.")
    p.add_argument('--add',action='append',dest='new_section',default=None,
                   help="Add a new section called NEW_SECTION to the config. "
                   "To add a new platform, use 'platform:NAME'. Multiple "
                   "--add options can be specified.")
    add_debug_option(p)
    # Deprecated options
    deprecated = p.add_argument_group('Deprecated/defunct options')
    deprecated.add_argument('--show',action='store_true',dest='show',
                            default=False,
                            help="Show the values of parameters and settings "
                            "(does nothing; use 'config' with no options to "
                            "display settings)")

def add_params_command(cmdparser):
    """
    Create a parser for the 'params' command
    """
    p  = cmdparser.add_command('params',
                               help="Query and change project parameters",
                               description="Query and change processing "
                               "parameters and settings for ANALYSIS_DIR.")
    p.add_argument('--set',action='append',dest='key_value',default=None,
                   help="Set the value of a parameter. KEY_VALUE should be "
                   "of the form '<param>=<value>'. Multiple --set options "
                   "can be specified.")
    add_debug_option(p)
    p.add_argument('analysis_dir',metavar="ANALYSIS_DIR",nargs="?",
                   help="auto_process analysis directory (optional: defaults "
                   "to the current directory)")

def add_metadata_command(cmdparser):
    """
    Create a parser for the 'metadata' command
    """
    p  = cmdparser.add_command('metadata',
                               help="Query and update analysis metadata",
                               description="Query and change metadata "
                               "associated with ANALYSIS_DIR.")
    p.add_argument('--set',action='append',dest='key_value',default=None,
                   help="Set the value of a metadata item. KEY_VALUE should "
                   "be of the form '<param>=<value>'. Multiple --set options "
                   "can be specified.")
    p.add_argument('--update',action='store_true',dest='update',default=False,
                   help="Automatically update metadata items where possible "
                   "(e.g. for older analyses which have old or missing "
                   "metadata files)")
    add_debug_option(p)
    p.add_argument('analysis_dir',metavar="ANALYSIS_DIR",nargs="?",
                   help="auto_process analysis directory (optional: defaults "
                   "to the current directory)")

def add_samplesheet_command(cmdparser):
    """
    Create a parser for the 'samplesheet' command
    """
    p = cmdparser.add_command('samplesheet',
                              help="Sample sheet manipulation",
                              description="Query and manipulate sample sheets")
    p.add_argument('-e','--edit',action='store_true',dest='edit',
                   default=False,
                   help="bring up sample sheet file in an editor to make "
                   "changes manually")
    add_debug_option(p)
    p.add_argument('analysis_dir',metavar="ANALYSIS_DIR",nargs="?",
                   help="auto_process analysis directory (optional: defaults "
                   "to the current directory)")

def add_make_fastqs_command(cmdparser):
    """
    Create a parser for the 'make_fastqs' command
    """
    p = cmdparser.add_command('make_fastqs',
                              help="Run Fastq generation",
                              description="Generate fastq files from raw "
                              "bcl files produced by Illumina sequencer.")
    # General options
    add_modulefiles_option(p)
    add_no_save_option(p)
    add_debug_option(p)
    # Primary data management
    primary_data = p.add_argument_group('Primary data management')
    primary_data.add_argument('--force-copy',action='store_true',
                              dest='force_copy',default=False,
                              help="force primary data to be copied (by "
                              "default only data on a remote system will be "
                              "copied; data on a local system will be "
                              "symlinked)")
    # General Fastq generation options
    fastq_generation = p.add_argument_group('General Fastq generation')
    fastq_generation.add_argument('--protocol',choices=PROTOCOLS,
                                  dest='protocol',default='standard',
                                  help="specify Fastq generation protocol "
                                  "depending on the data being processed "
                                  "(default: 'standard')")
    fastq_generation.add_argument('--sample-sheet',action="store",
                                  dest="sample_sheet",default=None,
                                  help="use an alternative sample sheet to "
                                  "the default 'custom_SampleSheet.csv' "
                                  "created on setup.")
    fastq_generation.add_argument('--lanes',action='append',
                                  dest='lanes',
                                  metavar="LANES[:OPTIONS]",
                                  help="define a set of lanes to group for "
                                  "processing. LANES can be a single lane "
                                  "(e.g. '1'), a list ('1,2,3,7'), a range "
                                  "('1-3'), or a combination ('1-3,7'). "
                                  "Specified lanes are processed together "
                                  "in a group, using OPTIONS (if supplied). "
                                  "OPTIONS takes the form "
                                  "'[PROTOCOL:][KEY=VALUE:[KEY=VALUE]...] "
                                  "(for example "
                                  "--lanes=1-4:standard:trim_adapters=no)")
    fastq_generation.add_argument('--output-dir',action='store',
                                  dest='out_dir',default=None,
                                  help="set the directory for the output "
                                  "Fastqs (default: 'bcl2fastq')")
    fastq_generation.add_argument('--platform',action="store",
                                  dest="platform",default=None,
                                  help="explicitly specify the sequencing "
                                  "platform. Only use this if the platform "
                                  "cannot be identified from the instrument "
                                  "name")
    fastq_generation.add_argument('--use-bases-mask',action="store",
                                  dest="bases_mask",default=None,
                                  help="explicitly set the bases-mask string "
                                  "to indicate how each cycle should be used "
                                  "in the bcl2fastq conversion (overrides "
                                  "default). Set to 'auto' to determine "
                                  "automatically")
    # Options to control bcl2fastq
    bcl_to_fastq = p.add_argument_group('Bcl2fastq options')
    # Use lane splitting
    # Determine defaults to report to user
    no_lane_splitting_platforms = []
    use_lane_splitting_platforms = []
    for platform in __settings.platform:
        if __settings.platform[platform].no_lane_splitting is not None:
            if __settings.platform[platform].no_lane_splitting:
                no_lane_splitting_platforms.append(platform)
            else:
                use_lane_splitting_platforms.append(platform)
    if __settings.bcl2fastq.no_lane_splitting:
        if use_lane_splitting_platforms:
            default_no_lane_splitting = \
                            "Used by default for all platforms except %s" % \
                            ', '.join(use_lane_splitting_platforms)
            default_use_lane_splitting = "Used by default for %s" % \
                            ', '.join(use_lane_splitting_platforms)
        else:
            default_no_lane_splitting = "Default for all platforms"
            default_use_lane_splitting = ""
    else:
        if no_lane_splitting_platforms:
            default_use_lane_splitting = \
                            "Used by default for all platforms except %s" % \
                            ', '.join(no_lane_splitting_platforms)
            default_no_lane_splitting = "Used by default for %s" % \
                            ', '.join(no_lane_splitting_platforms)
        else:
            default_no_lane_splitting = ""
            default_use_lane_splitting = "Default for all platforms"
    if default_use_lane_splitting:
        default_use_lane_splitting = ". "+default_use_lane_splitting
    if default_no_lane_splitting:
        default_no_lane_splitting = ". "+default_no_lane_splitting
    bcl_to_fastq.add_argument('--no-lane-splitting',action='store_true',
                              dest='no_lane_splitting',default=False,
                              help="don't split the output FASTQ files by "
                              "lane (bcl2fastq v2 only; turn off using "
                              "--use-lane-splitting)%s" %
                              default_no_lane_splitting)
    bcl_to_fastq.add_argument('--use-lane-splitting',action='store_true',
                              dest='use_lane_splitting',default=False,
                              help="split the output FASTQ files by lane "
                              "(bcl2fastq v2 only; turn off using "
                              "--no-lane-splitting)%s" %
                              default_use_lane_splitting)
    # Creation of empty fastqs
    # Determine defaults to report to user
    create_empty_fastqs_platforms = []
    no_create_empty_fastqs_platforms = []
    for platform in __settings.platform:
        if __settings.platform[platform].create_empty_fastqs is not None:
            if __settings.platform[platform].create_empty_fastqs:
                create_empty_fastqs_platforms.append(platform)
            else:
                no_create_empty_fastqs_platforms.append(platform)
    if __settings.bcl2fastq.create_empty_fastqs:
        if no_create_empty_fastqs_platforms:
            default_create_empty_fastqs = \
                    "Used by default for all platforms except %s" % \
                    ', '.join(no_create_empty_fastqs_platforms)
            default_no_create_empty_fastqs = \
                    "Used by default for %s" % \
                    ', '.join(no_create_empty_fastqs_platforms)
        else:
            default_create_empty_fastqs = "Default for all platforms"
            default_no_create_empty_fastqs = ""
    else:
        if create_empty_fastqs_platforms:
            default_no_create_empty_fastqs = \
                    "Used by default for all platforms except %s" % \
                    ', '.join(create_empty_fastqs_platforms)
            default_create_empty_fastqs = \
                    "Used by default for %s" % \
                    ', '.join(create_empty_fastqs_platforms)
        else:
            default_no_create_empty_fastqs = "Default for all platforms"
            default_create_empty_fastqs = ""
    if default_create_empty_fastqs:
        default_create_empty_fastqs = ". "+default_create_empty_fastqs
    if default_no_create_empty_fastqs:
        default_no_create_empty_fastqs = ". "+default_no_create_empty_fastqs
    bcl_to_fastq.add_argument('--create-empty-fastqs',action='store_true',
                              dest='create_empty_fastqs',
                              default=False,
                              help="create empty files as placeholders for "
                              "missing FASTQs from bcl2fastq demultiplexing "
                              "step (NB bcl2fastq must have finished without "
                              "an error for this option to be applied; turn "
                              "off using --no-create-empty-fastqs)%s"
                              % default_create_empty_fastqs)
    bcl_to_fastq.add_argument('--no-create-empty-fastqs',action='store_true',
                              dest='no_create_empty_fastqs',
                              default=False,
                              help="don't create empty files as placeholders "
                              "for missing FASTQs from bcl2fastq "
                              "demultiplexing step (turn off using "
                              "--create-empty-fastqs)%s"
                              % default_no_create_empty_fastqs)
    # Fastqs for index reads
    bcl_to_fastq.add_argument('--create-fastq-for-index-reads',
                              action='store_true',
                              dest='create_fastq_for_index_read',
                              default=False,
                              help="also create FASTQs for index reads")
    # Require specific version
    bcl_to_fastq.add_argument('--require-bcl2fastq-version',action='store',
                              dest='bcl2fastq_version',default=None,
                              help="explicitly specify version of bcl2fastq "
                              "software to use (e.g. '=1.8.4' or '>=2.0').")
    # Number of processors
    default_nprocessors = []
    for platform in __settings.platform:
        if __settings.platform[platform].nprocessors is not None:
            default_nprocessors.append(
                "%s: %s" % 
                (platform,
                 __settings.platform[platform].nprocessors))
    if default_nprocessors:
        default_nprocessors.append("other platforms: %s" %
                                   __settings.bcl2fastq.nprocessors)
    else:
        default_nprocessors.append("%s" % __settings.bcl2fastq.nprocessors)
    default_nprocessors = ', '.join(default_nprocessors)
    add_nprocessors_option(bcl_to_fastq,default_nprocessors,
                           default_display=default_nprocessors)
    add_runner_option(bcl_to_fastq)
    # Adapter trimming/masking options
    adapters = p.add_argument_group('Adapter trimming and masking')
    adapters.add_argument('--adapter',action="store",
                          dest="adapter_sequence",default=None,
                          help="sequence of adapter to be trimmed. "
                          "Specify multiple adapters by separating "
                          "them with plus sign (+). Only used for read "
                          "1 if --adapter-read2 is also specified "
                          "(default: use adapter sequence from sample "
                          "sheet)")
    adapters.add_argument('--adapter-read2',action="store",
                          dest="adapter_sequence_read2",default=None,
                          help="sequence of adapter to be trimmed in "
                          "read 2. Specify multiple adapters by separating "
                          "them with plus sign (+) (default: use adapter "
                          "sequence from sample sheet)")
    adapters.add_argument('--minimum-trimmed-read-length',action="store",
                          dest="minimum_trimmed_read_length",default=None,
                          help="Minimum read length after adapter "
                          "trimming. bcl2fastq trims the adapter from "
                          "the read down to this value; if there is more "
                          "adapter match below this length then those "
                          "bases are masked not trimmed (i.e. replaced "
                          "by N rather than removed) (default: %d)" %
                          BCL2FASTQ_DEFAULTS['minimum_trimmed_read_length'])
    adapters.add_argument('--mask-short-adapter-reads',action="store",
                          dest="mask_short_adapter_reads",default=None,
                          help="minimum length of unmasked bases that "
                          "a read can be after adapter trimming; reads "
                          "with fewer ACGT bases will be completely "
                          "masked with Ns (default: %d)" %
                          BCL2FASTQ_DEFAULTS['mask_short_adapter_reads'])
    adapters.add_argument('--no-adapter-trimming',action="store_true",
                          dest="no_adapter_trimming",default=False,
                          help="turn off adapter trimming even if "
                          "adapter sequences are supplied")
    # ICELL8 options
    icell8 = p.add_argument_group('ICELL8 options (ICELL8 data only)')
    icell8.add_argument("--well-list",
                        dest="icell8_well_list",
                        default=None,
                        help="specify ICELL8 well list file")
    icell8.add_argument('--swap-i1-and-i2',
                        action='store_true',
                        dest="icell8_swap_i1_and_i2",
                        help="swap supplied I1 and I2 Fastqs when matching "
                        "ATAC barcodes against well list")
    icell8.add_argument('--reverse-complement',
                        choices=['i1','i2','both'],
                        dest="icell8_reverse_complement",
                        default=None,
                        help="can be 'i1','i2', or 'both'; reverse complement "
                        "the specified indices from the well list when "
                        "matching ATAC barcodes against well list")
    # Cellranger (10xgenomics Chromium SC 3') options
    default_cellranger_localcores = __settings['10xgenomics'].\
                                    cellranger_localcores
    if not default_cellranger_localcores:
        default_cellranger_localcores = 1
    cellranger = p.add_argument_group('Cellranger options (10xGenomics '
                                      'Chromium SC 3\' and ATAC data only)')
    cellranger.add_argument("--jobmode",
                            dest="job_mode",
                            default=__settings['10xgenomics'].\
                            cellranger_jobmode,
                            help="job mode to run cellranger in (default: "
                            "'%s')"
                            % __settings['10xgenomics'].cellranger_jobmode)
    cellranger.add_argument("--mempercore",
                            dest="mem_per_core",
                            default=__settings['10xgenomics'].\
                            cellranger_mempercore,
                            help="memory assumed per core (in Gbs; default: "
                            "%s); NB only used if jobmode is not 'local'" %
                            __settings['10xgenomics'].cellranger_mempercore)
    cellranger.add_argument("--localcores",
                            dest="local_cores",
                            default=default_cellranger_localcores,
                            help="maximum cores cellranger can request at one"
                            "time for jobmode 'local' (ignored for other "
                            "jobmodes) (default: %s)" %
                            ("taken from job runner"
                             if not default_cellranger_localcores
                             else default_cellranger_localcores))
    cellranger.add_argument("--localmem",
                            dest="local_mem",
                            default=__settings['10xgenomics'].\
                            cellranger_localmem,
                            help="maximum total memory cellranger can "
                            "request at one time for jobmode 'local' "
                            "(ignored for other jobmodes) (in Gbs; default: "
                            "%s)" %
                            __settings['10xgenomics'].cellranger_localmem)
    cellranger.add_argument("--maxjobs",type=int,
                            dest="max_jobs",
                            default=__settings.general.max_concurrent_jobs,
                            help="maxiumum number of concurrent jobs to run "
                            "(default: %d)"
                            % __settings.general.max_concurrent_jobs)
    cellranger.add_argument("--jobinterval",type=int,
                            dest="job_interval",
                            default=__settings['10xgenomics'].\
                            cellranger_jobinterval,
                            help="how often jobs are submitted (in ms; "
                            "default: %d)" %
                            __settings['10xgenomics'].cellranger_jobinterval)
    cellranger.add_argument("--ignore-dual-index",
                            action="store_true",
                            dest="ignore_dual_index",
                            help="on a dual-indexed flowcell where the "
                            "second index was not used for the 10x sample, "
                            "ignore it")
    # Statistics
    statistics = p.add_argument_group('Statistics generation')
    statistics.add_argument('--stats-file',action='store',
                            dest='stats_file',default=None,
                            help="specify output file for fastq statistics")
    statistics.add_argument('--per-lane-stats-file',action='store',
                            dest='per_lane_stats_file',default=None,
                            help="specify output file for per-lane statistics")
    statistics.add_argument('--no-stats',action='store_true',
                            dest='no_stats',default=False,
                            help="don't generate statistics file; use "
                            "'update_fastq_stats' command to (re)generate "
                            "statistics")
    # Barcode analysis
    barcodes = p.add_argument_group('Barcode analysis')
    barcodes.add_argument('--barcode-analysis-dir',action="store",
                          dest="barcode_analysis_dir",default=None,
                          help="specify subdirectory where barcode analysis "
                          "will be performed and outputs will be written")
    barcodes.add_argument('--no-barcode-analysis',action='store_true',
                          dest='no_barcode_analysis',default=False,
                          help="don't perform barcode analysis; use "
                          "'analyse_barcodes' command to run barcode "
                          "analysis separately")
    p.add_argument('analysis_dir',metavar="ANALYSIS_DIR",nargs="?",
                   help="auto_process analysis directory (optional: defaults "
                   "to the current directory)")
    # Advanced options
    advanced = p.add_argument_group('Advanced/debugging options')
    advanced.add_argument('--verbose',action="store_true",
                          dest="verbose",default=False,
                          help="run pipeline in 'verbose' mode")
    advanced.add_argument('--work-dir',action="store",
                          dest="working_dir",default=None,
                          help="specify the working directory for the "
                          "pipeline operations")

def add_setup_analysis_dirs_command(cmdparser):
    """Create a parser for the 'setup_analysis_dirs' command
    """
    p = cmdparser.add_command('setup_analysis_dirs',
                              help="Create project subdirectories",
                              description="Create analysis subdirectories "
                              "for projects defined in projects.info file "
                              "in ANALYSIS_DIR.")
    p.add_argument('--ignore-missing-metadata',action='store_true',
                   dest='ignore_missing_metadata',default=False,
                   help="force creation of project directories even if "
                   "metadata is not set (default is to fail if metadata "
                   "is missing)")
    p.add_argument('--unaligned-dir',action='store',
                   dest='unaligned_dir',default=None,
                   help="explicitly specify the subdirectory with "
                   "output Fastqs")
    p.add_argument('--undetermined',action='store',
                   dest='undetermined',default=None,
                   help="explicitly specify name for project directory "
                   "with 'undetermined' fastqs")
    p.add_argument('--short-fastq-names',action='store_true',
                   dest='short_fastq_names',default=False,
                   help="shorten fastq file names when copying or linking "
                   "from project directory (default is to keep long names "
                   "from bcl2fastq)")
    p.add_argument('--link-to-fastqs',action='store_true',
                   dest='link_to_fastqs',default=False,
                   help="create symbolic links to original fastqs from "
                   "project directory (default is to make hard links)")
    add_debug_option(p)
    p.add_argument('analysis_dir',metavar="ANALYSIS_DIR",nargs="?",
                   help="auto_process analysis directory (optional: defaults "
                   "to the current directory)")

def add_run_qc_command(cmdparser):
    """
    Create a parser for the 'run_qc' command
    """
    p = cmdparser.add_command('run_qc',help="Run QC procedures",
                              description="Run QC procedures for sequencing "
                              "projects in ANALYSIS_DIR.")
    # Defaults
    default_nthreads = __settings.qc.nprocessors
    fastq_screen_subset = __settings.qc.fastq_screen_subset
    max_concurrent_jobs = __settings.general.max_concurrent_jobs
    # Build parser
    p.add_argument('--projects',action='store',
                   dest='project_pattern',default=None,
                   help="simple wildcard-based pattern specifying a "
                   "subset of projects and samples to run the QC on. "
                   "PROJECT_PATTERN should be of the form 'pname[/sname]', "
                   "where 'pname' specifies a project (or set of "
                   "projects) and 'sname' optionally specifies a sample "
                   "(or set of samples).")
    p.add_argument('--fastq_screen_subset',action='store',dest='subset',
                   type=int,default=fastq_screen_subset,
                   help="specify size of subset of total reads to use for "
                   "fastq_screen (i.e. --subset option); (default %d, set to "
                   "0 to use all reads)" % fastq_screen_subset)
    p.add_argument('-t','--threads',action='store',dest="nthreads",
                   type=int,default=default_nthreads,
                   help="number of threads to use for QC script "
                   "(default: %s)" % ('taken from job runner'
                                      if not default_nthreads
                                      else default_nthreads,))
    p.add_argument('--ungzip-fastqs',action='store_true',dest='ungzip_fastqs',
                   help="create decompressed copies of fastq.gz files")
    p.add_argument('--max-jobs',action='store',
                   dest='max_jobs',default=max_concurrent_jobs,type=int,
                   help="explicitly specify maximum number of concurrent QC "
                   "jobs to run (default %s, change in settings file)"
                   % max_concurrent_jobs)
    p.add_argument('--qc_dir',action='store',dest='qc_dir',default='qc',
                   help="explicitly specify QC output directory (nb if "
                   "supplied then the same QC_DIR will be used for each "
                   "project. Non-absolute paths are assumed to be relative to "
                   "the project directory). Default: 'qc'")
    p.add_argument('--fastq_dir',action='store',dest='fastq_dir',default=None,
                   help="explicitly specify subdirectory of DIR with "
                   "Fastq files to run the QC on.")
    p.add_argument("--10x_chemistry",
                   choices=sorted(CELLRANGER_ASSAY_CONFIGS.keys()),
                   dest="cellranger_chemistry",default="auto",
                   help="assay configuration for 10xGenomics scRNA-seq; if "
                   "set to 'auto' (the default) then cellranger will attempt "
                   "to determine this automatically")
    p.add_argument('--10x_transcriptome',action='append',
                   metavar='ORGANISM=REFERENCE',
                   dest='cellranger_transcriptomes',
                   help="specify cellranger transcriptome reference datasets "
                   "to associate with organisms (overrides references defined "
                   "in config file)")
    p.add_argument('--10x_premrna_reference',action='append',
                   metavar='ORGANISM=REFERENCE',
                   dest='cellranger_premrna_references',
                   help="specify cellranger pre-mRNA reference datasets "
                   "to associate with organisms (overrides references defined "
                   "in config file)")
    p.add_argument('--report',action='store',dest='html_file',default=None,
                   help="file name for output HTML QC report (default: "
                   "<QC_DIR>_report.html)")
    add_runner_option(p)
    add_modulefiles_option(p)
    # Advanced options
    advanced = p.add_argument_group('Advanced/debugging options')
    advanced.add_argument('--verbose',action="store_true",
                          dest="verbose",default=False,
                          help="run pipeline in 'verbose' mode")
    advanced.add_argument('--work-dir',action="store",
                          dest="working_dir",default=None,
                          help="specify the working directory for the "
                          "pipeline operations")
    add_debug_option(advanced)
    # Deprecated options
    deprecated = p.add_argument_group('Deprecated/defunct options')
    deprecated.add_argument('--no-ungzip-fastqs',action='store_true',
                            dest='no_ungzip_fastqs',
                            help="don't create uncompressed copies of "
                            "fastq.gz files (does nothing; this is now the "
                            "default, use --ungzip-fastqs to turn on "
                            "decompression)")
    p.add_argument('analysis_dir',metavar="ANALYSIS_DIR",nargs="?",
                   help="auto_process analysis directory (optional: defaults "
                   "to the current directory)")

def add_publish_qc_command(cmdparser):
    """
    Create a parser for the 'publish_qc' command
    """
    # Defaults
    default_use_hierarchy = ("yes" if __settings.qc_web_server.use_hierarchy
                             else "no")
    default_exclude_zips = ("yes" if __settings.qc_web_server.exclude_zip_files
                             else "no")
    p = cmdparser.add_command('publish_qc',
                              help="Copy QC reports to publication area",
                              description="Copy QC reports from ANALYSIS_DIR "
                              "to local or remote directory (e.g. web "
                              "server). By default existing QC reports will "
                              "be copied without further checking; if no "
                              "report is found then QC results will be "
                              "verified and a report generated first.")
    p.add_argument('--projects',action='store',
                   dest='project_pattern',default=None,
                   help="simple wildcard-based pattern specifying a subset "
                   "of projects and samples to publish the QC for. "
                   "PROJECT_PATTERN can specify a single project, or a set "
                   "of projects.")
    p.add_argument('--qc_dir',action='store',
                   dest='qc_dir',default=None,
                   help="specify target directory to copy QC reports to. "
                   "QC_DIR can be a local directory, or a remote location "
                   "in the form '[[user@]host:]directory'. Overrides the "
                   "default settings.")
    p.add_argument('--use-hierarchy',choices=["yes","no"],
                   dest='use_hierarchy',default=default_use_hierarchy,
                   help="use YEAR/PLATFORM hierarchy under QC_DIR; can be "
                   "'yes' or 'no' (default: %s)" % default_use_hierarchy)
    p.add_argument('--exclude-zip-files',choices=["yes","no"],
                   dest='exclude_zip_files',default=default_exclude_zips,
                   help="exclude ZIP archives from publication; can be 'yes' "
                   "or 'no' (default: %s)" % default_exclude_zips)
    p.add_argument('--ignore-missing-qc',action='store_true',
                   dest='ignore_missing_qc',default=False,
                   help="skip projects where QC results are missing or "
                   "can't be verified, or where reports can't be generated.")
    p.add_argument('--regenerate-reports',action='store_true',
                   dest='regenerate_reports',default=False,
                   help="attempt to regenerate existing QC reports")
    p.add_argument('--force',action='store_true',
                   dest='force',default=False,
                   help="force generation of QC reports for all projects "
                   "even if verification has failed")
    p.add_argument('--legacy',action='store_true',
                   dest='legacy_mode',default=False,
                   help="legacy mode: include links to MultiQC, cellranger "
                   "count and ICELL8 reports in the top-level index page")
    add_runner_option(p)
    add_debug_option(p)
    p.add_argument('analysis_dir',metavar="ANALYSIS_DIR",nargs="?",
                   help="auto_process analysis directory (optional: defaults "
                   "to the current directory)")

def add_archive_command(cmdparser):
    """
    Create a parser for the 'archive' command
    """
    # Defaults
    default_group = __settings.archive.group
    default_chmod = __settings.archive.chmod
    # Build parser
    p = cmdparser.add_command('archive',
                              help="Copy analyses to 'archive' area",
                              description="Copy sequencing analysis data "
                              "directory ANALYSIS_DIR to 'archive' "
                              "destination.")
    p.add_argument('--archive_dir',action='store',
                   dest='archive_dir',default=None,
                   help="specify top-level archive directory to copy data "
                   "under. ARCHIVE_DIR can be a local directory, or a "
                   "remote location in the form '[[user@]host:]directory'. "
                   "Overrides the default settings.")
    p.add_argument('--platform',action='store',
                   dest='platform',default=None,
                   help="specify the platform e.g. 'hiseq', 'miseq' etc "
                   "(overrides automatically determined platform, if any). "
                   "Use 'other' for cases where the platform is unknown.")
    p.add_argument('--year',action='store',
                   dest='year',default=None,
                   help="specify the year e.g. '2014' (default is the "
                   "current year)")
    p.add_argument('--group',action='store',dest='group',
                   default=default_group,
                   help="specify the name of group for the archived files "
                   "(default: %s)" % default_group)
    p.add_argument('--chmod',action='store',dest='permissions',
                   default=default_chmod,
                   help="specify permissions for the archived files. "
                   "PERMISSIONS should be a string recognised by the "
                   "'chmod' command (e.g. 'o-rwX') (default: %s)" %
                   default_chmod)
    p.add_argument('--final',action='store_true',dest='final',
                   default=False,
                   help="copy data to final archive location (default is "
                   "to copy to staging area)")
    p.add_argument('--force',action='store_true',dest='force',default=False,
                   help="attempt to complete archiving operations ignoring "
                   "any errors (e.g. key metadata items not set, unable to "
                   "set group etc)")
    add_dry_run_option(p)
    add_debug_option(p)
    p.add_argument('analysis_dir',metavar="ANALYSIS_DIR",nargs="?",
                   help="auto_process analysis directory (optional: defaults "
                   "to the current directory)")

def add_report_command(cmdparser):
    """Create a parser for the 'report' command
    """
    p  = cmdparser.add_command('report',
                               help="Generate reporting information",
                               description="Report information on analysis "
                               "in ANALYSIS_DIR.")
    p.add_argument('--logging',action='store_true',dest='logging',
                   default=False,
                   help="print short report suitable for logging file")
    p.add_argument('--summary',action='store_true',dest='summary',
                   default=False,
                   help="print full report suitable for bioinformaticians")
    p.add_argument('--projects',action='store_true',dest='projects',
                   default=False,
                   help="print tab-delimited line (one per project) "
                   "suitable for injection into a spreadsheet")
    p.add_argument('--fields',action='store',dest='fields',default=None,
                   help="fields to report")
    p.add_argument('--template',action='store',dest='template',
                   default=None,
                   help="name of template with fields to report (templates "
                   "should be defined in the config file)")
    p.add_argument('--file',action='store',dest='out_file',default=None,
                   help="write report to OUT_FILE; destination can be a "
                   "local file, or a remote file specified as "
                   "[[USER@]HOST:]PATH (default is to write to stdout)")
    add_debug_option(p)
    p.add_argument('analysis_dir',metavar="ANALYSIS_DIR",nargs="?",
                   help="auto_process analysis directory (optional: defaults "
                   "to the current directory)")

def add_readme_command(cmdparser):
    """
    Create a parser for the 'readme' command
    """
    p = cmdparser.add_command('readme',
                              help="Add or amend top-level README file",
                              description="Add or amend a README file in the "
                              "analysis directory DIR.")
    p.add_argument('--init',action='store_true',dest='init',
                   default=False,
                   help="create a new README file")
    p.add_argument('-V','--view',action='store_true',dest='view',
                   default=False,
                   help="display the contents of the README file")
    p.add_argument('-e','--edit',action='store_true',dest='edit',
                   default=False,
                   help="bring up README file in an editor to make changes")
    p.add_argument('-m','--message',action='store',dest='message',
                   default=None,
                   help="append MESSAGE text to the README file")
    add_debug_option(p)
    p.add_argument('analysis_dir',metavar="ANALYSIS_DIR",nargs="?",
                   help="auto_process analysis directory (optional: defaults "
                   "to the current directory)")

def add_clone_command(cmdparser):
    """
    Create a parser for the 'clone' command
    """
    p = cmdparser.add_command('clone',
                              help="Make a copy of an analysis directory",
                              description="Make a copy of an existing "
                              "directory DIR in a new directory CLONE_DIR.")
    p.add_argument('--copy-fastqs',action='store_true',dest='copy_fastqs',
                   default=False,
                   help="Copy fastq.gz files from DIR into CLONE_DIR "
                   "(default is to make a link to the bcl-to-fastq "
                   "directory)")
    p.add_argument('--exclude-projects',action='store_true',
                   dest='exclude_projects',default=False,
                   help="Exclude (i.e. don't copy) project directories "
                   "from DIR")
    add_debug_option(p)
    p.add_argument('analysis_dir',metavar="ANALYSIS_DIR",nargs='?',
                   help="existing auto_process analysis directory to "
                   "clone (optional: defaults to the current directory)")
    p.add_argument('clone_dir',metavar="CLONE_DIR",
                   help="path to cloned directory")

def add_analyse_barcodes_command(cmdparser):
    """
    Create a parser for the 'analyse_barcodes' command
    """
    p = cmdparser.add_command('analyse_barcodes',
                              help="Analyse index (barcode) sequences",
                              description="Analyse barcode sequences for "
                              "Fastq files in specified lanes in "
                              "ANALYSIS_DIR, and report the most "
                              "common barcodes found across all reads from "
                              "each lane.")
    p.add_argument('--unaligned-dir',action='store',
                   dest='unaligned_dir',default='bcl2fastq',
                   help="explicitly set the (sub)directory with bcl-to-fastq "
                   "outputs")
    p.add_argument('--lanes',action='store',
                   dest='lanes',default=None,
                   help="specify which lanes to analyse barcodes for "
                   "(default is to do analysis for all lanes).")
    p.add_argument('--mismatches',action='store',dest='mismatches',
                   default=None,type=int,
                   help="maximum number of mismatches to use when grouping "
                   "similar barcodes (default is to determine automatically "
                   "from the bases mask)")
    p.add_argument('--cutoff',action='store',dest='cutoff',
                   default=0.001,type=float,
                   help="exclude barcodes with a smaller fraction of "
                   "associated reads than CUTOFF, e.g. '0.01' excludes "
                   "barcodes with < 1%% of reads (default is 0.01%%)")
    p.add_argument('--sample-sheet',action="store",
                   dest="sample_sheet",default=None,
                   help="use an alternative sample sheet to the default "
                   "'custom_SampleSheet.csv' created on setup.")
    p.add_argument('--barcode-analysis-dir',action="store",
                   dest="barcode_analysis_dir",default=None,
                   help="specify subdirectory where barcode analysis will "
                   "be performed and outputs will be written")
    p.add_argument('--force',action='store_true',
                   dest="force",default=False,
                   help="discard and regenerate counts (by default existing "
                   "counts will be used)")
    add_runner_option(p)
    add_debug_option(p)
    p.add_argument('analysis_dir',metavar="ANALYSIS_DIR",nargs="?",
                   help="auto_process analysis directory (optional: defaults "
                   "to the current directory)")

def add_merge_fastq_dirs_command(cmdparser):
    """
    Create a parser for the 'merge_fastq_dirs' command
    """
    p = cmdparser.add_command('merge_fastq_dirs',
                              help="Combine bcl-to-fastq runs",
                              description="Automatically merge fastq "
                              "directories from multiple bcl-to-fastq runs "
                              "within ANALYSIS_DIR. Use this command if "
                              "'make_fastqs' step was run multiple times to "
                              "process subsets of lanes.")
    p.add_argument('--primary-unaligned-dir',action='store',
                   dest='unaligned_dir',default='bcl2fastq',
                   help="merge fastqs from additional bcl-to-fastq "
                   "directories into UNALIGNED_DIR. Original data will "
                   "be moved out of the way first. Defaults to 'bcl2fastq'.")
    p.add_argument('--output-dir',action='store',
                   dest='output_dir',default=None,
                   help="merge fastqs into OUTPUT_DIR (relative to "
                   "ANALYSIS_DIR). Defaults to UNALIGNED_DIR.")
    add_dry_run_option(p)
    add_debug_option(p)
    p.add_argument('analysis_dir',metavar="ANALYSIS_DIR",nargs="?",
                   help="auto_process analysis directory (optional: defaults "
                   "to the current directory)")

def add_import_project_command(cmdparser):
    """
    Create a parser for the 'import_project' command
    """
    p = cmdparser.add_command('import_project',
                              help="Import a project directory",
                              description="Copy a project directory "
                              "PROJECT_DIR into ANALYSIS_DIR.")
    add_debug_option(p)
    p.add_argument('analysis_dir',metavar="ANALYSIS_DIR",nargs="?",
                   help="auto_process analysis directory (optional: defaults "
                   "to the current directory)")
    p.add_argument('project_dir',metavar="PROJECT_DIR",
                   help="path to project directory to import")

def add_update_fastq_stats_command(cmdparser):
    """
    Create a parser for the 'update_fastq_stats' command
    """
    p = cmdparser.add_command('update_fastq_stats',
                              help="(Re)generate Fastq statistics",
                              description="(Re)generate statistics for fastq "
                              "files produced from 'make_fastqs'.")
    p.add_argument('--unaligned-dir',action='store',
                   dest='unaligned_dir',default='bcl2fastq',
                   help="explicitly set the (sub)directory with "
                   "bcl-to-fastq outputs")
    p.add_argument('--stats-file',action='store',
                   dest='stats_file',default=None,
                   help="specify output file for fastq statistics")
    p.add_argument('--per-lane-stats-file',action='store',
                   dest='per_lane_stats_file',default=None,
                   help="specify output file for per-lane statistics")
    p.add_argument('-a','--add',action="store_true",dest="add_data",
                   help="add new data from UNALIGNED_DIR to existing "
                   "statistics")
    p.add_argument('--force',action="store_true",dest="force",
                   help="force statistics to be regenerated even if "
                   "existing statistics files are newer than fastqs ")
    add_nprocessors_option(p,__settings.fastq_stats.nprocessors)
    add_runner_option(p)
    add_debug_option(p)
    p.add_argument('analysis_dir',metavar="ANALYSIS_DIR",nargs="?",
                   help="auto_process analysis directory (optional: defaults "
                   "to the current directory)")

# Commands

def config(args):
    """
    Implement functionality for 'config' command
    """
    if args.init:
        # Create new settings file and reload
        settings_file = locate_settings_file(create_from_sample=True)
        settings = Settings(settings_file=settings_file)
    else:
        settings = __settings
    if args.key_value or args.new_section:
        if args.new_section is not None:
            # Add new sections
            for new_section in args.new_section:
                try:
                    new_section = new_section.strip("'").strip('"')
                    print("Adding section '%s'" % new_section)
                    settings.add_section(new_section)
                except ValueError:
                    logging.error("Can't process '%s'" % args.section)
        if args.key_value is not None:
            # Update parameters
            for key_value in args.key_value:
                try:
                    i = key_value.index('=')
                    key = key_value[:i]
                    value = key_value[i+1:].strip("'").strip('"')
                    print("Setting '%s' to '%s'" % (key,value))
                    settings.set(key,value)
                except ValueError:
                    logging.error("Can't process '%s'" % args.key_value)
        # Save the updated settings to file
        settings.save()
    elif not args.init:
        # Report the current configuration settings
        paginate(settings.report_settings())

def setup(args):
    """
    Implement functionality for 'setup' command
    """
    d = AutoProcess()
    if args.unaligned_dir is None:
        d.setup(args.run_dir,
                analysis_dir=args.analysis_dir,
                sample_sheet=args.sample_sheet,
                extra_files=args.extra_files,
                unaligned_dir=args.unaligned_dir)

def params(args):
    """
    Implement functionality for 'params' command
    """
    analysis_dir = args.analysis_dir
    if not analysis_dir:
        analysis_dir = os.getcwd()
    d = AutoProcess(analysis_dir,allow_save_params=False)
    # Update the project-specific parameters
    if args.key_value is not None:
        for key_value in args.key_value:
            try:
                i = key_value.index('=')
                key = key_value[:i]
                value = key_value[i+1:].strip("'").strip('"')
                d.set_param(key,value)
            except ValueError:
                logging.error("Can't process '%s'" % args.key_value)
            d.save_parameters(force=True)
    else:
        d.print_params()

def metadata(args):
    """
    Implement functionality for 'metadata' command
    """
    analysis_dir = args.analysis_dir
    if not analysis_dir:
        analysis_dir = os.getcwd()
    d = AutoProcess(analysis_dir,allow_save_params=False)
    # Update the metadata associated with the analysis
    if args.update:
        d.update_metadata()
    if args.key_value is not None:
        for key_value in args.key_value:
            try:
                i = key_value.index('=')
                key = key_value[:i]
                value = key_value[i+1:].strip("'").strip('"')
                d.set_metadata(key,value)
            except ValueError:
                logging.error("Can't process '%s'" % args.key_value)
    if args.key_value or args.update:
        # Save updates to file
        d.save_metadata(force=True)
    else:
        # Only print values
        d.print_metadata()

def samplesheet(args):
    """
    Implement functionality for 'samplesheet' command
    """
    analysis_dir = args.analysis_dir
    if not analysis_dir:
        analysis_dir = os.getcwd()
    d = AutoProcess(analysis_dir)
    # Sample sheet operations
    if args.edit:
        # Manually edit the sample sheet
        d.edit_samplesheet()
    else:
        # Predict the outputs
        paginate(predict_outputs(
            sample_sheet_file=d.params.sample_sheet))

def make_fastqs(args):
    """
    Implement functionality for 'make_fastqs' command
    """
    analysis_dir = args.analysis_dir
    if not analysis_dir:
        analysis_dir = os.getcwd()
    d = AutoProcess(analysis_dir)
    # Deal with --no-lane-splitting
    if args.no_lane_splitting and args.use_lane_splitting:
        raise Exception("--no-lane-splitting and --use-lane-splitting "
                        "are mutually exclusive")
    elif args.no_lane_splitting:
        no_lane_splitting = True
    elif args.use_lane_splitting:
        no_lane_splitting = False
    else:
        no_lane_splitting = None
    # Deal with --[no-]create-empty-fastqs
    if args.create_empty_fastqs and args.no_create_empty_fastqs:
        raise Exception("--create-empty-fastqs and --no-create-empty-fastqs "
                        "are mutually exclusive")
    elif args.create_empty_fastqs:
        create_empty_fastqs = True
    elif args.no_create_empty_fastqs:
        create_empty_fastqs = False
    else:
        create_empty_fastqs = None
    # Deal with --lanes
    if args.lanes:
        # Build subsets
        only_include_lanes = list()
        lane_subsets = list()
        for lanes_spec in args.lanes:
            # Split the lane spec into components
            lanes_spec = lanes_spec.rstrip(':').split(':')
            print("Lanes spec: %s" % lanes_spec)
            # Extract lanes for this subset
            lanes_ = bcf_utils.parse_lanes(lanes_spec[0])
            # Add to set of lanes to include
            only_include_lanes.extend(lanes_)
            # Collect additional parameters
            subset_options = dict()
            # Check for initial protocol specification
            try:
                if '=' not in lanes_spec[1]:
                    subset_options['protocol'] = lanes_spec[1]
                    del(lanes_spec[1])
            except IndexError:
                pass
            # Remaining options should be 'key=value' pairs
            for key_value in lanes_spec[1:]:
                logging.debug("-- handling %s" % key_value)
                key,value = key_value.split('=')
                # Normalise the keys
                # Remove leading '-'s and replace other '-'s
                # with underscores
                # (e.g. '--no-trimming' => 'no_trimming')
                key = key.lstrip('-').replace('-','_')
                # Some command line options don't map
                # directly so update those
                if key == 'use_bases_mask':
                    key = 'bases_mask'
                elif key == 'well_list':
                    key = 'icell8_well_list'
                elif key in ('swap_i1_and_i2',
                             'reverse_complement'):
                    key = "icell8_atac_%s" % key
                # Deal with True/False values
                # Ignore case and also allow 'yes' and 'no'
                if value.lower() in ('true','yes'):
                    value = True
                elif value.lower() in ('false','no'):
                    value = False
                subset_options[key] = value
            # Add the subset to the list
            lane_subsets.append(subset(lanes_,**subset_options))
        logging.debug("Lane subsets: %s" % lane_subsets)
        only_include_lanes = sorted(only_include_lanes)
    else:
        # No lane subsets
        only_include_lanes = None
        lane_subsets = None
    # Handle job runner specification
    if args.runner is not None:
        runner = fetch_runner(args.runner)
    else:
        runner = None
    # Do the make_fastqs step
    d.make_fastqs(
        protocol=args.protocol,
        nprocessors=args.nprocessors,
        runner=runner,
        force_copy_of_primary_data=args.force_copy,
        generate_stats=(not args.no_stats),
        analyse_barcodes=(not args.no_barcode_analysis),
        require_bcl2fastq_version=args.bcl2fastq_version,
        unaligned_dir=args.out_dir,
        sample_sheet=args.sample_sheet,
        bases_mask=args.bases_mask,
        lanes=only_include_lanes,lane_subsets=lane_subsets,
        icell8_well_list=args.icell8_well_list,
        icell8_swap_i1_and_i2=args.icell8_swap_i1_and_i2,
        icell8_reverse_complement=args.icell8_reverse_complement,
        platform=args.platform,
        no_lane_splitting=no_lane_splitting,
        trim_adapters=(not args.no_adapter_trimming),
        minimum_trimmed_read_length=args.minimum_trimmed_read_length,
        mask_short_adapter_reads=args.mask_short_adapter_reads,
        adapter_sequence=args.adapter_sequence,
        adapter_sequence_read2=args.adapter_sequence_read2,
        create_fastq_for_index_read=args.create_fastq_for_index_read,
        stats_file=args.stats_file,
        per_lane_stats_file=args.per_lane_stats_file,
        barcode_analysis_dir=args.barcode_analysis_dir,
        create_empty_fastqs=create_empty_fastqs,
        cellranger_jobmode=args.job_mode,
        cellranger_mempercore=args.mem_per_core,
        cellranger_maxjobs=args.max_jobs,
        cellranger_jobinterval=args.job_interval,
        cellranger_localcores=args.local_cores,
        cellranger_localmem=args.local_mem,
        cellranger_ignore_dual_index=args.ignore_dual_index,
        working_dir=args.working_dir,
        verbose=args.verbose)

def setup_analysis_dirs(args):
    """
    Implement functionality for 'setup_analysis_dirs' command
    """
    analysis_dir = args.analysis_dir
    if not analysis_dir:
        analysis_dir = os.getcwd()
    d = AutoProcess(analysis_dir)
    d.setup_analysis_dirs(unaligned_dir=args.unaligned_dir,
                          ignore_missing_metadata=
                          args.ignore_missing_metadata,
                          undetermined_project=args.undetermined,
                          short_fastq_names=args.short_fastq_names,
                          link_to_fastqs=args.link_to_fastqs)

def run_qc(args):
    """
    Implement functionality for 'run_qc' command
    """
    analysis_dir = args.analysis_dir
    if not analysis_dir:
        analysis_dir = os.getcwd()
    d = AutoProcess(analysis_dir)
    # Handle 10x transcriptomes
    cellranger_transcriptomes = dict()
    if args.cellranger_transcriptomes:
        for transcriptome in args.cellranger_transcriptomes:
            organism,reference =  transcriptome.split('=')
            cellranger_transcriptomes[organism] = reference
    # Handle 10x pre-mRNA references
    cellranger_premrna_references = dict()
    if args.cellranger_premrna_references:
        for premrna_reference in args.cellranger_premrna_references:
            organism,reference =  premrna_reference.split('=')
            cellranger_premrna_references[organism] = reference
    # Handle job runner specification
    if args.runner is not None:
        runner = fetch_runner(args.runner)
    else:
        runner = None
    # Do the run_qc step
    retcode = d.run_qc(projects=args.project_pattern,
                       max_jobs=args.max_jobs,
                       ungzip_fastqs=args.ungzip_fastqs,
                       fastq_screen_subset=args.subset,
                       nthreads=args.nthreads,
                       fastq_dir=args.fastq_dir,
                       qc_dir=args.qc_dir,
                       cellranger_chemistry=
                       args.cellranger_chemistry,
                       cellranger_transcriptomes=
                       cellranger_transcriptomes,
                       cellranger_premrna_references=
                       cellranger_premrna_references,
                       report_html=args.html_file,
                       runner=runner,
                       working_dir=args.working_dir,
                       verbose=args.verbose)
    sys.exit(retcode)

def publish_qc(args):
    """
    Implement functionality for 'publish_qc' command
    """
    analysis_dir = args.analysis_dir
    if not analysis_dir:
        analysis_dir = os.getcwd()
    # Handle job runner specification
    if args.runner is not None:
        runner = fetch_runner(args.runner)
    else:
        runner = None
    # Do the publish_qc step
    d = AutoProcess(analysis_dir)
    d.publish_qc(projects=args.project_pattern,
                 location=args.qc_dir,
                 use_hierarchy=(args.use_hierarchy == 'yes'),
                 exclude_zip_files=(args.exclude_zip_files == 'yes'),
                 ignore_missing_qc=args.ignore_missing_qc,
                 regenerate_reports=args.regenerate_reports,
                 force=args.force,legacy=args.legacy_mode,
                 runner=runner)

def archive(args):
    """
    Implement functionality for 'archive' command
    """
    analysis_dir = args.analysis_dir
    if not analysis_dir:
        analysis_dir = os.getcwd()
    d = AutoProcess(analysis_dir)
    retcode = d.archive(archive_dir=args.archive_dir,
                        platform=args.platform,
                        year=args.year,
                        group=args.group,
                        perms=args.permissions,
                        final=args.final,
                        force=args.force,
                        dry_run=args.dry_run)
    sys.exit(retcode)

def report(args):
    """
    Implement functionality for 'report' command
    """
    analysis_dir = args.analysis_dir
    if not analysis_dir:
        analysis_dir = os.getcwd()
    d = AutoProcess(analysis_dir,allow_save_params=False)
    if args.logging:
        mode = ReportingMode.CONCISE
    elif args.summary:
        mode = ReportingMode.SUMMARY
    elif args.projects:
        mode = ReportingMode.PROJECTS
    else:
        mode = None
    if args.fields:
        fields = str(args.fields).split(',')
    elif args.template:
        try:
            fields = str(__settings.reporting_templates[args.template]).split(',')
        except KeyError:
            print("Template '%s' not in list of custom templates:" %
                  args.template)
            if __settings.reporting_templates:
                for template in __settings.reporting_templates:
                    print("- %s" % template)
                else:
                    print("- (no templates defined)")
                logging.critical("--template: unrecognised template '%s'"
                                 % args.template)
                sys.exit(1)
    else:
        fields = None
    d.report(mode=mode,fields=fields,out_file=args.out_file)

def readme(args):
    """
    Implement functionality for 'readme' command
    """
    analysis_dir = args.analysis_dir
    if not analysis_dir:
        analysis_dir = os.getcwd()
    d = AutoProcess(analysis_dir)
    if args.init:
        d.init_readme()
    else:
        if d.readme_file is None:
            print("No README file for %s" % d.analysis_dir)
        else:
            if args.message:
                message = '\n'.join(
                    bcf_utils.split_into_lines(
                        "[%s]\n%s" % (time.ctime(),
                                      str(args.message)),
                        80,sympathetic=True))
                with open(d.readme_file,'a') as fp:
                    fp.write("\n%s\n" % message)
            elif args.edit:
                d.edit_readme()
            else:
                print(d.readme_file)
                if args.view:
                    paginate(open(d.readme_file,'r').read())

def clone(args):
    """
    Implement functionality for 'clone' command
    """
    analysis_dir = args.analysis_dir
    if not analysis_dir:
        analysis_dir = os.getcwd()
    d = AutoProcess(analysis_dir,allow_save_params=False)
    d.clone(args.clone_dir,
            copy_fastqs=args.copy_fastqs,
            exclude_projects=args.exclude_projects)

def analyse_barcodes(args):
    """
    Implement functionality for 'analyse_barcodes' command
    """
    analysis_dir = args.analysis_dir
    if not analysis_dir:
        analysis_dir = os.getcwd()
    d = AutoProcess(analysis_dir)
    # Deal with --lanes
    if args.lanes is not None:
        lanes = bcf_utils.parse_lanes(args.lanes)
    else:
        lanes = None
    # Handle job runner specification
    if args.runner is not None:
        runner = fetch_runner(args.runner)
    else:
        runner = None
    # Do barcode analysis
    d.analyse_barcodes(unaligned_dir=args.unaligned_dir,
                       lanes=lanes,
                       mismatches=args.mismatches,
                       cutoff=args.cutoff,
                       sample_sheet=args.sample_sheet,
                       barcode_analysis_dir=args.barcode_analysis_dir,
                       runner=runner,
                       force=args.force)

def merge_fastq_dirs(args):
    """
    Implement functionality for 'merge_fastq_dirs' command
    """
    analysis_dir = args.analysis_dir
    if not analysis_dir:
        analysis_dir = os.getcwd()
    d = AutoProcess(analysis_dir)
    d.merge_fastq_dirs(args.unaligned_dir,
                       output_dir=args.output_dir,
                       dry_run=args.dry_run)

def import_project(args):
    """
    Implement functionality for 'import_project' command
    """
    analysis_dir = args.analysis_dir
    if not analysis_dir:
        analysis_dir = os.getcwd()
    d = AutoProcess(analysis_dir)
    d.import_project(args.project_dir)

def update_fastq_stats(args):
    """
    Implement functionality for 'update_fastq_stats' command
    """
    analysis_dir = args.analysis_dir
    if not analysis_dir:
        analysis_dir = os.getcwd()
    d = AutoProcess(analysis_dir)
    # Handle job runner specification
    if args.runner is not None:
        runner = fetch_runner(args.runner)
    else:
        runner = None
    # Do the update
    d.update_fastq_stats(
        unaligned_dir=args.unaligned_dir,
        stats_file=args.stats_file,
        per_lane_stats_file=args.per_lane_stats_file,
        add_data=args.add_data,
        force=args.force,
        nprocessors=args.nprocessors,
        runner=runner)

# Supporting functions

def add_modulefiles_option(p,dest='modulefiles',default=None):
    """
    Add a --modulefiles argument to a parser instance

    The values can be accessed via the 'modulefiles' property of the
    parser.
    """
    add_arg(p,'--modulefiles',action='store',
            dest=dest,default=default,
            help="Specify comma-separated list of environment "
            "modules to load before executing commands (overrides "
            "any modules specified in the global settings)")

def set_debug(debug_flag):
    """
    Turn on debug output
    """
    if debug_flag: logging.getLogger().setLevel(logging.DEBUG)

# Main function

def main():
    """
    """
    # Set up the command line parser
    p = CommandParser(
        description="Automated processing & QC for Illumina sequencing data",
        version="%prog "+__version__,
        subparser=argparse.ArgumentParser)

    # Add commands
    add_config_command(p)
    add_setup_command(p)
    add_params_command(p)
    add_metadata_command(p)
    add_samplesheet_command(p)
    add_make_fastqs_command(p)
    add_setup_analysis_dirs_command(p)
    add_run_qc_command(p)
    add_publish_qc_command(p)
    add_archive_command(p)
    add_report_command(p)
    add_readme_command(p)
    add_clone_command(p)
    add_analyse_barcodes_command(p)
    add_merge_fastq_dirs_command(p)
    add_import_project_command(p)
    add_update_fastq_stats_command(p)

    # Map commands to functions
    commands = {
        'config': config,
        'setup': setup,
        'params': params,
        'metadata': metadata,
        'samplesheet': samplesheet,
        'make_fastqs': make_fastqs,
        'setup_analysis_dirs': setup_analysis_dirs,
        'run_qc': run_qc,
        'publish_qc': publish_qc,
        'archive': archive,
        'report': report,
        'readme': readme,
        'clone': clone,
        'analyse_barcodes': analyse_barcodes,
        'merge_fastq_dirs': merge_fastq_dirs,
        'import_project': import_project,
        'update_fastq_stats': update_fastq_stats,
    }
    
    # Process command line
    cmd,args = p.parse_args()

    # Report name and version
    print("%s version %s" % (os.path.basename(sys.argv[0]),__version__))

    # Turn on debugging?
    set_debug(args.debug)

    # Allow saving of parameters?
    try:
        allow_save = not args.no_save
    except AttributeError:
        allow_save = True

    # Set up environment modules
    modulefiles = None
    try:
        modulefiles = args.modulefiles
    except AttributeError:
        pass
    if modulefiles is None:
        try:
            modulefiles = __settings.modulefiles[cmd]
        except KeyError:
            # No environment modules specified
            pass
    if modulefiles is not None:
        for modulefile in modulefiles.split(','):
            envmod_load(modulefile)

    # Locate and run the requested command
    try:
        commands[cmd](args)
    except KeyError:
        raise Exception("Unrecognised command: '%s'" % cmd)
