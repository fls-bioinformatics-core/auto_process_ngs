#!/usr/bin/env python
#
#     auto_process.py: automated processing of Illumina sequence data
#     Copyright (C) University of Manchester 2013-15 Peter Briggs
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
    config
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

Additional commands are available:

    clone
    analyse_barcodes
    merge_fastq_dirs
    update_fastq_stats

but these are not part of the standard workflow - they are used for
special cases and testing.

"""

#######################################################################
# Imports
#######################################################################

import logging
#logging.basicConfig(format='[%(asctime)s] %(levelname)8s: %(message)s')
logging.basicConfig(format='%(levelname) 8s: %(message)s')

import sys
import os
import subprocess
import optparse
import bcftbx.utils as bcf_utils
from bcftbx.cmdparse import CommandParser
from bcftbx.cmdparse import add_debug_option
from bcftbx.cmdparse import add_no_save_option
from bcftbx.cmdparse import add_dry_run_option
from bcftbx.cmdparse import add_nprocessors_option
from bcftbx.cmdparse import add_runner_option
import auto_process_ngs
import auto_process_ngs.settings
import auto_process_ngs.envmod as envmod
from auto_process_ngs.auto_processor import AutoProcess

__version__ = auto_process_ngs.get_version()
__settings = auto_process_ngs.settings.Settings()

#######################################################################
# Functions
#######################################################################

# Command line parsers

def add_setup_command(cmdparser):
    """Create a parser for the 'setup' command
    """
    p = cmdparser.add_command('setup',help="Set up a new analysis directory",
                              usage="%prog setup [OPTIONS] DIR",
                              description="Set up automatic processing of Illumina "
                              "sequencing data from DIR.")
    p.add_option('--sample-sheet',action='store',dest='sample_sheet',default=None,
                 help="Copy sample sheet file from name and location SAMPLE_SHEET "
                 "(default is to look for SampleSheet.csv inside DIR)")
    p.add_option('--fastq-dir',action='store',dest='fastq_dir',default=None,
                 help="Import fastq.gz files from FASTQ_DIR (which should be a "
                 "subdirectory of DIR with the same structure as that produced "
                 "by CASAVA/bcl2fastq i.e. 'Project_<name>/Sample_<name>/<fastq>')")
    p.add_option('--analysis-dir',action='store',dest='analysis_dir',default=None,
                 help="Make new directory called ANALYSIS_DIR (otherwise default is "
                 "'DIR_analysis')")
    add_debug_option(p)

def add_config_command(cmdparser):
    """Create a parser for the 'config' command
    """
    p  = cmdparser.add_command('config',help="Query and change global configuration",
                               usage="%prog config [OPTIONS] [ANALYSIS_DIR]",
                               description="Query and change global configuration.")
    p.add_option('--init',action='store_true',dest='init',default=False,
                 help="Create a new configuration file from the sample.")
    p.add_option('--set',action='append',dest='key_value',default=None,
                 help="Set the value of a parameter. KEY_VALUE should be of the form "
                 "'<param>=<value>'. Multiple --set options can be specified.")
    add_debug_option(p)
    # Deprecated options
    deprecated = optparse.OptionGroup(p,'Deprecated/defunct options')
    deprecated.add_option('--show',action='store_true',dest='show',default=False,
                          help="Show the values of parameters and settings (does "
                          "nothing; use 'config' with no options to display settings)")
    p.add_option_group(deprecated)

def add_params_command(cmdparser):
    """Create a parser for the 'params' command
    """
    p  = cmdparser.add_command('params',help="Query and change project parameters",
                               usage="%prog config [OPTIONS] [ANALYSIS_DIR]",
                               description="Query and change processing parameters "
                               "and settings for ANALYSIS_DIR.")
    p.add_option('--set',action='append',dest='key_value',default=None,
                 help="Set the value of a parameter. KEY_VALUE should be of the form "
                 "'<param>=<value>'. Multiple --set options can be specified.")
    add_debug_option(p)

def add_make_fastqs_command(cmdparser):
    """Create a parser for the 'make_fastqs' command
    """
    p = cmdparser.add_command('make_fastqs',help="Run Fastq generation",
                              usage="%prog make_fastqs [OPTIONS] [ANALYSIS_DIR]",
                              description="Generate fastq files from raw bcl files "
                              "produced by Illumina sequencer within ANALYSIS_DIR.")
    # General options
    add_no_save_option(p)
    add_debug_option(p)
    # Primary data management
    primary_data = optparse.OptionGroup(p,'Primary data management')
    primary_data.add_option('--only-fetch-primary-data',action='store_true',
                            dest='only_fetch_primary_data',default=False,
                            help="only fetch the primary data, don't perform any other "
                            "operations")
    primary_data.add_option('--skip-rsync',action='store_true',
                            dest='skip_rsync',default=False,
                            help="don't rsync the primary data at the beginning of processing")
    primary_data.add_option('--remove-primary-data',action='store_true',
                            dest='remove_primary_data',default=False,
                            help="Delete the primary data at the end of processing (default "
                            "is to keep data)")
    p.add_option_group(primary_data)
    # Options to control bcl2fastq
    bcl_to_fastq = optparse.OptionGroup(p,'Bcl-to-fastq options')
    bcl_to_fastq.add_option('--skip-bcl2fastq',action='store_true',
                            dest='skip_bcl2fastq',default=False,
                            help="don't run the Fastq generation step")
    bcl_to_fastq.add_option('--output-dir',action='store',
                            dest='unaligned_dir',default=None,
                            help="explicitly set the output (sub)directory for bcl-to-fastq "
                            "conversion (overrides default)")
    bcl_to_fastq.add_option('--use-bases-mask',action="store",
                            dest="bases_mask",default=None,
                            help="explicitly set the bases-mask string to indicate how each "
                            "cycle should be used in the bcl-to-fastq conversion (overrides "
                            "default)")
    bcl_to_fastq.add_option('--sample-sheet',action="store",
                            dest="sample_sheet",default=None,
                            help="use an alternative sample sheet to the default "
                            "'custom_SampleSheet.csv' created on setup.")
    bcl_to_fastq.add_option('--ignore-missing-bcl',action='store_true',
                            dest='ignore_missing_bcl',default=False,
                            help="use the --ignore-missing-bcl option for bcl2fastq (treat "
                            "missing bcl files as no call)")
    bcl_to_fastq.add_option('--ignore-missing-stats',action='store_true',
                            dest='ignore_missing_stats',default=False,
                            help="use the --ignore-missing-stats option for bcl2fastq (fill "
                            "in with zeroes when *.stats files are missing)")
    add_nprocessors_option(bcl_to_fastq,None,
                           default_display=__settings.bcl2fastq.nprocessors)
    add_runner_option(bcl_to_fastq)
    p.add_option_group(bcl_to_fastq)
    # Statistics
    statistics = optparse.OptionGroup(p,'Statistics generation')
    statistics.add_option('--stats-file',action='store',
                          dest='stats_file',default=None,
                          help="specify output file for fastq statistics")
    statistics.add_option('--no-stats',action='store_true',
                          dest='no_stats',default=False,
                          help="don't generate statistics file; use 'update_fastq_stats' "
                          "command to (re)generate statistics")
    p.add_option_group(statistics)
    # Deprecated options
    deprecated = optparse.OptionGroup(p,'Deprecated/defunct options')
    deprecated.add_option('--keep-primary-data',action='store_true',
                          dest='keep_primary_data',default=False,
                          help="don't delete the primary data at the end of processing "
                          "(does nothing; primary data is kept by default unless "
                          "--remove-primary-data is specified)")
    deprecated.add_option('--generate-stats',action='store_true',
                          dest='generate_stats',default=False,
                          help="(re)generate statistics for fastq files (does nothing; "
                          "statistics are generated by default unless suppressed by "
                          "--no-stats)")
    deprecated.add_option('--report-barcodes',action='store_true',
                          dest='report_barcodes',default=False,
                          help="analyse and report barcode indices for all lanes after "
                          "generating fastq files (deprecated: use the "
                          "'analyse_barcodes' command instead)")
    deprecated.add_option('--barcodes-file',action='store',
                          dest='barcodes_file',default=None,
                          help="specify output file for barcode analysis report "
                          "(deprecated: use the 'analyse_barcodes' command instead)")
    p.add_option_group(deprecated)

def add_setup_analysis_dirs_command(cmdparser):
    """Create a parser for the 'setup_analysis_dirs' command
    """
    p = cmdparser.add_command('setup_analysis_dirs',help="Create project subdirectories",
                              usage="%prog setup_analysis_dirs [OPTIONS] [ANALYSIS_DIR]",
                              description="Create analysis subdirectories for projects "
                              "defined in projects.info file in ANALYSIS_DIR.")
    p.add_option('--ignore-missing-metadata',action='store_true',
                 dest='ignore_missing_metadata',default=False,
                 help="force creation of project directories even if metadata is not "
                 "set (default is to fail if metadata is missing)")
    p.add_option('--short-fastq-names',action='store_true',
                 dest='short_fastq_names',default=False,
                 help="shorten fastq file names when copying or linking from project "
                 "directory (default is to keep long names from bcl2fastq)")
    p.add_option('--link-to-fastqs',action='store_true',
                 dest='link_to_fastqs',default=False,
                 help="create symbolic links to original fastqs from project directory "
                 "(default is to make hard links)")
    add_debug_option(p)

def add_run_qc_command(cmdparser):
    """Create a parser for the 'run_qc' command
    """
    p = cmdparser.add_command('run_qc',help="Run QC procedures",
                              usage="%prog run_qc [OPTIONS] [ANALYSIS_DIR]",
                              description="Run QC procedures for sequencing projects in "
                              "ANALYSIS_DIR.")
    max_concurrent_jobs = __settings.general.max_concurrent_jobs
    p.add_option('--projects',action='store',
                 dest='project_pattern',default=None,
                 help="simple wildcard-based pattern specifying a subset of projects "
                 "and samples to run the QC on. PROJECT_PATTERN should be of the form "
                 "'pname[/sname]', where 'pname' specifies a project (or set of "
                 "projects) and 'sname' optionally specifies a sample (or set of "
                 "samples).")
    p.add_option('--ungzip-fastqs',action='store_true',dest='ungzip_fastqs',
                 help="create decompressed copies of fastq.gz files")
    p.add_option('--max-jobs',action='store',
                 dest='max_jobs',default=max_concurrent_jobs,type='int',
                 help="explicitly specify maximum number of concurrent QC jobs to run "
                 "(default %s, change in settings file)" % max_concurrent_jobs)
    add_runner_option(p)
    add_debug_option(p)
    # Deprecated options
    deprecated = optparse.OptionGroup(p,'Deprecated/defunct options')
    deprecated.add_option('--no-ungzip-fastqs',action='store_true',dest='no_ungzip_fastqs',
                          help="don't create uncompressed copies of fastq.gz files "
                          "(does nothing; this is now the default, use --ungzip-fastqs "
                          "to turn on decompression)")
    p.add_option_group(deprecated)

def add_publish_qc_command(cmdparser):
    """Create a parser for the 'publish_qc' command
    """
    p = cmdparser.add_command('publish_qc',help="Copy QC reports to publication area",
                              usage="%prog publish_qc [OPTIONS] [ANALYSIS_DIR]",
                              description="Copy QC reports from ANALYSIS_DIR to local "
                              "or remote directory (e.g. web server). By default existing "
                              "QC reports will be copied without further checking; if no "
                              "report is found then QC results will be verified and a "
                              "report generated first.")
    p.add_option('--projects',action='store',
                 dest='project_pattern',default=None,
                 help="simple wildcard-based pattern specifying a subset of projects "
                 "and samples to publish the QC for. PROJECT_PATTERN can specify a "
                 "single project, or a set of projects.")
    p.add_option('--qc_dir',action='store',
                 dest='qc_dir',default=None,
                 help="specify target directory to copy QC reports to. QC_DIR can "
                 "be a local directory, or a remote location in the form "
                 "'[[user@]host:]directory'. Overrides the default settings.")
    p.add_option('--ignore-missing-qc',action='store_true',
                 dest='ignore_missing_qc',default=False,
                 help="skip projects where QC results are missing or can't be verified, "
                 "or where reports can't be generated.")
    p.add_option('--regenerate-reports',action='store_true',
                 dest='regenerate_reports',default=False,
                 help="attempt to regenerate existing QC reports")
    add_debug_option(p)

def add_archive_command(cmdparser):
    """Create a parser for the 'archive' command
    """
    p = cmdparser.add_command('archive',help="Copy analyses to 'archive' area",
                              usage="%prog archive [OPTIONS] [ANALYSIS_DIR]",
                              version="%prog "+__version__,
                              description="Copy sequencing analysis data directory "
                              "ANALYSIS_DIR to 'archive' destination.")
    p.add_option('--archive_dir',action='store',
                 dest='archive_dir',default=None,
                 help="specify top-level archive directory to copy data under. "
                 "ARCHIVE_DIR can be a local directory, or a remote location in the "
                 "form '[[user@]host:]directory'. Overrides the default settings.")
    p.add_option('--platform',action='store',
                 dest='platform',default=None,
                 help="specify the platform e.g. 'hiseq', 'miseq' etc (overrides "
                 "automatically determined platform, if any). Use 'other' for cases "
                 "where the platform is unknown.")
    p.add_option('--year',action='store',
                 dest='year',default=None,
                 help="specify the year e.g. '2014' (default is the current year)")
    default_group = __settings.archive.group
    p.add_option('--group',action='store',dest='group',default=default_group,
                 help="specify the name of group for the archived files NB only works "
                 "when the archive is a local directory (default: %s)" % default_group)
    default_chmod = __settings.archive.chmod
    p.add_option('--chmod',action='store',dest='chmod',default=default_chmod,
                 help="specify chmod operations for the archived files (default: "
                 "%s)" % default_chmod)
    p.add_option('--force',action='store_true',dest='force',default=False,
                 help="perform archiving operation even if key metadata items are "
                 "not set")
    add_dry_run_option(p)
    add_debug_option(p)

def add_report_command(cmdparser):
    """Create a parser for the 'report' command
    """
    p  = cmdparser.add_command('report',help="Generate reporting information",
                               usage="%prog report [OPTIONS] [ANALYSIS_DIR]",
                               description="Report information on processed Illumina "
                               "sequence data in ANALYSIS_DIR.")
    p.add_option('--logging',action='store_true',dest='logging',default=False,
                 help="print short report suitable for logging file")
    p.add_option('--summary',action='store_true',dest='summary',default=False,
                 help="print full report suitable for bioinformaticians")
    p.add_option('--projects',action='store_true',dest='projects',default=False,
                 help="print tab-delimited line (one per project) suitable for "
                 "injection into a spreadsheet")
    p.add_option('--full',action='store_true',dest='full',default=False,
                 help="print summary report suitable for record-keeping")
    add_debug_option(p)

def add_clone_command(cmdparser):
    """Create a parser for the 'clone' command
    """
    p = cmdparser.add_command('clone',help="Make a copy of an analysis directory",
                              usage="%prog clone [OPTIONS] DIR CLONE_DIR",
                              description="Make a copy of an existing auto_processed analysis "
                              "directory DIR, in a new directory CLONE_DIR. The clone will "
                              "not include any project directories, but will copy the "
                              "projects.info file.")
    p.add_option('--copy-fastqs',action='store_true',dest='copy_fastqs',default=False,
                 help="Copy fastq.gz files from DIR into DIR2 (default is to make a "
                 "link to the bcl-to-fastq directory)")
    add_debug_option(p)

def add_analyse_barcodes_command(cmdparser):
    """Create a parser for the 'analyse_barcodes' command
    """
    p = cmdparser.add_command('analyse_barcodes',help="Analyse index (barcode) sequences",
                              usage="%prog analyse_barcodes [OPTIONS] [ANALYSIS_DIR]",
                              version="%prog "+__version__,
                              description="Analyse barcode sequences for fastq files "
                              "in specified lanes in ANALYSIS_DIR, and report the most "
                              "common barcodes found across all reads from each lane.")
    p.add_option('--unaligned-dir',action='store',
                 dest='unaligned_dir',default='bcl2fastq',
                 help="explicitly set the (sub)directory with bcl-to-fastq outputs")
    p.add_option('--lanes',action='store',
                 dest='lanes',default=None,
                 help="specify which lanes to analyse barcodes for (default is to do "
                 "analysis for all lanes).")
    p.add_option('--truncate',action='store',
                 dest='length',default=None,type='int',
                 help="truncate sample sheet barcodes to LENGTH for barcode analysis "
                 "(default is to use full barcodes)")
    add_nprocessors_option(p,1)
    add_runner_option(p)
    add_debug_option(p)

def add_merge_fastq_dirs_command(cmdparser):
    """Create a parser for the 'merge_fastq_dirs' command
    """
    p = cmdparser.add_command('merge_fastq_dirs',help="Combine bcl-to-fastq runs",
                              usage="%prog merge_fastq_dirs [OPTIONS] [ANALYSIS_DIR]",
                              description="Automatically merge fastq directories from "
                              "multiple bcl-to-fastq runs within ANALYSIS_DIR. Use this "
                              "command if 'make_fastqs' step was run multiple times to "
                              "process subsets of lanes.")
    p.add_option('--primary-unaligned-dir',action='store',
                 dest='unaligned_dir',default='bcl2fastq',
                 help="merge fastqs from additional bcl-to-fastq directories into "
                 "UNALIGNED_DIR. Original data will be moved out of the way first. "
                 "Defaults to 'bcl2fastq'.")
    add_dry_run_option(p)
    add_debug_option(p)

def add_update_fastq_stats_command(cmdparser):
    """Create a parser for the 'update_fastq_stats' command
    """
    p = cmdparser.add_command('update_fastq_stats',help="(Re)generate Fastq statistics",
                              usage="%prog update_fastq_stats [OPTIONS] [ANALYSIS_DIR]",
                              description="(Re)generate statistics for fastq "
                              "files produced from 'make_fastqs'.")
    p.add_option('--unaligned-dir',action='store',
                 dest='unaligned_dir',default='bcl2fastq',
                 help="explicitly set the (sub)directory with bcl-to-fastq outputs")
    p.add_option('--stats-file',action='store',
                 dest='stats_file',default=None,
                 help="specify output file for fastq statistics")
    add_nprocessors_option(p,__settings.fastq_stats.nprocessors)
    add_runner_option(p)
    add_debug_option(p)

def set_debug(debug_flag):
    """Turn on debug output
    """
    if debug_flag: logging.getLogger().setLevel(logging.DEBUG)

#######################################################################
# Main program
#######################################################################

if __name__ == "__main__":

    # Set up the command line parser
    p = CommandParser(description="Automated processing & QC for Illumina sequence data.",
                      version="%prog "+__version__)
    # Add commands
    add_config_command(p)
    add_setup_command(p)
    add_params_command(p)
    add_make_fastqs_command(p)
    add_setup_analysis_dirs_command(p)
    add_run_qc_command(p)
    add_publish_qc_command(p)
    add_archive_command(p)
    add_report_command(p)
    add_clone_command(p)
    add_analyse_barcodes_command(p)
    add_merge_fastq_dirs_command(p)
    add_update_fastq_stats_command(p)
    
    # Process remaining command line
    cmd,options,args = p.parse_args()

    # Report name and version
    print "%s version %s" % (os.path.basename(sys.argv[0]),__version__)

    # Turn on debugging?
    set_debug(options.debug)

    # Allow saving of parameters?
    try:
        allow_save = not options.no_save
    except AttributeError:
        allow_save = True

    # Set up environment modules
    try:
        modulefiles = __settings.modulefiles[cmd]
        if modulefiles is not None:
            for modulefile in modulefiles.split(','):
                envmod.load(modulefile)
    except KeyError:
        # No environment modules specified
        pass

    # Setup the processing object and run the requested command
    if cmd == 'config':
        if options.init:
            # Create new settings file and reload
            auto_process_ngs.settings.locate_settings_file()
            __settings = auto_process_ngs.settings.Settings()
        if options.key_value is not None:
            # Update parameters
            for key_value in options.key_value:
                try:
                    i = key_value.index('=')
                    key = key_value[:i]
                    value = key_value[i+1:].strip("'").strip('"')
                    print "Setting '%s' to '%s'" % (key,value)
                    __settings.set(key,value)
                except ValueError:
                    logging.error("Can't process '%s'" % options.key_value)
            __settings.save()
        else:
            # Report the current configuration settings
            __settings.report_settings()
    elif cmd == 'setup':
        if len(args) != 1:
            sys.stderr.write("Need to supply a data source location\n")
            sys.exit(1)
        d = AutoProcess()
        if options.fastq_dir is None:
            d.setup(args[0],
                    analysis_dir=options.analysis_dir,
                    sample_sheet=options.sample_sheet)
        else:
            d.setup_from_fastq_dir(args[0],options.fastq_dir)
    elif cmd == 'clone':
        if len(args) != 2:
            sys.stderr.write("Need to supply an existing analysis dir and directory for "
                             "the copy\n")
            sys.exit(1)
        d = AutoProcess(args[0])
        d.clone(args[1],copy_fastqs=options.copy_fastqs)
    else:
        # For other options check if an analysis
        # directory was specified
        if len(args) > 0:
            analysis_dir = args[0]
        else:
            analysis_dir = os.getcwd()
        d = AutoProcess(analysis_dir,allow_save_params=allow_save)
        # Run the specified stage
        if cmd == 'make_fastqs':
            # Do the make_fastqs step
            d.make_fastqs(skip_rsync=options.skip_rsync,
                          nprocessors=options.nprocessors,
                          runner=options.runner,
                          remove_primary_data=options.remove_primary_data,
                          ignore_missing_bcl=options.ignore_missing_bcl,
                          ignore_missing_stats=options.ignore_missing_stats,
                          generate_stats=(not options.no_stats),
                          unaligned_dir=options.unaligned_dir,
                          sample_sheet=options.sample_sheet,
                          bases_mask=options.bases_mask,
                          stats_file=options.stats_file,
                          report_barcodes=options.report_barcodes,
                          barcodes_file=options.barcodes_file,
                          skip_bcl2fastq=options.skip_bcl2fastq,
                          only_fetch_primary_data=options.only_fetch_primary_data)
        elif cmd == 'merge_fastq_dirs':
            d.merge_fastq_dirs(options.unaligned_dir,
                               dry_run=options.dry_run)
        elif cmd == 'update_fastq_stats':
            d.generate_stats(unaligned_dir=options.unaligned_dir,
                             stats_file=options.stats_file,
                             nprocessors=options.nprocessors,
                             runner=options.runner)
        elif cmd == 'analyse_barcodes':
            if options.lanes is not None:
                lanes = options.lanes.split(',')
            else:
                lanes = None
            d.analyse_barcodes(unaligned_dir=options.unaligned_dir,
                               lanes=lanes,
                               truncate_barcodes=options.length,
                               nprocessors=options.nprocessors,
                               runner=options.runner)
        elif cmd == 'setup_analysis_dirs':
            d.setup_analysis_dirs(ignore_missing_metadata=
                                  options.ignore_missing_metadata,
                                  short_fastq_names=options.short_fastq_names,
                                  link_to_fastqs=options.link_to_fastqs)
        elif cmd == 'run_qc':
            # Do the make_fastqs step
            retcode = d.run_qc(projects=options.project_pattern,
                               max_jobs=options.max_jobs,
                               ungzip_fastqs=options.ungzip_fastqs,
                               runner=options.runner)
            sys.exit(retcode)
        elif cmd == 'params':
            # Update the project-specific parameters
            if options.key_value is not None:
                for key_value in options.key_value:
                    try:
                        i = key_value.index('=')
                        key = key_value[:i]
                        value = key_value[i+1:].strip("'").strip('"')
                        d.set_param(key,value)
                    except ValueError:
                        logging.error("Can't process '%s'" % options.key_value)
            else:
                d.show_settings()
        elif cmd == 'archive':
            retcode = d.copy_to_archive(archive_dir=options.archive_dir,
                                        platform=options.platform,
                                        year=options.year,
                                        dry_run=options.dry_run,
                                        group=options.group,
                                        chmod=options.chmod,
                                        force=options.force)
            sys.exit(retcode)
        elif cmd == 'publish_qc':
            d.publish_qc(projects=options.project_pattern,
                         location=options.qc_dir,
                         ignore_missing_qc=options.ignore_missing_qc,
                         regenerate_reports=options.regenerate_reports)
        elif cmd == 'report':
            d.report(logging=options.logging,
                     summary=options.summary,
                     projects=options.projects,
                     full=options.full)
