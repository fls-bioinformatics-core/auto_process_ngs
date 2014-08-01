#!/bin/env python
#
#     auto_process.py: automated processing of Illumina sequence data
#     Copyright (C) University of Manchester 2013-14 Peter Briggs
#
#########################################################################
#
# auto_process.py
#
#########################################################################

"""
First attempt at an automated data processing & QC pipeline in Python

Implements a program for automating stages of a standard protocol for
processing and QC'ing Illumina sequencing data.

The stages are:

    setup
    config
    make_fastqs
    setup_analysis_dirs
    run_qc
    archive
    publish_qc

The 'setup' stage creates an analysis directory and acquires the basic
data about the sequencing run from a source directory. Subsequent stages
should be run in sequence to create fastq files, set up analysis
directories for each project, and run QC scripts for each sample in
each project.

Additional commands are available:

    clone
    merge_fastq_dirs
    update_fastq_stats

but these are not part of the standard workflow - they are used for
special cases and testing.

"""

__version__ = "0.0.77"

#######################################################################
# Imports
#######################################################################

import sys
import os
import subprocess
import logging
import optparse
import bcf_utils
import auto_process_ngs.auto_process_settings
from auto_process_ngs.auto_processor import AutoProcess

#######################################################################
# Functions
#######################################################################

def make_custom_sample_sheet(input_sample_sheet,output_sample_sheet=None):
    # Read sample sheet info from input_sample_sheet
    # Do clean up
    # Write to output_sample_sheet (if specified)
    # Return CasavaSampleSheet object
    sample_sheet = IlluminaData.get_casava_sample_sheet(input_sample_sheet)
    for line in sample_sheet:
        if not line['SampleProject']:
            line['SampleProject'] = line['SampleID']
    sample_sheet.fix_illegal_names()
    sample_sheet.fix_duplicated_names()
    if output_sample_sheet is not None:
        sample_sheet.write(output_sample_sheet)
    return sample_sheet

def get_bases_mask(run_info_xml,sample_sheet_file):
    # Return bases mask string generated from data in RunInfo.xml and
    # sample sheet files
    # Get initial bases mask
    bases_mask = IlluminaData.IlluminaRunInfo(run_info_xml).bases_mask
    print "Bases mask: %s (from RunInfo.xml)" % bases_mask
    # Update bases mask from sample sheet
    example_barcode = IlluminaData.get_casava_sample_sheet(sample_sheet_file)[0]['Index']
    bases_mask = IlluminaData.fix_bases_mask(bases_mask,example_barcode)
    print "Bases mask: %s (updated for barcode sequence '%s')" % (bases_mask,
                                                                  example_barcode)
    return bases_mask

def list_available_commands(cmds):
    # Pretty-print available commands
        print ""
        print "Available commands are:"
        for cmd in cmds:
            print "\t%s" % cmd
        print ""

def set_debug(flag):
    # Turn on debug output
    if flag:
        logging.getLogger().setLevel(logging.DEBUG)

# Command line parsers

def setup_parser():
    p = optparse.OptionParser(usage="%prog setup [OPTIONS] DIR",
                              version="%prog "+__version__,
                              description="Automatically process Illumina sequencing "
                              "data from DIR.")
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
    return p

def clone_parser():
    p = optparse.OptionParser(usage="%prog clone [OPTIONS] DIR CLONE_DIR",
                              version="%prog "+__version__,
                              description="Make a copy of an existing auto_processed analysis "
                              "directory DIR, in a new directory CLONE_DIR. The clone will "
                              "not include any project directories, but will copy the "
                              "projects.info file.")
    p.add_option('--copy-fastqs',action='store_true',dest='copy_fastqs',default=False,
                 help="Copy fastq.gz files from DIR into DIR2 (default is to make a "
                 "link to the bcl-to-fastq directory)")
    return p

def config_parser():
    p  = optparse.OptionParser(usage="%prog config [OPTIONS] [ANALYSIS_DIR]",
                               version="%prog "+__version__,
                               description="Query and configure automatic processing "
                               "parameters and settings for ANALYSIS_DIR.")
    p.add_option('--set',action='append',dest='key_value',default=None,
                 help="Set the value of a parameter. KEY_VALUE should be of the form "
                 "'<param>=<value>'. Multiple --set options can be specified.")
    # Deprecated options
    deprecated = optparse.OptionGroup(p,'Deprecated/defunct options')
    deprecated.add_option('--show',action='store_true',dest='show',default=False,
                          help="Show the values of parameters and settings (does "
                          "nothing; use 'config' with no options to display settings)")
    p.add_option_group(deprecated)
    return p

def make_fastqs_parser():
    p = optparse.OptionParser(usage="%prog make_fastqs [OPTIONS] [ANALYSIS_DIR]",
                              version="%prog "+__version__,
                              description="Generate fastq files from raw bcl files "
                              "produced by Illumina sequencer within ANALYSIS_DIR.")
    # General options
    add_no_save_option(p)
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
    nprocessors = auto_process_ngs.auto_process_settings.bcl2fastq.nprocessors
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
    bcl_to_fastq.add_option('--nprocessors',action='store',
                            dest='nprocessors',default=nprocessors,
                            help="explicitly specify number of processors to use for "
                            "bclToFastq (default %s, change in settings file)" % nprocessors)
    bcl_to_fastq.add_option('--ignore-missing-bcl',action='store_true',
                            dest='ignore_missing_bcl',default=False,
                            help="use the --ignore-missing-bcl option for bcl2fastq (treat "
                            "missing bcl files as no call)")
    bcl_to_fastq.add_option('--ignore-missing-stats',action='store_true',
                            dest='ignore_missing_stats',default=False,
                            help="use the --ignore-missing-stats option for bcl2fastq (fill "
                            "in with zeroes when *.stats files are missing)")
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
    statistics.add_option('--report-barcodes',action='store_true',
                          dest='report_barcodes',default=False,
                          help="analyse and report barcode indices for all lanes after "
                          "generating fastq files")
    statistics.add_option('--barcodes-file',action='store',
                          dest='barcodes_file',default=None,
                          help="specify output file for barcode analysis report")
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
    p.add_option_group(deprecated)
    return p

def merge_fastq_dirs_parser():
    p = optparse.OptionParser(usage="%prog merge_fastq_dirs [OPTIONS] [ANALYSIS_DIR]",
                              version="%prog "+__version__,
                              description="Automatically merge fastq directories froms "
                              "multiple bcl-to-fastq runs within ANALYSIS_DIR. Use this "
                              "command if 'make_fastqs' step was run multiple times to "
                              "process subsets of lanes.")
    nprocessors = auto_process_ngs.auto_process_settings.bcl2fastq.nprocessors
    p.add_option('--primary-unaligned-dir',action='store',
                 dest='unaligned_dir',default='bcl2fastq',
                 help="merge fastqs from additional bcl-to-fastq directories into "
                 "UNALIGNED_DIR. Original data will be moved out of the way first. "
                 "Defaults to 'bcl2fastq'.")
    add_dry_run_option(p)
    return p

def setup_analysis_dirs_parser():
    p = optparse.OptionParser(usage="%prog setup_analysis_dirs [OPTIONS] [ANALYSIS_DIR]",
                              version="%prog "+__version__,
                              description="Create analysis subdirectories for projects "
                              "defined in projects.info file in ANALYSIS_DIR.")
    p.add_option('--ignore-missing-metadata',action='store_true',
                 dest='ignore_missing_metadata',default=False,
                 help="force creation of project directories even if metadata is not "
                 "set (default is to fail if metadata is missing)")
    return p

def run_qc_parser():
    p = optparse.OptionParser(usage="%prog run_qc [OPTIONS] [ANALYSIS_DIR]",
                              version="%prog "+__version__,
                              description="Automatically process Illumina sequence from "
                              "ANALYSIS_DIR.")
    max_concurrent_jobs = auto_process_ngs.auto_process_settings.general.max_concurrent_jobs
    p.add_option('--projects',action='store',
                 dest='project_pattern',default=None,
                 help="simple wildcard-based pattern specifying a subset of projects "
                 "and samples to run the QC on. PROJECT_PATTERN should be of the form "
                 "'pname[/sname]', where 'pname' specifies a project (or set of "
                 "projects) and 'sname' optionally specifies a sample (or set of "
                 "samples).")
    p.add_option('--max-jobs',action='store',
                 dest='max_jobs',default=max_concurrent_jobs,type='int',
                 help="explicitly specify maximum number of concurrent QC jobs to run "
                 "(default %s, change in settings file)" % max_concurrent_jobs)
    p.add_option('--no-ungzip-fastqs',action='store_true',dest='no_ungzip_fastqs',
                 help="don't create uncompressed copies of fastq.gz files")
    return p

def publish_qc_parser():
    p = optparse.OptionParser(usage="%prog publish_qc [OPTIONS] [ANALYSIS_DIR]",
                              version="%prog "+__version__,
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
    return p

def archive_parser():
    p = optparse.OptionParser(usage="%prog archive [OPTIONS] [ANALYSIS_DIR]",
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
    default_group = auto_process_ngs.auto_process_settings.archive.group
    p.add_option('--group',action='store',dest='group',default=default_group,
                 help="specify the name of group for the archived files NB only works "
                 "when the archive is a local directory (default: %s)" % default_group)
    default_chmod = auto_process_ngs.auto_process_settings.archive.chmod
    p.add_option('--chmod',action='store',dest='chmod',default=default_chmod,
                 help="specify chmod operations for the archived files (default: "
                 "%s)" % default_chmod)
    add_dry_run_option(p)
    return p

def report_parser():
    p  = optparse.OptionParser(usage="%prog report [OPTIONS] [ANALYSIS_DIR]",
                              version="%prog "+__version__,
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
    return p

def generic_parser(cmd,description=None):
    if description is None:
        description = "Automatically process Illumina sequence from ANALYSIS_DIR."
    p  = optparse.OptionParser(usage="%prog "+cmd+" [OPTIONS] [ANALYSIS_DIR]",
                              version="%prog "+__version__,
                              description=description)
    return p

def add_no_save_option(p):
    # Add --no-save option to a parser
    p.add_option('--no-save',action='store_true',dest='no_save',default=False,
                 help="Don't save parameter changes to the auto_process.info file")
    return p

def add_dry_run_option(p):
    # Add --dry-run option to a parser
    p.add_option('--dry-run',action='store_true',dest='dry_run',default=False,
                 help="Dry run i.e. report but don't perform any actions")
    return p

def add_debug_option(p):
    # Add debug option to a parser
    p.add_option('--debug',action='store_true',dest='debug',default=False,
                 help="Turn on debugging output from Python libraries")
    return p

#######################################################################
# Main program
#######################################################################

if __name__ == "__main__":

    # Available commands and corresponding
    cmd_parsers = bcf_utils.OrderedDictionary()
    cmd_parsers['setup'] = setup_parser()
    cmd_parsers['clone'] = clone_parser()
    cmd_parsers['config'] = config_parser()
    cmd_parsers['make_fastqs'] = make_fastqs_parser()
    cmd_parsers['merge_fastq_dirs'] = merge_fastq_dirs_parser()
    cmd_parsers['update_fastq_stats'] = generic_parser('update_fastq_stats',
                                                       "(Re)generate statistics for fastq "
                                                       "files produced from 'make_fastqs'.")
    cmd_parsers['setup_analysis_dirs'] = setup_analysis_dirs_parser()
    cmd_parsers['run_qc'] = run_qc_parser()
    cmd_parsers['archive'] = archive_parser()
    cmd_parsers['publish_qc'] = publish_qc_parser()
    cmd_parsers['report'] = report_parser()

    # Process major command
    try:
        cmd = sys.argv[1]
    except IndexError:
        cmd = "help"
    if cmd == "help" or cmd == "--help" or cmd == "-h":
        list_available_commands(cmd_parsers)
        sys.exit(0)
    else:
        err = None
        if cmd not in cmd_parsers:
            err = "Unrecognised command '%s'" % cmd
        elif len(sys.argv) < 2:
            err = "Need to supply a command"
        if err is not None:
            print err
            list_available_commands(cmd_parsers)
            sys.stderr.write("%s\n" % err)
            sys.exit(1)
        p = cmd_parsers[cmd]

    # Add debug options (available for all commands)
    add_debug_option(p)
    
    # Process remaining command line arguments
    options,args = p.parse_args(sys.argv[2:])

    # Report name and version
    print "%s version %s" % (os.path.basename(sys.argv[0]),__version__)

    # Turn on debugging?
    set_debug(options.debug)

    # Allow saving of parameters?
    try:
        allow_save = not options.no_save
    except AttributeError:
        allow_save = True

    # Setup the processing object and run the requested command
    if cmd == 'setup':
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
            d.make_fastqs(skip_rsync=options.skip_rsync,
                          nprocessors=options.nprocessors,
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
            d.generate_stats()
        elif cmd == 'setup_analysis_dirs':
            d.setup_analysis_dirs(ignore_missing_metadata=
                                  options.ignore_missing_metadata)
        elif cmd == 'run_qc':
            d.run_qc(projects=options.project_pattern,
                     max_jobs=options.max_jobs,
                     no_ungzip_fastqs=options.no_ungzip_fastqs)
        elif cmd == 'config':
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
            d.copy_to_archive(archive_dir=options.archive_dir,
                              platform=options.platform,
                              year=options.year,
                              dry_run=options.dry_run,
                              group=options.group,
                              chmod=options.chmod)
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
