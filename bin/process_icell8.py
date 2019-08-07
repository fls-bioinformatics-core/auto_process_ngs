#!/usr/bin/env python
#
#     process_icell8.py: perform processing of Wafergen iCell8 data
#     Copyright (C) University of Manchester 2017-2019 Peter Briggs
#
"""
process_icell8.py

Utility to peform initial processing of data from Wafergen iCell8
platform.
"""

######################################################################
# Imports
######################################################################

import os
import sys
import logging
import argparse
import shutil
from bcftbx.utils import mkdir
from bcftbx.IlluminaData import IlluminaData
from bcftbx.IlluminaData import IlluminaDataError
from bcftbx.JobRunner import fetch_runner
from auto_process_ngs.analysis import AnalysisFastq
from auto_process_ngs.analysis import AnalysisProject
from auto_process_ngs.icell8.utils import ICell8WellList
from auto_process_ngs.icell8.pipeline import ICell8QCFilter
from auto_process_ngs.icell8.pipeline import ICell8FinalReporting
from auto_process_ngs.qc.pipeline import QCPipeline
import auto_process_ngs.envmod as envmod

# Fetch configuration settings
import auto_process_ngs.settings
__settings = auto_process_ngs.settings.Settings()

# Module specific logger
logger = logging.getLogger(__name__)

######################################################################
# Main
######################################################################

if __name__ == "__main__":
    # Pipeline stages
    stages = ('default','contaminant_filter','qc','statistics','report')
    # Fetch defaults
    default_batch_size = __settings.icell8.batch_size
    default_aligner = __settings.icell8.aligner
    default_mammalian_conf = __settings.icell8.mammalian_conf_file
    default_contaminants_conf = __settings.icell8.contaminants_conf_file
    default_fastq_screen_subset = __settings.qc.fastq_screen_subset
    default_nprocessors = 1
    # Handle the command line
    p = argparse.ArgumentParser(
        description="Perform initial QC on FASTQs from Wafergen "
        "ICell8: assign to barcodes, filter on barcode & UMI quality, "
        "trim reads, perform contaminant filtering and split by "
        "barcode.")
    p.add_argument("well_list",metavar="WELL_LIST",help="Well list file")
    p.add_argument("fastqs",nargs='*',metavar="FASTQ_R1 FASTQ_R2",
                   help="FASTQ file pairs")
    p.add_argument("-u","--unaligned",
                   dest="unaligned_dir",default=None,
                   help="process FASTQs from 'unaligned' dir with output "
                   "from bcl2fastq (NB cannot be used with -p option)")
    p.add_argument("-p","--project",metavar="NAME",
                   dest="project",default=None,
                   help="process FASTQS from project directory NAME (NB "
                   "if -o not specified then this will also be used as "
                   "the output directory; cannot be used with -u option)")
    p.add_argument("-o","--outdir",
                   dest="outdir",default=None,
                   help="directory to write outputs to "
                   "(default: 'CWD/icell8', or project dir if -p "
                   "is specified)")
    p.add_argument("-m","--mammalian",
                   dest="mammalian_conf",
                   default=default_mammalian_conf,
                   help="fastq_screen 'conf' file with the "
                   "'mammalian' genome indices (default: %s)"
                   % default_mammalian_conf)
    p.add_argument("-c","--contaminants",
                   dest="contaminants_conf",
                   default=default_contaminants_conf,
                   help="fastq_screen 'conf' file with the "
                   "'contaminant' genome indices (default: "
                   "%s)" % default_contaminants_conf)
    p.add_argument("-q","--quality-filter",action='store_true',
                   dest="quality_filter",
                   help="filter out read pairs with low quality "
                   "barcode and UMI sequences (not recommended for "
                   "NextSeq data)")
    p.add_argument("-a","--aligner",
                   dest="aligner",default=None,
                   choices=["bowtie","bowtie2"],
                   help="aligner to use with fastq_screen (default: "
                   "don't specify the aligner)")
    p.add_argument("-r","--runner",metavar="STAGE=RUNNER",
                   action="append",dest="runners",default=list(),
                   help="explicitly specify runner definitions for "
                   "running pipeline jobs at each stage. STAGE "
                   "can be one of %s. If STAGE is not specified "
                   "then it is assumed to be 'default'. RUNNER "
                   "must be a valid job runner specification e.g. "
                   "'GEJobRunner(-j y)'. Multiple --runner arguments "
                   "can be specified (default: '%s')" %
                   (','.join(["'%s'" % s for s in stages]),
                    __settings.general.default_runner))
    p.add_argument("-n","--nprocessors",metavar='STAGE=N',
                   action="append",dest="nprocessors",default=list(),
                   help="specify number of processors to use at each "
                   "stage. STAGE can be one of %s. If STAGE is not "
                   "specified then it is assumed to be 'default'. "
                   "Multiple --nprocessors arguments can be "
                   "specified (default: %d)" %
                   (','.join(["'%s'" % s for s in stages]),
                    default_nprocessors))
    p.add_argument("-s","--size",type=int,
                   dest="batch_size",default=default_batch_size,
                   help="number of reads per batch when splitting "
                   "FASTQ files for processing (default: %s)" %
                   default_batch_size)
    p.add_argument("-j","--max-jobs",type=int,
                   dest="max_jobs",
                   default= __settings.general.max_concurrent_jobs,
                   help="maxiumum number of concurrent jobs to run "
                   "(default: %d)"
                   % __settings.general.max_concurrent_jobs)
    p.add_argument('--modulefiles',action='store',
                   dest='modulefiles',default=None,
                   help="comma-separated list of environment "
                   "modules to load before executing commands "
                   "(overrides any modules specified in the global "
                   "settings)")
    p.add_argument("--no-contaminant-filter",action='store_true',
                   dest='no_contaminant_filter',
                   help="don't perform contaminant filter step "
                   "(default is to do contaminant filtering)")
    p.add_argument("--no-cleanup",action='store_true',
                   dest="no_cleanup",
                   help="don't remove intermediate Fastq files "
                   "(default is to delete intermediate Fastqs once "
                   "no longer needed)")
    p.add_argument('--force',action='store_true',
                   dest='force',default=False,
                   help="force overwrite of existing outputs")
    p.add_argument("-v","--verbose",action='store_true',
                   dest="verbose",default=False,
                   help="produce verbose output for diagnostics")
    p.add_argument("--no-quality-filter",action='store_true',
                   dest="no_quality_filter",
                   help="deprecated: kept for backwards compatibility "
                   "only as barcode/UMI quality checks are now "
                   "disabled by default")
    p.add_argument("--threads",type=int,
                   dest="threads",default=None,
                   help="deprecated (use -n/--nprocessors option "
                   "instead): number of threads to use with multicore "
                   "tasks (e.g. 'contaminant_filter')")
    args = p.parse_args()

    # Deal with environment modules
    modulefiles = args.modulefiles
    if modulefiles is None:
        try:
            modulefiles = __settings.modulefiles['process_icell8']
        except KeyError:
            # No environment modules specified
            pass
    if modulefiles is not None:
        for modulefile in modulefiles.split(','):
            envmod.load(modulefile)

    # Deal with job runners
    runners = dict()
    for runner in args.runners:
        try:
            stage,runner_spec = runner.split('=')
        except ValueError: # too few values to unpack
            stage = 'default'
            runner_spec = runner
        if stage not in stages:
            logger.fatal("Bad stage for --runner option: %s" % stage)
            sys.exit(1)
        runners[stage] = fetch_runner(runner_spec)
    try:
        default_runner = runners['default']
    except KeyError:
        default_runner = __settings.runners.icell8
    for stage in stages:
        if stage not in runners:
            if stage == 'qc':
                # Use the general QC settings
                stage_runner = __settings.runners['qc']
            else:
                # Look for Icell8-specific runner
                try:
                    stage_runner = __settings.runners['icell8_%s' % stage]
                except KeyError:
                    stage_runner = default_runner
            runners[stage] = stage_runner

    # Deal with number of processors
    nprocessors = dict()
    for nprocs in args.nprocessors:
        try:
            stage,n = nprocs.split('=')
        except ValueError: # too few values to unpack
            stage = 'default'
            n = nprocs
        if stage not in stages:
            logger.fatal("Bad stage for --nprocessors option: %s" % stage)
            sys.exit(1)
        nprocessors[stage] = int(n)
    if args.threads is not None:
        for stage in ('contaminant_filter','statistics'):
            if stage not in nprocessors:
                logging.warning("Setting nprocessors for stage '%s' "
                                "from --threads option" % stage)
                nprocessors[stage] = args.threads
    try:
        default_nprocessors = nprocessors['default']
    except KeyError:
        default_nprocessors = 1
    for stage in stages:
        stage_nprocessors = default_nprocessors
        if stage not in nprocessors:
            if stage == 'qc':
                # Use the general QC settings
                try:
                    stage_nprocessors = __settings.qc.nprocessors
                except (KeyError,AttributeError):
                    pass
            else:
                # Look for Icell8-specific nprocessors
                try:
                    stage_nprocessors = __settings.icell8['nprocessors_%s'
                                                          % stage]
                except KeyError:
                    pass
            nprocessors[stage] = stage_nprocessors

    # Check for clashing -u/-p
    if args.project and args.unaligned_dir:
        logger.fatal("Cannot specify -u and -p together")
        sys.exit(1)

    # Check for contaminant filtering inputs
    if args.mammalian_conf is None or \
       not os.path.isfile(args.mammalian_conf):
        logging.fatal("Mammalian genome panel not specified "
                      "or doesn't exist (-m)")
        sys.exit(1)
    if args.contaminants_conf is None or \
       not os.path.isfile(args.contaminants_conf):
        logging.fatal("Contaminant genome panel not specified "
                      "or doesn't exist (-c)")
        sys.exit(1)

    # Output dir
    if args.outdir is None:
        if args.project:
            outdir = args.project
        else:
            outdir = "icell8"
    else:
        outdir = args.outdir

    # Other settings
    well_list = os.path.abspath(args.well_list)
    max_jobs = args.max_jobs
    do_quality_filter = args.quality_filter
    do_contaminant_filter = (not args.no_contaminant_filter)
    do_clean_up = (not args.no_cleanup)

    # Report settings
    print("Unaligned dir     : %s" % args.unaligned_dir)
    print("Project           : %s" % args.project)
    print("Well list file    : %s" % well_list)
    print("Output dir        : %s" % outdir)
    print("Batch size (reads): %s" % args.batch_size)
    print("Quality filter barcodes/UMIs: %s" %
          ('yes' if do_quality_filter else 'no'))
    print("Filter contaminants: %s" %
          ('yes' if do_contaminant_filter else 'no'))
    if do_contaminant_filter:
        print("Mammalian genome panel  : %s" % args.mammalian_conf)
        with open(args.mammalian_conf) as fp:
            for line in fp:
                if line.startswith("DATABASE"):
                    print("-- %s" % line.split('\t')[1])
        print("Contaminant genome panel: %s" % args.contaminants_conf)
        with open(args.contaminants_conf) as fp:
            for line in fp:
                if line.startswith("DATABASE"):
                    print("-- %s" % line.split('\t')[1])
            print("Fastq_screen aligner    : %s" % args.aligner)
    print("Maximum concurrent jobs : %s" % max_jobs)
    print("Stage specific settings :")
    for stage in stages:
        print("-- %s: %s (nprocs=%d)" % (stage,
                                         runners[stage],
                                         nprocessors[stage]))
    if modulefiles is not None:
        print("Environment modules:")
        for modulefile in modulefiles.split(','):
            print("-- %s" % modulefile)
    print("Clean-up intermediate Fastqs: %s" %
          ('yes' if do_clean_up else 'no'))

    # Check well list file
    try:
        ICell8WellList(well_list).barcodes()
    except Exception as ex:
        logger.fatal("Couldn't load data from well list file '%s'"
                     % well_list)
        sys.exit(1)

    # Get the input FASTQ file pairs
    fastqs = []
    # Collect files from command line
    for fq in args.fastqs:
        fastqs.append(os.path.abspath(fq))
    # Collect files from unaligned dir
    if fastqs and args.unaligned_dir is not None:
        logger.warning("Ignoring unaligned dir '%s'" %
                       args.unaligned_dir)
    elif args.unaligned_dir:
        try:
            illumina_data = IlluminaData(
                os.getcwd(),
                unaligned_dir=args.unaligned_dir)
            for project in illumina_data.projects:
                for sample in project.samples:
                    for fq in sample.fastq:
                        fastqs.append(os.path.join(sample.dirn,fq))
        except IlluminaDataError:
            logger.fatal("Couldn't find FASTQS in directory '%s'" %
                          args.unaligned_dir)
    # Collect files from project
    analysis_project = None
    if fastqs and args.project is not None:
        logger.warning("Ignoring project '%s'" % args.project)
    elif args.project:
        analysis_project = AnalysisProject(args.project,
                                           args.project,
                                           fastq_dir='fastqs')
        for sample in analysis_project.samples:
            for fq in sample.fastq:
                fastqs.append(os.path.join(
                    analysis_project.fastq_dir,
                    fq))
    if not fastqs:
        logger.fatal("No FASTQs found")
        sys.exit(1)

    # Make top-level output dirs
    icell8_dir = os.path.abspath(outdir)
    if os.path.exists(icell8_dir) and args.project is None:
        if not args.force:
            logger.fatal("Output destination '%s': already exists "
                         "(remove or use --force to overwrite)" %
                         icell8_dir)
            sys.exit(1)
        logger.warning("Removing existing output destination '%s'" %
                       icell8_dir)
        shutil.rmtree(icell8_dir)
    log_dir = os.path.join(icell8_dir,"logs")
    scripts_dir = os.path.join(icell8_dir,"scripts")
    for dirn in (icell8_dir,log_dir,scripts_dir):
        mkdir(dirn)

    # Copy well list file into output directory
    shutil.copy(well_list,outdir)
    well_list = os.path.join(outdir,os.path.basename(well_list))
    if analysis_project is not None:
        analysis_project.info['icell8_well_list'] = os.path.basename(well_list)
        analysis_project.info.save()

    # Set up pipelines
    pipelines = []

    # ICELL QC and filtering
    print("Setting up a pipeline for ICELL processing")
    pipelines.append(
        ICell8QCFilter(outdir,fastqs,well_list,
                       args.mammalian_conf,
                       args.contaminants_conf,
                       args.batch_size,
                       aligner=args.aligner,
                       do_contaminant_filter=do_contaminant_filter,
                       do_quality_filter=do_quality_filter,
                       do_clean_up=do_clean_up,
                       nprocessors=nprocessors))

    # Final reporting
    print("Setting up a pipeline for final reporting")
    pipelines.append(
        ICell8FinalReporting(outdir,
                             project=analysis_project))

    # Chain the pipelines
    print("Merging the pipelines")
    ppl = pipelines[0]
    for p in pipelines[1:]:
        ppl.append_pipeline(p)

    # Execute the pipelines
    print("Running the final pipeline")
    exit_status = ppl.run(log_dir=log_dir,scripts_dir=scripts_dir,
                          default_runner=default_runner,
                          runners=runners,
                          max_jobs=max_jobs,
                          verbose=args.verbose)
    if exit_status != 0:
        # Finished with error
        logger.critical("Pipeline failed: exit status %s" % exit_status)
        sys.exit(exit_status)

    # Run the QC
    print("Running the QC")
    runqc = QCPipeline()
    runqc.add_project(analysis_project,
                      qc_protocol="singlecell",
                      fastq_dir="fastqs.samples",
                      qc_dir="qc.samples",
                      multiqc=True)
    runqc.add_project(analysis_project,
                      qc_protocol="singlecell",
                      fastq_dir="fastqs.barcodes",
                      qc_dir="qc.barcodes",
                      multiqc=False)
    exit_status = runqc.run(max_jobs=max_jobs,
                            batch_size=25,
                            runners={
                                'qc_runner': runners['qc'],
                                'report_runner': default_runner,
                                'verify_runner': default_runner
                            },
                            fastq_strand_indexes=
                            __settings.fastq_strand_indexes,
                            nthreads=nprocessors['qc'],
                            default_runner=default_runner,
                            verbose=args.verbose)
    if exit_status != 0:
        # Finished with error
        logger.critical("QC failed: exit status %s" % exit_status)
        sys.exit(exit_status)

    # Finish
    print("All pipelines completed ok")
    sys.exit(0)
