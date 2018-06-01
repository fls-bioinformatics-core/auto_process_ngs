#!/usr/bin/env python
#
#     process_icell8.py: perform processing of Wafergen iCell8 data
#     Copyright (C) University of Manchester 2017 Peter Briggs
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
from auto_process_ngs.pipeliner import Pipeline
from auto_process_ngs.analysis import AnalysisFastq
from auto_process_ngs.analysis import AnalysisProject
from auto_process_ngs.icell8.utils import ICell8WellList
from auto_process_ngs.icell8.pipeline import *
from auto_process_ngs.qc.runqc import RunQC
from auto_process_ngs.qc.illumina_qc import IlluminaQC
from auto_process_ngs.qc.fastq_strand import build_fastq_strand_conf
import auto_process_ngs.envmod as envmod

# Fetch configuration settings
import auto_process_ngs.settings
__settings = auto_process_ngs.settings.Settings()

# Module specific logger
logger = logging.getLogger(__name__)

######################################################################
# Magic numbers
######################################################################

DEFAULT_BATCH_SIZE = __settings.icell8.batch_size

######################################################################
# Classes
######################################################################

# No classes defined

######################################################################
# Functions
######################################################################

# No functions defined

######################################################################
# Main
######################################################################

if __name__ == "__main__":
    # Pipeline stages
    stages = ('default','contaminant_filter','qc','statistics')
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
    print "Unaligned dir     : %s" % args.unaligned_dir
    print "Project           : %s" % args.project
    print "Well list file    : %s" % well_list
    print "Output dir        : %s" % outdir
    print "Batch size (reads): %s" % args.batch_size
    print "Quality filter barcodes/UMIs: %s" % \
        ('yes' if do_quality_filter else 'no')
    print "Filter contaminants: %s" % \
        ('yes' if do_contaminant_filter else 'no')
    if do_contaminant_filter:
        print "Mammalian genome panel  : %s" % args.mammalian_conf
        with open(args.mammalian_conf) as fp:
            for line in fp:
                if line.startswith("DATABASE"):
                    print "-- %s" % line.split('\t')[1]
        print "Contaminant genome panel: %s" % args.contaminants_conf
        with open(args.contaminants_conf) as fp:
            for line in fp:
                if line.startswith("DATABASE"):
                    print "-- %s" % line.split('\t')[1]
            print "Fastq_screen aligner    : %s" % args.aligner
    print "Maximum concurrent jobs : %s" % max_jobs
    print "Stage specific settings :"
    for stage in stages:
        print "-- %s: %s (nprocs=%d)" % (stage,
                                  runners[stage],
                                  nprocessors[stage])
    if modulefiles is not None:
        print "Environment modules:"
        for modulefile in modulefiles.split(','):
            print "-- %s" % modulefile
    print "Clean-up intermediate Fastqs: %s" % \
        ('yes' if do_clean_up else 'no')

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

    # Basename for output fastqs and job names etc
    basename = AnalysisFastq(fastqs[0]).sample_name

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
    stats_dir = os.path.join(icell8_dir,"stats")
    scripts_dir = os.path.join(icell8_dir,"scripts")
    for dirn in (icell8_dir,log_dir,stats_dir,scripts_dir):
        mkdir(dirn)

    # Copy well list file into output directory
    shutil.copy(well_list,outdir)
    well_list = os.path.join(outdir,os.path.basename(well_list))
    if analysis_project is not None:
        analysis_project.info['icell8_well_list'] = os.path.basename(well_list)
        analysis_project.info.save()

    # Final Fastq directories
    barcode_fastqs_dir = os.path.join(icell8_dir,"fastqs.barcodes")
    sample_fastqs_dir = os.path.join(icell8_dir,"fastqs.samples")

    # Set up pipelines
    pipelines = []

    # ICell8-specific QC and filtering tasks
    # Only run these stages if the final fastqs don't exist
    if not (os.path.exists(barcode_fastqs_dir) and os.path.exists(sample_fastqs_dir)):

        # Dedicated pipeline for these tasks
        print "Setting pipeline for ICell8 QC filter"
        ppl = Pipeline(name="ICell8: QC filter")

        # Initial stats
        initial_stats = GetICell8Stats("Initial statistics",
                                       fastqs,
                                       os.path.join(stats_dir,"icell8_stats.tsv"),
                                       well_list,
                                       unassigned=True,
                                       nprocs=nprocessors['statistics'])
        ppl.add_task(initial_stats,runner=runners['statistics'])

        # Split fastqs into batches
        batch_dir = os.path.join(icell8_dir,"_fastqs.batched")
        batch_fastqs = SplitFastqsIntoBatches("Batch Fastqs",fastqs,
                                              batch_dir,basename,
                                              batch_size=args.batch_size)
        ppl.add_task(batch_fastqs)
        collect_batch_fastqs = CollectFiles("Collect batched files",
                                            batch_dir,
                                            batch_fastqs.output().pattern)
        ppl.add_task(collect_batch_fastqs,requires=(batch_fastqs,))

        # Setup the filtering jobs as a group
        filter_dir = os.path.join(icell8_dir,"_fastqs.quality_filter")
        filter_fastqs = FilterICell8Fastqs("Filter Fastqs",
                                           collect_batch_fastqs.output(),
                                           filter_dir,
                                           well_list=well_list,
                                           mode='none',
                                           discard_unknown_barcodes=True,
                                           quality_filter=do_quality_filter)
        ppl.add_task(filter_fastqs,requires=(collect_batch_fastqs,))
        # Collect the files from the filtering jobs
        collect_filtered_fastqs = CollectFiles(
            "Collect filtered fastqs",
            filter_dir,
            filter_fastqs.output().patterns.assigned)
        collect_filtered_unassigned = CollectFiles(
            "Collect filtered fastqs (unassigned barcodes)",
            filter_dir,
            filter_fastqs.output().patterns.unassigned)
        collect_filtered_failed_barcodes = CollectFiles(
            "Collect filtered fastqs (failed barcodes)",
            filter_dir,
            filter_fastqs.output().patterns.failed_barcodes)
        collect_filtered_failed_umis = CollectFiles(
            "Collect filtered fastqs (failed UMIs)",
            filter_dir,
            filter_fastqs.output().patterns.failed_umis)
        for task in (collect_filtered_fastqs,
                     collect_filtered_unassigned,
                     collect_filtered_failed_barcodes,
                     collect_filtered_failed_umis):
            ppl.add_task(task,requires=(filter_fastqs,))
    
        # Post filtering stats
        filter_stats = GetICell8Stats("Post-filtering statistics",
                                      collect_filtered_fastqs.output(),
                                      initial_stats.output(),
                                      suffix="_filtered",
                                      append=True,
                                      nprocs=nprocessors['statistics'])
        ppl.add_task(filter_stats,requires=(initial_stats,
                                            collect_filtered_fastqs),
                     runner=runners['statistics'])

        # Use cutadapt to find reads with poly-G regions
        poly_g_dir = os.path.join(icell8_dir,"_fastqs.poly_g")
        get_poly_g_reads = GetReadsWithPolyGRegions(
            "Find reads with poly-G regions",
            collect_filtered_fastqs.output(),
            poly_g_dir)
        ppl.add_task(get_poly_g_reads,requires=(collect_filtered_fastqs,))
        collect_poly_g_fastqs = CollectFiles("Collect poly-G fastqs",
                                             poly_g_dir,
                                             get_poly_g_reads.output().pattern)
        ppl.add_task(collect_poly_g_fastqs,requires=(get_poly_g_reads,))
        poly_g_stats = GetICell8PolyGStats("Poly-G region statistics",
                                           collect_poly_g_fastqs.output(),
                                           initial_stats.output(),
                                           suffix="_poly_g",
                                           append=True,
                                           nprocs=nprocessors['statistics'])
        ppl.add_task(poly_g_stats,
                     requires=(collect_poly_g_fastqs,filter_stats),
                     runner=runners['statistics'])

        # Set up the cutadapt jobs as a group
        trim_dir = os.path.join(icell8_dir,"_fastqs.trim_reads")
        trim_reads = TrimReads("Read trimming",
                               collect_filtered_fastqs.output(),
                               trim_dir)
        ppl.add_task(trim_reads,requires=(collect_filtered_fastqs,))
        collect_trimmed_fastqs = CollectFiles("Collect trimmed fastqs",
                                              trim_dir,
                                              trim_reads.output().pattern)
        ppl.add_task(collect_trimmed_fastqs,requires=(trim_reads,))

        # Post read trimming stats
        trim_stats = GetICell8Stats("Post-trimming statistics",
                                    collect_trimmed_fastqs.output(),
                                    initial_stats.output(),
                                    suffix="_trimmed",
                                    append=True,
                                    nprocs=nprocessors['statistics'])
        ppl.add_task(trim_stats,requires=(collect_trimmed_fastqs,
                                          poly_g_stats),
                     runner=runners['statistics'])

        # Set up the contaminant filter jobs as a group
        if do_contaminant_filter:
            contaminant_filter_dir = os.path.join(
                icell8_dir,
                "_fastqs.contaminant_filter")
            contaminant_filter = FilterContaminatedReads(
                "Contaminant filtering",
                collect_trimmed_fastqs.output(),
                contaminant_filter_dir,
                args.mammalian_conf,
                args.contaminants_conf,
                aligner=args.aligner,
                threads=nprocessors['contaminant_filter'])
            ppl.add_task(contaminant_filter,
                         requires=(collect_trimmed_fastqs,),
                         runner=runners['contaminant_filter'])
            collect_contaminant_filtered = CollectFiles(
                "Collect contaminant-filtered fastqs",
                contaminant_filter_dir,
                contaminant_filter.output().pattern)
            ppl.add_task(collect_contaminant_filtered,
                         requires=(contaminant_filter,))

            # Post contaminant filter stats
            final_stats = GetICell8Stats(
                "Post-contaminant filter statistics",
                collect_contaminant_filtered.output(),
                initial_stats.output(),
                suffix="_contaminant_filtered",
                append=True,
                nprocs=nprocessors['statistics'])
            ppl.add_task(final_stats,
                         requires=(collect_contaminant_filtered,
                                   trim_stats),
                         runner=runners['statistics'])
            fastqs_in = collect_contaminant_filtered.output()
            split_barcodes_requires = (collect_contaminant_filtered,)
        else:
            fastqs_in = collect_trimmed_fastqs.output()
            split_barcodes_requires = (collect_trimmed_fastqs,)

        # Prepare for rebatching reads by barcode and sample by splitting
        # each batch by barcode
        split_barcoded_fastqs_dir = os.path.join(icell8_dir,"_fastqs.split_barcodes")
        split_barcodes = SplitByBarcodes("Split batches by barcode",
                                         fastqs_in,
                                         split_barcoded_fastqs_dir)
        ppl.add_task(split_barcodes,requires=split_barcodes_requires)
        collect_split_barcodes = CollectFiles("Collect barcode-split fastqs",
                                              split_barcoded_fastqs_dir,
                                              split_barcodes.output().pattern)
        ppl.add_task(collect_split_barcodes,requires=(split_barcodes,))
        # Merge (concat) fastqs into single pairs per barcode
        barcode_fastqs = MergeBarcodeFastqs(
            "Assemble reads by barcode",
            collect_split_barcodes.output(),
            collect_filtered_unassigned.output(),
            collect_filtered_failed_barcodes.output(),
            collect_filtered_failed_umis.output(),
            barcode_fastqs_dir,
            basename)
        ppl.add_task(barcode_fastqs,requires=(collect_split_barcodes,))
        collect_barcode_fastqs = CollectFiles(
            "Collect final barcoded fastqs",
            barcode_fastqs_dir,
            barcode_fastqs.output().patterns.assigned)
        ppl.add_task(collect_barcode_fastqs,requires=(barcode_fastqs,))
        # Merge (concat) fastqs into single pairs per barcode
        sample_fastqs_dir = os.path.join(icell8_dir,"fastqs.samples")
        sample_fastqs = MergeSampleFastqs("Assemble reads by sample",
                                          collect_split_barcodes.output(),
                                          well_list,
                                          sample_fastqs_dir)
        ppl.add_task(sample_fastqs,requires=(collect_split_barcodes,))

        # Final stats for verification
        final_barcode_stats = GetICell8Stats(
            "Post-barcode splitting and merging statistics",
            collect_barcode_fastqs.output(),
            initial_stats.output(),
            suffix="_final",
            append=True,
            nprocs=nprocessors['statistics'])
        if do_contaminant_filter:
            final_barcode_stats_requires = (collect_barcode_fastqs,
                                            final_stats,)
        else:
            final_barcode_stats_requires = (collect_barcode_fastqs,
                                            trim_stats,)
        ppl.add_task(final_barcode_stats,
                     requires=final_barcode_stats_requires,
                     runner=runners['statistics'])

        # Verify that barcodes are okay
        check_barcodes = CheckICell8Barcodes(
            "Verify barcodes are consistent",
            collect_barcode_fastqs.output())
        ppl.add_task(check_barcodes,requires=(collect_barcode_fastqs,))

        # Generate XLSX version of stats
        xlsx_stats = ConvertStatsToXLSX(
            "Convert statistics to XLSX",
            final_barcode_stats.output(),
            os.path.join(stats_dir,"icell8_stats.xlsx"))
        ppl.add_task(xlsx_stats,requires=(final_barcode_stats,))

        # Cleanup outputs
        cleanup_tasks = []
        cleanup_tasks.append(CleanupDirectory("Remove batched Fastqs",
                                              batch_dir))
        cleanup_tasks.append(CleanupDirectory("Remove filtered Fastqs",
                                              filter_dir))
        cleanup_tasks.append(CleanupDirectory("Remove poly-G region stats data",
                                              poly_g_dir))
        cleanup_tasks.append(CleanupDirectory("Remove trimmed Fastqs",
                                              trim_dir))
        cleanup_tasks.append(CleanupDirectory("remove barcode split Fastqs",
                                              split_barcoded_fastqs_dir))
        if do_contaminant_filter:
            cleanup_tasks.append(CleanupDirectory("Remove contaminant "
                                                  "filtered Fastqs",
                                                  contaminant_filter_dir))
            
        if do_clean_up:
            # Wait until all stages are finished before doing clean up
            cleanup_requirements = [collect_barcode_fastqs,
                                    sample_fastqs]
            if do_contaminant_filter:
                cleanup_requirements.append(final_stats)
            else:
                cleanup_requirements.append(trim_stats)
            for task in cleanup_tasks:
                ppl.add_task(task,requires=cleanup_requirements)

        # Add to list of pipelines
        pipelines.append(ppl)

    # Final reporting
    print "Setting up a pipeline for final reporting"
    ppl = Pipeline(name="ICell8: final reporting")
    # Reset primary fastq dir (if working in a project)
    if analysis_project is not None:
        update_project_data = UpdateProjectData(
            "Updating metadata associated with the project",
            icell8_dir,"fastqs.samples")
        ppl.add_task(update_project_data)
    # Final report
    final_report = ReportProcessing("Generate processing report",
                                    outdir)
    ppl.add_task(final_report)
    pipelines.append(ppl)

    # Execute the pipelines
    print "Running the pipelines"
    for ppl in pipelines:
        exit_status = ppl.run(log_dir=log_dir,scripts_dir=scripts_dir,
                              default_runner=runners['default'],
                              max_jobs=max_jobs)
        if exit_status != 0:
            # Finished with error
            logger.critical("Pipeline failed: exit status %s" % exit_status)
            sys.exit(exit_status)

    # Run the QC
    print "Running the QC"
    # Set up conf file for strandedness determination
    fastq_strand_conf = None
    if analysis_project.info.organism:
        print "Organisms: %s" % analysis_project.info.organism
        fastq_strand_indexes = build_fastq_strand_conf(
            analysis_project.info.organism.lower().split(','),
            __settings.fastq_strand_indexes)
        if fastq_strand_indexes:
            print "Setting up conf file for strandedness determination"
            fastq_strand_conf = os.path.join(analysis_project.dirn,
                                             "fastq_strand.conf")
            with open(fastq_strand_conf,'w') as fp:
                fp.write("%s\n" % fastq_strand_indexes)
        else:
            print "No matching indexes for strandedness determination"
    else:
        print "No organisms specified"
    # Set up the QC
    illumina_qc = IlluminaQC(
        nthreads=nprocessors['qc'],
        fastq_screen_subset=default_fastq_screen_subset,
        fastq_strand_conf=fastq_strand_conf)
    runqc = RunQC()
    runqc.add_project(analysis_project,
                      fastq_dir="fastqs.samples",
                      qc_dir="qc.samples",
                      illumina_qc=illumina_qc)
    runqc.add_project(analysis_project,
                      fastq_dir="fastqs.barcodes",
                      qc_dir="qc.barcodes",
                      illumina_qc=illumina_qc)
    exit_status = runqc.run(multiqc=True,
                            qc_runner=runners['qc'],
                            report_runner=default_runner,
                            max_jobs=max_jobs,
                            batch_size=25)
    if exit_status != 0:
        # Finished with error
        logger.critical("QC failed: exit status %s" % exit_status)
        sys.exit(exit_status)

    # Finish
    print "All pipelines completed ok"
    sys.exit(0)
