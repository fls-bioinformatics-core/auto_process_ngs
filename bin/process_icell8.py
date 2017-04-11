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
import glob
import shutil
from bcftbx.utils import mkdir
from bcftbx.utils import strip_ext
from bcftbx.IlluminaData import IlluminaData
from bcftbx.IlluminaData import IlluminaDataError
from bcftbx.JobRunner import fetch_runner
from auto_process_ngs.applications import Command
from auto_process_ngs.simple_scheduler import SimpleScheduler
from auto_process_ngs.simple_scheduler import SchedulerReporter
from auto_process_ngs.fastq_utils import pair_fastqs
from auto_process_ngs.utils import AnalysisFastq
import auto_process_ngs.envmod as envmod

# Fetch configuration settings
import auto_process_ngs.settings
__settings = auto_process_ngs.settings.Settings()

######################################################################
# Magic numbers
######################################################################

DEFAULT_BATCH_SIZE = 5000000

######################################################################
# Functions
######################################################################

def check_status(name,jobs,sched):
    """
    Check the exit status of a set of jobs

    This function is not called directly, instead it should be
    passed as part of a callback to check the exit status of a
    job or group of jobs from the scheduler.

    For example:

    >>> sched.callback("My job",check_status,wait_for(('my_job',))

    If any of the jobs have a non-zero exit status then
    the program is terminated.

    Arguments:
      name (str): name for the callback
      jobs (list): list of SchedulerJob instances
      sched (SimpleScheduler): scheduler instance

    """
    print "*** %s completed ***" % name
    for job in jobs:
        exit_code = job.exit_code
        if exit_code != 0:
            logging.critical("Job '%s' failed: exit code %s"
                             % (job.name,exit_code))
            sched.stop()
            sys.exit(1)
    print "%s: ok" % name

######################################################################
# Main
######################################################################

if __name__ == "__main__":
    # Handle the command line
    p = argparse.ArgumentParser(
        description="Perform initial QC on FASTQs from Wafergen "
        "ICell8: assign to barcodes, filter on barcode & UMI quality, "
        "trim reads, perform contaminant filtering and split by "
        "barcode.")
    p.add_argument("well_list",metavar="WELL_LIST",help="Well list file")
    p.add_argument("fastqs",nargs='*',metavar="FASTQ_R1 FASTQ_R2",
                   help="FASTQ file pairs")
    p.add_argument("--unaligned",
                   dest="unaligned_dir",default="bcl2fastq",
                   help="'unaligned' dir with output from "
                   "bcl2fastq")
    p.add_argument("-o","--outdir",
                   dest="outdir",default="icell8",
                   help="directory to write outputs to "
                   "(default: 'CWD/icell8')")
    p.add_argument("-p","--preferred",
                   dest="preferred_conf",
                   help="fastq_screen 'conf' file with the "
                   "'preferred' genome indices")
    p.add_argument("-c","--contaminants",
                   dest="contaminants_conf",
                   help="fastq_screen 'conf' file with the "
                   "'contaminant' genome indices")
    p.add_argument("-a","--aligner",
                   dest="aligner",default=None,
                   choices=["bowtie","bowtie2"],
                   help="aligner to use with fastq_screen (default: "
                   "don't specify the aligner)")
    p.add_argument("-n","--threads",type=int,
                   dest="threads",default=1,
                   help="number of threads to use with fastq_screen "
                   "(default: 1)")
    p.add_argument("-r","--runner",metavar="STAGE=RUNNER",
                   action="append",dest="runners",default=list(),
                   help="explicitly specify runner definitions for "
                   "running pipeline jobs at each stage. STAGE "
                   "can be one of 'default','contaminant_filter'. "
                   "RUNNER must be a valid job runner specification "
                   "e.g. 'GEJobRunner(-j y)'. Multiple --runner "
                   "arguments can be specified (default: "
                   "SimpleJobRunner)")
    p.add_argument("-s","--size",type=int,
                   dest="batch_size",default=DEFAULT_BATCH_SIZE,
                   help="number of reads per batch when splitting "
                   "FASTQ files for processing (default: %s)" %
                   DEFAULT_BATCH_SIZE)
    p.add_argument("-m","--max-jobs",type=int,
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
    p.add_argument('--force',action='store_true',
                   dest='force',default=False,
                   help="force overwrite of existing outputs")
    args = p.parse_args()

    # Deal with module files
    if args.modulefiles is not None:
        modulefiles = args.modulefiles.split(',')
        for modulefile in modulefiles:
            envmod.load(modulefile)

    # Deal with job runners
    stages = ('default','contaminant_filter')
    runners = dict()
    for runner in args.runners:
        try:
            stage,runner_spec = runner.split('=')
        except ValueError: # too few values to unpack
            stage = 'default'
            runner_spec = runner
        if stage not in stages:
            logging.fatal("Bad stage for --runner option: %s" % stage)
            sys.exit(1)
        runners[stage] = fetch_runner(runner_spec)
    try:
        default_runner = runners['default']
    except KeyError:
        default_runner = fetch_runner('SimpleJobRunner')
    for stage in stages:
        if stage not in runners:
            runners[stage] = default_runner

    # Other settings
    well_list = os.path.abspath(args.well_list)
    max_jobs = args.max_jobs

    # Report settings
    print "Unaligned dir     : %s" % args.unaligned_dir
    print "Well list file    : %s" % well_list
    print "Output dir        : %s" % args.outdir
    print "Batch size (reads): %s" % args.batch_size
    print "Preferred genomes screen: %s" % args.preferred_conf
    with open(args.preferred_conf) as fp:
        for line in fp:
            if line.startswith("DATABASE"):
                print "-- %s" % line.split('\t')[1]
    print "Contaminants screen     : %s" % args.contaminants_conf
    with open(args.contaminants_conf) as fp:
        for line in fp:
            if line.startswith("DATABASE"):
                print "-- %s" % line.split('\t')[1]
    print "Fastq_screen aligner    : %s" % args.aligner
    print "Fastq_screen threads    : %s" % args.threads
    print "Maximum concurrent jobs : %s" % max_jobs
    print "Job runners:"
    for stage in stages:
        print "-- %s: %s" % (stage,runners[stage])
    if args.modulefiles is not None:
        print "Environment modules:"
        for modulefile in modulefiles:
            print "-- %s" % modulefile

    # Get the input FASTQ file pairs
    fastqs = []
    try:
        illumina_data = IlluminaData(os.getcwd(),
                                     unaligned_dir=args.unaligned_dir)
        for project in illumina_data.projects:
            for sample in project.samples:
                for fq in sample.fastq:
                    fastqs.append(os.path.join(sample.dirn,fq))
    except IlluminaDataError:
        logging.warning("Couldn't find FASTQS in directory '%s'" %
                        args.unaligned_dir)
    for fq in args.fastqs:
        fastqs.append(os.path.abspath(fq))
    if not fastqs:
        logging.fatal("No FASTQs found")
        sys.exit(1)

    # Basename for output fastqs and job names etc
    basename = AnalysisFastq(fastqs[0]).sample_name

    # Set up a scheduler for running jobs
    sched_reporter = SchedulerReporter(
        job_start="Started  #%(job_number)d: %(job_name)s:\n-- %(command)s",
        job_end=  "Finished #%(job_number)d: %(job_name)s"
    )
    sched = SimpleScheduler(runner=runners['default'],
                            max_concurrent=max_jobs,
                            reporter=sched_reporter)
    sched.start()

    # Make top-level output dirs
    icell8_dir = os.path.abspath(args.outdir)
    if os.path.exists(icell8_dir):
        if not args.force:
            logging.fatal("Output destination '%s': already exists "
                          "(remove or use --force to overwrite)" %
                          icell8_dir)
            sys.exit(1)
        logging.warning("Removing existing output destination '%s'" %
                        icell8_dir)
        shutil.rmtree(icell8_dir)
    log_dir = os.path.join(icell8_dir,"logs")
    stats_dir = os.path.join(icell8_dir,"stats")
    scripts_dir = os.path.join(icell8_dir,"scripts")
    for dirn in (icell8_dir,log_dir,stats_dir,scripts_dir):
        mkdir(dirn)

    # Initial stats
    stats_file = os.path.join(stats_dir,"icell8_stats.tsv")
    icell8_stats_cmd = Command('icell8_stats.py',
                               '-f',stats_file,
                               '-w',well_list)
    icell8_stats_cmd.add_args(*fastqs)
    initial_stats = sched.submit(icell8_stats_cmd,
                                 wd=icell8_dir,
                                 name="initial_stats.%s" % basename,
                                 log_dir=log_dir)
    sched.callback("Initial statistics",
                   check_status,wait_for=(initial_stats.name,))

    # Split fastqs into batches
    batch_dir = os.path.join(icell8_dir,"_fastqs.batched")
    mkdir(batch_dir)
    batch_fastqs_cmd = Command('batch_fastqs.py',
                               '-s',args.batch_size,
                               '-o',batch_dir,
                               '-b',basename)
    batch_fastqs_cmd.add_args(*fastqs)
    batch_fastqs = sched.submit(batch_fastqs_cmd,
                                wd=icell8_dir,
                                name="batch_fastqs.%s" % basename,
                                log_dir=log_dir)
    # Wait for the batching job to complete
    sched.wait_for((batch_fastqs.name,))
    if batch_fastqs.exit_code != 0:
        logging.critical("Fastq batching stage failed (exit code %s)"
                         % batch_fastqs.exit_code)
        sys.exit(1)
    print "*** Fastq batching stage completed ***"

    # Collect the batched files for processing
    batched_fastqs = glob.glob(os.path.join(batch_dir,"*.B*.r*.fastq"))

    # Setup the quality filter jobs as a group
    fastq_pairs = pair_fastqs(batched_fastqs)[0]
    filter_dir = os.path.join(icell8_dir,"_fastqs.quality_filter")
    mkdir(filter_dir)
    quality_filter = sched.group("quality_filter.%s" % basename)
    for pair in fastq_pairs:
        print "-- %s\n   %s" % (pair[0],pair[1])
        batch_basename = os.path.basename(pair[0])[:-len(".r1.fastq")]
        print "   %s" % batch_basename
        filter_cmd = Command('split_icell8_fastqs.py',
                             '-w',well_list,
                             '-b',batch_basename,
                             '-o',filter_dir,
                             '-m','none',
                             '--filter')
        filter_cmd.add_args(*pair)
        # Submit the job
        job = quality_filter.add(filter_cmd,
                                 wd=filter_dir,
                                 name="quality_filter.%s" %
                                 batch_basename,
                                 log_dir=log_dir)
    quality_filter.close()
    sched.wait_for((quality_filter.name,))
    exit_code = max([j.exit_code for j in quality_filter.jobs])
    if exit_code != 0:
        logging.critical("Quality filtering stage failed (exit code %s)"
                         % exit_code)
        sys.exit(1)
    print "*** Quality filtering stage completed ***"

    # Collect the post-filtering files for stats and read trimming
    filtered_fastqs = glob.glob(os.path.join(filter_dir,
                                             "*.B*.filtered.r*.fastq"))

    # Post quality filter stats
    icell8_stats_cmd = Command('icell8_stats.py',
                               '-f',stats_file,
                               '--suffix','_quality_filtered',
                               '--append',)
    icell8_stats_cmd.add_args(*filtered_fastqs)
    icell8_stats_script = os.path.join(scripts_dir,
                                       "icell8_stats_quality_filtered.sh")
    icell8_stats_cmd.make_wrapper_script(filen=icell8_stats_script,
                                         shell="/bin/bash")
    filter_stats = sched.submit(Command('/bin/bash',icell8_stats_script),
                                wd=icell8_dir,
                                name="post_quality_filter_stats.%s" %
                                basename,
                                log_dir=log_dir,
                                wait_for=(initial_stats.name,))
    sched.callback("Post-quality filter statistics",
                   check_status,wait_for=(filter_stats.name,))
    
    # Set up the cutadapt jobs as a group
    fastq_pairs = pair_fastqs(filtered_fastqs)[0]
    trim_dir = os.path.join(icell8_dir,"_fastqs.trim_reads")
    mkdir(trim_dir)
    trim_reads = sched.group("trim_reads.%s" % basename)
    for pair in fastq_pairs:
        print "-- %s\n   %s" % (pair[0],pair[1])
        fqr1_in,fqr2_in = pair
        cutadapt_cmd = Command(
            'cutadapt',
            '-a','AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT',
            '-a','AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT',
            '-a','AAAAAAAA',
            '-a','TTTTTTTT',
            '-m',20,
            '--trim-n',
            '--max-n',0.7,
            '-q',25)
        # Output pair (nb R2 then R1)
        fqr1_out = strip_ext(os.path.basename(fqr1_in),'.fastq') \
                   + '.trimmed.fastq'
        fqr2_out = strip_ext(os.path.basename(fqr2_in),'.fastq') \
                   + '.trimmed.fastq'
        cutadapt_cmd.add_args('-o',fqr2_out,
                              '-p',fqr1_out)
        # Input pair (nb R2 then R1)
        cutadapt_cmd.add_args(fqr2_in,fqr1_in)
        # Submit the job
        job = trim_reads.add(cutadapt_cmd,
                             wd=trim_dir,
                             name="cutadapt.%s" % os.path.basename(fqr1_in),
                             log_dir=log_dir)
    trim_reads.close()
    sched.wait_for((trim_reads.name,))
    exit_code = max([j.exit_code for j in trim_reads.jobs])
    if exit_code != 0:
        logging.critical("Read trimming stage failed (exit code %s)"
                         % exit_code)
        sys.exit(1)
    print "*** Read trimming stage completed ***"

    # Collect the files for stats and contamination filter step
    trimmed_fastqs = glob.glob(os.path.join(trim_dir,"*.trimmed.fastq"))

    # Post read trimming stats
    icell8_stats_cmd = Command('icell8_stats.py',
                               '-f',stats_file,
                               '--suffix','_trimmed',
                               '--append',)
    icell8_stats_cmd.add_args(*trimmed_fastqs)
    icell8_stats_script = os.path.join(scripts_dir,
                                       "icell8_stats_trimmed.sh")
    icell8_stats_cmd.make_wrapper_script(filen=icell8_stats_script,
                                         shell="/bin/bash")
    trim_stats = sched.submit(Command('/bin/bash',icell8_stats_script),
                              wd=icell8_dir,
                              name="post_trimming_stats.%s" %
                              basename,
                              log_dir=log_dir,
                              wait_for=(filter_stats.name,))
    sched.callback("Post-trimming statistics",
                   check_status,wait_for=(trim_stats.name,))

    # Set up the contaminant filter jobs as a group
    fastq_pairs = pair_fastqs(trimmed_fastqs)[0]
    contaminant_filter_dir = os.path.join(icell8_dir,
                                          "_fastqs.contaminant_filter")
    mkdir(contaminant_filter_dir)
    contaminant_filter = sched.group("contaminant_filter.%s" % basename)
    for pair in fastq_pairs:
        print "-- %s\n   %s" % (pair[0],pair[1])
        fqr1_in,fqr2_in = pair
        contaminant_filter_cmd = Command(
            'icell8_contamination_filter.py',
            '-p',os.path.abspath(args.preferred_conf),
            '-c',os.path.abspath(args.contaminants_conf),
            '-o',contaminant_filter_dir,
            '-n',args.threads)
        if args.aligner is not None:
            contaminant_filter_cmd.add_args('-a',args.aligner)
        contaminant_filter_cmd.add_args(fqr1_in,fqr2_in)
        # Submit the job
        job = contaminant_filter.add(
            contaminant_filter_cmd,
            wd=contaminant_filter_dir,
            name="contaminant_filter.%s" %
            os.path.basename(fqr1_in),
            log_dir=log_dir,
            runner=runners['contaminant_filter'])
    contaminant_filter.close()
    sched.wait_for((contaminant_filter.name,))
    exit_code = max([j.exit_code for j in contaminant_filter.jobs])
    if exit_code != 0:
        logging.critical("Contaminant filtering stage failed (exit code %s)"
                         % exit_code)
        sys.exit(exit_code)
    print "*** Contaminant filtering stage completed ***"

    # Post contaminant filter stats
    filtered_fastqs = glob.glob(os.path.join(contaminant_filter_dir,
                                             "*.trimmed.filtered.fastq"))
    icell8_stats_cmd = Command('icell8_stats.py',
                               '-f',stats_file,
                               '--suffix','_contaminant_filtered',
                               '--append',)
    icell8_stats_cmd.add_args(*filtered_fastqs)
    icell8_stats_script = os.path.join(scripts_dir,
                                       "icell8_stats_contaminant_filtered.sh")
    icell8_stats_cmd.make_wrapper_script(filen=icell8_stats_script,
                                         shell="/bin/bash")
    contaminant_filter_stats = sched.submit(
        Command('/bin/bash',icell8_stats_script),
        wd=icell8_dir,
        name="post_contaminant_filter_stats.%s" %
        basename,
        log_dir=log_dir,
        wait_for=(trim_stats.name,))
    sched.callback("Post-contaminant filter statistics",
                   check_status,wait_for=(contaminant_filter_stats.name,))

    # Rebatch reads by barcode
    # First: split each batch by barcode
    fastq_pairs = pair_fastqs(filtered_fastqs)[0]
    barcoded_fastqs_dir = os.path.join(icell8_dir,"_fastqs.barcodes")
    mkdir(barcoded_fastqs_dir)
    split_barcodes = sched.group("split_barcodes.%s" % basename)
    for pair in fastq_pairs:
        print "-- %s\n   %s" % (pair[0],pair[1])
        batch_basename = os.path.basename(pair[0])[:-len(".r1.fastq")+1]
        print "   %s" % batch_basename
        split_cmd = Command('split_icell8_fastqs.py',
                            '-o',barcoded_fastqs_dir,
                            '-b',batch_basename,
                            '-m','barcodes')
        split_cmd.add_args(*pair)
        # Submit the job
        job = split_barcodes.add(split_cmd,
                                 wd=barcoded_fastqs_dir,
                                 name="split_barcodes.%s" %
                                 batch_basename,
                                 log_dir=log_dir)
    split_barcodes.close()
    sched.wait_for((split_barcodes.name,))
    exit_code = max([j.exit_code for j in split_barcodes.jobs])
    if exit_code != 0:
        logging.critical("Splitting by barcodes stage failed (exit code %s)"
                         % exit_code)
        sys.exit(1)
    print "*** Splitting by barcodes stage completed ***"
    # Second: merge across batches
    barcoded_fastqs = glob.glob(os.path.join(barcoded_fastqs_dir,
                                             "*.B*.*.fastq"))
    # Get barcodes
    barcodes = set()
    for fq in barcoded_fastqs:
        barcode = os.path.basename(fq).split('.')[-3]
        barcodes.add(barcode)
    barcodes = sorted(list(barcodes))
    # Group files by barcode
    fastq_groups = dict()
    for barcode in barcodes:
        fastqs = filter(lambda fq: (fq.endswith("%s.r1.fastq" % barcode) or
                                    fq.endswith("%s.r2.fastq" % barcode)),
                        barcoded_fastqs)
        fastq_groups[barcode] = fastqs
    # Group barcodes into subgroups
    GROUP_LEN = 25
    barcode_groups = [barcodes[i:i+GROUP_LEN]
                      for i in xrange(0,len(barcodes),GROUP_LEN)]
    # Merge (concat) fastqs into single pairs per barcode
    final_fastqs_dir = os.path.join(icell8_dir,"fastqs")
    mkdir(final_fastqs_dir)
    merge_fastqs = sched.group("merge_fastqs.%s" % basename)
    for i,barcode_group in enumerate(barcode_groups):
        group_name = "barcodes%06d" % i
        print "Barcode group: %s" % group_name
        merge_cmd = Command('split_icell8_fastqs.py',
                            '-o',final_fastqs_dir,
                            '-b',basename,
                            '-m','barcodes')
        for barcode in barcode_group:
            print "-- %s" % barcode
            merge_cmd.add_args(*fastq_groups[barcode])
        merge_fastqs_script = os.path.join(scripts_dir,
                                           "merge_fastqs.%s.sh" %
                                           group_name)
        merge_cmd.make_wrapper_script(filen=merge_fastqs_script,
                                      shell="/bin/bash")
        job = merge_fastqs.add(Command('/bin/bash',merge_fastqs_script),
                               wd=icell8_dir,
                               name="merge_fastqs.%s.%s" % (basename,
                                                            group_name),
                               log_dir=log_dir)
    # Deal with unassigned and failed quality reads
    for name in ("unassigned","failed_barcode","failed_umi"):
        fastqs = glob.glob(os.path.join(filter_dir,"*.%s.r*.fastq" % name))
        if not fastqs:
            continue
        merge_cmd = Command('batch_fastqs.py',
                            '-o',final_fastqs_dir,
                            '-b','%s.%s' % (basename,name),
                            '-m')
        merge_cmd.add_args(*fastqs)
        merge_fastqs_script = os.path.join(scripts_dir,
                                           "merge_fastqs.%s.sh" %
                                           name)
        merge_cmd.make_wrapper_script(filen=merge_fastqs_script,
                                      shell="/bin/bash")
        job = merge_fastqs.add(Command('/bin/bash',merge_fastqs_script),
                               wd=icell8_dir,
                               name="merge_fastqs.%s.%s" % (basename,
                                                            name),
                               log_dir=log_dir)
    # Close group and wait for merging barcodes to complete
    merge_fastqs.close()
    sched.wait_for((merge_fastqs.name,))
    exit_code = max([j.exit_code for j in merge_fastqs.jobs])
    if exit_code != 0:
        logging.critical("Merging barcoded FASTQs failed (exit code %s)"
                         % exit_code)
        sys.exit(exit_code)
    print "*** Merging barcoded FASTQs stage completed ***"

    # Finish
    sched.wait()
    print "All jobs completed"
    sched.stop()
