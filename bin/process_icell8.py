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
from bcftbx.utils import mkdir
from bcftbx.utils import strip_ext
from bcftbx.IlluminaData import IlluminaData
from bcftbx.JobRunner import fetch_runner
from auto_process_ngs.applications import Command
from auto_process_ngs.simple_scheduler import SimpleScheduler
from auto_process_ngs.simple_scheduler import SchedulerReporter
from auto_process_ngs.fastq_utils import pair_fastqs

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
    p = argparse.ArgumentParser()
    p.add_argument("WELL_LIST",help="Well list file")
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
                   dest="aligner",default="bowtie2",
                   choices=["bowtie","bowtie2"],
                   help="aligner to use with fastq_screen (default: "
                   "'bowtie2')")
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
    args = p.parse_args()

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
    max_jobs = args.max_jobs

    # Report settings
    print "Unaligned dir     : %s" % args.unaligned_dir
    print "Output dir        : %s" % args.outdir
    print "Batch size (reads): %s" % args.batch_size
    print "Preferred genomes screen: %s" % args.preferred_conf
    print "Contaminants screen     : %s" % args.contaminants_conf
    print "Fastq_screen aligner    : %s" % args.aligner
    print "Fastq_screen threads    : %s" % args.threads
    print "Job runners:"
    for stage in stages:
        print "-- %20s: %s" % (stage,runners[stage])
    print "Maximum concurrent jobs : %s" % max_jobs

    # Get the input FASTQ file pairs
    illumina_data = IlluminaData(os.getcwd(),
                                 unaligned_dir=args.unaligned_dir)
    fastqs = []
    for project in illumina_data.projects:
        for sample in project.samples:
            for fq in sample.fastq:
                fastqs.append(os.path.join(sample.dirn,fq))

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
    log_dir = os.path.join(icell8_dir,"logs")
    stats_dir = os.path.join(icell8_dir,"stats")
    scripts_dir = os.path.join(icell8_dir,"scripts")
    for dirn in (icell8_dir,log_dir,stats_dir,scripts_dir):
        mkdir(dirn)

    # Initial stats
    stats_file = os.path.join(stats_dir,"icell8_stats.tsv")
    icell8_stats_cmd = Command('icell8_stats.py',
                               '-f',stats_file,
                               '-w',os.path.abspath(args.WELL_LIST))
    icell8_stats_cmd.add_args(*fastqs)
    icell8_stats = sched.submit(icell8_stats_cmd,
                                wd=icell8_dir,
                                name="initial_stats.%s" %
                                os.path.basename(fastqs[0]),
                                log_dir=log_dir)
    sched.callback("Initial statistics",
                   check_status,wait_for=(icell8_stats.name,))

    # Setup the filter and splitting job
    split_dir = os.path.join(icell8_dir,"_fastqs.filter_and_split")
    mkdir(split_dir)
    filter_and_split_cmd = Command('split_icell8_fastqs.py',
                                   '-w',os.path.abspath(args.WELL_LIST),
                                   '-o',split_dir,
                                   '-m','batch',
                                   '-s',args.batch_size)
    filter_and_split_cmd.add_args(*fastqs)
    # Submit the job
    filter_and_split = sched.submit(filter_and_split_cmd,
                                    wd=split_dir,
                                    name="filter_and_split.%s" %
                                    os.path.basename(fastqs[0]),
                                    log_dir=log_dir)
    # Wait for the filtering job to complete
    # (necessary as we don't know ahead of time what the
    # names of the batched files will be)
    sched.wait_for((filter_and_split.name,))
    if filter_and_split.exit_code != 0:
        logging.critical("Filter/split stage failed (exit code %s)"
                         % filter_and_split.exit_code)
        sys.exit(1)
    print "*** Quality filter and splitting stage completed ***"

    # Collect the post-filtering files for stats and read trimming
    filtered_fastqs = glob.glob(os.path.join(split_dir,"*.B*.r*.fastq"))

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
    icell8_stats = sched.submit(Command('/bin/bash',icell8_stats_script),
                                wd=icell8_dir,
                                name="post_quality_filter_stats.%s" %
                                os.path.basename(fastqs[0]),
                                log_dir=log_dir)
    sched.callback("Post-quality filter statistics",
                   check_status,wait_for=(icell8_stats.name,))
    
    # Set up the cutadapt jobs as a group
    fastq_pairs = pair_fastqs(filtered_fastqs)[0]
    trim_dir = os.path.join(icell8_dir,"_fastqs.trim_reads")
    mkdir(trim_dir)
    trim_reads = sched.group("trim_reads.%s" % os.path.basename(fastqs[0]))
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
    icell8_stats = sched.submit(Command('/bin/bash',icell8_stats_script),
                                wd=icell8_dir,
                                name="post_trimming_stats.%s" %
                                os.path.basename(fastqs[0]),
                                log_dir=log_dir)
    sched.callback("Post-trimming statistics",
                   check_status,wait_for=(icell8_stats.name,))

    # Set up the contaminant filter jobs as a group
    fastq_pairs = pair_fastqs(trimmed_fastqs)[0]
    contaminant_filter_dir = os.path.join(icell8_dir,
                                          "_fastqs.contaminant_filter")
    mkdir(contaminant_filter_dir)
    contaminant_filter = sched.group("contaminant_filter.%s" %
                                     os.path.basename(fastqs[0]))
    for pair in fastq_pairs:
        print "-- %s\n   %s" % (pair[0],pair[1])
        fqr1_in,fqr2_in = pair
        contaminant_filter_cmd = Command(
            'icell8_contamination_filter.py',
            '-p',os.path.abspath(args.preferred_conf),
            '-c',os.path.abspath(args.contaminants_conf),
            '-o',contaminant_filter_dir,
            '-a',args.aligner,
            '-n',args.threads,
            fqr1_in,fqr2_in)
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
    icell8_stats = sched.submit(Command('/bin/bash',icell8_stats_script),
                                wd=icell8_dir,
                                name="post_contaminant_filter_stats.%s" %
                                os.path.basename(fastqs[0]),
                                log_dir=log_dir)
    sched.callback("Post-contaminant filter statistics",
                   check_status,wait_for=(icell8_stats.name,))

    # Rebatch reads by barcode
    barcoded_fastqs_dir = os.path.join(icell8_dir,"fastqs")
    mkdir(barcoded_fastqs_dir)
    split_barcodes_cmd = Command('split_icell8_fastqs.py',
                                 '-w',os.path.abspath(args.WELL_LIST),
                                 '-o',barcoded_fastqs_dir,
                                 '-m','barcodes',
                                 '--no-filter')
    split_barcodes_cmd.add_args(*filtered_fastqs)
    split_barcodes_script = os.path.join(scripts_dir,
                                         "icell8_split_barcodes.sh")
    split_barcodes_cmd.make_wrapper_script(filen=split_barcodes_script,
                                           shell="/bin/bash")
    split_barcodes = sched.submit(Command('/bin/bash',split_barcodes_script),
                                  wd=icell8_dir,
                                  name="split_barcodes.%s" %
                                  os.path.basename(fastqs[0]),
                                  log_dir=log_dir)
    sched.wait()
    if split_barcodes.exit_code != 0:
        logging.critical("Rebatching reads by barcode failed (exit code %s)"
                         % split_barcodes.exit_code)
        sys.exit(1)
    print "*** Rebatching reads by barcodes completed ***"

    # Make hard links to the unassigned and failed barcode/quality
    # fastq files from the filter step
    extra_fastqs = glob.glob(
        os.path.join(split_dir,"*.unassigned.r*.fastq"))
    extra_fastqs.extend(glob.glob(
        os.path.join(split_dir,"*.failed_barcode.r*.fastq")))
    extra_fastqs.extend(glob.glob(
        os.path.join(split_dir,"*.failed_umi.r*.fastq")))
    for fastq in extra_fastqs:
        print "Linking: %s" % (os.path.basename(fastq))
        fq = os.path.join(barcoded_fastqs_dir,
                          os.path.basename(fastq))
        os.link(fastq,fq)
    print "*** Linking unassigned and failed barcode/UMIs completed ***"

    # Finish
    sched.stop()
