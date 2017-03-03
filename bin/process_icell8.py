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

######################################################################
# Magic numbers
######################################################################

DEFAULT_BATCH_SIZE = 50000000

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
    p.add_argument("-r","--runner",
                   dest="runner",default="SimpleJobRunner",
                   help="explicitly specify a runner definition for "
                   "running pipeline jobs (e.g. 'GEJobRunner(-j y)')")
    p.add_argument("-s","--size",type=int,
                   dest="batch_size",default=DEFAULT_BATCH_SIZE,
                   help="number of reads per batch when splitting "
                   "FASTQ files for processing (default: %s)" %
                   DEFAULT_BATCH_SIZE)
    args = p.parse_args()

    # Get the input FASTQ file pairs
    illumina_data = IlluminaData(os.getcwd(),
                                 unaligned_dir=args.unaligned_dir)
    fastqs = []
    for project in illumina_data.projects:
        for sample in project.samples:
            for fq in sample.fastq:
                fastqs.append(os.path.join(sample.dirn,fq))

    # Set up a scheduler for running jobs
    icell8_qc_runner = runner = fetch_runner(args.runner)
    max_jobs = 4
    sched_reporter = SchedulerReporter(
        job_start="Started  #%(job_number)d: %(job_name)s:\n-- %(command)s",
        job_end=  "Finished #%(job_number)d: %(job_name)s"
    )
    sched = SimpleScheduler(runner=icell8_qc_runner,
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
    sched.wait()
    if icell8_stats.exit_code != 0:
        logging.critical("Initial stats stage failed (exit code %d)"
                         % icell8_stats.exit_code)
        sys.exit(1)
    print "*** Initial statistics stage completed ***"

    # Setup the filter and splitting job
    split_dir = os.path.join(icell8_dir,"filter_and_split")
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
    # Wait for the job to complete
    # (necessary as we don't know ahead of time what the
    # names of the batched files will be)
    sched.wait()
    if filter_and_split.exit_code != 0:
        logging.critical("Filter/split stage failed (exit code %d)"
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
    sched.wait()
    if icell8_stats.exit_code != 0:
        logging.critical("Post-quality filter stats stage failed (exit code %d)"
                         % icell8_stats.exit_code)
        sys.exit(1)
    print "*** Post-quality filter statistics stage completed ***"
    
    # Set up the cutadapt jobs as a group
    fastq_pairs = pair_fastqs(filtered_fastqs)[0]
    trim_reads = sched.group("trim_reads.%s" % os.path.basename(fastqs[0]))
    for pair in fastq_pairs:
        print "-- %s\n   %s" % (pair[0],pair[1])
        fqr1_in,fqr2_in = pair
        trim_dir = os.path.join(icell8_dir,"trim_reads")
        mkdir(trim_dir)
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
    sched.wait()
    exit_code = max([j.exit_code for j in trim_reads.jobs])
    if exit_code != 0:
        logging.critical("Read trimming stage failed (exit code %d)"
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
    sched.wait()
    if icell8_stats.exit_code != 0:
        logging.critical("Post-trimming stats stage failed (exit code %d)"
                         % icell8_stats.exit_code)
        sys.exit(1)
    print "*** Post-trimming statistics stage completed ***"

    # Set up the contaminant filter jobs as a group
    contaminant_filter = sched.group("contaminant_filter.%s" %
                                     os.path.basename(fastqs[0]))
    fastq_pairs = pair_fastqs(trimmed_fastqs)[0]
    for pair in fastq_pairs:
        print "-- %s\n   %s" % (pair[0],pair[1])
        fqr1_in,fqr2_in = pair
        contaminant_filter_dir = os.path.join(icell8_dir,
                                              "contaminant_filter")
        mkdir(contaminant_filter_dir)
        contaminant_filter_cmd = Command(
            'icell8_contamination_filter.py',
            '-p',os.path.abspath(args.preferred_conf),
            '-c',os.path.abspath(args.contaminants_conf),
            '-o',contaminant_filter_dir,
            '-a',args.aligner,
            fqr1_in,fqr2_in)
        # Submit the job
        job = contaminant_filter.add(
            contaminant_filter_cmd,
            wd=contaminant_filter_dir,
            name="contaminant_filter.%s" %
            os.path.basename(fqr1_in),
            log_dir=log_dir)
    contaminant_filter.close()
    sched.wait()
    exit_code = max([j.exit_code for j in contaminant_filter.jobs])
    if exit_code != 0:
        logging.critical("Contaminant filtering stage failed (exit code %d)"
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
    sched.wait()
    if icell8_stats.exit_code != 0:
        logging.critical("Post-contaminant filter stats stage failed (exit code %d)"
                         % icell8_stats.exit_code)
        sys.exit(1)
    print "*** Post-contaminant filter statistics stage completed ***"

    # Finish
    sched.stop()
