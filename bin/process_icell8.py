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
from bcftbx.IlluminaData import IlluminaData
from bcftbx.JobRunner import SimpleJobRunner
from auto_process_ngs.applications import Command
from auto_process_ngs.simple_scheduler import SimpleScheduler
from auto_process_ngs.fastq_utils import pair_fastqs

######################################################################
# Main
######################################################################

if __name__ == "__main__":
    # Handle the command line
    p = argparse.ArgumentParser()
    #p.add_argument("FQ_R1",help="R1 FASTQ file")
    #p.add_argument("FQ_R2",help="Matching R2 FASTQ file")
    p.add_argument("WELL_LIST",help="Well list file")
    p.add_argument("--unaligned",
                   dest="unaligned_dir",default="bcl2fastq",
                   help="'unaligned' dir with output from "
                   "bcl2fastq")
    p.add_argument("-p","--preferred",
                   dest="preferred_conf",
                   help="fastq_screen 'conf' file with the "
                   "'preferred' genome indices")
    p.add_argument("-c","--contaminants",
                   dest="contaminants_conf",
                   help="fastq_screen 'conf' file with the "
                   "'contaminant' genome indices")
    p.add_argument("-a","--aligner",
                   dest="aligner",default="barcodes",
                   choices=["bowtie","bowtie2"],
                   help="aligner to use with fastq_screen (default: "
                   "'bowtie2')")
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
    icell8_qc_runner = SimpleJobRunner()
    max_jobs = 4
    sched = SimpleScheduler(runner=icell8_qc_runner,
                            max_concurrent=max_jobs)
    sched.start()

    # Make top-level output dirs
    log_dir = os.path.join(os.getcwd(),"logs")
    icell8_dir = os.path.join(os.getcwd(),"icell8")
    stats_dir = os.path.join(os.getcwd(),"stats")
    for dirn in (log_dir,icell8_dir,stats_dir):
        mkdir(dirn)

    # Setup the filter and splitting job
    split_dir = os.path.join(icell8_dir,"filter_and_split")
    mkdir(split_dir)
    stats_file = os.path.join(stats_dir,"stats.filter.tsv")
    filter_and_split_cmd = Command('split_icell8_fastqs.py',
                                   '-w',os.path.abspath(args.WELL_LIST),
                                   '-o',split_dir,
                                   '-m','batch',
                                   '-s',5000,
                                   '-f',stats_file)
    filter_and_split_cmd.add_args(*fastqs)
    # Submit the job
    print "Running %s" % filter_and_split_cmd
    filter_and_split = sched.submit(filter_and_split_cmd,
                                    wd=split_dir,
                                    name="filter_and_split.%s" %
                                    os.path.basename(fastqs[0]),
                                    log_dir=log_dir)
    print "Job: %s" % filter_and_split
    # Wait for the job to complete
    # (necessary as we don't know ahead of time what the
    # names of the batched files will be)
    filter_and_split.wait()
    if filter_and_split.exit_code != 0:
        logging.critical("Filter/split stage failed (exit code %d)"
                         % filter_and_split.exit_code)
        sys.exit(1)
    print "*** Quality filter and splitting stage completed ***"

    # Collect the files for the read trimming step
    fastq_pairs = pair_fastqs(
        glob.glob(os.path.join(split_dir,"*.B*.r*.fastq")))[0]
    
    # Set up the cutadapt jobs as a group
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
        fqr1_out = '.'.join(os.path.basename(fqr1_in).split('.')[0:-1]) \
                   + '.trimmed.fastq'
        fqr2_out = '.'.join(os.path.basename(fqr2_in).split('.')[0:-1]) \
                   + '.trimmed.fastq'
        cutadapt_cmd.add_args('-o',fqr2_out,
                              '-p',fqr1_out)
        # Input pair (nb R2 then R1)
        cutadapt_cmd.add_args(fqr2_in,fqr1_in)
        # Submit the job
        print "Running %s" % cutadapt_cmd
        job = trim_reads.add(cutadapt_cmd,
                             wd=trim_dir,
                             name="cutadapt.%s" % os.path.basename(fqr1_in),
                             log_dir=log_dir)
        print "Job: %s" % job
    trim_reads.close()
    trim_reads.wait()
    exit_code = max([j.exit_code for j in trim_reads.jobs])
    if exit_code != 0:
        logging.critical("Read trimming stage failed (exit code %d)"
                         % exit_code)
        sys.exit(1)
    print "*** Read trimming stage completed ***"

    # Collect the files for the contamination filter step
    fastq_pairs = pair_fastqs(
        glob.glob(os.path.join(trim_dir,"*.r*.trimmed.fastq")))[0]

    # Set up the contaminant filter jobs as a group
    contaminant_filter = sched.group("contaminant_filter.%s" %
                                     os.path.basename(fastqs[0]))
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
        # Output pair (nb R2 then R1)
        ##fqr1_out = '.'.join(os.path.basename(fqr1_in).split('.')[0:-1]) \
        ##           + '.filtered.fastq'
        ##fqr2_out = '.'.join(os.path.basename(fqr2_in).split('.')[0:-1]) \
        ##           + '.filtered.fastq'
        ##contaminant_filter_cmd.add_args(fqr1_out,fqr1_out)
        # Submit the job
        print "Running %s" % contaminant_filter_cmd
        job = contaminant_filter.add(
            contaminant_filter_cmd,
            wd=contaminant_filter_dir,
            name="contaminant_filter.%s" %
            os.path.basename(fqr1_in),
            log_dir=log_dir)
        print "Job: %s" % job
    contaminant_filter.close()
    contaminant_filter.wait()
    exit_code = max([j.exit_code for j in contaminant_filter.jobs])
    if exit_code != 0:
        logging.critical("Contaminant filtering stage failed (exit code %d)"
                         % exit_code)
        sys.exit(exit_code)
    print "*** Contaminant filtering stage completed ***"

    # Finish
    sched.stop()
