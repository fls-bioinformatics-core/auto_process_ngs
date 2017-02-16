#!/usr/bin/env python
#
#     icell8_trim_reads.py: trims reads from Wafergen iCell8
#     Copyright (C) University of Manchester 2017 Peter Briggs
#
"""
trim_icell8_reads.py

Utility to process a set of FASTQ pairs from Wafergen iCell8 that
have been previous filtered on barcode id/quality and UMI quality,
and perform read trimming:

- Remove sequencing primers
- Remove poly-A/T and poly-N sequences
- Remove trailing bases below quality threshold

Post trimming short reads (<=20 bases) are removed.
"""

######################################################################
# Imports
######################################################################

import os
import sys
import logging
import argparse
from bcftbx.utils import mkdir
from bcftbx.FASTQFile import FastqIterator
from auto_process_ngs.applications import Command
from auto_process_ngs.simple_scheduler import SimpleScheduler
from auto_process_ngs import envmod

import logging
logging.basicConfig(format='%(levelname) 8s: %(message)s')

import auto_process_ngs.settings
__settings = auto_process_ngs.settings.Settings()

#try:
#    __modulefiles = __settings.modulefiles['run_qc']
#except KeyError:
#    # No environment modules specified
#    __modulefiles = None

######################################################################
# Functions
######################################################################

def pair_fastqs(fastqs):
    """
    Automagically pair up FASTQs files
    """
    fq_pairs = []
    seq_ids = {}
    for fq in [os.path.abspath(fq) for fq in fastqs]:
        # Get header from first read
        for r in FastqIterator(fq):
            seq_id = r.seqid
            break
        fq_pair = None
        for fq1 in seq_ids:
            if seq_id.is_pair_of(seq_ids[fq1]):
                # Found a pair
                if seq_id.pair_id == '1':
                    fq_pair = (fq,fq1)
                else:
                    fq_pair = (fq1,fq)
                fq_pairs.append(fq_pair)
                logging.debug("*** Paired: %s\n"
                              "          : %s" % fq_pair)
                # Remove paired fastq
                del(seq_ids[fq1])
                break
        if fq_pair is None:
            # Unable to pair, store for now
            logging.debug("Unpaired: %s" % fq)
            seq_ids[fq] = seq_id
    # Return paired and upaired fastqs
    return (fq_pairs,seq_ids.keys())

######################################################################
# Main
######################################################################

if __name__ == "__main__":
    # Handle the command line
    p = argparse.ArgumentParser()
    p.add_argument("fastq",metavar="FASTQ",nargs="+",
                   help="FASTQ file (R1/R2 pairs)")
    p.add_argument("-o","--outdir",
                   dest="out_dir",default=None,
                   help="directory to write output FASTQ files to "
                   "(default: current directory)")
    args = p.parse_args()

    # Get FASTQ pairs
    fq_pairs,unpaired_fqs = pair_fastqs(args.fastq)
    if unpaired_fqs:
        print "*** Unpaired FASTQs:"
        for fq in unpaired_fqs:
            print "-- %s" % fq
        logging.critical("Failed to pair all FASTQs")
        sys.exit(1)    
    # Report pairs
    for pair in fq_pairs:
        print "-- %s\n   %s" % (pair[0],pair[1])

    # Make output dir
    if args.out_dir is not None:
        out_dir = os.path.abspath(args.out_dir)
        mkdir(out_dir)
    else:
        out_dir = os.getcwd()

    # Schedule cutadapt jobs
    qc_runner = __settings.runners.qc
    max_jobs = __settings.general.max_concurrent_jobs
    sched = SimpleScheduler(runner=qc_runner,
                            max_concurrent=max_jobs)
    sched.start()
    for fq_pair in fq_pairs:
        fqr1_in,fqr2_in = fq_pair
        cutadapt_cmd = Command('cutadapt',
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
        job = sched.submit(cutadapt_cmd,
                           wd=out_dir,
                           name="cutadapt.%s" % os.path.basename(fqr1_in),
                           log_dir=out_dir)
        print "Job: %s" % job
    # Wait for the scheduler to run all jobs
    sched.wait()
    sched.stop()
