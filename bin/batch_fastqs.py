#!/usr/bin/env python
#
#     batch_fastqs.py: splits fastqs into batches
#     Copyright (C) University of Manchester 2017 Peter Briggs
#
"""
batch_fastqs.py

Utility to split FASTQs into a set of smaller "batch" FASTQ files
"""

######################################################################
# Imports
######################################################################

import os
import sys
import argparse
import time
import logging
from itertools import izip
from bcftbx.utils import mkdir
from bcftbx.FASTQFile import get_fastq_file_handle
from auto_process_ngs.utils import AnalysisFastq
from auto_process_ngs.fastq_utils import pair_fastqs

######################################################################
# Magic numbers
######################################################################

DEFAULT_BATCH_SIZE = 5000000

######################################################################
# Functions
######################################################################

def split_fastqs(fastq_pairs,basename,batch_size):
    start_time = time.time()
    print "Starting at %s" % time.ctime()
    # Initialise loop variables
    current_batch = None
    fpr1_out = None
    fpr2_out = None
    ntot = 1
    # Loop over file pairs
    for fastq_pair in fastq_pairs:
        # Open the input files
        fpr1 = get_fastq_file_handle(fastq_pair[0])
        fpr2 = get_fastq_file_handle(fastq_pair[1])
        # Iterate over file pair
        for i,pair in enumerate(izip(fpr1,fpr2)):
            # Read number
            n = i/4 + ntot
            # Batch number
            batch_number = (n-1)/batch_size
            # New batch?
            if current_batch != batch_number:
                print "-- Start of a new batch B%03d" % batch_number
                try:
                    fpr1_out.close()
                except AttributeError:
                    pass
                fpr1_out = open("%s.B%03d.r1.fastq" %
                                (basename,batch_number),'w')
                try:
                    fpr2_out.close()
                except AttributeError:
                    pass
                fpr2_out = open("%s.B%03d.r2.fastq" %
                                (basename,batch_number),'w')
                current_batch = batch_number
            if (n % 1000000) == 0 and (i % 4) == 0:
                print "-- Examining read pair #%d B%03d (%s)" % (n,
                                                                 batch_number,
                                                                 time.ctime())
            # Write the read pair
            fpr1_out.write(pair[0])
            fpr2_out.write(pair[1])
    fpr1_out.close()
    fpr2_out.close()
    print "Finished at %s" % time.ctime()
    print "Handled %d reads in %.0fs" % (n,
                                         time.time()-start_time)

def merge_fastqs(fastq_pairs,basename):
    start_time = time.time()
    print "Starting at %s" % time.ctime()
    ntot = 1
    # Open output files
    fpr1_out = open("%s.r1.fastq" % basename,'w')
    fpr2_out = open("%s.r2.fastq" % basename,'w')
    # Loop over file pairs
    for fastq_pair in fastq_pairs:
        # Open the input files
        fpr1 = get_fastq_file_handle(fastq_pair[0])
        fpr2 = get_fastq_file_handle(fastq_pair[1])
        for i,pair in enumerate(izip(fpr1,fpr2)):
            # Read number
            n = i/4 + ntot
            if (n % 1000000) == 0:
                print "-- Examining read pair #%d %s (%s)" % (n,
                                                              fastq_pair[0],
                                                              time.ctime())
            # Write the read pair
            fpr1_out.write(pair[0])
            fpr2_out.write(pair[1])
        ntot = n
    fpr1_out.close()
    fpr2_out.close()
    print "Finished at %s" % time.ctime()
    print "Handled %d reads in %.0fs" % (ntot,
                                         time.time()-start_time)

######################################################################
# Main
######################################################################

if __name__ == "__main__":
    # Handle the command line
    p = argparse.ArgumentParser()
    p.add_argument("fastqs",metavar="FASTQ_R1 FASTQ_R2",nargs="+",
                   help="FASTQ R1/R2 pair")
    p.add_argument("-s","--size",type=int,
                   dest="batch_size",default=DEFAULT_BATCH_SIZE,
                   help="number of reads per batch (default: %d)" %
                   DEFAULT_BATCH_SIZE)
    p.add_argument("-m","--merge",
                   dest="merge",action="store_true",
                   help="merge input pairs into a single pair")
    p.add_argument("-b","--basename",default=None,
                   help="basename for output FASTQ files (default: "
                   "sample name from FASTQ_R1)")
    p.add_argument("-o","--outdir",
                   dest="out_dir",default=None,
                   help="directory to write output FASTQ files to "
                   "(default: current directory)")
    args = p.parse_args()

    # Batch size
    batch_size = args.batch_size

    # Output directory
    if args.out_dir is not None:
        out_dir = os.path.abspath(args.out_dir)
        mkdir(out_dir)
    else:
        out_dir = os.getcwd()

    # Output file basename
    basename = args.basename
    if basename is None:
        basename = AnalysisFastq(args.fastq_r1).sample_name
    basename = os.path.join(out_dir,basename)

    # Sort input files into pairs
    pairs,unpaired = pair_fastqs(args.fastqs)
    if unpaired:
        print "Unpaired FASTQS:"
        for fq in unpaired:
            print "-- %s" % fq
        logging.critical("Couldn't pair all FASTQ files")
        sys.exit(1)

    # Do splitting or merging
    if args.merge:
        merge_fastqs(pairs,basename)
    else:
        split_fastqs(pairs,basename,batch_size)
