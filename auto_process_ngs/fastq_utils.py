#!/usr/bin/env python
#
#     fastq_utils.py: utility functions for operating on fastq files
#     Copyright (C) University of Manchester 2016 Peter Briggs
#
########################################################################
#
# fastq_utils.py
#
#########################################################################

"""
fastq_utils.py

Utility functions for operating on Fastq files:

- assign_barcodes_single_end: extract and assign inline barcodes

"""

#######################################################################
# Imports
#######################################################################

import gzip
from bcftbx.FASTQFile import FastqIterator

#######################################################################
# Functions
#######################################################################

def assign_barcodes_single_end(fastq_in,fastq_out,n=5):
    """
    Extract inline barcodes and assign to Fastq read headers

    Strips the first n bases from each read of the input
    FASTQ file and assigns it to the index sequence for that
    read in the output file.

    If the supplied output file name ends with '.gz' then it
    will be gzipped.

    Arguments:
      fastq_in (str): input FASTQ file (can be gzipped)
      fastq_out (str): output FASTQ file (will be gzipped if
        ending with '.gz')
      n (integer): number of bases to extract and assign as
        index sequence (default: 5)

    Returns:
      Integer: number of reads processed.

    """
    if fastq_out.endswith('.gz'):
        fp = gzip.GzipFile(filename=fastq_out,mode='wb')
    else:
        fp = open(fastq_out,'w')
    print "Processing reads from %s" % fastq_in
    nread = 0
    for read in FastqIterator(fastq_in):
        # Extract new barcode sequence
        barcode = read.sequence[:n]
        # Truncate sequence and quality accordingly
        sequence = read.sequence[n:]
        quality = read.quality[n:]
        # Assign new values and write to output
        read.seqid.index_sequence = barcode
        read.sequence = sequence
        read.quality = quality
        fp.write("%s\n" % read)
        nread += 1
    print "Finished (%d reads processed)" % nread
    return nread
