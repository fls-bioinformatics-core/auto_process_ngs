#!/usr/bin/env python
#
#     fastq_utils.py: utility functions for operating on fastq files
#     Copyright (C) University of Manchester 2016-17 Peter Briggs
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
- get_read_number: get the read number (1 or 2) from a Fastq file
- pair_fastqs: automagically pair up FASTQ files

"""

#######################################################################
# Imports
#######################################################################

import os
import gzip
import logging
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

def get_read_number(fastq):
    """
    Get the read number (1 or 2) from a Fastq file

    Arguments:
      fastq (str): path to a Fastq file

    Returns:
      Integer: read number (1 or 2) extracted from the first read.
    """
    for r in FastqIterator(fastq):
        seq_id = r.seqid
        break
    return int(seq_id.pair_id)

def pair_fastqs(fastqs):
    """
    Automagically pair up FASTQ files

    Given a list of FASTQ files, generate a list of R1/R2
    pairs by examining the header for the first read in
    each file.

    Arguments:
      fastqs (list): list of paths to FASTQ files which
        will be paired.

    Returns:
      Tuple: pair of lists of the form (paired,unpaired),
        where `paired` is a list of tuples consisting of
        FASTQ R1/R2 pairs and `unpaired` is a list of
        FASTQs which couldn't be paired.
    """
    fq_pairs = []
    seq_ids = {}
    bad_files = []
    for fq in [os.path.abspath(fq) for fq in fastqs]:
        # Get header from first read
        seq_id = None
        for r in FastqIterator(fq):
            seq_id = r.seqid
            break
        if seq_id is None:
            logging.debug("'Bad' file: %s" % fq)
            bad_files.append(fq)
            continue
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
    # Sort pairs into order
    fq_pairs = sorted(fq_pairs,lambda x,y: cmp(x[0],y[0]))
    unpaired = sorted(seq_ids.keys() + bad_files)
    # Return paired and upaired fastqs
    return (fq_pairs,unpaired)
