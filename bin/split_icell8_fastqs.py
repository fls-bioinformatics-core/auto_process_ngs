#!/usr/bin/env python
#
#     split_icell8_fastqs.py: splits fastq from Wafergen iCell8
#     Copyright (C) University of Manchester 2017 Peter Briggs
#
"""
split_icell8_fastqs.py

Utility to split FASTQ pair from Wafergen iCell8 into individual
FASTQ files based on the inline barcodes in read 1.

"""

######################################################################
# Imports
######################################################################

import sys
import logging
import argparse
from itertools import izip
from bcftbx.FASTQFile import FastqIterator

######################################################################
# Magic numbers
######################################################################

INLINE_BARCODE_LENGTH = 11
UMI_LENGTH = 10
INLINE_BARCODE_QUALITY_CUTOFF = 10
UMI_QUALITY_CUTOFF = 30

######################################################################
# Main
######################################################################

if __name__ == "__main__":
    # Handle the command line
    p = argparse.ArgumentParser()
    p.add_argument("FQ_R1",help="R1 FASTQ file")
    p.add_argument("FQ_R2",help="Matching R2 FASTQ file")
    p.add_argument("-l",
                   type=int,dest='inline_barcode_length',
                   default=INLINE_BARCODE_LENGTH,
                   help="length of inline barcodes (default: %d)" %
                   INLINE_BARCODE_LENGTH)
    p.add_argument("-b","--basename",
                   default="icell8",
                   help="basename for output FASTQ files (default: "
                   "'icell8'")
    args = p.parse_args()
    
    # Initialise positions for inline barcode and UMI
    bc_start = 0
    umi_start = args.inline_barcode_length + bc_start
    umi_end = umi_start + UMI_LENGTH

    # Convert quality cutoffs to character encoding
    barcode_quality_cutoff = chr(INLINE_BARCODE_QUALITY_CUTOFF + 33)
    umi_quality_cutoff = chr(UMI_QUALITY_CUTOFF + 33)

    # Count barcodes and rejections
    barcodes = {}
    nrejected = 0
    
    # Output Fastqs
    output_fqs = {}
    basename = args.basename

    # Iterate over read pairs from the Fastqs
    for n,read_pair in enumerate(izip(FastqIterator(args.FQ_R1),
                                      FastqIterator(args.FQ_R2))):
        # Check reads are a pair
        r1,r2 = read_pair
        if not r1.seqid.is_pair_of(r2.seqid):
            logging.critical("Read #%d: unpaired read headers detected"
                             % n)
            sys.exit(1)
        # Extract the inline barcode and UMI
        try:
            inline_barcode = r1.sequence[bc_start:umi_start]
            umi = r1.sequence[umi_start:umi_end]
        except IndexError:
            logging.critical("Read #%d: unable to extract inline "
                             "barcode or UMI" % n)
            sys.exit(1)
        # Count barcodes
        try:
            barcodes[inline_barcode] += 1
        except KeyError:
            print "New barcode: %s" % inline_barcode
            barcodes[inline_barcode] = 1
        # Apply quality filters
        min_barcode_quality = min(r1.quality[bc_start:umi_start])
        if min_barcode_quality < barcode_quality_cutoff:
            logging.warning("Read #%d: failed barcode quality filter"
                            % n)
            inline_barcode = 'rejected'
        min_umi_quality = min(r1.quality[umi_start:umi_end])
        if min_umi_quality < umi_quality_cutoff:
            logging.warning("Read #%d: failed UMI quality filter"
                            % n)
            inline_barcode = 'rejected'
        if inline_barcode == 'rejected':
            nrejected += 1
        # Write to appropriate output files
        fq_r1 = "%s_R1" % inline_barcode
        if fq_r1 not in output_fqs:
            output_fqs[fq_r1] = open("%s.%s.r1.fastq" %
                                     (basename,inline_barcode),'w')
        output_fqs[fq_r1].write("%s\n" % r1)
        fq_r2 = "%s_R2" % inline_barcode
        if fq_r2 not in output_fqs:
            output_fqs[fq_r2] = open("%s.%s.r2.fastq" %
                                     (basename,inline_barcode),'w')
        output_fqs[fq_r2].write("%s\n" % r2)

    # Close output files
    for fq in output_fqs:
        output_fqs[fq].close()
        
    # Report barcodes
    barcode_list = sorted(barcodes.keys(),
                          cmp=lambda x,y: cmp(barcodes[x],barcodes[y]),
                          reverse=True)
    for barcode in barcode_list:
        print "%s\t%d" % (barcode,barcodes[barcode])
    total_reads = sum([barcodes[b] for b in barcode_list])
    print "Total reads:\t%d" % total_reads
    print "# Rejected :\t%d" % nrejected
