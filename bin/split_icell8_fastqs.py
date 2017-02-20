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
from bcftbx.TabFile import TabFile
from bcftbx.utils import mkdir
from auto_process_ngs.icell8_utils import ICell8WellList
from auto_process_ngs.utils import OutputFiles

######################################################################
# Magic numbers
######################################################################

MAX_OPEN_FILES = 100
INLINE_BARCODE_LENGTH = 11
UMI_LENGTH = 10
INLINE_BARCODE_QUALITY_CUTOFF = 10
UMI_QUALITY_CUTOFF = 30
DEFAULT_BATCH_SIZE = 500

######################################################################
# Main
######################################################################

if __name__ == "__main__":
    # Handle the command line
    p = argparse.ArgumentParser()
    p.add_argument("FQ_R1",help="R1 FASTQ file")
    p.add_argument("FQ_R2",help="Matching R2 FASTQ file")
    p.add_argument("-w","--well-list",
                   dest="well_list_file",default=None,
                   help="iCell8 'well list' file")
    p.add_argument("-m","--mode",
                   dest="splitting_mode",default="barcodes",
                   choices=["barcodes","batch","none"],
                   help="how to split the input FASTQs (default: "
                   "'barcodes')")
    p.add_argument("-s","--size",type=int,
                   dest="batch_size",default=DEFAULT_BATCH_SIZE,
                   help="number of reads per batch in 'batch' mode "
                   "(default: %d)" % DEFAULT_BATCH_SIZE)
    p.add_argument("-b","--basename",
                   default="icell8",
                   help="basename for output FASTQ files (default: "
                   "'icell8')")
    p.add_argument("-o","--outdir",
                   dest="out_dir",default=None,
                   help="directory to write output FASTQ files to "
                   "(default: current directory)")
    p.add_argument("-l",
                   type=int,dest='inline_barcode_length',
                   default=INLINE_BARCODE_LENGTH,
                   help="length of inline barcodes (default: %d)" %
                   INLINE_BARCODE_LENGTH)
    args = p.parse_args()
    
    # Initialise positions for inline barcode and UMI
    bc_start = 0
    umi_start = args.inline_barcode_length + bc_start
    umi_end = umi_start + UMI_LENGTH

    # Convert quality cutoffs to character encoding
    barcode_quality_cutoff = chr(INLINE_BARCODE_QUALITY_CUTOFF + 33)
    umi_quality_cutoff = chr(UMI_QUALITY_CUTOFF + 33)

    # Get well list and expected barcodes
    well_list = ICell8WellList(args.well_list_file)
    expected_barcodes = well_list.barcodes()
    print "%d expected barcodes" % len(expected_barcodes)

    # Count barcodes and rejections
    barcodes = {}
    filtered = 0
    unassigned = 0
    nopenfiles = 0
    
    # Output Fastqs
    if args.out_dir is not None:
        mkdir(args.out_dir)
    output_fqs = OutputFiles(base_dir=args.out_dir)
    basename = args.basename
    batch_number = 0

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
        # Check that barcode is one we expected
        if inline_barcode not in expected_barcodes:
            assignment = "unassigned"
            unassigned += 1
        else:
            # Do filtering and splitting
            assignment = inline_barcode
            try:
                barcodes[inline_barcode]['unfiltered'] += 1
            except KeyError:
                logging.debug("New barcode: %s" % inline_barcode)
                barcodes[inline_barcode] = {
                    'unfiltered': 1,
                    'failed_barcode': 0,
                    'failed_UMI': 0,
                }
                # Barcode quality filter
                min_barcode_quality = min(r1.quality[bc_start:umi_start])
                if min_barcode_quality < barcode_quality_cutoff:
                    logging.debug("Read #%d: failed barcode quality filter"
                                  % n)
                    assignment = 'failed_barcode_quality'
                    barcodes[inline_barcode]['failed_barcode'] += 1
                # UMI quality filter
                min_umi_quality = min(r1.quality[umi_start:umi_end])
                if min_umi_quality < umi_quality_cutoff:
                    logging.debug("Read #%d: failed UMI quality filter"
                                  % n)
                    assignment = 'failed_UMI_quality'
                    barcodes[inline_barcode]['failed_UMI'] += 1
        # Assign read pair to appropriate output files
        if assignment == inline_barcode:
            filtered += 1
            if args.splitting_mode == "batch":
                # Output to a batch-specific file pair
                if (filtered - 1) % args.batch_size == 0:
                    batch_number += 1
                assignment = "B%03d" % batch_number
            elif args.splitting_mode == "none":
                # Output to a single file pair
                assignment = "filtered"
        # Write read pair
        fq_r1 = "%s_R1" % assignment
        fq_r2 = "%s_R2" % assignment
        if fq_r1 not in output_fqs:
            try:
                # Try to reopen file and append
                output_fqs.open(fq_r1,append=True)
            except KeyError:
                # Open new file
                output_fqs.open(fq_r1,
                                "%s.%s.r1.fastq" % (basename,assignment))
            nopenfiles += 1
            ##print "%d: %d files open..." % (n,nopenfiles)
        output_fqs.write(fq_r1,"%s" % r1)
        if fq_r2 not in output_fqs:
            try:
                # Try to reopen file and append
                output_fqs.open(fq_r2,append=True)
            except KeyError:
                # Open new file
                output_fqs.open(fq_r2,
                                "%s.%s.r2.fastq" % (basename,assignment))
            nopenfiles += 1
            ##print "%d: %d files open..." % (n,nopenfiles)
        output_fqs.write(fq_r2,"%s" % r2)
        # FIXME close the files if it looks like we have too
        # many open at once (to avoid IOError [Errno 24])
        if nopenfiles > MAX_OPEN_FILES:
            logging.debug("*** Closing output files ***")
            output_fqs.close()
            nopenfiles = 0
    # Close output files
    output_fqs.close()
    # Report barcode and filtering statistics
    barcode_list = sorted(barcodes.keys(),
                          cmp=lambda x,y: cmp(barcodes[x],barcodes[y]),
                          reverse=True)
    for barcode in barcode_list:
        line = [barcode,
                well_list.sample(barcode),
                barcodes[barcode]['unfiltered'],
                barcodes[barcode]['failed_barcode'],
                barcodes[barcode]['failed_UMI']
        ]
        print "%s" % '\t'.join([str(x) for x in line])
    assigned = sum([barcodes[b]['unfiltered'] for b in barcode_list])
    total_unfiltered = assigned + unassigned
    failed_barcode = sum([barcodes[b]['failed_barcode']
                          for b in barcode_list])
    failed_umi = sum([barcodes[b]['failed_UMI']
                      for b in barcode_list])
    print "For all barcodes:"
    print "-----------------"
    print "Number of barcodes found:\t%d/%d" % (len(barcode_list),
                                                len(expected_barcodes))
    print "Total reads (unfiltered):\t%d" % total_unfiltered
    print "Unassigned reads        :\t%d" % unassigned
    print "Failed barcode quality  :\t%d" % failed_barcode
    print "Failed UMI quality      :\t%d" % failed_umi
    print "Total reads (filtered)  :\t%d" % filtered
