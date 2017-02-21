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
from collections import Iterator
from bcftbx.FASTQFile import FastqIterator
from bcftbx.TabFile import TabFile
from bcftbx.utils import mkdir
from auto_process_ngs.icell8_utils import ICell8WellList
from auto_process_ngs.fastq_utils import pair_fastqs
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
# Classes
######################################################################

class ICell8ReadPair(object):
    def __init__(self,r1,r2):
        if not r1.seqid.is_pair_of(r2.seqid):
            raise Exception("Reads are not paired")
        self._r1 = r1
        self._r2 = r2
    @property
    def r1(self):
        return self._r1
    @property
    def r2(self):
        return self._r2
    @property
    def barcode(self):
        return self._r1.sequence[0:INLINE_BARCODE_LENGTH]
    @property
    def umi(self):
        return self._r1.sequence[INLINE_BARCODE_LENGTH:
                                 INLINE_BARCODE_LENGTH+UMI_LENGTH]
    @property
    def min_barcode_quality(self):
        return min(self._r1.quality[0:INLINE_BARCODE_LENGTH])
    @property
    def min_umi_quality(self):
        return min(self._r1.quality[INLINE_BARCODE_LENGTH:
                                    INLINE_BARCODE_LENGTH+UMI_LENGTH])

class ICell8FastqIterator(Iterator):
    def __init__(self,fqr1,fqr2):
        self._fqr1 = FastqIterator(fqr1)
        self._fqr2 = FastqIterator(fqr2)
    def next(self):
        return ICell8ReadPair(self._fqr1.next(),
                              self._fqr2.next())

######################################################################
# Main
######################################################################

if __name__ == "__main__":
    # Handle the command line
    p = argparse.ArgumentParser()
    p.add_argument("FQ_R1",help="R1 FASTQ file")
    p.add_argument("FQ_R2",help="Matching R2 FASTQ file")
    p.add_argument("FQ",nargs='*',help="Additional FASTQ file pairs")
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
    p.add_argument("-f","--stats_file",
                   dest="stats_file",default=None,
                   help="output file for statistics")
    args = p.parse_args()

    # Convert quality cutoffs to character encoding
    barcode_quality_cutoff = chr(INLINE_BARCODE_QUALITY_CUTOFF + 33)
    umi_quality_cutoff = chr(UMI_QUALITY_CUTOFF + 33)

    # Get well list and expected barcodes
    well_list = ICell8WellList(args.well_list_file)
    expected_barcodes = well_list.barcodes()
    print "%d expected barcodes" % len(expected_barcodes)

    # Count barcodes and rejections
    unassigned = 0
    filtered = 0
    unfiltered_counts = {}
    filtered_counts = {}
    unfiltered_umis = {}
    filtered_umis = {}

    # Input Fastqs
    fastqs = [args.FQ_R1,args.FQ_R2]
    for fq in args.FQ:
        fastqs.append(fq)
    fastqs = pair_fastqs(fastqs)[0]

    # Output Fastqs
    if args.out_dir is not None:
        mkdir(args.out_dir)
    output_fqs = OutputFiles(base_dir=args.out_dir)
    basename = args.basename

    # Iterate over pairs of Fastqs
    for fastq_pair in fastqs:
        # Iterate over read pairs from the Fastqs
        print "-- %s\n   %s" % fastq_pair
        for read_pair in ICell8FastqIterator(*fastq_pair):
            # Count the barcodes and unique UMIs
            inline_barcode = read_pair.barcode
            umi = read_pair.umi
            try:
                unfiltered_counts[inline_barcode] += 1
                if umi not in unfiltered_umis[inline_barcode]:
                    unfiltered_umis[inline_barcode].append(umi)
            except KeyError:
                unfiltered_counts[inline_barcode] = 1
                unfiltered_umis[inline_barcode] = []
            # Do filtering
            if inline_barcode not in expected_barcodes:
                assign_to = "unassigned"
                unassigned += 1
            elif read_pair.min_barcode_quality < barcode_quality_cutoff:
                assign_to = "failed_barcode"
            elif read_pair.min_umi_quality < umi_quality_cutoff:
                assign_to = "failed_umi"
            else:
                assign_to = inline_barcode
                filtered += 1
            logging.debug("%s" % '\t'.join([assign_to,
                                            inline_barcode,
                                            umi,
                                            read_pair.min_barcode_quality,
                                            read_pair.min_umi_quality]))
            # Post filtering counts
            if assign_to == inline_barcode:
                try:
                    filtered_counts[inline_barcode] += 1
                    if umi not in filtered_umis[inline_barcode]:
                        filtered_umis[inline_barcode].append(umi)
                except KeyError:
                    filtered_counts[inline_barcode] = 1
                    filtered_umis[inline_barcode] = []
                # Reassign read pair to appropriate output files
                if args.splitting_mode == "batch":
                    # Output to a batch-specific file pair
                    batch_number = filtered/args.batch_size
                    assign_to = "B%03d" % batch_number
                elif args.splitting_mode == "none":
                    # Output to a single file pair
                    assign_to = "filtered"
            # Write read pair
            fq_r1 = "%s_R1" % assign_to
            fq_r2 = "%s_R2" % assign_to
            if fq_r1 not in output_fqs:
                try:
                    # Try to reopen file and append
                    output_fqs.open(fq_r1,append=True)
                except KeyError:
                    # Open new file
                    output_fqs.open(fq_r1,
                                    "%s.%s.r1.fastq" % (basename,assign_to))
            output_fqs.write(fq_r1,"%s" % read_pair.r1)
            if fq_r2 not in output_fqs:
                try:
                    # Try to reopen file and append
                    output_fqs.open(fq_r2,append=True)
                except KeyError:
                    # Open new file
                    output_fqs.open(fq_r2,
                                    "%s.%s.r2.fastq" % (basename,assign_to))
            output_fqs.write(fq_r2,"%s" % read_pair.r2)
        # FIXME close the files if it looks like we have too
        # many open at once (to avoid IOError [Errno 24])
        if len(output_fqs) > MAX_OPEN_FILES:
            logging.debug("*** Closing output files ***")
            output_fqs.close()
    # Close output files
    output_fqs.close()

    # Report barcode and filtering statistics
    barcode_list = sorted(unfiltered_counts.keys(),
                          key=lambda bc: unfiltered_counts[bc],
                          reverse=True)
    if args.stats_file is not None:
        stats_file = TabFile(column_names=('Barcode',
                                           'Sample',
                                           'Nreads_pre_filter',
                                           'Unique_UMIs_pre_filter',
                                           'Nreads_post_filter',
                                           'Unique_UMIs_post_filter'))
        for barcode in barcode_list:
            if barcode in expected_barcodes:
                line = [barcode,
                        well_list.sample(barcode),
                        unfiltered_counts[barcode],
                        len(unfiltered_umis[barcode])]
            else:
                line = [barcode,
                        '',
                        unfiltered_counts[barcode],
                        len(unfiltered_umis[barcode])]
            if barcode in filtered_counts:
                line.extend([filtered_counts[barcode],
                             len(filtered_umis[barcode])])
            else:
                line.extend(['',''])
            stats_file.append(data=line)
        stats_file.write(args.stats_file,include_header=True)

    # Summary output to screen
    assigned = sum([unfiltered_counts[b] for b in barcode_list])
    total_reads = assigned + unassigned
    print "Summary:"
    print "--------"
    print "Number of barcodes      :\t%d" % len(barcode_list)
    print "Number of barcodes found:\t%d/%d" % (len(filtered_counts.keys()),
                                                len(expected_barcodes))
    print "Total reads (unfiltered):\t%d" % assigned
    print "Unassigned reads        :\t%d" % unassigned
    print "Total reads (filtered)  :\t%d" % filtered
