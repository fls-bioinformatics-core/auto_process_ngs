#!/bin/env python
#
#     count_barcodes.py: count index sequences (barcodes) in fastqs
#     Copyright (C) University of Manchester 2014 Peter Briggs
#
#########################################################################
#
# count_barcodes.py
#
#########################################################################

"""
count_barcodes.py

Counts the index sequences (i.e. barcodes) for all reads in one or more
Fastq files, and reports the most numerous.

"""

__version__ = "0.0.1"

import FASTQFile
import optparse
import logging
import sys
if __name__ == '__main__':
    # Handle command line
    p = optparse.OptionParser(usage="%prog FASTQ [FASTQ...]",
                              version=__version__,
                              description="Count the index/barcode sequences from one or "
                              "more Fastq files and report the most numerous. Sequences will "
                              "be pooled from all specified Fastqs before being analysed.")
    p.add_option('-n',action='store',dest='n',default=20,type='int',
                 help="report top N barcodes (default 20)")
    p.add_option('-o',action='store',dest='counts_file',default=None,
                 help="output all counts to tab-delimited file COUNTS_FILE")
    options,args = p.parse_args()
    if not args:
        p.error("Need to supply at least one input Fastq file")
    # Count the index sequences
    counts = dict()
    nreads = 0
    for fastq in args:
        print "Reading in data from %s" % fastq
        for read in FASTQFile.FastqIterator(fastq):
            nreads += 1
            seq = read.seqid.index_sequence
            if not seq:
                logging.error("No index sequence for read!")
                sys.exit(1)
            # Check if we've already encountered this sequence
            if seq in counts:
                # Already seen
                counts[seq] += 1
            else:
                # Novel sequence
                counts[seq] = 1
    # Finished processing files
    print "Number of reads: %d" % nreads
    print "Number of unique barcode sequences: %d" % len(counts.keys())
    # Sort barcodes into order of most common
    ordered_seqs = sorted(counts.keys(),cmp=lambda x,y: cmp(counts[y],counts[x]))
    # Output top barcodes
    print "Top %d barcode sequences" % options.n
    for i,seq in enumerate(ordered_seqs[:options.n]):
        print "% 8d\t%s\t%d" % (i+1,seq,counts[seq])
    print "Top %d barcode sequences (no 'N's)" % options.n
    for i,seq in enumerate(
            filter(lambda x: x if x.count('N') == 0 else None,ordered_seqs)[:options.n]):
        print "% 8d\t%s\t%d" % (i+1,seq,counts[seq])
    # Output all barcodes to file
    if options.counts_file is not None:
        print "Writing all counts to file '%s'" % options.counts_file
        fp = open(options.counts_file,'w')
        for i,seq in enumerate(ordered_seqs):
            fp.write("%d\t%s\t%d\n" % (i+1,seq,counts[seq]))
        fp.close()
    
    
