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

__version__ = "0.0.2"

import FASTQFile
import IlluminaData
import optparse
import logging
import sys
import os

def count_barcodes(fastqs):
    """Count the index sequences across one or more Fastq files

    Arguments:
      fastqs: list of one or more Fastq files to read barcodes from

    Returns:
      'counts' dictionary where counts[SEQ] holds the number of
      times index sequence SEQ occurs.

    """
    counts = dict()
    nreads = 0
    for fastq in fastqs:
        print "Reading in data from %s" % fastq
        for read in FASTQFile.FastqIterator(fastq):
            nreads += 1
            seq = read.seqid.index_sequence
            if not seq:
                raise ValueError,"No index sequence for read! %s" % read
            # Check if we've already encountered this sequence
            if seq in counts:
                # Already seen
                counts[seq] += 1
            else:
                # Novel sequence
                counts[seq] = 1
    # Return the counts dictionary
    return counts

def ordered_seqs(counts):
    """Return list of sequences sorted most to least numerous

    Arguments:
      counts: 'counts' dictionary returned by 'count_barcodes'

    """
    return sorted(counts.keys(),cmp=lambda x,y: cmp(counts[y],counts[x]))

def nreads(counts):
    """Return number of reads

    Arguments:
      counts: 'counts' dictionary returned by 'count_barcodes'

    """
    return sum([counts[x] for x in counts])

def report(counts,nseqs=20,exclude_ns=False):
    """Print a summary of the most numerous barcode sequences

    Arguments:
      counts: 'counts' dictionary returned by 'count_barcodes'
      nseqs: number of barcode sequences to report
      exclude_ns: if True then only report sequences which
        don't include any bases flagged as 'N'

    """
    seqs = ordered_seqs(counts)
    if exclude_ns:
        seqs = filter(lambda x: x if x.count('N') == 0 else None,seqs)
    n = nreads(counts)
    print "Number of reads: %d" % n
    print "Number of unique barcode sequences: %d" % len(seqs)
    # Output top barcodes
    print "Top %d barcode sequences%s" % (nseqs,
                                          (" (no 'N's)" if exclude_ns else ""))
    for i,seq in enumerate(seqs[:nseqs]):
        print "% 8d\t%s\t% 8d\t% 5.1f%%" % (i+1,seq,counts[seq],float(counts[seq])/n*100.0)

def output(counts,filen):
    """Write all sequences to a file

    Arguments:
      counts: 'counts' dictionary returned by 'count_barcodes'
      filen: name and path of output file to write counts to

    """
    # Output all barcodes to file
    print "Writing all counts to file '%s'" % counts_file
    fp = open(filen,'w')
    for i,seq in enumerate(ordered_seqs(counts)):
        fp.write("%d\t%s\t%d\n" % (i+1,seq,counts[seq]))
    fp.close()

if __name__ == '__main__':
    # Handle command line
    p = optparse.OptionParser(usage="\n\t%prog FASTQ [FASTQ...]\n\t%prog DIR",
                              version=__version__,
                              description="Count the index/barcode sequences from one or "
                              "more Fastq files and report the most numerous. Sequences will "
                              "be pooled from all specified Fastqs before being analysed. "
                              "Fastqs can either be specified explicitly, or if a single "
                              "directory is supplied instead then this will be assumed to "
                              "be an output directory from CASAVA or bclToFastq and files "
                              "will be processed on a per-lane basis.")
    p.add_option('-n',action='store',dest='n',default=20,type='int',
                 help="report top N barcodes (default 20)")
    p.add_option('-o',action='store',dest='counts_file',default=None,
                 help="output all counts to tab-delimited file COUNTS_FILE")
    options,args = p.parse_args()
    # Check arguments and acquire fastqs
    if not args:
        p.error("Need to supply at least one input Fastq file, or a bclToFastq output "
                "directory")
    if len(args) == 1 and os.path.isdir(args[0]):
        # Dealing with a bclToFastq output dir
        illumina_data = IlluminaData.IlluminaData(os.path.dirname(args[0]),
                                                  unaligned_dir=os.path.basename(args[0]))
        # Assign fastqs to lanes (R1 only)
        fastq_in_lane = dict()
        for p in illumina_data.projects:
            for s in p.samples:
                for f in s.fastq_subset(read_number=1,full_path=True):
                    lane = IlluminaData.IlluminaFastq(f).lane_number
                    if lane not in fastq_in_lane:
                        fastq_in_lane[lane] = []
                    fastq_in_lane[lane].append(f)
        if illumina_data.undetermined:
            for s in illumina_data.undetermined.samples:
                for f in s.fastq_subset(read_number=1,full_path=True):
                    lane = IlluminaData.IlluminaFastq(f).lane_number
                    if lane not in fastq_in_lane:
                        fastq_in_lane[lane] = []
                    fastq_in_lane[lane].append(f)
        lanes = fastq_in_lane.keys()
        lanes.sort()
        print "Lanes: %s" % lanes
        # Count barcodes in each lane
        for lane in lanes:
            print "*** Counting barcodes in lane %d ***" % lane
            counts = count_barcodes(fastq_in_lane[lane])
            report(counts,nseqs=options.n)
            if options.counts_file:
                counts_file = "%s.lane%s%s" % (os.path.splitext(options.counts_file)[0],
                                               lane,
                                               os.path.splitext(options.counts_file)[1])
                output(counts,counts_file)
    else:
        # One or more fastq files
        counts = count_barcodes(args)
        report(counts,nseqs=options.n)
        if options.counts_file:
            output(counts,options.counts_file)
    
    
