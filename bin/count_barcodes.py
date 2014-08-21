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

__version__ = "0.1.1"

import FASTQFile
import IlluminaData
import optparse
import logging
import sys
import os
from multiprocessing import Pool

#######################################################################
# Functions
#######################################################################

def count_barcodes_for_file(fastq):
    """Count the index sequences across a single Fastq file

    Arguments:
      fastq: Fastq file to read barcodes from

    Returns:
      'counts' dictionary where counts[SEQ] holds the number of
      times index sequence SEQ occurs.

    """
    counts = dict()
    nreads = 0
    print "Reading in data from %s" % fastq
    for read in FASTQFile.FastqIterator(fastq):
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

def count_barcodes(fastqs,n_processors=1):
    """Count the index sequences across one or more Fastq files

    Arguments:
      fastqs: list of one or more Fastq files to read barcodes from
      n_processors: number of processors to use (if 1>1 then uses
        the multiprocessing library to run the barcode counting
        using multiple cores).

    Returns:
      'counts' dictionary where counts[SEQ] holds the number of
      times index sequence SEQ occurs.

    """
    # Get counts for individual files
    print "Counting index sequences for each file"
    if n_processors > 1:
        # Multiple cores
        pool = Pool(n_processors)
        results = pool.map(count_barcodes_for_file,fastqs)
        pool.close()
        pool.join()
    else:
        # Single core
        results = map(count_barcodes_for_file,fastqs)
    # Merge the results
    print "Merging counts of index sequences"
    all_counts = dict()
    for counts in results:
        for seq in counts:
            if seq in all_counts:
                all_counts[seq] += counts[seq]
            else:
                all_counts[seq] = counts[seq]
    # Return all the counts
    return all_counts

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

def percent_reads(n,nreads):
    """Return n as a percentage of nreads

    """
    return float(n)/nreads*100.0

def nmismatches(seq1,seq2,truncate=True):
    """Return the number of mismatches between seq1 and seq2

    Return the number of mismatches between two sequences
    seq1 and seq2.

    If 'truncate' is True then the sequences are truncated
    before comparison, otherwise the mismatch counting is
    undefined and raises a ValueError.

    """
    if len(seq1) != len(seq2):
        if truncate:
            new_len = min(len(seq1),len(seq2))
            seq1 = seq1[:new_len]
            seq2 = seq2[:new_len]
        else:
            raise ValueError("mismatches undefined for sequences of unequal length")
    return sum(c1 != c2 for c1,c2 in zip(seq1,seq2))

def report(counts,nseqs=20,exclude_ns=False,fp=None,
           cutoff=None):
    """Print a summary of the most numerous barcode sequences

    Arguments:
      counts: 'counts' dictionary returned by 'count_barcodes'
      nseqs: number of barcode sequences to report
      exclude_ns: if True then only report sequences which
        don't include any bases flagged as 'N'
      cutoff: threshold fraction of reads below which not
        to report matches (default is to report everything)
      fp: specify output stream to write report to (default
        is to use stdout)

    """
    # Deal with cutoff
    if cutoff is not None:
        if cutoff < 1.0:
            cutoff = cutoff*100.0
    # Get sequences in frequency order
    seqs = ordered_seqs(counts)
    if exclude_ns:
        seqs = filter(lambda x: x if x.count('N') == 0 else None,seqs)
    n = nreads(counts)
    fp.write("Number of reads: %d\n" % n)
    fp.write("Number of unique barcode sequences: %d\n" % len(seqs))
    fp.write("Report up to %d sequences\n" % nseqs)
    if cutoff is not None:
        fp.write("Only report sequences appearing in >%.2f%% of reads\n"
                 % cutoff)
    if exclude_ns:
        fp.write("Don't include sequences that have N's\n")
    # Output top barcodes
    fp.write("Top barcode sequences\n")
    for i,seq in enumerate(seqs[:nseqs]):
        if cutoff is not None and percent_reads(counts[seq],n) < cutoff:
            # Stop reporting
            break
        fp.write("% 8d\t%s\t% 8d\t% 5.1f%%\n" % \
                 (i+1,seq,counts[seq],percent_reads(counts[seq],n)))

def get_barcodes_from_sample_sheet(sample_sheet,lanes,length=None):
    """Get a dictionary of barcode sequences from a sample sheet

    Arguments:
      sample_sheet: CasavaSampleSheet object
      lanes: list of lanes to get barcodes for
      length: (optional) if set then truncate barcodes to
        the specified length

    Returns:
      Dictionary with barcode sequences for keys, and sample names
      as the corresponding values.

    """
    # Filter by lane
    barcodes = dict()
    for line in sample_sheet:
        if line['Lane'] in lanes:
            # Get sample name
            sample_id = line['SampleID']
            # Fix dual index barcodes
            seq = ''.join(line['Index'].split('-'))
            # Truncate barcode
            if length is not None:
                seq = seq[:length]
            # Store sample name against index sequence
            barcodes[seq] = sample_id
    return barcodes

def match_barcodes(counts,samples,nseqs=20,fp=None,
                   max_mismatches=3,cutoff=None):
    """Find best matches for a set of supplied barcodes

    Arguments:
      counts: 'counts' dictionary returned by 'count_barcodes'
        nseqs: number of barcode sequences to report
      samples: dictionary with keys as barcodes and
        corresponding values as sample names
      nseqs: number of top barcode sequences to try and match
      max_mismatches: maximum number of mismatches to allow
        when reporting matches (will be corrected to allow
        for differing barcode lengths)
      cutoff: threshold fraction of reads below which not
        to report matches (default is to report everything)
      fp: specify output stream to write report to (default
        is to use stdout)

    """
    # Deal with cutoff
    if cutoff is not None:
        if cutoff < 1.0:
            cutoff = cutoff*100.0
    if fp is None: fp = sys.stdout
    fp.write("\nBest matches to supplied barcode sequences\n")
    fp.write("Report up to %d matches\n" % nseqs)
    if cutoff is not None:
        fp.write("Only report sequences appearing in >%.2f%% of reads\n"
                 % cutoff)
    fp.write("Allow up to %d mismatches\n" % max_mismatches)
    fp.write("NB all exact matches will be reported\n")
    # Get sequences in frequency order
    seqs = ordered_seqs(counts)[:nseqs]
    n = nreads(counts)
    # Find best matches for each sample barcode
    for barcode in samples:
        matches = {}
        for seq in seqs:
            matches[seq] = nmismatches(barcode,seq)
        # Sort into order of least to most mismatches
        ordered_matches = sorted(matches.keys(),
                                 cmp=lambda x,y: cmp(matches[x],matches[y]
                                                     if matches[x] != matches[y]
                                                     else cmp(counts[x],counts[y])))
        best_match = ordered_matches[0]
        # Maximum number of mismatches
        n_mismatches = max(max_mismatches + abs(len(barcode) - len(best_match)),
                           matches[best_match])
        # Report top matches i.e. lowest number of mismatches
        fp.write("%s\t%15s\t%s\t[%d]\t% 10d\t% 5.1f%%\n" % 
                 (samples[barcode],barcode,
                  best_match,
                  matches[best_match],
                  counts[best_match],
                  percent_reads(counts[best_match],n)))
        for seq in ordered_matches[1:]:
            if matches[seq] > 0 and \
               ((cutoff is not None and percent_reads(counts[seq],n) < cutoff) or \
               matches[seq] <= n_mismatches):
                # Stop reporting
                break
            fp.write("%s\t%15s\t%s\t[%d]\t% 10d\t% 5.1f%%\n" % 
                     ('','',
                      seq,
                      matches[seq],
                      counts[seq],
                      percent_reads(counts[seq],n)))

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

#######################################################################
# Main program
#######################################################################

if __name__ == '__main__':
    # Handle command line
    p = optparse.OptionParser(usage="\n\t%prog FASTQ [FASTQ...]\n\t%prog DIR\n\t%prog -c COUNTS_FILE",
                              version=__version__,
                              description="Report counts and statistics for Fastq index "
                              "sequences (aka barcodes). If from multiple Fastq files are "
                              "supplied then sequences will be pooled before being analysed. "
                              "If a single directory is supplied then this will be assumed to "
                              "be an output directory from CASAVA or bclToFastq and files "
                              "will be processed on a per-lane basis. If the -c option is "
                              "supplied then the input must be a file of barcode codes "
                              "generated previously using the -o option.")
    p.add_option('-c','--counts',action='store_true',dest='counts_file_in',
                 help="counts file generated by a previous run using -o")
    p.add_option('-n',action='store',dest='n',default=20,type='int',
                 help="report top N barcodes (default 20)")
    p.add_option('-o','--output',action='store',dest='counts_file',default=None,
                 help="output all counts to tab-delimited file COUNTS_FILE. This can be "
                 "used again in another run by specifying the '-c' option. NB if input "
                 "is a directory then this will be used as the base name for per-lane "
                 "count files")
    p.add_option('-l','--lanes',action='store',dest='lanes',default=None,
                 help="when processing a directory, specify one or more lane numbers to "
                 "report barcodes for (default is to process all lanes)")
    p.add_option('-r','--report',action='store',dest='report_file',default=None,
                 help="write report to REPORT_FILE (otherwise write to stdout)")
    p.add_option('-s','--sample-sheet',action='store',dest='sample_sheet',default=None,
                 help="report best matches against barcodes in SAMPLE_SHEET")
    p.add_option('-m','--mismatches',action='store',dest='mismatches',default=3,type='int',
                 help="maximum number of mismatches to use when reporting best matches "
                 "against barcodes from SAMPLE_SHEET (default is 3). Note that this is "
                 "corrected to allow for differences in sequence lengths.")
    p.add_option('-T','--truncate',action='store',dest='length',default=0,type='int',
                 help="truncate sample sheet sequences to LENGTH bases when reporting best "
                 "matches against barcodes from SAMPLE_SHEET.")
    p.add_option('-t','--threshold',action='store',dest='cutoff',default=None,type='float',
                 help="percentage or fraction of reads that a sequence must appear in "
                 "to be reported; sequences appearing fewer times will not be reported "
                 "(default is to report all barcodes).")
    p.add_option('-N','--nprocessors',action="store",dest="cores",default=1,type='int',
                 help="spread work across multiple processors/cores (default is 1)")
    options,args = p.parse_args()
    # Check arguments
    if not args:
        p.error("Need to supply at least one input Fastq file, a bclToFastq output "
                "directory, or a counts file from a previous run (if using -c)")
    if options.report_file is not None:
        print "Writing report to %s" % options.report_file
        fp = open(options.report_file,'w')
    else:
        fp = sys.stdout
    # Handle input sample sheet
    if options.sample_sheet is not None:
        print "Loading sample sheet data from %s" % options.sample_sheet
        sample_sheet = IlluminaData.get_casava_sample_sheet(options.sample_sheet)
    # Process according to inputs
    if options.counts_file_in:
        # Use counts from a previously generated file
        print "Loading counts from %s" % args[0]
        counts = dict()
        for line in open(args[0],'r'):
            seq = line.split('\t')[1]
            count = int(line.split('\t')[2])
            counts[seq] = count
        report(counts,nseqs=options.n,cutoff=options.cutoff,fp=fp)
        # Match barcodes to index sequences in sample sheet
        if options.sample_sheet:
            lanes = [int(lane) for lane in options.lanes.split(',')]
            barcodes = get_barcodes_from_sample_sheet(sample_sheet,
                                                      lanes=lanes,
                                                      length=options.length)
            match_barcodes(counts,barcodes,
                           nseqs=options.n,
                           max_mismatches=options.mismatches,
                           cutoff=options.cutoff,
                           fp=fp)
    elif len(args) == 1 and os.path.isdir(args[0]):
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
        print "Found lanes: %s" % lanes
        # Discard lanes not specified on the command line
        if options.lanes is not None:
            lanes = filter(lambda x: True if str(x) in options.lanes.split(',') else False,
                           lanes)
            print "Using lanes: %s" % lanes
        # Count barcodes in each lane
        for lane in lanes:
            fp.write("*** Counting barcodes in lane %d ***\n" % lane)
            for fq in fastq_in_lane[lane]:
                fp.write("Fastq: %s\n" % fq)
            counts = count_barcodes(fastq_in_lane[lane],n_processors=options.cores)
            report(counts,nseqs=options.n,cutoff=options.cutoff,fp=fp)
            if options.counts_file:
                counts_file = "%s.lane%s%s" % (os.path.splitext(options.counts_file)[0],
                                               lane,
                                               os.path.splitext(options.counts_file)[1])
                output(counts,counts_file)
            # Match barcodes to index sequences in sample sheet
            if options.sample_sheet:
                barcodes = get_barcodes_from_sample_sheet(sample_sheet,
                                                          lanes=lanes,
                                                          length=options.length)
                match_barcodes(counts,barcodes,
                               nseqs=options.n,
                               max_mismatches=options.mismatches,
                               cutoff=options.cutoff,
                               fp=fp)
    else:
        # One or more fastq files
        counts = count_barcodes(args)
        report(counts,nseqs=options.n,cutoff=options.cutoff,fp=fp)
        if options.counts_file:
            output(counts,options.counts_file)

