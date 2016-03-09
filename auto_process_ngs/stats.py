#!/usr/bin/env python
#
#     stats.py: utilities for generating run-related statistics
#     Copyright (C) University of Manchester 2016 Peter Briggs
#
#########################################################################

#######################################################################
# Imports
#######################################################################

import sys
import bcftbx.TabFile as tf
from bcftbx.IlluminaData import IlluminaFastq

#######################################################################
# Functions
#######################################################################

def report_per_lane_stats(stats_file,out_file=None):
    """
    Write per-lane statistics

    Analyse the 'statistics.info' format file from fastq_stats.py
    and summarise the number of reads per sample in each lane (and
    the percentage of reads that represents in the lane).

    Example output:

    Lane 1
    Total reads = 182851745
    - KatharineDibb/KD-K1	79888058	43.7%
    - KatharineDibb/KD-K3	97854292	53.5%
    - Undetermined_indices/lane1	5109395	2.8%
    ...

    Arguments:
      stats_file (str): path of input file with per-file
        statistics
      out_file (str): optional path to output file; if None then
        statistics will be written to stdout.

    """
    # Read in data
    stats = tf.TabFile(stats_file,first_line_is_header=True)
    data = dict()
    for line in stats:
        # Collect sample name, lane etc
        fq = IlluminaFastq(line['Fastq'])
        if fq.read_number != 1:
            # Only interested in R1 reads
            continue
        lane = fq.lane_number
        sample = "%s/%s" % (line['Project'],line['Sample'])
        nreads = line['Nreads']
        # Update information in dict
        if lane not in data:
            data[lane] = dict()
        try:
            data[lane][sample] += nreads
        except KeyError:
            data[lane][sample] = nreads
    # Get list of lanes
    lanes = [int(x) for x in data.keys()]
    lanes.sort()
    # Determine output stream
    if out_file is None:
        fp = sys.stdout
    else:
        fp = open(out_file,'w')
    # Report
    for lane in lanes:
        fp.write("\nLane %d\n" % lane)
        samples = data[lane].keys()
        samples.sort()
        total_reads = sum([data[lane][x] for x in samples])
        fp.write("Total reads = %d\n" % total_reads)
        for sample in samples:
            nreads = float(data[lane][sample])
            if total_reads > 0:
                fp.write("- %s\t%d\t%.1f%%\n" % (sample,nreads,
                                                 nreads/total_reads*100.0))
            else:
                fp.write("- %s\t%d\tn/a\n" % (sample,nreads))
    # Close file
    if out_file is not None:
        fp.close()
