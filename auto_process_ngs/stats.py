#!/usr/bin/env python
#
#     stats.py: utilities for generating run-related statistics
#     Copyright (C) University of Manchester 2016-17 Peter Briggs
#
#########################################################################

"""
stats.py

Classes and functions for collecting and reporting statistics for
a run:

- FastqReadCounter: implements various methods for counting reads
  in FASTQ files

"""

#######################################################################
# Imports
#######################################################################

import sys
import subprocess
import bcftbx.TabFile as tf
import bcftbx.FASTQFile as FASTQFile
from bcftbx.IlluminaData import IlluminaFastq

#######################################################################
# Classes
#######################################################################

class FastqReadCounter:
    """
    Implements various methods for counting reads in FASTQ file

    The methods are:

    - simple: a wrapper for the FASTQFile.nreads() function
    - fastqiterator: counts reads using FASTQFile.FastqIterator
    - zcat_wc: runs 'zcat | wc -l' in the shell
    - reads_per_lane: counts reads by lane using FastqIterator

    """
    @staticmethod
    def simple(fastq=None,fp=None):
        """
        Return number of reads in a FASTQ file

        Uses the FASTQFile.nreads function to do the counting.

        Arguments:
          fastq: fastq(.gz) file
          fp: open file descriptor for fastq file

        Returns:
          Number of reads

        """
        return FASTQFile.nreads(fastq=fastq,fp=fp)
    @staticmethod
    def fastqiterator(fastq=None,fp=None):
        """
        Return number of reads in a FASTQ file

        Uses the FASTQFile.FastqIterator class to do the
        counting.

        Arguments:
          fastq: fastq(.gz) file
          fp: open file descriptor for fastq file

        Returns:
          Number of reads

        """
        nreads = 0
        for r in FASTQFile.FastqIterator(fastq_file=fastq,fp=fp):
            nreads += 1
        return nreads
    @staticmethod
    def zcat_wc(fastq=None,fp=None):
        """
        Return number of reads in a FASTQ file

        Uses a system call to run 'zcat FASTQ | wc -l' to do
        the counting (or just 'wc -l' if not a gzipped FASTQ).

        Note that this can only operate on fastq files (not
        on streams provided via the 'fp' argument; this will
        raise an exception).

        Arguments:
          fastq: fastq(.gz) file
          fp: open file descriptor for fastq file

        Returns:
          Number of reads

        """
        if fastq is None:
            raise Exception("zcat_wc: can only operate on a file")
        if fastq.endswith(".gz"):
            cmd = "zcat %s | wc -l" % fastq
        else:
            cmd = "wc -l %s | cut -d' ' -f1" % fastq
        output = subprocess.check_output(cmd,shell=True)
        try:
            return int(output)/4
        except Exception,ex:
            raise Exception("zcat_wc returned: %s" % output)
    @staticmethod
    def reads_per_lane(fastq=None,fp=None):
        """
        Return counts of reads in each lane of FASTQ file

        Uses the FASTQFile.FastqIterator class to do the
        counting, with counts split by lane.

        Arguments:
          fastq: fastq(.gz) file
          fp: open file descriptor for fastq file

        Returns:
          Dictionary where keys are lane numbers (as integers)
            and values are number of reads in that lane.

        """
        nreads = {}
        for r in FASTQFile.FastqIterator(fastq_file=fastq,fp=fp):
            lane = int(r.seqid.flowcell_lane)
            try:
                nreads[lane] += 1
            except KeyError:
                nreads[lane] = 1
        return nreads

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
