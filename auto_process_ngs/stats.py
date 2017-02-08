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

- FastqStatistics: collects and reports stats on FASTQs from an
  Illumina sequencing run
- FastqStats: container for storing data about a FASTQ file
- FastqReadCounter: implements various methods for counting reads
  in FASTQ files
- collect_fastq_data: collect data from FASTQ file in a FastqStats
  instance

"""

#######################################################################
# Imports
#######################################################################

import sys
import os
import time
import subprocess
from multiprocessing import Pool
import bcftbx.FASTQFile as FASTQFile
import bcftbx.utils as bcf_utils
from bcftbx.IlluminaData import IlluminaFastq
from bcftbx.TabFile import TabFile

# Initialise logging
import logging
logger = logging.getLogger(__name__)

#######################################################################
# Classes
#######################################################################

class FastqStatistics:
    """
    Class for collecting and reporting stats on Illumina FASTQs

    Given a directory with fastq(.gz) files arranged in the same
    structure as the output from bcl2fastq or bcl2fastq2,
    collects statistics for each file and provides methods for
    reporting different aspects.

    Example usage:

    >>> from IlluminaData import IlluminaData
    >>> data = IlluminaData('120117_BLAH_JSHJHXXX','bcl2fastq')
    >>> stats = FastqStatistics(data)
    >>> stats.report_basic_stats('basic_stats.out')

    """
    def __init__(self,illumina_data,n_processors=1):
        """
        Create a new FastqStatistics instance

        Arguments:
          illumina_data: populated IlluminaData object describing the
            run.
          n_processors: number of processors to use (if >1 then uses
            the multiprocessing library to run the statistics gathering
            using multiple cores).
        """
        self._illumina_data = illumina_data
        self._n_processors = n_processors
        self._stats = None
        self._lane_names = []
        self._get_data()

    def _get_data(self):
        """
        Collect statistics for FASTQ outputs from an Illumina run
        """
        # Collect FASTQ files
        fastqstats = []
        for project in self._illumina_data.projects:
            for sample in project.samples:
                for fastq in sample.fastq:
                    fastqstats.append(
                        FastqStats(os.path.join(sample.dirn,fastq),
                                   project.name,
                                   sample.name))
        # Gather same information for undetermined reads (if present)
        if self._illumina_data.undetermined is not None:
            for lane in self._illumina_data.undetermined.samples:
                for fastq in lane.fastq:
                    fastqstats.append(
                        FastqStats(os.path.join(lane.dirn,fastq),
                                   self._illumina_data.undetermined.name,
                                   lane.name))
        # Collect the data for each file
        if self._n_processors > 1:
            # Multiple cores
            pool = Pool(self._n_processors)
            results = pool.map(collect_fastq_data,fastqstats)
            pool.close()
            pool.join()
        else:
            # Single core
            results = map(collect_fastq_data,fastqstats)
        # Set up class to hold all collected data
        self._stats = TabFile(column_names=('Project',
                                            'Sample',
                                            'Fastq',
                                            'Size',
                                            'Nreads',
                                            'Paired_end',
                                            'Read_number'))
        # Split result sets into R1 and R2
        results_r1 = filter(lambda f: f.read_number == 1,results)
        results_r2 = filter(lambda f: f.read_number == 2,results)
        # Determine which lanes are present and append
        # columns for each
        lanes = set()
        for fastq in results_r1:
            logger.debug("-- %s: lanes %s" %
                         (fastq.name,
                          ','.join([str(l) for l in fastq.lanes])))
            for lane in fastq.lanes:
                lanes.add(lane)
        self._lanes = sorted(list(lanes))
        logger.debug("Lanes found: %s" %
                     ','.join([str(l) for l in self._lanes]))
        for lane in self._lanes:
            self._stats.appendColumn("L%s" % lane)
        # Copy reads per lane from R1 FASTQs into R2
        for r2_fastq in results_r2:
            # Get corresponding R1 name
            logger.debug("-- Fastq R2: %s" % r2_fastq.name)
            r1_fastq_name = IlluminaFastq(r2_fastq.name)
            r1_fastq_name.read_number = 1
            r1_fastq_name = str(r1_fastq_name)
            logger.debug("--    -> R1: %s" % r1_fastq_name)
            # Locate corresponding data
            r1_fastq = filter(lambda f: f.name.startswith(r1_fastq_name),
                              results_r1)[0]
            r2_fastq.reads_by_lane = dict(r1_fastq.reads_by_lane)
        # Write the data into the tabfile
        paired_end = ('Y' if self._illumina_data.paired_end else 'N')
        for fastq in results:
            data = [fastq.project,
                    fastq.sample,
                    fastq.name,
                    bcf_utils.format_file_size(fastq.fsize),
                    fastq.nreads,
                    paired_end,
                    fastq.read_number]
            for lane in lanes:
                try:
                    data.append(fastq.reads_by_lane[lane])
                except:
                    data.append('')
            self._stats.append(data=data)

    @property
    def lane_names(self):
        """
        Return list of lane names (e.g. ['L1','L2',...])
        """
        return [("L%d" % l) for l in self._lanes]

    @property
    def raw(self):
        """
        Return the 'raw' statistics TabFile instance
        """
        return self._stats

    def report_full_stats(self,out_file=None,fp=None):
        """
        Report all statistics gathered for all FASTQs

        Essentially a dump of all the data.

        Arguments:
          out_file (str): name of file to write report
            to (used if 'fp' is not supplied)
          fp (File): File-like object open for writing
            (defaults to stdout if 'out_file' also not
            supplied)
        """
        # Determine output stream
        if fp is None:
            if out_file is None:
                fpp = sys.stdout
            else:
                fpp = open(out_file,'w')
        else:
            fpp = fp
        # Report
        self._stats.write(fp=fpp,include_header=True)
        # Close file
        if fp is None and out_file is not None:
            fpp.close()

    def report_basic_stats(self,out_file=None,fp=None):
        """
        Report the 'basic' statistics

        For each FASTQ file, report the following information:

        - Project name
        - Sample name
        - FASTQ file name (without leading directory)
        - Size (human-readable)
        - Nreads (number of reads)
        - Paired_end ('Y' for paired-end, 'N' for single-end)

        Arguments:
          out_file (str): name of file to write report
            to (used if 'fp' is not supplied)
          fp (File): File-like object open for writing
            (defaults to stdout if 'out_file' also not
            supplied)
        """
        # Determine output stream
        if fp is None:
            if out_file is None:
                fpp = sys.stdout
            else:
                fpp = open(out_file,'w')
        else:
            fpp = fp
        # Report
        stats = TabFile(column_names=('Project',
                                      'Sample',
                                      'Fastq',
                                      'Size',
                                      'Nreads',
                                      'Paired_end'))
        for line in self._stats:
            data = [line[c] for c in stats.header()]
            stats.append(data=data)
        stats.write(fp=fpp,include_header=True)
        # Close file
        if fp is None and out_file is not None:
            fpp.close()

    def report_per_lane_sample_stats(self,out_file=None,fp=None):
        """
        Report of reads per sample in each lane

        Reports the number of reads for each sample in each
        lane plus the total reads for each lane.

        Example output:

        Lane 1
        Total reads = 182851745
        - KatyDobbs/KD-K1      79888058        43.7%
        - KatyDobbs/KD-K3      97854292        53.5%
        - Undetermined_indices/lane1       5109395 2.8%
        ...

        Arguments:
          out_file (str): name of file to write report
            to (used if 'fp' is not supplied)
          fp (File): File-like object open for writing
            (defaults to stdout if 'out_file' also not
            supplied)
        """
        # Determine output stream
        if fp is None:
            if out_file is None:
                fpp = sys.stdout
            else:
                fpp = open(out_file,'w')
        else:
            fpp = fp
        # Report
        lanes = self.lane_names
        for lane in lanes:
            lane_number = int(lane[1:])
            samples = filter(lambda x:
                             x['Read_number'] == 1 and bool(x[lane]),
                             self._stats)
            total_reads = sum([s[lane] for s in samples])
            fpp.write("\nLane %d\n" % lane_number)
            fpp.write("Total reads = %d\n" % total_reads)
            for sample in samples:
                sample_name = "%s/%s" % (sample['Project'],
                                         sample['Sample'])
                nreads = float(sample[lane])
                fpp.write("- %s\t%d\t%.1f%%\n" % (sample_name,
                                                  nreads,
                                                  nreads/total_reads*100.0))
        # Close file
        if fp is None and out_file is not None:
            fpp.close()

    def report_per_lane_summary_stats(self,out_file=None,fp=None):
        """
        Report summary of total and unassigned reads per-lane

        Arguments:
          out_file (str): name of file to write report
            to (used if 'fp' is not supplied)
          fp (File): File-like object open for writing
            (defaults to stdout if 'out_file' also not
            supplied)
        """
        # Determine output stream
        if fp is None:
            if out_file is None:
                fpp = sys.stdout
            else:
                fpp = open(out_file,'w')
        else:
            fpp = fp
        # Set up TabFile to hold the data collected
        per_lane_stats = TabFile(column_names=('Lane',
                                               'Total reads',
                                               'Assigned reads',
                                               'Unassigned reads',
                                               '%assigned',
                                               '%unassigned'))
        assigned = {}
        unassigned = {}
        # Count assigned and unassigned (= undetermined) reads
        for line in filter(lambda x: x['Read_number'] == 1,
                           self._stats):
            if line['Project'] != 'Undetermined_indices':
                counts = assigned
            else:
                counts = unassigned
            for lane in self.lane_names:
                if line[lane]:
                    try:
                        counts[lane] += line[lane]
                    except KeyError:
                        counts[lane] = line[lane]
        # Write out data for each lane
        for lane in self.lane_names:
            lane_number = int(lane[1:])
            assigned_reads = assigned[lane]
            try:
                unassigned_reads = unassigned[lane]
            except KeyError:
                # lane doesn't have any unassigned reads
                unassigned_reads = 0
            total_reads = assigned_reads + unassigned_reads
            percent_assigned = float(assigned_reads)/float(total_reads)*100.0
            percent_unassigned = float(unassigned_reads)/float(total_reads)*100.0
            per_lane_stats.append(data=("Lane %d" % lane_number,
                                        total_reads,
                                        assigned_reads,
                                        unassigned_reads,
                                        "%.2f" % percent_assigned,
                                        "%.2f" % percent_unassigned))
        # Write to file
        per_lane_stats.write(fp=fpp,include_header=True)
        # Close file
        if fp is None and out_file is not None:
            fpp.close()

class FastqStats:
    """Container for storing data about a FASTQ file

    This is a convenience wrapper for holding together data
    for a FASTQ file (full path, associated project and sample
    names, number of reads and filesize).
    """
    def __init__(self,fastq,project,sample):
        """
        Create a new FastqStats instance

        Arguments:
          fastq (str): full path to FASTQ file
          project (str): project name associated
            with FASTQ file
          sample (str): sample name associated
            with FASTQ file
        """
        self.fastq = fastq
        self.project = project
        self.sample = sample
        self.nreads = None
        self.fsize = None
        self.reads_by_lane = {}
    @property
    def name(self):
        """
        FASTQ file name without leading directory
        """
        return os.path.basename(self.fastq)
    @property
    def lanes(self):
        """
        Lane numbers associated with the FASTQ file
        """
        return sorted(self.reads_by_lane.keys())
    @property
    def read_number(self):
        """
        Read number extracted from the FASTQ name
        """
        return IlluminaFastq(self.name).read_number

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

def collect_fastq_data(fqstats):
    """
    Collect data from FASTQ file in a FastqStats instance

    Given a FastqStats instance, collects and sets the
    following properties derived from the corresponding
    FASTQ file stored in that instance:

    - nreads: total number of reads
    - fsize: file size
    - reads_by_lane: (R1 FASTQs only) dictionary
      where keys are lane numbers and values are
      read counts

    Note that if the FASTQ file is an R2 file then the
    reads per lane will not be set.

    Arguments:
      fqstats (FastqStats): FastqStats instance

    Returns:
      FastqStats: input FastqStats instance with the
        appropriated properties updated.
    """
    fqs = fqstats
    fastq = fqs.fastq
    fastq_name = fqs.name
    print "* %s: starting" % fastq_name
    start_time = time.time()
    sys.stdout.flush()
    if fqs.read_number == 1:
        # Do full processing for R1 fastqs
        lane = IlluminaFastq(fastq_name).lane_number
        if lane is not None:
            # Lane number is in file name
            fqs.reads_by_lane[lane] = \
                FastqReadCounter.zcat_wc(fastq)
        else:
            # Need to get lane(s) from read headers
            fqs.reads_by_lane = \
                FastqReadCounter.reads_per_lane(fastq)
        # Store total reads
        fqs.nreads = sum([fqs.reads_by_lane[x]
                          for x in fqs.lanes])
    else:
        # Only get total reads for R2 fastqs
        fqs.nreads = FastqReadCounter.zcat_wc(fastq)
    fqs.fsize = os.path.getsize(fastq)
    print "- %s: finished" % fastq_name
    end_time = time.time()
    print "- %s: %d reads, %s" % (fastq_name,
                                  fqs.nreads,
                                  bcf_utils.format_file_size(fqs.fsize))
    print "- %s: took %.2fs" % (fastq_name,(end_time-start_time))
    return fqs
