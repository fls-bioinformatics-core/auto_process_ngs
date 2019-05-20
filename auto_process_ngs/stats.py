#!/usr/bin/env python
#
#     stats.py: utilities for generating run-related statistics
#     Copyright (C) University of Manchester 2016-19 Peter Briggs
#
#########################################################################

"""
stats.py

Classes and functions for collecting and reporting statistics for
a run:

- FastqStatistics: collects and reports stats on FASTQs from an
  Illumina sequencing run
- FastqStats: container for storing data about a FASTQ file
- collect_fastq_data: collect data from FASTQ file in a FastqStats
  instance

"""

#######################################################################
# Imports
#######################################################################

import sys
import os
import time
from multiprocessing import Pool
import bcftbx.utils as bcf_utils
from bcftbx.IlluminaData import IlluminaFastq
from bcftbx.IlluminaData import SampleSheet
from bcftbx.TabFile import TabFile
from bcftbx.TabFile import TabDataLine
from .fastq_utils import FastqReadCounter

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
    def __init__(self,illumina_data,n_processors=1,add_to=None):
        """
        Create a new FastqStatistics instance

        Arguments:
          illumina_data: populated IlluminaData object describing the
            run.
          n_processors: number of processors to use (if >1 then uses
            the multiprocessing library to run the statistics gathering
            using multiple cores).
          add_to: optional, add the data to that from an existing
            statistics file
        """
        self._illumina_data = illumina_data
        self._n_processors = n_processors
        self._stats = None
        self._lane_names = []
        self._get_data(filen=add_to)

    def _get_data(self,filen=None):
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
        # Set up tabfile to hold pre-existing data
        if filen is not None:
            existing_stats = TabFile(filen,first_line_is_header=True)
        else:
            existing_stats = None
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
        # Add lane numbers from pre-existing stats file
        if existing_stats is not None:
            for c in existing_stats.header():
                if c.startswith('L'):
                    lanes.add(int(c[1:]))
        self._lanes = sorted(list(lanes))
        logger.debug("Lanes found: %s" %
                     ','.join([str(l) for l in self._lanes]))
        for lane in self._lanes:
            self._stats.appendColumn("L%s" % lane)
        # Copy pre-existing stats into new tabfile
        if existing_stats:
            for line in existing_stats:
                data = [line['Project'],
                        line['Sample'],
                        line['Fastq'],
                        line['Size'],
                        line['Nreads'],
                        line['Paired_end'],
                        line['Read_number']]
                for lane in lanes:
                    try:
                        data.append(line["L%s" % lane])
                    except:
                        data.append('')
                self._stats.append(data=data)
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
            # Check for existing entry
            existing_entry = False
            for line in self._stats:
                if (line['Project'] == fastq.project and
                    line['Sample'] == fastq.sample and
                    line['Fastq'] == fastq.name):
                    # Overwrite the existing entry
                    existing_entry = True
                    break
            # Write the data
            if not existing_entry:
                # Append new entry
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
            else:
                # Overwrite existing entry
                logging.warning("Overwriting exisiting entry for "
                                "%s/%s/%s" % (fastq.project,
                                              fastq.sample,
                                              fastq.name))
                line['Size'] = bcf_utils.format_file_size(fastq.fsize)
                line['Nreads'] = fastq.nreads
                line['Paired_end'] = paired_end
                line['Read_number'] = fastq.read_number
                for lane in lanes:
                    lane_name = "L%d" % lane
                    try:
                        line[lane_name] = fastq.reads_by_lane[lane]
                    except:
                        line[lane_name] = ''

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

    def report_per_lane_sample_stats(self,out_file=None,fp=None,
                                     samplesheet=None):
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
          samplesheet (str): optional sample sheet file
            to get additional data from
        """
        # Determine output stream
        if fp is None:
            if out_file is None:
                fpp = sys.stdout
            else:
                fpp = open(out_file,'w')
        else:
            fpp = fp
        # Get data from samplesheet
        expected_samples = {}
        if samplesheet:
            s = SampleSheet(samplesheet)
            ncol = s.sample_id_column
            pcol = s.sample_project_column
            for data in s:
                if s.has_lanes:
                    lanes = ['L%d' % data['Lane']]
                else:
                    lanes = self.lane_names
                sample = {
                    'Project': data[pcol],
                    'Sample': data[ncol],
                }
                for lane in lanes:
                    try:
                        expected_samples[lane].append(sample)
                    except KeyError:
                        expected_samples[lane] = [sample,]
        # Report
        lanes = self.lane_names
        for lane in lanes:
            lane_number = int(lane[1:])
            samples = filter(lambda x:
                             x['Read_number'] == 1 and bool(x[lane]),
                             self._stats)
            # Additional samples from samplesheet
            if lane in expected_samples:
                for sample in expected_samples[lane]:
                    found_sample = False
                    for smpl in samples:
                        if smpl['Sample'] == sample['Sample'] and \
                           smpl['Project'] == sample['Project']:
                            found_sample = True
                            break
                    if not found_sample:
                        # Add the expected sample with zero reads
                        # for the lane being examined
                        samples.append(
                            TabDataLine(
                                line="%s\t%s\t0" % (sample['Project'],
                                                    sample['Sample']),
                                column_names=('Project','Sample',lane)))
                # Sort into order
                samples = sorted(samples,
                                 key=lambda x: (x['Project'],x['Sample']))
            try:
                total_reads = sum([int(s[lane]) for s in samples])
            except Exception as ex:
                for s in samples:
                    try:
                        int(s[lane])
                    except ValueError:
                        logging.critical("Bad value for read count in "
                                         "lane %s sample %s: '%s'" %
                                         (lane,s['Sample'],s[lane]))
                raise ex
            fpp.write("\nLane %d\n" % lane_number)
            fpp.write("Total reads = %d\n" % total_reads)
            for sample in samples:
                sample_name = "%s/%s" % (sample['Project'],
                                         sample['Sample'])
                nreads = float(sample[lane])
                if total_reads > 0:
                    frac_reads = "%.1f%%" % (nreads/total_reads*100.0)
                else:
                    frac_reads = "n/a"
                fpp.write("- %s\t%d\t%s\n" % (sample_name,
                                              nreads,
                                              frac_reads))
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
        # Initialise counts for each lane
        assigned = {}
        unassigned = {}
        for lane in self.lane_names:
            assigned[lane] = 0
            unassigned[lane] = 0
        # Count assigned and unassigned (= undetermined) reads
        for line in filter(lambda x: x['Read_number'] == 1 and
                           not IlluminaFastq(x['Fastq']).is_index_read,
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
            if total_reads > 0:
                percent_assigned = float(assigned_reads)/ \
                                   float(total_reads)*100.0
                percent_unassigned = float(unassigned_reads)/ \
                                     float(total_reads)*100.0
            else:
                percent_assigned = 0.0
                percent_unassigned = 0.0
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
