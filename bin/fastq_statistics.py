#!/usr/bin/env python
#
#     fastq_stats.py: get statistics for fastq files from Illumina run
#     Copyright (C) University of Manchester 2013-17 Peter Briggs
#
########################################################################
#
# fastq_stats.py
#
#########################################################################

"""fastq_stats.py

Generate statistics for fastq files for an Illumina sequencing run.

"""

#######################################################################
# Import modules that this module depends on
#######################################################################

import os
import sys
import optparse
import logging
import time
import subprocess

# Put .. onto Python search path for modules
SHARE_DIR = os.path.abspath(
    os.path.normpath(
        os.path.join(os.path.dirname(sys.argv[0]),'..')))
sys.path.append(SHARE_DIR)
import bcftbx.IlluminaData as IlluminaData
import bcftbx.TabFile as TabFile
import bcftbx.FASTQFile as FASTQFile 
import bcftbx.utils as bcf_utils
import auto_process_ngs.stats as auto_process_stats
from multiprocessing import Pool

from auto_process_ngs import get_version
__version__ = get_version()

#######################################################################
# Functions
#######################################################################

class FastqStats:
    """Container for storing data for a Fastq file

    This is a convenience wrapper for holding together data
    for a fastq file (full path, associated project and sample
    names, number of reads and filesize), for use with the
    'get_stats_for_file' and 'map' functions.

    """
    def __init__(self,fastq,project,sample):
        self.fastq = fastq
        self.project = project
        self.sample = sample
        self.nreads = None
        self.fsize = None
        self.reads_by_lane = {}
    @property
    def name(self):
        return os.path.basename(self.fastq)
    @property
    def lanes(self):
        return sorted(self.reads_by_lane.keys())
    @property
    def read_number(self):
        return IlluminaData.IlluminaFastq(self.name).read_number

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
        return FASTQFile.nreads
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
        the counting.

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
        output = subprocess.check_output(("zcat %s | wc -l" % fastq),
                                         shell=True)
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

def get_fastqs(illumina_data):
    """Get a list of fastq files from an Illumina run

    Arguments:
      illumina_data: populated IlluminaData object describing the
        run.

    Returns:
      List of populated FastqStats objects.

    """
    fastqs = []
    for project in illumina_data.projects:
        for sample in project.samples:
            for fastq in sample.fastq:
                fastqs.append(FastqStats(os.path.join(sample.dirn,fastq),
                                         project.name,
                                         sample.name))
    # Gather same information for undetermined reads (if present)
    if illumina_data.undetermined is not None:
        for lane in illumina_data.undetermined.samples:
            for fastq in lane.fastq:
                fastqs.append(FastqStats(os.path.join(lane.dirn,fastq),
                                         illumina_data.undetermined.name,
                                         lane.name))
    return fastqs

def get_stats_for_file(fq,read_counter=FastqReadCounter.zcat_wc):
    """Generate statistics for a single fastq file

    Given a FastqStats object, set the 'nreads' property to
    the number of reads and the 'fsize' property to the file
    size for the corresponding fastq file.

    Arguments:
      fq: FastqStats object with 'fastq' property set to the
        full path for a Fastq file
      read_counter: optional, specify function to use for
        counting reads in the fastq file

    Returns:
      Input FastqStats object with the 'nreads' and 'fsize'
      properties set.

    """
    print "* %s: starting" % fq.name
    start_time = time.time()
    sys.stdout.flush()
    fq.nreads = read_counter(fq.fastq)
    fq.fsize = os.path.getsize(fq.fastq)
    print "- %s: finished" % fq.name
    end_time = time.time()
    print "- %s: %d reads, %s" % (fq.name,
                                  fq.nreads,
                                  bcf_utils.format_file_size(fq.fsize))
    print "- %s: took %f.2s" % (fq.name,(end_time-start_time))
    return fq

def get_per_lane_stats_for_file(fq):
    """
    Collect statistics for a single fastq file

    Given a FastqStats object, sets the following properties
    for the corresponding FASTQ file:

    - nreads: total number of reads
    - fsize: file size
    - reads_by_lane: (R1 FASTQs only) dictionary where keys
      are lane numbers and values are read counts

    Arguments:
      fq: FastqStats object with 'fastq' property set to the
        full path for a Fastq file

    Returns:
      Input FastqStats object with the 'nreads', 'fsize'
      and 'reads_by_lane' (if R1) properties set.

    """
    print "* %s: starting" % fq.name
    start_time = time.time()
    sys.stdout.flush()
    if fq.read_number == 1:
        # Do full processing for R1 fastqs
        lane = IlluminaData.IlluminaFastq(fq.name).lane_number
        if lane is not None:
            # Lane number is in file name
            fq.reads_by_lane[lane] = FastqReadCounter.zcat_wc(fq.fastq)
        else:
            # Need to get lane(s) from read headers
            fq.reads_by_lane = FastqReadCounter.reads_per_lane(fq.fastq)
        # Store total reads
        fq.nreads = sum([fq.reads_by_lane[x] for x in fq.lanes])
    else:
        # Only get total reads for R2 fastqs
        fq.nreads = FastqReadCounter.zcat_wc(fq.fastq)
    fq.fsize = os.path.getsize(fq.fastq)
    print "- %s: finished" % fq.name
    end_time = time.time()
    print "- %s: %d reads, %s" % (fq.name,
                                  fq.nreads,
                                  bcf_utils.format_file_size(fq.fsize))
    print "- %s: took %f.2s" % (fq.name,(end_time-start_time))
    return fq

def fastq_statistics(illumina_data,n_processors=1):
    """Generate statistics for fastq outputs from an Illumina run

    Given a directory with fastq(.gz) files arranged in the same
    structure as the output from bcl2fastq or bcl2fastq2,
    generate statistics for each file.

    Arguments:
      illumina_data: populated IlluminaData object describing the
        run.
      n_processors: number of processors to use (if 1>1 then uses
        the multiprocessing library to run the statistics gathering
        using multiple cores).

    Returns:
      Populated TabFile object containing the statistics.

    """
    fastqs = get_fastqs(illumina_data)
    if n_processors > 1:
        # Multiple cores
        pool = Pool(n_processors)
        results = pool.map(get_per_lane_stats_for_file,fastqs)
        pool.close()
        pool.join()
    else:
        # Single core
        results = map(get_per_lane_stats_for_file,fastqs)
    # Per-file stats
    stats = TabFile.TabFile(column_names=('Project',
                                          'Sample',
                                          'Fastq',
                                          'Size',
                                          'Nreads',
                                          'Paired_end'))
    for fastq in results:
        stats.append(data=(fastq.project,
                           fastq.sample,
                           fastq.name,
                           bcf_utils.format_file_size(fastq.fsize),
                           fastq.nreads,
                           'Y' if illumina_data.paired_end else 'N'))
    # Per-lane stats
    stats_per_lane = TabFile.TabFile(column_names=('Project',
                                                   'Sample',
                                                   'Fastq',))
    results_r1 = filter(lambda f: f.read_number == 1,results)
    # Determine which lanes are present and append
    # columns for each
    lanes = set()
    for fastq in results_r1:
        print "%s: %s" % (fastq.name,fastq.lanes)
        for lane in fastq.lanes:
            lanes.add(lane)
    lanes = sorted(list(lanes))
    print "%s" % lanes
    for lane in lanes:
        stats_per_lane.appendColumn("L%s" % lane)
    # Write the number of reads in each lane
    for fastq in results_r1:
        data = [fastq.project,
                fastq.sample,
                fastq.name,]
        for lane in lanes:
            try:
                data.append(fastq.reads_by_lane[lane])
            except:
                data.append('')
        stats_per_lane.append(data=data)
    # Return the data
    return (stats,stats_per_lane)

def sequencer_stats(stats_per_lane):
    """Generate per-lane statistics for an Illumina run

    Given a populated TabFile object containing data about each
    fastq(.gz) file from an Illumina run (as generated by the
    'fastq_statistics' function), calculate numbers of assigned,
    unassigned and total reads for each lane.

    Arguments:
      stats: populated TabFile object with per file data from
        fastq_statistics.

    Returns:
      Populated TabFile object containing the per-lane statistics.

    """
    # Set up TabFile to hold the data collected
    sequencer_stats = TabFile.TabFile(column_names=('Lane',
                                                    'Total reads',
                                                    'Assigned reads',
                                                    'Unassigned reads'))
    # Figure out lanes
    lanes = stats_per_lane.header()[3:]
    assigned = {}
    unassigned = {}
    # Count assigned and unassigned (= undetermined) reads
    for line in stats_per_lane:
        if line['Project'] != 'Undetermined_indices':
            for lane in lanes:
                if line[lane]:
                    try:
                        assigned[lane] += line[lane]
                    except KeyError:
                        assigned[lane] = line[lane]
        else:
            for lane in lanes:
                if line[lane]:
                    try:
                        unassigned[lane] += line[lane]
                    except KeyError:
                        unassigned[lane] = line[lane]
    # Write out data for each lane
    for lane in lanes:
        lane_number = int(lane[1:])
        assigned_reads = assigned[lane]
        unassigned_reads = unassigned[lane]
        total_reads = assigned_reads + unassigned_reads
        sequencer_stats.append(data=("Lane %d" % lane_number,
                                     total_reads,
                                     assigned_reads,
                                     unassigned_reads))
    # Return the data
    return sequencer_stats

def report_sample_stats(stats_per_lane,out_file):
    """
    Report of reads per sample in each lane for an Illumina run

    Given a populated TabFile object containing per-lane data
    about each fastq(.gz) file from an Illumina run (as generated
    by the 'fastq_statistics' function), report the number of reads
    for each sample in each lane plus the total reads for each
    lane.

    Arguments:
      stats_per_lane: populated TabFile object with per-lane data
        for each FASTQ from fastq_statistics.
    """
    # Determine output stream
    if out_file is None:
        fp = sys.stdout
    else:
        fp = open(out_file,'w')
    # Report
    lanes = stats_per_lane.header()[3:]
    for lane in lanes:
        lane_number = int(lane[1:])
        samples = filter(lambda x: bool(x[lane]),stats_per_lane)
        total_reads = sum([s[lane] for s in samples])
        fp.write("\nLane %d\n" % lane_number)
        fp.write("Total reads = %d\n" % total_reads)
        for sample in samples:
            sample_name = "%s/%s" % (sample['Project'],
                                     sample['Sample'])
            nreads = float(sample[lane])
            fp.write("- %s\t%d\t%.1f%%\n" % (sample_name,
                                             nreads,
                                             nreads/total_reads*100.0))
    # Close file
    if out_file is not None:
        fp.close()

#######################################################################
# Main program
#######################################################################

if __name__ == '__main__':

    # Set up logger formatting
    logging.basicConfig(format='%(levelname)s %(message)s')

    # Process command line
    p = optparse.OptionParser(usage="%prog [OPTIONS] ILLUMINA_RUN_DIR",
                              version="%prog "+__version__,
                              description="Generate statistics for fastq files in "
                              "ILLUMINA_RUN_DIR (top-level directory of the Illumina "
                              "run to be processed, plus assigned and unassigned reads "
                              "per lane)")
    p.add_option('-o','--output',action="store",dest="stats_file",default='statistics.info',
                 help="name of output file for per-file statistics (default "
                 "is 'statistics.info')")
    p.add_option('-p','--per-lane-stats',action="store",
                 dest="per_lane_stats_file",default='per_lane_statistics.info',
                 help="name of output file for per-lane statistics (default "
                 "is 'per_lane_statistics.info')")
    p.add_option("--unaligned",action="store",dest="unaligned_dir",default="Unaligned",
                 help="specify an alternative name for the 'Unaligned' directory "
                 "containing the fastq.gz files")
    p.add_option("--nprocessors",action="store",dest="n",default=1,type='int',
                 help="spread work across N processors/cores (default is 1)")
    p.add_option("--force",action="store_true",dest="force",default=False,
                 help="force regeneration of statistics from fastq files")
    p.add_option("--debug",action="store_true",dest="debug",default=False,
                 help="turn on debugging output")
    options,args = p.parse_args()
    if len(args) != 1:
        p.error("expects a single argument (input directory)")

    # Report settings etc
    print "%s version %s" % (os.path.basename(sys.argv[0]),__version__)
    print "Source dir    : %s" % args[0]
    print "Unaligned dir : %s" % options.unaligned_dir
    print "Stats file    : %s" % options.stats_file
    print "Per-lane stats: %s" % options.per_lane_stats_file
    print "Nprocessors   : %s" % options.n
    print "Force?        : %s" % options.force
    print "Debug?        : %s" % options.debug

    # Handle debugging output if requested
    if options.debug:
        logging.getLogger().setLevel(logging.DEBUG)

    # Get the data from FASTQ files
    try:
        illumina_data = IlluminaData.IlluminaData(
            args[0],
            unaligned_dir=options.unaligned_dir)
    except IlluminaData.IlluminaDataError,ex:
        logging.error("Failed to get data from %s: %s" % (args[0],ex))
        sys.exit(1)
    # Generate statistics for fastq files
    stats,per_lane = fastq_statistics(illumina_data,
                                      n_processors=options.n)
    stats.write(options.stats_file,include_header=True)
    print "Statistics written to %s" % options.stats_file
    per_lane.write(options.per_lane_stats_file,include_header=True)
    print "Per-lane sequencer stats written to %s" % \
        options.per_lane_stats_file
    seqstats = sequencer_stats(per_lane)
    seq_stats_file="per_lane_statistics.summary"
    seqstats.write(seq_stats_file,include_header=True)
    print "Summary of per-lane sequencer stats written to %s" % \
        seq_stats_file
    report_sample_stats(per_lane,"sample_stats.info")
    sys.exit()
