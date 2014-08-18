#!/bin/env python
#
#     fastq_stats.py: get statistics for fastq files from Illumina run
#     Copyright (C) University of Manchester 2013 Peter Briggs
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
# Module metadata
#######################################################################

__version__ = "0.1.0"

#######################################################################
# Import modules that this module depends on
#######################################################################

import os
import sys
import optparse
import logging

# Put ../share onto Python search path for modules
SHARE_DIR = os.path.abspath(
    os.path.normpath(
        os.path.join(os.path.dirname(sys.argv[0]),'..','share')))
sys.path.append(SHARE_DIR)
import IlluminaData
import TabFile
import FASTQFile
import bcf_utils
from multiprocessing import Pool

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
    @property
    def name(self):
        return os.path.basename(self.fastq)

#######################################################################
# Functions
#######################################################################

def get_fastqs(illumina_data):
    """Get a list of fastq files from an Illumina run

    Given a directory with fastq(.gz) files arranged in the same
    structure as the output from bcl2fastq (i.e. subdirectory
    'Unaligned', then project directories within this called
    'Project_<NAME>', each containing sample directories called
    'Sample_<NAME>', and each of these containing fastq files),
    return a list of FastqStats objects representing each of the
    fastq files.

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

def get_stats_for_file(fq):
    """Generate statistics for a single fastq file

    Given a FastqStats object, set the 'nreads' property to
    the number of reads and the 'fsize' property to the file
    size for the corresponding fastq file.

    Arguments:
      fq: FastqStats object with 'fastq' property set to the
        full path for a Fastq file

    Returns:
      Input FastqStats object with the 'nreads' and 'fsize'
      properties set.

    """
    print "* %s" % fq.name
    fq.nreads = FASTQFile.nreads(fq.fastq)
    fq.fsize = os.path.getsize(fq.fastq)
    return fq

def fastq_statistics(illumina_data,n_processors=1):
    """Generate statistics for fastq outputs from an Illumina run

    Given a directory with fastq(.gz) files arranged in the same
    structure as the output from bcl2fastq (i.e. subdirectory
    'Unaligned', then project directories within this called
    'Project_<NAME>', each containing sample directories called
    'Sample_<NAME>', and each of these containing fastq files),
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
    stats = TabFile.TabFile(column_names=('Project',
                                          'Sample',
                                          'Fastq',
                                          'Size',
                                          'Nreads',
                                          'Paired_end'))
    fastqs = get_fastqs(illumina_data)
    if n_processors > 1:
        # Multiple cores
        pool = Pool(n_processors)
        results = pool.map(get_stats_for_file,fastqs)
        pool.close()
        pool.join()
    else:
        # Single core
        results = map(get_stats_for_file,fastqs)
    for fastq in results:
        stats.append(data=(fastq.project,
                           fastq.sample,
                           fastq.name,
                           bcf_utils.format_file_size(fastq.fsize),
                           fastq.nreads,
                           'Y' if illumina_data.paired_end else 'N'))
    return stats

def sequencer_stats(stats):
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
    # Get the data for each lane
    for lane in xrange(1,9):
        data = []
        pattern = "_L00%d_R1_001.fastq.gz" % lane
        total_reads = 0
        assigned_reads = 0
        unassigned_reads = 0
        for line in stats:
            if line['Fastq'].endswith(pattern):
                data.append(line)
                nreads = line['Nreads']
                total_reads += nreads
                if line['Project'] != 'Undetermined_indices':
                    assigned_reads += nreads
                else:
                    unassigned_reads += nreads
        if not data:
            continue
        sequencer_stats.append(data=("Lane %d" % lane,
                                     total_reads,
                                     assigned_reads,
                                     unassigned_reads))
    return sequencer_stats

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
                 help="name of output file (default is 'statistics.info')")
    p.add_option("--unaligned",action="store",dest="unaligned_dir",default="Unaligned",
                 help="specify an alternative name for the 'Unaligned' directory "
                 "containing the fastq.gz files")
    p.add_option("--nprocessors",action="store",dest="n",default=1,type='int',
                 help="spread work across N processors/cores (default is 1)")
    p.add_option("--force",action="store_true",dest="force",default=False,
                 help="force regeneration of statistics from fastq files")
    options,args = p.parse_args()
    if len(args) != 1:
        p.error("expects a single argument (input directory)")

    # Get the raw data
    stats = None
    filein = os.path.join(args[0],'statistics.info')
    if os.path.isfile(filein):
        # Existing statistics file detected, read data from here
        logging.warning("Existing statistics file detected")
        if not options.force:
            stats = TabFile.TabFile(filen=filein,first_line_is_header=True)
        else:
            logging.warning("--force specified, statistics will be regenerated")
    if stats is None:
        # (Re)create statistics from raw files
        try:
            illumina_data = IlluminaData.IlluminaData(args[0],
                                                      unaligned_dir=options.unaligned_dir)
        except IlluminaData.IlluminaDataError,ex:
            logging.error("Failed to get data from %s: %s" % (args[0],ex))
            sys.exit(1)
        # Generate statistics for fastq files
        stats = fastq_statistics(illumina_data,n_processors=options.n)
        stats.write(options.stats_file,include_header=True)
    print "Statistics written to %s" % options.stats_file
    # Per-lane sequencer statistics
    seq_stats_file = 'per_lane_stats.info'
    seqstats = sequencer_stats(stats)
    seqstats.write(seq_stats_file,include_header=True)
    print "Per-lane sequencer stats written to %s" % seq_stats_file
