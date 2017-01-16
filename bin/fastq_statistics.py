#!/usr/bin/env python
#
#     fastq_statisticss.py: get stats for fastq files from Illumina run
#     Copyright (C) University of Manchester 2013-17 Peter Briggs
#
########################################################################
#
# fastq_statistics.py
#
#########################################################################

"""
fastq_statistics.py

Collect and report various statistics for FASTQ files produced from
an Illumina sequencing run:

- basic statistics: file size and number of reads for each FASTQ
- full statistics: file size, number of reads and reads per lane
  for each FASTQ
- per-lane summary statistics: total numbers of assigned and
  unassigned reads for each lane across all FASTQs
- per-lane, per-sample statistics: number of reads assigned to each
  sample in each lane

"""

#######################################################################
# Import modules that this module depends on
#######################################################################

import os
import sys
import optparse
import logging
import bcftbx.IlluminaData as IlluminaData
import bcftbx.TabFile as TabFile
import bcftbx.utils as bcf_utils
from auto_process_ngs.stats import FastqStatistics
from auto_process_ngs.stats import FastqStats

from auto_process_ngs import get_version
__version__ = get_version()

#######################################################################
# Main program
#######################################################################

if __name__ == '__main__':

    # Set up logger formatting
    logging.basicConfig(format='%(levelname)s %(message)s')

    # Process command line
    p = optparse.OptionParser(usage="%prog [OPTIONS] ILLUMINA_RUN_DIR",
                              version="%prog "+__version__,
                              description="Generate statistics for FASTQ "
                              "files in ILLUMINA_RUN_DIR (top-level "
                              "directory of a processed Illumina run)")
    p.add_option('-o','--output',action="store",dest="stats_file",
                 default='statistics.info',
                 help="name of output file for per-file statistics (default "
                 "is 'statistics.info')")
    p.add_option('-p','--per-lane-stats',action="store",
                 dest="per_lane_stats_file",default='per_lane_statistics.info',
                 help="name of output file for per-lane statistics (default "
                 "is 'per_lane_statistics.info')")
    p.add_option("--unaligned",action="store",dest="unaligned_dir",
                 default="Unaligned",
                 help="specify an alternative name for the 'Unaligned' "
                 "directory containing the fastq.gz files")
    p.add_option("--nprocessors",action="store",dest="n",
                 default=1,type='int',
                 help="spread work across N processors/cores (default is 1)")
    p.add_option("--force",action="store_true",dest="force",default=False,
                 help="force regeneration of statistics from fastq files")
    p.add_option("--debug",action="store_true",dest="debug",default=False,
                 help="turn on debugging output")
    options,args = p.parse_args()
    if len(args) != 1:
        p.error("expects a single argument (input directory)")

    # Report settings etc
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
    stats = FastqStatistics(illumina_data,
                            n_processors=options.n)
    stats.report_basic_stats(options.stats_file)
    stats.report_per_lane_sample_stats("per_lane_sample_stats.info")
    stats.report_per_lane_summary_stats(options.per_lane_stats_file)
    stats.report_full_stats("statistics_full.info")
    print "Basic statistics written to %s" % options.stats_file
    print "Per-lane summary statistics written to %s" % \
        options.per_lane_stats_file
