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
import bcftbx.IlluminaData as IlluminaData
import bcftbx.TabFile as TabFile
import bcftbx.utils as bcf_utils
from auto_process_ngs.stats import FastqStatistics
from auto_process_ngs.stats import FastqStats

from auto_process_ngs import get_version
__version__ = get_version()

# Initialise logging
import logging
logger = logging.getLogger("fastq_statistics")

# Format logging
logging.basicConfig(format='%(levelname)s %(name)s: %(message)s')

#######################################################################
# Main program
#######################################################################

if __name__ == '__main__':

    # Process command line
    p = optparse.OptionParser(usage="%prog [OPTIONS] ILLUMINA_RUN_DIR",
                              version="%prog "+__version__,
                              description="Generate statistics for FASTQ "
                              "files in ILLUMINA_RUN_DIR (top-level "
                              "directory of a processed Illumina run)")
    p.add_option("--unaligned",action="store",dest="unaligned_dir",
                 default="Unaligned",
                 help="specify an alternative name for the 'Unaligned' "
                 "directory containing the fastq.gz files")
    p.add_option('-o','--output',action="store",dest="stats_file",
                 default='statistics.info',
                 help="name of output file for per-file statistics (default "
                 "is 'statistics.info')")
    p.add_option('-p','--per-lane-stats',action="store",
                 dest="per_lane_stats_file",
                 default='per_lane_statistics.info',
                 help="name of output file for per-lane statistics (default "
                 "is 'per_lane_statistics.info')")
    p.add_option('-s','--per-lane-sample-stats',action="store",
                 dest="per_lane_sample_stats_file",
                 default='per_lane_sample_stats.info',
                 help="name of output file for per-lane statistics (default "
                 "is 'per_lane_sample_stats.info')")
    p.add_option('-f','--full-stats',action="store",
                 dest="full_stats_file",
                 default='statistics_full.info',
                 help="name of output file for full statistics (default "
                 "is 'statistics_full.info')")
    p.add_option('-u','--update',action="store_true",dest="update",
                 help="update existing full statistics file with stats for "
                 "additional files")
    p.add_option('-n',"--nprocessors",action="store",dest="n",
                 default=1,type='int',
                 help="spread work across N processors/cores (default is 1)")
    p.add_option("--debug",action="store_true",dest="debug",default=False,
                 help="turn on debugging output")
    # Deprecated options
    deprecated = optparse.OptionGroup(p,'Deprecated/defunct options')
    deprecated.add_option("--force",action="store_true",dest="force",
                          default=False,
                          help="does nothing: retained for backwards "
                          "compatibility")
    p.add_option_group(deprecated)
    options,args = p.parse_args()
    if len(args) != 1:
        p.error("expects a single argument (input directory)")

    # Report settings etc
    print "Source directory      : %s" % args[0]
    print "Unaligned subdirectory: %s" % options.unaligned_dir
    print "Basic stats file      : %s" % options.stats_file
    print "Full stats file       : %s" % options.full_stats_file
    print "Update existing stats?: %s" % ('yes' if options.update else 'no')
    print "Per-lane summary stats: %s" % options.per_lane_stats_file
    print "Per-lane sample stats : %s" % options.per_lane_sample_stats_file
    print "Number of processors  : %s" % options.n
    print "Debug?                : %s" % ('yes' if options.debug else 'no')

    # Check for existing stats
    if options.update:
        existing_stats_file = os.path.abspath(options.full_stats_file)
        if not os.path.exists(existing_stats_file):
            logging.fatal("No file '%s': cannot update" %
                          existing_stats_file)
            sys.exit(1)
    else:
        existing_stats_file = None

    # Ignore 'force'
    if options.force:
        logger.warn("ignoring deprecated option '--force'")

    # Handle debugging output if requested
    if options.debug:
        logging.getLogger("auto_process_ngs").setLevel(logging.DEBUG)

    # Get the data from FASTQ files
    try:
        illumina_data = IlluminaData.IlluminaData(
            args[0],
            unaligned_dir=options.unaligned_dir)
    except IlluminaData.IlluminaDataError,ex:
        logger.critical("failed to get data from %s: %s" % (args[0],ex))
        sys.exit(1)
    # Generate statistics for fastq files
    stats = FastqStatistics(illumina_data,
                            n_processors=options.n,
                            add_to=existing_stats_file)
    stats.report_full_stats(options.full_stats_file)
    print "Full statistics written to %s" % options.full_stats_file
    stats.report_basic_stats(options.stats_file)
    print "Basic statistics written to %s" % options.stats_file
    stats.report_per_lane_sample_stats(options.per_lane_sample_stats_file)
    print "Per-lane sample statistics written to %s" % \
        options.per_lane_sample_stats_file
    stats.report_per_lane_summary_stats(options.per_lane_stats_file)
    print "Per-lane summary statistics written to %s" % \
        options.per_lane_stats_file
