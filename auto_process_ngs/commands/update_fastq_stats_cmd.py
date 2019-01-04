#!/usr/bin/env python
#
#     update_fastq_stats_cmd.py: implement update_fastq_stats command
#     Copyright (C) University of Manchester 2019 Peter Briggs
#
#########################################################################

#######################################################################
# Imports
#######################################################################

import logging
from .make_fastqs_cmd import fastq_statistics

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Command functions
#######################################################################

def update_fastq_stats(ap,stats_file=None,per_lane_stats_file=None,
                       unaligned_dir=None,add_data=False,nprocessors=None,
                       runner=None):
    """Update statistics for Fastq files

    Updates the statistics for all Fastq files found in the
    'unaligned' directory, by running the 'fastq_statistics.py'
    program.

    Arguments
      ap (AutoProcessor): autoprocessor pointing to the analysis
        directory to create Fastqs for
      stats_file (str): path of a non-default file to write the
        statistics to (defaults to 'statistics.info' unless
        over-ridden by local settings)
      per_lane_stats_file (str): path for per-lane statistics
        output file (defaults to 'per_lane_statistics.info'
        unless over-ridden by local settings)
      unaligned_dir (str): output directory for bcl-to-fastq
        conversion
      add_data (bool): if True then add stats to the existing
        stats files (default is to overwrite existing stats
        files)
      nprocessors (int): number of cores to use when running
        'fastq_statistics.py'
      runner (JobRunner): (optional) specify a non-default job
        runner to use for running 'fastq_statistics.py'
    """
    ap.set_log_dir(ap.get_log_subdir("update_fastq_stats"))
    fastq_statistics(unaligned_dir=options.unaligned_dir,
                     stats_file=options.stats_file,
                     per_lane_stats_file=options.per_lane_stats_file,
                     add_data=options.add_data,
                     nprocessors=options.nprocessors,
                     runner=options.runner)
