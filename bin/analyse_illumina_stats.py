#!/usr/bin/env python
#
#     analyse_illumina_stats.py: summarise per lane/per sample stats
#     Copyright (C) University of Manchester 2015 Peter Briggs
#
#########################################################################
#
# analyse_illumina_stats.py
#
#########################################################################

"""analyse_illumina_stats

Utility to analyse the 'statistics.info' file from auto_process.py
and summarise the number of reads per sample in each lane (and the
percentage of reads that represents in the lane).

Example output:

Lane 1
Total reads = 182851745
- KatharineDibb/KD-K1	79888058	43.7%
- KatharineDibb/KD-K3	97854292	53.5%
- Undetermined_indices/lane1	5109395	2.8%
...

"""

#######################################################################
# Module metadata
#######################################################################

__version__ = "0.0.1"

#######################################################################
# Imports
#######################################################################

import optparse
from auto_process_ngs import stats

#######################################################################
# Main script
#######################################################################

if __name__ == '__main__':

    # Process command line
    p = optparse.OptionParser(usage="%prog [OPTIONS] STATS_INFO_FILE",
                              version="%prog "+__version__,
                              description="Summarise the per sample stats for each lane "
                              "for an Illumina sequencing run, using the data from "
                              "STATS_INFO_FILE (typically called 'statistics.info' and "
                              "found in the top-level directory of sequencing runs "
                              "processed using auto_process.py.")
    options,args = p.parse_args()
    if len(args) == 0:
        stats_file = "statistics.info"
    elif len(args) > 1:
        p.error("expects a single argument (statistics file)")
    else:
        stats_file = args[0]
    stats.report_per_lane_stats(stats_file)
