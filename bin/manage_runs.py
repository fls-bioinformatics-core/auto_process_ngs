#!/bin/env python
#
#     manage_runs.py: utility for managing Illumina sequencing run directories
#     Copyright (C) University of Manchester 2014 Peter Briggs
#
#########################################################################
#
# manage_runs.py
#
#########################################################################

"""
Utility for managing Illumina sequencing run directories

"""

#######################################################################
# Imports
#######################################################################

import optparse
import os
import sys
import logging
import auto_process_ngs.settings
from auto_process_ngs.utils import AnalysisDir

__version__ = auto_process_ngs.settings.version

#######################################################################
# Main program
#######################################################################

if __name__ == "__main__":
    # Process command line
    p = optparse.OptionParser(usage="%prog [OPTIONS] DIR",
                              version="%prog "+__version__,
                              description="Report information on processed Illumina "
                              "sequence data in ANALYSIS_DIR.")
    p.add_option('--debug',action='store_true',dest='debug',default=False,
                 help="Turn on debugging output")
    options,args = p.parse_args()
    if len(args) != 1:
        p.error("Need to supply a directory to examine")
    if options.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    # Acquire data for putative analysis directory
    dirn = os.path.abspath(args[0])
    print "Examining %s" % dirn
    analysis_dir = AnalysisDir(dirn)
    if not (analysis_dir.n_projects and analysis_dir.n_sequencing_data):
        logging.warning("Not an analysis directory")
        sys.exit()
    # Report projects from sequencing data
    if analysis_dir.n_sequencing_data:
        for data in analysis_dir.sequencing_data:
            print "Sequencing data projects found in '%s':" % \
                os.path.basename(data.unaligned_dir)
            for project in data.projects:
                print "%s" % project.name
