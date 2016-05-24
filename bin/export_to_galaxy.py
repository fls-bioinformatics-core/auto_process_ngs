#!/usr/bin/env python
#
#     export_to_galaxy.py: export fastqs to Galaxy data libraries
#     Copyright (C) University of Manchester 2016 Peter Briggs
#

"""
Automatically export Fastqs to Galaxy data libraries

"""

#######################################################################
# Imports
#######################################################################

import logging
logging.basicConfig(format='%(levelname) 8s: %(message)s')

import optparse
from bcftbx.IlluminaData import split_run_name
from auto_process_ngs.utils import AnalysisDir
from auto_process_ngs.galaxy import build_library_directory
from auto_process_ngs.galaxy import create_data_library

#######################################################################
# Main program
#######################################################################

if __name__ == "__main__":

    # Parse command line
    p = optparse.OptionParser(usage="%prog GALAXY LIBRARY DEST RUNDIR")
    p.add_option('-n','--no-verify',action='store_true',dest='no_verify',
                 default=False,help="don't verify HTTPS connections")
    opts,args = p.parse_args()

    # Collect arguments
    if len(args) != 4:
        p.error("Wrong args")
    galaxy_url = args[0]
    library_name = args[1]
    dest = args[2]
    rundir = args[3]

    # Turn off SSL certificate verification?
    if opts.no_verify:
        logging.warning("SSL certificate verification disabled")
        turn_off_urllib3_warnings()

    # Get analysis dir
    analysis_dir = AnalysisDir(rundir)
    print "Run name: %s" % analysis_dir.run_name
    print "Run number: %s" % analysis_dir.metadata.run_number
    print "Date stamp: %s" % analysis_dir.date_stamp
    try:
        # Create the directory structure on the server
        build_library_directory(analysis_dir,dest)
        # Upload to the Galaxy instance
        create_data_library(galaxy_url,library_name,analysis_dir,dest,
                            no_verify=opts.no_verify)
    except ex:
        logging.critical("Failed to create data library: %s" % ex)
        sys.exit(1)
