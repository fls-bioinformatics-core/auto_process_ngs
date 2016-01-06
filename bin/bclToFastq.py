#!/usr/bin/env python
#
#     bclToFastq.py: wrapper for converting Illumina bcl files to fastq
#     Copyright (C) University of Manchester 2013-2016 Peter Briggs
#
########################################################################
#
# bclToFastq.py
#
#########################################################################

"""bclToFastq.py

Wrapper for running bcl to fastq conversion using the Illumina
CASAVA or bcl2fastq software packages.

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
import re
import time
import logging
import optparse
import subprocess
import bcftbx.IlluminaData as IlluminaData
import auto_process_ngs.utils as utils
from auto_process_ngs.bcl2fastq_utils import run_bcl2fastq_1_8

#######################################################################
# Functions
#######################################################################

def run_bcl_to_fastq(basecalls_dir,sample_sheet,output_dir="Unaligned",
                     mismatches=None,
                     bases_mask=None,
                     nprocessors=None,
                     force=False,
                     ignore_missing_bcl=False,
                     ignore_missing_stats=False,
                     ignore_missing_control=False):
    """
    Run the bcl-to-fastq conversion software

    This wrapper function determines which version of the
    bcl-to-fastq conversion software is available and
    invokes the appropriate library function to execute it.

    """
    # Determine bcl2fastq version
    bcl2fastq_info = utils.bcl_to_fastq_info()
    bcl2fastq_version = bcl2fastq_info[2]
    # Return with error code if no version detected
    if bcl2fastq_version is None:
        logging.error("Cannot determine bcl2fastq software version")
        return 1
    # Run the appropriate pipeline based on the version
    if bcl2fastq_version.startswith('1.8.'):
        # 1.8.* pipeline
        return run_bcl2fastq_1_8(basecalls_dir,
                                 sample_sheet,
                                 output_dir=output_dir,
                                 mismatches=mismatches,
                                 bases_mask=bases_mask,
                                 nprocessors=nprocessors,
                                 force=force,
                                 ignore_missing_bcl=ignore_missing_bcl,
                                 ignore_missing_stats=ignore_missing_stats,
                                 ignore_missing_control=ignore_missing_control)
    else:
        # Unimplemented version
        logging.error("Don't know how to run bcl2fastq version %s" %
                      bcl2fastq_version)
        return 1

#######################################################################
# Main program
#######################################################################

if __name__ == '__main__':
    p = optparse.OptionParser(usage="%prog [OPTIONS] ILLUMINA_RUN_DIR OUTPUT_DIR [ SAMPLE_SHEET ]",
                              version="%prog "+__version__,
                              description="Wrapper to automate the Illumina bcl to fastq "
                              "conversion process: runs configureBclToFastq.pl followed "
                              "by make step. ILLUMINA_RUN_DIR is the top-level directory "
                              "of the Illumina run to be processed; output will be written "
                              "to OUTPUT_DIR. Optionally a SAMPLE_SHEET  file can also be "
                              "specified, otherwise the SampleSheet.csv file in the BaseCalls "
                              "directory will be used (if present).")
    p.add_option('--nmismatches',action="store",dest="nmismatches",default=None,
                 help="set number of mismatches to allow; recommended values are 0 "
                 "for samples without multiplexing, 1 for multiplexed samples with tags "
                 "of length 6 or longer (see the CASAVA user guide for details of "
                 "the --nmismatches option)")
    p.add_option('--use-bases-mask',action="store",dest="bases_mask",default=None,
                 help="specify a bases-mask string to tell CASAVA how to use each cycle. "
                 "The supplied value is passed directly to configureBcltoFastq.pl "
                 "(see the CASAVA user guide for details of how --use-bases-mask "
                 "works)")
    p.add_option('--nprocessors',action="store",dest="nprocessors",default=None,
                 help="set the number of processors to use (defaults to 1). "
                 "This is passed to the -j option of the 'make' step after running "
                 "configureBcltoFastq.pl (see the CASAVA user guide for details of "
                 "how -j works)")
    p.add_option('--ignore-missing-bcl',action="store_true",
                 dest="ignore_missing_bcl",default=False,
                 help="interpret missing bcl files as no call "
                 "(see the CASAVA user guide for details of how --ignore-missing-bcl "
                 "works)")
    p.add_option('--ignore-missing-stats',action="store_true",
                 dest="ignore_missing_stats",default=False,
                 help="fill in with zeroes when *.stats files are missing "
                 "(see the CASAVA user guide for details of how --ignore-missing-stats "
                 "works)")
    p.add_option('--ignore-missing-control',action="store_true",
                 dest="ignore_missing_control",default=False,
                 help="interpret missing control files as not-set control bits "
                 "(see the CASAVA user guide for details of how --ignore-missing-control "
                 "works)")
        
    options,args = p.parse_args()
    if not (2 <= len(args) <=3):
        p.error("input is an input directory, output directory and an optional sample sheet")
    # Locate run directory (and strip any trailing slash)
    illumina_run_dir = os.path.abspath(args[0].rstrip(os.sep))
    if not os.path.isdir(illumina_run_dir):
        logging.error("%s: doesn't exist or is not a directory" % illumina_run_dir)
        sys.exit(1)
    illumina_run = IlluminaData.IlluminaRun(illumina_run_dir)
    # Output directory
    output_dir = os.path.abspath(args[1].rstrip(os.sep))
    # Sample sheet
    if len(args) == 3:
        sample_sheet = os.path.abspath(args[2])
    else:
        sample_sheet = illumina_run.sample_sheet_csv
    # Bases mask
    if options.bases_mask is not None:
        bases_mask = options.bases_mask
    else:
        bases_mask = IlluminaData.IlluminaRunInfo(illumina_run.runinfo_xml).bases_mask
    # Report settings
    print "Illumina run directory: %s" % illumina_run.run_dir
    print "Platform              : %s" % illumina_run.platform
    print "Basecalls directory   : %s" % illumina_run.basecalls_dir
    print "Bcl file extension    : %s" % illumina_run.bcl_extension
    print "SampleSheet.csv file  : %s" % sample_sheet
    print "Output dir            : %s" % output_dir
    print "Nmismatches           : %s" % options.nmismatches
    print "Bases mask            : %s" % bases_mask
    print "Nprocessors           : %s" % options.nprocessors
    print "Ignore missing bcl    : %s" % options.ignore_missing_bcl
    print "Ignore missing stats  : %s" % options.ignore_missing_stats
    print "Ignore missing control: %s" % options.ignore_missing_control
    # Run bclToFastq conversion
    status = run_bcl_to_fastq(illumina_run.basecalls_dir,
                              sample_sheet,
                              output_dir=output_dir,
                              mismatches=options.nmismatches,
                              bases_mask=bases_mask,
                              force=True,
                              nprocessors=options.nprocessors,
                              ignore_missing_bcl=options.ignore_missing_bcl,
                              ignore_missing_stats=options.ignore_missing_stats,
                              ignore_missing_control=options.ignore_missing_control)
    print "bclToFastq returncode: %s" % status
    if status != 0:
        logging.error("bclToFastq failure")
    sys.exit(status)
