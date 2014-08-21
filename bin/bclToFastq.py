#!/bin/env python
#
#     bclToFastq.py: wrapper for converting Illumina bcl files to fastq
#     Copyright (C) University of Manchester 2013 Peter Briggs
#
########################################################################
#
# bclToFastq.py
#
#########################################################################

"""bclToFastq.py

Wrapper for running bcl to fastq conversion using the Illumina
configureBclToFastq.pl pipeline.

"""

#######################################################################
# Module metadata
#######################################################################

__version__ = "0.0.6"

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
import applications

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
    """Wrapper for running the CASAVA bcl to fastq pipeline

    Runs the CASAVA bcl to fastq pipeline, specifically:
    
    1. Executes the 'configureBclToFastq.pl' script to generate a
       Makefile to perform the conversion, then
    2. Executes 'make' using this Makefile to generate fastq files
       from bcl files.

    Arguments:
      basecalls_dir: path to the top-level directory holding the bcl
        files (typically 'Data/Intensities/Basecalls/' subdirectory)
      sample_sheet: path to the sample sheet file to use
      output_dir: optional, path to the output directory. Defaults to
        'Unaligned'. If this directory already exists then the
        conversion will fail unless the force option is set to True
      mismatches: optional, specify maximum number of mismatched bases
        allowed for matching index sequences during multiplexing.
        Recommended values are zero for indexes shorter than 6 base
        pairs, 1 for indexes of 6 or longer
        (If not specified and bases_mask is supplied then mismatches
        will be derived automatically from the bases mask string)
      bases_mask: optional, specify string indicating how to treat
        each cycle within each read e.g. 'y101,I6,y101'
      nprocessors: optional, number of processors to use when running
        'make' step
      force: optional, if True then force overwrite of an existing
        output directory (default is False).
      ignore_missing_bcl: optional, if True then interpret missing bcl
        files as no call (default is False)
      ignore_missing_stats: optional, if True then fill in with zeroes
        when *.stats files are missing (default is False)
      ignore_missing_control: optional, if True then interpret missing
        control files as not-set control bits (default is False)

    Returns:
      0 on success; if a problem is encountered then returns -1 for
      errors within the function (e.g. missing Makefile) or the exit
      code from the failed program.

    """
    # Set up and run configureBclToFastq
    configure_cmd = applications.bcl2fastq.configureBclToFastq(
        basecalls_dir,
        sample_sheet,
        output_dir=output_dir,
        mismatches=mismatches,
        bases_mask=bases_mask,
        force=force,
        ignore_missing_bcl=ignore_missing_bcl,
        ignore_missing_stats=ignore_missing_stats,
        ignore_missing_control=ignore_missing_control
    )
    if not configure_cmd.has_exe:
        logging.error("'%s' missing, cannot run" % configure_cmd.command)
        return -1
    print "Running command: %s" % configure_cmd
    returncode = configure_cmd.run_subprocess()
    # Check returncode
    if returncode != 0:
        logging.error("configureToBclFastq.pl returned %s" % returncode)
        return returncode
    # Check outputs (directory and makefile)
    if not os.path.isdir(output_dir):
        logging.error("Output directory '%s' not found" % output_dir)
        return -1
    makefile = os.path.join(output_dir,'Makefile')
    if not os.path.isfile(makefile):
        logging.error("Makefile not found in %s" % output_dir)
        return -1
    # Set up and run make command
    make_cmd = applications.general.make(makefile=makefile,
                                         working_dir=output_dir,
                                         nprocessors=nprocessors)
    if not make_cmd.has_exe:
        logging.error("'%s' missing, cannot run" % make_cmd.command)
        return -1
    print "Running command: %s" % make_cmd
    returncode = make_cmd.run_subprocess()
    # Check returncode
    if returncode != 0:
        logging.error("make returned %s" % returncode)
    return returncode

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
