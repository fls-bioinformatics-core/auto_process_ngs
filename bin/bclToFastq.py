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

__version__ = "0.2.0"

# Set of versions that this is known to work with
BCL2FASTQ_VERSIONS = ('1.8','2.17',)

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
from auto_process_ngs.bcl2fastq_utils import available_bcl2fastq_versions
from auto_process_ngs.bcl2fastq_utils import bcl_to_fastq_info
from auto_process_ngs.bcl2fastq_utils import run_bcl2fastq_1_8
from auto_process_ngs.bcl2fastq_utils import run_bcl2fastq_2_17

#######################################################################
# Functions
#######################################################################

def main():
    p = optparse.OptionParser(
        usage="%prog [OPTIONS] ILLUMINA_RUN_DIR OUTPUT_DIR [ SAMPLE_SHEET ]",
        version="%prog "+__version__,
        description="Wrapper to automate the Illumina bcl to fastq "
        "conversion process. It will either run the CASAVA/bcl2fastq v1.8 "
        "configureBclToFastq.pl/make pipeline or bcl2fastq v2 directly, "
        "depending on which software package is detected. ILLUMINA_RUN_DIR "
        "is the top-level directory of the Illumina run to be processed; "
        "output will be written to OUTPUT_DIR. Optionally a SAMPLE_SHEET "
        "file can also be specified, otherwise the SampleSheet.csv file in "
        "the BaseCalls directory will be used (if present).")
    # Options common to both bcl2fastq/bcl2fastq v2
    p.add_option('--nmismatches',action="store",dest="nmismatches",
                 default=None,
                 help="set number of mismatches to allow; recommended "
                 "values are 0 for samples without multiplexing, 1 for "
                 "multiplexed samples with tags of length 6 or longer "
                 "(CASAVA/bcl2fastq v1.8 --mismatches option, bcl2fastq "
                 "v2 --barcode-mismatches option)")
    p.add_option('--use-bases-mask',action="store",dest="bases_mask",
                 default=None,
                 help="specify a bases-mask string to tell CASAVA how "
                 "to use each cycle (the supplied value is passed "
                "to the --use-bases-mask option)")
    p.add_option('--nprocessors',action="store",dest="nprocessors",
                 default=None,
                 help="set the number of processors to use (defaults to "
                 "1; for CASAVA/bcl2fastq v1.8 this is passed to the "
                 "-j option of the 'make' step after running "
                 "configureBcltoFastq.pl, for bcl2fastq v2 this is "
                 "the maximum number of CPUs that should be used by "
                 "the -r, -d, -p and -w options)")
    p.add_option('--ignore-missing-bcl',action="store_true",
                 dest="ignore_missing_bcl",default=False,
                 help="interpret missing bcl files as no call "
                 "(CASAVA/bcl2fastq v1.8 --ignore-missing-bcl option, "
                 "bcl2fastq v2 --ignore-missing-bcls option)")
    p.add_option('--bcl2fastq_path',action="store",
                 dest="bcl2fastq_path",default=None,
                 help="explicitly specify the path to the CASAVA or "
                 "bcl2fastq software to use.")
    # CASAVA/bcl2fastq 1.8.* only
    casava = optparse.OptionGroup(p,'CASAVA/bcl2fastq v1.8 only')
    casava.add_option('--ignore-missing-stats',action="store_true",
                      dest="ignore_missing_stats",default=False,
                      help="fill in with zeroes when *.stats files are missing "
                      "(see the CASAVA user guide for details of how "
                      "--ignore-missing-stats works)")
    casava.add_option('--ignore-missing-control',action="store_true",
                      dest="ignore_missing_control",default=False,
                 help="interpret missing control files as not-set control "
                      "bits (see the CASAVA user guide for details of how "
                      "--ignore-missing-control works)")
    p.add_option_group(casava)
    # bcl2fastq 2 only
    bcl2fastq2 = optparse.OptionGroup(p,'bcl2fastq v2 only')
    bcl2fastq2.add_option('--no-lane-splitting',action="store_true",
                          dest="no_lane_splitting",default=False,
                          help="Don't split output FASTQ files by lane")
    p.add_option_group(bcl2fastq2)

    options,args = p.parse_args()
    if not (2 <= len(args) <=3):
        p.error("input is an input directory, output directory and an "
                "optional sample sheet")
    # Acquire bcl2fastq software
    bcl2fastq = available_bcl2fastq_versions(paths=(options.bcl2fastq_path,))
    if not bcl2fastq:
        logging.error("No bcl2fastq software found")
        return 1
    else:
        bcl2fastq_exe = bcl2fastq[0]
    # Determine bcl2fastq version
    bcl2fastq_info = bcl_to_fastq_info(bcl2fastq_exe)
    if bcl2fastq_info[0] is None:
        logging.error("No bcl2fastq software found")
        return 1
    print "Using conversion software from %s" % os.path.dirname(
        bcl2fastq_info[0])
    # Return with error code if no version detected
    bcl2fastq_package = bcl2fastq_info[1]
    bcl2fastq_version = bcl2fastq_info[2]
    if bcl2fastq_version is None:
        logging.error("Cannot determine bcl2fastq software version")
        return 1
    print "Package: %s" % bcl2fastq_package
    print "Version: %s" % bcl2fastq_version
    known_version = None
    for version in BCL2FASTQ_VERSIONS:
        if bcl2fastq_version.startswith("%s." % version):
            known_version = version
            break
    if known_version is None:
        # Unimplemented version
        logging.error("Don't know how to run bcl2fastq version %s" %
                      bcl2fastq_version)
        return 1
    # Locate run directory (and strip any trailing slash)
    illumina_run_dir = os.path.abspath(args[0].rstrip(os.sep))
    if not os.path.isdir(illumina_run_dir):
        logging.error("%s: doesn't exist or is not a directory" %
                      illumina_run_dir)
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
        bases_mask = IlluminaData.IlluminaRunInfo(
            illumina_run.runinfo_xml).bases_mask
    # Report settings
    print "Illumina run directory: %s" % illumina_run.run_dir
    print "Basecalls directory   : %s" % illumina_run.basecalls_dir
    print "Platform              : %s" % illumina_run.platform
    print "Bcl file extension    : %s" % illumina_run.bcl_extension
    print "SampleSheet.csv file  : %s" % sample_sheet
    print "Output dir            : %s" % output_dir
    print "Nmismatches           : %s" % options.nmismatches
    print "Bases mask            : %s" % bases_mask
    print "Nprocessors           : %s" % options.nprocessors
    print "Ignore missing bcl    : %s" % options.ignore_missing_bcl
    if known_version == '1.8':
        print "Ignore missing stats  : %s" % options.ignore_missing_stats
        print "Ignore missing control: %s" % options.ignore_missing_control
    elif known_version == '2.17':
        print "No lane splitting     : %s" % options.no_lane_splitting
    # Run bclToFastq conversion based on the version
    if bcl2fastq_version.startswith('1.8.'):
        # 1.8.* pipeline
        status = run_bcl2fastq_1_8(
            illumina_run.basecalls_dir,
            sample_sheet,
            output_dir=output_dir,
            mismatches=options.nmismatches,
            bases_mask=options.bases_mask,
            nprocessors=options.nprocessors,
            force=True,
            ignore_missing_bcl=options.ignore_missing_bcl,
            ignore_missing_stats=options.ignore_missing_stats,
            ignore_missing_control=options.ignore_missing_control
        )
    elif bcl2fastq_version.startswith('2.17.'):
        # bcl2fastq 2.17.*
        if options.nprocessors is not None:
            # Explicitly set number of threads for each stage
            nprocessors=int(options.nprocessors)
            loading_threads=min(4,nprocessors)
            writing_threads=min(4,nprocessors)
            demultiplexing_threads=max(int(float(nprocessors)*0.2),
                                       nprocessors)
            processing_threads=nprocessors
            print "Explicitly setting number of threads for each stage:"
            print "Loading (-r)       : %d" % loading_threads
            print "Demultiplexing (-d): %d" % demultiplexing_threads
            print "Processing (-p)    : %d" % processing_threads
            print "Writing (-w)       : %d" % writing_threads
        else:
            # Use the defaults
            loading_threads = None
            demultiplexing_threads = None
            processing_threads = None
            writing_threads = None
        # Run the bcl to fastq conversion
        status = run_bcl2fastq_2_17(
            illumina_run.run_dir,
            sample_sheet,
            output_dir=output_dir,
            mismatches=options.nmismatches,
            bases_mask=options.bases_mask,
            ignore_missing_bcl=options.ignore_missing_bcl,
            no_lane_splitting=options.no_lane_splitting,
            loading_threads=loading_threads,
            demultiplexing_threads=demultiplexing_threads,
            processing_threads=processing_threads,
            writing_threads=writing_threads
        )
    print "bclToFastq returncode: %s" % status
    if status != 0:
        logging.error("bclToFastq failure")
    return status

#######################################################################
# Main program
#######################################################################

if __name__ == '__main__':
    sys.exit(main())
