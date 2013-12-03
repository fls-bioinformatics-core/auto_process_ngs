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

__version__ = "0.0.1"

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

#######################################################################
# Functions
#######################################################################

def fastq_statistics(illumina_data):
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

    Returns:
      Populated TabFile object containing the statistics.

    """
    # Set up TabFile to hold the data collected
    stats = TabFile.TabFile(column_names=('Project',
                                          'Sample',
                                          'Fastq',
                                          'Size',
                                          'Nreads',
                                          'Paired_end'))
    # Get number of reads and file size for each file across all
    # projects and samples
    for project in illumina_data.projects:
        for sample in project.samples:
            for fastq in sample.fastq:
                print "* %s" % fastq
                fq = os.path.join(sample.dirn,fastq)
                nreads = FASTQFile.nreads(fq)
                fsize = os.path.getsize(fq)
                stats.append(data=(project.name,
                                   sample,fastq,
                                   bcf_utils.format_file_size(fsize),
                                   nreads,
                                   'Y' if sample.paired_end else 'N'))
    # Gather same information for undetermined reads (if present)
    if illumina_data.undetermined is not None:
        for sample in illumina_data.undetermined.samples:
            for fastq in sample.fastq:
                fq = os.path.join(sample.dirn,fastq)
                nreads = FASTQFile.nreads(fq)
                fsize = os.path.getsize(fq)
                stats.append(data=(illumina_data.undetermined.name,
                                   sample,fastq,
                                   bcf_utils.format_file_size(fsize),
                                   nreads,
                                   'Y' if sample.paired_end else 'N'))
    # Return the data
    return stats

#######################################################################
# Main program
#######################################################################

if __name__ == '__main__':

    # Process command line
    p = optparse.OptionParser(usage="%prog [OPTIONS] ILLUMINA_RUN_DIR",
                              version="%prog "+__version__,
                              description="Generate statistics for fastq files in "
                              "ILLUMINA_RUN_DIR (top-level directory of the Illumina "
                              "run to be processed.")
    p.add_option('-o','--output',action="store",dest="stats_file",default='statistics.info',
                 help="name of output file (default is 'statistics.info')")
    p.add_option("--unaligned",action="store",dest="unaligned_dir",default="Unaligned",
                 help="specify an alternative name for the 'Unaligned' directory "
                 "containing the fastq.gz files")
    options,args = p.parse_args()
    if len(args) != 1:
        p.error("expects a single argument (input directory)")

    # Get the data
    try:
        illumina_data = IlluminaData.IlluminaData(args[0],
                                                  unaligned_dir=options.unaligned_dir)
    except IlluminaData.IlluminaDataError,ex:
        logging.error("Failed to get data from %s: %s" % (args[0],ex))
        sys.exit(1)

    # Generate statistics for fastq files
    stats = fastq_statistics(illumina_data)
    stats.write(options.stats_file,include_header=True)
    print "Statistics written to %s" % options.stats_file
