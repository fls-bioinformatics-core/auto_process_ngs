#!/usr/bin/env python
#
#     bcl2fastq_utils.py: utility functions for bcl2fastq conversion
#     Copyright (C) University of Manchester 2013-2016 Peter Briggs
#
########################################################################
#
# bclToFastq.py
#
#########################################################################

"""
bcl2fastq_utils.py

Utility functions for bcl to fastq conversion operations:

make_custom_sample_sheet: create a fixed copy of a sample sheet file
get_nmismatches: determine number of mismatches from bases mask
get_bases_mask: get a bases mask string
run_bcl2fastq_1_8: run bcl-to-fastq conversion from CASAVA/bcl2fastq 1.8.*
run_bcl2fastq_2_17: run bcl-to-fastq conversion from bcl2fastq 2.17.*

"""

#######################################################################
# Imports
#######################################################################
import os
import logging
import auto_process_ngs.applications as applications
import bcftbx.IlluminaData as IlluminaData

#######################################################################
# Functions
#######################################################################

def make_custom_sample_sheet(input_sample_sheet,output_sample_sheet=None,
                             fmt=None):
    """
    Creates a corrected copy of a sample sheet file

    Creates and returns a SampleSheet object with a copy of the
    input sample sheet, with any illegal or duplicated names fixed.
    Optionally also writes the updated sample sheet data to a new
    file.

    Arguments:
      input_sample_sheet (str): name and path of the original sample
        sheet file
      output_sample_sheet (str): (optional) name and path to write
        updated sample sheet to, or `None`
      fmt (str): (optional) format for the output sample sheet,
        either 'CASAVA' or 'IEM'; if this `None` then the format of
        the original file will be used

    Returns:
      SampleSheet object with the data for the corrected sample
      sheet.

    """
    # Load the sample sheet data
    sample_sheet = IlluminaData.SampleSheet(input_sample_sheet)
    # Determine the column names for this format
    if sample_sheet.format == 'CASAVA':
        sample_col = 'SampleID'
        project_col = 'SampleProject'
    elif sample_sheet.format == 'IEM':
        sample_col = 'Sample_ID'
        project_col = 'Sample_Project'
    else:
        raise Exception("Unknown sample sheet format: %s" %
                        sample_sheet.format)
    # Add project names if not supplied
    for line in sample_sheet:
        if not line[project_col]:
            line[project_col] = line[sample_col]
    # Fix other problems
    sample_sheet.fix_illegal_names()
    sample_sheet.fix_duplicated_names()
    # Write out new sample sheet
    if output_sample_sheet is not None:
        sample_sheet.write(output_sample_sheet,fmt=fmt)
    return sample_sheet

def get_bases_mask(run_info_xml,sample_sheet_file):
    """
    Get bases mask string

    Generates initial bases mask based on data in RunInfo.xml (which
    says how many reads there are, how many cycles in each read, and
    which are index reads). Then updates this using the barcode
    information in the sample sheet file.

    Arguments:
      run_info_xml: name and path of RunInfo.xml file from the
        sequencing run
      sample_sheet_file: name and path of sample sheet file.

    Returns:
      Bases mask string e.g. 'y101,I6'. 

    """
    # Get initial bases mask
    bases_mask = IlluminaData.IlluminaRunInfo(run_info_xml).bases_mask
    print "Bases mask: %s (from RunInfo.xml)" % bases_mask
    # Update bases mask from sample sheet
    example_barcode = IlluminaData.get_casava_sample_sheet(sample_sheet_file)[0]['Index']
    bases_mask = IlluminaData.fix_bases_mask(bases_mask,example_barcode)
    print "Bases mask: %s (updated for barcode sequence '%s')" % (bases_mask,
                                                                  example_barcode)
    return bases_mask

def get_nmismatches(bases_mask):
    """
    Determine number of mismatches from bases mask

    Automatically determines the maximum number of mismatches that shoud
    be allowed for a bcl to fastq conversion run, based on the tag
    length i.e. the length of the index barcode sequences.

    Tag lengths of 6 or more use 1 mismatch, otherwise use zero
    mismatches.

    The number of mismatches should be supplied to the bclToFastq
    conversion process.

    Arguments:
      bases_mask: bases mask string of the form e.g. 'y101,I6,y101'

    Returns:
      Integer value of number of mismatches. (If the bases mask doesn't
      contain any index reads then returns zero.)

    """
    for read in bases_mask.split(','):
        if read.startswith('I'):
            try:
                i = read.index('n')
                read = read[:i]
            except ValueError:
                pass
            index_length = int(read[1:].rstrip('n'))
            if index_length >= 6:
                return 1
            else:
                return 0
    # Failed to find any indexed reads
    return 0

def run_bcl2fastq_1_8(basecalls_dir,sample_sheet,
                      output_dir="Unaligned",
                      mismatches=None,
                      bases_mask=None,
                      nprocessors=None,
                      force=False,
                      ignore_missing_bcl=False,
                      ignore_missing_stats=False,
                      ignore_missing_control=False):
    """
    Wrapper for running the CASAVA/bcl2fastq 1.8.* pipeline

    Runs the CASAVA-style bcl to fastq pipeline, specifically:

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
    # Check the executable exists
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

def run_bcl2fastq_2_17(basecalls_dir,sample_sheet,
                       output_dir="Unaligned",
                       mismatches=None,
                       bases_mask=None,
                       nprocessors=None,
                       force=False,
                       ignore_missing_bcl=False,
                       no_lane_splitting=False):
    """
    Wrapper for running bcl2fastq 2.17.*

    Runs the bcl2fastq 2.17.* software to generate fastq files
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
      no_lane_splitting: optional, if True then output FASTQ files
        will not be split by lane (default is False, i.e. do split
        FASTQs by lane)

    Returns:
      0 on success; if a problem is encountered then returns -1 for
      errors within the function (e.g. missing Makefile) or the exit
      code from the failed program.

    """
    # Set up and run bcl2fastq2
    bcl2fastq2_cmd = applications.bcl2fastq.bcl2fastq2(
        basecalls_dir,
        sample_sheet,
        output_dir=output_dir,
        mismatches=mismatches,
        bases_mask=bases_mask,
        force=force,
        ignore_missing_bcl=ignore_missing_bcl,
        no_lane_splitting=no_lane_splitting
    )
    # Check the executable exists
    if not bcl2fastq2_cmd.has_exe:
        logging.error("'%s' missing, cannot run" % bcl2fastq2_cmd.command)
        return -1
    print "Running command: %s" % bcl2fastq2_cmd
    returncode = bcl2fastq2_cmd.run_subprocess()
    # Check returncode
    if returncode != 0:
        logging.error("bcl2fastq returned %s" % returncode)
        return returncode
    # Check outputs (directory and makefile)
    if not os.path.isdir(output_dir):
        logging.error("Output directory '%s' not found" % output_dir)
        return -1
    return returncode
