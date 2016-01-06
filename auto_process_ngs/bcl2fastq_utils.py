#!/usr/bin/env python
#
#     bcl2fastq_utils.py: utility functions for bcl2fastq conversion
#     Copyright (C) University of Manchester 2013-2014 Peter Briggs
#
########################################################################
#
# bclToFastq.py
#
#########################################################################

"""bcl2fastq_utils.py

Utility functions for bcl to fastq conversion operations:

make_custom_sample_sheet: create a fixed copy of a sample sheet file
get_nmismatches: determine number of mismatches from bases mask
get_bases_mask: get a bases mask string

"""

#######################################################################
# Imports
#######################################################################
import bcftbx.IlluminaData as IlluminaData

#######################################################################
# Functions
#######################################################################

def make_custom_sample_sheet(input_sample_sheet,output_sample_sheet=None,
                             fmt=None):
    """Creates a corrected copy of a sample sheet file

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
    """Get bases mask string

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
    """Determine number of mismatches from bases mask

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
