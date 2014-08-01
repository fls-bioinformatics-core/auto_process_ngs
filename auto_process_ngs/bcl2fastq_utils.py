#!/bin/env python
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

Utility functions for bcl to fastq conversion.

get_nmismatches: determine number of mismatches from bases mask

"""

#######################################################################
# Functions
#######################################################################

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
