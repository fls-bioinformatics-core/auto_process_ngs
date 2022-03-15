#!/usr/bin/env python
#
#     icell8.constants.py: constants for ICELL8 data handling
#     Copyright (C) University of Manchester 2018 Peter Briggs
#

"""
Constants for ICELL8 data handling.
"""
######################################################################
# ICELL8 magic numbers
######################################################################

INLINE_BARCODE_LENGTH = 11
UMI_LENGTH = 14

######################################################################
# Other constants
######################################################################

SAMPLENAME_ILLEGAL_CHARS = "?()[]/\\=+<>:;\"',*^|& \t"

MAXIMUM_BATCH_SIZE = 100000000
DEFAULT_BATCH_SIZE = 5000000
