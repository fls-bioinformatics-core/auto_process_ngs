#!/usr/bin/env python
#
#     qc.constants.py: constants for QC pipeline
#     Copyright (C) University of Manchester 2019 Peter Briggs
#

"""
Provides the following constants for the QC pipeline:

- FASTQ_SCREENS: tuple of screen names
- PROTOCOLS: tuple of QC protocol names
"""
######################################################################
# QC constants
######################################################################

# QC protocols
PROTOCOLS = ('standardPE',
             'standardSE',
             'singlecell',
             '10x_scATAC')

# Screen names
FASTQ_SCREENS = ('model_organisms',
                 'other_organisms',
                 'rRNA',)
