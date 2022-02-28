#!/usr/bin/env python
#
#     qc.constants.py: constants for QC pipeline
#     Copyright (C) University of Manchester 2019-2022 Peter Briggs
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
from .protocols import QC_PROTOCOLS
PROTOCOLS = tuple(QC_PROTOCOLS.keys())

# Screen names
FASTQ_SCREENS = ('model_organisms',
                 'other_organisms',
                 'rRNA',)

# Sequence deduplication thresholds
SEQUENCE_DEDUP_CUTOFFS = {
    "warn": 0.3,
    "fail": 0.2,
}
