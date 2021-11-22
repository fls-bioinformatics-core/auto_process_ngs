#!/usr/bin/env python
#
#     qc.constants.py: constants for QC pipeline
#     Copyright (C) University of Manchester 2019,2021 Peter Briggs
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
             '10x_scRNAseq',
             '10x_snRNAseq',
             '10x_scATAC',
             '10x_Visium',
             '10x_Multiome_GEX',
             '10x_Multiome_ATAC',
             '10x_CellPlex',
             'ICELL8_scATAC',)

# Screen names
FASTQ_SCREENS = ('model_organisms',
                 'other_organisms',
                 'rRNA',)

# Sequence deduplication thresholds
SEQUENCE_DEDUP_CUTOFFS = {
    "warn": 0.3,
    "fail": 0.2,
}
