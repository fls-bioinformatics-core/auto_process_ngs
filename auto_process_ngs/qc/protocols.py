#!/usr/bin/env python
#
#     protocols: define and handle QC protocols
#     Copyright (C) University of Manchester 2022 Peter Briggs
#

"""
QC protocol definitions and utility functions for handling protocols
and modules for QC of analysis projects.

QC protocols are defined with the ``QC_PROTOCOLS`` dictionary, where
each key is a protocol name and the corresponding values are
dictionaries which specify the reads used for sequence data ("data reads")
and for QC ("QC reads"), along with a list of QC modules that comprise
the protocol.

For example:

::

    "ExampleProtocol": {
        "reads": { "data": ('r1','r2'), "qc": ('r2') },
        "qc_modules": ['fastqc','fastq_screen','sequence_lengths']
    }

The available QC modules are those supported by the ``QCVerifier``
class in the ``check_outputs`` method; new modules must be added
there before they can be specified in protocol definitions.

Optional modifiers can also be added to QC module specifications,
using the format ``NAME(KEY=VALUE;...)``.

For example:

::

    cellranger_count(cellranger_version=*;cellranger_refdata=*)

The available modifiers are the same as the parameter list for the
``check_outputs`` in the ``QCVerifier`` class.

This module also provides the following functions:

- fetch_protocol_definition: get the definition for a QC protocol
"""

#######################################################################
# Imports
#######################################################################

from bcftbx.utils import AttributeDictionary

#######################################################################
# Data
#######################################################################

# QC protocol definitions

QC_PROTOCOLS = {

    "standardSE": {
        "reads": {
            "data": ('r1',),
            "qc": ('r1',)
        },
        "qc_modules": ['fastqc',
                       'fastq_screen',
                       'sequence_lengths',
                       'strandedness',
                       'multiqc']
    },

    "standardPE": {
        "reads": {
            "data": ('r1','r2',),
            "qc": ('r1','r2',)
        },
        "qc_modules": ['fastqc',
                       'fastq_screen',
                       'sequence_lengths',
                       'strandedness',
                       'multiqc']
    },

    "singlecell": {
        "reads": {
        "data": ('r2',),
            "qc": ('r1','r2',)
        },
        "qc_modules": ['fastqc',
                       'fastq_screen',
                       'sequence_lengths',
                       'strandedness',
                       'multiqc']
    },

    "10x_scRNAseq": {
        "reads": {
            "data": ('r2',),
            "qc": ('r1','r2',)
        },
        "qc_modules": ['fastqc',
                       'fastq_screen',
                       'sequence_lengths',
                       'strandedness',
                       'multiqc',
                       'cellranger_count']
    },

    "10x_snRNAseq": {
        "reads": {
            "data": ('r2',),
            "qc": ('r1','r2',)
        },
        "qc_modules": ['fastqc',
                       'fastq_screen',
                       'sequence_lengths',
                       'strandedness',
                       'multiqc',
                       'cellranger_count']
    },

    "10x_scATAC": {
        "reads": {
            "data": ('r1','r3',),
            "qc": ('r1','r3',)
        },
        "qc_modules": ['fastqc',
                       'fastq_screen',
                       'sequence_lengths',
                       'strandedness',
                       'multiqc',
                       'cellranger-atac_count']
    },

    "10x_Multiome_GEX": {
        "reads": {
            "data": ('r2',),
            "qc": ('r1','r2',)
        },
        "qc_modules": ['fastqc',
                       'fastq_screen',
                       'sequence_lengths',
                       'strandedness',
                       'multiqc',
                       'cellranger_count(cellranger_version=*;'
                       'cellranger_refdata=*)',
                       'cellranger-arc_count']
    },

    "10x_Multiome_ATAC": {
        "reads": {
            "data": ('r1','r3',),
            "qc": ('r1','r3',)
        },
        "qc_modules": ['fastqc',
                       'fastq_screen',
                       'sequence_lengths',
                       'strandedness',
                       'multiqc',
                       'cellranger-arc_count',
                       'cellranger-atac_count(cellranger_version=*;'
                       'cellranger_refdata=*)']
    },
    "10x_CellPlex": {
        "reads": {
            "data": ('r2',),
            "qc": ('r1','r2',)
        },
        "qc_modules": ['fastqc',
                       'fastq_screen',
                       'sequence_lengths',
                       'strandedness',
                       'multiqc',
                       'cellranger_count',
                       'cellranger_multi']
    },

    "10x_Visium": {
        "reads": {
            "data": ('r1','r2',),
            "qc": ('r1','r2',)
        },
        "qc_modules": ['fastqc',
                       'fastq_screen',
                       'sequence_lengths',
                       'strandedness',
                       'multiqc']
    }
}

#######################################################################
# Functions
#######################################################################

def fetch_protocol_definition(name):
    """
    Return the definition for a QC protocol

    Arguments:
      name (str): name of the QC protocol

    Returns:
      Tuple: definition as a tuple of the form
        (reads,qc_modules) where 'reads' is an
        AttributeDictionary with elements 'data'
        and 'qc' (listing data and QC reads
        respectively) and 'qc_modules' is a list
        of QC module definitions.
    """
    if name not in QC_PROTOCOLS:
        raise KeyError("%s: undefined QC protocol" % name)
    protocol_defn = QC_PROTOCOLS[name]
    reads = AttributeDictionary()
    try:
        reads['data'] = list(protocol_defn['reads']['data'])
        reads['qc'] = list(protocol_defn['reads']['qc'])
        qc_modules = [m for m in protocol_defn['qc_modules']]
    except KeyError as ex:
        raise Exception("%s: exception loading QC protocol "
                        "definition: %s" % (name,ex))
    return (reads,qc_modules)
