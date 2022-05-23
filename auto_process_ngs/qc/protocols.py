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
dictionaries which specify the reads used for sequence data and for
index information, along with a list of QC modules that comprise
the protocol.

For example:

::

    "ExampleProtocol": {
        "reads": { "seq_data": ('r1','r3'), "index": ('r2') },
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

- determine_qc_protocol: get QC protocol for a project
- fetch_protocol_definition: get the definition for a QC protocol
- get_read_numbers: get the read numbers associated with a QC protocol
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
            "seq_data": ('r1',),
            "index": ()
        },
        "qc_modules": [
            'fastqc',
            'fastq_screen',
            'sequence_lengths',
            'strandedness'
        ]
    },

    "standardPE": {
        "reads": {
            "seq_data": ('r1','r2',),
            "index": ()
        },
        "qc_modules": [
            'fastqc',
            'fastq_screen',
            'sequence_lengths',
            'strandedness'
        ]
    },

    "singlecell": {
        "reads": {
            "seq_data": ('r2',),
            "index": ('r1',)
        },
        "qc_modules": [
            'fastqc',
            'fastq_screen',
            'sequence_lengths',
            'strandedness'
        ]
    },

    "10x_scRNAseq": {
        "reads": {
            "seq_data": ('r2',),
            "index": ('r1',)
        },
        "qc_modules": [
            'fastqc',
            'fastq_screen',
            'sequence_lengths',
            'strandedness',
            'cellranger_count'
        ]
    },

    "10x_snRNAseq": {
        "reads": {
            "seq_data": ('r2',),
            "index": ('r1',)
        },
        "qc_modules": [
            'fastqc',
            'fastq_screen',
            'sequence_lengths',
            'strandedness',
            'cellranger_count'
        ]
    },

    "10x_scATAC": {
        "reads": {
            "seq_data": ('r1','r3',),
            "index": ()
        },
        "qc_modules": [
            'fastqc',
            'fastq_screen',
            'sequence_lengths',
            'strandedness',
            'cellranger-atac_count'
        ]
    },

    "10x_Multiome_GEX": {
        "reads": {
            "seq_data": ('r2',),
            "index": ('r1',)
        },
        "qc_modules": [
            'fastqc',
            'fastq_screen',
            'sequence_lengths',
            'strandedness',
            'cellranger_count(cellranger_version=*;'
                             'cellranger_refdata=*)',
            'cellranger-arc_count'
        ]
    },

    "10x_Multiome_ATAC": {
        "reads": {
            "seq_data": ('r1','r3',),
            "index": ()
        },
        "qc_modules": [
            'fastqc',
            'fastq_screen',
            'sequence_lengths',
            'strandedness',
            'cellranger-arc_count',
            'cellranger-atac_count(cellranger_version=*;'
                                  'cellranger_refdata=*)'
        ]
    },

    "10x_CellPlex": {
        "reads": {
            "seq_data": ('r2',),
            "index": ('r1',)
        },
        "qc_modules": [
            'fastqc',
            'fastq_screen',
            'sequence_lengths',
            'strandedness',
            'cellranger_count(cellranger_use_multi_config=True)',
            'cellranger_multi'
        ]
    },

    "10x_Visium": {
        "reads": {
            "seq_data": ('r2',),
            "index": ('r1',)
        },
        "qc_modules": [
            'fastqc',
            'fastq_screen',
            'sequence_lengths',
            'strandedness'
        ]
    },

    "ParseEvercode": {
        "reads": {
            "seq_data": ('r1',),
            "index": ('r2',)
        },
        "qc_modules": [
            'fastqc',
            'fastq_screen',
            'sequence_lengths',
            'strandedness'
        ]
    },

    "ICELL8_scATAC": {
        "reads": {
            "seq_data": ('r1','r2',),
            "index": ()
        },
        "qc_modules": [
            'fastqc',
            'fastq_screen',
            'sequence_lengths',
            'strandedness'
        ]
    },
}

#######################################################################
# Functions
#######################################################################

def determine_qc_protocol(project):
    """
    Determine the QC protocol for a project

    Arguments:
      project (AnalysisProject): project instance

    Return:
      String: QC protocol for the project
    """
    # Standard protocols
    if project.info.paired_end:
        protocol = "standardPE"
    else:
        protocol = "standardSE"
    # Single cell protocols
    if project.info.single_cell_platform is not None:
        # Default
        protocol = "singlecell"
        single_cell_platform = project.info.single_cell_platform
        library_type = project.info.library_type
        if single_cell_platform.startswith('10xGenomics Chromium 3\''):
            if library_type == "scRNA-seq":
                # 10xGenomics scATAC-seq
                protocol = "10x_scRNAseq"
            elif library_type == "snRNA-seq":
                # 10xGenomics snRNA-seq
                protocol = "10x_snRNAseq"
            elif library_type in ("CellPlex",
                                  "CellPlex scRNA-seq",
                                  "CellPlex snRNA-seq"):
                # 10xGenomics CellPlex (cell multiplexing)
                protocol = "10x_CellPlex"
        elif single_cell_platform == 'Parse Evercode':
            if library_type == "scRNA-seq":
                # Parse Evercode snRNAseq
                protocol = "ParseEvercode"
        elif library_type in ("scATAC-seq",
                              "snATAC-seq",):
            if single_cell_platform == "10xGenomics Single Cell ATAC":
                # 10xGenomics scATAC-seq
                protocol = "10x_scATAC"
            elif single_cell_platform == "ICELL8":
                # ICELL8 scATAC-seq
                protocol = "ICELL8_scATAC"
    # Spatial RNA-seq
    if project.info.single_cell_platform == "10xGenomics Visium":
        # 10xGenomics Visium spatial transcriptomics
        protocol = "10x_Visium"
    # Multiome ATAC+GEX
    if project.info.single_cell_platform == "10xGenomics Single Cell Multiome":
        if library_type == "ATAC":
            # 10xGenomics single cell Multiome ATAC
            protocol = "10x_Multiome_ATAC"
        elif library_type == "GEX":
            # 10xGenomics single cell Multiome gene expression
            protocol = "10x_Multiome_GEX"
    return protocol

def fetch_protocol_definition(name):
    """
    Return the definition for a QC protocol

    Arguments:
      name (str): name of the QC protocol

    Returns:
      Tuple: definition as a tuple of the form
        (reads,qc_modules) where 'reads' is an
        AttributeDictionary with elements 'seq_data',
        'index', and 'qc' (listing sequence data,
        index reads, and all reads for QC,
        respectively) and 'qc_modules' is a list
        of QC module definitions.
    """
    if name not in QC_PROTOCOLS:
        raise KeyError("%s: undefined QC protocol" % name)
    protocol_defn = QC_PROTOCOLS[name]
    reads = AttributeDictionary()
    try:
        reads['seq_data'] = list(protocol_defn['reads']['seq_data'])
        reads['index'] = list(protocol_defn['reads']['index'])
        reads['qc'] = sorted(reads.seq_data + reads.index)
        qc_modules = [m for m in protocol_defn['qc_modules']]
    except KeyError as ex:
        raise Exception("%s: exception loading QC protocol "
                        "definition: %s" % (name,ex))
    return (reads,qc_modules)

def get_read_numbers(protocol):
    """
    Return the read numbers for a QC protocol definition

    Given a QC protocol, returns the integer read numbers
    for sequence data reads, index reads, and QC reads
    associated with that protocol.

    Arguments:
      protocol (str): name of the QC protocol

    Returns:
      AttributeDictionary: dictionary with keys 'seq_data',
        'index' and 'qc', mapping to lists of read numbers
        for each type of read data.
    """
    reads,qc_modules = fetch_protocol_definition(protocol)
    read_numbers = AttributeDictionary(seq_data=[],
                                       index=[],
                                       qc=[])
    for name in ('seq_data','index','qc'):
        for read in reads[name]:
            if read.startswith('r'):
                read_numbers[name].append(int(read[1:]))
    return read_numbers
