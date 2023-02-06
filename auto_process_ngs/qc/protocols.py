#!/usr/bin/env python
#
#     protocols: define and handle QC protocols
#     Copyright (C) University of Manchester 2022-2023 Peter Briggs
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
        "description": "Example QC protocol"
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

This module also provides the following classes and functions:

- QCProtocol: class representing a QC protocol
- determine_qc_protocol: get QC protocol for a project
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
        "description": "Standard single-end data (R1 Fastqs only)",
        "reads": {
            "seq_data": ('r1',),
            "index": ()
        },
        "qc_modules": [
            'fastqc',
            'fastq_screen',
            'sequence_lengths',
            'strandedness',
            'rseqc_genebody_coverage',
            'qualimap_rnaseq'
        ]
    },

    "standardPE": {
        "description": "Standard paired-end data (R1/R2 Fastq pairs)",
        "reads": {
            "seq_data": ('r1','r2',),
            "index": ()
        },
        "qc_modules": [
            'fastqc',
            'fastq_screen',
            'sequence_lengths',
            'strandedness',
            'picard_insert_size_metrics',
            'rseqc_genebody_coverage',
            'qualimap_rnaseq'
        ]
    },

    "singlecell": {
        "description": "ICELL8 single cell RNA-seq",
        "reads": {
            "seq_data": ('r2',),
            "index": ('r1',)
        },
        "qc_modules": [
            'fastqc',
            'fastq_screen',
            'sequence_lengths',
            'strandedness',
            'rseqc_genebody_coverage',
            'qualimap_rnaseq'
        ]
    },

    "10x_scRNAseq": {
        "description": "10xGenomics single cell RNA-seq",
        "reads": {
            "seq_data": ('r2',),
            "index": ('r1',)
        },
        "qc_modules": [
            'fastqc',
            'fastq_screen',
            'sequence_lengths',
            'strandedness',
            'rseqc_genebody_coverage',
            'qualimap_rnaseq',
            'cellranger_count'
        ]
    },

    "10x_snRNAseq": {
        "description": "10xGenomics single nuclei RNA-seq",
        "reads": {
            "seq_data": ('r2',),
            "index": ('r1',)
        },
        "qc_modules": [
            'fastqc',
            'fastq_screen',
            'sequence_lengths',
            'strandedness',
            'rseqc_genebody_coverage',
            'qualimap_rnaseq',
            'cellranger_count'
        ]
    },

    "10x_scATAC": {
        "description": "10xGenomics single cell ATAC-seq",
        "reads": {
            "seq_data": ('r1','r3',),
            "index": ()
        },
        "qc_modules": [
            'fastqc',
            'fastq_screen',
            'sequence_lengths',
            'strandedness',
            'picard_insert_size_metrics',
            'rseqc_genebody_coverage',
            'qualimap_rnaseq',
            'cellranger-atac_count'
        ]
    },

    "10x_Multiome_GEX": {
        "description": "10xGenomics single cell multiome gene expression data",
        "reads": {
            "seq_data": ('r2',),
            "index": ('r1',)
        },
        "qc_modules": [
            'fastqc',
            'fastq_screen',
            'sequence_lengths',
            'strandedness',
            'rseqc_genebody_coverage',
            'qualimap_rnaseq',
            'cellranger-arc_count',
            # Also run Cellranger
            # Set the chemistry to 'ARC-v1' and library to 'snRNA-seq'
            # See https://kb.10xgenomics.com/hc/en-us/articles/360059656912
            'cellranger_count(chemistry=ARC-v1;'
                             'library=snRNA-seq;'
                             'cellranger_version=*;'
                             'cellranger_refdata=*;'
                             'set_cell_count=false;'
                             'set_metadata=False)'
        ]
    },

    "10x_Multiome_ATAC": {
        "description": "10xGenomics single cell multiome ATAC-seq data",
        "reads": {
            "seq_data": ('r1','r3',),
            "index": ()
        },
        "qc_modules": [
            'fastqc',
            'fastq_screen',
            'sequence_lengths',
            'strandedness',
            'picard_insert_size_metrics',
            'rseqc_genebody_coverage',
            'qualimap_rnaseq',
            'cellranger-arc_count',
            # Also run Cellranger ATAC
            # Set the chemistry to 'ARC-v1' and library to 'scATAC-seq'
            # See https://kb.10xgenomics.com/hc/en-us/articles/360061165691
            'cellranger-atac_count(chemistry=ARC-v1;'
                                  'library=scATAC-seq;'
                                  'cellranger_version=*;'
                                  'cellranger_refdata=*;'
                                  'set_cell_count=false;'
                                  'set_metadata=False)'
        ]
    },

    "10x_CellPlex": {
        "description": "10xGenomics CellPlex cell multiplexing data",
        "reads": {
            "seq_data": ('r2',),
            "index": ('r1',)
        },
        "qc_modules": [
            'fastqc',
            'fastq_screen',
            'sequence_lengths',
            'strandedness',
            'rseqc_genebody_coverage',
            'qualimap_rnaseq',
            'cellranger_count(cellranger_use_multi_config=True;'
                             'set_cell_count=false;'
                             'set_metadata=False)',
            'cellranger_multi'
        ]
    },

    "10x_Flex": {
        "description": "10xGenomics fixed RNA profiling (Flex) data",
        "reads": {
            "seq_data": ('r2',),
            "index": ('r1',)
        },
        "qc_modules": [
            'fastqc',
            'fastq_screen',
            'sequence_lengths',
            'strandedness',
            'rseqc_genebody_coverage',
            'qualimap_rnaseq',
            'cellranger_multi'
        ]
    },

    "10x_Visium": {
        "description": "10xGenomics Visium spatial RNA-seq",
        "reads": {
            "seq_data": ('r2',),
            "index": ('r1',)
        },
        "qc_modules": [
            'fastqc',
            'fastq_screen',
            'sequence_lengths',
            'strandedness',
            'rseqc_genebody_coverage',
            'qualimap_rnaseq'
        ]
    },

    "ParseEvercode": {
        "description": "Parse Biosciences Evercode data",
        "reads": {
            "seq_data": ('r1',),
            "index": ('r2',)
        },
        "qc_modules": [
            'fastqc',
            'fastq_screen',
            'sequence_lengths',
            'strandedness',
            'strandedness',
            'rseqc_genebody_coverage',
            'qualimap_rnaseq'
        ]
    },

    "ICELL8_scATAC": {
        "description": "ICELL8 single cell ATAC-seq",
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
# Classes
#######################################################################

class QCProtocol:
    """
    Class defining a QC protocol

    Properties:

    - name: protocol name
    - description: text description
    - reads: AttributeDictionary with elements 'seq_data',
      'index', and 'qc' (listing sequence data, index reads,
      and all read for QC, respectively)
    - read_numbers: AttributeDictionary with the same
      elements as 'reads', listing non-index read numbers
    - read_range: AttributeDictionary with normalised read
      names as elements and range of bases (as a tuple) as
      the values
    - qc_modules: list of QC module definitions

    Reads are supplied as 'r1', 'i2' etc; read numbers are
    integers without the leading 'r' (NB index reads are
    not included).

    Read ranges define subsequences within each read which
    contain the biologically significant data, and can be
    appended to the supplied reads using the syntax:

    READ[:[START]-[END]]

    For example: 'r1:1-50'. These can ben accessed via the
    'read_range' property as e.g. 'read_range.r1' and will
    be returned as either 'None' (if no range was supplied),
    or a tuple '(START,END)' (where either 'START' or 'END'
    will be 'None' if no limit was supplied).

    Arguments:
      name (str): name of the protocol
      description (str): protocol description
      seq_data_reads (list): read names associated
        with sequence data
      index_reads (list): read names associated with
        index data
      qc_modules (list): list of names of associated
        QC modules
    """
    def __init__(self,name,description,seq_data_reads,index_reads,
                 qc_modules):
        # Store name, description and modules
        self.name = str(name)
        self.description = (str(description) if description is not None
                            else "")
        self.qc_modules = (sorted([m for m in qc_modules])
                           if qc_modules is not None else [])
        # Normalise and store supplied read names
        self.reads = AttributeDictionary(
            seq_data=self.__reads(seq_data_reads),
            index=self.__reads(index_reads))
        # Generate and store QC read names
        self.reads['qc'] = self.__reads(self.reads.seq_data +
                                        self.reads.index)
        # Extract and store read numbers
        self.read_numbers = AttributeDictionary(
            seq_data=self.__read_numbers(self.reads.seq_data),
            index=self.__read_numbers(self.reads.index),
            qc=self.__read_numbers(self.reads.qc))
        # Extract and store sequence ranges for reads
        self.read_range = AttributeDictionary()
        if not seq_data_reads:
            seq_data_reads = tuple()
        if not index_reads:
            index_reads = tuple()
        for r in list(seq_data_reads) + list(index_reads):
            rd,rng = self.__parse_read_defn(r)
            self.read_range[rd] = rng

    def __parse_read_defn(self,read):
        # Internal: process a read definition string of the
        # form 'READ[:[START]-[END]]' and return a tuple
        # of (READ,(START,END)) or (READ,None) (if no range
        # was supplied)
        # Extract read and range
        try:
            rd,rng = read.split(':')
        except ValueError:
            rd,rng = (read,None)
        # Convert read to lower case
        rd = rd.lower()
        # Deal with range
        if rng:
            rng = tuple([int(s) for s in rng.split('-')])
        else:
            rng = None
        return (rd,rng)

    def __reads(self,reads):
        # Internal: normalise read names (remove range spec,
        # then convert to lowercase and sort)
        if not reads:
            return tuple()
        return tuple(sorted([self.__parse_read_defn(r)[0]
                             for r in reads]))

    def __read_numbers(self,reads):
        # Internal: extract read numbers (discard leading
        # character from read names and convert to integer)
        read_numbers = []
        for r in reads:
            if r.startswith('r'):
                read_numbers.append(int(r[1:]))
        return tuple(sorted(read_numbers))

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
            elif library_type == "Flex":
                # 10xGenomics Flex (fixed RNA profiling)
                protocol = "10x_Flex"
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
      QCProtocol: QCProtocol object representing
        the requested protocol

    Raises:
      KeyError: when the requested protocol isn't
        defined.
    """
    if name not in QC_PROTOCOLS:
        raise KeyError("%s: undefined QC protocol" % name)
    protocol_defn = QC_PROTOCOLS[name]
    return QCProtocol(name=name,
                      description=protocol_defn['description'],
                      seq_data_reads=protocol_defn['reads']['seq_data'],
                      index_reads=protocol_defn['reads']['index'],
                      qc_modules=protocol_defn['qc_modules'])
