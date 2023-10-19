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
- determine_qc_protocol: determine built-in protocol for a project
- fetch_protocol_definition: get the definition for a QC protocol
- parse_protocol_repr: get a QCProtocol object from a string
"""

#######################################################################
# Imports
#######################################################################

import logging
from bcftbx.utils import AttributeDictionary

# Module specific logger
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

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
            "seq_data": ('r2:1-50',),
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

    "10x_Visium_FFPE": {
        "description": "10xGenomics Visium FFPE spatial RNA-seq",
        "reads": {
            "seq_data": ('r2:1-50',),
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
      and all reads for QC, respectively)
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

    For example: 'r1:1-50'. These can be accessed via the
    'read_range' property as e.g. 'read_range.r1' and will
    be returned as either 'None' (if no range was supplied),
    or a tuple '(START,END)' (where either 'START' or 'END'
    will be 'None' if no limit was supplied).

    The 'seq_data_reads' and 'index_reads' properties
    store the original read specifications.

    QCProtocol instances can also be created directly from
    protocol specification strings using the
    'from_specification' class method. (Specification
    strings are returned from 'repr' on existing QCProtocol
    instances, or can alternatively be constructed
    manually.)

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
        # Store supplied reads
        self.seq_data_reads = [str(r).lower() for r in seq_data_reads] \
                              if seq_data_reads else []
        self.index_reads = [str(r).lower() for r in index_reads] \
                           if index_reads else []
        # Normalise and store supplied read names
        self.reads = AttributeDictionary(
            seq_data=self.__reads(self.seq_data_reads),
            index=self.__reads(self.index_reads))
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

    @classmethod
    def from_specification(cls,s):
        """
        Create new QCProtocol instance from specification

        Given a specification string (such as that
        returned by 'repr(...)'), create a new
        QCProtocol instance initialised from that
        specification.

        Example usage:

        >>> p = QCProtocol.from_specification("custom:...")

        Arguments:
          s (str): QC protocol specification string
        """
        return cls(**parse_protocol_spec(s))

    def summarise(self):
        """
        Summarise protocol

        Generate plain-text description of the protocol
        """
        summary = []
        if self.reads.seq_data:
            summary.append("biological data in %s" %
                           self.__summarise_reads(self.reads.seq_data))
        else:
            summary.append("no reads explicitly assigned as biological data")
        if self.reads.index:
            summary.append("index data in %s" %
                           self.__summarise_reads(self.reads.index))
        else:
            summary.append("no reads explicitly assigned as index data")
        if self.__mapped_metrics():
            has_ranges = False
            if self.reads.seq_data:
                for rd in self.reads.seq_data:
                    if self.read_range[rd]:
                        has_ranges = True
                        break
            if has_ranges:
                summary.append("mapped metrics generated using only "
                               "subsequences of biological data reads")
            else:
                summary.append("mapped metrics generated using only "
                               "biological data reads")
        summary = '; '.join(summary)
        if self.name:
            summary = "'%s' protocol: %s" % (self.name,summary)
        return summary

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

    def __summarise_reads(self,rds):
        # Internal: generate a plain text description of the
        # sequence data and index reads and ranges
        reads = []
        for r in rds:
            rd = r.upper()
            if r in self.read_range:
                rng = self.read_range[r]
                if rng and (rng[0] or rng[1]):
                    if not rng[0]:
                        rd += " (%s first %sbp)" % (rd,rng[1])
                    elif not rng[1]:
                        rd += " (%s from %s base)" % (rd,rng[0])
                    else:
                        rd += " (%s bases %s to %s)" % (rd,rng[0],rng[1])
                reads.append(rd)
        if len(reads) == 1:
            reads[0] += " only"
        elif len(reads) > 1:
            reads[-2] += " and %s" % reads[-1]
            reads = reads[0:-1]
        return ', '.join(reads)

    def __mapped_metrics(self):
        # Internal: return list of QC module names within
        # the protocol which produce mapped metrics
        mapped_metrics = []
        for qc_defn in self.qc_modules:
            qc_module = qc_defn.split('(')[0]
            if qc_module in ('fastq_screen',
                             'strandedness',
                             'picard_insert_size_metrics',
                             'rseqc_genebody_coverage',
                             'qualimap_rnaseq'):
                mapped_metrics.append(qc_module)
        return mapped_metrics

    def __repr_reads(self,rds):
        # Internal: get string representation of reads
        # for __repr__
        reads = []
        for rd in rds:
            rng = self.read_range[rd]
            if rng and (rng[0] or rng[1]):
                if not rng[0]:
                    reads.append("%s:1-%s" % (rd,rng[1]))
                elif not rng[1]:
                    reads.append("%s:%s-" % (rd,rng[0]))
                else:
                    reads.append("%s:%s-%s" % (rd,rng[0],rng[1]))
            else:
                reads.append(rd)
        return "[%s]" % ','.join(r for r in reads)

    def __repr__(self):
        # Generate string representation
        qc_modules = "[%s]" % ','.join(m for m in self.qc_modules)
        return \
            "{name}:'{descr}':seq_reads={seq_reads}:"\
            "index_reads={index_reads}:"\
            "qc_modules={qc_modules}".format(
                name=self.name,
                descr=self.description,
                seq_reads=self.__repr_reads(self.reads.seq_data),
                index_reads=self.__repr_reads(self.reads.index),
                qc_modules=qc_modules)

    def __eq__(self,p):
        # Check if another object is equal to this one
        return (isinstance(p,QCProtocol) and
                repr(self) == repr(p))

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
            if library_type in ("scRNA-seq",
                                "snRNA-seq"):
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
    if project.info.single_cell_platform in ("10xGenomics Visium",
                                             "10xGenomics CytAssist Visium"):
        # 10xGenomics Visium spatial transcriptomics
        if project.info.library_type == "FFPE Spatial RNA-seq":
            protocol = "10x_Visium_FFPE"
        else:
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

def fetch_protocol_definition(p):
    """
    Return the definition for a QC protocol

    Fetches a QCProtocol instance for the supplied
    protocol, which can either be a protocol
    definition string or the name of a built-in
    QC protocol.

    Arguments:
      p (str): name or specification for the QC
       protocol

    Returns:
      QCProtocol: QCProtocol object representing
        the requested protocol.

    Raises:
      KeyError: when a QCProtocol instance can't
        be returned for the requested protocol.
    """
    # Try protocol specification
    if ':' in p:
        return QCProtocol.from_specification(p)
    # Try looking up a built-in protocol
    if p not in QC_PROTOCOLS:
        raise KeyError("%s: not built-in QC protocol" % p)
    protocol_defn = QC_PROTOCOLS[p]
    return QCProtocol(name=p,
                      description=protocol_defn['description'],
                      seq_data_reads=protocol_defn['reads']['seq_data'],
                      index_reads=protocol_defn['reads']['index'],
                      qc_modules=protocol_defn['qc_modules'])

def parse_protocol_spec(s):
    """
    Parse QC protocol specification string

    Parses a QC protocol specification string (such as
    one returned by the '__repr__' built-in of an
    existing QCProtocol instance) and returns an
    AttributeDictionary with the following elements
    extracted from the specification:

    - name
    - description
    - seq_reads
    - index_reads
    - qc_modules

    These can then be used to create a new QCProtocol
    instance which matches the specification using
    e.g.

    >>> p = QCProtocol(**parse_protocol_spec("..."))

    Arguments:
      s (string): QC protocol specification string

    Returns:
      AttributeDictionary: AttributeDictionary with
        keys mapped to values from the supplied
        specification.

    Raises:
      QCProtocolParseSpecError: if the specification
        string cannot be parsed correctly.
    """
    # Initialise
    seq_data_reads = None
    index_reads = None
    qc_modules = None
    # Extract name
    try:
        ix = s.index(':')
        name = s[:ix]
        s = s[ix+1:]
    except ValueError:
        name = s
        s = ""
    logger.debug("Name: %s" % name)
    logger.debug("(Remaining string: %s)" % s)
    if not s:
        raise QCProtocolParseSpecError("Incomplete specification: "
                                       "nothing after 'name'?")
    # Extract description
    if s.startswith('"') or s.startswith("'"):
        quote = s[0]
        s = s[1:]
        try:
            ix = s.index(quote)
            if s[ix+1] != ':':
                raise QCProtocolParseSpecError("Unable to parse "
                                               "description: continues "
                                               "after closing quote?")
            description = s[:ix]
            s = s[ix+2:]
        except ValueError:
            raise QCProtocolParseSpecError("Unable to parse description: "
                                           "no closing quote?")
    else:
        try:
            ix = s.index(':')
            description = s[:ix]
            s = s[ix:]
        except ValueError:
            description = s
            s = ""
    logger.debug("Description: %s" % description)
    logger.debug("(Remaining string: %s)" % s)
    if not s:
        raise QCProtocolParseSpecError("Incomplete specification: "
                                       "nothing after description?")
    # Basic checks on remainder of string
    if s.count("[") != s.count("]"):
        raise QCProtocolParseSpecError("Specification syntax error: "
                                       "unmatched braces?")
    # Break up remainder of the string
    while s:
        if s.startswith(':'):
            # Strip leading colon
            s = s[1:]
        if s.startswith("seq_reads=["):
            # Sequencing data reads
            ix = s.index(']')
            seq_data_reads = list(
                filter(lambda r: r != '',
                       s[len("seq_reads=["):ix].split(',')))
            s = s[ix+1:]
            logger.debug("Seq data reads: %s" % seq_data_reads)
            logger.debug("(Remaining string: %s)" % s)
        elif s.startswith("index_reads=["):
            # Index reads
            try:
                ix = s.index(']')
                index_reads = list(
                    filter(lambda r: r != '',
                           s[len("index_reads=["):ix].split(',')))
                s = s[ix+1:]
            except ValueError:
                raise QCProtocolParseSpecError("Unable to parse index "
                                               "reads: no closing brace?")
            logger.debug("Index reads: %s" % index_reads)
            logger.debug("Remaining string: %s" % s)
        elif s.startswith("qc_modules=["):
            # QC modules
            try:
                ix = s.index(']')
                qc_modules = s[len("qc_modules=["):ix].split(',')
                s = s[ix+1:]
            except ValueError:
                qc_modules = s[len("qc_modules=["):].split(',')
                s = ""
            qc_modules = list(filter(lambda m: m != '',qc_modules))
            logger.debug("QC modules: %s" % qc_modules)
            logger.debug("Remaining string: %s" % s)
        elif s:
            # Unknown component
            raise QCProtocolParseSpecError("Unable to parse section "
                                           "starting '%s...': "
                                           "unknown specification "
                                           "element?" % s[:5])
    # Return mapping of extracted components
    return AttributeDictionary(name=name,
                               description=description,
                               seq_data_reads=seq_data_reads,
                               index_reads=index_reads,
                               qc_modules=qc_modules)

######################################################################
# Custom exceptions
######################################################################

class QCProtocolError(Exception):
    """
    Base class for QC protocol-specific exceptions
    """

class QCProtocolParseSpecError(QCProtocolError):
    """
    Exception raised when a protocol specificaion
    can't be parsed

    Arguments:
      message (str): error message
    """
    def __init__(self,message=None):
        if message is None:
            message = "Error parsing QC protocol specification"
        self.message = message
        QCProtocolError.__init__(self,self.message)
