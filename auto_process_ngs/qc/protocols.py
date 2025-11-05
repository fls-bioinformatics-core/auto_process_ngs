#!/usr/bin/env python
#
#     protocols: define and handle QC protocols
#     Copyright (C) University of Manchester 2022-2025 Peter Briggs
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
- determine_qc_protocol_from_metadata: determine built-in QC protocol
- determine_qc_protocol: determine built-in protocol for a project
- fetch_protocol_definition: get the definition for a QC protocol
- parse_protocol_repr: get a QCProtocol object from a string
- parse_qc_module_spec: process QC module specification string
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

# QC modules
from .qc_modules import QC_MODULES
from .qc_modules import QC_MODULE_NAMES

# QC protocol definitions

QC_PROTOCOLS = {

    "minimal": {
        "description": "Minimal generic QC (all non-index Fastqs)",
        "reads": {
            "seq_data": ('r*',),
            "index": ()
        },
        "qc_modules": [
            'fastqc',
            'fastq_screen',
            'sequence_lengths'
        ]
    },

    "standard": {
        "description": "Standard QC for single and paired-end data",
        "reads": {
            "seq_data": ('r*',),
            "index": ()
        },
        "qc_modules": [
            'fastqc',
            'fastq_screen',
            'sequence_lengths',
            'picard_insert_size_metrics',
            'rseqc_genebody_coverage',
            'rseqc_infer_experiment',
            'qualimap_rnaseq'
        ]
    },

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
            'rseqc_genebody_coverage',
            'rseqc_infer_experiment',
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
            'picard_insert_size_metrics',
            'rseqc_genebody_coverage',
            'rseqc_infer_experiment',
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
            'rseqc_genebody_coverage',
            'rseqc_infer_experiment',
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
            'rseqc_genebody_coverage',
            'rseqc_infer_experiment',
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
            'picard_insert_size_metrics',
            'rseqc_genebody_coverage',
            'rseqc_infer_experiment',
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
            'rseqc_genebody_coverage',
            'rseqc_infer_experiment',
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
            'picard_insert_size_metrics',
            'rseqc_genebody_coverage',
            'rseqc_infer_experiment',
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
            'rseqc_genebody_coverage',
            'rseqc_infer_experiment',
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
            'rseqc_genebody_coverage',
            'rseqc_infer_experiment',
            'qualimap_rnaseq',
            'cellranger_multi'
        ]
    },

    "10x_ImmuneProfiling": {
        "description": "10xGenomics single cell immune profiling data",
        "reads": {
            "seq_data": ('r2',),
            "index": ('r1',)
        },
        "qc_modules": [
            'fastqc',
            'fastq_screen',
            'sequence_lengths',
            'rseqc_genebody_coverage',
            'rseqc_infer_experiment',
            'qualimap_rnaseq',
            'cellranger_count(cellranger_use_multi_config=True;'
                             'set_cell_count=false;'
                             'set_metadata=False)',
            'cellranger_multi(cellranger_required_version=">=9")'
        ]
    },

    "10x_Visium_GEX": {
        "description": "10xGenomics Visium spatial gene expression",
        "reads": {
            # 50bp insert in R2
            "seq_data": ('r2:1-50',),
            "index": ('r1',)
        },
        "qc_modules": [
            'fastqc',
            'fastq_screen',
            'sequence_lengths',
            'rseqc_genebody_coverage',
            'rseqc_infer_experiment',
            'qualimap_rnaseq'
        ]
    },

    "10x_Visium_GEX_90bp_insert": {
        "description": "10xGenomics Visium spatial gene expression "\
        "(90bp insert R2)",
        "reads": {
            # 90bp insert in R2
            "seq_data": ('r2:1-90',),
            "index": ('r1',)
        },
        "qc_modules": [
            'fastqc',
            'fastq_screen',
            'sequence_lengths',
            'rseqc_genebody_coverage',
            'rseqc_infer_experiment',
            'qualimap_rnaseq'
        ]
    },

    "10x_Visium_PEX": {
        "description": "10xGenomics Visium spatial protein expression",
        "reads": {
            # 50bp insert in R2
            "seq_data": ('r2:1-50',),
            "index": ('r1',)
        },
        "qc_modules": [
            'fastqc',
            'sequence_lengths'
        ]
    },

    "10x_Visium_legacy": {
        "description": "10xGenomics Visium spatial RNA-seq (legacy QC)",
        "reads": {
            "seq_data": ('r2',),
            "index": ('r1',)
        },
        "qc_modules": [
            'fastqc',
            'fastq_screen',
            'sequence_lengths',
            'rseqc_genebody_coverage',
            'rseqc_infer_experiment',
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
            'rseqc_genebody_coverage',
            'rseqc_infer_experiment',
            'qualimap_rnaseq'
        ]
    },

    "BioRad_ddSEQ_ATAC": {
        "description": "Bio-Rad ddSEQ Single Cell ATAC data",
        "reads": {
            "seq_data": ('r1:79-118', 'r2',),
            "index": ()
        },
        "qc_modules": [
            'fastqc',
            'fastq_screen',
            'sequence_lengths',
            'rseqc_genebody_coverage',
            'rseqc_infer_experiment',
            'qualimap_rnaseq'
        ]
    },

    "singlecell": {
        "description": "Generic single cell RNA-seq",
        "reads": {
            "seq_data": ('r2',),
            "index": ('r1',)
        },
        "qc_modules": [
            'fastqc',
            'fastq_screen',
            'sequence_lengths',
            'rseqc_genebody_coverage',
            'qualimap_rnaseq'
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
        self.seq_data_reads = []
        self.index_reads = []
        self.qc_modules = []
        self.update(seq_data_reads=seq_data_reads,
                    index_reads=index_reads,
                    qc_modules=qc_modules)

    def update(self,seq_data_reads=None,index_reads=None,
                 qc_modules=None):
        """
        Update the reads and QC modules for the protocol

        Allows the sequence data reads, index reads and QC
        modules associated with the protocol to be updated.
        Checks that values are valid and that internal data
        is also correctly updated.

        Arguments:
          seq_data_reads (list): read names associated
            with sequence data (if not supplied then
            existing read data will be kept)
          index_reads (list): read names associated
            with index data (if not supplied then existing
            read data will be kept)
          qc_modules (list): list of names of associated
            QC modules (if not supplied then existing
            modules data will be kept)
        """
        # Store QC modules
        if qc_modules is not None:
            self.qc_modules = (sorted([m for m in qc_modules])
                               if qc_modules is not None else [])
            # Check QC modules are valid
            for m in self.qc_modules:
                name = m.split('(')[0]
                if name not in QC_MODULE_NAMES:
                    raise QCProtocolError("'%s': unrecognised QC module"
                                          % name)
        # Store supplied reads
        if seq_data_reads is not None:
            self.seq_data_reads = [str(r).lower() for r in seq_data_reads] \
                                  if seq_data_reads else []
        if index_reads is not None:
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
        if not self.seq_data_reads:
            seq_data_reads = tuple()
        else:
            seq_data_reads = self.seq_data_reads
        if not index_reads:
            index_reads = tuple()
        else:
            index_reads = self.index_reads
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

    @property
    def expected_outputs(self):
        """
        Return a list of the expected QC outputs

        The expected outputs are based on the
        QC modules plus the reads associated
        with the protocol.
        """
        expected = []
        for qc_module in self.qc_module_names:
            if qc_module == "fastqc":
                for r in self.reads.qc:
                    expected.append("fastqc_%s" % r)
            elif qc_module == "fastq_screen":
                for r in self.reads.seq_data:
                    expected.append("screens_%s" % r)
            else:
                expected.append(qc_module)
        return sorted(expected)

    @property
    def qc_module_names(self):
        """
        Return list of QC module names without parameters

        Returns a list of the QC modules associated with
        the protocol, with any parameter lists (i.e.
        trailing '(...)') removed.

        For example the QC module list:

        ['cellranger(use_10x_multi_config=true)','fastqc']

        will be returned as:

        ['cellranger','fastqc']
        """
        return [parse_qc_module_spec(m)[0] for m in self.qc_modules]

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
                try:
                    read_numbers.append(int(r[1:]))
                except ValueError:
                    # Check for wildcard read number
                    if r[1:] == "*":
                        read_numbers.append("*")
                    else:
                        raise QCProtocolError(f"'{r}': invalid read "
                                              "specification")
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
        for m in QC_MODULES:
            if m.name in self.qc_module_names and m.mapped_metrics:
                mapped_metrics.append(m.name)
        return sorted(mapped_metrics)

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

def determine_qc_protocol_from_metadata(library_type,
                                        single_cell_platform,
                                        paired_end):
    """
    Determine the QC protocol from metadata values

    Arguments:
      library_type (str): library or application
      single_cell_platform (str): single cell platform (or None)
      paired_end (bool): whether data are paired end

    Return:
      String: QC protocol for the project
    """
    # Standard protocols
    if library_type is None or not library_type:
        protocol = "minimal"
    elif library_type in ("DNA-seq", "WGS") or \
         library_type.startswith("CRISPR"):
        protocol = "minimal"
    elif paired_end:
        protocol = "standardPE"
    else:
        protocol = "standardSE"
    # Single cell protocols
    if single_cell_platform is not None:
        # 10x Genomics
        if single_cell_platform.startswith('10xGenomics Chromium 3\'') or \
           single_cell_platform.startswith('10xGenomics Chromium GEM-X') or \
           single_cell_platform.startswith('10xGenomics Chromium Next GEM'):
            # Default for 10x
            protocol = "singlecell"
            # 10xGenomics scRNA-seq
            if library_type == "scRNA-seq":
                protocol = "10x_scRNAseq"
            # 10xGenomics snRNA-seq
            elif library_type == "snRNA-seq":
                protocol = "10x_snRNAseq"
            # 10xGenomics CellPlex (cell multiplexing)
            elif library_type in ("CellPlex",
                                  "CellPlex scRNA-seq",
                                  "CellPlex snRNA-seq"):
                protocol = "10x_CellPlex"
            # 10xGenomics Flex (fixed RNA profiling)
            elif library_type == "Flex":
                protocol = "10x_Flex"
        # 10x single cell immune profiling
        elif single_cell_platform.startswith('10xGenomics Chromium 5\''):
            if library_type == "Single Cell Immune Profiling":
                protocol = "10x_ImmuneProfiling"
        # 10x ATAC
        elif single_cell_platform == "10xGenomics Single Cell ATAC":
            if library_type in ("scATAC-seq",
                                "snATAC-seq",):
                # 10xGenomics scATAC-seq
                protocol = "10x_scATAC"
        # 10x Multiome ATAC+GEX
        elif single_cell_platform == "10xGenomics Single Cell Multiome":
            if library_type == "ATAC":
                # 10xGenomics single cell Multiome ATAC
                protocol = "10x_Multiome_ATAC"
            elif library_type == "GEX":
                # 10xGenomics single cell Multiome gene expression
                protocol = "10x_Multiome_GEX"
        # Parse Evercode
        elif single_cell_platform == 'Parse Evercode':
            if library_type in ("scRNA-seq",
                                "snRNA-seq",
                                "TCR",
                                "TCR scRNA-seq",
                                "WT",
                                "WT scRNA-seq"):
                # Parse Evercode snRNAseq
                protocol = "ParseEvercode"
        # Bio-Rad RNA-Seq
        elif single_cell_platform == "Bio-Rad ddSEQ Single Cell 3' RNA-Seq":
            if library_type in ("scRNA-seq",
                                "snRNA-seq",):
                # Bio-Rad ddSeq RNA-Seq
                protocol = "minimal"
        # Bio-Rad ATAC
        elif single_cell_platform == "Bio-Rad ddSEQ Single Cell ATAC":
            if library_type in ("scATAC-seq",
                                "snATAC-seq",):
                # Bio-Rad ddSeq ATAC
                protocol = "BioRad_ddSEQ_ATAC"
        # Visium/spatial data
        elif single_cell_platform in ("10xGenomics Visium",
                                      "10xGenomics Visium (CytAssist)",
                                      "10xGenomics CytAssist Visium"):
            if library_type.lower() in ("ffpe spatial pex",
                                        "ffpe spatial protein expression"):
                # Spatial protein expression
                protocol = "10x_Visium_PEX"
            elif library_type.lower() in \
                 ("fresh frozen spatial gene expression",
                  "fresh frozen spatial gex") \
                  and single_cell_platform == "10xGenomics Visium":
                # Special case (GEX with 90bp insert in R2)
                protocol = "10x_Visium_GEX_90bp_insert"
            elif library_type.lower() == "spatial rna-seq" and \
                 single_cell_platform in ("10xGenomics Visium",
                                          "10xGenomics CytAssist Visium"):
                # Legacy spatial RNA-seq
                protocol = "10x_Visium_legacy"
            else:
                # Default (GEX with 50bp insert in R2)
                protocol = "10x_Visium_GEX"
    return protocol

def determine_qc_protocol(project):
    """
    Determine the QC protocol for a project

    Arguments:
      project (AnalysisProject): project instance

    Return:
      String: QC protocol for the project
    """
    return determine_qc_protocol_from_metadata(
        project.info.library_type,
        project.info.single_cell_platform,
        project.info.paired_end)

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

def parse_qc_module_spec(module_spec):
    """
    Parse QC module spec into name and parameters

    Parse a QC module specification of the form
    ``NAME`` or ``NAME(KEY=VALUE;...)`` and return
    the module name and any additional parameters
    in the form of a dictionary.

    For example:

    >>> parse_qc_module_spec('NAME')
    ('NAME', {})
    >>> parse_qc_module_spec('NAME(K1=V1;K2=V2)')
    ('NAME', { 'K1':'V1', 'K2':'V2' })

    By default values are returned as strings (with
    surrounding single or double quotes removed);
    however basic type conversion is also applied to
    certain values:

    - True/true and False/false are returned as the
      appropriate boolean value

    Arguments:
      module_spec (str): QC module specification

    Returns:
      Tuple: tuple of the form (name,params) where
        'name' is the QC module name and 'params'
        is a dictionary with the extracted key-value
        pairs.
    """
    # Handle module specification string of the form
    # 'NAME[(KEY=VALUE;...)]'
    items = module_spec.split('(')
    # Extract the module name and associated parameter list
    name = items[0]
    params = {}
    # Extract the key-value pairs from the parameters
    try:
        current_key = None
        current_value = None
        param_str = items[1].rstrip(')')
        value_quote = None
        # Need to handle string character-by-character
        for c in param_str:
            if current_key is None:
                # Start of a new key
                current_key = c
            elif current_value is None:
                if c == "=":
                    # End of key, start of a new value
                    current_value = ""
                else:
                    # Still processing the key
                    current_key += c
            elif not current_value and not value_quote:
                # First character in new value
                if c in ("'", "\""):
                    # Value is quoted
                    value_quote = c
                # Store full value (including quotes)
                current_value = c
            else:
                # Inside a value string
                if value_quote:
                    # Inside a quoted string
                    if c == value_quote:
                        # End of the quoted string
                        value_quote = None
                    # Append to current value (including quotes)
                    current_value += c
                elif c == ";":
                    # End of value, store extracted
                    # key-value pair
                    params[current_key] = current_value
                    # Reset for next pair
                    current_key = None
                    current_value = None
                else:
                    # Still part of current value
                    current_value += c
        # Store any final key-value pair
        if current_key and current_value:
            params[current_key] = current_value
        else:
            # Encountered an error
            raise Exception(f"Error parsing parameters: {param_str}")
    except IndexError:
        # No closing ")"
        pass
    # Post-process values
    for key in params:
        value = params[key]
        if value[0] in ("'", "\"") and value[0] == value[-1]:
            # Strip quotes
            value = value[1:-1]
        elif value in ("True", "true"):
            # Boolean true
            value = True
        elif value in ("False", "false"):
            # Boolean false
            value = False
        params[key] = value
    return (name, params)

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
