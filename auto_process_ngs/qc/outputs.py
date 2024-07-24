#!/usr/bin/env python
#
#     qc/outputs: utilities to predict and check QC pipeline outputs
#     Copyright (C) University of Manchester 2019-2024 Peter Briggs
#
"""
Provides utility classes and functions for QC outputs.

Provides the following classes:

- QCOutputs: detect and characterise QC outputs
- ExtraOutputs: helper class for reading 'extra_outputs.tsv' file
"""

#######################################################################
# Imports
#######################################################################

import os
import logging
from bcftbx.qc.report import strip_ngs_extensions
from bcftbx.utils import AttributeDictionary
from ..analysis import AnalysisFastq
from ..analysis import split_sample_name
from ..fastq_utils import group_fastqs_by_name
from ..fastq_utils import remove_index_fastqs
from ..metadata import AnalysisProjectQCDirInfo
from ..tenx.cellplex import CellrangerMultiConfigCsv
from .modules import QCDir
from .modules.cellranger_arc_count import CellrangerArcCount
from .modules.cellranger_atac_count import CellrangerAtacCount
from .modules.cellranger_count import CellrangerCount
from .modules.cellranger_multi import CellrangerMulti
from .modules.fastqc import Fastqc
from .modules.fastq_screen import FastqScreen
from .modules.multiqc import Multiqc
from .modules.picard_insert_size_metrics import PicardInsertSizeMetrics
from .modules.qualimap_rnaseq import QualimapRnaseq
from .modules.rseqc_genebody_coverage import RseqcGenebodyCoverage
from .modules.rseqc_infer_experiment import RseqcInferExperiment
from .modules.sequence_lengths import SequenceLengths
from .modules.strandedness import Strandedness

# Module specific logger
logger = logging.getLogger(__name__)

######################################################################
# Supported QC module classes
######################################################################

QC_MODULES = (CellrangerCount,
              CellrangerAtacCount,
              CellrangerArcCount,
              CellrangerMulti,
              Fastqc,
              FastqScreen,
              Multiqc,
              PicardInsertSizeMetrics,
              QualimapRnaseq,
              RseqcGenebodyCoverage,
              RseqcInferExperiment,
              SequenceLengths,
              Strandedness)

#######################################################################
# Classes
#######################################################################

class QCOutputs:
    """
    Class to detect and characterise QC outputs

    On instantiation this class scans the supplied
    directory to identify and classify various
    artefacts which are typically produced by the
    QC pipeline.

    The following attributes are available:

    - fastqs: sorted list of Fastq names
    - reads: list of reads (e.g. 'r1', 'r2', 'i1' etc)
    - samples: sorted list of sample names extracted
      from Fastqs
    - seq_data_samples: sorted list of samples with
      biological data (rather than e.g. feature
      barcodes)
    - bams: sorted list of BAM file names
    - organisms: sorted list of organism names
    - fastq_screens: sorted list of screen names
    - cellranger_references: sorted list of reference
      datasets used with 10x pipelines
    - cellranger_probe_sets: sorted list of probe set
      files used with 10x pipelines
    - multiplexed_samples: sorted list of sample names
      for multiplexed samples (e.g. 10x CellPlex)
    - outputs: list of QC output categories detected (see
      below for valid values)
    - output_files: list of absolute paths to QC output
      files
    - software: dictionary with information on the
      QC software packages
    - stats: AttrtibuteDictionary with useful stats from
      across the project
    - config_files: list of QC configuration files found
      in the QC directory (see below for valid values)

    The following are valid values for the 'outputs'
    property:

    - 'fastqc_[r1...]'
    - 'screens_[r1...]'
    - 'strandedness'
    - 'sequence_lengths'
    - 'picard_insert_size_metrics'
    - 'rseqc_genebody_coverage'
    - 'rseqc_infer_experiment'
    - 'qualimap_rnaseq'
    - 'icell8_stats'
    - 'icell8_report'
    - 'cellranger_count'
    - 'cellranger_multi'
    - 'cellranger-atac_count'
    - 'cellranger-arc_count'
    - 'multiqc'
    - 'extra_outputs'

    The following are valid values for the 'config_files'
    property:

    - fastq_strand.conf
    - 10x_multi_config.csv
    - 10x_multi_config.<SAMPLE>.csv
    - libraries.<SAMPLE>.csv

    Arguments:
      qc_dir (str): path to directory to examine
      fastq_attrs (BaseFastqAttrs): (optional) class for
        extracting data from Fastq names
    """
    def __init__(self,qc_dir,fastq_attrs=None):
        # Directory to examine
        self.qc_dir = os.path.abspath(qc_dir)
        # How to handle Fastq names
        if fastq_attrs is not None:
            self.fastq_attrs = fastq_attrs
        else:
            self.fastq_attrs = AnalysisFastq
        # QC metadata
        self.qc_info = None
        qc_info_file = os.path.join(self.qc_dir,"qc.info")
        if os.path.exists(qc_info_file):
            self.qc_info = AnalysisProjectQCDirInfo(qc_info_file)
        # Properties
        self.fastqs = set()
        self.samples = []
        self.seq_data_samples = []
        self.bams = set()
        self.reads = []
        self.organisms = set()
        self.stats = AttributeDictionary(
            max_seqs=None,
            min_sequence_length_read={},
            max_sequence_length_read={},
            min_sequence_length=None,
            max_sequence_length=None,
        )
        self.software = {}
        self.outputs = set()
        self.output_files = []
        self.config_files = set()
        # Raw QC output data
        self._qc_data = {}
        # Scan the directory for QC outputs
        self._collect_qc_outputs()
        # Report
        print("Collected %d output files" % len(self.output_files))
        print("Outputs:")
        for output in self.outputs:
            print("- %s" % output)
        print("Fastqs:")
        for fq in self.fastqs:
            print("- %s" % fq)
        print("BAMs:")
        for b in self.bams:
            print("- %s" % b)
        print("Samples:")
        for s in self.samples:
            print("- %s" % s)
        print("Reads:")
        for r in self.reads:
            print("- %s" % r)
        print("Organisms:")
        for o in self.organisms:
            print("- %s" % o)
        print("Biological samples:")
        for s in self.seq_data_samples:
            print("- %s" % s)
        print("Software:")
        for sw in self.software:
            print("- %s: %s" % (sw,','.join(self.software[sw])))
        print("Screens:")
        for scrn in self.fastq_screens:
            print("- %s" % scrn)
        print("QC config files:")
        for cfg in self.config_files:
            print("- %s" % cfg)

    def data(self,name):
        """
        Return the 'raw' data associated with a QC output

        Arguments:
          name (str): name identifier for a QC output
            (e.g. 'fastqc')

        Returns:
          AttributeDictionary containing the raw data
            associated with the named QC output.

        Raises:
          KeyError: if the name doesn't match a stored
            QC output.
        """
        return self._qc_data[name]

    def _collect_qc_outputs(self):
        """
        Internal: determine which QC outputs are present
        """
        # QC directory object
        qcdir = QCDir(self.qc_dir,fastq_attrs=self.fastq_attrs)
        # Get list of files
        print("Scanning contents of %s" % self.qc_dir)
        files = qcdir.file_list
        print("\t- %d objects found" % len(files))
        logger.debug("files: %s" % files)
        # Collect QC outputs from modules
        for m in QC_MODULES:
            self._add_qc_outputs(m.collect_qc_outputs(qcdir))
        # Non-QC module outputs
        for qc_data in (
                self._collect_icell8(self.qc_dir),
                self._collect_extra_outputs(self.qc_dir),):
            self._add_qc_outputs(qc_data)
        # Fastq screens
        fastq_screen = self.data('fastq_screen')
        self.fastq_screens = sorted(fastq_screen.screen_names)
        # Sequence lengths
        seq_lengths = self.data('sequence_lengths')
        if "sequence_lengths" in seq_lengths.tags:
            # Maximum number of sequences
            self.stats['max_seqs'] = seq_lengths.max_seqs
            # Maximum and minimum sequence lengths for each read
            for read in seq_lengths.reads:
                self.stats.min_sequence_length_read[read] = \
                    seq_lengths.min_seq_length[read]
                self.stats.max_sequence_length_read[read] = \
                    seq_lengths.max_seq_length[read]
            # Maximum and minimum sequence lengths for all reads
            min_lens = [seq_lengths.min_seq_length[r]
                        for r in seq_lengths.min_seq_length]
            if min_lens:
                self.stats['min_sequence_length'] = min(min_lens)
            max_lens = [seq_lengths.max_seq_length[r]
                        for r in seq_lengths.max_seq_length]
            if max_lens:
                self.stats['max_sequence_length'] = max(max_lens)
        # Single library analyses reference data
        cellranger_references = []
        for name in ('cellranger_count',
                     'cellranger-atac_count',
                     'cellranger-arc_count'):
            cellranger_count = self.data(name)
            cellranger_references.extend([ref for ref in
                                          cellranger_count.references])
        self.cellranger_references = sorted(list(set(cellranger_references)))
        # Cellranger multi
        cellranger_multi = self.data('cellranger_multi')
        # Multiplexed samples
        self.multiplexed_samples = sorted(cellranger_multi.multiplexed_samples,
                                          key=lambda s: split_sample_name(s))
        for reference_data in cellranger_multi.references:
            if reference_data not in self.cellranger_references:
                self.cellranger_references.append(reference_data)
        self.cellranger_references = sorted(self.cellranger_references)
        self.cellranger_probe_sets = cellranger_multi.probe_sets
        # Fastqs
        self.fastqs = sorted(list(self.fastqs))
        # BAMs
        self.bams = sorted(list(self.bams))
        # Configs
        self.config_files = sorted(list(self.config_files))
        # Samples
        samples = set([self.fastq_attrs(fq).sample_name
                       for fq in self.fastqs])
        for bam in self.bams:
            samples.add(self.fastq_attrs(bam).sample_name)
        for s in cellranger_count.samples:
            samples.add(s)
        self.samples = sorted(list(samples),
                              key=lambda s: split_sample_name(s))
        # Samples with biological sequence data
        if self.qc_info is not None and self.qc_info.seq_data_samples:
            # Explicitly defined in QC metadata
            self.seq_data_samples.extend(
                [s for s in self.qc_info.seq_data_samples.split(',')
                 if s in samples])
        for cf in self.config_files:
            # Implicitly defined in 10x multi config file(s) as GEX data
            if cf.startswith('10x_multi_config.'):
                cf = CellrangerMultiConfigCsv(os.path.join(self.qc_dir,cf))
                self.seq_data_samples.extend([s for s in cf.gex_libraries
                                              if s in samples])
        if not self.seq_data_samples:
            # Nothing defined - assume all samples are biological data
            self.seq_data_samples = [s for s in self.samples]
        self.seq_data_samples = sorted(list(set(self.seq_data_samples)))
        # Determine reads
        reads = set()
        for fastq in self.fastqs:
            fq = self.fastq_attrs(fastq)
            if fq.is_index_read:
                reads.add("i%s" % (fq.read_number
                                   if fq.read_number is not None else '1'))
            else:
                reads.add("r%s" % (fq.read_number
                                   if fq.read_number is not None else '1'))
        self.reads = sorted(list(reads))
        # Organisms
        self.organisms = sorted(list(self.organisms))
        # Sort outputs
        self.outputs = sorted(self.outputs)
        # Sort output files
        self.output_files = sorted(self.output_files)

    def _add_qc_outputs(self,data):
        """
        Internal: store 'raw' data for a specific QC output

        The data will be an AttributeDictionary instance
        returned from one of the ``_collect_...`` methods,
        which should contain at minimum the following
        attributes:

        - name: identifier for the QC output (e.g. 'fastqc')
        - fastqs: list of Fastq names associated with the
          output
        - software: a dictionary compromising the names of
          associated software packages as keys, and a list
          of versions as the values
        - output_files: a list of associated output files
        - tags: a list of output 'classes' associated with
          the QC output (e.g. 'fastqc_r1')

        The data can be accessed using the name identifier,
        via the 'data' method; the information about Fastqs,
        software, output files and tags are added to the
        global properties.
        """
        print("- adding data for '%s'" % data.name)
        if data.name in self._qc_data:
            raise KeyError("'%s': data already stored against "
                           "this name" % data.name)
        ##print(data)
        # Store raw data
        self._qc_data[data.name] = data
        # Output files
        for f in data.output_files:
            self.output_files.append(f)
        # Fastqs
        try:
            for fq in data.fastqs:
                self.fastqs.add(fq)
        except AttributeError:
            pass
        # BAMs
        try:
            for bam in data.bam_files:
                self.bams.add(bam)
        except AttributeError:
            pass
        # Organisms
        try:
            for organism in data.organisms:
                self.organisms.add(organism)
        except AttributeError:
            pass
        # Config files
        try:
            for cf in data.config_files:
                self.config_files.add(cf)
        except AttributeError:
            pass
        # Software versions
        for sw in data.software.keys():
            self.software[sw] = [str(v) if v else '?'
                                 for v in data.software[sw]]
        # Tags
        for tag in data.tags:
            self.outputs.add(tag)

    def _collect_icell8(self,qc_dir):
        """
        Collect information on ICell8 reports

        Returns an AttributeDictionary with the following
        attributes:

        - name: set to 'icell8'
        - software: dictionary of software and versions
        - fastqs: list of associated Fastq names
        - output_files: list of associated output files
        - tags: list of associated output classes

        Arguments:
          qc_dir (str): top-level directory to look under.
        """
        output_files = list()
        tags = set()
        # Look for ICELL8 outputs
        icell8_top_dir = os.path.dirname(qc_dir)
        print("Checking for ICELL8 reports in %s/stats" %
              icell8_top_dir)
        icell8_stats_xlsx = os.path.join(icell8_top_dir,
                                         "stats",
                                         "icell8_stats.xlsx")
        if os.path.exists(icell8_stats_xlsx):
            tags.add("icell8_stats")
            output_files.append(icell8_stats_xlsx)
        icell8_report_html = os.path.join(icell8_top_dir,
                                          "icell8_processing.html")
        if os.path.exists(icell8_report_html):
            tags.add("icell8_report")
            output_files.append(icell8_report_html)
        # Return collected information
        return AttributeDictionary(
            name='icell8',
            software={},
            fastqs=[],
            output_files=output_files,
            tags=sorted(list(tags))
        )

    def _collect_extra_outputs(self,qc_dir):
        """
        Collect information on additional outputs

        Returns an AttributeDictionary with the following
        attributes:

        - name: set to 'extra_outputs'
        - software: dictionary of software and versions
        - output_files: list of associated output files
        - tags: list of associated output classes

        Arguments:
          qc_dir (str): top-level directory to look under.
        """
        version = None
        software = dict()
        output_files = set()
        tags = set()
        # Look for extra_outputs.tsv
        extra_outputs_tsv = os.path.join(qc_dir,"extra_outputs.tsv")
        for output in ExtraOutputs(extra_outputs_tsv).outputs:
            output_files.add(os.path.join(qc_dir,output.file_path))
            if output.additional_files:
                # Add in any additional files or dirs
                for f in output.additional_files:
                    output_files.add(os.path.join(qc_dir,f))
        if output_files:
            tags.add("extra_outputs")
        # Return collected information
        return AttributeDictionary(
            name='extra_outputs',
            software=software,
            fastqs=[],
            output_files=sorted(list(output_files)),
            tags=sorted(list(tags))
        )

class ExtraOutputs:
    """
    Class for handling files specifying external QC outputs

    Reads data from the supplied tab-delimited (TSV) file
    specifying one or more external QC output files.

    Each line in the file should have up to three items
    separated by tabs:

    - file or directory (relative to the qc dir)
    - text description (used in HTML)
    - optionally, comma-separated list of additional files
      or directories to include in the final ZIP archive
      (relative to the qc dir)

    Blank lines and lines starting with the '#' comment
    character are ignored.

    The data from each line of the file is then available
    via the 'outputs' attribute, which provides a list
    of 'AttributeDictionary' objects with the following
    properties:

    - 'file_path': relative path to the output file
    - 'description': associated description
    - 'additional_files': list of the associated files

    Arguments:
      tsv_file (str): path to the input TSV file
    """
    def __init__(self,tsv_file):
        self.tsv_file = os.path.abspath(tsv_file)
        self.outputs = list()
        self.__load()

    def __load(self):
        """
        Internal: loads data from the supplied TSV file
        """
        if not os.path.exists(self.tsv_file):
            return
        with open(self.tsv_file,'rt') as fp:
            for line in fp:
                line = line.strip()
                if line.startswith('#') or not line:
                    continue
                items = line.split('\t')
                if len(items) > 3:
                    raise IndexError("Bad line (too many items): '%s'" %
                                     line)
                try:
                    file_path = items[0]
                    description = items[1]
                except IndexError:
                    raise IndexError("Bad line (not enough items): '%s'" %
                                     line)
                try:
                    additional_files = items[2].split(',')
                except IndexError:
                    additional_files = None
                self.outputs.append(
                    AttributeDictionary(file_path=file_path,
                                        description=description,
                                        additional_files=additional_files))
