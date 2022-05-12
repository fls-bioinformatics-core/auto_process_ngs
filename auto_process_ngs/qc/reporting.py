#!/usr/bin/env python
#
#     reporting: report QC from analysis projects
#     Copyright (C) University of Manchester 2018-2022 Peter Briggs
#

"""
Utilities for reporting QC pipeline outputs.

Provides the following core class:

- QCReport: create QC report document for one or more projects

In addition there are a number of supporting classes:

- QCProject: gather information about the QC associated with a project
- SampleQCReporter: reports the QC for a sample
- FastqGroupQCReporter: reports the QC for a group of Fastqs
- FastqQCReporter: interface to QC outputs for a single Fastq

There are also a number of utility functions:

- report: report the QC for a project
- pretty_print_reads: print number of reads with commas at each thousand
- sanitize_name: replace 'unsafe' characters in HTML link targets

Overview
--------

The ``SampleQCReporter``, ``FastqGroupQCReporter`` and
``FastqQCReporter`` classes are used by the top-level ``QCProject``
class to report QC outputs at the level of samples (for example,
single library analyses), groups of Fastqs (for example, strandedness),
and individual Fastqs (for example, FastQC or screen data).

Adding support for new metadata
-------------------------------

Support for new metadata items should be implemented within the
``_init_metadata_table`` method of ``QCReport``. Descriptions of
new items should also be added to the ``METADATA_FIELD_DESCRIPTIONS``
module constant.

Adding support for new QC outputs
---------------------------------

When adding reporting of new QC outputs it is recommended first to
ensure that they are detected by the ``QCOutputs`` class (in the
``outputs`` module); then the relevant reporter class should be
extended depending on the level that the QC outputs are associated
with (i.e. project, sample, Fastq group or individual Fastq).

Typically this is done by adding support for new summary table fields,
which can implemented within the ``get_value`` or ``get_10x_value``
methods in ``SampleQCreporter`` (for sample-level QC) or the ``get_value``
method of ``FastqGroupQCReporter`` (for Fastq-group level QC). (This
may require additional internal functionality to be implemented
within the relevant class.)

Descriptions for new fields also need to be added to the
``SUMMARY_FIELD_DESCRIPTIONS`` module constant. Additionally: if the
QC outputs are produced using new software packages then these should
be added to the ``SOFTWARE_PACKAGE_NAMES`` module constant; as long
as these are reported by ``QCOutputs`` then they will also be listed
automatically within the QC report.
"""

#######################################################################
# Imports
#######################################################################

import sys
import os
import re
import logging
import time
import ast
import shutil
import textwrap
from collections import defaultdict
from bcftbx.IlluminaData import IlluminaFastq
from bcftbx.TabFile import TabFile
from bcftbx.qc.report import strip_ngs_extensions
from bcftbx.utils import AttributeDictionary
from bcftbx.utils import extract_prefix
from bcftbx.utils import extract_index
from bcftbx.utils import walk
from ..analysis import AnalysisFastq
from ..analysis import run_reference_id
from ..analysis import split_sample_name
from ..analysis import split_sample_reference
from ..docwriter import Document
from ..docwriter import Section
from ..docwriter import Table
from ..docwriter import Img
from ..docwriter import Link
from ..docwriter import Target
from ..docwriter import List
from ..docwriter import Para
from ..docwriter import WarningIcon
from ..docwriter import DocumentIcon
from ..metadata import AnalysisDirMetadata
from ..metadata import AnalysisProjectQCDirInfo
from ..fastq_utils import group_fastqs_by_name
from .fastqc import Fastqc
from .fastq_screen import Fastqscreen
from .fastq_strand import Fastqstrand
from .cellranger import CellrangerCount
from .cellranger import CellrangerMulti
from .outputs import fastqc_output
from .outputs import fastq_screen_output
from .outputs import fastq_strand_output
from .outputs import QCOutputs
from .picard import CollectInsertSizeMetrics
from .plots import RGB_COLORS
from .plots import Plot
from .plots import useqlenplot
from .plots import ureadcountplot
from .plots import uscreenplot
from .plots import ufastqcplot
from .plots import uboxplot
from .plots import ustrandplot
from .plots import uduplicationplot
from .plots import uadapterplot
from .plots import encode_png
from .qualimap import QualimapRnaseq
from .rseqc import InferExperiment
from .seqlens import SeqLens
from ..tenx_genomics_utils import MultiomeLibraries
from ..utils import ZipArchive
from .. import get_version

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Data
#######################################################################

from .constants import FASTQ_SCREENS
from .constants import QC_REPORT_CSS_STYLES

# Metadata field descriptions
METADATA_FIELD_DESCRIPTIONS = {
    'project_id': 'Project ID',
    'run_id': 'Run ID',
    'run': 'Run name',
    'user': 'User',
    'PI': 'PI',
    'library_type': 'Library type',
    'sequencer_model': 'Sequencer model',
    'single_cell_platform': 'Single cell preparation platform',
    'number_of_cells': 'Number of cells',
    'organism': 'Organism',
    'protocol': 'QC protocol',
    'cellranger_reference': 'Cellranger reference datasets',
    'multiqc': 'MultiQC report',
    'icell8_stats': 'ICELL8 statistics',
    'icell8_report': 'ICELL8 processing report',
}

# Software package names
SOFTWARE_PACKAGE_NAMES = {
    'bcl2fastq': 'Bcl2fastq',
    'bcl-convert': 'BCL Convert',
    'cellranger': 'Cellranger',
    'cellranger-atac': 'Cellranger ATAC',
    'cellranger-arc': 'Cellranger ARC',
    'fastqc': 'FastQC',
    'fastq_screen': 'FastqScreen',
    'fastq_strand': 'FastqStrand',
    'multiqc': 'MultiQC',
    'picard': 'Picard Tools',
    'qualimap': 'Qualimap',
    'rseqc:genebody_coverage': 'RSeQC geneBody_coverage.py',
    'spaceranger': 'Spaceranger',
}

# Field descriptions for summary tables
# Dictionary with keys matching field labels, and values
# consisting of tuples with short and long descriptions
SUMMARY_FIELD_DESCRIPTIONS = {
    'sample': ('Sample','Sample name'),
    'fastq' : ('Fastq','Fastq file'),
    'fastqs': ('Fastqs','Fastq files in each sample'),
    'bam_file': ('BAM','BAM file'),
    'reads': ('#reads','Number of reads/read pairs'),
    'read_lengths': ('Lengths',
                     'Mean sequence length and range'),
    'read_counts': ('Counts',
                    'Relative number of total reads, and proportions '
                    'of masked and padded reads in each Fastq'),
    'sequence_duplication': ('Dup%',
                             'Fraction of reads with duplicated '
                             'sequences in each Fastq'),
    'adapter_content': ('Adapters',
                        'Fraction of data containing adapter sequences '
                        'in each Fastq'),
    'read_lengths_dist_r1': ('Dist[R1]',
                             'Distributions of R1 sequence lengths'),
    'fastqc_r1': ('FastQC[R1]',
                  'Summary of FastQC metrics for R1'),
    'boxplot_r1': ('Quality[R1]',
                   'Per base sequence quality for R1'),
    'screens_r1': ('Screens[R1]',
                   'Outputs from FastqScreen running R1 against multiple '
                   'panels'),
    'read_lengths_dist_r2': ('Dist[R2]',
                             'Distributions of R2 sequence lengths'),
    'fastqc_r2': ('FastQC[R2]',
                  'Summary of FastQC metrics for R2'),
    'boxplot_r2': ('Quality[R2]',
                   'Per base sequence quality for R2'),
    'screens_r2': ('Screens[R2]',
                   'Outputs from FastqScreen running R2 against multiple '
                   'panels'),
    'fastqc_r3': ('FastQC[R3]',
                  'Summary of FastQC metrics for R3'),
    'read_lengths_dist_r3': ('Dist[R3]',
                             'Distributions of R3 sequence lengths'),
    'boxplot_r3': ('Quality[R3]',
                   'Per base sequence quality for R3'),
    'screens_r3': ('Screens[R3]',
                   'Outputs from FastqScreen running R3 against multiple '
                   'panels'),
    'strandedness': ('Strand',
                     'Proportions of reads mapping to forward and reverse '
                     'strands'),
    'strandedness_.*': ('Strand',
                        'Fraction of reads mapping to forward and reverse '
                        'strands (from RSeQC infer_experiment.py)'),
    'endedness_.*' : ('Endedness',
                      'Whether data is Paired- or single-ended (from RSeQC '
                      'infer_experiment.py'),
    'insert_size_metrics_.*': ('Insert size (bp)',
                               'Picard insert size: mean (SD)'),
    'insert_size_histogram_.*': ('Insert size histogram',
                                 'Picard insert size histogram'),
    'coverage_profile_along_genes_.*': ('Gene coverage',
                                        'Qualimap rna-seq gene coverage '
                                        'profile'),
    'reads_genomic_origin_.*': ('Genomic origin',
                                'Qualimap rna-seq genomic origin of '
                                'reads'),
    'qualimap_rnaseq_report_.*': ('Qualimap report',
                                  'Qualimap rna-seq HTML report'),
    'cellranger_count': ('Single library analyses',
                         'Web summary from Cellranger* single library '
                         'analysis for this sample'),
    '10x_cells': ('#cells','Number of cells'),
    '10x_reads_per_cell': ('#reads/cell','Average reads per cell'),
    '10x_genes_per_cell': ('#genes/cell','Median genes per cell'),
    '10x_frac_reads_in_cell': ('%reads in cells',
                               'Fraction of reads in cells'),
    '10x_fragments_per_cell': ('#fragments/cell',
                               'Median fragments per cell'),
    '10x_fragments_overlapping_targets': ('%fragments overlapping targets',
                                          'Fraction of fragments '
                                          'overlapping targets'),
    '10x_fragments_overlapping_peaks': ('%fragments overlapping peaks',
                                        'Fraction of fragments overlapping '
                                        'peaks'),
    '10x_tss_enrichment_score': ('TSS enrichment score',
                                 'TSS enrichment score'),
    '10x_atac_fragments_per_cell': ('#ATAC fragments/cell',
                                    'ATAC Median high-quality fragments '
                                    'per cell'),
    '10x_gex_genes_per_cell': ('#GEX genes/cell','GEX Median genes per cell'),
    '10x_genes_detected': ('#genes','Total genes detected'),
    '10x_umis_per_cell': ('#UMIs/cell','Median UMI counts per cell'),
    '10x_pipeline': ('Pipeline','Name of the 10x Genomics pipeline used'),
    '10x_reference': ('Reference dataset',
                      'Reference dataset used for the analysis'),
    '10x_web_summary': ('HTML report','Link to the web_summary.html report'),
    'linked_sample': ('Linked sample',
                      'Corresponding sample for single cell multiome analysis')
}

#######################################################################
# Classes
#######################################################################

class QCReport(Document):
    """
    Create a QC report document for one or more projects

    Example usage:

    >>> report = QCReport(project)
    >>> report.write("qc_report.html")

    To control the fields written to the summary table, specify
    a list of field names via the 'summary_fields' argument.
    Valid field names are:

    - sample: sample name
    - fastq: Fastq name
    - fastqs: Fastq R1/R2 names
    - reads: number of reads
    - read_lengths: length of reads
    - read_lengths_distributions: mini-plots of read length distributions
    - read_counts: mini-plots of fractions of masked/padded/total reads
    - adapter_content: mini-plots summarising adapter content for all reads
    - read_length_dist_[read]: length dist mini-plot for [read] (r1,r2,...)
    - fastqc_[read]: FastQC mini-plot for [read]
    - boxplot_[read]: FastQC per-base-quality mini-boxplot' for [read]
    - screens_[read]: FastQScreen mini-plots for [read]
    - strandedness: 'forward', 'reverse' or 'unstranded' for pair
    - cellranger_count: 'cellranger count' outputs for each sample

    To control the elements written to the reports for each Fastq
    pair, specify a list of element names via the 'report_attrs'
    argument. Valid element names are:

    - fastqc: FastQC report
    - fastq_screen: FastQCScreen report
    - program_versions: program versions

    Arguments:
      projects (AnalysisProject): list of projects to report
        QC for
      title (str): title for the report (defaults to
        "QC report: <PROJECT_NAME>")
      qc_dir (str): path to the QC output dir; relative path
        will be treated as a subdirectory of the project
      report_attrs (list): list of elements to report for
        each Fastq pair
      summary_fields (list): list of fields to report for each
        sample in the summary table
      relpath (str): if set then make link paths relative to
        'relpath'
      data_dir (str): if set then copy external data files to
        this directory and make link paths to these copies;
        relative path will be treated as a subdirectory of the
        project
    """
    # Field descriptions for summary table
    field_descriptions = SUMMARY_FIELD_DESCRIPTIONS
    # Titles for metadata items
    metadata_titles = METADATA_FIELD_DESCRIPTIONS
    # Software packages and names
    software_names = SOFTWARE_PACKAGE_NAMES
    software_packages = sorted(list(software_names.keys()))

    def __init__(self,projects,title=None,qc_dir=None,report_attrs=None,
                 summary_fields=None,relpath=None,data_dir=None):
        """
        Create a new QCReport instance
        """
        # Convert projects to QCProjects
        projects = [QCProject(p,qc_dir=qc_dir) for p in projects]
        # Flag to indicate if report references multiple projects
        self.multi_project = (len(projects) > 1)
        # Flag to indicate if any project has single cell data
        self.has_single_cell = any([p.is_single_cell for p in projects])
        # Flag to indicate if we should make a separate single library
        # analysis table
        self.use_single_library_table = True
        # Primary project
        primary_project = projects[0]
        # Set up title
        if title is None:
            if self.multi_project:
                title = "QC report: %s &amp; %s" % (
                    ', '.join([p.id for p in projects[0:-1]]),
                    projects[-1].id)
            else:
                title = "QC report: %s" % primary_project.name
        # Initialise superclass
        Document.__init__(self,title)
        # Collect initial QC data across all projects
        self.output_files = []
        self.outputs = set()
        for project in projects:
            # Output files
            if project.output_files:
                self.output_files.extend([f for f in project.output_files])
            # Outputs
            for output in project.outputs:
                self.outputs.add(output)
        self.outputs = sorted(list(self.outputs))
        # Software packages
        self.software = []
        for pkg in self.software_packages:
            for project in projects:
                if project.software_info(pkg,exclude_processing=True):
                    self.software.append(pkg)
                    break
        # Use data directory?
        if data_dir is not None:
            if not os.path.abspath(data_dir):
                data_dir = os.path.join(primary_project.dirn,data_dir)
            data_dir = os.path.normpath(data_dir)
        self.data_dir = data_dir
        # Relative paths
        if relpath is not None:
            relpath = os.path.normpath(os.path.abspath(relpath))
        self.relpath = relpath
        # Status of report
        self.status = True
        # Initialise tables
        self._init_metadata_table(projects)
        self._init_processing_software_table()
        self._init_qc_software_table()
        # Initialise report sections
        self.preamble = self._init_preamble_section()
        self.warnings = self._init_warnings_section()
        self.summary = self._init_summary_section(project)
        # Initialise data directory
        self._init_data_dir(projects)
        # Report each project
        for project in projects:
            # Report QC data for project
            print("Project: %s" % project.name)
            if project.samples:
                print("Samples found:")
                for sample in project.samples:
                    print("\t- %s" % sample)
            else:
                logger.warning("%s: no samples found" % project.name)
            if project.multiplexed_samples:
                print("Multiplexed samples found:")
                for sample in project.multiplexed_samples:
                    print("\t- %s" % sample)
            if project.fastqs:
                print("Fastqs referenced:")
                for fastq in project.fastqs:
                    print("\t- %s" % fastq)
                print("Reads found:")
                for read in project.reads:
                    print("\t- %s" % read)
            else:
                logger.warning("%s: no Fastqs referenced" % project.name)
            if project.outputs:
                print("Available QC outputs:")
                for output in project.outputs:
                    print("\t- %s" % output)
            else:
                logger.warning("%s: no QC outputs found" % project.name)
            if project.fastq_screens:
                print("Fastq screens:")
                for screen in project.fastq_screens:
                    print("\t- %s" % screen)
            if project.stats.max_seqs:
                print("Maximum number of sequences: %d" %
                      project.stats.max_seqs)
            if project.stats.min_sequence_length:
                print("Minimum sequence length    : %d" %
                      project.stats.min_sequence_length)
            for read in project.reads:
                if read in project.stats.min_sequence_length_read:
                    print("\t- %s: %s" %
                          (read,project.stats.min_sequence_length_read[read]))
            if project.stats.max_sequence_length:
                print("Maximum sequence length    : %d" %
                      project.stats.max_sequence_length)
            for read in project.reads:
                if read in project.stats.max_sequence_length_read:
                    print("\t- %s: %s" %
                          (read,project.stats.max_sequence_length_read[read]))
            if project.software:
                print("Software versions:")
                for package in project.software:
                    print("\t- %s: %s" %
                          (package,
                           ', '.join(project.software[package])))
            # Fields to report in summary table
            if not summary_fields:
                if len(project.reads) > 1:
                    summary_fields_ = ['sample',
                                       'fastqs',
                                       'reads',
                                       'read_counts',
                                       'read_lengths',
                                       'strandedness',
                                       'sequence_duplication',
                                       'adapter_content']
                else:
                    summary_fields_ = ['sample',
                                       'fastq',
                                       'reads',
                                       'read_counts',
                                       'read_lengths',
                                       'strandedness',
                                       'sequence_duplication',
                                       'adapter_content',]
                if 'strandedness' not in project.outputs:
                    try:
                        summary_fields_.remove('strandedness')
                    except ValueError:
                        pass
                if 'sequence_lengths' in project.outputs:
                    for read in project.reads:
                        summary_fields_.append('read_lengths_dist_%s' % read)
                else:
                    try:
                        summary_fields_.remove('read_counts')
                    except ValueError:
                        pass
                for read in project.reads:
                    if ('fastqc_%s' % read) in project.outputs:
                        summary_fields_.append('fastqc_%s' % read)
                for read in project.reads:
                    if ('fastqc_%s' % read) in project.outputs:
                        summary_fields_.append('boxplot_%s' % read)
                for read in project.reads:
                    if ('screens_%s' % read) in project.outputs:
                        summary_fields_.append('screens_%s' % read)
                if 'cellranger_count' in project.outputs and \
                   not self.use_single_library_table:
                    summary_fields_.append('cellranger_count')
            # Attributes to report for each sample
            if report_attrs is None:
                report_attrs_ = ['fastqc',
                                 'fastq_screen',]
                if 'strandedness' in project.outputs:
                    report_attrs_.append('strandedness')
            # Add data for this project to the report
            print("Adding project '%s' to the report..." % project.name)
            self.report_metadata(project)
            self.report_processing_software(project)
            self.report_qc_software(project)
            self.report_comments(project)
            # Create a summary subsection for multi-project reporting
            if self.multi_project:
                project_summary = self.summary.add_subsection(
                    project.id,
                    name=sanitize_name(project.id))
            else:
                project_summary = self.summary
            # Create a new summary table
            project_summary.add(
                "%d sample%s | %d fastq%s" % (
                    len(project.samples),
                    ('s' if len(project.samples) != 1 else ''),
                    len(project.fastqs),
                    ('s' if len(project.fastqs) != 1 else ''))
            )
            summary_table = self.add_summary_table(project,
                                                   summary_fields_,
                                                   section=project_summary)
            # Report each sample
            for sample in project.samples:
                self.report_sample(project,sample,report_attrs_,
                                   summary_table,summary_fields_)
            # Report single library analyses
            if self.use_single_library_table:
                for single_library in ('cellranger_count',
                                       'cellranger-atac_count',
                                       'cellranger-arc_count'):
                    if single_library not in project.outputs:
                        # Skip missing analysis
                        continue
                    # Set up fields for reporting
                    if single_library == 'cellranger_count':
                        pkg = 'cellranger'
                        single_library_fields = ['sample',
                                                 '10x_cells',
                                                 '10x_frac_reads_in_cell',
                                                 '10x_reads_per_cell',
                                                 '10x_genes_per_cell']
                    elif single_library == 'cellranger-atac_count':
                        pkg = 'cellranger-atac'
                        single_library_fields = ['sample',
                                                 '10x_cells',
                                                 '10x_fragments_per_cell',
                                                 '10x_tss_enrichment_score']
                        for v in project.software[pkg]:
                            # Add version specific fields to summary table
                            v = v.split('.')
                            if v[0] == '2':
                                extra_fields = \
                                    ['10x_fragments_overlapping_peaks']
                            else:
                                extra_fields = \
                                    ['10x_fragments_overlapping_targets']
                            for f in extra_fields:
                                if f not in single_library_fields:
                                    single_library_fields.append(f)
                    elif single_library == 'cellranger-arc_count':
                        pkg = 'cellranger-arc'
                        single_library_fields = ['sample',
                                                 'linked_sample',
                                                 '10x_cells',
                                                 '10x_atac_fragments_per_cell',
                                                 '10x_gex_genes_per_cell']
                    # Add column for multiple versions
                    if len(project.software[pkg]) > 1:
                        single_library_fields.append('10x_pipeline')
                    # Add column for multiple reference datasets
                    if len(project.cellranger_references) > 1:
                        single_library_fields.append('10x_reference')
                    # Always link to web summary
                    single_library_fields.append('10x_web_summary')
                    # Create a new table
                    single_library_analysis_table = \
                        self.add_single_library_analysis_table(
                            pkg,
                            project,
                            single_library_fields,
                            section=project_summary)
                    # Report analyses for each sample
                    for sample in project.samples:
                        self.report_single_library_analyses(
                            pkg,
                            project,
                            sample,
                            single_library_analysis_table,
                            single_library_fields)
            # Report 10x multiplexing analyses
            if 'cellranger_multi' in project.outputs:
                # Set up fields for reporting
                pkg = 'cellranger'
                multiplex_analysis_fields = ['sample',
                                             '10x_cells',
                                             '10x_reads_per_cell',
                                             '10x_genes_per_cell',
                                             '10x_genes_detected',
                                             '10x_umis_per_cell']
                # Add column for multiple versions
                if len(project.software[pkg]) > 1:
                    multiplex_analysis_fields.append('10x_pipeline')
                # Add column for multiple reference datasets
                if len(project.cellranger_references) > 1:
                    multiplex_analysis_fields.append('10x_reference')
                # Always link to web summary
                multiplex_analysis_fields.append('10x_web_summary')
                # Create a new table
                multiplex_analysis_table = \
                    self.add_multiplex_analysis_table(
                        project,
                        multiplex_analysis_fields,
                        section=project_summary)
                # Report analyses for each multiplexed sample
                for sample in project.multiplexed_samples:
                    self.report_multiplexing_analyses(
                        project,
                        sample,
                        multiplex_analysis_table,
                        multiplex_analysis_fields)
            # Report extended QC metrics
            include_extended_metrics = \
                ('rseqc_genebody_coverage' in project.outputs) or \
                ('rseqc_infer_experiment' in project.outputs) or \
                ('picard_insert_size_metrics' in project.outputs) or \
                ('qualimap_rnaseq' in project.outputs)
            if include_extended_metrics:
                self.report_extended_qc_metrics(
                    project,
                    section=project_summary)
        # Report the status
        self.report_status()

    def _init_metadata_table(self,projects):
        """
        Internal: set up a table for project metadata

        Associated CSS class is 'metadata'
        """
        # Identify metadata items
        metadata_items = ['run_id',
                          'run',
                          'sequencer_model',
                          'user',
                          'PI',
                          'library_type',
                          'organism',
                          'protocol',]
        if self.has_single_cell:
            for item in ('single_cell_platform',
                         'number_of_cells',):
                metadata_items.insert(metadata_items.index('organism'),
                                      item)
        if 'cellranger_count' in self.outputs:
            metadata_items.append('cellranger_reference')
        if 'multiqc' in self.outputs:
            metadata_items.append('multiqc')
        if 'icell8_stats' in self.outputs:
            metadata_items.append('icell8_stats')
        if 'icell8_report' in self.outputs:
            metadata_items.append('icell8_report')
        if self.multi_project:
            metadata_items[metadata_items.index('run_id')] = 'project_id'
        # Make table with one column per project
        columns = ['item']
        for project in projects:
            columns.append(project.id)
        metadata_table = Table(columns)
        metadata_table.no_header()
        metadata_table.add_css_classes('metadata')
        # Add rows for metadata items
        for item in metadata_items:
            metadata_table.add_row(item=self.metadata_titles[item])
        # Store table and metadata items as attribute
        self.metadata_table = metadata_table
        self.metadata_items = metadata_items

    def _init_processing_software_table(self):
        """
        Internal: set up a table for processing software information

        Associated CSS class is 'metadata'
        """
        # Determine what software packages were used across
        # all projects
        software_table = Table(('program','version',))
        software_table.no_header()
        software_table.add_css_classes('metadata')
        self.processing_software_table = software_table

    def _init_qc_software_table(self):
        """
        Internal: set up a table for QC software information

        Associated CSS class is 'metadata'
        """
        # Determine what software packages were used across
        # all projects
        software_table = Table(('program','version',))
        software_table.no_header()
        software_table.add_css_classes('metadata')
        for pkg in self.software:
            software_table.add_row(
                program=self.software_names[pkg],
                version=None)
        self.qc_software_table = software_table

    def _init_preamble_section(self):
        """
        Internal: set up a "preamble" section

        Associated name is 'preamble'
        """
        preamble = self.add_section(name='preamble')
        preamble.add("Report generated by auto_process %s on %s" %
                     (get_version(),time.asctime()))
        return preamble

    def _init_summary_section(self,project):
        """
        Internal: set up a summary section for the report

        Associated name is 'summary'
        """
        summary = self.add_section("Summary",name='summary',
                                   css_classes=("summary",))
        # Add subsections inside a container
        info = summary.add_subsection()
        general_info = info.add_subsection("General information",
                                           css_classes=("info",))
        general_info.add(self.metadata_table)
        processing_software_info = info.add_subsection(
            "Processing software",
            css_classes=("info",))
        processing_software_info.add(self.processing_software_table)
        qc_software_info = info.add_subsection("QC software",
                                               css_classes=("info",))
        qc_software_info.add(self.qc_software_table)
        # Add an empty section to clear HTML floats
        clear = summary.add_subsection(css_classes=("clear",))
        # Add additional subsections for comments etc
        info = summary.add_subsection()
        # Add comments section
        self.comments = info.add_subsection("Comments",
                                            css_classes=("info",))
        # Add an empty section to clear HTML floats
        clear = summary.add_subsection(css_classes=("clear",))
        return summary

    def _init_warnings_section(self):
        """
        Internal: creates a section which is displayed if there
        are warnings
        """
        warnings = self.add_section(name='warnings',
                                    css_classes=("warnings",))
        warnings.add(Para(WarningIcon(size=50),
                          "There are issues with this project"))
        return warnings

    def _init_data_dir(self,projects):
        """
        Internal: initialise the data directory

        Creates directory, copies QC artefacts and updates the
        list of output files associated with the report
        """
        if not self.data_dir:
            return
        # Remove pre-existing data directory
        if os.path.exists(self.data_dir):
            print("Removing existing data directory: %s" %
                  self.data_dir)
            shutil.rmtree(self.data_dir)
        # (Re)create directory
        print("Creating data directory: %s" % self.data_dir)
        os.mkdir(self.data_dir)
        # Copy QC artefacts from each project
        output_files = []
        for project in projects:
            # Set up a project-specific top-level directory
            project_data_dir = os.path.join(self.data_dir,
                                            sanitize_name(project.id))
            if not os.path.exists(project_data_dir):
                logging.debug("Creating subdirectory %s" %
                              project_data_dir)
                os.makedirs(project_data_dir)
            for output_file in project.output_files:
                # Convert file paths to be relative to project dir
                output_file = os.path.relpath(output_file,
                                              project.dirn)
                # Make a copy of the file in the data dir
                subdir = os.path.dirname(output_file)
                # File includes a leading directory path so
                # build a copy of this first
                if subdir:
                    subdir = os.path.join(project_data_dir,subdir)
                    if not os.path.exists(subdir):
                        os.makedirs(subdir)
                # Copy the file
                src = os.path.join(project.dirn,output_file)
                dest = os.path.join(project_data_dir,output_file)
                logging.debug("-- copying %s to %s" % (output_file,
                                                       dest))
                shutil.copyfile(src,dest)
                # Add to the updated list of output files
                output_files.append(dest)
        # Update the list of files
        self.output_files = output_files

    def add_summary_table(self,project,fields,section):
        """
        Create a new table for summarising samples from a project

        Associated CSS classes are 'summary' and 'fastq_summary'
        """
        # Generate column headers for table
        tbl_headers = {}
        for field in fields:
            # Look for matching field
            th = None
            for f in self.field_descriptions.keys():
                if re.match('^%s$' % f,field):
                    # Construct the header
                    th = "<span title=\"%s\">%s</span>" % \
                                       ('\n'.join(textwrap.wrap(
                                           self.field_descriptions[f][1],
                                           width=20)),
                                       self.field_descriptions[f][0])
                    break
            if th is None:
                # Generic header for unrecognised field
                th = "<span title=\"%s\">%s</span>" % (field,field)
            tbl_headers[field] = th
        # Create the table
        summary_tbl = Table(fields,**tbl_headers)
        summary_tbl.add_css_classes('summary','fastq_summary')
        if "cellranger_count" in fields:
            summary_tbl.add_css_classes('single_library_analyses',
                                        column='cellranger_count')
        section.add(summary_tbl)
        return summary_tbl

    def add_single_library_analysis_table(self,package,project,fields,
                                          section):
        """
        Create a new table for summarising 10x single library analyses
        """
        # Add title
        section = section.add_subsection(
            "Single library analysis (%s)" % package,
            name="single_library_analysis_%s_%s" % (package,
                                                    sanitize_name(project.id)),
            css_classes=('single_library_summary',))
        # Generate headers for table
        tbl_headers = {
            f: "<space title=\"%s\">%s</span>" %
            ('\n'.join(textwrap.wrap(self.field_descriptions[f][1],width=20)),
             self.field_descriptions[f][0])
            for f in self.field_descriptions.keys()
        }
        # Create the table
        single_library_tbl = Table(fields,**tbl_headers)
        single_library_tbl.add_css_classes('summary','single_library_summary')
        # Append to the summary section
        section.add(single_library_tbl)
        return single_library_tbl

    def add_multiplex_analysis_table(self,project,fields,section):
        """
        Create a new table for summarising 10x multiplexing analyses
        """
        # Add title
        section = section.add_subsection(
            "Multiplexing analysis",
            name="multiplexing_analysis_%s" % sanitize_name(project.id),
            css_classes=('single_library_summary',))
        # Generate headers for table
        tbl_headers = {
            f: "<space title=\"%s\">%s</span>" %
            ('\n'.join(textwrap.wrap(self.field_descriptions[f][1],width=20)),
             self.field_descriptions[f][0])
            for f in self.field_descriptions.keys()
        }
        # Create the table
        multiplexing_tbl = Table(fields,**tbl_headers)
        multiplexing_tbl.add_css_classes('summary','single_library_summary')
        # Append to the summary section
        section.add(multiplexing_tbl)
        return multiplexing_tbl

    def report_metadata(self,project):
        """
        Report the project metadata

        Adds entries for the project metadata to the "metadata"
        table in the report
        """
        # Determine the root directory for QC outputs
        if self.data_dir:
            project_data_dir = os.path.join(self.data_dir,
                                            sanitize_name(project.id))
        else:
            project_data_dir = project.dirn
        if self.relpath:
            project_data_dir = os.path.relpath(project_data_dir,
                                               self.relpath)
        # Add metadata items
        for idx,item in enumerate(self.metadata_items):
            # Try to acquire the value from QC metadata
            try:
                value = project.qc_info[item]
            except KeyError:
                # Fall back to project metadata
                try:
                    if project.info[item]:
                        value = project.info[item]
                    else:
                        # No value set, skip this item
                        continue
                except KeyError:
                    # Additional non-metadata items, or items
                    # requiring additional processing
                    if item == 'project_id':
                        value = Link(project.id,
                                     "#%s" % sanitize_name(project.id))
                    elif item == 'run_id':
                        value = project.run_id
                        if value is None:
                            # Unable to determine run id
                            continue
                    elif item == 'cellranger_reference':
                        if len(project.cellranger_references) == 1:
                            # Single reference dataset
                            value = os.path.basename(
                                project.cellranger_references[0])
                        elif len(project.cellranger_references) > 1:
                            # Many reference datasets
                            value = List()
                            for ref in project.cellranger_references:
                                value.add_item(os.path.basename(ref))
                        else:
                            # No reference datasets
                            continue
                    elif item == 'multiqc':
                        multiqc_report = "multi%s_report.html" \
                                         % os.path.basename(project.qc_dir)
                        value = Link(multiqc_report)
                    elif item == 'icell8_stats':
                        value = Link("icell8_stats.xlsx",
                                     os.path.join(project_data_dir,
                                                  "stats",
                                                  "icell8_stats.xlsx"))
                    elif item == 'icell8_report':
                        value = Link("icell8_processing.html",
                                     os.path.join(project_data_dir,
                                                  "icell8_processing.html"))
                    else:
                        raise Exception("Unrecognised item to report: '%s'"
                                        % item)
            # Update the value in the metadata table
            self.metadata_table.set_value(idx,
                                          project.id,
                                          value)

    def report_processing_software(self,project):
        """
        Report the software versions used in the processing

        Adds entries for the software versions to the
        "processing software" table in the report
        """
        # Report processing software packages
        for pkg in project.processing_software:
            # Acquire the value
            try:
                program = self.software_names[pkg]
            except KeyError:
                # Unknown software package
                program = pkg
            value = project.processing_software[pkg][2]
            self.processing_software_table.add_row(
                program=program,
                version=value)

    def report_qc_software(self,project):
        """
        Report the software versions used in the QC

        Adds entries for the software versions to the
        "qc_software" table in the report
        """
        # Report software packages
        for pkg in self.software_packages:
            # Acquire the value
            value = project.software_info(pkg,
                                          exclude_processing=True)
            if value is None:
                continue
            # Get row index in the metadata table
            idx = self.software.index(pkg)
            self.qc_software_table.set_value(idx,
                                             "version",
                                             value)

    def report_comments(self,project):
        """
        Report the comments associated with a project

        Adds the comments from the project metadata as
        a list to the comments section.

        Arguments:
          project (QCProject): project to report
        """
        comments_list = List()
        try:
            if project.comments:
                for comment in project.comments.split(';'):
                    comments_list.add_item(comment.strip())
            else:
                # Drop out with exception
                raise AttributeError
        except AttributeError:
            comments_list.add_item("N/A")
        if self.multi_project:
            self.comments.add("<p><b>%s</b></p>" % project.id)
        self.comments.add(comments_list)

    def report_sample(self,project,sample,report_attrs,summary_table,
                      summary_fields):
        """
        Report the QC for a sample

        Reports the QC for the sample and Fastqs to
        the summary table and appends a section with
        detailed reports to the document.

        Arguments:
          project (QCProject): project to report
          sample (str): name of sample to report
          report_attrs (list): list of elements to report for
            each set of Fastqs
          summary_table (Table): summary table to report
            each sample in
          summary_fields (list): list of fields to report
            for each sample in the summary table
        """
        # Create a unique name and title
        if self.multi_project:
            sample_name = "sample_%s_%s" % (sanitize_name(project.id),
                                            sample)
            sample_title = "Sample: %s/%s" % (project.id,sample)
        else:
            sample_name = "sample_%s" % sample
            sample_title = "Sample: %s" % sample
        # Determine location of QC artefacts
        qc_dir = self.fetch_qc_dir(project)
        # Create a new section
        sample_report = self.add_section(
            sample_title,
            name=sample_name,
            css_classes=('sample',))
        reporter = SampleQCReporter(project,
                                    sample,
                                    qc_dir=qc_dir,
                                    fastq_attrs=project.fastq_attrs)
        reads = reporter.reads
        n_fastq_groups = len(reporter.fastq_groups)
        if len(reads) == 1:
            sample_report.add("%d %s Fastq%s" %
                              (n_fastq_groups,
                               reads[0].upper(),
                               's' if n_fastq_groups > 1 else ''))
        elif len(reads) == 2:
            sample_report.add("%d %s Fastq pair%s" %
                              (n_fastq_groups,
                               '/'.join([r.upper() for r in reads]),
                              's' if n_fastq_groups > 1 else ''))
        else:
            sample_report.add("%d %s Fastq group%s" %
                              (n_fastq_groups,
                               '/'.join([r.upper() for r in reads]),
                              's' if n_fastq_groups > 1 else ''))
        # Report the Fastq groups within the sample
        reporter.report_fastq_groups(sample_report,
                                     attrs=report_attrs,
                                     relpath=self.relpath)
        # Update the summary table
        status = reporter.update_summary_table(summary_table,
                                               sample_report=sample_report,
                                               fields=summary_fields,
                                               relpath=self.relpath)
        if not status:
            # Update flag to indicate problems with the
            # report
            self.status = False

    def report_single_library_analyses(self,package,project,sample,
                                       single_library_analysis_table,
                                       fields):
        """
        Report the single library analyses for a sample

        Writes lines to the single library analysis summary
        table for each analysis found that is associated
        with the specified sample.

        Arguments:
          package (str): name of 10x package to report
          project (QCProject): project to report
          sample (str): name of sample to report
          single_library_analysis_table (Table): summary table
            to report each analysis in
          fields (list): list of fields to report for each
            analysis in the summary table
        """
        # Determine location of QC artefacts
        qc_dir = self.fetch_qc_dir(project)
        # Get a reporter for the sample
        reporter = SampleQCReporter(project,
                                    sample,
                                    qc_dir=qc_dir,
                                    fastq_attrs=project.fastq_attrs)
        # Update the single library analysis table
        status = reporter.update_single_library_table(
            package,
            single_library_analysis_table,
            fields,
            relpath=self.relpath)
        if not status:
            # Update flag to indicate problems with the
            # report
            self.status = False

    def report_multiplexing_analyses(self,project,sample,
                                     multiplexing_analysis_table,
                                     fields):
        """
        Report the multiplexing analyses for a sample

        Writes lines to the multiplexing analysis summary
        table for each analysis found that is associated
        with the specified sample.

        Arguments:
          project (QCProject): project to report
          sample (str): name of multiplexed sample to report
          multiplexing_analysis_table (Table): summary table
            to report each analysis in
          fields (list): list of fields to report for each
            analysis in the summary table
        """
        # Determine location of QC artefacts
        qc_dir = self.fetch_qc_dir(project)
        # Get a reporter for the sample
        reporter = SampleQCReporter(project,
                                    sample,
                                    qc_dir=qc_dir)
        # Update the single library analysis table
        status = reporter.update_multiplexing_analysis_table(
            multiplexing_analysis_table,
            fields,
            relpath=self.relpath)
        if not status:
            # Update flag to indicate problems with the
            # report
            self.status = False

    def report_extended_qc_metrics(self,project,section):
        """
        Report the extended QC metrics

        Creates a new subsection in the report document and
        adds content for reporting on:

        - RSeQC gene body coverage
        - Picard insert size metrics
        - Qualimap RNA-seq metrics

        Arguments:
          project (QCProject): project to report
          section (Section): document section to add
            the report to
        """
        # Create top-level container for extended metrics
        extended_metrics = section.add_subsection(
            "Extended QC metrics (experimental)",
            name='extended_qc_metrics')
        # Make a subsection for each organism
        for organism in project.organisms:
            extended_metrics_subsection = extended_metrics.add_subsection(
                "%s" % str(organism).title(),
                name='extended_qc_metrics_%s' % organism,
                css_classes=('info',))
            # Project-wide metrics
            # - RSeQC gene body coverage
            if 'rseqc_genebody_coverage' in project.outputs:
                rseqc_genebody_coverage = self.report_rseqc_genebody_coverage(
                    project,
                    organism,
                    extended_metrics_subsection)
            # Sample/Fastq group-level metrics
            # - Picard insert sizes
            # - Qualimap RNA-seq metrics
            fields = []
            if 'rseqc_infer_experiment' in project.outputs:
                fields.extend(['endedness_%s' % organism,
                               'strandedness_%s' % organism])
            if 'picard_insert_size_metrics' in project.outputs:
                fields.extend(['insert_size_metrics_%s' % organism,
                               'insert_size_histogram_%s' % organism])
            if 'qualimap_rnaseq' in project.outputs:
                fields.extend(['coverage_profile_along_genes_%s' % organism,
                               'reads_genomic_origin_%s' % organism,
                               'qualimap_rnaseq_report_%s' % organism])
            if fields:
                # Create a new summary table
                fields.insert(0,'sample')
                fields.insert(1,'bam_file')
                extended_metrics_table = self.add_summary_table(
                    project,
                    fields,
                    section=extended_metrics_subsection)
                # Add metrics for each sample
                for sample in project.samples:
                    reporter = SampleQCReporter(project,
                                                sample,
                                                qc_dir=project.qc_dir,
                                                fastq_attrs=project.fastq_attrs)
                    reporter.update_summary_table(extended_metrics_table,
                                                  fields=fields,
                                                  relpath=self.relpath)
        # Add an empty section to clear HTML floats
        clear = section.add_subsection(css_classes=("clear",))

    def report_rseqc_genebody_coverage(self,project,organism,section):
        """
        Add report of RSeQC gene body coverage to a document section

        Arguments:
          project (QCProject): parent project
          organism (str): name of organism to report coverage for
          section (Section): section to add the report to
        """
        # Determine location of QC artefacts
        qc_dir = self.fetch_qc_dir(project)
        # Top-level gene body coverage output dir
        genebody_coverage_dir = os.path.join(qc_dir,
                                             "rseqc_genebody_coverage",
                                             organism)
        if not os.path.exists(genebody_coverage_dir):
            # No outputs for specifed organism
            return None
        # Create a container for the outputs
        coverage = section.add_subsection(
            "Gene Body Coverage",
            name='rseqc_genebody_coverage_%s' %
            organism)
        # Acquire plot PNG
        png = os.path.join(genebody_coverage_dir,
                           "%s.geneBodyCoverage.curves.png"
                           % project.name)
        if not os.path.exists(png):
            # Missing PNG
            coverage.add(WarningIcon(),"No RSeQC gene body coverage plot")
        else:
            # Embed image in report
            if self.relpath:
                # Convert to relative path
                png = os.path.relpath(png,self.relpath)
            coverage.add("Gene body coverage from RSeQC for %s mapped "
                         "to %s" % (project.name,organism),
                         Img(png,
                             href=png,
                             title="Gene body coverage from RSeQC "
                             "(mapped to %s); click for PNG" %
                             organism,
                             name="gene_body_coverage_%s" % organism))
        # Return the subsection
        return coverage

    def fetch_qc_dir(self,project):
        """
        Return path to QC dir for reporting

        If a 'data directory' has been defined for this
        report then QC artefacts will have been copied
        to a project-specific subdirectory of that
        directory; otherwise QC artefacts will be in the
        QC directory of the original project.
        """
        # Determine location of QC artefacts
        if self.data_dir:
            return os.path.join(self.data_dir,
                                sanitize_name(project.id),
                                os.path.basename(project.qc_dir))
        else:
            return project.qc_dir

    def report_status(self):
        """
        Set the visibility of the "warnings" section
        """
        if self.status:
            # Turn off display of warnings section
            self.warnings.add_css_classes("hide")

class QCProject:
    """
    Gather information about the QC associated with a project

    Collects data about the QC for an AnalysisProject and
    makes it available via the following properties:

    - project: the AnalysisProject instance
    - qc_dir: the directory to examine for QC outputs
    - run_metadata: AnalysisDirMetadata instance with
      metadata from the parent run (if present)
    - processing_software: dictionary with information on
      software used to process the initial data

    Properties that shortcut to properties of the parent
    AnalysisProject:

    - name: project name
    - dirn: path to associated directory
    - comments: comments associated with the project
    - info: shortcut to the project's AnalysisProjectMetadata
      instance
    - qc_info: shortcut to the QC directory's
      AnalysisProjectQCDirInfo instance
    - fastq_attrs: class to use for extracting data from
      Fastq file names

    Properties based on artefacts detected in the QC
    directory:

    - fastqs: sorted list of Fastq names
    - reads: list of reads (e.g. 'r1', 'r2', 'i1' etc)
    - samples: sorted list of sample names extracted
      from Fastqs
    - bams: sorted list of BAM file names
    - multiplexed_samples: sorted list of sample names
      for multiplexed samples (e.g. 10x CellPlex)
    - organisms: sorted list of organism names
    - outputs: list of QC output categories detected (see
      below for valid values)
    - output_files: list of absolute paths to QC output
      files
    - software: dictionary with information on the
      QC software packages
    - stats: AttrtibuteDictionary with useful stats from
      across the project

    The 'stats' property has the following attributes:

    - max_seqs: maximum number of sequences across all
      Fastq files
    - min_sequence_length: minimum sequence length across
      all Fastq files
    - max_sequence_length: maximum sequence length across
      all Fastq files
    - max_sequence_length_read[READ]: maximum sequence
      length across all READ Fastqs (where READ is 'r1',
      'r2' etc)
    - min_sequence_length_read[READ]: minimum sequence
      length across all READ Fastqs (where READ is 'r1',
      'r2' etc)

    Valid values of the 'outputs' property are taken from
    the QCOutputs class.

    General properties about the project:

    - is_single_cell: True if the project has single cell
        data (10xGenomics, ICELL8 etc)

    Arguments:
      project (AnalysisProject): project to report QC for
      qc_dir (str): path to the QC output dir; relative
        path will be treated as a subdirectory of the
        project
    """
    def __init__(self,project,qc_dir=None):
        """
        Create a new QCProject instance
        """
        logger.debug("QCProject: project         : %s" % project.name)
        logger.debug("QCProject: project dir     : %s" % project.dirn)
        logger.debug("QCProject: qc_dir (initial): %s" % qc_dir)
        # Store project
        self.project = project
        # Sort out target QC dir
        if qc_dir is None:
            qc_dir = self.project.qc_dir
        else:
            if not os.path.isabs(qc_dir):
                qc_dir = os.path.join(self.project.dirn,
                                      qc_dir)
        self.qc_dir = qc_dir
        logger.debug("QCProject: qc_dir (final): %s" % self.qc_dir)
        # How to handle Fastq names
        self.fastq_attrs = self.project.fastq_attrs
        # Additional metrics
        self.stats = AttributeDictionary(
            max_seqs=None,
            min_sequence_length=None,
            max_sequence_length=None,
            min_sequence_length_read={},
            max_sequence_length_read={},
        )
        # Detect outputs
        self._detect_outputs()
        # Expose metadata from parent run
        self.run_metadata = AnalysisDirMetadata()
        run_metadata_file = os.path.join(
            os.path.dirname(os.path.abspath(self.project.dirn)),
            "metadata.info")
        if os.path.exists(run_metadata_file):
            print("Loading run metadata from %s" % run_metadata_file)
            self.run_metadata.load(run_metadata_file)
        else:
            logger.warning("Run metadata file '%s' not found"
                           % run_metadata_file)
        # Collect processing software metadata
        try:
            self.processing_software = ast.literal_eval(
                self.run_metadata.processing_software)
        except ValueError:
            self.processing_software = dict()
        if not self.processing_software:
            # Fallback to legacy metadata items
            try:
                self.processing_software['bcl2fastq'] = ast.literal_eval(
                    self.run_metadata.bcl2fastq_software)
            except ValueError:
                pass
            try:
                self.processing_software['cellranger'] = ast.literal_eval(
                    self.run_metadata.cellranger_software)
            except ValueError:
                pass
        # Expose project metadata
        self.info = self.project.info
        # Expose QC info
        self.qc_info = self.project.qc_info(self.qc_dir)

    @property
    def id(self):
        """
        Identifier for the project

        This is a string of the form:

        <RUN_ID>:<PROJECT_NAME>

        e.g. ``MINISEQ_201120#22:PJB``

        If the run id can't be determined then the name
        of the parent directory is used instead.
        """
        run_id = self.run_id
        if not run_id:
            run_id = os.path.basename(os.path.dirname(self.dirn))
        return "%s:%s" % (run_id,self.name)

    @property
    def name(self):
        """
        Name of project
        """
        return self.project.name

    @property
    def dirn(self):
        """
        Path to project directory
        """
        return self.project.dirn

    @property
    def comments(self):
        """
        Comments associated with the project
        """
        return self.project.info.comments

    @property
    def run_id(self):
        """
        Identifier for parent run

        This is the standard identifier constructed
        from the platform, datestamp and facility
        run number (e.g. ``MINISEQ_201120#22``).

        If an identifier can't be constructed then
        ``None`` is returned.
        """
        try:
            return run_reference_id(self.info['run'],
                                    platform=self.info['platform'],
                                    facility_run_number=
                                    self.run_metadata['run_number'])
        except (AttributeError,TypeError) as ex:
            logger.warning("Run reference ID can't be "
                           "determined: %s (ignored)" % ex)
            return None

    @property
    def is_single_cell(self):
        """
        Check whether project has single cell data
        """
        return (self.info['single_cell_platform'] is not None)

    def software_info(self,pkg,exclude_processing=False):
        """
        Get information on software package

        Arguments:
          pkg (str): name of software package to get
            information about
          exclude_processing (bool): if True then don't
            fall back to processing software information
            if package is not found (default: False, do
            fall back to checking processing software)

        Returns:
          String: software version information, or
            None if no information is stored.
        """
        # Acquire the value
        try:
            if self.software[pkg]:
                return ', '.join(self.software[pkg])
        except KeyError:
            if exclude_processing:
                # Don't check processing software
                return None
            try:
                return self.processing_software[pkg][2]
            except KeyError:
                # Not processing software
                return None
            except TypeError:
                # Not valid data
                return None

    def _detect_outputs(self):
        """
        Internal: determine which QC outputs are present
        """
        # Outsource to external class
        qc_outputs = QCOutputs(self.qc_dir,
                               fastq_attrs=self.fastq_attrs)
        # Fastqs
        self.fastqs = qc_outputs.fastqs
        # BAMs
        self.bams = qc_outputs.bams
        # Reads
        self.reads = qc_outputs.reads
        # Samples
        self.samples = sorted(list(
            set(qc_outputs.samples +
                [s.name for s in self.project.samples])),
                key=lambda s: split_sample_name(s))
        # Organisms
        self.organisms = qc_outputs.organisms
        # Fastq screens
        self.fastq_screens = qc_outputs.fastq_screens
        # Single library analyses reference data
        self.cellranger_references = qc_outputs.cellranger_references
        # Multiplexed samples
        self.multiplexed_samples = qc_outputs.multiplexed_samples
        # QC outputs
        self.outputs = qc_outputs.outputs
        # Software versions
        self.software = qc_outputs.software
        # Output files
        self.output_files = qc_outputs.output_files
        # Sequence length stats
        self.stats = AttributeDictionary(**qc_outputs.stats)

class SampleQCReporter:
    """
    Utility class for reporting the QC for a sample

    Provides the following properties:

    sample: name of the sample
    fastqs: list of the Fastqs associated with the sample
    reads: list of read ids e.g. ['r1','r2']
    fastq_groups: list of FastqGroupQCReporter instances from
      grouped Fastqs associated with the sample
    cellranger_count: list of CellrangerCount instances
      associated with the sample
    multiome_libraries: MultiomeLibraries instance with
      data for 10xGenomics libraries for the project (or
      None, if the libraries file doesn't exist)

    Provides the following methods:

    report_fastq_groups: add reports for Fastq groups in
      the sample to a document section
    update_summary_table: add lines to summary table for
      the sample
    update_single_library_table: add lines to the single
      library analysis summary table
    update_multiplexing_analysis_table: add lines to the
      multiplexing analysis summary table
    """
    # Summary fields provided at sample level
    sample_summary_fields = (
        'sample',
        'cellranger_count',
    )
    def __init__(self,project,sample,qc_dir=None,
                 fastq_attrs=AnalysisFastq):
        """
        Create a new SampleQCReporter

        Arguments:
          project (QCProject): project to report
          sample (str): name of sample to report
          qc_dir (str): path to the directory holding the
            QC artefacts
          fastq_attrs (BaseFastqAttrs): class for extracting
            data from Fastq names
        """
        self.sample = str(sample)
        self.fastq_groups = []
        self.cellranger_count = []
        self.cellranger_multi = []
        self.multiome_libraries = None
        # Location for QC artefacts
        if qc_dir:
            if not os.path.isabs(qc_dir):
                qc_dir = os.path.join(project.dirn)
        else:
            qc_dir = project.qc_dir
        # Group Fastqs associated with this sample
        self.fastqs = sorted(list(
            filter(lambda fq:
                   project.fastq_attrs(fq).sample_name == sample,
                   project.fastqs)))
        for fqs in group_fastqs_by_name(self.fastqs,fastq_attrs):
            self.fastq_groups.append(FastqGroupQCReporter(
                fqs,
                qc_dir=qc_dir,
                project=project,
                project_id=project.id,
                fastq_attrs=fastq_attrs))
        # Reads associated with the sample
        if self.fastq_groups:
            self.reads = [r for r in self.fastq_groups[0].reads]
        else:
            self.reads = []
        # 10x single library analyses
        cellranger_count_dir = os.path.join(qc_dir,
                                            "cellranger_count")
        if os.path.isdir(cellranger_count_dir):
            for dirn in walk(cellranger_count_dir):
                if os.path.isdir(dirn) and \
                   os.path.basename(dirn) == sample:
                    if not os.path.join(cellranger_count_dir,
                                        sample) == dirn:
                        # Extract version number and reference from
                        # the path
                        reference = dirn.split(os.sep)[-2]
                        version = dirn.split(os.sep)[-3]
                    else:
                        # No version or reference in path
                        reference = None
                        version = None
                    self.cellranger_count.append(
                        CellrangerCount(dirn,
                                        version=version,
                                        reference_data=reference))
        # 10x multiplexing analyses
        cellranger_multi_dir = os.path.join(qc_dir,
                                            "cellranger_multi")
        if os.path.isdir(cellranger_multi_dir):
            # Descend into version and reference subdirs
            for version in os.listdir(cellranger_multi_dir):
                for reference in os.listdir(
                        os.path.join(cellranger_multi_dir,version)):
                    # Add the cellranger multi information
                    self.cellranger_multi.append(
                        CellrangerMulti(os.path.join(cellranger_multi_dir,
                                                     version,
                                                     reference),
                                        version=version,
                                        reference_data=reference))
        # 10x multiome libraries
        multiome_libraries_file = os.path.join(
            project.dirn,
            "10x_multiome_libraries.info")
        if os.path.exists(multiome_libraries_file):
            self.multiome_libraries = MultiomeLibraries(
                multiome_libraries_file)

    def report_fastq_groups(self,sample_report,attrs,relpath=None):
        """
        Add reports for all Fastq groups to a document section

        Creates new subsections in 'sample_report' for each
        Fastq group, to report data on each Fastq file in the
        group.

        The following 'attributes' that can be reported for
        each Fastq are those available for the 'report' method
        of the 'FastqGroupQCReporter' class.

        By default all attributes are reported.

        Arguments:
          sample_report (Section): section to add the reports
            to
          attrs (list): optional list of custom 'attributes'
            to report
          relpath (str): if set then make link paths
            relative to 'relpath'
        """
        # Report each Fastq group
        for fastq_group in self.fastq_groups:
            # Report Fastq group in a subsection
            if attrs:
                fastq_group.report(sample_report,
                                   attrs=attrs,
                                   relpath=relpath)

    def update_summary_table(self,summary_table,fields=None,
                             sample_report=None,relpath=None):
        """
        Add lines to a summary table reporting a sample

        Creates new lines in 'summary_table' for the sample
        (one line per Fastq group), adding content for each
        specified field.

        See the 'get_value' method for a list of valid fields.

        Arguments:
          summary_table (Table): table to update
          fields (list): list of custom fields to report
          relpath (str): if set then make link paths
            relative to 'relpath'

        Returns:
          Boolean: True if report didn't contain any issues,
            False otherwise.
        """
        # Flag to indicate whether we're writing to the first
        # line of the summary table for this sample, and only
        # report name and sample-level metrics on the first line
        first_line = True
        # Reporting status
        has_problems = False
        # Report each Fastq group
        for fastq_group in self.fastq_groups:
            # Add line in summary table
            if first_line:
                # Only display sample name on first line
                if sample_report:
                    # Link to the sample report section
                    idx = summary_table.add_row(
                        sample=Link(self.sample,
                                    sample_report))
                else:
                    # No link
                    idx = summary_table.add_row(
                        sample=self.sample)
                # Add sample-level metrics
                for field in fields:
                    if field == "sample":
                        # Ignore: was already set
                        continue
                    elif field not in self.sample_summary_fields:
                        # Ignore non-sample level fields
                        continue
                    try:
                        value = self.get_value(field,
                                               relpath=relpath)
                        summary_table.set_value(idx,field,value)
                    except KeyError as ex:
                        # Ignore and carry on
                        logger.warning("Exception setting '%s' in "
                                       "summary table "
                                       "for sample %s: %s" %
                                       (field,self.sample,ex))
            else:
                # Don't report sample name on subsequent
                # lines in summary table
                idx = summary_table.add_row(sample="&nbsp;")
            status = fastq_group.update_summary_table(
                summary_table,idx=idx,
                fields=fields,
                relpath=relpath)
            if not status:
                has_problems = True
            # Update flag to indicate no longer on
            # first line for this sample
            first_line = False
        return (not has_problems)

    def update_single_library_table(self,package,single_library_table,
                                    fields=None,relpath=None):
        """
        Add lines to a table reporting single library analyses

        Creates new lines in 'single_library_table' for the
        sample (one line per single library analysis group),
        adding content for each specified field.

        Valid fields are any supported by the 'get_10x_value'
        method.

        Arguments:
          package (str): 10x package to report on
          summary_table (Table): table to update
          fields (list): list of custom fields to report
          relpath (str): if set then make link paths
            relative to 'relpath'

        Returns:
          Boolean: True if report didn't contain any issues,
            False otherwise.
        """
        # Flag to indicate whether we're writing to the first
        # line of the single library table for this sample, and
        # only report name on the first line
        first_line = True
        # Report each single library analysis
        for cellranger_count in self.cellranger_count:
            # Check package
            if cellranger_count.pipeline_name != package:
                # Skip
                continue
            # Shortcut to metrics
            metrics = cellranger_count.metrics
            web_summary = cellranger_count.web_summary
            # Add line in summary table
            if first_line:
                # Only display sample name on first line
                idx = single_library_table.add_row(
                    sample=self.sample)
                first_line = False
            else:
                idx = single_library_table.add_row()
            # Set values for fields
            for field in fields:
                if field == "sample":
                    continue
                value = self.get_10x_value(field,
                                           cellranger_count,
                                           metrics,
                                           web_summary,
                                           relpath=relpath)
                single_library_table.set_value(idx,field,value)
        # Finished
        return True

    def update_multiplexing_analysis_table(self,multiplexing_analysis_table,
                                           fields=None,relpath=None):
        """
        Add lines to a table reporting multiplexing analyses

        Creates new lines in 'multiplexing_analysis_table' for the
        sample (one line per multiplexed analysis), adding content
        for each specified field.

        Valid fields are any supported by the 'get_10x_value'
        method.

        Arguments:
          summary_table (Table): table to update
          fields (list): list of custom fields to report
          relpath (str): if set then make link paths
            relative to 'relpath'

        Returns:
          Boolean: True if report didn't contain any issues,
            False otherwise.
        """
        # Flag to indicate whether we're writing to the first
        # line of the multiplexing analysis table for this sample, and
        # only report name on the first line
        first_line = True
        # Report each multiplexing analysis
        for cellranger_multi in self.cellranger_multi:
            # Shortcut to metrics
            metrics = cellranger_multi.metrics(self.sample)
            web_summary = cellranger_multi.web_summary(self.sample)
            # Add line in summary table
            if first_line:
                # Only display sample name on first line
                idx = multiplexing_analysis_table.add_row(
                    sample=self.sample)
                first_line = False
            else:
                idx = multiplexing_analysis_table.add_row()
            # Set values for fields
            for field in fields:
                if field == "sample":
                    continue
                value = self.get_10x_value(field,
                                           cellranger_multi,
                                           metrics,
                                           web_summary,
                                           relpath=relpath)
                multiplexing_analysis_table.set_value(idx,field,value)
        # Finished
        return True

    def get_10x_value(self,field,cellranger_data,metrics,web_summary,
                      relpath=None):
        """
        Return the value for the specified field

        The following fields can be reported for each sample:

        Valid fields are:

        - 10x_cells
        - 10x_genes_per_cell
        - 10x_frac_reads_in_cell
        - 10x_fragments_per_cell
        - 10x_fragments_overlapping_targets
        - 10x_fragments_overlapping_peaks
        - 10x_tss_enrichment_score
        - 10x_atac_fragments_per_cell
        - 10x_gex_genes_per_cell
        - 10x_genes_detected
        - 10x_umis_per_cell
        - 10x_pipeline
        - 10x_reference
        - 10x_web_summary
        - linked_sample

        Arguments:
          field (str): name of the field to report; if the
            field is not recognised then KeyError is raised
          cellranger_data (object): parent CellrangerCount
            object for the sample
          metrics (MetricsSummary): summary metrics for the
            sample
          web_summary (str): path to the web_summary.html
            report for the sample
          relpath (str): if set then make link paths
            relative to 'relpath'
        """
        if field == "sample":
            return
        if field == "linked_sample":
            try:
                value = ','.join(
                    [split_sample_reference(s)[2]
                     for s in self.multiome_libraries.linked_samples(
                             self.sample)])
            except AttributeError:
                value = '?'
        elif field == "10x_cells":
            try:
                value = metrics.estimated_number_of_cells
            except AttributeError:
                try:
                    value = metrics.annotated_cells
                except AttributeError:
                    value = metrics.cells
            value = pretty_print_reads(value)
        elif field == "10x_reads_per_cell":
            try:
                value = metrics.mean_reads_per_cell
            except AttributeError:
                value = metrics.median_reads_per_cell
            value = pretty_print_reads(value)
        elif field == "10x_genes_per_cell":
            value = pretty_print_reads(metrics.median_genes_per_cell)
        elif field == "10x_frac_reads_in_cell":
            value = metrics.frac_reads_in_cells
        elif field == "10x_fragments_per_cell":
            value = pretty_print_reads(
                metrics.median_fragments_per_cell)
        elif field == "10x_fragments_overlapping_targets":
            try:
                value = "%.1f%%" % \
                        (metrics.frac_fragments_overlapping_targets
                         *100.0,)
            except AttributeError:
                value = 'N/A'
        elif field == "10x_fragments_overlapping_peaks":
            value = "%.1f%%" % \
                    (metrics.frac_fragments_overlapping_peaks*100.0,)
        elif field == "10x_tss_enrichment_score":
            try:
                value = "%.2f" % (metrics.tss_enrichment_score,)
            except AttributeError:
                value = 'N/A'
        elif field == "10x_atac_fragments_per_cell":
            value = pretty_print_reads(
                metrics.atac_median_high_quality_fragments_per_cell)
        elif field == "10x_gex_genes_per_cell":
            value = pretty_print_reads(
                metrics.gex_median_genes_per_cell)
        elif field == "10x_pipeline":
            if cellranger_data.version:
                value = "%s %s" % (cellranger_data.pipeline_name,
                                   cellranger_data.version)
            else:
                value = cellranger_data.pipeline_name
        elif field == "10x_genes_detected":
            value = pretty_print_reads(metrics.total_genes_detected)
        elif field == "10x_umis_per_cell":
            value = pretty_print_reads(metrics.median_umi_counts_per_cell)
        elif field == "10x_reference":
            value = os.path.basename(cellranger_data.\
                                     reference_data)
        elif field == "10x_web_summary":
            if relpath:
                web_summary = os.path.relpath(web_summary,
                                              relpath)
            value = Link(self.sample,web_summary)
        else:
            raise KeyError("'%s': unrecognised field for single "
                           "library analysis table" % field)
        return value

    def get_value(self,field,relpath=None):
        """
        Return the value for the specified field

        The following fields can be reported for each Fastq
        pair:

        - sample
        - cellranger_count

        Arguments:
          field (str): name of the field to report; if the
            field is not recognised then KeyError is raised
          relpath (str): if set then make link paths
            relative to 'relpath'
        """
        if field == "sample":
            value = self.sample
        elif field == "cellranger_count":
            if len(self.cellranger_count) == 1:
                # Report single summary
                cellranger_count = self.cellranger_count[0]
                web_summary = cellranger_count.web_summary
                if relpath:
                    web_summary = os.path.relpath(web_summary,
                                                  relpath)
                    value = Link(cellranger_count.sample_name,
                                 web_summary)
            else:
                # Report multiple single library outputs
                value = List()
                for cellranger_count in self.cellranger_count:
                    web_summary = cellranger_count.web_summary
                    if relpath:
                        web_summary = os.path.relpath(web_summary,
                                                      relpath)
                    pipeline_name = cellranger_count.pipeline_name
                    version = cellranger_count.version
                    reference = os.path.basename(
                        cellranger_count.reference_data)
                    cellranger = "%s%s" % (pipeline_name,
                                           (' %s' % version
                                            if version else ''))
                    value.add_item(Link("%s %s" % (cellranger,
                                                   reference),
                                        web_summary))
        else:
            raise KeyError("'%s': unrecognised field for summary "
                           "table" % field)
        return value

class FastqGroupQCReporter:
    """
    Utility class for reporting the QC for a Fastq group

    Provides the following properties:

    - reads: list of read ids e.g. ['r1','r2']
    - fastqs: dictionary mapping read ids to Fastq paths
    - reporters: dictionary mapping read ids to FastqQCReporter
      instances
    - paired_end: whether FastqGroup is paired end
    - bam: associated BAM file name
    - fastq_strand_txt: location of associated Fastq_strand
      output

    Provides the following methods:

    - infer_experiment: fetch data from RSeQC
      'infer_experiment.py'
    - ustrandednessplot: return mini-plot of strandedness
      data from RSeQC 'infer_experiment.py'
    - insert_size_metrics: fetch insert size metrics
    - uinsertsizeplot: return mini-plot of insert size
      histogram
    - qualimap_rnaseq: fetch Qualimap 'rnaseq' metrics
    - ucoverageprofileplot: return mini-plot of gene
      coverage profile
    - ugenomicoriginplot: return mini-plot of genomic origin
      of reads data
    - strandedness: fetch strandedness data for this group
    - ustrandplot: return mini-strand stats summary plot
    - report_strandedness: write report for strandedness
    - report: write report for the group
    - update_summary_table: add line to summary table for
      the group

    Arguments:
      fastqs (list): list of paths to Fastqs in the group
      qc_dir (str): path to the QC output dir; relative
        path will be treated as a subdirectory of the
        project
      project (QCProject): parent project
      project_id (str): identifier for the project
      fastq_attrs (BaseFastqAttrs): class for extracting
        data from Fastq names
    """
    def __init__(self,fastqs,qc_dir,project,project_id=None,
                 fastq_attrs=AnalysisFastq):
        """
        Create a new FastqGroupQCReporter
        """
        self.qc_dir = qc_dir
        self.project = project
        self.project_id = project_id
        self.fastq_attrs = fastq_attrs
        # Assign fastqs to reads
        self.fastqs = defaultdict(lambda: None)
        self.reporters = defaultdict(lambda: None)
        self.reads = set()
        for fastq in fastqs:
            fq = self.fastq_attrs(fastq)
            if fq.is_index_read:
                read = 'i'
            else:
                read = 'r'
            read = "%s%d" % (read,
                             fq.read_number
                             if fq.read_number is not None else 1)
            self.fastqs[read] = fastq
            self.reporters[read] = FastqQCReporter(fastq,
                                                   self.qc_dir,
                                                   self.project_id,
                                                   fastq_attrs=
                                                   self.fastq_attrs)
            self.reads.add(read)
        self.reads = sorted(list(self.reads))
        # Locate matching BAM file
        self.bam = None
        for fastq in fastqs:
            if fastq in self.project.bams:
                self.bam = fastq

    @property
    def paired_end(self):
        """
        True if pair consists of R1/R2 files
        """
        return ('r1' in self.fastqs and 'r2' in self.fastqs)

    @property
    def fastq_strand_txt(self):
        """
        Locate output from fastq_strand (None if not found)
        """
        fastq_strand_txt = None
        for fq in (self.fastqs['r1'],self.fastqs['r2']):
            if fq is not None:
                fastq_strand_txt = os.path.join(
                    self.qc_dir,
                    fastq_strand_output(fq))
                if os.path.isfile(fastq_strand_txt):
                    return fastq_strand_txt
        return None

    def infer_experiment(self,organism):
        """
        Return RSeQC infer_experiment.py data for organism
        """
        if self.bam:
            infer_experiment_log = os.path.join(
                self.qc_dir,
                "rseqc_infer_experiment",
                organism,
                "%s.infer_experiment.log" % self.bam)
            if os.path.isfile(infer_experiment_log):
                return InferExperiment(infer_experiment_log)
        return None

    def ustrandednessplot(self,organism,width=50,height=30):
        """
        Return a mini-plot for RSeQC strandness data
        """
        # Fetch data
        infer_experiment = self.infer_experiment(organism)
        fwd_data = [ self.infer_experiment(organism).forward,
                     1.0 - self.infer_experiment(organism).forward ]
        rvs_data = [ self.infer_experiment(organism).reverse,
                     1.0 - self.infer_experiment(organism).reverse ]
        # Create plot with two horizontal bars
        p = Plot(width,height)
        bar_length = width - 4
        bar_height = height/2 - 4
        # Start and end on x-axis (same for both)
        x1 = (width - bar_length)/2
        x2 = (width + bar_length)/2
        # Plot the 'reverse' bar
        y1 = (height/4 - bar_height/2)
        y2 = (height/4 + bar_height/2)
        p.bar(rvs_data,
              (x1,y1),(x2,y2),
              (RGB_COLORS.black,RGB_COLORS.lightgrey))
        # Plot the 'forward' bar
        y1 = y1 + height/2
        y2 = y2 + height/2
        p.bar(fwd_data,
              (x1,y1),(x2,y2),
              (RGB_COLORS.black,RGB_COLORS.lightgrey))
        return p.encoded_png()

    def insert_size_metrics(self,organism):
        """
        Return Picard insert size metrics for specified organism
        """
        if self.bam:
            insert_size_metrics_txt = os.path.join(
                self.qc_dir,
                "picard",
                organism,
                "%s.insert_size_metrics.txt" % self.bam)
            if os.path.isfile(insert_size_metrics_txt):
                return CollectInsertSizeMetrics(insert_size_metrics_txt)
        return None

    def uinsertsizeplot(self,organism):
        """
        Return a mini-plot with the Picard insert size histogram
        """
        insert_size_metrics = self.insert_size_metrics(organism)
        data = insert_size_metrics.histogram
        p = Plot(50,40)
        p.stripe(RGB_COLORS.lightgrey,RGB_COLORS.white)
        p.plot(data,RGB_COLORS.red,fill=True)
        return p.encoded_png()

    def qualimap_rnaseq(self,organism):
        """
        Return Qualimap 'rnaseq' outputs instance for specified organism
        """
        if self.bam:
            qualimap_rnaseq_dir = os.path.join(
                self.qc_dir,
                "qualimap-rnaseq",
                organism,
                self.bam)
            if os.path.isdir(qualimap_rnaseq_dir):
                return QualimapRnaseq(qualimap_rnaseq_dir)
        return None

    def ucoverageprofileplot(self,organism):
        """
        Return a mini-plot of the gene coverage profile
        """
        data = self.qualimap_rnaseq(organism).\
               raw_coverage_profile_along_genes_total
        p = Plot(50,40)
        p.stripe(RGB_COLORS.gainsboro,RGB_COLORS.white)
        p.plot(data,RGB_COLORS.red,interpolation="minmax")
        return p.encoded_png()

    def ugenomicoriginplot(self,organism,width=100,height=40):
        """
        Return a mini-barplot with the genomic origin of reads data
        """
        # Internal constants
        names = ['exonic',
                 'intronic',
                 'intergenic',
                 'overlapping exon',
                 'rRNA']
        colors = {
            'exonic': RGB_COLORS.maroon,
            'intronic': RGB_COLORS.green,
            'intergenic': RGB_COLORS.navyblue,
            'overlapping exon': RGB_COLORS.orange,
            'rRNA': RGB_COLORS.grey
        }
        # Fetch data
        data = self.qualimap_rnaseq(organism).reads_genomic_origin
        bar_length = width - 4
        bar_height = height - 10
        # Determine limits of the bar
        x1 = (width - bar_length)/2
        y1 = (height - bar_height)/2
        x2 = (width - bar_length)/2 + bar_length
        y2 = (height - bar_height)/2 + bar_height
        # Create plot
        p = Plot(width,height)
        p.bar([data[name][1] for name in names],
              (x1,y1),
              (x2,y2),
              [colors[name] for name in names])
        return p.encoded_png()

    def strandedness(self):
        """
        Return strandedness from fastq_strand.py
        """
        # Locate strandedness data
        txt = self.fastq_strand_txt
        # No file found
        if txt is None:
            return "!!!No strandedness data!!!"
        # Return stats
        strandedness = Fastqstrand(txt)
        output = []
        for genome in strandedness.genomes:
            try:
                ratio = "%.2f" % strandedness.stats[genome].ratio
            except TypeError:
                ratio = "NaN"
            output.append("%s: F%.2f%% | R%.2f%% (%s)" %
                          (genome,
                           strandedness.stats[genome].forward,
                           strandedness.stats[genome].reverse,
                           ratio))
        return "\n".join(output)

    def ustrandplot(self):
        """
        Return a mini-strand stats summary plot
        """
        txt = self.fastq_strand_txt
        return ustrandplot(txt,inline=True,dynamic=True)

    def report_strandedness(self,document):
        """
        Report the strandedness outputs to a document

        Creates a new subsection called "Strandedness"
        with a table of the strandedness determination
        outputs.

        Arguments:
          document (Section): section to add report to
        """
        rd = self.reads[0]
        strandedness_report = document.add_subsection(
            "Strandedness",
            name="strandedness_%s" % self.reporters[rd].safe_name)
        strandedness_report.add_css_classes("strandedness")
        # Locate strandedness
        txt = self.fastq_strand_txt
        # No file found
        if txt is None:
            strandedness_report.add(WarningIcon(),
                                    "No strandedness data available")
        else:
            strandedness = Fastqstrand(txt)
            strandedness_report.add("Strandedness statistics from "
                                    "<tt>fastq_strand %s</tt>:"
                                    % strandedness.version)
            strandedness_tbl = Table(("genome",
                                      "forward",
                                      "reverse",
                                      "ratio",
                                      "strand"),
                                     genome="Genome",
                                     forward="Forward %",
                                     reverse="Reverse %",
                                     ratio="Ratio (forward/reverse)",
                                     strand="Strandedness *")
            strandedness_tbl.add_css_classes("summary")
            for genome in strandedness.genomes:
                try:
                    ratio = "%.2f" % strandedness.stats[genome].ratio
                except TypeError:
                    ratio = "NaN"
                strandedness_tbl.add_row(
                    genome=genome,
                    forward=strandedness.stats[genome].forward,
                    reverse=strandedness.stats[genome].reverse,
                    ratio=ratio,
                    strand=strandedness.stats[genome].strandedness)
            strandedness_report.add(strandedness_tbl)
            strandedness_report.add("* Strandedness is 'reverse' if "
                                    "ratio &lt;0.2, 'forward' if ratio "
                                    "is &gt;5, 'unstranded' otherwise "
                                    "(or undetermined if 'NaN')")
        return strandedness_report

    def report(self,sample_report,attrs=None,relpath=None):
        """
        Add report for Fastq group to a document section

        Creates a new subsection in 'sample_report' for the
        Fastq pair, within which are additional subsections
        for each Fastq file.

        The following 'attributes' can be reported for each
        Fastq:

        - fastqc
        - fastq_screen
        - program versions

        By default all attributes are reported.

        Arguments:
          sample_report (Section): section to add the report
            to
          attrs (list): optional list of custom 'attributes'
            to report
          relpath (str): if set then make link paths
            relative to 'relpath'
        """
        # Attributes to report
        if attrs is None:
            attrs = ('fastqc',
                     'fastq_screen',
                     'program_versions',
                     'strandedness')
        # Add container section for Fastq pair
        fastqs_report = sample_report.add_subsection(css_classes=('fastqs',))
        # Create sections for individual Fastqs
        for read in self.reads:
            fq = self.reporters[read]
            if fq is None:
                continue
            fq_report = fastqs_report.add_subsection(fq.name,
                                                     name=fq.safe_name,
                                                     css_classes=('fastq',))
            # Add reports for each requested 'attribute'
            for attr in attrs:
                if attr == "fastqc":
                    # FastQC outputs
                    fq.report_fastqc(fq_report,relpath=relpath)
                elif attr == "fastq_screen":
                    # FastQScreens
                    fq.report_fastq_screens(fq_report,relpath=relpath)
                elif attr == "program_versions":
                    # Versions of programs used
                    fq.report_program_versions(fq_report)
                elif attr == "strandedness":
                    # Strandedness - handle separately
                    pass
                else:
                    raise KeyError("'%s': unrecognised reporting element "
                                   % attr)
        # Sections for pairwise data
        if "strandedness" in attrs:
            # Strandedness
            self.report_strandedness(fastqs_report)
        # Add an empty section to clear HTML floats
        clear = fastqs_report.add_subsection(css_classes=("clear",))

    def update_summary_table(self,summary_table,idx=None,fields=None,
                             relpath=None):
        """
        Add a line to a summary table reporting a Fastq group

        Creates a new line in 'summary_table' (or updates an
        existing line) for the Fastq pair, adding content for
        each specified field.

        See the 'get_value' method for a list of valid fields.

        Arguments:
          summary_table (Table): table to add the summary to
          idx (int): if supplied then indicates which existing
            table row to update (if None then a new row is
            appended)
          fields (list): list of custom fields to report
          relpath (str): if set then make link paths
            relative to 'relpath'

        Returns:
          Boolean: True if report didn't contain any issues,
            False otherwise.
        """
        # Flag indicating issues
        has_problems = False
        # Fields to report
        if fields is None:
            if self.paired_end:
                fields = ('fastqs',
                          'reads',
                          'read_lengths',
                          'read_lengths_dist_r1',
                          'read_lengths_dist_r2',
                          'boxplot_r1','boxplot_r2',
                          'fastqc_r1','fastqc_r2',
                          'screens_r1','screens_r2')
            else:
                fields = ('fastq',
                          'reads',
                          'read_lengths',
                          'read_lengths_dist_r1',
                          'boxplot_r1',
                          'fastqc_r1',
                          'screens_r1')
        else:
            # Drop fields for reads that aren't present
            # For example if there are a mixture of
            # single and paired-end Fastqs
            updated_fields = []
            for field in fields:
                if field[:-1].endswith("_r"):
                    read = field.split("_")[-1]
                    if read not in self.reads:
                        logging.warning("Dropping '%s' from summary table"
                                        % field)
                        continue
                updated_fields.append(field)
            fields = updated_fields
        # Add row to summary table
        if idx is None:
            idx = summary_table.add_row()
        # Populate with data
        for field in fields:
            if field in SampleQCReporter.sample_summary_fields:
                # Ignore sample-level fields
                continue
            try:
                # Populate fields in the table
                value = self.get_value(field,relpath=relpath)
                summary_table.set_value(idx,field,value)
            except KeyError as ex:
                # Assume this is an unrecognised field
                logger.error("%s" % ex)
                raise ex
            except Exception as ex:
                # Encountered an exception trying to acquire the value
                # for the field
                fastqs = ", ".join([self.reporters[r].name
                                    for r in self.reads])
                logger.warning("Exception setting '%s' in summary table "
                              "for Fastq group { %s }: %s"
                               % (field,
                                  fastqs,
                                  ex))
                # Put error value into the table
                summary_table.set_value(idx,field,
                                        WarningIcon("Unable to get value for "
                                                    "%s for %s" % (field,
                                                                   fastqs),
                                                     size=25))
                # Update flag
                has_problems = True
        # Return the status
        return (not has_problems)

    def get_value(self,field,relpath=None):
        """
        Return the value for the specified field

        The following fields can be reported for each Fastq
        group:

        - fastqs (if paired-end)
        - fastq (if single-end)
        - bam_file
        - reads
        - read_lengths
        - read_counts
        - sequence_duplication
        - adapter_content
        - read_lengths_dist_r1
        - read_lengths_dist_r2
        - read_lengths_dist_r3
        - boxplot_r1
        - boxplot_r2
        - boxplot_r3
        - fastqc_r1
        - fastqc_r2
        - fastqc_r3
        - screens_r1
        - screens_r2
        - screens_r3
        - strandedness
        - endedness_<organism>
        - strandedness_<organism>
        - insert_size_metrics_<organism>
        - insert_size_histogram_<organism>
        - reads_genomic_origin_<organism>
        - coverage_profile_along_genes_<organism>
        - qualimap_rnaseq_report_<organism>

        Arguments:
          field (str): name of the field to report; if the
            field is not recognised then KeyError is raised
          relpath (str): if set then make link paths
            relative to 'relpath'
        """
        if field == "fastq" or field == "fastqs":
            value = []
            for read in self.reads:
                value.append(Link(self.reporters[read].name,
                                  "#%s" % self.reporters[read].safe_name))
            value = "<br />".join([str(x) for x in value])
        elif field == "bam_file":
            value = self.bam
        elif field == "reads":
            if self.reporters[self.reads[0]].sequence_lengths:
                value = pretty_print_reads(
                    self.reporters[self.reads[0]].sequence_lengths.nreads)
            else:
                value = pretty_print_reads(
                    self.reporters[self.reads[0]].fastqc.data.basic_statistics(
                    'Total Sequences'))
        elif field == "read_counts":
            value = []
            title = []
            for read in self.reads:
                value.append(
                    Img(self.reporters[read].ureadcountplot(
                        max_reads=self.project.stats.max_seqs)))
                title.append(
                    "%s: %s reads\n* %s masked (%.1f%%)\n* %s padded (%.1f%%)"
                    % (read.upper(),
                       pretty_print_reads(
                           self.reporters[read].sequence_lengths.nreads),
                       pretty_print_reads(
                           self.reporters[read].sequence_lengths.nmasked),
                       self.reporters[read].sequence_lengths.frac_masked,
                       pretty_print_reads(
                           self.reporters[read].sequence_lengths.npadded),
                       self.reporters[read].sequence_lengths.frac_padded))
            value = "<div title=\"%s\">%s</div>" % ("\n\n".join(title),
                                                    "<br />".join(
                                                        [str(x)
                                                         for x in value]))
        elif field == "read_lengths":
            value = []
            for read in self.reads:
                if self.reporters[read].sequence_lengths:
                    value.append(
                        "%.1f&nbsp;(%s)" %
                        (self.reporters[read].sequence_lengths.mean,
                         self.reporters[read].sequence_lengths.range.replace('-','&#8209;')))
                else:
                    value.append(Link(
                        self.reporters[read].fastqc.data.basic_statistics(
                            'Sequence length'),
                        self.reporters[read].fastqc.summary.link_to_module(
                            'Sequence Length Distribution',
                            relpath=relpath)))
            value = "<br />".join([str(x) for x in value])
        elif field.startswith("read_lengths_dist_"):
            read = field.split('_')[-1]
            min_seq_len = self.project.stats.min_sequence_length_read[read]
            max_seq_len = self.project.stats.max_sequence_length_read[read]
            if (max_seq_len - min_seq_len) < 50:
                min_seq_len = max(0,max_seq_len - 50)
            if min_seq_len == max_seq_len:
                length_range = ""
            else:
                length_range = " (%s-%s)" % (min_seq_len,max_seq_len)
            value = Img(
                self.reporters[read].useqlenplot(
                    min_len=min_seq_len,
                    max_len=max_seq_len),
                href=self.reporters[read].fastqc.summary.link_to_module(
                    'Sequence Length Distribution',
                    relpath=relpath),
                title="%s: sequence length distribution%s\n(click for "
                "FastQC plot)" % (read.upper(),length_range))
        elif field == "sequence_duplication":
            value = []
            for read in self.reads:
                value.append(
                    Img(self.reporters[read].uduplicationplot(mode='dup'),
                        href=self.reporters[read].fastqc.summary.link_to_module(
                            'Sequence Duplication Levels',
                            relpath=relpath),
                        title="%s: percentage of sequences removed by "
                        "deduplication: %.2f%%\n(click for FastQC Sequence "
                        "Duplication Levels plot)" %
                        (read.upper(),
                         (100.0 - self.reporters[read].\
                          sequence_deduplication_percentage))))
            value = '<br />'.join([str(x) for x in value])
        elif field == "adapter_content":
            value = []
            for read in self.reads:
                value.append(
                    Img(self.reporters[read].uadapterplot(),
                        href=self.reporters[read].fastqc.summary.link_to_module(
                            'Adapter Content',
                            relpath=relpath),
                        title="%s fraction adapter content:\n%s" %
                        (read.upper(),
                         '\n'.join(
                             ["* %s: %.2f" %
                              (adapter,
                               self.reporters[read].adapters_summary[adapter])
                              for adapter in self.reporters[read].adapters]))))
            value = ''.join([str(x) for x in value])
        elif field.startswith("boxplot_"):
            read = field.split('_')[-1]
            value = Img(self.reporters[read].uboxplot(),
                        href="#boxplot_%s" %
                        self.reporters[read].safe_name,
                        title="%s: sequence quality distribution (click for "
                        "FastQC plot)" % read.upper())
        elif field.startswith("fastqc_"):
            read = field.split('_')[-1]
            value = Img(self.reporters[read].ufastqcplot(),
                        href="#fastqc_%s" %
                        self.reporters[read].safe_name,
                        title=self.reporters[read].fastqc_summary())
        elif field.startswith("screens_"):
            read = field.split('_')[-1]
            value = Img(self.reporters[read].uscreenplot(),
                        href="#fastq_screens_%s" %
                        self.reporters[read].safe_name)
        elif field == "strandedness":
            value = Img(self.ustrandplot(),
                        href="#strandedness_%s" %
                        self.reporters[self.reads[0]].safe_name,
                        title=self.strandedness())
        elif field.startswith("endedness_"):
            organism = field[len("endedness_"):]
            infer_experiment = self.infer_experiment(organism)
            if infer_experiment:
                value = ("PE" if infer_experiment.paired_end else "SE")
        elif field.startswith("strandedness_"):
            organism = field[len("strandedness_"):]
            infer_experiment = self.infer_experiment(organism)
            if infer_experiment:
                value = Img(self.ustrandednessplot(organism),
                            title="%s: F%.2f%% | R%.2f%% (U%.2f%%)"
                            % (self.bam,
                               infer_experiment.forward*100.0,
                               infer_experiment.reverse*100.0,
                               infer_experiment.unstranded*100.0))
        elif field.startswith("insert_size_metrics_"):
            organism = field[len("insert_size_metrics_"):]
            insert_size_metrics = self.insert_size_metrics(organism)
            if insert_size_metrics:
                insert_size = \
                    "%.2f bp (%.3f)" % \
                    (insert_size_metrics.metrics['MEAN_INSERT_SIZE'],
                     insert_size_metrics.metrics['STANDARD_DEVIATION'])
                insert_size_file = insert_size_metrics.metrics_file
                if relpath:
                    insert_size_file =  os.path.relpath(insert_size_file,
                                                        relpath)
                value = Link(insert_size,insert_size_file)
        elif field.startswith("insert_size_histogram_"):
            organism = field[len("insert_size_histogram_"):]
            histogram_pdf = os.path.join(
                os.path.dirname(
                    self.insert_size_metrics(organism).metrics_file),
                "%s.insert_size_histogram.pdf" % self.bam)
            if relpath:
                histogram_pdf = os.path.relpath(histogram_pdf,
                                                relpath)
            value = Img(self.uinsertsizeplot(organism),
                        title="%s: insert size histogram (click for PDF)"
                        % self.bam,
                        href=histogram_pdf)
        elif field.startswith("qualimap_rnaseq_report_"):
            organism = field[len("qualimap_rnaseq_report_"):]
            qualimap_outputs = self.qualimap_rnaseq(organism)
            if qualimap_outputs:
                qualimap_report = qualimap_outputs.html_report
                if relpath:
                    qualimap_report =  os.path.relpath(qualimap_report,
                                                       relpath)
                value = Link(DocumentIcon(title="%s: Qualimap rnaseq report"
                                          % self.bam,
                                          size=40),
                             target=qualimap_report)
        elif field.startswith("reads_genomic_origin_"):
            organism = field[len("reads_genomic_origin_"):]
            title = ["Genomic origin of reads:"]
            genomic_origin = \
                self.qualimap_rnaseq(organism).reads_genomic_origin
            for name in genomic_origin:
                title.append("- %s: %s (%s%%)" % (name,
                                                  genomic_origin[name][0],
                                                  genomic_origin[name][1]))
            link = self.qualimap_rnaseq(organism).\
                   link_to_output("Reads Genomic Origin",
                                  relpath=relpath)
            value = Img(self.ugenomicoriginplot(organism),
                        title='\n'.join(title),
                        href=link)
        elif field.startswith("coverage_profile_along_genes_"):
            organism = field[len("coverage_profile_along_genes_"):]
            link = self.qualimap_rnaseq(organism).\
                   link_to_output("Coverage Profile Along Genes (Total)",
                                  relpath=relpath)
            value = Img(self.ucoverageprofileplot(organism),
                        title="%s: Qualimap gene coverage profile" %
                        self.bam,
                        href=link)
        else:
            raise KeyError("'%s': unrecognised field for summary "
                           "table" % field)
        return value

class FastqQCReporter:
    """
    Provides interface to QC outputs for Fastq file

    Provides the following attributes:

    - name: basename of the Fastq
    - path: path to the Fastq
    - safe_name: name suitable for use in HTML links etc
    - sample_name: sample name derived from the Fastq basename
    - sequence_lengths: SeqLens instance
    - fastqc: Fastqc instance
    - fastq_screen.names: list of FastQScreen names
    - fastq_screen.SCREEN.description: description of SCREEN
    - fastq_screen.SCREEN.png: associated PNG file for SCREEN
    - fastq_screen.SCREEN.txt: associated TXT file for SCREEN
    - fastq_screen.SCREEN.version: associated version for SCREEN
    - program_versions.NAME: version of package NAME
    - adapters: list of adapters from Fastqc
    - adapters_summary: dictionary summarising adapter content

    Provides the following methods:

    - report_fastqc
    - report_fastq_screens
    - report_program_versions
    - useqlenplot
    - ureadcountplot
    - uboxplot
    - ufastqcplot
    - useqduplicationplot
    - uadapterplot
    - uscreenplot

    Arguments:
      fastq (str): path to Fastq file
      qc_dir (str): path to QC directory
      project_id (str): identifier for the parent project
      fastq_attrs (BaseFastqAttrs): class for extracting
        data from Fastq names
    """
    def __init__(self,fastq,qc_dir,project_id=None,
                 fastq_attrs=AnalysisFastq):
        """
        Create a new FastqQCReporter instance
        """
        # Source data
        logging.debug("******** FASTQQCREPORTER ************")
        self.name = os.path.basename(fastq)
        self.path = os.path.abspath(fastq)
        self.project_id = project_id
        self.safe_name = strip_ngs_extensions(self.name)
        if self.project_id:
            self.safe_name = "%s_%s" % (sanitize_name(self.project_id),
                                        self.safe_name)
        self.fastq_attrs = fastq_attrs
        # Sample name
        self.sample_name = fastq_attrs(self.name).sample_name
        logging.debug("Name  : %s" % self.name)
        logging.debug("Path  : %s" % self.path)
        logging.debug("ID    : %s" % self.project_id)
        logging.debug("Nsafe : %s" % self.safe_name)
        logging.debug("QCDir : %s" % qc_dir)
        # FastQC
        logging.debug("Fastqc: %s" %
                      os.path.join(qc_dir,fastqc_output(fastq)[0]))
        try:
            self.fastqc = Fastqc(os.path.join(
                qc_dir,fastqc_output(fastq)[0]))
        except Exception as ex:
            self.fastqc = None
        # Fastqscreen
        self.fastq_screen = AttributeDictionary()
        self.fastq_screen['names'] = list()
        fastq_base = self.fastq_attrs(fastq).basename
        # Legacy screens
        for screen_name in FASTQ_SCREENS:
            for f in list(filter(lambda f:
                                 f.startswith(fastq_base) and
                                 f.endswith("_%s_screen.txt" % screen_name),
                                 os.listdir(qc_dir))):
                png,txt = fastq_screen_output(fastq,screen_name,
                                              legacy=True)
                png = os.path.join(qc_dir,png)
                txt = os.path.join(qc_dir,txt)
                if os.path.exists(png) and os.path.exists(txt):
                    self.fastq_screen['names'].append(screen_name)
                    self.fastq_screen[screen_name] = AttributeDictionary()
                    self.fastq_screen[screen_name]["description"] = \
                                        screen_name.replace('_',' ').title()
                    self.fastq_screen[screen_name]["png"] = png
                    self.fastq_screen[screen_name]["txt"] = txt
                    self.fastq_screen[screen_name]["version"] = \
                                        Fastqscreen(txt).version
        # New-style screens (FASTQ_screen_SCREENNAME)
        for f in list(filter(lambda f:
                             f.startswith(fastq_base) and
                             f.endswith(".txt") and
                             "_screen_" in f,
                             os.listdir(qc_dir))):
            screen_name = f[:-len(".txt")].split("_screen_")[1]
            png,txt = fastq_screen_output(fastq,screen_name)
            png = os.path.join(qc_dir,png)
            txt = os.path.join(qc_dir,txt)
            if os.path.exists(png) and os.path.exists(txt):
                self.fastq_screen['names'].append(screen_name)
                self.fastq_screen[screen_name] = AttributeDictionary()
                self.fastq_screen[screen_name]["description"] = \
                                        screen_name.replace('_',' ').title()
                self.fastq_screen[screen_name]["png"] = png
                self.fastq_screen[screen_name]["txt"] = txt
                self.fastq_screen[screen_name]["version"] = \
                                        Fastqscreen(txt).version
        self.fastq_screen['names'] = sorted(self.fastq_screen['names'])
        # Sequence lengths
        try:
            self.sequence_lengths = SeqLens(
                os.path.join(qc_dir,
                             "%s_seqlens.json" % self.fastq_attrs(fastq)))
        except Exception as ex:
            self.sequence_lengths = None
        # Sequence deduplication percentage
        if self.fastqc is not None:
            self.sequence_deduplication_percentage = \
                self.fastqc.data.sequence_deduplication_percentage()
        else:
            self.sequence_deduplication_percentage = None
        # Adapters
        if self.fastqc is not None:
            self.adapters_summary = \
                self.fastqc.data.adapter_content_summary()
            self.adapters = self.adapters_summary.keys()
        else:
            self.adapters = None
        # Program versions
        self.program_versions = AttributeDictionary()
        if self.fastqc is not None:
            fastqc_version = self.fastqc.version
        else:
            fastqc_version = WarningIcon(size=20)
        self.program_versions['fastqc'] = fastqc_version
        fastq_screen_versions = list(
            set([self.fastq_screen[s].version
                 for s in self.fastq_screen.names]))
        if fastq_screen_versions:
            fastq_screen_versions = ','.join(sorted(fastq_screen_versions))
        else:
            fastq_screen_versions = WarningIcon(size=20)
        self.program_versions['fastq_screen'] = fastq_screen_versions

    def report_fastqc(self,document,relpath=None):
        """
        Report the FastQC outputs to a document

        Creates a new subsection called "FastQC" with
        a copy of the FastQC sequence quality boxplot and
        a summary table of the results from each FastQC
        module.

        Arguments:
          document (Section): section to add report to
          relpath (str): if set then make link paths
            relative to 'relpath'
        """
        fastqc_report = document.add_subsection("FastQC")
        if self.fastqc:
            # FastQC quality boxplot
            fastqc_report.add("Per base sequence quality boxplot:")
            boxplot = Img(self.fastqc.quality_boxplot(inline=True),
                          height=250,
                          width=480,
                          href=self.fastqc.summary.link_to_module(
                              'Per base sequence quality',
                              relpath=relpath),
                          name="boxplot_%s" % self.safe_name)
            fastqc_report.add(boxplot)
            # FastQC summary table
            fastqc_report.add("FastQC summary:")
            fastqc_tbl = Target("fastqc_%s" % self.safe_name)
            fastqc_report.add(fastqc_tbl,
                              self.fastqc.summary.html_table(relpath=relpath))
            if relpath:
                fastqc_html_report = os.path.relpath(self.fastqc.html_report,
                                                     relpath)
            else:
                fastqc_html_report = self.fastqc.html_report
            fastqc_report.add("%s for %s" % (Link("Full FastQC report",
                                                  fastqc_html_report),
                                             self.name))
        else:
            fastqc_report.add(WarningIcon(),"No FastQC data available")
        return fastqc_report

    def report_fastq_screens(self,document,relpath=None):
        """
        Report the FastQScreen outputs to a document

        Creates a new subsection called "Screens" with
        copies of the screen plots for each screen and
        links to the "raw" text files.

        Arguments:
          document (Section): section to add report to
          relpath (str): if set then make link paths
            relative to 'relpath'
        """
        screens_report = document.add_subsection("Screens",
                                                 name="fastq_screens_%s" %
                                                 self.safe_name)
        if not self.fastq_screen.names:
            screens_report.add(WarningIcon(),"No screens found")
            return screens_report
        raw_data = list()
        for name in self.fastq_screen.names:
            # Title for the screen
            screens_report.add(self.fastq_screen[name].description)
            # Gather data for screen
            png = self.fastq_screen[name].png
            txt = self.fastq_screen[name].txt
            if relpath:
                png_href = os.path.relpath(png,relpath)
                txt_href = os.path.relpath(txt,relpath)
            else:
                png_href = png
                txt_href = txt
            # Screen plot (PNG)
            if os.path.exists(png):
                screens_report.add(Img(encode_png(png),
                                       height=250,
                                       href=png_href))
            else:
                screens_report.add("!!!No FastqScreen plot available!!!")
            # Link to raw data (TXT)
            if os.path.exists(txt):
                raw_data.append(
                    Link(self.fastq_screen[name].description,
                         txt_href).html())
            else:
                screens_report.add("!!!No FastqScreen data available!!!")
        # Add links to raw data (TXT files)
        screens_report.add("Raw screen data: " +
                           "| ".join(raw_data))
        return screens_report

    def report_program_versions(self,document):
        """
        Report the program versions to a document

        Creates a new subsection called "Program versions"
        with a table listing the versions of the QC
        programs.

        Arguments:
          document (Section): section to add report to
          relpath (str): if set then make link paths
            relative to 'relpath'
        """
        versions = document.add_subsection("Program versions")
        programs = Table(("Program","Version"))
        programs.add_css_classes("programs","summary")
        programs.add_row(Program='fastqc',
                         Version=self.program_versions.fastqc)
        programs.add_row(Program='fastq_screen',
                         Version=self.program_versions.fastq_screen)
        versions.add(programs)
        return versions

    def fastqc_summary(self):
        """
        Return plaintext version of the FastQC summary
        """
        output = []
        for m in self.fastqc.summary.modules:
            output.append("%s: %s" %
                          (self.fastqc.summary.status(m),m))
        return "\n".join(output)

    def useqlenplot(self,max_len=None,min_len=None,height=None,
                    inline=True):
        """
        Return a mini-sequence length distribution plot

        Arguments:
          max_len (int): set the upper limit of the x-axis
          min_len (int): set the lower limit of the x-axis
          height (int): optionally set the plot height in pixels
          inline (bool): if True then return plot in format for
            inlining in HTML document
        """
        return useqlenplot(self.sequence_lengths.dist,
                           self.sequence_lengths.masked_dist,
                           min_len=min_len,
                           max_len=max_len,
                           height=height,
                           inline=inline)

    def ureadcountplot(self,max_reads=None,inline=True):
        """
        Return a mini-sequence composition plot

        Arguments:
          max_reads (int): if set then scale the reads for this
            Fastq against this value in the plot
          inline (bool): if True then return plot in format for
            inlining in HTML document
        """
        return ureadcountplot(self.sequence_lengths.nreads,
                              self.sequence_lengths.nmasked,
                              self.sequence_lengths.npadded,
                              max_reads=max_reads,
                              inline=inline)

    def uboxplot(self,inline=True):
        """
        Return a mini-sequence quality boxplot

        Arguments:
          inline (bool): if True then return plot in format for
            inlining in HTML document
        """
        return uboxplot(self.fastqc.data.path,inline=inline)

    def uduplicationplot(self,mode='dup',inline=True):
        """
        Return a mini-sequence duplication plot

        Arguments:
          mode (str): either 'dup' or 'dedup'
          inline (bool): if True then return plot in format for
            inlining in HTML document
        """
        return uduplicationplot(self.sequence_deduplication_percentage,
                                  mode=mode,
                                  inline=inline)

    def uadapterplot(self,height=40,multi_bar=False,inline=True):
        """
        Return a mini-adapter content summary plot

        Arguments:
          height (int): optionally set the plot height in pixels
          multi_bar (bool): if True then create a plot with one
            bar per adapter class (otherwise summarise adapter
            content on one bar in the plot)
          inline (bool): if True then return plot in format for
            inlining in HTML document
        """
        return uadapterplot(self.adapters_summary,self.adapters,
                            height=height,multi_bar=multi_bar,
                            inline=inline)

    def ufastqcplot(self,inline=True):
        """
        Return a mini-FastQC summary plot

        Arguments:
          inline (bool): if True then return plot in format for
            inlining in HTML document
        """
        return ufastqcplot(self.fastqc.summary.path,inline=inline)

    def uscreenplot(self,inline=True):
        """
        Return a mini-FastQScreen summary plot

        Arguments:
          inline (bool): if True then return plot in format for
            inlining in HTML document
        """
        screen_files = list()
        for name in self.fastq_screen.names:
            screen_files.append(self.fastq_screen[name].txt)
        return uscreenplot(screen_files,screen_width=30,inline=inline)

#######################################################################
# Functions
#######################################################################

def report(projects,title=None,filename=None,qc_dir=None,
           report_attrs=None,summary_fields=None,
           relative_links=False,use_data_dir=False,
           make_zip=False):
    """
    Report the QC for a project

    Arguments:
      projects (list): AnalysisProject instances to report
        QC for
      title (str): optional, specify title for the report
        (defaults to '<PROJECT_NAME>: QC report')
      filename (str): optional, specify path and name for
        the output report file (defaults to
        '<PROJECT_NAME>.qc_report.html')
      qc_dir (str): path to the QC output dir
      report_attrs (list): optional, list of elements to
        report for each Fastq pair
      summary_fields (list): optional, list of fields to
        report for each sample in the summary table
        relative_links (boolean): optional, if set to True
        then use relative paths for links in the report
        (default is to use absolute paths)
      use_data_dir (boolean): if True then copy QC
        artefacts to a data directory parallel to the
        output report file
      make_zip (boolean): if True then also create a ZIP
        archive of the QC report and outputs (default is
        not to create the ZIP archive)

    Returns:
      String: filename of the output HTML report.
    """
    # Primary project
    project = projects[0]
    # Set output destination
    if filename is None:
        filename = "%s.qc_report.html" % project.name
    # Use data directory
    if use_data_dir:
        data_dir = "%s_data" % os.path.splitext(filename)[0]
    else:
        data_dir = None
    # Use relative paths for links
    if relative_links:
        relpath = os.path.dirname(filename)
    else:
        relpath = None
    # Initialise report
    report = QCReport(projects,
                      title=title,
                      qc_dir=qc_dir,
                      report_attrs=report_attrs,
                      summary_fields=summary_fields,
                      data_dir=data_dir,
                      relpath=relpath)
    # Styles
    report.add_css_rule(QC_REPORT_CSS_STYLES)
    # Write the report
    report.write(filename)
    # Make ZIP file
    if make_zip:
        # Name for ZIP file
        out_dir = os.path.dirname(filename)
        basename = os.path.splitext(os.path.basename(filename))[0]
        run_name = project.info.run
        zip_prefix = "%s.%s%s" % (basename,
                                  project.name,
                                  '.%s' % run_name if run_name else '')
        logging.debug("ZIP prefix: %s" % zip_prefix)
        zip_name = os.path.join(out_dir,"%s.zip" % zip_prefix)
        logging.debug("ZIP file: %s" % zip_name)
        # Directory with QC artefacts
        if data_dir:
            qc_dir = data_dir
        elif qc_dir is None:
            qc_dir = project.qc_dir
        else:
            if not os.path.isabs(qc_dir):
                qc_dir = os.path.join(project.dirn,
                                      qc_dir)
        logging.debug("report: qc_dir  = %s" % qc_dir)
        # Create ZIP archive
        try:
            tmp_zip = "%s.tmp" % zip_name
            zip_file = ZipArchive(tmp_zip,
                                  relpath=os.path.dirname(qc_dir),
                                  prefix=zip_prefix)
            # Add the report
            zip_file.add_file(filename)
            # Add the QC outputs
            logging.debug("Adding QC outputs to ZIP archive '%s'" % zip_name)
            for f in report.output_files:
                # Output files are absolute paths
                logging.debug("-- Adding %s" % f)
                if os.path.exists(f):
                    zip_file.add(f)
                else:
                    raise Exception("Unable to add missing file '%s'" % f)
            zip_file.close()
            # Rename temporary ZIP to final location
            os.rename(tmp_zip,zip_name)
        except Exception as ex:
            # Remove temporary ZIP
            os.remove(tmp_zip)
            raise Exception("Failed to create ZIP archive '%s': %s" %
                            (zip_name,ex))
    # Return the output filename
    return filename

def pretty_print_reads(n):
    """
    Print the number of reads with commas at each thousand

    For example:

    >>> pretty_print_reads(10409789)
    10,409,789

    Arguments:
      n (int): number of reads

    Returns:
      String: representation with commas for every thousand.
    """
    n = str(int(n))[::-1]
    n0 = []
    while len(n) >= 3:
        n0.append(n[0:3])
        n = n[3:]
    if n: n0.append(n)
    return (','.join(n0))[::-1]

def sanitize_name(s,new_char='_'):
    """
    Replace 'unsafe' characters in HTML link targets

    Replaces 'unsafe' characters in a string used as a
    link target in an HTML document with underscore
    characters.

    Arguments:
      s (str): string to sanitize
      replace_with (str): string to replace 'unsafe'
        characters with (default: '_')
    """
    for c in ('#',':','/'):
        s = s.replace(c,new_char)
    return s
