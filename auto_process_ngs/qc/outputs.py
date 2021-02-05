#!/usr/bin/env python
#
#     qc/outputs: utilities to predict and check QC pipeline outputs
#     Copyright (C) University of Manchester 2019-2021 Peter Briggs
#
"""
Provides utility functions for QC outputs.

Provides the following functions:

- fastq_screen_output: get names for fastq_screen outputs
- fastqc_output: get names for FastQC outputs
- fastq_strand_output: get name for fastq_strand.py output
- cellranger_count_output: get names for cellranger count output
- cellranger_atac_count_output: get names for cellranger-atac count output
- cellranger_arc_count_output: get names for cellranger-arc count output
- check_illumina_qc_outputs: fetch Fastqs without illumina_qc.sh outputs
- check_fastq_strand_outputs: fetch Fastqs without fastq_strand.py outputs
- check_cellranger_count_outputs: fetch sample names without cellranger
  count outputs
- check_cellranger_atac_count_outputs: fetch sample names without
  cellranger-atac count outputs
- check_cellranger_arc_count_outputs: fetch sample names without
  cellranger-arc count outputs
- expected_outputs: return expected QC outputs for a project
"""

#######################################################################
# Imports
#######################################################################

import os
from bcftbx.qc.report import strip_ngs_extensions
import logging
from .constants import FASTQ_SCREENS
from ..fastq_utils import group_fastqs_by_name
from ..fastq_utils import remove_index_fastqs

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Functions
#######################################################################

def fastq_screen_output(fastq,screen_name):
    """
    Generate name of fastq_screen output files

    Given a Fastq file name and a screen name, the outputs from
    fastq_screen will look like:

    - {FASTQ}_{SCREEN_NAME}_screen.png
    - {FASTQ}_{SCREEN_NAME}_screen.txt

    Arguments:
       fastq (str): name of Fastq file
       screen_name (str): name of screen

    Returns:
       tuple: fastq_screen output names (without leading path)

    """
    base_name = "%s_%s_screen" % (strip_ngs_extensions(os.path.basename(fastq)),
                                  str(screen_name))

    return (base_name+'.png',base_name+'.txt')

def fastqc_output(fastq):
    """
    Generate name of FastQC outputs

    Given a Fastq file name, the outputs from FastQC will look
    like:

    - {FASTQ}_fastqc/
    - {FASTQ}_fastqc.html
    - {FASTQ}_fastqc.zip

    Arguments:
       fastq (str): name of Fastq file

    Returns:
       tuple: FastQC outputs (without leading paths)

    """
    base_name = "%s_fastqc" % strip_ngs_extensions(os.path.basename(fastq))
    return (base_name,base_name+'.html',base_name+'.zip')

def fastq_strand_output(fastq):
    """
    Generate name for fastq_strand.py output

    Given a Fastq file name, the output from fastq_strand.py
    will look like:

    - {FASTQ}_fastq_strand.txt

    Arguments:
       fastq (str): name of Fastq file

    Returns:
       tuple: fastq_strand.py output (without leading paths)

    """
    return "%s_fastq_strand.txt" % strip_ngs_extensions(
        os.path.basename(fastq))

def cellranger_count_output(project,sample_name=None,
                            prefix="cellranger_count"):
    """
    Generate list of 'cellranger count' outputs

    Given an AnalysisProject, the outputs from 'cellranger
    count' will look like:

    - {PREFIX}/{SAMPLE_n}/outs/metrics_summary.csv
    - {PREFIX}/{SAMPLE_n}/outs/web_summary.html

    for each SAMPLE_n in the project.

    If a sample name is supplied then outputs are limited
    to those for that sample

    Arguments:
      project (AnalysisProject): project to generate
        output names for
      sample_name (str): sample to limit outputs to
      prefix (str): directory for outputs (defaults
        to "cellranger_count")

    Returns:
       tuple: cellranger count outputs (without leading paths)
    """
    outputs = []
    # Metrics and web summary files
    for sample in project.samples:
        if sample_name and sample_name != sample.name:
            continue
        sample_count_dir = os.path.join(prefix,
                                        sample.name)
        for f in ("metrics_summary.csv",
                  "web_summary.html"):
            outputs.append(os.path.join(sample_count_dir,
                                        "outs",f))
    return tuple(outputs)

def cellranger_atac_count_output(project,sample_name=None,
                                 prefix="cellranger_count"):
    """
    Generate list of 'cellranger-atac count' outputs

    Given an AnalysisProject, the outputs from 'cellranger-atac
    count' will look like:

    - {PREFIX}/{SAMPLE_n}/outs/summary.csv
    - {PREFIX}/{SAMPLE_n}/outs/web_summary.html

    for each SAMPLE_n in the project.

    If a sample name is supplied then outputs are limited
    to those for that sample

    Arguments:
      project (AnalysisProject): project to generate
        output names for
      sample_name (str): sample to limit outputs to
      prefix (str): directory for outputs (defaults
        to "cellranger_count")

    Returns:
       tuple: cellranger count outputs (without leading paths)
    """
    outputs = []
    # Metrics and web summary files
    for sample in project.samples:
        if sample_name and sample_name != sample.name:
            continue
        sample_count_dir = os.path.join(prefix,
                                        sample.name)
        for f in ("summary.csv",
                  "web_summary.html"):
            outputs.append(os.path.join(sample_count_dir,
                                        "outs",f))
    return tuple(outputs)

def cellranger_arc_count_output(project,sample_name=None,
                                prefix="cellranger_count"):
    """
    Generate list of 'cellranger-arc count' outputs

    Given an AnalysisProject, the outputs from 'cellranger-arc
    count' will look like:

    - {PREFIX}/{SAMPLE_n}/outs/summary.csv
    - {PREFIX}/{SAMPLE_n}/outs/web_summary.html

    for each SAMPLE_n in the project.

    If a sample name is supplied then outputs are limited
    to those for that sample

    Arguments:
      project (AnalysisProject): project to generate
        output names for
      sample_name (str): sample to limit outputs to
      prefix (str): directory for outputs (defaults
        to "cellranger_count")

    Returns:
       tuple: cellranger count outputs (without leading paths)
    """
    outputs = []
    # Metrics and web summary files
    for sample in project.samples:
        if sample_name and sample_name != sample.name:
            continue
        sample_count_dir = os.path.join(prefix,
                                        sample.name)
        for f in ("summary.csv",
                  "web_summary.html"):
            outputs.append(os.path.join(sample_count_dir,
                                        "outs",f))
    return tuple(outputs)

def check_illumina_qc_outputs(project,qc_dir,qc_protocol=None):
    """
    Return Fastqs missing QC outputs from illumina_qc.sh

    Returns a list of the Fastqs from a project for which
    one or more associated outputs from `illumina_qc.sh`
    don't exist in the specified QC directory.

    Arguments:
      project (AnalysisProject): project to check the
        QC outputs for
      qc_dir (str): path to the QC directory (relative
        path is assumed to be a subdirectory of the
        project)
      qc_protocol (str): QC protocol to predict outputs
        for; if not set then defaults to standard QC
        based on ended-ness

    Returns:
      List: list of Fastq files with missing outputs.
    """
    if not os.path.isabs(qc_dir):
        qc_dir = os.path.join(project.dirn,qc_dir)
    fastqs = set()
    for fastq in remove_index_fastqs(project.fastqs,
                                     project.fastq_attrs):
        if qc_protocol in ('10x_scATAC',
                           '10x_Multiome_ATAC',):
            if project.fastq_attrs(fastq).read_number == 2:
                # Ignore the R2 reads for 10x single-cell ATAC
                continue
        # FastQC
        for output in [os.path.join(qc_dir,f)
                       for f in fastqc_output(fastq)]:
            if not os.path.exists(output):
                fastqs.add(fastq)
        # Fastq_screen
        if qc_protocol in ('singlecell',
                           '10x_scRNAseq',
                           '10x_snRNAseq',
                           '10x_Visium',
                           '10x_Multiome_GEX',):
            if project.fastq_attrs(fastq).read_number == 1:
                # No screens for R1 for single cell
                continue
        for screen in FASTQ_SCREENS:
            for output in [os.path.join(qc_dir,f)
                        for f in fastq_screen_output(fastq,screen)]:
                if not os.path.exists(output):
                    fastqs.add(fastq)
    return sorted(list(fastqs))

def check_fastq_strand_outputs(project,qc_dir,fastq_strand_conf,
                               qc_protocol=None):
    """
    Return Fastqs missing QC outputs from fastq_strand.py

    Returns a list of the Fastqs from a project for which
    one or more associated outputs from `fastq_strand.py`
    don't exist in the specified QC directory.

    Arguments:
      project (AnalysisProject): project to check the
        QC outputs for
      qc_dir (str): path to the QC directory (relative
        path is assumed to be a subdirectory of the
        project)
      fastq_strand_conf (str): path to a fastq_strand
        config file; strandedness QC outputs will be
        included unless the path is `None` or the
        config file doesn't exist. Relative path is
        assumed to be a subdirectory of the project
      qc_protocol (str): QC protocol to predict outputs
        for; if not set then defaults to standard QC
        based on ended-ness

    Returns:
      List: list of Fastq file "pairs" with missing
        outputs; pairs are (R1,R2) tuples, with 'R2'
        missing if only one Fastq is used for the
        strandedness determination.
    """
    # Sort out QC directory
    if not os.path.isabs(qc_dir):
        qc_dir = os.path.join(project.dirn,qc_dir)
    # Sort out fastq_strand config file
    if fastq_strand_conf is not None:
        if not os.path.isabs(fastq_strand_conf):
            fastq_strand_conf = os.path.join(project.dirn,
                                             fastq_strand_conf)
    if not os.path.exists(fastq_strand_conf):
        # No conf file, nothing to check
        return list()
    fastq_pairs = set()
    for fq_group in group_fastqs_by_name(
            remove_index_fastqs(project.fastqs,
                                project.fastq_attrs),
            fastq_attrs=project.fastq_attrs):
        # Strand stats output
        if qc_protocol in ('10x_scATAC',
                           '10x_Multiome_ATAC',):
            # Strand stats output based on R1/R3 pair
            fq_pair = (fq_group[0],fq_group[2])
        elif qc_protocol in ('singlecell',
                             '10x_scRNAseq',
                             '10x_snRNAseq',
                             '10x_Visium',
                             '10x_Multiome_GEX',):
            # Strand stats output based on R2
            fq_pair = (fq_group[1],)
        else:
            # All other protocols use R1 (single-end)
            # or R1/R2 (paired-end)
            if len(fq_group) > 1:
                fq_pair = (fq_group[0],fq_group[1])
            else:
                fq_pair = (fq_group[0],)
        output = os.path.join(qc_dir,
                              fastq_strand_output(fq_pair[0]))
        if not os.path.exists(output):
            fastq_pairs.add(fq_pair)
    return sorted(list(fastq_pairs))

def check_cellranger_count_outputs(project,qc_dir=None,
                                   prefix="cellranger_count"):
    """
    Return samples missing QC outputs from 'cellranger count'

    Returns a list of the samples from a project for which
    one or more associated outputs from `cellranger count`
    don't exist in the specified QC directory.

    Arguments:
      project (AnalysisProject): project to check the
        QC outputs for
      qc_dir (str): path to QC directory (if not the default
        QC directory for the project)
      prefix (str): directory for outputs (defaults
        to "cellranger_count")

    Returns:
      List: list of sample names with missing outputs
    """
    if qc_dir is None:
        qc_dir = project.qc_dir
    qc_dir = os.path.abspath(qc_dir)
    samples = set()
    for sample in project.samples:
        for output in cellranger_count_output(project,
                                              sample.name,
                                              prefix):
            if not os.path.exists(os.path.join(qc_dir,output)):
                samples.add(sample.name)
    return sorted(list(samples))

def check_cellranger_atac_count_outputs(project,qc_dir=None,
                                        prefix="cellranger_count"):
    """
    Return samples missing QC outputs from 'cellranger-atac count'

    Returns a list of the samples from a project for which
    one or more associated outputs from `cellranger-atac count`
    don't exist in the specified QC directory.

    Arguments:
      project (AnalysisProject): project to check the
        QC outputs for
      qc_dir (str): path to QC directory (if not the default
        QC directory for the project)
      prefix (str): directory for outputs (defaults
        to "cellranger_count")

    Returns:
      List: list of sample names with missing outputs
    """
    if qc_dir is None:
        qc_dir = project.qc_dir
    qc_dir = os.path.abspath(qc_dir)
    samples = set()
    for sample in project.samples:
        for output in cellranger_atac_count_output(project,
                                                   sample.name,
                                                   prefix):
            if not os.path.exists(os.path.join(qc_dir,output)):
                samples.add(sample.name)
    return sorted(list(samples))

def check_cellranger_arc_count_outputs(project,qc_dir=None,
                                       prefix="cellranger_count"):
    """
    Return samples missing QC outputs from 'cellranger-arc count'

    Returns a list of the samples from a project for which
    one or more associated outputs from `cellranger-arc count`
    don't exist in the specified QC directory.

    Arguments:
      project (AnalysisProject): project to check the
        QC outputs for
      qc_dir (str): path to QC directory (if not the default
        QC directory for the project)
      prefix (str): directory for outputs (defaults
        to "cellranger_count")

    Returns:
      List: list of sample names with missing outputs
    """
    if qc_dir is None:
        qc_dir = project.qc_dir
    qc_dir = os.path.abspath(qc_dir)
    samples = set()
    for sample in project.samples:
        if not os.path.exists(os.path.join(qc_dir,
                                           "libraries.%s.csv"
                                           % sample.name)):
            # Skip if there is no libraries.csv for the sample
            continue
        for output in cellranger_arc_count_output(project,
                                                  sample.name,
                                                  prefix):
            if not os.path.exists(os.path.join(qc_dir,output)):
                samples.add(sample.name)
    return sorted(list(samples))

def expected_outputs(project,qc_dir,fastq_strand_conf=None,
                     cellranger_version=None,cellranger_refdata=None,
                     qc_protocol=None):
    """
    Return expected QC outputs for a project

    Arguments:
      project (AnalysisProject): project to predict the
        QC outputs for
      qc_dir (str): path to the QC directory (relative
        path is assumed to be a subdirectory of the
        project)
      fastq_strand_conf (str): path to a fastq_strand
        config file; strandedness QC outputs will be
        included unless the path is `None` or the
        config file doesn't exist. Relative path is
        assumed to be a subdirectory of the project
      cellranger_version (str): version of the cellranger
        package used for single library analysis
      cellranger_refdata (str): path to a cellranger
        reference dataset; cellranger count outputs will
        be included for 10x protocols unless this path
        is set to `None`
      qc_protocol (str): QC protocol to predict outputs
        for; if not set then defaults to standard QC
        based on ended-ness

    Returns:
      List: list of full paths to the expected QC
        outputs from the project.
    """
    # Sort out QC directory
    if not os.path.isabs(qc_dir):
        qc_dir = os.path.join(project.dirn,qc_dir)
    # Sort out fastq_strand config file
    if fastq_strand_conf is not None:
        if not os.path.isabs(fastq_strand_conf):
            fastq_strand_conf = os.path.join(project.dirn,
                                             fastq_strand_conf)
    # Assemble the expected outputs
    outputs = set()
    for fastq in remove_index_fastqs(project.fastqs,
                                     project.fastq_attrs):
        if qc_protocol in ('10x_scATAC','10x_Multiome_ATAC') and \
           project.fastq_attrs(fastq).read_number == 2:
            # No outputs for R2 for 10x single cell ATAC-seq
            continue
        # FastQC
        for output in [os.path.join(qc_dir,f)
                       for f in fastqc_output(fastq)]:
            outputs.add(output)
        # Fastq_screen
        if (qc_protocol in ('singlecell',
                            '10x_scRNAseq',
                            '10x_snRNAseq',
                            '10x_Visium',
                            '10x_Multiome_GEX',)) \
            and project.fastq_attrs(fastq).read_number == 1:
            # No screens for R1 for single cell, Visium or
            # multiome GEX
            continue
        for screen in FASTQ_SCREENS:
            for output in [os.path.join(qc_dir,f)
                        for f in fastq_screen_output(fastq,screen)]:
                outputs.add(output)
    if fastq_strand_conf and os.path.exists(fastq_strand_conf):
        for fq_group in group_fastqs_by_name(
                remove_index_fastqs(project.fastqs,
                                    project.fastq_attrs),
                fastq_attrs=project.fastq_attrs):
            # Strand stats output
            if qc_protocol in ('singlecell',
                               '10x_scRNAseq',
                               '10x_snRNAseq',
                               '10x_Visium',
                               '10x_Multiome_GEX',):
                # Strand stats output based on R2
                output = os.path.join(qc_dir,
                                      fastq_strand_output(fq_group[1]))
            else:
                # Strand stats output based on R1
                output = os.path.join(qc_dir,
                                      fastq_strand_output(fq_group[0]))
            outputs.add(output)
    # Cellranger count output
    if cellranger_refdata is not None:
        prefix = "cellranger_count"
        if cellranger_version:
            prefix = os.path.join(prefix,
                                  cellranger_version,
                                  os.path.basename(cellranger_refdata))
        if qc_protocol in ('10x_scRNAseq','10x_snRNAseq',):
            for output in cellranger_count_output(project,
                                                  prefix=prefix):
                outputs.add(os.path.join(qc_dir,output))
        elif qc_protocol == '10x_scATAC':
            for output in cellranger_atac_count_output(project,
                                                       prefix=prefix):
                outputs.add(os.path.join(qc_dir,output))
        elif qc_protocol in ('10x_Multiome_ATAC','10x_Multiome_GEX',):
            # Only expect single library analysis output for 10x
            # single cell multiome data if there is a
            # 'libraries.<SAMPLE>.csv' file linking GEX and ATAC
            # datasets
            for sample in project.samples:
                if os.path.exists(os.path.join(qc_dir,
                                               "libraries.%s.csv"
                                               % sample.name)):
                    for output in cellranger_arc_count_output(
                            project,
                            sample_name=sample.name,
                            prefix=prefix):
                        outputs.add(os.path.join(qc_dir,output))
    return sorted(list(outputs))
