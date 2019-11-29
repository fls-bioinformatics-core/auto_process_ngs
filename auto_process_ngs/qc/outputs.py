#!/usr/bin/env python
#
#     qc/outputs: utilities to predict and check QC pipeline outputs
#     Copyright (C) University of Manchester 2019 Peter Briggs
#
"""
Provides utility functions for QC outputs.

Provides the following functions:

- fastq_screen_output: get names for fastq_screen outputs
- fastqc_output: get names for FastQC outputs
- fastq_strand_output: get name for fastq_strand.py output
- rseqc_gene_body_coverage_output: get names for RSeQC
  geneBody_coverage.py outputs
- rseqc_inner_distance_output: get names for RSeQC inner_distance.py
  outputs
- cellranger_count_output: get names for cellranger count output
- cellranger_atac_count_output: get names for cellranger-atac count output
- check_illumina_qc_outputs: fetch Fastqs without illumina_qc.sh outputs
- check_fastq_strand_outputs: fetch Fastqs without fastq_strand.py outputs
- check_rseqc_gene_body_coverage_outputs: fetch Fastqs without RSeQC
  geneBody_coverage.py outputs
- check_rseqc_inner_distance_outputs: fetch Fastqs without RSeQC
  inner_distance.py outputs
- check_cellranger_count_outputs: fetch sample names without cellranger
  count outputs
- check_cellranger_atac_count_outputs: fetch sample names without
  cellranger-atac count outputs
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
from ..utils import get_organism_list

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

def rseqc_gene_body_coverage_output(project,organism=None,
                                    include_heatmap=False):
    """
    Generate list of RSeQC geneBody_coverage.py outputs

    Given an AnalysisProject, the outputs from RSeQC's
    geneBody_coverage.py utility will look like:

    - {PREFIX}.geneBodyCoverage.txt
    - {PREFIX}.geneBodyCoverage.r
    - {PREFIX}.geneBodyCoverage.curves.png

    Optionally a heatmap will also be produced, if there are
    more than 3 input Fastqs/Fastq pairs:

    - {PREFIX}.geneBodyCoverage.heatMap.png

    'PREFIX' is generated as the project name, with the
    organism appended with a dot (if supplied), e.g.
    'PJB' or 'PJB.mouse'.

    Arguments:
      project (AnalysisProject): project to generate
        output names for
      organism (str): optional organism name
      include_heatmap (bool): if true then include the
        heatmap in the list of outputs

    Returns:
       tuple: RSeQC geneBody_coverage.py outputs (without
         leading paths)

    """
    prefix = project.name
    if organism:
        prefix += "." + organism
    outputs = ["%s.geneBodyCoverage.txt" % prefix,
               "%s.geneBodyCoverage.r" % prefix,
               "%s.geneBodyCoverage.curves.png" % prefix]
    if include_heatmap:
        outputs.append("%s.geneBodyCoverage.heatMap.png" % prefix)
    return tuple(outputs)

def rseqc_inner_distance_output(fastq,organism=None):
    """
    Generate list of RSeQC inner_distance.py outputs

    Given a Fastq file name, the output from inner_distance.py
    will look like:

    - {PREFIX}_inner_distance.txt
    - {PREFIX}_inner_distance_plot.r
    - {PREFIX}_inner_distance_plot_png.r
    - {PREFIX}_inner_distance_freq.txt
    - {PREFIX}_inner_distance_plot.pdf
    - {PREFIX}_inner_distance_plot.png

    'PREFIX' is generated as the basename for the Fastq file,
    with the organism appended with a dot (if supplied), e.g.
    'PJB_S1_R1_000' or 'PJB_S1_R1_000.mouse'.

    Arguments:
      fastq (str): name of Fastq file
      organism (str): optional organism name

    Returns:
       tuple: fastq_strand.py output (without leading paths)

    """
    prefix = os.path.basename(fastq)
    while prefix.split('.')[-1] in ('fastq','gz'):
        prefix = '.'.join(prefix.split('.')[:-1])
    if organism:
        prefix += "." + organism
    outputs = ["%s.inner_distance.txt" % prefix,
               "%s.inner_distance_plot.r" % prefix,
               "%s.inner_distance_plot_png.r" % prefix,
               "%s.inner_distance_freq.txt" % prefix,
               "%s.inner_distance_plot.pdf" % prefix,
               "%s.inner_distance_plot.png" % prefix]
    return tuple(outputs)

def cellranger_count_output(project,sample_name=None):
    """
    Generate list of 'cellranger count' outputs

    Given an AnalysisProject, the outputs from 'cellranger
    count' will look like:

    - cellranger_count/{SAMPLE_n}/outs/metrics_summary.csv
    - cellranger_count/{SAMPLE_n}/outs/web_summary.html

    for each SAMPLE_n in the project.

    If a sample name is supplied then outputs are limited
    to those for that sample

    Arguments:
      project (AnalysisProject): project to generate
        output names for
      sample_name (str): sample to limit outputs to

    Returns:
       tuple: cellranger count outputs (without leading paths)
    """
    outputs = []
    # Metrics and web summary files
    for sample in project.samples:
        if sample_name and sample_name != sample.name:
            continue
        sample_count_dir = os.path.join("cellranger_count",
                                        sample.name)
        for f in ("metrics_summary.csv",
                  "web_summary.html"):
            outputs.append(os.path.join(sample_count_dir,
                                        "outs",f))
    return tuple(outputs)

def cellranger_atac_count_output(project,sample_name=None):
    """
    Generate list of 'cellranger-atac count' outputs

    Given an AnalysisProject, the outputs from 'cellranger-atac
    count' will look like:

    - cellranger_count/{SAMPLE_n}/outs/summary.csv
    - cellranger_count/{SAMPLE_n}/outs/web_summary.html

    for each SAMPLE_n in the project.

    If a sample name is supplied then outputs are limited
    to those for that sample

    Arguments:
      project (AnalysisProject): project to generate
        output names for
      sample_name (str): sample to limit outputs to

    Returns:
       tuple: cellranger count outputs (without leading paths)
    """
    outputs = []
    # Metrics and web summary files
    for sample in project.samples:
        if sample_name and sample_name != sample.name:
            continue
        sample_count_dir = os.path.join("cellranger_count",
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
        if qc_protocol == '10x_scATAC':
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
                           '10x_snRNAseq',):
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
        if qc_protocol == '10x_scATAC':
            # Strand stats output based on R1/R3 pair
            fq_pair = (fq_group[0],fq_group[2])
        elif qc_protocol in ('singlecell',
                             '10x_scRNAseq',
                             '10x_snRNAseq',):
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

def check_rseqc_gene_body_coverage_outputs(project,qc_dir,
                                           organism=None):
    """
    Return Fastqs missing QC outputs from geneBody_coverage.py

    Returns a list of the Fastqs from a project for which
    one or more associated outputs from RSeqC's
    `geneBody_coverage.py` utility don't exist in the specified
    QC directory.

    Arguments:
      project (AnalysisProject): project to check the
        QC outputs for
      qc_dir (str): path to the QC directory (relative
        path is assumed to be a subdirectory of the
        project)
      organism (str): optional organism name

    Returns:
      List: list of Fastq files with missing outputs.

    """
    if not os.path.isabs(qc_dir):
        qc_dir = os.path.join(project.dirn,qc_dir)
    for output in rseqc_gene_body_coverage_output(project,
                                                  organism=organism):
        if not os.path.exists(os.path.join(qc_dir,output)):
            return project.fastqs
    return list()

def check_rseqc_inner_distance_outputs(project,qc_dir=None,organism=None):
    """
    Return Fastqs missing QC outputs from inner_distance.py

    Returns a list of the Fastqs from a project for which
    one or more associated outputs from RSeqC's
    `inner_distance.py` utility don't exist in the specified
    QC directory.

    Arguments:
      project (AnalysisProject): project to check the
        QC outputs for
      qc_dir (str): path to the QC directory (relative
        path is assumed to be a subdirectory of the
        project)
      organism (str): optional organism name

    Returns:
      List: list of Fastq (R1,R2) file pairs which have
        missing outputs.

    """
    if not os.path.isabs(qc_dir):
        qc_dir = os.path.join(project.dirn,qc_dir)
    fastq_pairs = list()
    for fq_group in group_fastqs_by_name(
            remove_index_fastqs(project.fastqs,
                                project.fastq_attrs),
            fastq_attrs=project.fastq_attrs):
        for output in rseqc_inner_distance_output(fq_group[0],
                                                  organism=organism):
            if not os.path.exists(os.path.join(qc_dir,output)):
                fastq_pairs.append(tuple(fq_group))
                break
    return sorted(list(fastq_pairs))

def check_cellranger_count_outputs(project,qc_dir=None):
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

    Returns:
      List: list of sample names with missing outputs
    """
    if qc_dir is None:
        qc_dir = project.qc_dir
    qc_dir = os.path.abspath(qc_dir)
    samples = set()
    for sample in project.samples:
        for output in cellranger_count_output(project,sample.name):
            if not os.path.exists(os.path.join(qc_dir,output)):
                samples.add(sample.name)
    return sorted(list(samples))

def check_cellranger_atac_count_outputs(project,qc_dir=None):
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

    Returns:
      List: list of sample names with missing outputs
    """
    if qc_dir is None:
        qc_dir = project.qc_dir
    qc_dir = os.path.abspath(qc_dir)
    samples = set()
    for sample in project.samples:
        for output in cellranger_atac_count_output(project,sample.name):
            if not os.path.exists(os.path.join(qc_dir,output)):
                samples.add(sample.name)
    return sorted(list(samples))

def expected_outputs(project,qc_dir,fastq_strand_conf=None,
                     cellranger_refdata=None,qc_protocol=None):
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
        if qc_protocol == '10x_scATAC' and \
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
                            '10x_snRNAseq',)) \
            and project.fastq_attrs(fastq).read_number == 1:
            # No screens for R1 for single cell
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
                               '10x_snRNAseq',):
                # Strand stats output based on R2
                output = os.path.join(qc_dir,
                                      fastq_strand_output(fq_group[1]))
            else:
                # Strand stats output based on R1
                output = os.path.join(qc_dir,
                                      fastq_strand_output(fq_group[0]))
            outputs.add(output)
    # RSeQC
    if project.info.library_type == 'RNA-seq':
        for organism in get_organism_list(project.info.organism):
            # Gene body coverage
            for output in rseqc_gene_body_coverage_output(
                    project,
                    organism):
                outputs.add(os.path.join(qc_dir,output))
            # Inner distance
            for fq_group in group_fastqs_by_name(
                    remove_index_fastqs(project.fastqs,
                                        project.fastq_attrs),
                    fastq_attrs=project.fastq_attrs):
                for output in rseqc_inner_distance_output(
                        fq_group[0],
                        organism):
                    outputs.add(os.path.join(qc_dir,output))
    # Cellranger count output
    if qc_protocol in ('10x_scRNAseq','10x_snRNAseq',) and \
       cellranger_refdata is not None:
        for output in cellranger_count_output(project):
            outputs.add(os.path.join(qc_dir,output))
    elif qc_protocol == '10x_scATAC':
        for output in cellranger_atac_count_output(project):
            outputs.add(os.path.join(qc_dir,output))
    return sorted(list(outputs))
