#!/usr/bin/env python
#
#     run_qc_cmd.py: implement auto process run_qc command
#     Copyright (C) University of Manchester 2018-2022 Peter Briggs
#
#########################################################################

#######################################################################
# Imports
#######################################################################

import os
import logging
from ..analysis import AnalysisProject
from ..command import Command
from ..qc.pipeline import QCPipeline
from ..qc.fastq_strand import build_fastq_strand_conf
from ..qc.protocols import determine_qc_protocol
from ..settings import fetch_reference_data
from ..utils import get_organism_list

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Command functions
#######################################################################

def run_qc(ap,projects=None,fastq_screens=None,
           fastq_subset=100000,nthreads=None,
           runner=None,fastq_dir=None,qc_dir=None,
           cellranger_exe=None,
           cellranger_chemistry='auto',
           cellranger_force_cells=None,
           cellranger_transcriptomes=None,
           cellranger_premrna_references=None,
           cellranger_extra_project_dirs=None,
           report_html=None,run_multiqc=True,
           working_dir=None,verbose=None,
           max_jobs=None,max_cores=None,
           batch_limit=None,enable_conda=None,
           conda_env_dir=None,poll_interval=None):
    """Run QC pipeline script for projects

    Run the illumina_qc.sh script to perform QC on projects.

    Note that if all QC outputs already exist for a project then
    the QC will *not* be run for that project.

    A subset of projects can be selected for QC by setting the
    'projects' argument to a name or pattern, only matching
    projects will be examined.

    Arguments:
      projects (str): specify a pattern to match one or more
        projects to run the QC for (default is to run QC for all
        projects)
      fastq_screens (dict): mapping of Fastq screen names to
        corresponding conf files, to use for contaminant screens
      fastq_subset (int): maximum size of subset of reads to use
        for FastQScreen, BAM file generation etc; set to zero or
        None to use all reads (default: 100000)
      nthreads (int): specify number of threads to run the QC jobs
        with (default: 1)
      runner (JobRunner): specify a non-default job runner to use
        for the QC jobs
      fastq_dir (str): specify the subdirectory to take the
        Fastq files from; will be used for all projects that are
        processed (default: 'fastqs')
      qc_dir (str): specify a non-standard directory to write the
        QC outputs to; will be used for all projects that are
        processed (default: 'qc')
      cellranger_exe (str): explicitly specify path to cellranger
        executable to use for 10xGenomics projects (default:
        determine appropriate executable automatically)
      cellranger_chemistry (str): assay configuration for
        10xGenomics scRNA-seq data (set to 'auto' to let cellranger
        determine this automatically; default: 'auto')
      cellranger_force_cells (int): override cell detection
        algorithm and set number of cells in 'cellranger' and
        'cellranger-atac' (set to 'None' to use built-in cell
        detection; default: 'None')
      cellranger_transcriptomes (dict): mapping of organism names
        to cellranger transcriptome reference data
      cellranger_premrna_references (dict): mapping of organism
        names to cellranger pre-mRNA reference data
      cellranger_extra_project_dirs (str): optional list of
        additional project dirs to use in single library analyses
      report_html (str): specify the name for the output HTML QC
        report (default: '<QC_DIR>_report.html')
      run_multiqc (bool): if True then run MultiQC at the end of
        the QC run (default)
      working_dir (str): path to a working directory (defaults to
         temporary directory in the current directory)
      verbose (bool): if True then report additional information
         for pipeline diagnostics
      max_jobs (int): maximum number of jobs that will be
        scheduled to run at one time (passed to the scheduler;
        default: no limit)
      max_cores (int): maximum number of cores available to
        the scheduler (default: no limit)
      batch_limit (int): if set then run commands in each task in
         batches, with the batch size set dyanmically so as not to
         exceed this limit
      enable_conda (bool): if True then use conda to resolve
        dependencies declared on tasks in the pipeline
      conda_env_dir (str): path to non-default directory for conda
        environments
      poll_interval (float): specifies non-default polling
        interval for scheduler used for running QC

    Returns:
      Integer: UNIX-style integer returncode where 0 = successful
        termination, non-zero indicates an error.
    """
    # Process project pattern matching
    if projects is None:
        project_pattern = '*'
        sample_pattern = '*'
    else:
        project_pattern = projects.split('/')[0]
        try:
            sample_pattern = projects.split('/')[1]
        except IndexError:
            sample_pattern = '*'
    # Get project dir data
    projects = ap.get_analysis_projects(project_pattern)
    # Check we have projects
    if len(projects) == 0:
        logger.warning("No projects found for QC analysis")
        return 1
    # Set up dictionaries for indices and reference data
    # STAR indexes
    star_indexes = fetch_reference_data(ap.settings,'star_index')
    # Annotation BEDs
    annotation_bed_files = fetch_reference_data(ap.settings,
                                                'annotation_bed')
    # Annotation GTFs
    annotation_gtf_files = fetch_reference_data(ap.settings,
                                                'annotation_gtf')
    # Set 10x cellranger reference data
    cellranger_transcriptomes_ = fetch_reference_data(
        ap.settings,
        'cellranger_reference')
    cellranger_premrna_references_ = fetch_reference_data(
        ap.settings,
        'cellranger_premrna_reference')
    cellranger_atac_references = fetch_reference_data(
        ap.settings,
        'cellranger_atac_reference')
    cellranger_multiome_references = fetch_reference_data(
        ap.settings,
        'cellranger_arc_reference')
    # Overload 10x transcriptome and pre-mRNA references
    if cellranger_transcriptomes:
        for organism in cellranger_transcriptomes:
            cellranger_transcriptomes_[organism] = \
                cellranger_transcriptomes[organism]
    cellranger_transcriptomes = cellranger_transcriptomes_
    if cellranger_premrna_references:
        for organism in cellranger_premrna_references:
             cellranger_premrna_references_[organism] = \
                cellranger_premrna_references_[organism]
    cellranger_premrna_references = cellranger_premrna_references_
    # Extra cellranger projects
    if cellranger_extra_project_dirs:
        cellranger_extra_projects = [AnalysisProject(d.strip())
                                     for d in
                                     cellranger_extra_project_dirs.split(',')]
    else:
        cellranger_extra_projects = None
    # Legacy FastqScreen naming convention
    legacy_screens = bool(ap.settings.qc.use_legacy_screen_names)
    # Set up runners
    if runner is None:
        default_runner = ap.settings.general.default_runner
        runners={
            'cellranger_count_runner': ap.settings.runners.cellranger_count,
            'cellranger_multi_runner': ap.settings.runners.cellranger_multi,
            'fastqc_runner': (ap.settings.runners.fastqc
                              if ap.settings.runners.fastqc
                              else ap.settings.runners.qc),
            'fastq_screen_runner': (ap.settings.runners.fastq_screen
                                    if ap.settings.runners.fastq_screen
                                    else ap.settings.runners.qc),
            'qualimap_runner': ap.settings.runners.qualimap,
            'rseqc_runner': ap.settings.runners.rseqc,
            'star_runner': ap.settings.runners.star,
            'verify_runner': default_runner,
            'report_runner': default_runner,
        }
    else:
        default_runner = runner
        runners={
            'cellranger_count_runner': runner,
            'cellranger_multi_runner': runner,
            'fastqc_runner': runner,
            'fastq_screen_runner': runner,
            'qualimap_runner': runner,
            'rseqc_runner': runner,
            'star_runner': runner,
            'verify_runner': runner,
            'report_runner': runner,
        }
    # Get environment modules
    envmodules = dict()
    for name in ('fastqc',
                 'fastq_screen',
                 'fastq_strand',
                 'cellranger',
                 'report_qc',):
        try:
            envmodules[name] = ap.settings.modulefiles[name]
        except KeyError:
            envmodules[name] = None
    # Conda dependency resolution
    if enable_conda is None:
        enable_conda = ap.settings.conda.enable_conda
    if conda_env_dir is None:
        conda_env_dir = ap.settings.conda.env_dir
    # Set scheduler parameters
    if poll_interval is None:
        poll_interval = ap.settings.general.poll_interval
    # Set up a master log directory and file
    ap.set_log_dir(ap.get_log_subdir('run_qc'))
    log_file = os.path.join(ap.log_dir,"run_qc.log")
    # Set up the QC for each project
    runqc = QCPipeline()
    for project in projects:
        # Determine the QC protocol
        protocol = determine_qc_protocol(project)
        runqc.add_project(project,
                          qc_dir=qc_dir,
                          fastq_dir=fastq_dir,
                          organism=project.info.organism,
                          qc_protocol=protocol,
                          sample_pattern=sample_pattern,
                          multiqc=True)
    # Collect the cellranger data and parameters
    cellranger_settings = ap.settings['10xgenomics']
    cellranger_jobmode = cellranger_settings.cellranger_jobmode
    cellranger_maxjobs = cellranger_settings.cellranger_maxjobs
    cellranger_mempercore = cellranger_settings.cellranger_mempercore
    cellranger_jobinterval = cellranger_settings.cellranger_jobinterval
    cellranger_localcores = cellranger_settings.cellranger_localcores
    cellranger_localmem = cellranger_settings.cellranger_localmem
    # Run the QC
    status = runqc.run(nthreads=nthreads,
                       fastq_screens=fastq_screens,
                       fastq_subset=fastq_subset,
                       star_indexes=star_indexes,
                       annotation_bed_files=annotation_bed_files,
                       annotation_gtf_files=annotation_gtf_files,
                       cellranger_transcriptomes=cellranger_transcriptomes,
                       cellranger_premrna_references=\
                       cellranger_premrna_references,
                       cellranger_atac_references=cellranger_atac_references,
                       cellranger_arc_references=cellranger_multiome_references,
                       cellranger_chemistry=cellranger_chemistry,
                       cellranger_force_cells=cellranger_force_cells,
                       cellranger_jobmode=cellranger_jobmode,
                       cellranger_maxjobs=cellranger_maxjobs,
                       cellranger_mempercore=cellranger_mempercore,
                       cellranger_jobinterval=cellranger_jobinterval,
                       cellranger_localcores=cellranger_localcores,
                       cellranger_localmem=cellranger_localmem,
                       cellranger_extra_projects=cellranger_extra_projects,
                       cellranger_exe=cellranger_exe,
                       log_file=log_file,
                       poll_interval=poll_interval,
                       max_jobs=max_jobs,
                       max_slots=max_cores,
                       batch_limit=batch_limit,
                       runners=runners,
                       default_runner=default_runner,
                       enable_conda=enable_conda,
                       conda_env_dir=conda_env_dir,
                       envmodules=envmodules,
                       working_dir=working_dir,
                       legacy_screens=legacy_screens,
                       verbose=verbose)
    return status
