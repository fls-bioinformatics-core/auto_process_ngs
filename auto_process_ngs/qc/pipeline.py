#!/usr/bin/env python
#
#     qc.pipeline.py: pipelines for running QC
#     Copyright (C) University of Manchester 2019-2022 Peter Briggs
#

"""
Pipeline components for running the QC pipeline.

Pipeline classes:

- QCPipeline

Pipeline task classes:

- SetupQCDirs
- GetSeqLengthStats
- CheckFastqScreenOutputs
- RunFastqScreen
- CheckFastQCOutputs
- RunFastQC
- SetupFastqStrandConf
- CheckFastqStrandOutputs
- RunFastqStrand
- DetermineRequired10xPackage
- GetCellrangerReferenceData
- MakeCellrangerArcCountLibraries
- GetCellrangerMultiConfig
- CheckCellrangerCountOutputs
- RunCellrangerCount
- RunCellrangerMulti
- SetCellCountFromCellrangerCount
- GetReferenceDataset
- GetBAMFiles
- RunRSeQCGenebodyCoverage
- RunPicardCollectInsertSizeMetrics
- CollateInsertSizes
- RunRSeQCInferExperiment
- RunQualimapRnaseq
- ReportQC

Also imports the following pipeline tasks:

- Get10xPackage

"""

######################################################################
# Imports
######################################################################

import os
import logging
import tempfile
import shutil
import random
from bcftbx.JobRunner import SimpleJobRunner
from bcftbx.TabFile import TabFile
from bcftbx.utils import mkdir
from bcftbx.utils import mkdirs
from bcftbx.utils import find_program
from bcftbx.ngsutils import getreads
from bcftbx.ngsutils import getreads_subset
from ..analysis import AnalysisFastq
from ..analysis import copy_analysis_project
from ..bcl2fastq.pipeline import Get10xPackage
from ..bcl2fastq.pipeline import FunctionParam
from ..command import Command
from ..fastq_utils import pair_fastqs_by_name
from ..fastq_utils import group_fastqs_by_name
from ..fastq_utils import remove_index_fastqs
from ..pipeliner import Pipeline
from ..pipeliner import PipelineTask
from ..pipeliner import PipelineFunctionTask
from ..pipeliner import PipelineCommandWrapper
from ..pipeliner import PipelineParam as Param
from ..pipeliner import ListParam
from ..pipeliner import PipelineFailure
from ..tenx_genomics_utils import CellrangerMultiConfigCsv
from ..tenx_genomics_utils import MultiomeLibraries
from ..tenx_genomics_utils import add_cellranger_args
from ..utils import get_organism_list
from .outputs import fastq_screen_output
from .outputs import fastq_strand_output
from .outputs import picard_collect_insert_size_metrics_output
from .outputs import rseqc_genebody_coverage_output
from .outputs import qualimap_rnaseq_output
from .outputs import check_fastq_screen_outputs
from .outputs import check_fastqc_outputs
from .outputs import check_fastq_strand_outputs
from .outputs import check_cellranger_count_outputs
from .outputs import check_cellranger_atac_count_outputs
from .outputs import check_cellranger_arc_count_outputs
from .picard import CollectInsertSizeMetrics
from .protocols import determine_qc_protocol
from .protocols import get_read_numbers
from .rseqc import InferExperiment
from .utils import get_bam_basename
from .utils import set_cell_count_for_project
from .verification import verify_project
from .fastq_strand import build_fastq_strand_conf
from .seqlens import get_sequence_lengths

# Module specific logger
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

######################################################################
# Pipeline classes
######################################################################

class QCPipeline(Pipeline):
    """
    Run the QC pipeline on one or more projects

    Pipeline to run QC on multiple projects.

    Example usage:

    >>> qc = QCPipeline()
    >>> qc.add_project(AnalysisProject("AB","./AB")
    >>> qc.add_project(AnalysisProject("CDE","./CDE")
    >>> qc.run()
    """
    def __init__(self):
        """
        Create a new QCPipeline instance
        """
        # Initialise the pipeline superclass
        Pipeline.__init__(self,name="QC")

        # Define parameters
        self.add_param('nthreads',type=int)
        self.add_param('annotation_bed_files',type=dict)
        self.add_param('annotation_gtf_files',type=dict)
        self.add_param('fastq_subset',type=int)
        self.add_param('cellranger_exe',type=str)
        self.add_param('cellranger_reference_dataset',type=str)
        self.add_param('cellranger_out_dir',type=str)
        self.add_param('cellranger_chemistry',type=str)
        self.add_param('cellranger_force_cells',type=int)
        self.add_param('cellranger_extra_projects',type=list)
        self.add_param('cellranger_transcriptomes',type=dict)
        self.add_param('cellranger_premrna_references',type=dict)
        self.add_param('cellranger_atac_references',type=dict)
        self.add_param('cellranger_arc_references',type=dict)
        self.add_param('cellranger_jobmode',type=str,value='local')
        self.add_param('cellranger_maxjobs',type=int)
        self.add_param('cellranger_mempercore',type=int)
        self.add_param('cellranger_jobinterval',type=int)
        self.add_param('cellranger_localcores',type=int)
        self.add_param('cellranger_localmem',type=int)
        self.add_param('fastq_screens',type=dict)
        self.add_param('star_indexes',type=dict)
        self.add_param('legacy_screens',type=bool,value=False)

        # Define runners
        self.add_runner('verify_runner')
        self.add_runner('fastq_screen_runner')
        self.add_runner('fastqc_runner')
        self.add_runner('star_runner')
        self.add_runner('qualimap_runner')
        self.add_runner('cellranger_runner')
        self.add_runner('report_runner')

        # Define module environment modules
        self.add_envmodules('fastqc')
        self.add_envmodules('fastq_screen')
        self.add_envmodules('fastq_strand')
        self.add_envmodules('cellranger')
        self.add_envmodules('report_qc')

    def add_project(self,project,qc_dir=None,organism=None,fastq_dir=None,
                    qc_protocol=None,report_html=None,multiqc=False,
                    sample_pattern=None,log_dir=None,
                    include_extended_metrics=False):
        """
        Add a project to the QC pipeline

        Arguments:
          project (AnalysisProject): project to run
            QC for
          qc_dir (str): directory for QC outputs (defaults
            to subdirectory 'qc' of project directory)
          organism (str): organism(s) for project
            (defaults to organism defined in project
            metadata)
          fastq_dir (str): directory holding Fastq files
            (defaults to primary fastq_dir in project)
          qc_protocol (str): QC protocol to use
          multiqc (bool): if True then also run MultiQC
            (default is not to run MultiQC)
          sample_pattern (str): glob-style pattern to
            match a subset of projects and samples (not
            implemented)
          log_dir (str): directory to write log files to
            (defaults to 'logs' subdirectory of the QC
            directory)
          include_extended_metrics (bool): if True then
            include generation of extended QC metrics in
            the pipeline (experimental feature)
        """
        ###################
        # Do internal setup
        ###################

        self.report("Adding project: %s" % project.name)

        # Determine QC protocol (if not set)
        if qc_protocol is None:
            qc_protocol = determine_qc_protocol(project)

        # Clone the supplied project
        project = copy_analysis_project(project,fastq_dir=fastq_dir)

        # Handle sample subsetting
        if sample_pattern is None:
            sample_pattern = '*'
        if sample_pattern != '*':
            raise NotImplementedError("Sample subsetting not "
                                      "supported")

        # Sort out the QC dir
        if qc_dir is None:
            qc_dir = 'qc'
        if not os.path.isabs(qc_dir):
            qc_dir = os.path.join(project.dirn,qc_dir)

        # Sort out other parameters
        if organism is None:
            organism = project.info.organism

        # Report details
        self.report("-- Protocol  : %s" % qc_protocol)
        self.report("-- Directory : %s" % project.dirn)
        self.report("-- Fastqs dir: %s" % project.fastq_dir)
        self.report("-- QC dir    : %s" % qc_dir)
        self.report("-- Library   : %s" % project.info.library_type)
        self.report("-- Organism  : %s" % organism)
        self.report("-- Report    : %s" % report_html)
        self.report("-- ExtendedQC: %s" % include_extended_metrics)

        ####################
        # Build the pipeline
        ####################

        if log_dir is None:
            log_dir = os.path.join(qc_dir,'logs')
        else:
            log_dir = os.path.abspath(log_dir)

        project_name = "%s%s" % (project.name,
                                 ":%s" % os.path.basename(fastq_dir)
                                 if fastq_dir is not None
                                 else '')

        # Set up QC dirs
        setup_qc_dirs = SetupQCDirs(
            "%s: set up QC directories" % project_name,
            project,
            qc_dir,
            log_dir=log_dir,
            qc_protocol=qc_protocol
        )
        self.add_task(setup_qc_dirs,
                      log_dir=log_dir)

        # Build a dictionary of QC metadata items to
        # update
        qc_metadata = dict(protocol=qc_protocol,
                           organism=organism,
                           fastq_dir=project.fastq_dir)

        # Update QC metadata
        update_qc_metadata = UpdateQCMetadata(
            "%s: update QC metadata" % project_name,
            project,
            qc_dir,
            qc_metadata,
            legacy_screens=self.params.legacy_screens)
        self.add_task(update_qc_metadata,
                      requires=(setup_qc_dirs,))

        # Verify QC
        verify_qc = VerifyQC(
            "%s: verify QC outputs" % project_name,
            project,
            qc_dir)
        self.add_task(verify_qc,
                      requires=(update_qc_metadata,),
                      runner=self.runners['verify_runner'],
                      log_dir=log_dir)

        # Make QC report
        report_qc = ReportQC(
            "%s: make QC report" % project_name,
            project,
            qc_dir,
            report_html=report_html,
            multiqc=multiqc,
            force=True
        )
        self.add_task(report_qc,
                      requires=(verify_qc,),
                      runner=self.runners['report_runner'],
                      envmodules=self.envmodules['report_qc'],
                      log_dir=log_dir)

        # Get Fastq sequence length statistics
        get_seq_lengths = GetSeqLengthStats(
            "%s: get sequence length statistics" %
            project_name,
            project,
            qc_dir,
            qc_protocol=qc_protocol,
            fastq_attrs=project.fastq_attrs)
        self.add_task(get_seq_lengths,
                      requires=(setup_qc_dirs,),
                      runner=self.runners['fastqc_runner'],
                      log_dir=log_dir)
        verify_qc.requires(get_seq_lengths)

        # Check outputs for FastqScreen
        check_fastq_screen = CheckFastqScreenOutputs(
            "%s: check FastqScreen outputs" %
            project_name,
            project,
            qc_dir,
            self.params.fastq_screens,
            qc_protocol=qc_protocol,
            legacy=self.params.legacy_screens,
            verbose=self.params.VERBOSE
        )
        self.add_task(check_fastq_screen,
                      requires=(setup_qc_dirs,),
                      runner=self.runners['verify_runner'],
                      log_dir=log_dir)

        # Run FastqScreen
        run_fastq_screen = RunFastqScreen(
            "%s: FastqScreen" % project_name,
            check_fastq_screen.output.fastqs,
            qc_dir,
            self.params.fastq_screens,
            subset=self.params.fastq_subset,
            nthreads=self.params.nthreads,
            qc_protocol=qc_protocol,
            fastq_attrs=project.fastq_attrs,
            legacy=self.params.legacy_screens
        )
        self.add_task(run_fastq_screen,
                      requires=(check_fastq_screen,),
                      runner=self.runners['fastq_screen_runner'],
                      envmodules=self.envmodules['fastq_screen'],
                      log_dir=log_dir)
        qc_metadata['fastq_screens'] = self.params.fastq_screens
        qc_metadata['legacy_screens'] = self.params.legacy_screens
        verify_qc.requires(run_fastq_screen)

        # Check outputs for FastQC
        check_fastqc = CheckFastQCOutputs(
            "%s: check FastQC outputs" %
            project_name,
            project,
            qc_dir,
            qc_protocol=qc_protocol,
            verbose=self.params.VERBOSE
        )
        self.add_task(check_fastqc,
                      requires=(setup_qc_dirs,),
                      runner=self.runners['verify_runner'],
                      log_dir=log_dir)

        # Run FastqQC
        run_fastqc = RunFastQC(
            "%s: FastQC" % project_name,
            check_fastqc.output.fastqs,
            qc_dir,
            nthreads=self.params.nthreads
        )
        self.add_task(run_fastqc,
                      requires=(check_fastqc,),
                      runner=self.runners['fastqc_runner'],
                      envmodules=self.envmodules['fastqc'],
                      log_dir=log_dir)
        verify_qc.requires(run_fastqc)

        # Set up fastq_strand.conf file
        setup_fastq_strand_conf = SetupFastqStrandConf(
            "%s: conf file for strandedness (fastq_strand)" %
            project_name,
            project,
            qc_dir=qc_dir,
            organism=organism,
            star_indexes=self.params.star_indexes
        )
        self.add_task(setup_fastq_strand_conf,
                      requires=(setup_qc_dirs,),
                      log_dir=log_dir)

        # Check outputs for fastq_strand.py
        check_fastq_strand = CheckFastqStrandOutputs(
            "%s: check strandedness outputs (fastq_strand)" %
            project_name,
            project,
            qc_dir,
            setup_fastq_strand_conf.output.fastq_strand_conf,
            qc_protocol=qc_protocol,
            verbose=self.params.VERBOSE
        )
        self.add_task(check_fastq_strand,
                      requires=(setup_fastq_strand_conf,),
                      runner=self.runners['verify_runner'],
                      log_dir=log_dir)

        # Run fastq_strand.py
        run_fastq_strand = RunFastqStrand(
            "%s: strandedness (fastq_strand.py)" %
            project_name,
            check_fastq_strand.output.fastq_pairs,
            qc_dir,
            setup_fastq_strand_conf.output.fastq_strand_conf,
            fastq_strand_subset=self.params.fastq_subset,
            nthreads=self.params.nthreads,
            qc_protocol=qc_protocol
        )
        self.add_task(run_fastq_strand,
                      requires=(check_fastq_strand,),
                      runner=self.runners['star_runner'],
                      envmodules=self.envmodules['fastq_strand'],
                      log_dir=log_dir)
        verify_qc.requires(run_fastq_strand)

        if qc_protocol in ("10x_scRNAseq",
                           "10x_snRNAseq",
                           "10x_scATAC",
                           "10x_Multiome_ATAC",
                           "10x_Multiome_GEX",):
            # Run cellranger* count
            run_cellranger_count = self.add_cellranger_count(
                project_name,
                project,
                qc_dir,
                organism,
                fastq_dir,
                qc_protocol,
                chemistry=self.params.cellranger_chemistry,
                force_cells=self.params.cellranger_force_cells,
                reference_dataset=self.params.cellranger_reference_dataset,
                extra_projects=self.params.cellranger_extra_projects,
                log_dir=log_dir,
                required_tasks=(setup_qc_dirs,))

            # Update metadata
            qc_metadata['cellranger_version'] = \
                    run_cellranger_count.output.cellranger_version
            qc_metadata['cellranger_refdata'] = \
                    run_cellranger_count.output.cellranger_refdata
            update_qc_metadata.requires(run_cellranger_count)

            # Set cell count
            set_cellranger_cell_count = SetCellCountFromCellrangerCount(
                "%s: set cell count from single library analysis" %
                project_name,
                project,
                qc_dir
            )
            self.add_task(set_cellranger_cell_count,
                          requires=(run_cellranger_count,
                                    update_qc_metadata),)
            verify_qc.requires(set_cellranger_cell_count)

            # Extra protocols for multiome
            if qc_protocol == "10x_Multiome_ATAC":
                # See https://kb.10xgenomics.com/hc/en-us/articles/360061165691
                # Need to set the chemistry to indicate it's
                # multiome ATAC data
                run_cellranger_count = self.add_cellranger_count(
                    project_name,
                    project,
                    qc_dir,
                    organism,
                    fastq_dir,
                    qc_protocol="10x_scATAC",
                    chemistry="ARC-v1",
                    force_cells=self.params.cellranger_force_cells,
                    reference_dataset=\
                    self.params.cellranger_reference_dataset,
                    log_dir=log_dir,
                    required_tasks=(setup_qc_dirs,))
                verify_qc.requires(run_cellranger_count)

            elif qc_protocol == "10x_Multiome_GEX":
                # See https://kb.10xgenomics.com/hc/en-us/articles/360059656912
                # Need to set the chemistry to indicate it's
                # multiome snRNA-seq data
                run_cellranger_count = self.add_cellranger_count(
                    project_name,
                    project,
                    qc_dir,
                    organism,
                    fastq_dir,
                    qc_protocol="10x_snRNAseq",
                    chemistry="ARC-v1",
                    force_cells=self.params.cellranger_force_cells,
                    reference_dataset=\
                    self.params.cellranger_reference_dataset,
                    log_dir=log_dir,
                    required_tasks=(setup_qc_dirs,))
                verify_qc.requires(run_cellranger_count)

        elif qc_protocol in ("10x_CellPlex",):

            # Tasks 'check_cellranger_multi' depends on
            check_cellranger_multi_requires = []

            # Locate cellranger
            required_cellranger = DetermineRequired10xPackage(
                "%s: determine required 'cellranger' package" %
                project_name,
                qc_protocol,
                self.params.cellranger_exe)
            self.add_task(required_cellranger)

            get_cellranger = Get10xPackage(
                "%s: get information on cellranger" % project_name,
                require_package=\
                required_cellranger.output.require_cellranger)
            self.add_task(get_cellranger,
                          requires=(required_cellranger,),
                          envmodules=self.envmodules['cellranger'])
            check_cellranger_multi_requires.append(get_cellranger)
            qc_metadata['cellranger_version'] = \
                    get_cellranger.output.package_version
            update_qc_metadata.requires(get_cellranger)

            # Locate config.csv file for 'cellranger multi'
            get_cellranger_multi_config = GetCellrangerMultiConfig(
                "%s: get config file for 'cellranger multi'" %
                project_name,
                project,
                qc_dir
            )
            self.add_task(get_cellranger_multi_config,
                          requires=(setup_qc_dirs,),
                          log_dir=log_dir)
            check_cellranger_multi_requires.append(
                get_cellranger_multi_config)
            qc_metadata['cellranger_refdata'] = \
                    get_cellranger_multi_config.output.reference_data_path
            update_qc_metadata.requires(get_cellranger_multi_config)

            # Parent directory for cellranger multi outputs
            # Set to project directory unless the 'cellranger_out_dir'
            # parameter is set
            cellranger_out_dir = FunctionParam(
                lambda out_dir,project_dir:
                out_dir if out_dir is not None else project_dir,
                self.params.cellranger_out_dir,
                project.dirn)

            # Run cellranger multi
            run_cellranger_multi = RunCellrangerMulti(
                "%s: analyse cell multiplexing data (cellranger multi)" %
                project_name,
                project,
                get_cellranger_multi_config.output.config_csv,
                get_cellranger_multi_config.output.samples,
                get_cellranger_multi_config.output.reference_data_path,
                cellranger_out_dir,
                qc_dir=qc_dir,
                working_dir=self.params.WORKING_DIR,
                cellranger_exe=get_cellranger.output.package_exe,
                cellranger_version=get_cellranger.output.package_version,
                cellranger_jobmode=self.params.cellranger_jobmode,
                cellranger_maxjobs=self.params.cellranger_maxjobs,
                cellranger_mempercore=self.params.cellranger_mempercore,
                cellranger_jobinterval=self.params.cellranger_jobinterval,
                cellranger_localcores=self.params.cellranger_localcores,
                cellranger_localmem=self.params.cellranger_localmem,
                qc_protocol=qc_protocol
            )
            self.add_task(run_cellranger_multi,
                          requires=(get_cellranger,
                                    get_cellranger_multi_config,),
                          runner=self.runners['cellranger_runner'],
                          envmodules=self.envmodules['cellranger'],
                          log_dir=log_dir)

            # Set cell count
            set_cellranger_cell_count = SetCellCountFromCellrangerCount(
                "%s: set cell count from cell multiplexing analysis" %
                project_name,
                project,
                qc_dir
            )
            self.add_task(set_cellranger_cell_count,
                          requires=(run_cellranger_multi,),)
            verify_qc.requires(set_cellranger_cell_count)

            # Extra protocol for cellplex data
            # (NB only run on the GEX "samples")
            run_cellranger_count = self.add_cellranger_count(
                project_name,
                project,
                qc_dir,
                organism,
                fastq_dir,
                qc_protocol="10x_scRNAseq",
                chemistry=self.params.cellranger_chemistry,
                force_cells=self.params.cellranger_force_cells,
                log_dir=log_dir,
                samples=get_cellranger_multi_config.output.gex_libraries,
                fastq_dirs=get_cellranger_multi_config.output.fastq_dirs,
                reference_dataset=\
                get_cellranger_multi_config.output.reference_data_path,
                required_tasks=(setup_qc_dirs,))
            verify_qc.requires(run_cellranger_count)

        # Optional additional QC metrics
        if include_extended_metrics and organism:
            self.add_extended_metrics(project,
                                      qc_protocol,
                                      qc_dir,
                                      organism,
                                      log_dir,
                                      pre_tasks=(setup_qc_dirs,),
                                      post_tasks=(verify_qc,))

    def add_cellranger_count(self,project_name,project,qc_dir,
                             organism,fastq_dir,qc_protocol,chemistry,
                             force_cells,log_dir,samples=None,
                             fastq_dirs=None,reference_dataset=None,
                             extra_projects=None,required_tasks=None):
        """
        Add tasks to pipeline to run 'cellranger* count'

        Arguments:
          project_name (str): name to associate with project for
            reporting tasks
          project (AnalysisProject): project to run 10x
            cellranger pipeline within
          qc_dir (str): directory for QC outputs (defaults
            to subdirectory 'qc' of project directory)
          organism (str): organism for pipeline
          fastq_dir (str): directory holding Fastq files
          qc_protocol (str): QC protocol to use
          chemistry (str): chemistry to use in single
            library analysis
          force_cells (int): if set then bypasses
            the cell detection algorithm in 'cellranger'
            and 'cellranger-atac' using the '--force-cells'
            option (does nothing for 'cellranger-arc')
          log_dir (str): directory to write log files to
          samples (list): optional, list of samples to
            restrict single library analyses to (or None
            to use all samples in project)
          fastq_dirs (dict): optional, a dictionary mapping
            sample names to Fastq directories which will
            be used to override the paths set by the
            'fastq_dirs' argument
          reference_dataset (str): optional, path to
            reference dataset (otherwise will be determined
            automatically based on organism)
          extra_projects (list): optional list of extra
            AnalysisProjects to include Fastqs from when
            running cellranger pipeline
          required_tasks (list): list of tasks that the
            cellranger pipeline should wait for
        """
        # Tasks 'check_cellranger_count' depends on
        check_cellranger_count_requires = []

        # Locate cellranger
        required_cellranger = DetermineRequired10xPackage(
            "%s: determine required 10x pipeline package (%s)" %
            (project_name,qc_protocol),
            qc_protocol,
            self.params.cellranger_exe)
        self.add_task(required_cellranger)

        get_cellranger = Get10xPackage(
            "%s: get information on 10x pipeline package (%s)" %
            (project_name,qc_protocol),
            require_package=\
            required_cellranger.output.require_cellranger)
        self.add_task(get_cellranger,
                      requires=(required_cellranger,),
                      envmodules=self.envmodules['cellranger'])
        check_cellranger_count_requires.append(get_cellranger)

        # Get reference data for cellranger
        get_cellranger_reference_data = GetCellrangerReferenceData(
            "%s: get single library analysis reference data (%s)" %
            (project_name,qc_protocol),
            project,
            organism=organism,
            transcriptomes=self.params.cellranger_transcriptomes,
            premrna_references=self.params.cellranger_premrna_references,
            atac_references=self.params.cellranger_atac_references,
            multiome_references=self.params.cellranger_arc_references,
            cellranger_exe=get_cellranger.output.package_exe,
            cellranger_version=get_cellranger.output.package_version,
            qc_protocol=qc_protocol,
            force_reference_data=reference_dataset
        )
        self.add_task(get_cellranger_reference_data,
                      requires=(get_cellranger,),
                      runner=self.runners['verify_runner'],
                      log_dir=log_dir)
        check_cellranger_count_requires.append(get_cellranger_reference_data)

        # Make libraries.csv files (cellranger-arc only)
        if qc_protocol in ("10x_Multiome_ATAC",
                           "10x_Multiome_GEX",):
            make_cellranger_libraries = MakeCellrangerArcCountLibraries(
                "%s: make libraries files for 'cellranger-arc count'" %
                project_name,
                project,
                qc_dir
            )
            self.add_task(make_cellranger_libraries,
                          requires=required_tasks,
                          log_dir=log_dir)
            check_cellranger_count_requires.append(
                make_cellranger_libraries)

        # Check QC outputs for cellranger count
        check_cellranger_count = CheckCellrangerCountOutputs(
            "%s: check for single library analysis outputs (%s)" %
            (project_name,qc_protocol),
            project,
            fastq_dir=fastq_dir,
            samples=samples,
            qc_dir=qc_dir,
            qc_protocol=qc_protocol,
            extra_projects=extra_projects,
            cellranger_version=get_cellranger.output.package_version,
            cellranger_ref_data=\
            get_cellranger_reference_data.output.reference_data_path,
            verbose=self.params.VERBOSE
        )
        self.add_task(check_cellranger_count,
                      requires=check_cellranger_count_requires,
                      runner=self.runners['verify_runner'],
                      log_dir=log_dir)

        # Parent directory for cellranger count outputs
        # Set to project directory unless the 'cellranger_out_dir'
        # parameter is set
        cellranger_out_dir = FunctionParam(
            lambda out_dir,project_dir:
            out_dir if out_dir is not None else project_dir,
            self.params.cellranger_out_dir,
            project.dirn)

        # Run cellranger count
        run_cellranger_count = RunCellrangerCount(
            "%s: run single library analysis (%s)" %
            (project_name,qc_protocol),
            check_cellranger_count.output.samples,
            check_cellranger_count.output.fastq_dir,
            get_cellranger_reference_data.output.reference_data_path,
            cellranger_out_dir,
            qc_dir=qc_dir,
            cellranger_exe=get_cellranger.output.package_exe,
            cellranger_version=get_cellranger.output.package_version,
            chemistry=chemistry,
            fastq_dirs=fastq_dirs,
            force_cells=force_cells,
            cellranger_jobmode=self.params.cellranger_jobmode,
            cellranger_maxjobs=self.params.cellranger_maxjobs,
            cellranger_mempercore=self.params.cellranger_mempercore,
            cellranger_jobinterval=self.params.cellranger_jobinterval,
            cellranger_localcores=self.params.cellranger_localcores,
            cellranger_localmem=self.params.cellranger_localmem,
            qc_protocol=qc_protocol
        )
        self.add_task(run_cellranger_count,
                      requires=(get_cellranger,
                                get_cellranger_reference_data,
                                check_cellranger_count,),
                      runner=self.runners['cellranger_runner'],
                      envmodules=self.envmodules['cellranger'],
                      log_dir=log_dir)

        # Return the 'cellranger count' task
        return run_cellranger_count

    def add_extended_metrics(self,project,qc_protocol,qc_dir,
                             organism,log_dir,pre_tasks,post_tasks):
        """
        Add additional optional QC metrics for a project

        Arguments:
          project (AnalysisProject): project to run additional
            QC metrics for
          qc_protocol (str): QC protocol to use
          qc_dir (str): directory for QC outputs
          organism (str): organism(s) associated with the
            run
          log_dir (str): directory for log files (defaults
            to 'logs' subdirectory of the QC directory
          pre_tasks (list): list of tasks that must complete
            before the extended metrics tasks can run
          post_tasks (list): list of tasks that depend on
            for the extended metrics to complete
        """
        # Read numbers with sequence data (based on protocol)
        # Used to set which Fastqs should be mapped
        read_numbers = get_read_numbers(qc_protocol).seq_data

        # Sanitise organism name
        organism_name = str(organism).\
                        strip().\
                        lower().\
                        replace(' ','_')

        # Indicate if data are paired
        paired = (len(read_numbers) > 1)

        # Set up tasks for BAM file generation
        get_star_index = GetReferenceDataset(
            "%s: get STAR index for '%s'" % (project.name,
                                             organism),
            organism,
            self.params.star_indexes)
        self.add_task(get_star_index)
        get_bam_files = GetBAMFiles(
            "%s: get BAM files" % project.name,
            project.fastqs,
            get_star_index.output.reference_dataset,
            os.path.join(qc_dir,'bam_files',organism_name),
            self.params.fastq_subset,
            self.params.nthreads,
            reads=read_numbers,
            verbose=self.params.VERBOSE)
        self.add_task(get_bam_files,
                      requires=pre_tasks,
                      runner=self.runners['star_runner'],
                      log_dir=log_dir)

        # Get reference gene model for RSeQC
        get_reference_gene_model = GetReferenceDataset(
            "%s: get RSeQC reference gene model for '%s'" % (project.name,
                                                             organism),
            organism,
            self.params.annotation_bed_files)
        self.add_task(get_reference_gene_model,
                      log_dir=log_dir)

        # Run RSeQC infer experiment
        rseqc_infer_experiment = RunRSeQCInferExperiment(
            "%s: infer experiment from BAM files (RSeQC)" % project.name,
            get_bam_files.output.bam_files,
            get_reference_gene_model.output.reference_dataset,
            os.path.join(qc_dir,'rseqc_infer_experiment',organism_name))
        self.add_task(rseqc_infer_experiment,
                      log_dir=log_dir)

        # Run RSeQC gene body coverage
        rseqc_gene_body_coverage = RunRSeQCGenebodyCoverage(
            "%s: calculate gene body coverage (RSeQC)" % project.name,
            get_bam_files.output.bam_files,
            get_reference_gene_model.output.reference_dataset,
            os.path.join(qc_dir,'rseqc_genebody_coverage',organism_name),
            name=project.name)
        self.add_task(rseqc_gene_body_coverage,
                      log_dir=log_dir)
        for task in post_tasks:
            rseqc_gene_body_coverage.required_by(task)

        # Run Picard's CollectInsertSizeMetrics
        insert_size_metrics = RunPicardCollectInsertSizeMetrics(
            "%s: Picard: collect insert size metrics" % project.name,
            get_bam_files.output.bam_files,
            os.path.join(qc_dir,'picard',organism_name),
            bam_properties=rseqc_infer_experiment.output.experiments)
        self.add_task(insert_size_metrics,
                      log_dir=log_dir)

        collate_insert_sizes = CollateInsertSizes(
            "%s: collate insert size data" % project.name,
            get_bam_files.output.bam_files,
            os.path.join(qc_dir,'picard',organism_name),
            os.path.join(qc_dir,'insert_sizes.%s.tsv' % organism_name))
        self.add_task(collate_insert_sizes,
                      requires=(insert_size_metrics,),
                      log_dir=log_dir)
        for task in post_tasks:
            collate_insert_sizes.required_by(task)

        # Get reference gene model for Qualimap
        get_annotation_gtf = GetReferenceDataset(
            "%s: get GTF annotation for '%s'" % (project.name,
                                                 organism),
            organism,
            self.params.annotation_gtf_files)
        self.add_task(get_annotation_gtf,
                      log_dir=log_dir)

        # Run Qualimap RNA-seq analysis
        qualimap_rnaseq = RunQualimapRnaseq(
            "%s: Qualimap rnaseq: RNA-seq metrics" % project.name,
            get_bam_files.output.bam_files,
            get_annotation_gtf.output.reference_dataset,
            os.path.join(qc_dir,'qualimap-rnaseq',organism_name),
            bam_properties=rseqc_infer_experiment.output.experiments)
        self.add_task(qualimap_rnaseq,
                      runner=self.runners['qualimap_runner'],
                      log_dir=log_dir)
        for task in post_tasks:
            qualimap_rnaseq.required_by(task)

    def run(self,nthreads=None,fastq_screens=None,star_indexes=None,
            annotation_bed_files=None,annotation_gtf_files=None,
            fastq_subset=None,cellranger_chemistry='auto',
            cellranger_force_cells=None,cellranger_transcriptomes=None,
            cellranger_premrna_references=None,
            cellranger_atac_references=None,
            cellranger_arc_references=None,cellranger_jobmode='local',
            cellranger_maxjobs=None,cellranger_mempercore=None,
            cellranger_jobinterval=None,cellranger_localcores=None,
            cellranger_localmem=None,cellranger_exe=None,
            cellranger_extra_projects=None,cellranger_reference_dataset=None,
            cellranger_out_dir=None,working_dir=None,log_file=None,
            batch_size=None,batch_limit=None,max_jobs=1,max_slots=None,
            poll_interval=5,runners=None,default_runner=None,
            enable_conda=False,conda=None,conda_env_dir=None,
            envmodules=None,legacy_screens=False,verbose=False):
        """
        Run the tasks in the pipeline

        Arguments:
          nthreads (int): number of threads/processors to
            use for QC jobs (defaults to number of slots set
            in job runners)
          fastq_screens (dict): mapping of screen IDs to
            FastqScreen conf files
          star_indexes (dict): mapping of organism IDs to
            directories with STAR indexes
          annotation_bed_files (dict): mapping of organism
            IDs to BED files with annotation data
          annotation_gtf_files (dict): mapping of organism
            IDs to GTF files with annotation data
          fastq_subset (int): explicitly specify
            the subset size for subsetting running Fastqs
          cellranger_chemistry (str): explicitly specify
            the assay configuration (set to 'auto' to let
            cellranger determine this automatically; ignored
            if not scRNA-seq)
          force_cells (int): explicitly specify number of
            cells for 'cellranger' and 'cellranger-atac'
            (set to 'None' to use the cell detection
            algorithm; ignored for 'cellranger-arc')
          cellranger_transcriptomes (mapping): mapping of
            organism names to reference transcriptome data
            for cellranger
          cellranger_premrna_references (mapping):
            mapping of organism names to "pre-mRNA"
            reference data for cellranger
          cellranger_atac_references (mapping): mapping of
            organism names to ATAC-seq reference genome data
            for cellranger-atac
          cellranger_arc_references (mapping): mapping of
            organism names to multiome reference datasets
            for cellranger-arc
          cellranger_jobmode (str): specify the job mode to
            pass to cellranger (default: "local")
          cellranger_maxjobs (int): specify the maximum
            number of jobs to pass to cellranger (default:
            None)
          cellranger_mempercore (int): specify the memory
            per core (in Gb) to pass to cellranger (default:
            None)
          cellranger_jobinterval (int): specify the interval
            between launching jobs (in ms) to pass to
            cellranger (default: None)
          cellranger_localcores (int): maximum number of cores
            cellranger can request in jobmode 'local'
            (default: None)
          cellranger_localmem (int): maximum memory cellranger
            can request in jobmode 'local' (default: None)
          cellranger_exe (str): optional, explicitly specify
            the cellranger executable to use for single
            library analysis (default: cellranger executable
            is determined automatically)
          cellranger_extra_projects (list): optional list of
            extra AnalysisProjects to include Fastqs from
            when running cellranger pipeline
          cellranger_reference_dataset (str): optional,
            explicitly specify the path to the reference
            dataset to use for single library analysis
            (default: reference dataset is determined
            automatically)
          cellranger_out_dir (str): specify directory to
            put full cellranger outputs into (default:
            project directory)
          working_dir (str): optional path to a working
            directory (defaults to temporary directory in
            the current directory)
          log_dir (str): path of directory where log files
            will be written to
          batch_size (int): if set then run commands in
            each task in batches, with each batch running
            this many commands at a time (default is to run
            one command per job)
          batch_limit (int): if set then run commands in
            each task in batches, with the batch size set
            dyanmically so as not to exceed this limit
            (default is to use fixed batch sizes)
          max_jobs (int): optional maximum number of
            concurrent jobs in scheduler (defaults to 1)
          max_slots (int): optional maximum number of 'slots'
            (i.e. concurrent threads or maximum number of
            CPUs) available to the scheduler (defaults to
            no limit)
          poll_interval (float): optional polling interval
            (seconds) to set in scheduler (defaults to 5s)
          runners (dict): mapping of names to JobRunner
            instances; valid names are 'fastqc_runner',
            'fastq_screen_runner','star_runner',
            'report_runner','cellranger_runner',
            'verify_runner', and 'default'
          enable_conda (bool): if True then enable use of
            conda environments to satisfy task dependencies
          conda (str): path to conda
          conda_env_dir (str): path to non-default
            directory for conda environments
          envmodules (mapping): mapping of names to
            environment module file lists; valid names are
            'fastqc','fastq_screen','fastq_strand','cellranger',
            'report_qc'
          default_runner (JobRunner): optional default
            job runner to use
          legacy_screens (bool): if True then use 'legacy'
            naming convention for FastqScreen outputs
          verbose (bool): if True then report additional
            information for diagnostics
        """
        # Working directory
        clean_up_on_completion = False
        if working_dir is None:
            working_dir = tempfile.mkdtemp(prefix="__qc.",
                                           suffix=".tmp",
                                           dir=os.getcwd())
            clean_up_on_completion = True
        working_dir = os.path.abspath(working_dir)
        if not os.path.exists(working_dir):
            mkdir(working_dir)

        # Log and script directories
        tasks_work_dir = os.path.join(working_dir,"work")
        log_dir = os.path.join(working_dir,"logs")
        scripts_dir = os.path.join(working_dir,"scripts")

        # Runners
        if runners is None:
            runners = dict()

        # Execute the pipeline
        status = Pipeline.run(self,
                              working_dir=working_dir,
                              tasks_work_dir=tasks_work_dir,
                              log_dir=log_dir,
                              scripts_dir=scripts_dir,
                              log_file=log_file,
                              enable_conda=enable_conda,
                              conda=conda,
                              conda_env_dir=conda_env_dir,
                              batch_size=batch_size,
                              batch_limit=batch_limit,
                              exit_on_failure=PipelineFailure.DEFERRED,
                              params={
                                  'nthreads': nthreads,
                                  'annotation_bed_files': annotation_bed_files,
                                  'annotation_gtf_files': annotation_gtf_files,
                                  'fastq_subset': fastq_subset,
                                  'cellranger_chemistry': cellranger_chemistry,
                                  'cellranger_force_cells':
                                  cellranger_force_cells,
                                  'cellranger_transcriptomes':
                                  cellranger_transcriptomes,
                                  'cellranger_premrna_references':
                                  cellranger_premrna_references,
                                  'cellranger_atac_references':
                                  cellranger_atac_references,
                                  'cellranger_arc_references':
                                  cellranger_arc_references,
                                  'cellranger_jobmode': cellranger_jobmode,
                                  'cellranger_maxjobs': cellranger_maxjobs,
                                  'cellranger_mempercore': cellranger_mempercore,
                                  'cellranger_jobinterval': cellranger_jobinterval,
                                  'cellranger_localcores': cellranger_localcores,
                                  'cellranger_localmem': cellranger_localmem,
                                  'cellranger_exe': cellranger_exe,
                                  'cellranger_reference_dataset':
                                  cellranger_reference_dataset,
                                  'cellranger_out_dir': cellranger_out_dir,
                                  'cellranger_extra_projects':
                                  cellranger_extra_projects,
                                  'fastq_screens': fastq_screens,
                                  'star_indexes': star_indexes,
                                  'legacy_screens': legacy_screens,
                              },
                              poll_interval=poll_interval,
                              max_jobs=max_jobs,
                              max_slots=max_slots,
                              runners=runners,
                              default_runner=default_runner,
                              envmodules=envmodules,
                              verbose=verbose)

        # Clean up working dir
        if status == 0 and clean_up_on_completion:
            shutil.rmtree(working_dir)

        # Return pipeline status
        return status

######################################################################
# Pipeline task classes
######################################################################

class SetupQCDirs(PipelineFunctionTask):
    """
    Set up the directories for the QC run
    """
    def init(self,project,qc_dir,log_dir=None,qc_protocol=None):
        """
        Initialise the SetupQCDirs task

        Arguments:
          project (AnalysisProject): project to run
            QC for
          qc_dir (str): directory for QC outputs (defaults
            to subdirectory 'qc' of project directory)
          log_dir (str): directory for log files (defaults
            to 'logs' subdirectory of the QC directory
          qc_protocol (str): QC protocol to use
        """
        pass
    def setup(self):
        # Check the QC protocol
        qc_info = self.args.project.qc_info(self.args.qc_dir)
        stored_protocol = qc_info.protocol
        if stored_protocol is not None and \
           stored_protocol != self.args.qc_protocol:
            logger.warning("QC protocol mismatch for %s: "
                           "'%s' stored, '%s' specified"
                           % (self.args.project.name,
                              stored_protocol,
                              self.args.qc_protocol))
            logger.warning("Stored protocol will be ignored")
        # Set up QC dir
        if not os.path.exists(self.args.qc_dir):
            mkdir(self.args.qc_dir)
        # Set up log dir
        if self.args.log_dir is None:
            log_dir = os.path.join(self.args.qc_dir,'logs')
        else:
            log_dir = self.args.log_dir
        if not os.path.exists(log_dir):
            mkdir(log_dir)

class UpdateQCMetadata(PipelineTask):
    """
    Update the metadata stored for this QC run
    """
    def init(self,project,qc_dir,metadata,legacy_screens=False):
        """
        Initialise the UpdateQCMetadata task

        Arguments:
          project (AnalysisProject): project to run
            QC for
          qc_dir (str): directory for QC outputs (defaults
            to subdirectory 'qc' of project directory)
          metadata (dict): mapping of metadata items to
            values
          legacy_screens (bool): if True then 'legacy'
            naming convention was used for FastqScreen
            outputs
        """
        pass
    def setup(self):
        # Resolve the supplied metadata values
        metadata = {}
        for item in self.args.metadata:
            try:
                value = self.args.metadata[item].value
            except AttributeError:
                value = self.args.metadata[item]
            metadata[item] = value
        # Deal with FastqScreen metadata
        fastq_screens = (metadata['fastq_screens']
                         if 'fastq_screens' in metadata else None)
        if fastq_screens and not self.args.legacy_screens:
            # Collapse list into a string
            metadata['fastq_screens'] = ','.join([s for s
                                                  in fastq_screens])
        else:
            metadata['fastq_screens'] = None
        # Store the QC metadata
        qc_info = self.args.project.qc_info(self.args.qc_dir)
        for item in metadata:
            print("-- %s: %s" % (item,metadata[item]))
            qc_info[item] = metadata[item]
        qc_info.save()

class GetSeqLengthStats(PipelineFunctionTask):
    """
    Get data on sequence lengths, masking and padding
    for Fastqs in a project, and write the data to
    JSON files.
    """
    def init(self,project,qc_dir,qc_protocol=None,fastq_attrs=None):
        """
        Initialise the GetSeqLengthStats task

        Arguments:
          project (AnalysisProject): project with Fastqs
            to get the sequence length data from
          qc_dir (str): directory for QC outputs (defaults
            to subdirectory 'qc' of project directory)
          qc_protocol (str): QC protocol being used
          fastq_attrs (BaseFastqAttrs): class to use for
            extracting data from Fastq names
        """
        self._fastqs = list()
    def setup(self):
        # Remove index Fastqs
        self._fastqs = remove_index_fastqs(
            self.args.project.fastqs,
            fastq_attrs=self.args.fastq_attrs)
        # Get sequence length data for Fastqs
        read_numbers = get_read_numbers(self.args.qc_protocol).qc
        for fastq in self._fastqs:
            if self.args.fastq_attrs(fastq).read_number not in read_numbers:
                continue
            outfile = os.path.join(self.args.qc_dir,
                                   "%s_seqlens.json" %
                                   self.args.fastq_attrs(fastq))
            if os.path.exists(outfile):
                continue
            self.add_call(
                "Get read lengths for %s" % os.path.basename(fastq),
                get_sequence_lengths,
                fastq,
                outfile=outfile)
    def finish(self):
        for result in self.result():
            # Fastq name
            fastq = result['fastq']
            # Check output file exists
            outfile = os.path.join(self.args.qc_dir,
                                   "%s_seqlens.json" %
                                   self.args.fastq_attrs(fastq))
            if not os.path.exists(outfile):
                raise Exception("Missing sequence length file: %s"  %
                                outfile)

class CheckFastqScreenOutputs(PipelineFunctionTask):
    """
    Check the outputs from FastqScreen
    """
    def init(self,project,qc_dir,screens,qc_protocol=None,
             legacy=False,verbose=False):
        """
        Initialise the CheckFastqScreenOutputs task.

        Arguments:
          project (AnalysisProject): project to run
            QC for
          qc_dir (str): directory for QC outputs (defaults
            to subdirectory 'qc' of project directory)
          screens (mapping): mapping of screen names to
            FastqScreen conf files
          qc_protocol (str): QC protocol to use
          legacy (bool): if True then use 'legacy' naming
            convention for output files (default is to
            use new format)
          verbose (bool): if True then print additional
            information from the task

        Outputs:
          fastqs (list): list of Fastqs that have
            missing FastqScreen outputs under the specified
            QC protocol
        """
        self.add_output('fastqs',list())
    def setup(self):
        if self.args.screens is None:
            print("No screens supplied: no missing QC outputs from "
                  "FastqScreen:")
            return
        # Report if legacy naming convention will be used
        if self.args.legacy:
            print("Checking for legacy FastqScreen output names")
        for screen in self.args.screens:
            self.add_call("Check FastqScreen outputs for %s (%s)"
                          % (self.args.project.name,screen),
                          check_fastq_screen_outputs,
                          self.args.project,
                          self.args.qc_dir,
                          screen,
                          self.args.qc_protocol,
                          legacy=self.args.legacy)
    def finish(self):
        fastqs = set()
        for result in self.result():
            for fq in result:
                fastqs.add(fq)
        self.output.fastqs.extend(list(fastqs))
        if self.output.fastqs:
            if self.args.verbose:
                print("Fastqs with missing QC outputs from "
                      "FastqScreen:")
                for fq in self.output.fastqs:
                    print("-- %s" % fq)
            else:
                print("%s Fastqs with missing QC outputs from "
                      "FastqScreen" % len(self.output.fastqs))
        else:
            print("No Fastqs with missing QC outputs from "
                  "FastqScreen")

class RunFastqScreen(PipelineTask):
    """
    Run FastqScreen
    """
    def init(self,fastqs,qc_dir,screens,subset=None,nthreads=None,
             qc_protocol=None,fastq_attrs=None,legacy=False):
        """
        Initialise the RunIlluminaQC task.

        Arguments:
          fastqs (list): list of paths to Fastq files to
            run Fastq Screen on (it is expected that this
            list will come from the CheckIlluminaQCOutputs
            task)
          qc_dir (str): directory for QC outputs (defaults
            to subdirectory 'qc' of project directory)
          screens (mapping): mapping of screen names to
            FastqScreen conf files
          subset (int): explicitly specify the subset size
            for running Fastq_screen
          nthreads (int): number of threads/processors to
            use (defaults to number of slots set in runner)
          qc_protocol (str): QC protocol to use
          fastq_attrs (BaseFastqAttrs): class to use for
            extracting data from Fastq names
          legacy (bool): if True then use 'legacy' naming
            convention for output files (default is to
            use new format)
        """
        self.conda("fastq-screen=0.14.0",
                   "bowtie=1.2.3")
        # Also need to specify tbb=2020.2 for bowtie
        # See https://www.biostars.org/p/494922/
        self.conda("tbb=2020.2")
    def setup(self):
        if not self.args.fastqs:
            print("Nothing to do")
            return
        # Report if legacy naming convention will be used
        if self.args.legacy:
            print("Using legacy FastqScreen output names")
        # Set up the FastqScreen runs for each Fastq
        read_numbers = get_read_numbers(self.args.qc_protocol).seq_data
        for fastq in self.args.fastqs:
            if self.args.fastq_attrs(fastq).read_number not in read_numbers:
                continue
            # Base name for Fastq file
            fastq_basename = os.path.basename(fastq)
            while fastq_basename.split('.')[-1] in ('fastq','gz'):
                fastq_basename = '.'.join(fastq_basename.split('.')[:-1])
            # Run FastqScreen for each screen
            for screen in self.args.screens:
                # Locate the FastqScreen conf file
                fastq_screen_conf = os.path.abspath(self.args.screens[screen])
                if not os.path.isfile(fastq_screen_conf):
                    raise Exception("%s: conf file not found (or is not "
                                    "a file)" % fastq_screen_conf)
                # Determine base name for outputs
                if self.args.legacy:
                    screen_basename = "%s_%s_screen" % (fastq_basename,
                                                        screen)
                else:
                    screen_basename = "%s_screen_%s" % (fastq_basename,
                                                        screen)
                # Set parameters
                if self.args.nthreads:
                    nthreads = self.args.nthreads
                else:
                    nthreads = self.runner_nslots
                if self.args.subset is not None:
                    subset_option = "--subset %d" %  self.args.subset
                else:
                    subset_option = ""
                # Run FastqScreen
                self.add_cmd(
                    "Run FastqScreen '%s' on %s" % (screen,
                                                    os.path.basename(fastq)),
                    """
                    # Make temporary working directory
                    WORKDIR=$(mktemp -d --tmpdir=.)
                    # Run FastqScreen: {screen_name}
                    fastq_screen \\
                    --conf {fastq_screen_conf} \\
                    --threads {nthreads} {subset_option} \\
                    --outdir $WORKDIR --force \\
                    {fastq}
                    # Rename and move outputs to final location
                    out_txt=$WORKDIR/{fastq_basename}_screen.txt
                    if [ -e $out_txt ] ; then
                        /bin/mv $out_txt {qc_dir}/{screen_basename}.txt
                    else
                       echo "ERROR missing $out_txt" >&2
                       exit 1
                    fi
                    out_png=$WORKDIR/{fastq_basename}_screen.png
                    if [ -e $out_png ] ; then
                        /bin/mv $out_png {qc_dir}/{screen_basename}.png
                    else
                       echo "ERROR missing $out_png output" >&2
                       exit 1
                    fi
                    """.format(screen_name=screen,
                               fastq=fastq,
                               qc_dir=self.args.qc_dir,
                               fastq_screen_conf=fastq_screen_conf,
                               nthreads=nthreads,
                               subset_option=subset_option,
                               fastq_basename=fastq_basename,
                               screen_basename=screen_basename))

class CheckFastQCOutputs(PipelineFunctionTask):
    """
    Check the outputs from FastQC
    """
    def init(self,project,qc_dir,qc_protocol,verbose=False):
        """
        Initialise the CheckFastQCOutputs task.

        Arguments:
          project (AnalysisProject): project to run
            QC for
          qc_dir (str): directory for QC outputs (defaults
            to subdirectory 'qc' of project directory)
          qc_protocol (str): QC protocol to use
          verbose (bool): if True then print additional
            information from the task

        Outputs:
          fastqs (list): list of Fastqs that have
            missing FastQC outputs under the specified
            QC protocol
        """
        self.add_output('fastqs',list())
    def setup(self):
        self.add_call("Check FastQC outputs for %s"
                      % self.args.project.name,
                      check_fastqc_outputs,
                      self.args.project,
                      self.args.qc_dir,
                      self.args.qc_protocol)
    def finish(self):
        fastqs = set()
        for result in self.result():
            for fq in result:
                fastqs.add(fq)
        self.output.fastqs.extend(list(fastqs))
        if self.output.fastqs:
            if self.args.verbose:
                print("Fastqs with missing QC outputs from "
                      "FastQC:")
                for fq in self.output.fastqs:
                    print("-- %s" % fq)
            else:
                print("%s Fastqs with missing QC outputs from "
                      "FastQC" % len(self.output.fastqs))
        else:
            print("No Fastqs with missing QC outputs from "
                  "FastQC")

class RunFastQC(PipelineTask):
    """
    Run FastQC
    """
    def init(self,fastqs,qc_dir,nthreads=None):
        """
        Initialise the RunIlluminaQC task.

        Arguments:
          fastqs (list): list of paths to Fastq files to
            run Fastq Screen on (it is expected that this
            list will come from the CheckIlluminaQCOutputs
            task)
          qc_dir (str): directory for QC outputs (defaults
            to subdirectory 'qc' of project directory)
          nthreads (int): number of threads/processors to
            use (defaults to number of slots set in runner)
        """
        self.conda("fastqc=0.11.3")
    def setup(self):
        if not self.args.fastqs:
            print("Nothing to do")
            return
        # Set number of threads
        if self.args.nthreads:
            nthreads = self.args.nthreads
        else:
            nthreads = self.runner_nslots
        # Run FastQC for each Fastq
        for fastq in self.args.fastqs:
            self.add_cmd(
                "Run FastQC on %s" % os.path.basename(fastq),
                """
                # Run FastQC
                fastqc \\
                --outdir {qc_dir} \\
                --threads {nthreads} \\
                --nogroup --extract \\
                {fastq}
                """.format(qc_dir=self.args.qc_dir,
                           nthreads=nthreads,
                           fastq=fastq))

class SetupFastqStrandConf(PipelineFunctionTask):
    """
    Set up a fastq_strand.conf file
    """
    def init(self,project,qc_dir=None,organism=None,star_indexes=None):
        """
        Initialise the SetupFastqStrandConf task.

        Arguments:
          project (AnalysisProject): project to run
            QC for
          qc_dir (str): if supplied then points to directory
            for QC outputs (defaults to subdirectory 'qc'
            of project directory)
          organism (str): if supplied then must be a
            string with the names of one or more organisms,
            with multiple organisms separated by spaces
            (defaults to the organisms associated with
            the project)
          star_indexes (dict): dictionary mapping
            normalised organism names to STAR indexes

        Outputs:
          fastq_strand_conf (PipelineParam): pipeline
            parameter instance that resolves to a string
            with the path to the generated config file.
        """
        self.add_output('fastq_strand_conf',Param(type=str))
    def setup(self):
        qc_dir = self.args.qc_dir
        if qc_dir is None:
            qc_dir = 'qc'
        if not os.path.isabs(qc_dir):
            qc_dir = os.path.join(self.args.project.dirn,qc_dir)
        self.fastq_strand_conf = os.path.join(qc_dir,
                                              "fastq_strand.conf")
        self.add_call("Setup fastq_strand.conf for %s"
                      % self.args.project.name,
                      self.setup_fastq_strand_conf,
                      self.args.project,
                      self.fastq_strand_conf,
                      organism=self.args.organism,
                      star_indexes=self.args.star_indexes)
    def setup_fastq_strand_conf(self,project,fastq_strand_conf,
                                organism=None,star_indexes=None):
        # Remove existing fastq_strand.conf file
        if os.path.exists(fastq_strand_conf):
            print("Removing existing conf file: %s" % fastq_strand_conf)
            os.remove(fastq_strand_conf)
        # Sort out organism(s)
        if not organism:
            organism = project.info.organism
        if not organism:
            print("No organisms specified")
            return
        print("Organism(s): %s" % organism)
        # Fetch matching STAR indexes
        if not star_indexes:
            print("No STAR indexes available")
            return
        print("STAR indexes: %s" % star_indexes)
        star_indexes = build_fastq_strand_conf(
            get_organism_list(organism),
            star_indexes)
        if not star_indexes:
            print("No matching indexes for strandedness determination")
            return
        # Create the conf file
        print("Writing conf file: %s" % fastq_strand_conf)
        with open(fastq_strand_conf,'wt') as fp:
            fp.write("%s\n" % star_indexes)
    def finish(self):
        if os.path.exists(self.fastq_strand_conf):
            self.output.fastq_strand_conf.set(
                self.fastq_strand_conf)

class CheckFastqStrandOutputs(PipelineFunctionTask):
    """
    Check the outputs from the fastq_strand.py utility
    """
    def init(self,project,qc_dir,fastq_strand_conf,
             qc_protocol=None,verbose=False):
        """
        Initialise the CheckFastqStrandOutputs task.

        Arguments:
          project (AnalysisProject): project to run
            QC for
          qc_dir (str): directory for QC outputs (defaults
            to subdirectory 'qc' of project directory)
          fastq_strand_conf (str): path to the fastq_strand
            config file
          qc_protocol (str): QC protocol to use
          verbose (bool): if True then print additional
            information from the task

        Outputs:
          fastq_pairs (list): list of tuples with Fastq
            "pairs" that have missing outputs from
            fastq_strand.py under the specified QC
            protocol. A "pair" may be an (R1,R2) tuple,
            or a single Fastq (e.g. (fq,)).
        """
        self.add_output('fastq_pairs',list())
    def setup(self):
        fastq_strand_conf = self.args.fastq_strand_conf
        if not fastq_strand_conf or not os.path.exists(fastq_strand_conf):
            print("No conf file, nothing to check")
            return
        self.add_call("Check fastq_strand.py outputs for %s"
                      % self.args.project.name,
                      check_fastq_strand_outputs,
                      self.args.project,
                      self.args.qc_dir,
                      fastq_strand_conf,
                      qc_protocol=self.args.qc_protocol)
    def finish(self):
        for result in self.result():
            self.output.fastq_pairs.extend(result)
        if self.output.fastq_pairs:
            if self.args.verbose:
                print("Fastq pairs with missing QC outputs from "
                      "fastq_strand.py:")
                for fq_pair in self.output.fastq_pairs:
                    print("-- %s" % (fq_pair,))
            else:
                print("%s Fastq pairs with missing QC outputs from "
                      "fastq_strand.py" % len(self.output.fastq_pairs))
        else:
            print("No Fastqs with missing QC outputs from "
                  "fastq_strand.py")

class RunFastqStrand(PipelineTask):
    """
    Run the fastq_strand.py utility
    """
    def init(self,fastq_pairs,qc_dir,fastq_strand_conf,
             fastq_strand_subset=None,nthreads=None,
             qc_protocol=None):
        """
        Initialise the RunFastqStrand task.

        Arguments:
          fastq_pairs (list): list of tuples with "pairs"
            of Fastq files to run fastq_strand.py on (it is
            expected that this list will come from the
            CheckFastqStrandOutputs task)
          qc_dir (str): directory for QC outputs (defaults
            to subdirectory 'qc' of project directory)
          fastq_strand_conf (str): path to the fastq_strand
            config file to use
          fastq_strand_subset (int): explicitly specify
            the subset size for running fastq_strand
          nthreads (int): number of threads/processors to
            use (defaults to number of slots set in job
            runner)
          qc_protocol (str): QC protocol to use
        """
        self.conda("star=2.7.7a",
                   "future")
    def setup(self):
        if not self.args.fastq_pairs:
            print("No Fastqs: nothing to do")
            return
        elif not self.args.fastq_strand_conf or \
             not os.path.exists(self.args.fastq_strand_conf):
            print("No conf file: nothing to do")
            return
        for fastq_pair in self.args.fastq_pairs:
            cmd = PipelineCommandWrapper(
                "Run fastq_strand.py for %s" %
                os.path.basename(fastq_pair[0]),
                'fastq_strand.py',
                '--conf',self.args.fastq_strand_conf,
                '--outdir',
                os.path.abspath(self.args.qc_dir))
            if self.args.fastq_strand_subset:
                cmd.add_args('--subset',
                             self.args.fastq_strand_subset)
            if self.args.nthreads:
                cmd.add_args('-n',self.args.nthreads)
            else:
                cmd.add_args('-n',self.runner_nslots)
            cmd.add_args(*fastq_pair)
            # Add the command
            self.add_cmd(cmd)

class DetermineRequired10xPackage(PipelineTask):
    """
    Determine which 10xGenomics software package is required

    By default determines the package name based on the
    supplied QC protocol, but this can be overridden by
    explicitly supplying a required package (which can
    also be a path to an executable).

    The output 'require_cellranger' parameter should be
    supplied to the 'Get10xPackage' task, which will
    do the job of actually locating an executable.
    """
    def init(self,qc_protocol,require_cellranger=None):
        """
        Initialise the DetermineRequired10xPackage task

        Argument:
          qc_protocol (str): QC protocol to use
          require_cellranger (str): optional package name
            or path to an executable; if supplied then
            overrides the automatic package determination

        Outputs:
          require_cellranger (pipelineParam): the 10xGenomics
            software package name or path to use
        """
        self.add_output('require_cellranger',Param(type=str))
    def setup(self):
        protocols = {
            "10x_scRNAseq": "cellranger",
            "10x_snRNAseq": "cellranger",
            "10x_scATAC": "cellranger-atac",
            "10x_Multiome_ATAC": "cellranger-arc",
            "10x_Multiome_GEX": "cellranger-arc",
            "10x_CellPlex": "cellranger",
        }
        require_cellranger = self.args.require_cellranger
        if require_cellranger is None:
            try:
                require_cellranger = protocols[self.args.qc_protocol]
            except KeyError:
                raise Exception("Can't identify 10xGenomics package "
                                "required for protocol '%s'" %
                                self.args.qc_protocol)
        print("Required 10x package: %s" % require_cellranger)
        self.output.require_cellranger.set(require_cellranger)

class GetCellrangerReferenceData(PipelineFunctionTask):
    """
    """
    def init(self,project,organism=None,transcriptomes=None,
             premrna_references=None,atac_references=None,
             multiome_references=None,qc_protocol=None,
             cellranger_exe=None,cellranger_version=None,
             force_reference_data=None):
        """
        Initialise the GetCellrangerReferenceData task

        Arguments:
          project (AnalysisProject): project to run
            QC for
          organism (str): if supplied then must be a
            string with the names of one or more organisms,
            with multiple organisms separated by spaces
            (defaults to the organisms associated with
            the project)
          transcriptomes (mapping): mapping of organism names
            to reference transcriptome data for cellranger
          premrna_references (mapping): mapping of organism
            names to "pre-mRNA" reference data for cellranger
          atac_references (mapping): mapping of organism names
            to reference genome data for cellranger-atac
          multiome_references (mapping): mapping of organism
            names to reference datasets for cellranger-arc
          cellranger_exe (str): the path to the Cellranger
            software package to use (e.g. 'cellranger',
            'cellranger-atac', 'spaceranger')
          cellranger_version (str): the version string for the
            Cellranger package
          qc_protocol (str): QC protocol to use
          force_reference_data (str): if supplied then
            will be used as the reference dataset, instead of
            trying to locate appropriate reference data
            automatically

        Outputs:
          reference_data_path (PipelineParam): pipeline
            parameter instance which resolves to a string
            with the path to the reference data set
            corresponding to the supplied organism.
        """
        self.add_output('reference_data_path',Param(type=str))
    def setup(self):
        # Report inputs
        print("Cellranger : %s" % self.args.cellranger_exe)
        print("Version    : %s" % self.args.cellranger_version)
        print("QC protocol: %s" % self.args.qc_protocol)
        print("Organism(s): %s" % self.args.organism)
        # Check pre-determined reference dataset
        if self.args.force_reference_data:
            reference_data = self.args.force_reference_data
            print("Using pre-determined dataset: %s" %
                  reference_data)
            self.output.reference_data_path.set(reference_data)
            return
        # Check that a cellranger exe was supplied
        if self.args.cellranger_exe is None:
            print("No cellranger package supplied")
            return
        # Set the class of references we're going to use,
        # based on cellranger package and QC protocol (where
        # appropriate)
        cellranger_pkg = os.path.basename(self.args.cellranger_exe)
        if cellranger_pkg == "cellranger":
            references = self.args.transcriptomes
            if self.args.qc_protocol == "10x_snRNAseq" and \
               int(self.args.cellranger_version.split('.')[0]) < 5:
                references = self.args.premrna_references
        elif cellranger_pkg == "cellranger-atac":
            references = self.args.atac_references
        elif cellranger_pkg == "cellranger-arc":
            references = self.args.multiome_references
        else:
            self.fail(message="Don't know which reference "
                      "dataset to use for '%s' and protocol '%s'" %
                      (self.args.cellranger_exe,
                       self.args.qc_protocol))
            return
        if references is None:
            references = {}
        # Get organism name
        organism = self.args.organism
        if organism is None:
            organism = self.args.project.info.organism
        if not organism:
            self.fail(message="No organism specified")
            return
        organisms = get_organism_list(organism)
        if len(organisms) > 1:
            self.fail(message="Can't handle multiple organisms")
            return
        # Look up reference
        try:
            self.output.reference_data_path.set(
                references[organisms[0]])
            print("Reference dataset: %s" % references[organisms[0]])
        except KeyError:
            # No reference data available
            print("No reference data available for '%s'" % organism)
        except Exception as ex:
            # Some other problem
            self.fail(message="Failed to get reference data for "
                      "'%s': %s" % (organism,ex))

class MakeCellrangerArcCountLibraries(PipelineFunctionTask):
    """
    Make 'libraries.csv' files for cellranger-arc count
    """
    def init(self,project,qc_dir):
        """
        Initialise the MakeCellrangerArcCountLibraries task.

        Arguments:
          project (AnalysisProject): project to run
            QC for
          qc_dir (str): top-level QC directory to put
            'libraries.csv' files
        """
        pass
    def setup(self):
        # Check for top-level libraries file
        libraries_file = os.path.join(self.args.project.dirn,
                                      "10x_multiome_libraries.info")
        if not os.path.exists(libraries_file):
            # Nothing to do
            print("No 10x multiome libraries file '%s': nothing to "
                  "do " % libraries_file)
            return
        # Read the file and get sample associations between
        # 'local' and 'remote' datasets
        libraries = MultiomeLibraries(libraries_file)
        # Set up libraries.csv files for each sample for
        # cellranger-arc count
        for local_sample in libraries.local_samples:
            libraries_csv = os.path.join(self.args.qc_dir,
                                         "libraries.%s.csv" % local_sample)
            print("Sample '%s': making %s" % (local_sample,
                                              libraries_csv))
            libraries.write_libraries_csv(
                local_sample,
                self.args.project.fastq_dir,
                self.args.project.info.library_type,
                filen=libraries_csv)
            print("Generated %s" % libraries_csv)

class GetCellrangerMultiConfig(PipelineFunctionTask):
    """
    Locate 'config.csv' file for cellranger multi
    """
    def init(self,project,qc_dir):
        """
        Initialise the GetCellrangerMultiConfig task.

        Arguments:
          project (AnalysisProject): project to run
            QC for
          qc_dir (str): top-level QC directory to put
            'config.csv' files
        """
        self.add_output('config_csv',Param(type=str))
        self.add_output('samples',ListParam())
        self.add_output('gex_libraries',ListParam())
        self.add_output('fastq_dirs',dict())
        self.add_output('reference_data_path',Param(type=str))
    def setup(self):
        # Check for top-level libraries file
        config_file = os.path.join(self.args.project.dirn,
                                   "10x_multi_config.csv")
        if not os.path.exists(config_file):
            # Nothing to do
            print("No 10x multi config file '%s': nothing to do" %
                  config_file)
            return
        # Extract information from config.csv file
        print("Reading config.csv file")
        config_csv = CellrangerMultiConfigCsv(config_file)
        samples = config_csv.sample_names
        gex_libraries = config_csv.gex_libraries
        reference_data_path = config_csv.reference_data_path
        fastq_dirs = config_csv.fastq_dirs
        print("Samples:")
        for sample in samples:
            print("- %s" % sample)
        print("GEX libraries:")
        for library in gex_libraries:
            print("- %s" % library)
        print("Reference dataset: %s" % reference_data_path)
        # Copy config file to QC dir
        print("Copy '%s' into %s" % (config_file,self.args.qc_dir))
        shutil.copy(config_file,self.args.qc_dir)
        # Set outputs
        self.output.config_csv.set(os.path.join(self.args.qc_dir,
                                                os.path.basename(config_file)))
        self.output.samples.extend(samples)
        self.output.gex_libraries.extend(gex_libraries)
        self.output.reference_data_path.set(reference_data_path)
        for sample in fastq_dirs:
            self.output.fastq_dirs[sample] = fastq_dirs[sample]

class CheckCellrangerCountOutputs(PipelineFunctionTask):
    """
    Check the outputs from cellranger(-atac) count
    """
    def init(self,project,fastq_dir=None,samples=None,qc_dir=None,
             qc_protocol=None,extra_projects=None,
             cellranger_version=None,cellranger_ref_data=None,
             verbose=False):
        """
        Initialise the CheckCellrangerCountOutputs task.

        Arguments:
          project (AnalysisProject): project to run
            QC for
          fastq_dir (str): directory holding Fastq files
          samples (list): list of samples to restrict
            checks to (all samples in project are checked
            by default)
          qc_dir (str): top-level QC directory to look
            for 'count' QC outputs (e.g. metrics CSV and
            summary HTML files)
          qc_protocol (str): QC protocol to use
          extra_projects (list): optional list of extra
            AnalysisProjects to include Fastqs from when
            running cellranger pipeline
          cellranger_version (str): version number of
            10xGenomics package
          cellranger_ref_data (str): name or path to
            reference dataset for single library analysis
          verbose (bool): if True then print additional
            information from the task

        Outputs:
          fastq_dir (PipelineParam): pipeline parameter
            instance that resolves to a string with the
            path to directory with Fastq files
          samples (list): list of sample names that have
            missing outputs from 'cellranger count'
        """
        self.add_output('fastq_dir',Param(type=str))
        self.add_output('samples',list())
    def setup(self):
        if not self.args.cellranger_ref_data:
            # No reference data, nothing to check
            return
        # Determine which checking function to use
        if self.args.qc_protocol in ("10x_scRNAseq",
                                     "10x_snRNAseq",):
            check_outputs = check_cellranger_count_outputs
        elif self.args.qc_protocol == "10x_scATAC":
            check_outputs = check_cellranger_atac_count_outputs
        elif self.args.qc_protocol in ("10x_Multiome_ATAC",
                                       "10x_Multiome_GEX",):
            check_outputs = check_cellranger_arc_count_outputs
        # Set the prefix for cellranger/10x outputs
        prefix = os.path.join("cellranger_count",
                              self.args.cellranger_version,
                              os.path.basename(self.args.cellranger_ref_data))
        # Check if the outputs exist
        projects = [self.args.project]
        if self.args.extra_projects:
            projects.extend(self.args.extra_projects)
        for project in projects:
            self.add_call("Check cellranger count outputs from %s"
                          % self.args.project.name,
                          check_outputs,
                          project,
                          self.args.qc_dir,
                          prefix=prefix)
    def finish(self):
        if not self.args.cellranger_ref_data:
            # No reference data, nothing to check
            print("No reference data to check against")
            return
        # Collect the sample names with missing outputs
        samples = set()
        for result in self.result():
            for smpl in result:
                if self.args.samples is None or \
                   smpl in self.args.samples:
                    samples.add(smpl)
        self.output.samples.extend(sorted(list(samples)))
        if self.output.samples:
            if self.args.verbose:
                print("Samples with missing outputs from "
                      "cellranger count:")
                for sample in self.output.samples:
                    print("-- %s" % sample)
            else:
                print("%s samples with missing outputs from "
                      "cellranger count" %
                      len(self.output.samples))
        else:
            print("No samples with missing outputs from "
                  "cellranger count")
        # Set the fastq_dir that these are found in
        fastq_dir = self.args.project.fastq_dir
        if self.args.extra_projects:
            for project in self.args.extra_projects:
                fastq_dir += ",%s" % project.fastq_dir
        self.output.fastq_dir.set(fastq_dir)

class RunCellrangerCount(PipelineTask):
    """
    Run 'cellranger count'
    """
    def init(self,samples,fastq_dir,reference_data_path,out_dir,
             qc_dir=None,cellranger_exe=None,cellranger_version=None,
             chemistry='auto',fastq_dirs=None,force_cells=None,
             cellranger_jobmode='local',cellranger_maxjobs=None,
             cellranger_mempercore=None,cellranger_jobinterval=None,
             cellranger_localcores=None,cellranger_localmem=None,
             qc_protocol=None):
        """
        Initialise the RunCellrangerCount task.

        Arguments:
          samples (list): list of sample names to run
            cellranger count on (it is expected that this
            list will come from the
            CheckCellrangerCountsOutputs task)
          fastq_dir (str): path to directory holding the
            Fastq files
          reference_data_path (str): path to the cellranger
            compatible reference dataset
          out_dir (str): top-level directory to copy all
            final 'count' outputs into. Outputs won't be
            copied if no value is supplied
          qc_dir (str): top-level QC directory to put
            'count' QC outputs (e.g. metrics CSV and summary
            HTML files) into. Outputs won't be copied if
            no value is supplied
          cellranger_exe (str): the path to the Cellranger
            software package to use (e.g. 'cellranger',
            'cellranger-atac', 'spaceranger')
          cellranger_version (str): the version string for
            the Cellranger package
          fastq_dirs (dict): optional, a dictionary mapping
            sample names to Fastq directories which will
            be used to override the paths set by the
            'fastq_dirs' argument
          force_cells (int): optional, if set then bypasses
            the cell detection algorithm in 'cellranger'
            and 'cellranger-atac' using the '--force-cells'
            option (does nothing for 'cellranger-arc')
          chemistry (str): assay configuration (set to
            'auto' to let cellranger determine this
            automatically; ignored if not scRNA-seq)
          cellranger_jobmode (str): specify the job mode to
            pass to cellranger (default: "local")
          cellranger_maxjobs (int): specify the maximum
            number of jobs to pass to cellranger (default:
            None)
          cellranger_mempercore (int): specify the memory
            per core (in Gb) to pass to cellranger (default:
            None)
          cellranger_jobinterval (int): specify the interval
            between launching jobs (in ms) to pass to
            cellranger (default: None)
          cellranger_localcores (int): maximum number of cores
            cellranger can request in jobmode 'local'
            (defaults to number of slots set in runner)
          cellranger_localmem (int): maximum memory cellranger
            can request in jobmode 'local' (default: None)
          qc_protocol (str): QC protocol to use
        """
        # Add outputs
        self.add_output('cellranger_version',Param(type=str))
        self.add_output('cellranger_refdata',Param(type=str))
        self.add_output('cellranger_exe',Param(type=str))
    def setup(self):
        # Check if there's anything to do
        if not self.args.samples:
            print("No samples: nothing to do")
            return
        if not self.args.reference_data_path:
            print("No reference data available: nothing to do")
            return
        if not self.args.cellranger_exe:
            raise Exception("No cellranger executable provided")
        if not (self.args.out_dir or self.args.qc_dir):
            raise Exception("Need to provide at least one of "
                            "output directory and QC directory")
        # Cellranger details
        cellranger_exe = self.args.cellranger_exe
        cellranger_package = os.path.basename(cellranger_exe)
        cellranger_version = self.args.cellranger_version
        cellranger_major_version = int(cellranger_version.split('.')[0])
        # Expected outputs from cellranger
        self._top_level_files = ("_cmdline",)
        if cellranger_package in ("cellranger",):
            self._outs_files = ("web_summary.html",
                                "metrics_summary.csv")
        elif cellranger_package in ("cellranger-atac",
                                    "cellranger-arc",):
            self._outs_files = ("web_summary.html",
                                "summary.csv")
        # Check output for each sample, and run cellranger if required
        for sample in self.args.samples:
            # Check final outputs
            run_cellranger_count = False
            counts_dir = os.path.abspath(
                os.path.join(self.args.out_dir,
                             "cellranger_count",
                             cellranger_version,
                             os.path.basename(self.args.reference_data_path),
                             sample))
            outs_dir = os.path.join(counts_dir,"outs")
            for f in self._outs_files:
                path = os.path.join(outs_dir,f)
                if not os.path.exists(path):
                    run_cellranger_count = True
                    break
            for f in self._top_level_files:
                path = os.path.join(counts_dir,f)
                if not os.path.exists(path):
                    run_cellranger_count = True
                    break
            if not run_cellranger_count:
                print("Sample '%s': found existing outputs" % sample)
                continue
            # Create a working directory for this sample
            work_dir = "tmp.count.%s" % sample
            # Get the path(s) to the Fastq directory(ies)
            if self.args.fastq_dirs and sample in self.args.fastq_dirs:
                fastq_dir = self.args.fastq_dirs[sample]
            else:
                fastq_dir = self.args.fastq_dir
            # Build cellranger command
            cmd = Command(cellranger_exe,
                          "count",
                          "--id",sample)
            # Add package-specific options
            if cellranger_package == "cellranger":
                # Cellranger (gene expression data)
                cmd.add_args("--fastqs",self.args.fastq_dir,
                             "--sample",sample,
                             "--transcriptome",
                             self.args.reference_data_path,
                             "--chemistry",self.args.chemistry)
                # Force cells
                if self.args.force_cells:
                    cmd.add_args("--force-cells",
                                 self.args.force_cells)
                # Additional options for cellranger 5.0+
                if cellranger_major_version >= 5:
                    # Hard-trim the input R1 sequence to 26bp
                    cmd.add_args("--r1-length=26")
                    if self.args.qc_protocol == "10x_snRNAseq":
                        # For single nuclei RNA-seq specify the
                        # --include-introns for cellranger 5.0+
                        if cellranger_major_version == 7:
                            cmd.add_args("--include-introns",
                                         "true")
                        else:
                            cmd.add_args("--include-introns")
            elif cellranger_package == "cellranger-atac":
                # Cellranger-ATAC
                cmd.add_args("--fastqs",fastq_dir,
                             "--sample",sample,
                             "--reference",
                             self.args.reference_data_path)
                # Force cells
                if self.args.force_cells:
                    cmd.add_args("--force-cells",
                                 self.args.force_cells)
                # Additional options for cellranger-atac 2+
                if cellranger_major_version >= 2:
                    # Enable chemistry to be specified
                    if self.args.chemistry == "auto":
                        # 'auto' not recognised by cellranger-atac 2.0.0
                        print("Dropping \"--chemistry='auto'\" from "
                              "cellranger-atac 2+ command line")
                    else:
                        cmd.add_args("--chemistry",self.args.chemistry)
            elif cellranger_package == "cellranger-arc":
                # Cellranger-ARC (multiome GEX + ATAC data)
                cmd.add_args("--reference",
                             self.args.reference_data_path,
                             "--libraries",
                             os.path.join(self.args.qc_dir,
                                          "libraries.%s.csv"
                                          % sample))
            else:
                # Unimplemented package
                raise Exception("Don't know how to run 'count' "
                                "for %s" % cellranger_package)
            # Set number of local cores
            if self.args.cellranger_localcores:
                localcores = self.args.cellranger_localcores
            elif self.args.cellranger_jobmode == "local":
                # Get number of local cores from runner
                localcores = self.runner_nslots
            else:
                # Not in jobmode 'local'
                localcores = None
            add_cellranger_args(cmd,
                                jobmode=self.args.cellranger_jobmode,
                                mempercore=self.args.cellranger_mempercore,
                                maxjobs=self.args.cellranger_maxjobs,
                                jobinterval=self.args.cellranger_jobinterval,
                                localcores=localcores,
                                localmem=self.args.cellranger_localmem)
            # Append to command in task
            print("Running %s" % cmd)
            self.add_cmd(PipelineCommandWrapper(
                "Run %s count for %s" % (cellranger_exe,sample),
                "mkdir","-p",work_dir,
                "&&",
                "cd",work_dir,
                "&&",
                *cmd.command_line))
    def finish(self):
        # If no reference data then ignore and return
        if not self.args.reference_data_path:
            print("No reference data: single library analysis was "
                  "skipped")
            return
        # Set outputs
        self.output.cellranger_exe.set(self.args.cellranger_exe)
        self.output.cellranger_refdata.set(self.args.reference_data_path)
        self.output.cellranger_version.set(self.args.cellranger_version)
        # Handle outputs from cellranger count
        has_errors = False
        for sample in self.args.samples:
            # Check outputs
            top_dir = os.path.join("tmp.count.%s" % sample,
                                   sample)
            print("Sample: %s" % sample)
            if not os.path.exists(top_dir):
                # Cellranger count wasn't run for this sample
                print("'cellranger count' not run for this sample?")
                continue
            outs_dir = os.path.join(top_dir,"outs")
            missing_files = []
            for f in self._outs_files:
                path = os.path.join(outs_dir,f)
                if not os.path.exists(path):
                    print("Missing: %s" % path)
                    missing_files.append(path)
            for f in self._top_level_files:
                path = os.path.join(top_dir,f)
                if not os.path.exists(path):
                    print("Missing: %s" % path)
                    missing_files.append(path)
            if missing_files:
                # Skip this sample
                print("Some files missing for this sample, skipping")
                has_errors = True
            else:
                # Move count outputs to final destination
                count_dir = os.path.abspath(
                    os.path.join(self.args.out_dir,
                                 "cellranger_count",
                                 self.args.cellranger_version,
                                 os.path.basename(
                                     self.args.reference_data_path)
                                 ))
                print("Moving %s to %s" % (top_dir,count_dir))
                if not os.path.exists(count_dir):
                    mkdirs(count_dir)
                shutil.move(top_dir,count_dir)
        # Also copy outputs to QC directory
        if self.args.qc_dir and self.args.samples:
            print("Copying outputs to QC directory")
            # Top level output directory
            count_dir = os.path.abspath(
                os.path.join(self.args.out_dir,
                             "cellranger_count",
                             self.args.cellranger_version,
                             os.path.basename(
                                 self.args.reference_data_path)))
            for sample in self.args.samples:
                print("Sample: %s" % sample)
                # Location of outputs
                top_dir = os.path.join(count_dir,sample)
                outs_dir = os.path.join(top_dir,"outs")
                # Set location to copy QC outputs to
                qc_dir = os.path.abspath(
                    os.path.join(self.args.qc_dir,
                                 "cellranger_count",
                                 self.args.cellranger_version,
                                 os.path.basename(
                                     self.args.reference_data_path),
                                 sample))
                qc_outs_dir = os.path.join(qc_dir,"outs")
                # Make directories and copy the files
                mkdirs(qc_outs_dir)
                for f in self._outs_files:
                    path = os.path.join(outs_dir,f)
                    print("Copying %s from %s to %s" % (f,
                                                        outs_dir,
                                                        qc_outs_dir))
                    shutil.copy(path,qc_outs_dir)
                for f in self._top_level_files:
                    path = os.path.join(top_dir,f)
                    print("Copying %s from %s to %s" % (f,
                                                        top_dir,
                                                        qc_dir))
                    shutil.copy(path,qc_dir)
        # Delayed task failure from earlier errors
        if has_errors:
            self.fail(message="Some outputs missing from cellranger "
                      "count")
            return

class RunCellrangerMulti(PipelineTask):
    """
    Run 'cellranger multi'
    """
    def init(self,project,config_csv,samples,reference_data_path,out_dir,
             qc_dir=None,cellranger_exe=None,cellranger_version=None,
             cellranger_jobmode='local',cellranger_maxjobs=None,
             cellranger_mempercore=None,cellranger_jobinterval=None,
             cellranger_localcores=None,cellranger_localmem=None,
             qc_protocol=None,working_dir=None):
        """
        Initialise the RunCellrangerMulti task.

        Arguments:
          project (AnalysisProject): project to run
            QC for
          config_csv (str): path to 'cellranger multi'
            configuration file
          samples (list): list of sample names from the
            config.csv file
          reference_data_path (str): path to the cellranger
            compatible reference dataset from the config.csv
            file
          out_dir (str): top-level directory to copy all
            final 'count' outputs into. Outputs won't be
            copied if no value is supplied
          qc_dir (str): top-level QC directory to put
            'count' QC outputs (e.g. metrics CSV and summary
            HTML files) into. Outputs won't be copied if
            no value is supplied
          cellranger_exe (str): the path to the Cellranger
            software package to use (e.g. 'cellranger',
            'cellranger-atac', 'spaceranger')
          cellranger_version (str): the version string for
            the Cellranger package
          cellranger_jobmode (str): specify the job mode to
            pass to cellranger (default: "local")
          cellranger_maxjobs (int): specify the maximum
            number of jobs to pass to cellranger (default:
            None)
          cellranger_mempercore (int): specify the memory
            per core (in Gb) to pass to cellranger (default:
            None)
          cellranger_jobinterval (int): specify the interval
            between launching jobs (in ms) to pass to
            cellranger (default: None)
          cellranger_localcores (int): maximum number of cores
            cellranger can request in jobmode 'local'
            (defaults to number of slots set in runner)
          cellranger_localmem (int): maximum memory cellranger
            can request in jobmode 'local' (default: None)
          qc_protocol (str): QC protocol to use
        """
        # Internal: top-level working directory
        self._working_dir = None
        # Samples from config.csv
        self._samples = []
        # Whether to run cellranger multi
        self.run_cellranger_multi = False
    def setup(self):
        # Check if there's anything to do
        if not self.args.config_csv:
            print("No config file: nothing to do")
            return
        if not self.args.cellranger_exe:
            raise Exception("No cellranger executable provided")
        if not (self.args.out_dir or self.args.qc_dir):
            raise Exception("Need to provide at least one of "
                            "output directory and QC directory")
        # Top-level working directory
        self._working_dir = self.args.working_dir
        # Cellranger details
        cellranger_exe = self.args.cellranger_exe
        cellranger_package = os.path.basename(cellranger_exe)
        cellranger_version = self.args.cellranger_version
        cellranger_major_version = int(cellranger_version.split('.')[0])
        # Expected outputs from cellranger multi
        self._top_level_files = ("_cmdline",)
        self._outs_files_multi_analysis = ("tag_calls_summary.csv",)
        self._outs_files_per_sample = ("web_summary.html",
                                       "metrics_summary.csv")
        # Check outputs and run cellranger if required
        multi_dir = os.path.abspath(
            os.path.join(self.args.out_dir,
                         "cellranger_multi",
                         cellranger_version,
                         os.path.basename(self.args.reference_data_path)))
        outs_dir = os.path.join(multi_dir,"outs")
        # Per sample outputs
        for sample in self.args.samples:
            for f in self._outs_files_per_sample:
                path = os.path.join(outs_dir,
                                    "per_sample_outs",
                                    sample,f)
                if not os.path.exists(path):
                    self.run_cellranger_multi = True
                    break
        # Multiplexing analysis outputs
        for f in self._outs_files_multi_analysis:
            path = os.path.join(outs_dir,
                                "multi",
                                "multiplexing_analysis",f)
            if not os.path.exists(path):
                self.run_cellranger_multi = True
                break
        # Top level outputs
        for f in self._top_level_files:
            path = os.path.join(multi_dir,f)
            if not os.path.exists(path):
                self.run_cellranger_multi = True
                break
        if not self.run_cellranger_multi:
            print("Found existing outputs")
            return
        # Create a working directory for this sample
        work_dir = os.path.join(self._working_dir,
                                "tmp.cellranger_multi.%s" %
                                self.args.project.name)
        # Build cellranger command
        cmd = Command(cellranger_exe,
                      "multi",
                      "--id",self.args.project.name,
                      "--csv",self.args.config_csv)
        # Set number of local cores
        if self.args.cellranger_localcores:
            localcores = self.args.cellranger_localcores
        elif self.args.cellranger_jobmode == "local":
            # Get number of local cores from runner
            localcores = self.runner_nslots
        else:
            # Not in jobmode 'local'
            localcores = None
        add_cellranger_args(cmd,
                            jobmode=self.args.cellranger_jobmode,
                            mempercore=self.args.cellranger_mempercore,
                            maxjobs=self.args.cellranger_maxjobs,
                            jobinterval=self.args.cellranger_jobinterval,
                            localcores=localcores,
                            localmem=self.args.cellranger_localmem)
        # Add command to task
        print("Running %s" % cmd)
        self.add_cmd(PipelineCommandWrapper(
            "Run %s multi" % cellranger_exe,
            "mkdir","-p",work_dir,
            "&&",
            "cd",work_dir,
            "&&",
            *cmd.command_line))
    def finish(self):
        # If no config.csv then ignore and return
        if not self.args.config_csv:
            print("No config file: cell multiplexing analysis was skipped")
            return
        # Handle outputs from cellranger multi
        has_errors = False
        # Check outputs
        if self.run_cellranger_multi:
            top_dir = os.path.join(self._working_dir,
                                   "tmp.cellranger_multi.%s" %
                                   self.args.project.name,
                                   self.args.project.name)
        else:
            top_dir = os.path.abspath(
                os.path.join(self.args.out_dir,
                             "cellranger_multi",
                             self.args.cellranger_version,
                             os.path.basename(
                                 self.args.reference_data_path)
                ))
        outs_dir = os.path.join(top_dir,"outs")
        missing_files = []
        # Per sample outputs
        for sample in self.args.samples:
            for f in self._outs_files_per_sample:
                path = os.path.join(outs_dir,
                                    "per_sample_outs",
                                    sample,f)
                if not os.path.exists(path):
                    print("Missing: %s" % path)
                    missing_files.append(path)
        # Multiplexing analysis outputs
        for f in self._outs_files_multi_analysis:
            path = os.path.join(outs_dir,
                                "multi",
                                "multiplexing_analysis",f)
            if not os.path.exists(path):
                print("Missing: %s" % path)
                missing_files.append(path)
        # Top level outputs
        for f in self._top_level_files:
            path = os.path.join(top_dir,f)
            if not os.path.exists(path):
                print("Missing: %s" % path)
                missing_files.append(path)
        if missing_files:
            print("Some output files missing from multiplexing analysis")
            has_errors = True
        elif self.run_cellranger_multi:
            # Move multi outputs to final destination
            multi_dir = os.path.abspath(
                os.path.join(self.args.out_dir,
                             "cellranger_multi",
                             self.args.cellranger_version,
                             os.path.basename(
                                 self.args.reference_data_path)
                ))
            print("Moving contents of %s to %s" % (top_dir,multi_dir))
            if not os.path.exists(multi_dir):
                mkdirs(multi_dir)
            for d in os.listdir(top_dir):
                shutil.move(os.path.join(top_dir,d),
                            multi_dir)
        # Also copy outputs to QC directory
        if self.args.qc_dir:
            print("Copying outputs to QC directory")
            # Top level output directory
            top_dir = os.path.abspath(
                os.path.join(self.args.out_dir,
                             "cellranger_multi",
                             self.args.cellranger_version,
                             os.path.basename(
                                 self.args.reference_data_path)))
            # Per sample outputs
            for sample in self.args.samples:
                print("Sample: %s" % sample)
                # Location of outputs
                outs_dir = os.path.join(top_dir,
                                        "outs",
                                        "per_sample_outs",
                                        sample)
                # Set location to copy QC outputs to
                qc_outs_dir = os.path.abspath(
                    os.path.join(self.args.qc_dir,
                                 "cellranger_multi",
                                 self.args.cellranger_version,
                                 os.path.basename(
                                     self.args.reference_data_path),
                                 "outs",
                                 "per_sample_outs",
                                 sample))
                # Make directories and copy the files
                mkdirs(qc_outs_dir)
                for f in self._outs_files_per_sample:
                    path = os.path.join(outs_dir,f)
                    print("Copying %s from %s to %s" % (f,
                                                        outs_dir,
                                                        qc_outs_dir))
                    shutil.copy(path,qc_outs_dir)
            # Multiplexing analysis outputs
            outs_dir = os.path.join(top_dir,
                                    "outs",
                                    "multi",
                                    "multiplexing_analysis")
            qc_outs_dir = os.path.abspath(
                    os.path.join(self.args.qc_dir,
                                 "cellranger_multi",
                                 self.args.cellranger_version,
                                 os.path.basename(
                                     self.args.reference_data_path),
                                 "outs",
                                 "multi",
                                 "multiplexing_analysis"))
            mkdirs(qc_outs_dir)
            for f in self._outs_files_multi_analysis:
                path = os.path.join(outs_dir,f)
                print("Copying %s from %s to %s" % (f,
                                                    outs_dir,
                                                    qc_outs_dir))
                shutil.copy(path,qc_outs_dir)
            # Top level outputs
            qc_top_dir = os.path.abspath(
                    os.path.join(self.args.qc_dir,
                                 "cellranger_multi",
                                 self.args.cellranger_version,
                                 os.path.basename(
                                     self.args.reference_data_path)))
            mkdirs(qc_top_dir)
            for f in self._top_level_files:
                path = os.path.join(top_dir,f)
                print("Copying %s from %s to %s" % (f,
                                                    top_dir,
                                                    qc_top_dir))
                shutil.copy(path,qc_top_dir)
        # Delayed task failure from earlier errors
        if has_errors:
            self.fail(message="Some outputs missing from cellranger multi")
            return

class SetCellCountFromCellrangerCount(PipelineTask):
    """
    Update the number of cells in the project metadata from
    'cellranger count' output
    """
    def init(self,project,qc_dir=None):
        """
        Initialise the SetCellCountFromCellrangerCount task.

        Arguments:
          project (AnalysisProject): project to update the
            number of cells for
          qc_dir (str): directory for QC outputs (defaults
            to subdirectory 'qc' of project directory)
        """
        pass
    def setup(self):
        # Extract and store the cell count from the cellranger
        # metric file
        try:
            set_cell_count_for_project(self.args.project.dirn,
                                       self.args.qc_dir)
        except Exception as ex:
            print("Failed to set the cell count: %s" % ex)

class GetReferenceDataset(PipelineTask):
    """
    Acquire reference data for an organism from mapping

    Generic lookup task which attempts to locate the matching
    reference dataset from a mapping/dictionary.
    """
    def init(self,organism,references):
        """
        Initialise the GetReferenceDataset task

        Arguments:
          organism (str): name of the organism
          references (mapping): mapping with organism names
            as keys and reference datasets as corresponding
            values

        Outputs:
          reference_dataset: reference dataset (set to None
            if no dataset could be located)
        """
        self.add_output('reference_dataset',Param(type='str'))
    def setup(self):
        organism = str(self.args.organism).lower()
        try:
            self.output.reference_dataset.set(self.args.references[organism])
            print("%s: located %s" % (organism,
                                      self.output.reference_dataset.value))
        except Exception as ex:
            print("Unable to locate reference data for organism '%s'"
                  % organism)

class GetBAMFiles(PipelineFunctionTask):
    """
    Create BAM files from Fastqs using STAR

    Runs STAR to generate BAM files from Fastq files. The
    BAMs are then sorted and indexed using samtools.
    """
    def init(self,fastqs,star_index,out_dir,subset_size=None,
             nthreads=None,reads=None,fastq_attrs=None,
             verbose=False):
        """
        Initialise the GetBamFiles task

        Arguments:
          fastqs (list): list of Fastq files to generate
            BAM files from
          star_index (str): path to STAR index to use
          out_dir (str): path to directory to write final
            BAM files to
          subset_size (int): specify size of a random subset
            of reads to use in BAM file generation
          nthreads (int): number of cores for STAR
            to use
          reads (list): optional, list of read numbers to
            include (e.g. [1,2], [2] etc)
          fastq_attrs (IlluminaFastq): optional, class to
            use for extracting information from Fastq file
            names
          verbose (bool): if True then print additional
            information from the task

        Outputs:
          bam_files: list of sorted BAM files
        """
        # Conda dependencies
        self.conda("star=2.4.2a",
                   "samtools")
        # Internal variables
        self._bam_files = list()
        # Outputs
        self.add_output('bam_files',ListParam())
    def setup(self):
        # Check for STAR index
        if self.args.star_index:
            print("STAR index: %s" % self.args.star_index)
        else:
            self.report("STAR index is not set, cannot generate BAM files")
            return
        # Remove index reads and group Fastqs
        fq_pairs = group_fastqs_by_name(
            remove_index_fastqs(self.args.fastqs))
        # Filter reads
        if self.args.fastq_attrs:
            fastq_attrs = self.args.fastq_attrs
        else:
            fastq_attrs = AnalysisFastq
        if self.args.reads:
            filtered_pairs = []
            for fq_pair in fq_pairs:
                filtered_pair = []
                for fq in fq_pair:
                    if fastq_attrs(fq).read_number in self.args.reads:
                        filtered_pair.append(fq)
                filtered_pairs.append(filtered_pair)
            fq_pairs = filtered_pairs
        # Deal with threads
        if self.args.nthreads:
            nthreads = self.args.nthreads
        else:
            nthreads = self.runner_nslots
        # Check for output directory
        if not os.path.exists(self.args.out_dir):
            print("Making output dir: %s" % self.args.out_dir)
            os.makedirs(self.args.out_dir)
        # Check each Fastq pair to see if a corresponding
        # BAM file already exists
        for fq_pair in fq_pairs:
            if self.args.verbose:
                print("-- Fastq pair: %s" % fq_pair)
            bam_file = os.path.join(self.args.out_dir,
                                    "%s.bam" % get_bam_basename(
                                        fq_pair[0],
                                        self.args.fastq_attrs))
            if self.args.verbose:
                print("   BAM file: %s" % bam_file)
            # Add to the list of expected outputs
            self._bam_files.append(bam_file)
            if not os.path.exists(bam_file) and self.args.star_index:
                # Generate BAM file
                self.add_call("Make BAM file",
                              self.make_bam_file,
                              fq_pair,
                              self.args.star_index,
                              bam_file,
                              size=self.args.subset_size,
                              nthreads=nthreads)
    def make_bam_file(self,fastqs,genomedir,bam_file,size=None,
                      nthreads=None):
        ##############################
        # Generate and sort a BAM file
        ##############################
        # Basename and prefix for output files
        basename = os.path.basename(bam_file[:-4])
        prefix = "%s_" % basename
        # Make a temporary directory for fastqs
        if size:
            fastqs_dir = os.path.abspath("__%s_subset%d" % (basename,size))
        else:
            fastqs_dir = os.path.abspath("__%s" % basename)
        print("Creating directory for Fastq subsetting: %s" % fastqs_dir)
        mkdir(fastqs_dir)
        # Generate subset of input reads
        nreads = sum(1 for i in getreads(os.path.abspath(fastqs[0])))
        if not size:
            print("Using all reads/read pairs in Fastq file(s)")
            size = nreads
        elif size > nreads:
            print("Number of reads smaller than requested subset")
            size = nreads
        else:
            print("Using random subset of %d reads/read pairs"
                  % size)
        # Generate subset indices to extract
        if size == nreads:
            subset_indices = [i for i in range(nreads)]
        else:
            subset_indices = random.sample(range(nreads),size)
        # Do the subsetting
        fqs_in = filter(lambda fq: fq is not None,fastqs)
        fastqs = []
        for fq in fqs_in:
            fq_subset = os.path.join(fastqs_dir,os.path.basename(fq))
            if fq_subset.endswith(".gz"):
                fq_subset = '.'.join(fq_subset.split('.')[:-1])
            with open(fq_subset,'w') as fp:
                for read in getreads_subset(os.path.abspath(fq),
                                            subset_indices):
                    fp.write('\n'.join(read) + '\n')
                fastqs.append(fq_subset)
        # Make a temporary directory for STAR
        star_dir = os.path.abspath("__%s_STAR" % basename)
        print("Creating directory for STAR subsetting: %s" % star_dir)
        mkdir(star_dir)
        # Output BAM file name will be "<prefix>Aligned.out.bam"
        star_bam_file = os.path.join(star_dir,
                                     "%sAligned.out.bam" % prefix)
        # Build the STAR command line for mapping
        star_cmd = Command('STAR',
                           '--runMode','alignReads',
                           '--genomeLoad','NoSharedMemory',
                           '--genomeDir',os.path.abspath(genomedir),
                           '--readFilesIn',fastqs[0])
        if len(fastqs) > 1:
            star_cmd.add_args(fastqs[1])
        star_cmd.add_args('--outSAMtype','BAM','Unsorted',
                          '--outSAMstrandField','intronMotif',
                          '--outFileNamePrefix',prefix,
                          '--runThreadN',nthreads)
        print("Running %s" % star_cmd)
        status = star_cmd.run_subprocess(working_dir=star_dir)
        if status != 0:
            raise Exception("STAR returned non-zero exit code: %s"
                            % status)
        # Make a temporary directory for sorting the BAM file
        sort_dir = os.path.abspath("__%s_sort" % basename)
        print("Creating directory for sorting BAM file: %s" % sort_dir)
        mkdir(sort_dir)
        # Sort the BAM file
        sorted_bam_file = os.path.join(sort_dir,
                                       "%s.sorted.bam" %
                                       os.path.basename(star_bam_file)[:-4])
        # Run the sorting
        samtools_sort_cmd = Command('samtools',
                                    'sort',
                                    '-o',
                                    sorted_bam_file,
                                    star_bam_file)
        print("Running %s" % samtools_sort_cmd)
        status = samtools_sort_cmd.run_subprocess(working_dir=sort_dir)
        if status != 0:
            raise Exception("samtools sort returned non-zero exit code: %s" %
                            status)
        # Index the sorted BAM file (makes BAI file)
        sorted_bam_file_index = "%s.bai" % sorted_bam_file
        samtools_index_cmd = Command('samtools',
                                     'index',
                                     sorted_bam_file,
                                     sorted_bam_file_index)
        print("Running %s" % samtools_index_cmd)
        status = samtools_index_cmd.run_subprocess(working_dir=sort_dir)
        if status != 0:
            raise Exception("samtools index returned non-zero exit code: %s" %
                            status)
        # Move the BAM and BAI files to final location
        os.rename(sorted_bam_file,bam_file)
        os.rename(sorted_bam_file_index,"%s.bai" % bam_file)
        # Remove the temporary working directories
        for dirn in (fastqs_dir,star_dir,sort_dir):
            print("Removing %s" % dirn)
            shutil.rmtree(dirn)
        # Return the BAM file name
        return bam_file
    def finish(self):
        for bam_file in self._bam_files:
            if os.path.exists(bam_file):
                self.output.bam_files.append(bam_file)

class RunRSeQCInferExperiment(PipelineTask):
    """
    Run RSeQC's 'infer_experiment.py' on BAM files

    Given a list of BAM files, for each file runs the
    RSeQC 'infer_experiment.py' utility
    (http://rseqc.sourceforge.net/#infer-experiment-py).

    The log for each run is written to a file called
    '<BASENAME>.infer_experiment.log'; the data are
    also extracted and put into an output parameter
    for direct consumption by downstream tasks.
    """
    def init(self,bam_files,reference_gene_model,out_dir):
        """
        Initialise the RunRSeQCInferExperiment task

        Arguments:
          bam_files (list): list of paths to BAM files
            to run infer_experiment.py on
          reference_gene_model (str): path to BED file
            with the reference gene model data
          out_dir (str): path to a directory where the
            output files will be written

        Outputs:
          experiments: a dictionary with BAM files as
            keys; each value is another dictionary with
            keys 'paired_end' (True for paired-end data,
            False for single-end), 'reverse', 'forward'
            and 'unstranded' (fractions of reads mapped
            in each configuration).
        """
        # Conda dependencies
        self.conda("rseqc=4.0.0",
                   "r-base=4")
        # Outputs
        self.add_output('experiments',Param())
    def setup(self):
        # Check for reference gene model
        if self.args.reference_gene_model:
            print("Reference gene model: %s" %
                  self.args.reference_gene_model)
        else:
            print("Reference gene model is not set, cannot "
                  "run RSeQC infer_experiment.py")
            return
        # Set up command to run infer_experiment.py
        for bam_file in self.args.bam_files:
            if not os.path.exists(os.path.join(
                    self.args.out_dir,
                    "%s.infer_experiment.log" %
                    os.path.basename(bam_file)[:-4])):
                self.add_cmd("Run RSeQC infer_experiment.py",
                             """
                             infer_experiment.py \\
                             -r {reference_gene_model} \\
                             -i {bam_file} >{basename}.infer_experiment.log
                             """.format(
                                 reference_gene_model=\
                                 self.args.reference_gene_model,
                                 bam_file=bam_file,
                                 basename=os.path.basename(bam_file)[:-4]))
    def finish(self):
        if not self.args.reference_gene_model:
            return
        outputs = dict()
        for bam_file in self.args.bam_files:
            infer_expt_log = os.path.join(self._working_dir,
                                          "%s.infer_experiment.log" %
                                          os.path.basename(bam_file)[:-4])
            # Copy to final destination
            if os.path.exists(infer_expt_log):
                if not os.path.exists(self.args.out_dir):
                    os.makedirs(self.args.out_dir)
                shutil.copy(infer_expt_log,self.args.out_dir)
            # Load data from log file
            infer_expt_log = os.path.join(self.args.out_dir,
                                          os.path.basename(infer_expt_log))
            infer_expt = InferExperiment(infer_expt_log)
            # Dump data associated with previous BAM
            outputs[bam_file] = {
                'paired_end': infer_expt.paired_end,
                'unstranded': infer_expt.unstranded,
                'forward': infer_expt.forward,
                'reverse': infer_expt.reverse,
            }
        # Set output
        self.output.experiments.set(outputs)

class RunRSeQCGenebodyCoverage(PipelineTask):
    """
    Run RSeQC's 'genebody_coverage.py' on BAM files

    Given a collection of BAM files, runs the RSeQC
    'genebody_coverage.py' utility
    (http://rseqc.sourceforge.net/#genebody-coverage-py).
    """
    def init(self,bam_files,reference_gene_model,out_dir,name="rseqc"):
        """
        Initialise the RunRSeQCGenebodyCoverage task

        Arguments:
          bam_files (list): list of paths to BAM files
            to run genebody_coverage.py on
          reference_gene_model (str): path to BED file
            with the reference gene model data
          out_dir (str): path to a directory where the
            output files will be written
          name (str): optional basename for the output
            files (defaults to 'rseqc')
        """
        # Conda dependencies
        self.conda("rseqc=4.0.0",
                   "r-base=4")
    def setup(self):
        # Check if outputs already exist
        outputs_exist = True
        for f in rseqc_genebody_coverage_output(self.args.name,
                                                self.args.out_dir):
            outputs_exist = (outputs_exist and os.path.exists(f))
        if outputs_exist:
            print("All outputs exist already, nothing to do")
            return
        # Check for reference gene model
        if self.args.reference_gene_model:
            print("Reference gene model: %s" %
                  self.args.reference_gene_model)
        else:
            print("Reference gene model is not set, cannot run RSeQC "
                  "genebody_coverage.py")
            return
        # Set up command to run genebody_coverage.py
        self.add_cmd("Run RSeQC genebody_coverage.py",
                     """
                     geneBody_coverage.py \\
                         -r {reference_gene_model} \\
                         -i {bam_files} \\
                         -f png \\
                         -o {basename}
                     """.format(
                         reference_gene_model=self.args.reference_gene_model,
                         bam_files=','.join(self.args.bam_files),
                         basename=self.args.name))
    def finish(self):
        if not self.args.reference_gene_model:
            return
        # Copy outputs to final location
        if not os.path.exists(self.args.out_dir):
            os.makedirs(self.args.out_dir,exist_ok=True)
        for f in rseqc_genebody_coverage_output(self.args.name,
                                                self.args.out_dir):
            if not os.path.exists(f):
                # Copy new version to ouput location
                shutil.copy(os.path.basename(f),self.args.out_dir)

class RunPicardCollectInsertSizeMetrics(PipelineTask):
    """
    Run Picard 'CollectInsertSizeMetrics' on BAM files

    Given a list of BAM files, for each file first runs
    the Picard 'CleanSam' utility (to remove alignments
    that would otherwise cause problems for the insert
    size calculations) and then 'CollectInsertSizeMetrics'
    to generate the insert size metrics.

    Note that this task should only be run on BAM files
    with paired-end data.
    """
    def init(self,bam_files,out_dir,bam_properties):
        """
        Initialise the RunPicardCollectInsertSizeMetrics
        task

        Arguments:
          bam_files (list): list of paths to BAM files
            to run CollectInsertSizeMetrics on
          out_dir (str): path to a directory where the
            output files will be written
          bam_properties (mapping): properties for each
            BAM file from RSeQC 'infer_experiment.py'
            (used to determine if BAM is paired and
            what the strand-specificity is)R
        """
        # Conda dependencies
        self.conda("picard=2.27.1",
                   "r-base=4")
    def setup(self):
        # Filter list of BAM files down to those which have
        # associated properties, and which are paired-end
        if self.args.bam_properties:
            self.bam_files = list(
                filter(lambda f: f in self.args.bam_properties and
                       self.args.bam_properties[f]['paired_end'],
                       self.args.bam_files))
            if not self.bam_files:
                print("No paired-end BAM files found")
        else:
            self.bam_files = list()
            print("No properties for BAM files, cannot run "
                  "CollectInsertSizeMetrics")
        # Set up commands to run CleanSam and
        # CollectInsertSizeMetrics for each BAM file
        for bam in self.bam_files:
            # Check if outputs already exist
            outputs_exist = True
            for f in picard_collect_insert_size_metrics_output(
                    bam,
                    self.args.out_dir):
                outputs_exist = (outputs_exist and os.path.exists(f))
            if outputs_exist:
                # Skip this BAM
                continue
            # Add command to get insert sizes
            self.add_cmd("%s: collect insert size metrics" %
                         os.path.basename(bam),
                         """
                         mkdir tmp
                         picard CleanSam \\
                             -I {bam} \\
                             -O tmp/{basename}.bam \\
                             -XX:ActiveProcessorCount={nslots} \\
                         && \\
                         picard CollectInsertSizeMetrics \\
                             -I tmp/{basename}.bam \\
                             -O {basename}.insert_size_metrics.txt \\
                             -H {basename}.insert_size_histogram.pdf \\
                             -XX:ActiveProcessorCount={nslots}
                         """.format(bam=bam,
                                    basename=os.path.basename(bam)[:-4],
                                    nslots=self.runner_nslots))
    def finish(self):
        # Check if there were any files processed
        if not self.bam_files:
            return
        # Copy outputs to final location
        if not os.path.exists(self.args.out_dir):
            print("Creating output dir '%s'" % self.args.out_dir)
            os.makedirs(self.args.out_dir)
        for bam in self.bam_files:
            for f in picard_collect_insert_size_metrics_output(
                    bam,
                    self.args.out_dir):
                if not os.path.exists(f):
                    # Copy new version to ouput location
                    shutil.copy(os.path.basename(f),self.args.out_dir)

class CollateInsertSizes(PipelineTask):
    """
    Collate insert size metrics data from multiple BAMs

    Gathers together the Picard insert size data from a
    set of BAM files and puts them into a single TSV
    file.
    """
    def init(self,bam_files,picard_out_dir,out_file,delimiter='\t'):
        """
        Initialise the CollateInsertSizes task

        Arguments:
          bam_files (list): list of paths to BAM files
            to get associated insert size data for
          picard_out_dir (str): path to the directory
            containing the Picard CollectInsertSizeMetrics
            output files
          out_file (str): path to the output TSV file
          delimiter (str): specify the delimiter to use
            in the output file
        """
        pass
    def setup(self):
        # Set up a TabFile instance for the collated data
        tf = TabFile(column_names=("Bam file",
                                   "Mean insert size",
                                   "Standard deviation",
                                   "Median insert size",
                                   "Median absolute deviation"))
        metrics_files = []
        for bam in self.args.bam_files:
            # Get metrics file associated with this BAM file
            outputs = list(filter(lambda f:
                                  f.endswith('.txt') and
                                  os.path.exists(f),
                                  picard_collect_insert_size_metrics_output(
                                      os.path.basename(bam)[:-4],
                                      prefix=self.args.picard_out_dir)))
            if not outputs:
                # No metrics located
                print("%s: no associated Picard insert size "
                      "metrics file found in %s" %
                      (bam,
                       self.args.picard_out_dir))
                continue
            metrics_files.append(outputs[0])
        # Check there is data to collate
        if not metrics_files:
            print("no insert size metrics files recovered")
            return
        # Set up a TabFile instance for the collated data
        tf = TabFile(column_names=("Bam file",
                                   "Mean insert size",
                                   "Standard deviation",
                                   "Median insert size",
                                   "Median absolute deviation"))
        for metrics_file in metrics_files:
            # Get mean and median insert sizes
            insert_size_metrics = CollectInsertSizeMetrics(metrics_file)
            tf.append(data=(
                os.path.basename(bam),
                insert_size_metrics.metrics['MEAN_INSERT_SIZE'],
                insert_size_metrics.metrics['STANDARD_DEVIATION'],
                insert_size_metrics.metrics['MEDIAN_INSERT_SIZE'],
                insert_size_metrics.metrics['MEDIAN_ABSOLUTE_DEVIATION']
            ))
        # Output to file
        print("Writing to %s" % self.args.out_file)
        tf.write(self.args.out_file,
                 include_header=True,
                 delimiter=self.args.delimiter)

class RunQualimapRnaseq(PipelineTask):
    """
    Run Qualimap's 'rnaseq' module on BAM files

    Given a list of BAM files, for each file runs the
    Qualimap 'rnaseq' module
    (http://qualimap.conesalab.org/doc_html/command_line.html#rna-seq-qc)
    """
    def init(self,bam_files,feature_file,out_dir,
             bam_properties):
        """
        Initialise the RunQualimapRnaseq task

        Arguments:
          bam_files (list): list of paths to BAM files
            to run Qualimap rnaseq on
          feature_file (str): path to GTF file with the
            reference annotation data
          out_dir (str): path to a directory where the
            output files will be written
          bam_properties (mapping): properties for each
            BAM file from RSeQC 'infer_experiment.py'
            (used to determine if BAM is paired and
            what the strand-specificity is)
        """
        self.conda("qualimap=2.2")
        self.java_mem_size = '8G'
    def setup(self):
        # Check for feature file
        if self.args.feature_file:
            print("Feature file: %s" % self.args.feature_file)
        else:
            print("Feature file is not set, cannot run Qualimap rnaseq")
            return
        # Check for BAM file properties
        if self.args.bam_properties:
            self.bam_files = self.args.bam_files
        else:
            print("No properties for BAM files, cannot run Qualimap "
                  "rnaseq")
            return
        # Set up Qualimap rnaseq for each BAM file
        for bam in self.args.bam_files:
            # Output directory for individual BAM file
            bam_name = os.path.basename(bam)[:-4]
            out_dir = os.path.join(self.args.out_dir,bam_name)
            # Check for existing outputs
            outputs_exist = True
            for f in qualimap_rnaseq_output(out_dir):
                outputs_exist = (outputs_exist and os.path.exists(f))
            if outputs_exist:
                # Skip running Qualimap for this BAM
                continue
            # Get properties for BAM file
            paired_end = self.args.bam_properties[bam]['paired_end']
            unstranded = self.args.bam_properties[bam]['unstranded']
            forward = self.args.bam_properties[bam]['forward']
            reverse = self.args.bam_properties[bam]['reverse']
            # Set sequencing protocol (aka strand specificity)
            # Qualimap sequencing protocol can be
            # 'strand-specific-forward', 'strand-specific-reverse',
            # or 'non-strand-specific'
            if reverse > forward and reverse > unstranded:
                seq_protocol = "strand-specific-reverse"
            elif forward > reverse and forward > unstranded:
                seq_protocol = "strand-specific-forward"
            else:
                seq_protocol = "non-strand-specific"
            # Run Qualimap
            self.add_cmd("Run qualimap rnaseq on %s" %
                         os.path.basename(bam),
                         """
                         export _JAVA_OPTIONS="-XX:ParallelGCThreads={nthreads} -Xmx{java_mem_size}"
                         qualimap rnaseq \\
                             -bam {bam} \\
                             -gtf {feature_file} \\
                             -p {sequencing_protocol} \\
                             {paired} \\
                             -outdir {out_dir} \\
                             -outformat HTML \\
                             --java-mem-size={java_mem_size}
                         """.format(
                             bam=bam,
                             feature_file=self.args.feature_file,
                             sequencing_protocol=seq_protocol,
                             paired=('-pe' if paired_end else ''),
                             out_dir=bam_name,
                             nthreads=self.runner_nslots,
                             java_mem_size=self.java_mem_size))
    def finish(self):
        if not self.args.feature_file:
            return
        if not self.args.bam_properties:
            return
        for bam in self.args.bam_files:
            # Check outputs for each BAM
            bam_name = os.path.basename(bam)[:-4]
            out_dir = os.path.join(self.args.out_dir,bam_name)
            outputs_exist = True
            for f in qualimap_rnaseq_output(out_dir):
                outputs_exist = (outputs_exist and os.path.exists(f))
            if outputs_exist:
                print("outputs already exist for %s" % bam_name)
                continue
            os.makedirs(self.args.out_dir,exist_ok=True)
            print("copying outputs for %s" % bam_name)
            if os.path.exists(out_dir):
                # Remove existing (incomplete) outputs
                shutil.rmtree(out_dir)
            shutil.copytree(bam_name,out_dir)

class VerifyQC(PipelineFunctionTask):
    """
    Verify outputs from the QC pipeline
    """
    def init(self,project,qc_dir):
        """
        Initialise the VerifyQC task.

        Arguments:
          project (AnalysisProject): project to update the
            number of cells for
          qc_dir (str): directory for QC outputs (defaults
            to subdirectory 'qc' of project directory)
        """
        pass
    def setup(self):
        # Run the 'verify_project' function
        self.add_call(
            "Verify QC outputs for %s" % self.args.project.name,
            verify_project,
            self.args.project,
            self.args.qc_dir)
    def finish(self):
        # Report the verification output
        for line in self.stdout.split('\n'):
            if not line.startswith("#### "):
                print(line)
        # Check the status
        verified = self.result()[0]
        if not verified:
            self.fail(message="Failed to verify QC outputs")
        print("Verified QC outputs")

class ReportQC(PipelineTask):
    """
    Generate the QC report
    """
    def init(self,project,qc_dir,report_html=None,fastq_dir=None,
             multiqc=False,force=False,zip_outputs=True):
        """
        Initialise the ReportQC task.

        Arguments:
          project (AnalysisProject): project to generate
            QC report for
          qc_dir (str): directory for QC outputs (defaults
            to subdirectory 'qc' of project directory)
          report_html (str): set the name of the output
            HTML file for the report
          fastq_dir (str): directory holding Fastq files
            (defaults to current fastq_dir in project)
          multiqc (bool): if True then also generate
            MultiQC report (default: don't run MultiQC)
          force (bool): if True then force HTML report to
            be generated even if QC outputs fail
            verification (default: don't write report)
          zip_outputs (bool): if True then also generate
            a ZIP archive of the QC reports
        """
        self.conda("multiqc=1.8",
                   "pillow")
        # Specify Python version to use to avoid similar
        # issue as reported here:
        # https://github.com/ewels/MultiQC/issues/1413
        self.conda("python=3.8")
    def setup(self):
        # Check for 10x multiome libraries file with linked projects
        libraries_file = os.path.join(self.args.project.dirn,
                                      "10x_multiome_libraries.info")
        if os.path.exists(libraries_file):
            # Read the file and get sample associations between
            # 'local' and 'remote' datasets
            libraries = MultiomeLibraries(libraries_file)
            linked_projects = libraries.linked_projects()
        else:
            linked_projects = None
        # Build the reportqc.py command
        project = self.args.project
        fastq_dir = self.args.fastq_dir
        qc_base = os.path.basename(self.args.qc_dir)
        if self.args.report_html is None:
            out_file = '%s_report.html' % qc_base
        else:
            out_file = self.args.report_html
        if not os.path.isabs(out_file):
            out_file = os.path.join(project.dirn,out_file)
        if project.info.run is None:
            title = "%s" % project.name
        else:
            title = "%s/%s" % (project.info.run,
                               project.name)
        if fastq_dir is not None:
            title = "%s (%s)" % (title,fastq_dir)
        title = "%s: QC report" % title
        cmd = PipelineCommandWrapper(
            "Generate QC report for %s" % project.name,
            "reportqc.py",
            "--filename",out_file,
            "--title",title)
        if self.args.force:
            cmd.add_args("--force")
        if self.args.qc_dir is not None:
            if not os.path.isabs(self.args.qc_dir):
                # QC dir specified as a relative path
                cmd.add_args("--qc_dir",self.args.qc_dir)
        if fastq_dir is not None:
            cmd.add_args("--fastq_dir",fastq_dir)
        if self.args.multiqc:
            cmd.add_args("--multiqc")
        if self.args.zip_outputs:
            cmd.add_args("--zip")
        # Add the primary project/QC directory
        if self.args.qc_dir:
            if os.path.isabs(self.args.qc_dir):
                # QC dir specified as an absolute path
                # Use this directly
                cmd.add_args(self.args.qc_dir)
            else:
                # QC dir specified as a relative path
                # Use the project directory
                cmd.add_args(project.dirn)
        # Append the additional linked projects
        if linked_projects:
            cmd.add_args(*[p.dirn for p in linked_projects])
        self.add_cmd(cmd)
