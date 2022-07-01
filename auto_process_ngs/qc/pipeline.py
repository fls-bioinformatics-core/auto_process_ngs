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
from bcftbx.JobRunner import SimpleJobRunner
from bcftbx.utils import mkdir
from bcftbx.utils import mkdirs
from bcftbx.utils import find_program
from ..analysis import copy_analysis_project
from ..bcl2fastq.pipeline import Get10xPackage
from ..bcl2fastq.pipeline import FunctionParam
from ..command import Command
from ..fastq_utils import pair_fastqs_by_name
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
from .outputs import check_fastq_screen_outputs
from .outputs import check_fastqc_outputs
from .outputs import check_fastq_strand_outputs
from .outputs import check_cellranger_count_outputs
from .outputs import check_cellranger_atac_count_outputs
from .outputs import check_cellranger_arc_count_outputs
from .protocols import determine_qc_protocol
from .protocols import get_read_numbers
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
        self.add_param('fastq_subset',type=int)
        self.add_param('cellranger_exe',type=str)
        self.add_param('cellranger_reference_dataset',type=str)
        self.add_param('cellranger_out_dir',type=str)
        self.add_param('cellranger_chemistry',type=str)
        self.add_param('cellranger_force_cells',type=int)
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
                    sample_pattern=None,log_dir=None):
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

    def add_cellranger_count(self,project_name,project,qc_dir,
                             organism,fastq_dir,qc_protocol,chemistry,
                             force_cells,log_dir,samples=None,
                             fastq_dirs=None,reference_dataset=None,
                             required_tasks=None):
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

    def run(self,nthreads=None,fastq_screens=None,star_indexes=None,
            fastq_subset=None,cellranger_chemistry='auto',
            cellranger_force_cells=None,cellranger_transcriptomes=None,
            cellranger_premrna_references=None,
            cellranger_atac_references=None,
            cellranger_arc_references=None,cellranger_jobmode='local',
            cellranger_maxjobs=None,cellranger_mempercore=None,
            cellranger_jobinterval=None,cellranger_localcores=None,
            cellranger_localmem=None,cellranger_exe=None,
            cellranger_reference_dataset=None,
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
        self.conda("star=2.4.2a",
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
             qc_protocol=None,cellranger_version=None,
             cellranger_ref_data=None,verbose=False):
        """
        Initialise the CheckCellrangerCountOutputs task.

        Arguments:
          project (AnalysisProject): project to run
            QC for
          fastq_dir (str): directory holding Fastq files
            (defaults to current fastq_dir in project)
          samples (list): list of samples to restrict
            checks to (all samples in project are checked
            by default)
          qc_dir (str): top-level QC directory to look
            for 'count' QC outputs (e.g. metrics CSV and
            summary HTML files)
          qc_protocol (str): QC protocol to use
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
        self.add_call("Check cellranger count outputs for %s"
                      % self.args.project.name,
                      check_outputs,
                      self.args.project,
                      self.args.qc_dir,
                      prefix=prefix)
    def finish(self):
        if not self.args.cellranger_ref_data:
            # No reference data, nothing to check
            print("No reference data to check against")
            return
        # Collect the sample names with missing outputs
        for result in self.result():
            for smpl in result:
                if not self.args.samples or smpl in self.args.samples:
                    self.output.samples.append(smpl)
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
        self.output.fastq_dir.set(self.args.project.fastq_dir)

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
