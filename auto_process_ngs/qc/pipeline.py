#!/usr/bin/env python
#
#     qc.pipeline.py: pipelines for running QC
#     Copyright (C) University of Manchester 2019-2024 Peter Briggs
#

"""
Pipeline components for running the QC pipeline.

Pipeline classes:

- QCPipeline

Pipeline task classes:

- SetupQCDirs
- SplitFastqsByLane
- GetSequenceDataSamples
- GetSequenceDataFastqs
- UpdateQCMetadata
- VerifyFastqs
- SetCellCountFromCellranger
- GetReferenceDataset
- GetBAMFile
- ConvertGTFToBed
- ReportQC
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
from bcftbx.FASTQFile import FastqIterator
from bcftbx.TabFile import TabFile
from bcftbx.utils import mkdir
from bcftbx.utils import mkdirs
from bcftbx.utils import find_program
from bcftbx.ngsutils import getreads
from bcftbx.ngsutils import getreads_subset
from ..analysis import AnalysisFastq
from ..analysis import copy_analysis_project
from ..command import Command
from ..fastq_utils import group_fastqs_by_name
from ..fastq_utils import remove_index_fastqs
from ..pipeliner import Pipeline
from ..pipeliner import PipelineTask
from ..pipeliner import PipelineFunctionTask
from ..pipeliner import PipelineCommandWrapper
from ..pipeliner import PipelineParam as Param
from ..pipeliner import ListParam
from ..pipeliner import PipelineFailure
from ..tenx.cellplex import CellrangerMultiConfigCsv
from ..tenx.multiome import MultiomeLibraries
from ..tenx.utils import add_cellranger_args
from ..utils import get_organism_list
from .modules.cellranger_atac_count import CellrangerAtacCount
from .modules.cellranger_arc_count import CellrangerArcCount
from .modules.cellranger_count import CellrangerCount
from .modules.cellranger_multi import CellrangerMulti
from .modules.cellranger_multi import GetCellrangerMultiConfig
from .modules.fastqc import Fastqc
from .modules.fastq_screen import FastqScreen
from .modules.picard_insert_size_metrics import PicardInsertSizeMetrics
from .modules.qualimap_rnaseq import QualimapRnaseq
from .modules.rseqc_genebody_coverage import RseqcGenebodyCoverage
from .modules.rseqc_infer_experiment import RseqcInferExperiment
from .modules.sequence_lengths import SequenceLengths
from .modules.strandedness import Strandedness
from .protocols import determine_qc_protocol
from .protocols import fetch_protocol_definition
from .utils import get_bam_basename
from .utils import get_seq_data_samples
from .utils import set_cell_count_for_project
from .verification import parse_qc_module_spec
from .verification import verify_project

# Module specific logger
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

######################################################################
# Supported QC module classes
######################################################################

from .qc_modules import QC_MODULES

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

        # Default log directory
        self._default_log_dir = None

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
        self.add_param('force_star_index',type=str)
        self.add_param('force_gtf_annotation',type=str)
        self.add_param('legacy_screens',type=bool,value=False)

        # Define runners
        self.add_runner('verify_runner')
        self.add_runner('fastq_screen_runner')
        self.add_runner('fastqc_runner')
        self.add_runner('star_runner')
        self.add_runner('picard_runner')
        self.add_runner('qualimap_runner')
        self.add_runner('rseqc_runner')
        self.add_runner('cellranger_count_runner')
        self.add_runner('cellranger_multi_runner')
        self.add_runner('report_runner')

        # Define module environment modules
        self.add_envmodules('fastqc')
        self.add_envmodules('fastq_screen')
        self.add_envmodules('fastq_strand')
        self.add_envmodules('cellranger')
        self.add_envmodules('report_qc')

    def add_task(self,task,requires=(),**kws):
        """
        Override base class method

        Automatically set log dir when tasks are added
        """
        updated_kws = { x: kws[x] for x in kws }
        if 'log_dir' not in updated_kws and self.default_log_dir:
            kws['log_dir'] = self.default_log_dir
        return Pipeline.add_task(self,
                                 task,
                                 requires=requires,
                                 **kws)

    def set_default_log_dir(self,log_dir):
        """
        Set the default log directory for tasks
        """
        self._default_log_dir = log_dir

    @property
    def default_log_dir(self):
        """
        Return current value of default log dir
        """
        return self._default_log_dir

    def add_project(self,project,protocol,qc_dir=None,organism=None,
                    fastq_dir=None,report_html=None,multiqc=False,
                    sample_pattern=None,log_dir=None,convert_gtf=True,
                    verify_fastqs=False,split_fastqs_by_lane=False):
        """
        Add a project to the QC pipeline

        Arguments:
          project (AnalysisProject): project to run
            QC for
          protocol (QCProtocol): QC protocol to use
          qc_dir (str): directory for QC outputs (defaults
            to subdirectory 'qc' of project directory)
          organism (str): organism(s) for project
            (defaults to organism defined in project
            metadata)
          fastq_dir (str): directory holding Fastq files
            (defaults to primary fastq_dir in project)
          multiqc (bool): if True then also run MultiQC
            (default is not to run MultiQC)
          sample_pattern (str): glob-style pattern to
            match a subset of projects and samples (not
            implemented)
          log_dir (str): directory to write log files to
            (defaults to 'logs' subdirectory of the QC
            directory)
          convert_gtf (bool): if True then convert
            GTF files to BED for 'infer_experiment.py'
            (default; otherwise only use the explicitly
            defined BED files)
          verify_fastqs (bool): if True then verify
            Fastq integrity as part of the pipeline
            (default: False, skip verification)
          split_fastqs_by_lanes (bool): if True then
            split input Fastqs into lanes and run QC
            as per-lane (default: False, don't split
            QC by lanes)
        """
        ###################
        # Do internal setup
        ###################

        self.report("Adding project: %s" % project.name)

        # QC modules
        qc_modules = protocol.qc_modules

        # Read numbers for sequence data and QC
        read_numbers = protocol.read_numbers

        # Determine whether sequence data are paired
        paired = (len(read_numbers.seq_data) > 1)

        # Determine if BAM files are required
        require_bam_files = False
        for m in QC_MODULES:
            if m.name in protocol.qc_module_names:
                require_bam_files = (m.require_bam_files or require_bam_files)
                if require_bam_files:
                    break

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

        # Sanitised organism name
        organism_name = str(organism).\
                        strip().\
                        lower().\
                        replace(' ','_')

        # Report details
        self.report("-- Protocol   : %s" % protocol.name)
        self.report("-- Directory  : %s" % project.dirn)
        self.report("-- Fastqs dir : %s" % project.fastq_dir)
        self.report("-- QC dir     : %s" % qc_dir)
        self.report("-- Library    : %s" % project.info.library_type)
        self.report("-- SC platform: %s" % project.info.single_cell_platform)
        self.report("-- Organism   : %s" % organism)
        self.report("-- Report     : %s" % report_html)
        self.report("Reads")
        self.report("-- Seq data   : %s" % protocol.seq_data_reads)
        self.report("-- Index      : %s" % protocol.index_reads)
        self.report("QC modules")
        for qc_module in qc_modules:
            self.report("-- %s" % qc_module)

        # Sort out the log directory for this project
        if log_dir is None:
            log_dir = os.path.join(qc_dir,'logs')
        else:
            log_dir = os.path.abspath(log_dir)

        # Set as pipeline default (will be changed for
        # each project that is added)
        self.set_default_log_dir(log_dir)

        ####################
        # Build the pipeline
        ####################

        project_name = "%s%s" % (project.name,
                                 ":%s" % os.path.basename(fastq_dir)
                                 if fastq_dir is not None
                                 else '')

        startup_tasks = []

        # Set up QC dirs
        setup_qc_dirs = SetupQCDirs(
            "%s: set up QC directories" % project_name,
            project,
            qc_dir,
            log_dir=log_dir,
            protocol=protocol
        )
        self.add_task(setup_qc_dirs)
        startup_tasks.append(setup_qc_dirs)

        # Characterise samples
        get_seq_data = GetSequenceDataSamples(
            "%s: identify sequence data samples" % project.name,
            project,
            fastq_attrs=project.fastq_attrs)
        self.add_task(get_seq_data)

        # Build a dictionary of QC metadata items to
        # update
        qc_metadata = dict(protocol=protocol.name,
                           protocol_summary=protocol.summarise(),
                           protocol_specification=repr(protocol),
                           organism=organism,
                           seq_data_samples=\
                           get_seq_data.output.seq_data_samples,
                           fastq_dir=project.fastq_dir,
                           fastqs_split_by_lane=split_fastqs_by_lane)

        # Verify Fastqs
        if verify_fastqs:
            verify_fqs = VerifyFastqs(
                "%s: verify Fastqs" % project_name,
                project,
                fastq_attrs=project.fastq_attrs)
            self.add_task(verify_fqs,
                          requires=(setup_qc_dirs,),
                          runner=self.runners['verify_runner'])
            startup_tasks.append(verify_fqs)

        # Update QC metadata
        update_qc_metadata = UpdateQCMetadata(
            "%s: update QC metadata" % project_name,
            project,
            qc_dir,
            qc_metadata,
            legacy_screens=self.params.legacy_screens)
        self.add_task(update_qc_metadata,
                      requires=(setup_qc_dirs,
                                get_seq_data))

        # Split Fastqs by lane for QC?
        if split_fastqs_by_lane:
            split_fastqs = SplitFastqsByLane(
                "%s: split Fastqs by lane" % project_name,
                project,
                os.path.join(qc_dir,'__fastqs.split'))
            self.add_task(split_fastqs,
                          requires=(setup_qc_dirs,))
            # Subsequent start-up tasks should wait on this task
            startup_tasks.append(split_fastqs)
            # Ensure QC metadata also waits for this task
            update_qc_metadata.requires(split_fastqs)
            fastqs_in = split_fastqs.output.fastqs
        else:
            fastqs_in = project.fastqs
        qc_metadata['fastqs'] = fastqs_in

        # Verify QC
        verify_qc = VerifyQC(
            "%s: verify QC outputs" % project_name,
            project,
            qc_dir,
            str(protocol),
            fastqs=fastqs_in)
        self.add_task(verify_qc,
                      requires=startup_tasks,
                      runner=self.runners['verify_runner'])

        # Make QC report
        report_qc = ReportQC(
            "%s: make QC report" % project_name,
            project,
            qc_dir,
            report_html=report_html,
            force=True
        )
        self.add_task(report_qc,
                      requires=(verify_qc,),
                      runner=self.runners['report_runner'],
                      envmodules=self.envmodules['report_qc'])

        # Run MultiQC
        if multiqc:
            run_multiqc = RunMultiQC(
                f"{project_name}: make MultiQC report",
                project,
                qc_dir
            )
            self.add_task(run_multiqc,
                          requires=(verify_qc,),
                          runner=self.runners['report_runner'])
            report_qc.requires(run_multiqc)

        # Get sequence data Fastqs
        get_seq_fastqs = GetSequenceDataFastqs(
            "%s: get sequence data Fastqs" % project.name,
            project,
            os.path.join(qc_dir,'__fastqs'),
            read_range=protocol.read_range,
            samples=get_seq_data.output.seq_data_samples,
            fastq_attrs=project.fastq_attrs,
            fastqs=fastqs_in)
        self.add_task(get_seq_fastqs,
                      requires=startup_tasks)

        # Set up tasks to generate and characterise BAM files
        if require_bam_files:

            # Get STAR index
            get_star_index = GetReferenceDataset(
                "%s: get STAR index for '%s'" % (project.name,
                                                 organism),
                organism,
                self.params.star_indexes,
                force_reference=self.params.force_star_index)
            self.add_task(get_star_index)
            qc_metadata['star_index'] = \
                get_star_index.output.reference_dataset

            # Fetch BAM files
            get_bam_files = GetBAMFiles(
                "%s: get BAM files" % project.name,
                get_seq_fastqs.output.fastqs,
                get_star_index.output.reference_dataset,
                os.path.join(qc_dir,'__bam_files',organism_name),
                self.params.fastq_subset,
                self.params.nthreads,
                reads=read_numbers.seq_data,
                include_samples=get_seq_data.output.seq_data_samples,
                fastq_attrs=project.fastq_attrs,
                verbose=self.params.VERBOSE)
            self.add_task(get_bam_files,
                          requires=startup_tasks,
                          runner=self.runners['star_runner'])

            # Get reference gene model for RSeQC
            get_reference_gene_model = GetReferenceDataset(
                "%s: get RSeQC reference gene model for '%s'" %
                (project.name,
                 organism),
                organism,
                self.params.annotation_bed_files)
            self.add_task(get_reference_gene_model)
            qc_metadata['annotation_bed'] = \
                get_reference_gene_model.output.reference_dataset

            # Get GTF annotation
            get_annotation_gtf = GetReferenceDataset(
                "%s: get GTF annotation for '%s'" % (project.name,
                                                     organism),
                organism,
                self.params.annotation_gtf_files,
                force_reference=self.params.force_gtf_annotation)
            self.add_task(get_annotation_gtf)
            qc_metadata['annotation_gtf'] = \
                get_annotation_gtf.output.reference_dataset

            # BED annotation for infer experiment
            if convert_gtf:
                # Convert GTF annotation to BED
                get_bed_annotation_from_gtf = ConvertGTFToBed(
                    "%s: convert GTF annotation to BED for '%s'" %
                    (project.name,
                     organism),
                    get_annotation_gtf.output.reference_dataset,
                    os.path.join(qc_dir,'%s.annotation.bed' %
                                 organism_name))
                self.add_task(get_bed_annotation_from_gtf,
                              requires=startup_tasks)
                reference_gene_model_file = get_bed_annotation_from_gtf.\
                                            output.bed_file
            else:
                # Use pre-defined BED file
                reference_gene_model_file = get_reference_gene_model.\
                                            output.reference_dataset

            # Run RSeQC infer experiment
            rseqc_infer_experiment = RseqcInferExperiment.add_to_pipeline(
                self,
                project_name,
                qc_dir,
                get_bam_files.output.bam_files,
                reference_gene_model_file,
                organism_name,
                rseqc_runner=self.runners['rseqc_runner'])
            self.add_task(rseqc_infer_experiment)
            verify_qc.requires(rseqc_infer_experiment)

        ################
        # Add QC modules
        ################

        for qc_module in qc_modules:

            qc_module_name,qc_module_params = parse_qc_module_spec(qc_module)

            ##################################
            # Fastq sequence length statistics
            ##################################
            if qc_module_name == 'sequence_lengths':
                get_seq_lengths = SequenceLengths.add_to_pipeline(
                    self,
                    project_name,
                    project,
                    qc_dir,
                    read_numbers=read_numbers.qc,
                    fastqs=fastqs_in,
                    require_tasks=startup_tasks,
                    compute_runner=self.runners['fastqc_runner'])
                verify_qc.requires(get_seq_lengths)

            #############
            # FastqScreen
            #############
            if qc_module_name == 'fastq_screen':
                run_fastq_screen = FastqScreen.add_to_pipeline(
                    self,
                    project_name,
                    project,
                    qc_dir,
                    self.params.fastq_screens,
                    get_seq_fastqs.output.fastqs,
                    read_numbers.seq_data,
                    include_samples=get_seq_data.output.seq_data_samples,
                    nthreads=self.params.nthreads,
                    fastq_subset=self.params.fastq_subset,
                    legacy=self.params.legacy_screens,
                    verbose=self.params.VERBOSE,
                    requires_tasks=startup_tasks,
                    verify_runner=self.runners['verify_runner'],
                    compute_runner=self.runners['fastq_screen_runner'],
                    envmodules=self.envmodules['fastq_screen'])
                qc_metadata['fastq_screens'] = self.params.fastq_screens
                qc_metadata['legacy_screens'] = self.params.legacy_screens
                verify_qc.requires(run_fastq_screen)

            ########
            # FastQC
            ########
            if qc_module_name == 'fastqc':
                run_fastqc = Fastqc.add_to_pipeline(
                    self,
                    project_name,
                    project,
                    qc_dir,
                    read_numbers.qc,
                    fastqs_in,
                    verbose=self.params.VERBOSE,
                    nthreads=self.params.nthreads,
                    require_tasks=startup_tasks,
                    verify_runner=self.runners['verify_runner'],
                    compute_runner=self.runners['fastqc_runner'],
                    envmodules=self.envmodules['fastqc'])
                verify_qc.requires(run_fastqc)

            ##############
            # Fastq_strand
            ##############
            if qc_module_name == 'strandedness':
                run_fastq_strand = Strandedness.add_to_pipeline(
                    self,
                    project_name,
                    project,
                    qc_dir,
                    organism=organism,
                    read_numbers=read_numbers.seq_data,
                    fastqs=get_seq_fastqs.output.fastqs,
                    star_indexes=self.params.star_indexes,
                    include_samples=get_seq_data.output.seq_data_samples,
                    nthreads=self.params.nthreads,
                    fastq_subset=self.params.fastq_subset,
                    require_tasks=startup_tasks,
                    verify_runner=self.runners['verify_runner'],
                    compute_runner=self.runners['star_runner'],
                    envmodules=self.envmodules['fastq_strand'],
                    verbose=self.params.VERBOSE
                )
                verify_qc.requires(run_fastq_strand)

            #############################
            # 10x single library analysis
            #############################
            if qc_module_name in ('cellranger_count',
                                  'cellranger-atac_count',
                                  'cellranger-arc_count',):
                # Set base QC module
                if qc_module_name == "cellranger_count":
                    count = CellrangerCount
                elif qc_module_name == "cellranger-atac_count":
                    count = CellrangerAtacCount
                elif qc_module_name == "cellranger-arc_count":
                    count = CellrangerArcCount

                # Set library type
                try:
                    library_type = qc_module_params['library']
                except KeyError:
                    library_type = project.info.library_type

                # Set chemistry
                try:
                    chemistry = qc_module_params['chemistry']
                except KeyError:
                    chemistry = self.params.cellranger_chemistry

                # Locate 'cellranger multi' config.csv file
                # if required
                try:
                    cellranger_use_multi_config = \
                        qc_module_params['cellranger_use_multi_config']
                except KeyError:
                    cellranger_use_multi_config = False
                if cellranger_use_multi_config:
                    get_cellranger_multi_config = GetCellrangerMultiConfig(
                        "%s: get config file for 'cellranger multi'" %
                        project_name,
                        project,
                        qc_dir
                    )
                    self.add_task(get_cellranger_multi_config,
                                  requires=startup_tasks)
                    samples = get_cellranger_multi_config.output.gex_libraries
                    fastq_dirs = get_cellranger_multi_config.output.fastq_dirs
                    reference_dataset = \
                        get_cellranger_multi_config.output.reference_data_path
                else:
                    samples = None
                    fastq_dirs = None
                    reference_dataset = \
                        self.params.cellranger_reference_dataset

                # Whether to set metadata
                try:
                    set_metadata = qc_module_params['set_metadata']
                except KeyError:
                    set_metadata = True

                # Whether to set cell count
                try:
                    set_cell_count = (set_metadata and
                                      qc_module_params['set_cell_count'])
                except KeyError:
                    set_cell_count = set_metadata

                # Run cellranger* count
                run_cellranger_count = count.add_to_pipeline(
                    self,
                    project_name,
                    project,
                    qc_dir,
                    organism,
                    fastq_dir,
                    qc_module_name,
                    library_type=library_type,
                    transcriptome_references=\
                    self.params.cellranger_transcriptomes,
                    premrna_references=\
                    self.params.cellranger_premrna_references,
                    atac_references=self.params.cellranger_atac_references,
                    multiome_references=self.params.cellranger_arc_references,
                    chemistry=chemistry,
                    force_cells=self.params.cellranger_force_cells,
                    reference_dataset=reference_dataset,
                    samples=samples,
                    fastq_dirs=fastq_dirs,
                    cellranger_exe=self.params.cellranger_exe,
                    extra_projects=self.params.cellranger_extra_projects,
                    cellranger_out_dir=self.params.cellranger_out_dir,
                    cellranger_jobmode=self.params.cellranger_jobmode,
                    cellranger_maxjobs=self.params.cellranger_maxjobs,
                    cellranger_mempercore=self.params.cellranger_mempercore,
                    cellranger_jobinterval=self.params.cellranger_jobinterval,
                    cellranger_localcores=self.params.cellranger_localcores,
                    cellranger_localmem=self.params.cellranger_localmem,
                    required_tasks=startup_tasks,
                    verify_runner=self.runners['verify_runner'],
                    cellranger_runner=self.runners['cellranger_count_runner'],
                    envmodules=self.envmodules['cellranger'],
                    verbose=self.params.VERBOSE)
                verify_qc.requires(run_cellranger_count)

                # Update metadata
                if set_metadata:
                    qc_metadata['cellranger_version'] = \
                        run_cellranger_count.output.cellranger_version
                    qc_metadata['cellranger_refdata'] = \
                        run_cellranger_count.output.cellranger_refdata
                    update_qc_metadata.requires(run_cellranger_count)

                # Set cell count
                if set_cell_count:
                    set_cellranger_cell_count = \
                        SetCellCountFromCellranger(
                            "%s: set cell count from single library analysis" %
                            project_name,
                            project,
                            qc_dir,
                            source="count"
                        )
                    self.add_task(set_cellranger_cell_count,
                                  requires=(run_cellranger_count,
                                            update_qc_metadata),)
                    verify_qc.requires(set_cellranger_cell_count)

            ################################
            # 10x cell multiplexing analysis
            ################################
            if qc_module_name == "cellranger_multi":

                run_cellranger_multi = CellrangerMulti.add_to_pipeline(
                    self,
                    project_name,
                    project,
                    qc_dir,
                    qc_module_name,
                    cellranger_exe=self.params.cellranger_exe,
                    cellranger_out_dir=self.params.cellranger_out_dir,
                    cellranger_jobmode=self.params.cellranger_jobmode,
                    cellranger_maxjobs=self.params.cellranger_maxjobs,
                    cellranger_mempercore=self.params.cellranger_mempercore,
                    cellranger_jobinterval=self.params.cellranger_jobinterval,
                    cellranger_localcores=self.params.cellranger_localcores,
                    cellranger_localmem=self.params.cellranger_localmem,
                    required_tasks=startup_tasks,
                    cellranger_runner=self.runners['cellranger_multi_runner'],
                    envmodules=self.envmodules['cellranger'],
                    working_dir=self.params.WORKING_DIR
                )
                verify_qc.requires(run_cellranger_multi)

                # Update metadata
                qc_metadata['cellranger_version'] = \
                    run_cellranger_multi.output.cellranger_version
                qc_metadata['cellranger_refdata'] = \
                    run_cellranger_multi.output.cellranger_refdata
                qc_metadata['cellranger_probeset'] = \
                    run_cellranger_multi.output.cellranger_probeset
                update_qc_metadata.requires(run_cellranger_multi)

                # Set cell count
                set_cellranger_cell_count = SetCellCountFromCellranger(
                    "%s: set cell count from cell multiplexing analysis" %
                    project_name,
                    project,
                    qc_dir,
                    source="multi"
                )
                self.add_task(set_cellranger_cell_count,
                              requires=(run_cellranger_multi,),)
                verify_qc.requires(set_cellranger_cell_count)

            ##########################
            # RSeQC gene body coverage
            ##########################
            if qc_module_name == 'rseqc_genebody_coverage':
                rseqc_gene_body_coverage = RseqcGenebodyCoverage.add_to_pipeline(
                    self,
                    project_name,
                    qc_dir,
                    get_bam_files.output.bam_files,
                    get_reference_gene_model.output.reference_dataset,
                    organism_name,
                    required_tasks=startup_tasks,
                    rseqc_runner=self.runners['rseqc_runner'])
                verify_qc.requires(rseqc_gene_body_coverage)

            ##########################
            # Picard insert sizes
            ##########################
            if qc_module_name == 'picard_insert_size_metrics' and paired:
                collate_insert_sizes = PicardInsertSizeMetrics.add_to_pipeline(
                    self,
                    project_name,
                    project,
                    qc_dir,
                    get_bam_files.output.bam_files,
                    organism_name,
                    required_tasks=startup_tasks,
                    compute_runner=self.runners['picard_runner']
                )
                verify_qc.requires(collate_insert_sizes)

            #################
            # Qualimap RNASEQ
            #################
            if qc_module_name == 'qualimap_rnaseq':
                # Run Qualimap RNA-seq analysis
                qualimap_rnaseq = QualimapRnaseq.add_to_pipeline(
                    self,
                    project_name,
                    qc_dir,
                    get_bam_files.output.bam_files,
                    get_annotation_gtf.output.reference_dataset,
                    organism_name,
                    rseqc_infer_experiment_outputs=\
                    rseqc_infer_experiment.output.experiments,
                    required_tasks=startup_tasks,
                    qualimap_runner=self.runners['qualimap_runner']
                )
                verify_qc.requires(qualimap_rnaseq)

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
            cellranger_out_dir=None,force_star_index=None,
            force_gtf_annotation=None,working_dir=None,
            log_file=None,batch_size=None,batch_limit=None,max_jobs=1,
            max_slots=None,poll_interval=5,runners=None,default_runner=None,
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
          force_star_index (str): explicitly specify STAR
            index to use (default: index is determined
            automatically)
          force_gtf_annotation (str): explicitly specify
            GTF annotation to use (default: annotation
            file is determined automatically)
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
            'fastq_screen_runner','star_runner','rseqc_runner',
            'qualimap_runner','cellranger_count_runner',
            'cellranger_multi_runner','report_runner',
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
                                  'force_star_index': force_star_index,
                                  'force_gtf_annotation': force_gtf_annotation,
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
    def init(self,project,qc_dir,log_dir=None,protocol=None):
        """
        Initialise the SetupQCDirs task

        Arguments:
          project (AnalysisProject): project to run
            QC for
          qc_dir (str): directory for QC outputs (defaults
            to subdirectory 'qc' of project directory)
          log_dir (str): directory for log files (defaults
            to 'logs' subdirectory of the QC directory
          protocol (QCProject): QC protocol being used
        """
        pass
    def setup(self):
        # Get the existing QC metadata
        qc_info = self.args.project.qc_info(self.args.qc_dir)
        # Check the QC protocol
        stored_protocol_name = qc_info.protocol
        stored_protocol_spec = qc_info.protocol_specification
        if (stored_protocol_spec is not None and
            stored_protocol_spec != repr(self.args.protocol)) or \
            (stored_protocol_name is not None and \
             stored_protocol_name != self.args.protocol.name):
            logger.warning("QC protocol mismatch between stored and "
                           "supplied protocol information for %s"
                           % self.args.project.name)
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

class SplitFastqsByLane(PipelineTask):
    """
    Split reads into multiple Fastqs according to lane
    """
    def init(self,project,out_dir):
        """
        Initialise the SplitFastqsByLane task

        Arguments:
          project (AnalysisProject): project with source
            Fastqs to split by lane
          out_dir (str): path to directory where split
            Fastqs will be written

        Outputs:
          fastqs (list): list of paths to output Fastqs
            split by lanes
        """
        self.add_output('fastqs',list())
    def setup(self):
        pre_split_fastqs = set()
        # Make output directory
        if not os.path.exists(self.args.out_dir):
            print("Making output dir: %s" % self.args.out_dir)
            os.makedirs(self.args.out_dir)
        else:
            # Identify Fastqs that have already been split
            for fq in os.listdir(self.args.out_dir):
                print("Checking existing Fastq: %s" %
                      os.path.basename(fq))
                fqname = self.args.project.fastq_attrs(fq)
                fqname.lane_number = None
                pre_split_fastqs.add(str(fqname))
        print("Pre-split Fastqs: %s" % pre_split_fastqs)
        # Make copies of Fastqs split by lane (if not already
        # present)
        for fq in self.args.project.fastqs:
            fqname = self.args.project.fastq_attrs(fq)
            fqname.lane_number = None
            fqname = str(fqname)
            if fqname in pre_split_fastqs:
                print("%s: already split" % os.path.basename(fq))
                continue
            self.add_cmd("Split %s by lane" % os.path.basename(fq),
                         """
                         echo "Making temp dir"
                         tmp_dir=$(mktemp -d --tmpdir=.)
                         cd $tmp_dir
                         echo "Moved to $(pwd)"
                         split_fastq.py {fastq}
                         for f in $(ls *.fastq) ; do
                            echo "Compressing $f"
                            gzip $f
                         done
                         echo "Moving .fastq.gz files to final dir"
                         mv -f *.fastq.gz {out_dir}
                         """.format(fastq=fq,
                                    out_dir=self.args.out_dir))
    def finish(self):
        # Collect split files
        self.output.fastqs.extend(
            sorted([os.path.join(self.args.out_dir,fq)
                    for fq in os.listdir(self.args.out_dir)
                    if fq.endswith(".fastq.gz")]))

class GetSequenceDataSamples(PipelineTask):
    """
    Identify samples with sequence (i.e. biological) data
    """
    def init(self,project,fastq_attrs):
        """
        Initialise the GetSequenceDataSamples task

        Arguments:
          project (AnalysisProject): project to get
            samples for
          fastq_attrs (BaseFastqAttrs): class to use for
            extracting data from Fastq names

        Outputs:
          seq_data_samples (list): list of samples with
            biological data
        """
        self.add_output('seq_data_samples',list())
    def setup(self):
        # Identify samples with biological data
        seq_data_samples = get_seq_data_samples(self.args.project.dirn,
                                                self.args.fastq_attrs)
        # Report
        print("Samples with sequence data:")
        for sample in seq_data_samples:
            print("- %s" % sample)
        # Set outputs
        self.output.seq_data_samples.extend(seq_data_samples)

class GetSequenceDataFastqs(PipelineTask):
    """
    Set up Fastqs with sequence (i.e. biological) data
    """
    def init(self,project,out_dir,read_range,samples,
             fastq_attrs,fastqs=None):
        """
        Initialise the GetSequenceDataFastqs task

        Arguments:
          project (AnalysisProject): project to get
            Fastqs for
          out_dir (str): path to directory to write final
            Fastq files to
          read_range (dict): mapping of read names to
            tuples of subsequence ranges
          samples (list): list of samples with sequence
            data
          fastq_attrs (BaseFastqAttrs): class to use for
            extracting data from Fastq names
          fastqs (list): optional, list of Fastq files
            (overrides Fastqs in project)

        Outputs:
          fastqs (list): list of Fastqs with biological
            data
        """
        self.conda("seqtk=1.3")
        self.add_output('fastqs',ListParam())
    def setup(self):
        # Report what will happen
        print("Read ranges:")
        for rd in sorted(list(self.args.read_range.keys())):
            rng = self.args.read_range[rd]
            if rng is None:
                rng = ""
            else:
                rng = "%s-%s" % (rng[0] if rng[0] else "",
                                 rng[1] if rng[1] else "")
            print("- %s: %s" % (rd.upper(),rng))
        # Get input Fastqs
        if self.args.fastqs:
            fastqs = self.args.fastqs
        else:
            fastqs = self.args.project.fastqs
        # Remove Fastqs not in listed sample names
        if self.args.fastq_attrs:
            fastq_attrs = self.args.fastq_attrs
        else:
            fastq_attrs = AnalysisFastq
        fastqs = [fq for fq in fastqs
                      if fastq_attrs(fq).sample_name in
                      self.args.samples]
        # Check for output directory
        if not os.path.exists(self.args.out_dir):
            print("Making output dir: %s" % self.args.out_dir)
            os.makedirs(self.args.out_dir)
        # Get read ranges
        read_range = { int(r[1:]): self.args.read_range[r]
                       for r in self.args.read_range }
        # Store example command lines
        examples = {}
        # Process Fastqs
        for fq in fastqs:
            # Build path for final Fastq file
            ffq = os.path.join(self.args.out_dir,
                               os.path.basename(fq))
            self.output.fastqs.append(ffq)
            # Remove existing symlinks in case
            # original files have since moved
            if os.path.islink(ffq):
                os.remove(ffq)
            elif os.path.exists(ffq):
                # Don't regenerate existing files
                continue
            # Always symlink to index reads
            if fastq_attrs(fq).is_index_read:
                os.symlink(fq,ffq)
                continue
            # Get read ranges
            read_number = fastq_attrs(fq).read_number
            try:
                rng = read_range[read_number]
            except KeyError:
                rng = None
            if rng is None:
                # Make a symlink
                os.symlink(fq,ffq)
            else:
                # Generate a new Fastq using seqtk
                # Get trimming limits
                if rng[0] and rng[0] > 1:
                    trim_leading = rng[0] - 1
                else:
                    trim_leading = None
                length = rng[1]
                # Build the command
                cmd = ["zcat {fq_in}".format(fq_in=fq)]
                if length:
                    cmd.append("seqtk trimfq -L %d -" % length)
                if trim_leading:
                    cmd.append("seqtk trimfq -b %d -" % trim_leading)
                cmd.append("gzip -{compression} >{fq_out}".\
                           format(fq_out=ffq,compression=2))
                self.add_cmd(
                    "Run SeqTK 'trimfq' on '%s'" % os.path.basename(fq),
                    """
                    {cmd}
                    """.format(cmd=' | '.join(cmd)))
                if read_number not in examples:
                    examples[read_number] = ' | '.join(cmd)
        # Print examples
        if examples:
            print("Example commands:")
            for rd in sorted(list(examples.keys())):
                print("- %s" % examples[rd ])

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
        # Strip leading paths from Fastqs
        fastqs = ([os.path.basename(fq) for fq in metadata['fastqs']]
                   if 'fastqs' in metadata else None)
        if fastqs:
            # Collapse list into a string
            metadata['fastqs'] = ','.join(fastqs)
        else:
            metadata['fastqs'] = None
        # Deal with sequence data (biological) samples
        seq_data_samples = (metadata['seq_data_samples']
                            if 'seq_data_samples' in metadata else None)
        if seq_data_samples:
            # Collapse list into a string
            metadata['seq_data_samples'] = ','.join([s for s
                                                     in seq_data_samples])
        else:
            metadata['seq_data_samples'] = None
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

class VerifyFastqs(PipelineFunctionTask):
    """
    Check Fastqs are valid
    """
    def init(self,project,fastq_attrs=None):
        """
        Initialise the VerifyFastqs task

        Arguments:
          project (AnalysisProject): project with Fastqs
            to check
          fastq_attrs (BaseFastqAttrs): class to use for
            extracting data from Fastq names
        """
        pass
    def setup(self):
        # Remove index Fastqs
        fastqs = remove_index_fastqs(
            self.args.project.fastqs,
            fastq_attrs=self.args.fastq_attrs)
        # Check each Fastq is readable
        for fq in fastqs:
            self.add_call(
                "Check %s can be read" % os.path.basename(fq),
                self.read_fastq,
                fq)
    def read_fastq(self,fastq):
        # Iterate through Fastq file
        fq = os.path.basename(fastq)
        try:
            for r in FastqIterator(fastq_file=fastq):
                continue
        except Exception as ex:
            print("%s...FAILED: '%s'" % (fq,ex))
            return False
        print("%s...PASSED" % fq)
        return True
    def finish(self):
        # Report the output (will contain any errors)
        for line in self.stdout.split('\n'):
            if not line.startswith("#### "):
                print(line)
        # Check the verification status
        verified = all(r for r in self.result())
        if not verified:
            self.fail(message="Failed to verify Fastq files")
        else:
            print("Verified Fastq files")

class SetCellCountFromCellranger(PipelineTask):
    """
    Update the number of cells in the project metadata from
    'cellranger count' or 'cellranger multi' output
    """
    def init(self,project,qc_dir=None,source="count"):
        """
        Initialise the SetCellCountFromCellranger task.

        Arguments:
          project (AnalysisProject): project to update the
            number of cells for
          qc_dir (str): directory for QC outputs (defaults
            to subdirectory 'qc' of project directory)
          source (str): either 'count' (the default) or
            'multi'
        """
        pass
    def setup(self):
        # Extract and store the cell count from the cellranger
        # metric file
        try:
            set_cell_count_for_project(self.args.project.dirn,
                                       self.args.qc_dir,
                                       source=self.args.source)
        except Exception as ex:
            print("Failed to set the cell count: %s" % ex)

class GetReferenceDataset(PipelineTask):
    """
    Acquire reference data for an organism from mapping

    Generic lookup task which attempts to locate the matching
    reference dataset from a mapping/dictionary.
    """
    def init(self,organism,references,force_reference=None):
        """
        Initialise the GetReferenceDataset task

        Arguments:
          organism (str): name of the organism
          references (mapping): mapping with organism names
            as keys and reference datasets as corresponding
            values
          force_reference (str): if specified then return
            the supplied value instead of determining from
            the organism

        Outputs:
          reference_dataset: reference dataset (set to None
            if no dataset could be located)
        """
        self.add_output('reference_dataset',Param(type='str'))
    def setup(self):
        if self.args.force_reference:
            print("Using supplied dataset: %s" %
                  self.args.force_reference)
            self.output.reference_dataset.set(
                self.args.force_reference)
        elif self.args.organism:
            organism = str(self.args.organism).lower()
            try:
                self.output.reference_dataset.set(
                    self.args.references[organism])
                print("%s: located %s" %
                      (organism,
                       self.output.reference_dataset.value))
            except Exception as ex:
                print("Unable to locate reference data for organism '%s'"
                      % organism)
        else:
            print("No organism supplied, unable to look up reference data")

class GetBAMFiles(PipelineFunctionTask):
    """
    Create BAM files from Fastqs using STAR

    Runs STAR to generate BAM files from Fastq files. The
    BAMs are then sorted and indexed using samtools.
    """
    def init(self,fastqs,star_index,out_dir,subset_size=None,
             nthreads=None,reads=None,include_samples=None,
             fastq_attrs=None,verbose=False):
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
          include_samples (list): optional, list of sample
            names to include
          fastq_attrs (IlluminaFastq): optional, class to
            use for extracting information from Fastq file
            names
          verbose (bool): if True then print additional
            information from the task

        Outputs:
          bam_files: list of sorted BAM files
        """
        # Conda dependencies
        self.conda("star=2.7.7a",
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
        # Remove Fastqs not in listed sample names
        if self.args.fastq_attrs:
            fastq_attrs = self.args.fastq_attrs
        else:
            fastq_attrs = AnalysisFastq
        if self.args.include_samples:
            fastqs = [fq for fq in self.args.fastqs
                      if fastq_attrs(fq).sample_name in
                      self.args.include_samples]
        else:
            fastqs = self.args.fastqs
        # Remove index reads and group Fastqs
        fq_pairs = group_fastqs_by_name(
            remove_index_fastqs(fastqs,
                                fastq_attrs=fastq_attrs),
            fastq_attrs=fastq_attrs)
        # Filter reads
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
        get_versions = False
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
                get_versions = True
        # Get versions of STAR and samtools
        if get_versions:
            version_file = os.path.join(self.args.out_dir,
                                        "_versions")
            self.add_call("Get STAR and samtools versions",
                          self.get_versions,
                          version_file)
    def get_versions(self,version_file):
        # Get STAR version
        star_cmd = Command('STAR','--version')
        status = star_cmd.run_subprocess(log="__version_STAR")
        if status != 0:
            raise Exception("STAR returned non-zero exit code: %s"
                            % status)
        star_version = None
        with open("__version_STAR",'rt') as fp:
            for line in fp:
                if line.startswith("STAR_"):
                    # Example: STAR_2.4.2a
                    star_version = '_'.join(line.strip().split('_')[1:])
                    break
        # Get samtools version
        samtools_cmd = Command('samtools','--version')
        status = samtools_cmd.run_subprocess(log="__version_samtools")
        if status != 0:
            raise Exception("samtools returned non-zero exit code: %s"
                            % status)
        samtools_version = None
        with open("__version_samtools",'rt') as fp:
            for line in fp:
                if line.startswith("samtools"):
                    # Example: samtools 1.15.1
                    samtools_version = ' '.join(line.strip().split()[1:])
                    break
        # Write to output file
        with open(version_file,'wt') as fp:
            if star_version:
                fp.write("star\t%s\n" % star_version)
            if samtools_version:
                fp.write("samtools\t%s\n" % samtools_version)
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
        shutil.move(sorted_bam_file,bam_file)
        shutil.move(sorted_bam_file_index,"%s.bai" % bam_file)
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

class ConvertGTFToBed(PipelineTask):
    """
    Convert a GTF file to a BED file using BEDOPS 'gtf2bed'
    """
    def init(self,gtf_in,bed_out):
        """
        Initialise the ConvertGTFToBed task

        Arguments:
          gtf_in (str): path to the input GTF file
          bed_out (str): path to the output BED file
        """
        # Conda dependencies
        self.conda("bedops=2.4.41")
        # Outputs
        self.add_output('bed_file',Param())
    def setup(self):
        # Check for output BED
        if os.path.exists(self.args.bed_out):
            print("Output BED file aleady exists")
            return
        # Check for input GTF
        if self.args.gtf_in:
            print("Input GTF file: %s" % self.args.gtf_in)
        else:
            print("Input GTF file not supplied")
            return
        # Set up command to run BEDOPS gtf2bed
        self.add_cmd("Run RSeQC infer_experiment.py",
                     """
                     # Get version
                     echo gtf2bed $(gtf2bed --version 2>&1 | grep version: | cut -d: -f2) >_versions
                     # Convert GTF to BED
                     gtf2bed <{gtf_in} >out.bed
                     if [ $? -ne 0 ] ; then
                       echo "GTF to BED conversion failed"
                       exit 1
                     fi
                     """.format(gtf_in=self.args.gtf_in))
    def finish(self):
        # No output expected
        if not self.args.gtf_in:
            return
        if not os.path.exists(self.args.bed_out):
            # Copy converted BED file to final location
            bed_out = os.path.join(self._working_dir,"out.bed")
            if os.path.exists(bed_out):
                print("Copy BED file to %s" % self.args.bed_out)
                shutil.copy(bed_out,self.args.bed_out)
            else:
                raise Exception("failed to generate BED file")
        # Set output
        self.output.bed_file.set(self.args.bed_out)

class VerifyQC(PipelineFunctionTask):
    """
    Verify outputs from the QC pipeline
    """
    def init(self,project,qc_dir,protocol,fastqs):
        """
        Initialise the VerifyQC task.

        Arguments:
          project (AnalysisProject): project to update the
            number of cells for
          qc_dir (str): directory for QC outputs (defaults
            to subdirectory 'qc' of project directory)
          protocol (QCProtocl): QC protocol to verify against
          fastqs (list): Fastqs to include in the
            verification
        """
        pass
    def setup(self):
        # Call the verification function
        self.add_call(
            "Verify QC outputs for %s" % self.args.project.name,
            verify_project,
            self.args.project,
            self.args.qc_dir,
            qc_protocol=self.args.protocol,
            fastqs=self.args.fastqs)
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

class RunMultiQC(PipelineTask):
    """
    Run MultiQC
    """
    def init(self,project,qc_dir,fastq_dir=None):
        """
        Initialise the RunMultiQC task.

        Arguments:
          project (AnalysisProject): project to generate
            QC report for
          qc_dir (str): directory for QC outputs (defaults
            to subdirectory 'qc' of project directory)
          fastq_dir (str): directory holding Fastq files
            (defaults to current fastq_dir in project)
        """
        self.conda("multiqc=1.24")
        # Specify Python version to use to avoid similar
        # issue as reported here:
        # https://github.com/ewels/MultiQC/issues/1413
        self.conda("python=3.12")
    def setup(self):
        # MultiQC report file
        project = self.args.project
        qc_base = os.path.basename(self.args.qc_dir)
        self.multiqc_report = os.path.join(project.dirn,
                                           "multi%s_report.html" %
                                           qc_base)
        # Report title
        if project.info.run is None:
            title = "%s" % project.name
        else:
            title = "%s/%s" % (project.info.run,
                               project.name)
        if self.args.fastq_dir is not None:
            title = "%s (%s)" % (title,self.args.fastq_dir)
        # Add the command
        self.add_cmd(PipelineCommandWrapper(
            "Run MultiQC",
            "multiqc",
            "--title",title,
            "--filename",self.multiqc_report,
            "--force",
            self.args.qc_dir))
    def finish(self):
        # Check report was generated
        if not os.path.exists(self.multiqc_report):
            self.fail(message="Failed to create MultiQC report: %s" %
                      self.multiqc_report)

class ReportQC(PipelineTask):
    """
    Generate the QC report
    """
    def init(self,project,qc_dir,report_html=None,fastq_dir=None,
             force=False,zip_outputs=True):
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
          force (bool): if True then force HTML report to
            be generated even if QC outputs fail
            verification (default: don't write report)
          zip_outputs (bool): if True then also generate
            a ZIP archive of the QC reports
        """
        pass
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
