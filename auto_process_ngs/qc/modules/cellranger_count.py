#!/usr/bin/env python3
#
#     cellranger_count: implements 'cellranger_count' QC module
#     Copyright (C) University of Manchester 2024 Peter Briggs

"""
Implements the 'cellranger_count' QC module:

* CellrangerCount: core QCModule class

Pipeline task classes:

* DetermineRequired10xPackage
* GetCellrangerReferenceData
* MakeCellrangerArcCountLibraries
* CheckCellrangerCountOutputs
* RunCellrangerCount

Pipeline helper functions:

* add_cellranger_count: adds tasks to pipeline to run 'cellranger* count'
* filter_10x_pipelines: filters a list of 10x pipelines
* verify_10x_pipleine: check for and verify outputs for a 10x package

Also imports the following pipeline tasks:

- Get10xPackage

Additional helper functions:

* check_cellranger_count_outputs: fetch names of samples with
  with missing 'cellranger count' outputs
* check_cellranger_atac_count_outputs: fetch names of samples
  with missing 'cellranger-atac count' outputs
* check_cellranger_arc_count_outputs: fetch names of samples
  with missing 'cellranger-arc count' outputs
"""

#######################################################################
# Imports
#######################################################################

import os
import shutil
import logging
from bcftbx.utils import AttributeDictionary
from bcftbx.utils import mkdirs
from . import QCModule
from ..cellranger import CellrangerCount as CellrangerCountOutputs
from ..cellranger import cellranger_count_output
from ..cellranger import cellranger_atac_count_output
from ..cellranger import cellranger_arc_count_output
from ...bcl2fastq.pipeline import Get10xPackage
from ...command import Command
from ...pipeliner import PipelineTask
from ...pipeliner import PipelineFunctionTask
from ...pipeliner import PipelineParam as Param
from ...pipeliner import FunctionParam
from ...tenx.cellplex import CellrangerMultiConfigCsv
from ...tenx.multiome import MultiomeLibraries
from ...tenx.utils import add_cellranger_args
from ...utils import get_organism_list

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Core class
#######################################################################

class CellrangerCount(QCModule):
    """
    Class for handling the 'cellranger_count' QC module
    """
    name = "cellranger_count"
    mapped_metrics = False
    runners = ("verify_runner",
               "cellranger_count_runner")
    envmodules = ("cellranger",)
    
    def __init__(self):
        QCModule.__init__(self)

    @classmethod
    def collect_qc_outputs(self,qc_dir):
        """
        Collect information on Cellranger count outputs

        Returns an AttributeDictionary with the following
        attributes:

        - name: set to 'cellranger_count'
        - software: dictionary of software and versions
        - references: list of associated reference datasets
        - fastqs: list of associated Fastq names
        - samples: list of associated sample names
        - pipelines: list of tuples defining 10x pipelines
          in the form (name,version,reference)
        - samples_by_pipeline: dictionary with lists of
          sample names associated with each 10x pipeline
          tuple
        - output_files: list of associated output files
        - tags: list of associated output classes

        Arguments:
          qc_dir (QCDir): QC directory to examine
        """
        software = {}
        output_files = list()
        cellranger_samples = []
        cellranger_references = set()
        samples_by_pipeline = dict()
        tags = set()
        # Look for cellranger_count outputs
        cellranger_count_dir = os.path.join(qc_dir.path,
                                            "cellranger_count")
        ##print("Checking for cellranger* count outputs under %s" %
        ##      cellranger_count_dir)
        cellranger_versioned_samples = {}
        if os.path.isdir(cellranger_count_dir):
            cellranger_name = None
            versions = set()
            # Old-style (unversioned)
            for d in filter(
                    lambda f:
                    os.path.isdir(os.path.join(cellranger_count_dir,f)),
                    os.listdir(cellranger_count_dir)):
                sample_dir = os.path.join(cellranger_count_dir,d)
                try:
                    cellranger = CellrangerCountOutputs(sample_dir)
                    output_files.append(cellranger.web_summary)
                    output_files.append(cellranger.metrics_csv)
                    output_files.append(cellranger.cmdline_file)
                    cellranger_samples.append(d)
                    cellranger_name = cellranger.pipeline_name
                    cellranger_references.add(cellranger.reference_data)
                    # Store as version '?'
                    ref = os.path.basename(cellranger.reference_data)
                    if cellranger_name not in cellranger_versioned_samples:
                        cellranger_versioned_samples[cellranger_name] = {}
                    if '?' not in cellranger_versioned_samples[cellranger_name]:
                        cellranger_versioned_samples[cellranger_name]['?'] = {}
                    if ref not in \
                       cellranger_versioned_samples[cellranger_name]['?']:
                        cellranger_versioned_samples[cellranger_name]['?'][ref] = []
                    cellranger_versioned_samples[cellranger_name]['?'][ref].append(d)
                    versions.add('?')
                except OSError:
                    pass
            if cellranger_samples:
                tags.add("%s_count" % cellranger_name)
            # New-style (versioned)
            cellranger_name = None
            for ver in filter(
                    lambda f:
                    os.path.isdir(os.path.join(cellranger_count_dir,f)),
                    os.listdir(cellranger_count_dir)):
                # Check putative version numbers
                for ref in filter(
                        lambda f:
                        os.path.isdir(os.path.join(cellranger_count_dir,ver,f)),
                        os.listdir(os.path.join(cellranger_count_dir,ver))):
                    # Check putative reference dataset names
                    samples = []
                    for smpl in filter(
                            lambda f:
                            os.path.isdir(os.path.join(cellranger_count_dir,
                                                       ver,ref,f)),
                            os.listdir(os.path.join(cellranger_count_dir,
                                                    ver,ref))):
                        sample_dir = os.path.join(cellranger_count_dir,
                                                  ver,ref,smpl)
                        cellranger_name = None
                        try:
                            cellranger = CellrangerCountOutputs(sample_dir)
                            output_files.append(cellranger.web_summary)
                            output_files.append(cellranger.metrics_csv)
                            output_files.append(cellranger.cmdline_file)
                            samples.append(smpl)
                            cellranger_name = cellranger.pipeline_name
                            cellranger_references.add(
                                cellranger.reference_data)
                        except OSError:
                            pass
                    # Add outputs, samples and version
                    if samples:
                        tags.add("%s_count" % cellranger_name)
                        if cellranger_name not in cellranger_versioned_samples:
                            cellranger_versioned_samples[cellranger_name] = {}
                        if ver not in cellranger_versioned_samples[cellranger_name]:
                            cellranger_versioned_samples[cellranger_name][ver] = {}
                        cellranger_versioned_samples[cellranger_name][ver][ref] = samples
                        versions.add(ver)
                        for smpl in cellranger_versioned_samples[cellranger_name][ver][ref]:
                            if smpl not in cellranger_samples:
                                cellranger_samples.append(smpl)
            # Store cellranger versions
            for cellranger_name in cellranger_versioned_samples:
                software[cellranger_name] = sorted(list(cellranger_versioned_samples[cellranger_name].keys()))
        # Store sample lists associated with pipeline,
        # version and reference dataset
        for name in cellranger_versioned_samples:
            for version in cellranger_versioned_samples[name]:
                for reference in cellranger_versioned_samples[name][version]:
                    pipeline_key = (name,version,reference)
                    samples_by_pipeline[pipeline_key] = \
                        [s for s in
                         cellranger_versioned_samples[name][version][reference]]
        # Return collected information
        return AttributeDictionary(
            name=self.name,
            software=software,
            references=sorted(list(cellranger_references)),
            fastqs=[],
            samples=cellranger_samples,
            pipelines=sorted([p for p in samples_by_pipeline]),
            samples_by_pipeline=samples_by_pipeline,
            output_files=output_files,
            tags=sorted(list(tags))
        )

    @classmethod
    def verify(self,params,qc_outputs):
        """
        Verify 'cellranger_count' QC module against outputs

        Returns one of 3 values:

        - True: outputs verified ok
        - False: outputs failed to verify
        - None: verification not possible

        Arguments:
          params (AttributeDictionary): values of parameters
            used as inputs
          qc_outputs (AttributeDictionary): QC outputs returned
            from the 'collect_qc_outputs' method
        """
        if params.cellranger_use_multi_config:
            # Take parameters from 10x_multi_config.csv
            cf_file = os.path.join(params.qc_dir,
                                   "10x_multi_config.csv")
            if not os.path.exists(cf_file):
                # No multi config file so no outputs expected
                return True
            # Get GEX sample names and reference dataset from
            # multi config file
            cf = CellrangerMultiConfigCsv(cf_file)
            samples = cf.gex_libraries
            cellranger_refdata = cf.reference_data_path
        else:
            # Use supplied parameters
            samples = params.samples
            cellranger_refdata = params.cellranger_refdata
        if not samples:
            # No samples so cellranger outputs not expected
            return True
        if cellranger_refdata is None:
            # No reference data so cellranger outputs not expected
            return True
        # Check expected samples against actual samples
        # associated with specified version and dataset
        return verify_10x_pipeline(('cellranger',
                                    params.cellranger_version,
                                    cellranger_refdata),
                                   samples,
                                   qc_outputs)

    @classmethod
    def add_to_pipeline(self,*args,**kws):
        """
        Adds tasks for 'cellranger_count' module to pipeline

        Wrapper for the 'add_cellranger_count' function
        """
        return add_cellranger_count(*args,**kws)

#######################################################################
# Pipeline tasks
#######################################################################

class DetermineRequired10xPackage(PipelineTask):
    """
    Determine which 10xGenomics software package is required

    By default determines the package name based on the
    supplied QC module, but this can be overridden by
    explicitly supplying a required package (which can
    also be a path to an executable).

    The output 'require_cellranger' parameter should be
    supplied to the 'Get10xPackage' task, which will
    do the job of actually locating an executable.
    """
    def init(self,qc_module,require_cellranger=None):
        """
        Initialise the DetermineRequired10xPackage task

        Argument:
          qc_module (str): QC module being used
          require_cellranger (str): optional package name
            or path to an executable; if supplied then
            overrides the automatic package determination

        Outputs:
          require_cellranger (pipelineParam): the 10xGenomics
            software package name or path to use
        """
        self.add_output('require_cellranger',Param(type=str))
    def setup(self):
        qc_modules = {
            "cellranger_count": "cellranger",
            "cellranger-atac_count": "cellranger-atac",
            "cellranger-arc_count": "cellranger-arc",
            "cellranger_multi": "cellranger",
        }
        require_cellranger = self.args.require_cellranger
        if require_cellranger is None:
            try:
                require_cellranger = qc_modules[self.args.qc_module]
            except KeyError:
                raise Exception("Can't identify 10xGenomics package "
                                "required for QC module '%s'" %
                                self.args.qc_module)
        print("Required 10x package: %s" % require_cellranger)
        self.output.require_cellranger.set(require_cellranger)

class GetCellrangerReferenceData(PipelineFunctionTask):
    """
    """
    def init(self,project,organism=None,transcriptomes=None,
             premrna_references=None,atac_references=None,
             multiome_references=None,cellranger_exe=None,
             cellranger_version=None,force_reference_data=None):
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
            if self.args.project.info.library_type == "snRNA-seq" and \
               int(self.args.cellranger_version.split('.')[0]) < 5:
                references = self.args.premrna_references
        elif cellranger_pkg == "cellranger-atac":
            references = self.args.atac_references
        elif cellranger_pkg == "cellranger-arc":
            references = self.args.multiome_references
        else:
            self.fail(message="Don't know which reference "
                      "dataset to use for '%s' and library '%s'" %
                      (self.args.cellranger_exe,
                       self.args.project.info.library_type))
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

class CheckCellrangerCountOutputs(PipelineFunctionTask):
    """
    Check the outputs from cellranger(-atac) count
    """
    def init(self,project,fastq_dir=None,samples=None,qc_dir=None,
             qc_module=None,extra_projects=None,
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
          qc_module (str): QC protocol being used
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
        if self.args.qc_module == "cellranger_count":
            check_outputs = check_cellranger_count_outputs
        elif self.args.qc_module == "cellranger-atac_count":
            check_outputs = check_cellranger_atac_count_outputs
        elif self.args.qc_module == "cellranger-arc_count":
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
    Run 'cellranger* count'
    """
    def init(self,samples,fastq_dir,reference_data_path,library_type,
             out_dir,qc_dir=None,cellranger_exe=None,
             cellranger_version=None,chemistry='auto',fastq_dirs=None,
             force_cells=None,cellranger_jobmode='local',
             cellranger_maxjobs=None,cellranger_mempercore=None,
             cellranger_jobinterval=None,cellranger_localcores=None,
             cellranger_localmem=None):
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
          library_type (str): type of data being analysed
            (e.g. 'scRNA-seq')
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
            count_dir = os.path.abspath(
                os.path.join(self.args.out_dir,
                             "cellranger_count",
                             cellranger_version,
                             os.path.basename(self.args.reference_data_path),
                             sample))
            outs_dir = os.path.join(count_dir,"outs")
            for f in self._outs_files:
                path = os.path.join(outs_dir,f)
                if not os.path.exists(path):
                    run_cellranger_count = True
                    break
            for f in self._top_level_files:
                path = os.path.join(count_dir,f)
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
                    if self.args.library_type == "snRNA-seq":
                        # For single nuclei RNA-seq specify the
                        # --include-introns for cellranger 5.0+
                        if cellranger_major_version in (7,8):
                            cmd.add_args("--include-introns",
                                         "true")
                        else:
                            cmd.add_args("--include-introns")
                # Additional options for cellranger 8.0+
                if cellranger_major_version >= 8:
                    # --create-bam is compulsory
                    # Recommended to set to 'true'
                    cmd.add_args("--create-bam","true")
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
            self.add_cmd("Run %s count for %s" % (cellranger_exe,sample),
                         """
                         # Create working dir
                         mkdir -p {work_dir} && cd {work_dir}
                         # Run single library analysis
                         {cellranger_count}
                         if [ $? -ne 0 ] ; then
                           echo "{sample}: single library analysis failed"
                           exit 1
                         fi
                         # Check expected outputs
                         for f in {top_level_files} ; do
                           if [ ! -e {sample}/$f ] ; then
                             echo "{sample}: missing top-level file $f"
                             exit 1
                           fi
                         done
                         for f in {outs_files} ; do
                           if [ ! -e {sample}/outs/$f ] ; then
                             echo "{sample}: missing outs file $f"
                             exit 1
                           fi
                         done
                         # Move outputs to final location
                         mkdir -p {dest_dir}
                         mv {sample} {dest_dir}
                         """.format(cellranger_count=str(cmd),
                                    outs_files=' '.join(self._outs_files),
                                    top_level_files=' '.join(
                                        self._top_level_files),
                                    sample=sample,
                                    work_dir=work_dir,
                                    dest_dir=os.path.dirname(count_dir)))
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
        # Copy outputs to QC directory
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

#######################################################################
# Helper functions
#######################################################################

def add_cellranger_count(p,project_name,project,qc_dir,
                         organism,fastq_dir,qc_module_name,
                         library_type,chemistry,
                         transcriptome_references,premrna_references,
                         atac_references,multiome_references,
                         force_cells,samples=None,fastq_dirs=None,
                         cellranger_exe=None,reference_dataset=None,
                         extra_projects=None,cellranger_out_dir=None,
                         cellranger_jobmode=None,cellranger_maxjobs=None,
                         cellranger_mempercore=None,
                         cellranger_jobinterval=None,
                         cellranger_localcores=None,
                         cellranger_localmem=None,required_tasks=None,
                         verify_runner=None,cellranger_runner=None,
                         envmodules=None,verbose=False):
    """
    Add tasks to pipeline to run 'cellranger* count'

    Arguments:
      p (Pipeline): pipeline to extend
      project_name (str): name to associate with project for
        reporting tasks
      project (AnalysisProject): project to run 10x
        cellranger pipeline within
      qc_dir (str): directory for QC outputs (defaults
        to subdirectory 'qc' of project directory)
      organism (str): organism for pipeline
      fastq_dir (str): directory holding Fastq files
      qc_module (str): QC module being used
      library_type (str): type of data being analysed (e.g.
        'scRNA-seq')
      chemistry (str): chemistry to use in single
        library analysis
      transcriptome_references (mapping): mapping of
        organism names to reference transcriptome data
        for cellranger
      premrna_references (mapping): mapping of organism
        names to "pre-mRNA" reference data for cellranger
      atac_references (mapping): mapping of
        organism names to ATAC-seq reference genome data
        for cellranger-atac
      multiome_references (mapping): mapping of
        organism names to multiome reference datasets
        for cellranger-arc
      force_cells (int): if set then bypasses
        the cell detection algorithm in 'cellranger'
        and 'cellranger-atac' using the '--force-cells'
        option (does nothing for 'cellranger-arc')
      samples (list): optional, list of samples to
        restrict single library analyses to (or None
        to use all samples in project)
      fastq_dirs (dict): optional, a dictionary mapping
        sample names to Fastq directories which will
        be used to override the paths set by the
        'fastq_dirs' argument
      cellranger_exe (str): optional, explicitly specify
        the cellranger executable to use for single
        library analysis (default: cellranger executable
        is determined automatically)
      reference_dataset (str): optional, path to
        reference dataset (otherwise will be determined
        automatically based on organism)
      extra_projects (list): optional list of extra
        AnalysisProjects to include Fastqs from when
        running cellranger pipeline
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
      required_tasks (list): list of tasks that the
        cellranger pipeline should wait for
      verify_runner (JobRunner): runner to use for checks
      cellranger_runner (JobRunner): runner to use for
        running 'cellranger* count'
      envmodules (list): environment module names to
        load for running Cellranger
      verbose (bool): enable verbose output
    """
    # Tasks 'check_cellranger_count' depends on
    check_cellranger_count_requires = []

    # Locate cellranger
    required_cellranger = DetermineRequired10xPackage(
        "%s: determine required 10x pipeline package (%s)" %
        (project_name,qc_module_name),
        qc_module_name,
        cellranger_exe)
    p.add_task(required_cellranger)

    get_cellranger = Get10xPackage(
        "%s: get information on 10x pipeline package (%s)" %
        (project_name,qc_module_name),
        require_package=\
        required_cellranger.output.require_cellranger)
    p.add_task(get_cellranger,
               requires=(required_cellranger,),
               envmodules=envmodules)
    check_cellranger_count_requires.append(get_cellranger)

    # Get reference data for cellranger
    get_cellranger_reference_data = GetCellrangerReferenceData(
        "%s: get single library analysis reference data (%s)" %
        (project_name,qc_module_name),
        project,
        organism=organism,
        transcriptomes=transcriptome_references,
        premrna_references=premrna_references,
        atac_references=atac_references,
        multiome_references=multiome_references,
        cellranger_exe=get_cellranger.output.package_exe,
        cellranger_version=get_cellranger.output.package_version,
        force_reference_data=reference_dataset
    )
    p.add_task(get_cellranger_reference_data,
               requires=(get_cellranger,),
               runner=verify_runner)
    check_cellranger_count_requires.append(get_cellranger_reference_data)

    # Make libraries.csv files (cellranger-arc only)
    if qc_module_name == "cellranger-arc_count":
        make_cellranger_libraries = MakeCellrangerArcCountLibraries(
            "%s: make libraries files for 'cellranger-arc count'" %
            project_name,
            project,
            qc_dir
        )
        p.add_task(make_cellranger_libraries,
                   requires=required_tasks)
        check_cellranger_count_requires.append(
            make_cellranger_libraries)

    # Check QC outputs for cellranger count
    check_cellranger_count = CheckCellrangerCountOutputs(
        "%s: check for single library analysis outputs (%s)" %
        (project_name,qc_module_name),
        project,
        fastq_dir=fastq_dir,
        samples=samples,
        qc_dir=qc_dir,
        qc_module=qc_module_name,
        extra_projects=extra_projects,
        cellranger_version=get_cellranger.output.package_version,
        cellranger_ref_data=\
        get_cellranger_reference_data.output.reference_data_path,
        verbose=verbose
    )
    p.add_task(check_cellranger_count,
               requires=check_cellranger_count_requires,
               runner=verify_runner)

    # Parent directory for cellranger count outputs
    # Set to project directory unless the 'cellranger_out_dir'
    # parameter is set
    cellranger_out_dir = FunctionParam(
        lambda out_dir,project_dir:
        out_dir if out_dir is not None else project_dir,
        cellranger_out_dir,
        project.dirn)

    # Run cellranger count
    run_cellranger_count = RunCellrangerCount(
        "%s: run single library analysis (%s)" %
        (project_name,project.info.library_type),
        check_cellranger_count.output.samples,
        check_cellranger_count.output.fastq_dir,
        get_cellranger_reference_data.output.reference_data_path,
        library_type=library_type,
        out_dir=cellranger_out_dir,
        qc_dir=qc_dir,
        cellranger_exe=get_cellranger.output.package_exe,
        cellranger_version=get_cellranger.output.package_version,
        chemistry=chemistry,
        fastq_dirs=fastq_dirs,
        force_cells=force_cells,
        cellranger_jobmode=cellranger_jobmode,
        cellranger_maxjobs=cellranger_maxjobs,
        cellranger_mempercore=cellranger_mempercore,
        cellranger_jobinterval=cellranger_jobinterval,
        cellranger_localcores=cellranger_localcores,
        cellranger_localmem=cellranger_localmem
    )
    p.add_task(run_cellranger_count,
               requires=(get_cellranger,
                         get_cellranger_reference_data,
                         check_cellranger_count,),
               runner=cellranger_runner,
               envmodules=envmodules)

    # Return the 'cellranger count' task
    return run_cellranger_count

def filter_10x_pipelines(p,pipelines):
    """
    Filter list of 10x pipelines

    Pipelines are described using tuples of the form:

    (NAME,VERSION,REFERENCE)

    for example:

    ('cellranger','6.1.2','refdata-gex-2020')

    Only pipelines matching the specified name, version
    and reference data will be included in the returned
    list.

    Where the supplied version or reference dataset name
    are either None or '*', these will match any version
    and/or reference dataset.

    Arguments:
      p (tuple): tuple specifying pipeline(s) to match
        against
      pipelines (list): list of pipeline tuples to filter

    Returns:
      List: list of matching 10x pipeline tuples.
    """
    # Extract elements from pipeline pattern
    name = p[0]
    version = p[1]
    refdata = p[2]
    # Check for wildcard versions and reference data
    if version == "*":
        version = None
    if refdata == "*":
        refdata = None
    # Normalise reference dataset name
    refdata = (os.path.basename(refdata) if refdata else None)
    # Find all matching pipelines
    matching_pipelines = list()
    for pipeline in pipelines:
        if pipeline[0] != name:
            # Wrong 10x package name
            continue
        if version and pipeline[1] != version:
            # Wrong version
            continue
        if refdata and pipeline[2] != refdata:
            # Wrong reference dataset
            continue
        # Passed all filters
        matching_pipelines.append(pipeline)
    return matching_pipelines

def verify_10x_pipeline(pipeline,samples,qc_outputs):
    """
    Check for and verify outputs for 10x package

    Arguments:
      pipeline (tuple): tuple specifying pipeline(s) to
        verify
      samples (list): list of sample names to verify
      qc_outputs (AttributeDictionary): QC outputs
        returned from the 'collect_qc_outputs' method

    Returns:
      Boolean: True if at least one set of valid outputs
        exist for the specified pipeline and sample list,
        False otherwise.
    """
    pipelines = filter_10x_pipelines(pipeline,
                                     qc_outputs.pipelines)
    for pipeline in pipelines:
        verified_pipeline = True
        for sample in samples:
            if sample not in qc_outputs.samples_by_pipeline[pipeline]:
                # At least one sample missing outputs from
                # this pipeline, so move on to the next
                verified_pipeline = False
                break
        if verified_pipeline:
            # At least one matching pipeline has
            # verified
            return True
    # No matching outputs from cellranger count
    return False

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
