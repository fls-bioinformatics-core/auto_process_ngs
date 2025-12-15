#!/usr/bin/env python3
#
#     cellranger_multi: implements 'cellranger_multi' QC module
#     Copyright (C) University of Manchester 2024-2025 Peter Briggs

"""
Implements the 'cellranger_multi' QC module:

* CellrangerMulti: core QCModule class
* GetCellrangerMultiConfig: pipeline task to acquire multi config file
* RunCellrangerMulti: pipeline task to run 'cellranger multi'
* expected_outputs: helper function for handling 'cellranger multi' outputs

Also imports the following pipeline tasks:

- Get10xPackage
- DetermineRequired10xPackage
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
from .cellranger_count import verify_10x_pipeline
from .cellranger_count import DetermineRequired10xPackage
from ..apps.cellranger import CellrangerMulti as CellrangerMultiOutputs
from ..apps.cellranger import extract_path_data
from ..apps.cellranger import fetch_cellranger_multi_output_dirs
from ...bcl2fastq.pipeline import Get10xPackage
from ...command import Command
from ...tenx.cellplex import CellrangerMultiConfigCsv
from ...tenx.utils import add_cellranger_args
from ...pipeliner import PipelineTask
from ...pipeliner import PipelineFunctionTask
from ...pipeliner import PipelineParam as Param
from ...pipeliner import FunctionParam
from ...pipeliner import ListParam
from ...utils import check_required_version
from ...utils import parse_version

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Core class
#######################################################################

class CellrangerMulti(QCModule):
    """
    Class for handling the 'cellranger_multi' QC module
    """
    name = "cellranger_multi"
    mapped_metrics = False
    runners = ("cellranger_multi_runner",)
    envmodules = ("cellranger",)

    def __init__(self):
        QCModule.__init__(self)

    @classmethod
    def collect_qc_outputs(self,qc_dir):
        """
        Collect information on Cellranger multi outputs

        Returns an AttributeDictionary with the following
        attributes:

        - name: set to 'cellranger_multi'
        - software: dictionary of software and versions
        - references: list of associated reference datasets
        - probe_sets: list of associated probe sets
        - fastqs: list of associated Fastq names
        - multiplexed_samples: list of associated multiplexed
          sample names
        - pipelines: list of tuples defining 10x pipelines
          in the form (name,version,reference)
        - samples_by_pipeline: dictionary with lists of
          multiplexed sample names associated with each 10x
          pipeline tuple
        - config_files: list of associated config files
          ('10x_multi_config[.<SAMPLE>].csv')
        - output_files: list of associated output files
        - tags: list of associated output classes

        Arguments:
          qc_dir (QCDir): QC directory to examine
        """
        software = {}
        output_files = list()
        multiplexed_samples = set()
        physical_samples = set()
        cellranger_references = set()
        cellranger_probe_sets = set()
        samples_by_pipeline = dict()
        tags = set()
        # Look for cellranger multi configs
        config_files = list(
            filter(lambda f:
                   f.endswith(".csv") and f.startswith("10x_multi_config."),
                   [os.path.basename(f) for f in qc_dir.file_list]))
        # Look for cellranger multi outputs
        cellranger_multi_dir = os.path.join(qc_dir.path,
                                            "cellranger_multi")
        ##print("Checking for cellranger multi outputs under %s" %
        ##      cellranger_multi_dir)
        versions = set()
        cellranger_name = None
        cellranger_multi_samples = {}
        for multi_output_dir in fetch_cellranger_multi_output_dirs(
                cellranger_multi_dir):
            # Extract version, reference and physical sample
            # data from intermediate dir names
            version, refdata, psample = extract_path_data(
                multi_output_dir,
                cellranger_multi_dir)
            # Store version
            versions.add(version)
            # Store physical sample
            if psample is not None:
                physical_samples.add(psample)
            # Extend data structures
            if version not in cellranger_multi_samples:
                cellranger_multi_samples[version] = {}
            if refdata not in cellranger_multi_samples[version]:
                cellranger_multi_samples[version][refdata] = []
            # Load data
            cellranger_multi = CellrangerMultiOutputs(multi_output_dir)
            # Try to set software name
            cellranger_name = cellranger_multi.pipeline_name
            if cellranger_name is None:
                # Default if name can't be explicitly extracted
                cellranger_name = "cellranger"
            # Reference data
            if cellranger_multi.reference_data:
                cellranger_references.add(
                    cellranger_multi.reference_data)
            # Probeset
            if cellranger_multi.probe_set:
                cellranger_probe_sets.add(cellranger_multi.probe_set)
            # Loop over multiplexed samples
            for smpl in cellranger_multi.sample_names:
                # Add sample
                multiplexed_samples.add(smpl)
                # Add outputs
                cellranger_multi_samples[version][refdata].append(smpl)
                output_files.extend(
                    [cellranger_multi.web_summary(smpl),
                     cellranger_multi.metrics_csv(smpl)])
                # Store sample lists associated with pipeline,
                # version and reference data
                pipeline_key = (cellranger_name, version, refdata)
                samples_by_pipeline[pipeline_key] = \
                    [s for s in
                     cellranger_multi_samples[version][refdata]]
        # If any outputs were found then add tag
        if multiplexed_samples:
            tags.add("cellranger_multi")
        # Store cellranger versions
        if cellranger_name and versions:
            if cellranger_name not in software:
                software[cellranger_name] = list(versions)
            else:
                for version in list(versions):
                    if version not in software[cellranger_name]:
                        software[cellranger_name].append(version)
            software[cellranger_name] = sorted(software[cellranger_name])
        # Return collected information
        # - dictionary of 'software' (keys are cellranger names,
        #   values are sorted lists of matching versions)
        # - list of references across all outputs
        # - list of probesets across all outputs
        # - list of unique multiplexed samplenames across all outputs
        # - list of "pipelines" (a "pipeline" is defined by a tuple
        #   of (cellranger_name, version, refdata))
        # - list of (multiplexed) samples for each pipeline
        # - list of output files (i.e. web_summary and metrics
        #   across all outputs)
        # - list of config files
        return AttributeDictionary(
            name='cellranger_multi',
            software=software,
            references=sorted(list(cellranger_references)),
            probe_sets=sorted(list(cellranger_probe_sets)),
            fastqs=[],
            physical_samples=sorted(list(physical_samples)),
            multiplexed_samples=sorted(list(multiplexed_samples)),
            pipelines=sorted([p for p in samples_by_pipeline]),
            samples_by_pipeline=samples_by_pipeline,
            output_files=output_files,
            config_files=sorted(config_files),
            tags=sorted(list(tags))
        )

    @classmethod
    def verify(self,params,qc_outputs):
        """
        Verify 'cellranger_multi' QC module against outputs

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
        # Check for config files
        config_files = [os.path.join(params.qc_dir, f)
                        for f in os.listdir(params.qc_dir)
                        if f.startswith("10x_multi_config.")
                        and f.endswith(".csv")]
        if not config_files:
            # No multi config files so no outputs expected
            return None
        # Verification paramters
        cellranger_version = params.cellranger_version
        cellranger_refdata = params.cellranger_refdata
        if cellranger_version is None and cellranger_refdata is None:
            # If no version and no refdata then verification is
            # not possible
            return None
        # Update parameters
        if cellranger_version == "*":
            cellranger_version = None
        if cellranger_refdata == "*":
            cellranger_refdata = None
        if cellranger_refdata:
            cellranger_refdata = os.path.basename(cellranger_refdata)
        # Locate and filter cellranger multi output directories based on
        # version and refdata specification
        cellranger_multi_dir = os.path.join(params.qc_dir, "cellranger_multi")
        multi_dirs = {}
        for multi_dir in fetch_cellranger_multi_output_dirs(
                cellranger_multi_dir):
            version, refdata, psample = extract_path_data(
                multi_dir,
                cellranger_multi_dir)
            if (not cellranger_version or version == cellranger_version) and \
               (not cellranger_refdata or refdata == cellranger_refdata):
                try:
                    multi_dirs[psample].append(multi_dir)
                except KeyError:
                    multi_dirs[psample] = [multi_dir]
        # Verify against each config file
        for cf_file in config_files:
            config = CellrangerMultiConfigCsv(cf_file)
            psample = config.physical_sample
            multiplexed_samples = config.sample_names
            if psample not in multi_dirs:
                # Verification failure (no matching outputs)
                return False
            # Look for at least one valid output for the config
            verified_config = False
            for multi_dir in multi_dirs[psample]:
                outputs = CellrangerMultiOutputs(multi_dir)
                output_samples = outputs.sample_names
                missing_sample = False
                if multiplexed_samples:
                    # Check that output exists for each of the
                    # multiplexed samples in the config
                    for s in multiplexed_samples:
                        if s not in output_samples:
                            missing_sample = True
                            break
                    if not missing_sample:
                        verified_config = True
                        break
                else:
                    # No multiplexed samples - check
                    # for a single output
                    if len(output_samples) == 1:
                        verified_config = True
                        break
            if not verified_config:
                # Verification failure (no outputs with
                # all multiplexed samples)
                return False
        # No verification failures
        return True

    @classmethod
    def add_to_pipeline(self,p,project_name,project,qc_dir,
                        qc_module_name,cellranger_exe=None,
                        cellranger_out_dir=None,cellranger_jobmode=None,
                        cellranger_maxjobs=None,
                        cellranger_mempercore=None,
                        cellranger_jobinterval=None,
                        cellranger_localcores=None,
                        cellranger_localmem=None,
                        cellranger_required_version=None,
                        required_tasks=None,
                        cellranger_runner=None,
                        envmodules=None,working_dir=None):
        """
        Adds tasks for 'cellranger_multi' module to pipeline

        Arguments:
        p (Pipeline): pipeline to extend
        project_name (str): name to associate with project for
          reporting tasks
        project (AnalysisProject): project to run 10x
          cellranger pipeline within
        qc_dir (str): directory for QC outputs (defaults
          to subdirectory 'qc' of project directory)
        qc_module_name (str): QC module being used
        cellranger_exe (str): optional, explicitly specify
          the cellranger executable to use (default: cellranger
          executable is determined automatically)
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
        cellranger_runner (JobRunner): runner to use for
          running 'cellranger multi'
        envmodules (list): environment module names to
          load for running Cellranger
        working_dir (str): explicitly specify path to working
          directory
        """
        # Tasks 'check_cellranger_multi' depends on
        check_cellranger_multi_requires = []

        # Locate cellranger
        required_cellranger = DetermineRequired10xPackage(
            "%s: determine required 'cellranger' package" %
            project_name,
            qc_module_name,
            cellranger_exe)
        p.add_task(required_cellranger)

        get_cellranger = Get10xPackage(
            "%s: get information on cellranger" % project_name,
            require_package=\
            required_cellranger.output.require_cellranger)
        p.add_task(get_cellranger,
                   requires=(required_cellranger,),
                   envmodules=envmodules)
        check_cellranger_multi_requires.append(get_cellranger)

        # Locate config.csv files for 'cellranger multi'
        get_cellranger_multi_configs = GetCellrangerMultiConfigs(
            "%s: get config files for 'cellranger multi'" %
            project_name,
            project,
            qc_dir
        )
        p.add_task(get_cellranger_multi_configs,
                   requires=required_tasks)
        check_cellranger_multi_requires.append(
            get_cellranger_multi_configs)

        # Parent directory for cellranger multi outputs
        # Set to project directory unless the 'cellranger_out_dir'
        # parameter is set
        cellranger_out_dir = FunctionParam(
            lambda out_dir,project_dir:
            out_dir if out_dir is not None else project_dir,
            cellranger_out_dir,
            project.dirn)

        # Run cellranger multi
        run_cellranger_multi = RunCellrangerMulti(
            "%s: run cellranger multi (%s %s)" %
            (project_name,
             project.info.single_cell_platform,
             project.info.library_type),
            project,
            get_cellranger_multi_configs.output.config_csvs,
            get_cellranger_multi_configs.output.samples,
            get_cellranger_multi_configs.output.reference_data_path,
            get_cellranger_multi_configs.output.probe_set_path,
            cellranger_out_dir,
            qc_dir=qc_dir,
            working_dir=working_dir,
            cellranger_exe=get_cellranger.output.package_exe,
            cellranger_version=get_cellranger.output.package_version,
            cellranger_jobmode=cellranger_jobmode,
            cellranger_maxjobs=cellranger_maxjobs,
            cellranger_mempercore=cellranger_mempercore,
            cellranger_jobinterval=cellranger_jobinterval,
            cellranger_localcores=cellranger_localcores,
            cellranger_localmem=cellranger_localmem,
            cellranger_required_version=cellranger_required_version
        )
        p.add_task(run_cellranger_multi,
                   requires=(get_cellranger,
                             get_cellranger_multi_configs,),
                   runner=cellranger_runner,
                   envmodules=envmodules)
        return run_cellranger_multi

#######################################################################
# Pipeline tasks
#######################################################################

class GetCellrangerMultiConfigs(PipelineFunctionTask):
    """
    Locate 'config.csv' files for cellranger multi
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
        self.add_output('config_csvs',ListParam())
        self.add_output('samples',ListParam())
        self.add_output('gex_libraries',ListParam())
        self.add_output('fastq_dirs',dict())
        self.add_output('reference_data_path',Param(type=str))
        self.add_output('probe_set_path',Param(type=str))
    def setup(self):
        # Check for top-level multi config files
        config_files = sorted([os.path.join(self.args.project.dirn,f)
                               for f in os.listdir(self.args.project.dirn)
                               if f.startswith("10x_multi_config.")
                               and f.endswith(".csv")])
        if not config_files:
            # No configs found
            print("No 10x multi config files found: nothing to do")
            return
        samples = list()
        gex_libraries = list()
        fastq_dirs = dict()
        for config_file in config_files:
            # Extract information from each config.csv file
            print("Reading '%s'" % os.path.basename(config_file))
            config_csv = CellrangerMultiConfigCsv(config_file,
                                                  strict=False)
            if not config_csv.is_valid:
                print("Errors found in config.csv file:")
                for err in config_csv.get_errors():
                    print(f"- {err}")
                self.fail(message=f"problems with 10x multi config "
                          f"file {config_file}")
                return
            reference_data_path = config_csv.reference_data_path
            probe_set_path = config_csv.probe_set_path
            samples.extend(config_csv.sample_names)
            gex_libraries.extend(config_csv.gex_libraries)
            for sample in config_csv.fastq_dirs:
                fastq_dirs[sample] = config_csv.fastq_dirs[sample]
        print("Samples:")
        for sample in samples:
            print("- %s" % sample)
        print("GEX libraries:")
        for library in gex_libraries:
            print("- %s" % library)
        print("Reference dataset: %s" % reference_data_path)
        print("Probe set        : %s" % probe_set_path)
        # Copy config files to QC dir
        for config_file in config_files:
            print("Copy '%s' into %s" % (config_file,self.args.qc_dir))
            shutil.copy(config_file,self.args.qc_dir)
        # Set outputs
        self.output.config_csvs.extend(
            [os.path.join(self.args.qc_dir,os.path.basename(cf))
             for cf in config_files])
        self.output.samples.extend(samples)
        self.output.gex_libraries.extend(gex_libraries)
        self.output.reference_data_path.set(reference_data_path)
        self.output.probe_set_path.set(probe_set_path)
        for sample in fastq_dirs:
            self.output.fastq_dirs[sample] = fastq_dirs[sample]

class RunCellrangerMulti(PipelineTask):
    """
    Run 'cellranger multi'
    """
    def init(self,project,config_csvs,samples,reference_data_path,
             probe_set_path,out_dir,qc_dir=None,cellranger_exe=None,
             cellranger_version=None,cellranger_jobmode='local',
             cellranger_maxjobs=None,cellranger_mempercore=None,
             cellranger_jobinterval=None,cellranger_localcores=None,
             cellranger_localmem=None,cellranger_required_version=None,
             working_dir=None):
        """
        Initialise the RunCellrangerMulti task.

        Arguments:
          project (AnalysisProject): project to run
            QC for
          config_csvs (list): list of paths to
            'cellranger multi' configuration files
          samples (list): list of sample names from the
            config.csv file
          reference_data_path (str): path to the cellranger
            compatible reference dataset from the config.csv
            file
          probe_set_path (str): path to the probe set
            reference dataset from the config.csv file
          out_dir (str): top-level directory to copy all
            final 'multi' outputs into. Outputs won't be
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
          cellranger_required_version (str): string specifying
            the required Cellranger version (default: None)
        """
        # Internal: top-level working directory
        self._working_dir = None
        # Samples from config.csv
        self._samples = []
        # Add outputs
        self.add_output('cellranger_version',Param(type=str))
        self.add_output('cellranger_refdata',Param(type=str))
        self.add_output('cellranger_probeset',Param(type=str))
        self.add_output('cellranger_exe',Param(type=str))
        self.add_output('cellranger_package',Param(type=str))
    def setup(self):
        # Check if there's anything to do
        if not self.args.config_csvs:
            print("No config file: nothing to do")
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
        cellranger_major_version = int(parse_version(cellranger_version)[0])
        # Check required version, if specified
        if self.args.cellranger_required_version is not None:
            required_version = self.args.cellranger_required_version
            if not check_required_version(cellranger_version,
                                          required_version):
                print(f"Cellranger version {cellranger_version}: doesn't "
                      f"meet version requirement ({required_version})")
                # Don't run cellranger multi
                return
            print(f"Cellranger version {cellranger_version} meets version "
                  f"requirements ({required_version})")
        # Top-level working directory
        self._working_dir = self.args.working_dir
        # Top-level output directory
        multi_dir = os.path.abspath(
            os.path.join(self.args.out_dir,
                         "cellranger_multi",
                         cellranger_version,
                         os.path.basename(self.args.reference_data_path)))
        # Check final outputs from cellranger multi for each config file
        config_files = []
        for config_csv in self.args.config_csvs:
            # Load config data
            config = CellrangerMultiConfigCsv(config_csv)
            # Set prefix
            prefix = multi_dir
            if config.physical_sample:
                prefix = os.path.join(prefix, config.physical_sample)
            # Set ID for cellranger multi
            multi_id = self.args.project.name
            if config.physical_sample:
                multi_id = config.physical_sample
            # Check whether outputs already exist
            for f in expected_outputs(config_csv, multi_id, prefix=prefix):
                if not os.path.exists(f):
                    # At least one expected file is missing so run
                    # 'cellranger multi' for this config file
                    config_files.append(config_csv)
                    break
        if not config_files:
            # Nothing to do
            print("Found existing outputs")
            return
        # Run 'cellranger multi' with each config file
        for config_csv in config_files:
            # Create a working directory for this sample
            work_dir = f"tmp.cellranger_multi.{self.args.project.name}"
            config = CellrangerMultiConfigCsv(config_csv)
            if config.physical_sample:
                work_dir = f"{work_dir}.{config.physical_sample}"
            work_dir = os.path.join(self._working_dir, work_dir)
            # Set ID for cellranger multi run
            multi_id = self.args.project.name
            if config.physical_sample:
                multi_id = config.physical_sample
            # Build cellranger command
            cmd = Command(cellranger_exe,
                          "multi",
                          "--id", multi_id,
                          "--csv", config_csv)
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
            # Generate relative paths to expected outputs
            expected_files = expected_outputs(config_csv, multi_id)
            # Destination for outputs
            dest_dir = multi_dir
            if config.physical_sample:
                dest_dir = os.path.join(dest_dir, config.physical_sample)
            # Add command to run "cellranger multi" in the working
            # directory and copy outputs to the final destination
            print("Running %s" % cmd)
            self.add_cmd(
                "Run %s multi" % cellranger_exe,
                """
                # Create working dir
                mkdir -p {work_dir} && cd {work_dir}
                # Run multi analysis
                {cellranger_multi}
                if [ $? -ne 0 ] ; then
                  echo "Multi library analysis failed"
                  exit 1
                fi
                # Check expected outputs
                ls -ltrh
                for f in {top_level_files} ; do
                  if [ ! -e {multi_id}/$f ] ; then
                    echo "Missing top-level file $f"
                    exit 1
                  fi
                done
                for s in {samples} ; do
                  for f in web_summary.html metrics_summary.csv ; do
                    if [ ! -e {multi_id}/outs/per_sample_outs/$s/$f ] ; then
                      echo "$s: missing outs file $f"
                      exit 1
                    fi
                  done
                done
                # Move outputs to final location
                mkdir -p {dest_dir}
                mv {multi_id}/* {dest_dir}
                """.format(work_dir=work_dir,
                           cellranger_multi=str(cmd),
                           multi_id=multi_id,
                           samples=' '.join(config.sample_names),
                           top_level_files=' '.join(expected_files),
                           dest_dir=dest_dir))
    def finish(self):
        # If no config.csv then ignore and return
        if not self.args.config_csvs:
            print("No config file: cellranger multi analysis was skipped")
            return
        # If version doesn't meet requirements then ignore and return
        if self.args.cellranger_required_version:
            if not check_required_version(
                    self.args.cellranger_version,
                    self.args.cellranger_required_version):
                print(f"Cellranger version not suitable: cellranger multi "
                      "analysis was skipped")
                return
        # Location of final outputs from cellranger multi
        multi_dir = os.path.abspath(
            os.path.join(self.args.out_dir,
                         "cellranger_multi",
                         self.args.cellranger_version,
                         os.path.basename(
                             self.args.reference_data_path)
            ))
        # Top-level destination for QC outputs
        qc_multi_dir = os.path.abspath(
            os.path.join(self.args.qc_dir,
                         "cellranger_multi",
                         self.args.cellranger_version,
                         os.path.basename(
                             self.args.reference_data_path)))
        # Copy subset of outputs to QC directory for each config file
        if self.args.qc_dir:
            for config_csv in self.args.config_csvs:
                print(f"Copying outputs from {os.path.basename(config_csv)} "
                      f"to QC directory")
                # Load data for config file
                config = CellrangerMultiConfigCsv(config_csv)
                # Determine ID used for cellranger multi
                multi_id = self.args.project.name
                if config.physical_sample:
                    multi_id = config.physical_sample
                # Build list of QC outputs to copy
                qc_files = [f for f in expected_outputs(config_csv, multi_id)]
                qc_files.append(os.path.join("outs",
                                             "config.csv"))
                qc_files.append(os.path.join("outs",
                                             "multi",
                                             "multiplexing_analysis",
                                             "tag_calls_summary.csv"))
                # Prefix dirs
                prefix = multi_dir
                if config.physical_sample:
                    prefix = os.path.join(prefix, config.physical_sample)
                qc_prefix = qc_multi_dir
                if config.physical_sample:
                    qc_prefix = os.path.join(qc_prefix,
                                             config.physical_sample)
                # Copy the files
                for path in qc_files:
                    src = os.path.join(prefix, path)
                    if not os.path.exists(src):
                        # File not found, note and skip
                        print("INFO: no file '%s' found" % path)
                        continue
                    # Determine destination and copy
                    dst = os.path.normpath(os.path.join(qc_prefix,
                                                        os.path.dirname(path)))
                    print("Copying file '%s'" % path)
                    if not os.path.exists(dst):
                        mkdirs(dst)
                    shutil.copy(src, dst)
        # Set outputs
        self.output.cellranger_exe.set(self.args.cellranger_exe)
        self.output.cellranger_refdata.set(self.args.reference_data_path)
        self.output.cellranger_probeset.set(self.args.probe_set_path)
        self.output.cellranger_version.set(self.args.cellranger_version)
        self.output.cellranger_package.set(os.path.basename(self.args.cellranger_exe))

#######################################################################
# Helper functions
#######################################################################

def expected_outputs(config_csv, multi_id=None, prefix=None):
    """
    Generate expected output file paths from 10x multi config

    Arguments:
      config_csv (str): path to the 10x multi config file
        to generate the output file names for
      multi_id (str): optional, the ID of the multi run
        (supplied via the --id argument), used if the config
        file doesn't define multiplexed samples
      prefix (str): optional path to prepend to the
        expected file paths

    Returns:
      List: list of paths to expected output files.
    """
    # Load config file
    config = CellrangerMultiConfigCsv(config_csv)
    # Generate list of expected files
    expected_files = ["_cmdline",]
    sample_names = config.sample_names
    if not sample_names:
        # No multiplexed samples
        if multi_id:
            sample_names = [multi_id]
    for sample in sample_names:
        for f in ("web_summary.html",
                  "metrics_summary.csv"):
            # Per-sample outputs
            expected_files.append(os.path.join("outs",
                                               "per_sample_outs",
                                               sample,
                                               f))
    # Prepend prefix if supplied
    if prefix:
        expected_files = [os.path.join(prefix, f)
                          for f in expected_files]
    return expected_files
