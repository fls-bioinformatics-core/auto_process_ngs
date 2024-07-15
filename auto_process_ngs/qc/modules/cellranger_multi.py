#!/usr/bin/env python3
#
#     cellranger_multi: implements 'cellranger_multi' QC module
#     Copyright (C) University of Manchester 2024 Peter Briggs

"""
Implements the 'cellranger_multi' QC module:

* CellrangerMulti: core QCModule class
* GetCellrangerMultiConfig: pipeline task to acquire multi config file
* RunCellrangerMulti: pipeline task to run 'cellranger multi'

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
from ..cellranger import CellrangerMulti as CellrangerMultiOutputs
from ...bcl2fastq.pipeline import Get10xPackage
from ...command import Command
from ...tenx.cellplex import CellrangerMultiConfigCsv
from ...tenx.utils import add_cellranger_args
from ...pipeliner import PipelineTask
from ...pipeliner import PipelineFunctionTask
from ...pipeliner import PipelineParam as Param
from ...pipeliner import FunctionParam
from ...pipeliner import ListParam

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
        - output_files: list of associated output files
        - tags: list of associated output classes

        Arguments:
          qc_dir (QCDir): QC directory to examine
        """
        software = {}
        output_files = list()
        multiplexed_samples = set()
        cellranger_references = set()
        cellranger_probe_sets = set()
        samples_by_pipeline = dict()
        tags = set()
        # Look for cellranger multi outputs
        cellranger_multi_dir = os.path.join(qc_dir.path,
                                            "cellranger_multi")
        ##print("Checking for cellranger multi outputs under %s" %
        ##      cellranger_multi_dir)
        if os.path.isdir(cellranger_multi_dir):
            cellranger_name = None
            versions = set()
            cellranger_multi_samples = {}
            for ver in filter(
                    lambda f:
                    os.path.isdir(os.path.join(cellranger_multi_dir,f)),
                    os.listdir(cellranger_multi_dir)):
                cellranger_multi_samples[ver] = {}
                for ref in filter(
                        lambda f:
                        os.path.isdir(os.path.join(cellranger_multi_dir,ver,f)),
                        os.listdir(os.path.join(cellranger_multi_dir,ver))):
                    # Check putative reference dataset names
                    cellranger_multi_samples[ver][ref] = []
                    cellranger_multi = CellrangerMultiOutputs(
                        os.path.join(
                            cellranger_multi_dir,
                            ver,
                            ref))
                    for smpl in cellranger_multi.sample_names:
                        cellranger_multi_samples[ver][ref].append(smpl)
                        try:
                            output_files.append(cellranger_multi.
                                                web_summary(smpl))
                            output_files.append(cellranger_multi.
                                                metrics_csv(smpl))
                            cellranger_name = cellranger_multi.pipeline_name
                            if cellranger_name is None:
                                cellranger_name = 'cellranger'
                            if cellranger_multi.reference_data:
                                cellranger_references.add(
                                    cellranger_multi.reference_data)
                            if cellranger_multi.probe_set:
                                cellranger_probe_sets.add(
                                    cellranger_multi.probe_set)
                        except OSError:
                            pass
                    # Add outputs, samples and version
                    if cellranger_multi_samples[ver][ref]:
                        tags.add("cellranger_multi")
                        versions.add(ver)
                    for smpl in cellranger_multi_samples[ver][ref]:
                        multiplexed_samples.add(smpl)
            # Store sample lists associated with pipeline,
            # version and reference dataset
            for version in cellranger_multi_samples:
                for reference in cellranger_multi_samples[version]:
                    pipeline_key = (cellranger_name,version,reference)
                    samples_by_pipeline[pipeline_key] = \
                        [s for s in
                         cellranger_multi_samples[version][reference]]
            # Store cellranger versions
            if cellranger_name and versions:
                if cellranger_name not in software:
                    software[cellranger_name] = list(versions)
                else:
                    for v in list(versions):
                        if v not in software[cellranger_name]:
                            software[cellranger_name].append(v)
                software[cellranger_name] = sorted(software[cellranger_name])
        # Return collected information
        return AttributeDictionary(
            name='cellranger_multi',
            software=software,
            references=sorted(list(cellranger_references)),
            probe_sets=sorted(list(cellranger_probe_sets)),
            fastqs=[],
            multiplexed_samples=sorted(list(multiplexed_samples)),
            pipelines=sorted([p for p in samples_by_pipeline]),
            samples_by_pipeline=samples_by_pipeline,
            output_files=output_files,
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
        # Check for config file
        cf_file = os.path.join(params.qc_dir,
                               "10x_multi_config.csv")
        if not os.path.exists(cf_file):
            # No multi config file so no outputs expected
            return True
        # Get expected multiplexed sample names
        # from config file
        multiplexed_samples = CellrangerMultiConfigCsv(cf_file).\
                              sample_names
        if not multiplexed_samples:
            # No samples to check outputs for
            return True
        # Check against actual multiplexed samples
        # associated with specified version and dataset
        return verify_10x_pipeline(('cellranger',
                                    params.cellranger_version,
                                    params.cellranger_refdata),
                                   multiplexed_samples,
                                   qc_outputs)

    @classmethod
    def add_to_pipeline(self,p,project_name,project,qc_dir,
                        qc_module_name,cellranger_exe=None,
                        cellranger_out_dir=None,cellranger_jobmode=None,
                        cellranger_maxjobs=None,
                        cellranger_mempercore=None,
                        cellranger_jobinterval=None,
                        cellranger_localcores=None,
                        cellranger_localmem=None,required_tasks=None,
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

        # Locate config.csv file for 'cellranger multi'
        get_cellranger_multi_config = GetCellrangerMultiConfig(
            "%s: get config file for 'cellranger multi'" %
            project_name,
            project,
            qc_dir
        )
        p.add_task(get_cellranger_multi_config,
                   requires=required_tasks)
        check_cellranger_multi_requires.append(
            get_cellranger_multi_config)

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
            "%s: run cellranger multi (%s)" %
            (project_name,
             project.info.library_type),
            project,
            get_cellranger_multi_config.output.config_csvs,
            get_cellranger_multi_config.output.samples,
            get_cellranger_multi_config.output.reference_data_path,
            get_cellranger_multi_config.output.probe_set_path,
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
        )
        p.add_task(run_cellranger_multi,
                   requires=(get_cellranger,
                             get_cellranger_multi_config,),
                   runner=cellranger_runner,
                   envmodules=envmodules)
        return run_cellranger_multi

#######################################################################
# Pipeline tasks
#######################################################################

class GetCellrangerMultiConfig(PipelineFunctionTask):
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
        config_files = [os.path.join(self.args.project.dirn,f)
                        for f in os.listdir(self.args.project.dirn)
                        if f.startswith("10x_multi_config.")
                        and f.endswith(".csv")]
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
            config_csv = CellrangerMultiConfigCsv(config_file)
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
             cellranger_localmem=None,working_dir=None):
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
        """
        # Internal: top-level working directory
        self._working_dir = None
        # Samples from config.csv
        self._samples = []
        # Whether to run cellranger multi
        self.run_cellranger_multi = False
        # Add outputs
        self.add_output('cellranger_version',Param(type=str))
        self.add_output('cellranger_refdata',Param(type=str))
        self.add_output('cellranger_probeset',Param(type=str))
        self.add_output('cellranger_exe',Param(type=str))
    def setup(self):
        # Check if there's anything to do
        if not self.args.config_csvs:
            print("No config file: nothing to do")
            return
        if len(self.args.config_csvs) > 1:
            print("Too many config files: skipping multi analysis")
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
        self._expected_files = ["_cmdline",]
        for sample in self.args.samples:
            for f in ("web_summary.html",
                      "metrics_summary.csv"):
                # Per-sample outputs
                self._expected_files.append(os.path.join("outs",
                                                         "per_sample_outs",
                                                         sample,
                                                         f))
        # Check outputs and run cellranger if required
        multi_dir = os.path.abspath(
            os.path.join(self.args.out_dir,
                         "cellranger_multi",
                         cellranger_version,
                         os.path.basename(self.args.reference_data_path)))
        for path in self._expected_files:
            if not os.path.exists(os.path.join(multi_dir,path)):
                # At least one expected file is missing
                # so we need to run 'multi'
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
                      "--csv",self.args.config_csvs[0])
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
              if [ ! -e {project}/$f ] ; then
                echo "Missing top-level file $f"
                exit 1
              fi
            done
            for s in {samples} ; do
              for f in web_summary.html metrics_summary.csv ; do
                if [ ! -e {project}/outs/per_sample_outs/$s/$f ] ; then
                  echo "$s: missing outs file $f"
                  exit 1
                fi
              done
            done
            # Move outputs to final location
            mkdir -p {dest_dir}
            mv {project}/* {dest_dir}
            """.format(work_dir=work_dir,
                       cellranger_multi=str(cmd),
                       project=self.args.project.name,
                       samples=' '.join(self.args.samples),
                       top_level_files=' '.join(self._expected_files),
                       dest_dir=multi_dir))
    def finish(self):
        # If no config.csv then ignore and return
        if not self.args.config_csvs:
            print("No config file: cellranger multi analysis was skipped")
            return
        elif len(self.args.config_csvs) > 1:
            print("Too many config files: cellranger multi analysis was "
                  "skipped")
            return
        # Location of final outputs
        multi_dir = os.path.abspath(
            os.path.join(self.args.out_dir,
                         "cellranger_multi",
                         self.args.cellranger_version,
                         os.path.basename(
                             self.args.reference_data_path)
            ))
        # Copy subset of outputs to QC directory
        if self.args.qc_dir:
            print("Copying outputs to QC directory")
            # Build list of QC outputs to copy
            qc_files = [f for f in self._expected_files]
            qc_files.append(os.path.join("outs",
                                         "config.csv"))
            qc_files.append(os.path.join("outs",
                                         "multi",
                                         "multiplexing_analysis",
                                         "tag_calls_summary.csv"))
            # Final destination for QC outputs
            qc_multi_dir = os.path.abspath(
                os.path.join(self.args.qc_dir,
                             "cellranger_multi",
                             self.args.cellranger_version,
                             os.path.basename(
                                 self.args.reference_data_path)))
            # Copy the files
            for path in qc_files:
                src = os.path.join(multi_dir,path)
                if not os.path.exists(src):
                    # File not found, note and skip
                    print("INFO: no file '%s' found" % path)
                    continue
                # Determine destination and copy
                dst = os.path.normpath(os.path.join(qc_multi_dir,
                                                    os.path.dirname(path)))
                print("Copying file '%s'" % path)
                if not os.path.exists(dst):
                    mkdirs(dst)
                shutil.copy(src,dst)
        # Set outputs
        self.output.cellranger_exe.set(self.args.cellranger_exe)
        self.output.cellranger_refdata.set(self.args.reference_data_path)
        self.output.cellranger_probeset.set(self.args.probe_set_path)
        self.output.cellranger_version.set(self.args.cellranger_version)
