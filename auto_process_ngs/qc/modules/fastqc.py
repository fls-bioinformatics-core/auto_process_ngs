#!/usr/bin/env python3
#
#     fastqc: implements 'fastqc' QC module
#     Copyright (C) University of Manchester 2024 Peter Briggs

"""
Implements the 'fastqc' QC module:

* Fastqc: core QCModule class
* CheckFastqcOutputs: pipeline task to check outputs
* RunFastqc: pipeline task to run Fastqc
* check_fastqc_outputs: helper function for checking outputs
"""

#######################################################################
# Imports
#######################################################################

import os
import logging
from bcftbx.utils import AttributeDictionary
from . import QCModule
from ..fastqc import Fastqc as FastqcOutput
from ..fastqc import fastqc_output_files
from ..utils import filter_fastqs
from ...fastq_utils import remove_index_fastqs
from ...pipeliner import PipelineTask
from ...pipeliner import PipelineFunctionTask

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Core class
#######################################################################

class Fastqc(QCModule):
    """
    Class for handling the 'fastqc' QC module
    """
    name = "fastqc"
    mapped_metrics = False
    runners = ("verify_runner",
               "fastqc_runner")
    envmodules = ("fastqc",)
    
    def __init__(self):
        QCModule.__init__(self)

    @classmethod
    def collect_qc_outputs(self,qc_dir):
        """
        Collect information on FastQC outputs

        Returns an AttributeDictionary with the following
        attributes:

        - name: set to 'fastqc'
        - software: dictionary of software and versions
        - fastqs: list of associated Fastq names
        - output_files: list of associated output files
        - tags: list of associated output classes

        Arguments:
          qc_dir (QCDir): QC directory to examine
        """
        versions = set()
        output_files = list()
        fastqs = set()
        tags = set()
        # Look for fastqc outputs
        fastqcs = list(
            filter(lambda f:
                   f.endswith("_fastqc") and
                   os.path.exists("%s.html" % f) and
                   os.path.exists(os.path.join(f,"summary.txt")) and
                   os.path.exists(os.path.join(f,"fastqc_data.txt")),
                   qc_dir.file_list))
        logger.debug("Fastqc: %s" % fastqcs)
        print("\t- %d fastqc files" % len(fastqcs))
        if fastqcs:
            versions = set()
            # Pull out the Fastq names from the Fastqc files
            for fastqc in fastqcs:
                f = os.path.basename(fastqc)[:-len("_fastqc")]
                fastqs.add(f)
                fq = qc_dir.fastq_attrs(f)
                tags.add("fastqc_%s%s" %
                         (('i' if fq.is_index_read else 'r'),
                          (fq.read_number
                           if fq.read_number is not None else '1')))
                versions.add(FastqcOutput(fastqc).version)
            # Store the fastqc files needed for reporting
            output_files.extend(["%s.html" % f for f in fastqcs])
            output_files.extend([os.path.join(f,"summary.txt")
                                 for f in fastqcs])
            output_files.extend([os.path.join(f,"fastqc_data.txt")
                                 for f in fastqcs])
            # Fastqc plot images
            for png in ("per_base_quality",):
                output_files.extend([os.path.join(f,
                                                  "Images",
                                                  "%s.png" % png)
                                     for f in fastqcs])
        # Return collected information
        if versions:
            software = { 'fastqc': sorted(list(versions)) }
        else:
            software = {}
        return AttributeDictionary(
            name=self.name,
            software=software,
            fastqs=sorted(list(fastqs)),
            output_files=output_files,
            tags=sorted(list(tags))
        )

    @classmethod
    def verify(self,params,qc_outputs):
        """
        Verify 'fastqc' QC module against outputs

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
        if not params.fastqs:
            # Nothing to check
            return None
        try:
            # Filter Fastq names
            fastqs = filter_fastqs(params.qc_reads,
                                   params.fastqs)
            # Check that outputs exist for every Fastq
            for fq in fastqs:
                if fq not in qc_outputs.fastqs:
                    return False
            return True
        except KeyError:
            # No Fastqc outputs present
            return False

    @classmethod
    def add_to_pipeline(self,p,project_name,project,qc_dir,
                        read_numbers,fastqs,verbose=True,
                        nthreads=None,require_tasks=[],
                        verify_runner=None,compute_runner=None,
                        envmodules=None):
        """
        Adds tasks for 'fastqc' module to pipeline

        Arguments:
          p (Pipeline): pipeline to extend
          project_name (str): name of project
          project (AnalysisProject): project to run module on
          qc_dir (str): path to QC directory
          read_numbers (list): read numbers to include
          fastqs (list): Fastqs to run the module on
          verbose (bool): enable verbose output
          nthreads (int): number of threads (if not set then
            will be taken from the runner)
          require_tasks (list): list of tasks that the module
            needs to wait for
          verify_runner (JobRunner): runner to use for checks
          compute_runner (JobRunner): runner to use for
            computation
        """
        # Check outputs for FastQC
        check_fastqc = CheckFastQCOutputs(
            "%s: check FastQC outputs" %
            project_name,
            project,
            qc_dir,
            read_numbers=read_numbers,
            fastqs=fastqs,
            verbose=verbose
        )
        p.add_task(check_fastqc,
                   requires=require_tasks,
                   runner=verify_runner)
        # Run FastqQC
        run_fastqc = RunFastQC(
            "%s: FastQC" % project_name,
            check_fastqc.output.fastqs,
            qc_dir,
            nthreads=nthreads
        )
        p.add_task(run_fastqc,
                   requires=(check_fastqc,),
                   runner=compute_runner,
                   envmodules=envmodules)
        return run_fastqc

#######################################################################
# Pipeline tasks
#######################################################################

class CheckFastQCOutputs(PipelineFunctionTask):
    """
    Check the outputs from FastQC
    """
    def init(self,project,qc_dir,read_numbers,fastqs=None,
             verbose=False):
        """
        Initialise the CheckFastQCOutputs task.

        Arguments:
          project (AnalysisProject): project to run
            QC for
          qc_dir (str): directory for QC outputs (defaults
            to subdirectory 'qc' of project directory)
          read_numbers (list): list of read numbers to
            include
          fastqs (list): optional, list of Fastq files
            (overrides Fastqs in project)
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
                      fastqs=self.args.fastqs,
                      read_numbers=self.args.read_numbers)
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
        self.conda("fastqc=0.12.1")
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

#######################################################################
# Helper functions
#######################################################################

def check_fastqc_outputs(project,qc_dir,fastqs=None,
                         read_numbers=None):
    """
    Return Fastqs missing QC outputs from FastQC

    Returns a list of the Fastqs from a project for which
    one or more associated outputs from FastQC don't exist
    in the specified QC directory.

    Arguments:
      project (AnalysisProject): project to check the
        QC outputs for
      qc_dir (str): path to the QC directory (relative
        path is assumed to be a subdirectory of the
        project)
      fastqs (list): optional list of Fastqs to check
        against (defaults to Fastqs from the project)
      read_numbers (list): read numbers to predict
        outputs for

    Returns:
      List: list of Fastq files with missing outputs.
    """
    if not os.path.isabs(qc_dir):
        qc_dir = os.path.join(project.dirn,qc_dir)
    if not fastqs:
        fastqs_in = project.fastqs
    else:
        fastqs_in = fastqs
    fastqs = set()
    for fastq in remove_index_fastqs(fastqs_in,
                                     project.fastq_attrs):
        if read_numbers and \
           project.fastq_attrs(fastq).read_number not in read_numbers:
            # Ignore non-QC reads
            continue
        # FastQC outputs
        for output in [os.path.join(qc_dir,f)
                       for f in fastqc_output_files(fastq)]:
            if not os.path.exists(output):
                fastqs.add(fastq)
    return sorted(list(fastqs))
