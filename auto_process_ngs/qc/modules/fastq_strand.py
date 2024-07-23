#!/usr/bin/env python3
#
#     fastq_strand: implements 'fastq_strand' QC module
#     Copyright (C) University of Manchester 2024 Peter Briggs

#FIXME there is some ambiguity elsewhere in the code as to
#FIXME this should be referred to as 'strandedness' or as
#FIXME 'fastq_strand'?

"""
Implements the 'fastq_strand' QC module:

* FastqStrand: core QCModule class
* SetupFastqStrandConf: pipeline task to set up conf file
* CheckFastqStrandOutputs: pipeline task to check outputs
* RunFastqStrand: pipeline task to run 'fastq_strand'
* check_fastq_strand_outputs: helper function for checking output files
"""

#######################################################################
# Imports
#######################################################################

import os
import logging
from bcftbx.utils import AttributeDictionary
from . import QCModule
from ..fastq_strand import Fastqstrand as FastqstrandOutput
from ..fastq_strand import build_fastq_strand_conf
from ..fastq_strand import fastq_strand_output
from ..utils import filter_fastqs
from ...fastq_utils import group_fastqs_by_name
from ...fastq_utils import remove_index_fastqs
from ...pipeliner import PipelineTask
from ...pipeliner import PipelineFunctionTask
#FIXME can we avoid using PipelineCommandWrapper?
from ...pipeliner import PipelineCommandWrapper
from ...pipeliner import PipelineParam as Param
from ...utils import get_organism_list

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Core class
#######################################################################

class FastqStrand(QCModule):
    """
    Class for handling the 'fastq_strand' QC module
    """
    name = "strandedness"
    mapped_metrics = True
    runners = ("star_runner",
               "verify_runner")
    envmodules = ("star")
    
    def __init__(self):
        QCModule.__init__(self)

    @classmethod
    def collect_qc_outputs(self,qc_dir):
        """
        Collect information on FastqStrand outputs

        Returns an AttributeDictionary with the following
        attributes:

        - name: set to 'fastq_strand'
        - software: dictionary of software and versions
        - fastqs: list of associated Fastq names
        - config_files: list of associated config files
          ('fastq_strand.conf')
        - output_files: list of associated output files
        - tags: list of associated output classes

        Arguments:
          qc_dir (QCDir): QC directory to examine
        """
        versions = set()
        output_files = list()
        fastqs = set()
        tags = set()
        # Look for fastq_strand config file
        if os.path.join(qc_dir.path,"fastq_strand.conf") in qc_dir.file_list:
            config_files = ["fastq_strand.conf"]
        else:
            config_files = []
        # Look for fastq_strand outputs
        fastq_strand = list(filter(lambda f:
                                   f.endswith("_fastq_strand.txt"),
                                   qc_dir.file_list))
        logger.debug("fastq_strand: %s" % fastq_strand)
        print("\t- %d fastq_strand files" % len(fastq_strand))
        if fastq_strand:
            tags.add("strandedness")
            for f in fastq_strand:
                fq = qc_dir.fastq_attrs(os.path.splitext(f)[0])
                fastqs.add(
                    os.path.basename(
                        os.path.splitext(f)[0])[:-len("_fastq_strand")])
                versions.add(FastqstrandOutput(f).version)
            # Store the fastq_strand files
            output_files.extend(fastq_strand)
        # Return collected information
        if versions:
            software = { 'fastq_strand': sorted(list(versions)) }
        else:
            software = {}
        return AttributeDictionary(
            name=self.name,
            software=software,
            fastqs=sorted(list(fastqs)),
            output_files=output_files,
            config_files=sorted(config_files),
            tags=sorted(list(tags))
        )

    @classmethod
    def verify(self,params,qc_outputs):
        """
        Verify 'fastq_strand' QC module against outputs

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
        if not params.seq_data_fastqs or \
           "fastq_strand.conf" not in qc_outputs.config_files:
            # No Fastqs or no conf file so strandedness
            # outputs not expected
            return None
        # Filter Fastq names
        fastqs = filter_fastqs(params.seq_data_reads[:1],
                               params.seq_data_fastqs)
        # Check that outputs exist for every Fastq
        for fq in fastqs:
            if fq not in qc_outputs.fastqs:
                return False
        return True

    @classmethod
    def add_to_pipeline(self,p,project_name,project,qc_dir,
                        organism,read_numbers,fastqs,star_indexes,
                        include_samples=None,nthreads=None,
                        fastq_subset=None,require_tasks=[],
                        verify_runner=None,compute_runner=None,
                        envmodules=None,verbose=False):
        """
        Adds tasks for 'cellranger_count' module to pipeline

        Arguments:
          p (Pipeline): pipeline to extend
          project_name (str): name of project
          project (AnalysisProject): project to run module on
          qc_dir (str): path to QC directory
          organism (str): name of organism(s)
          read_numbers (list): read numbers to include
          fastqs (list): Fastqs to run the module on
          star_indexes (mapping) map organism names to
            associated STAR indexes
          include_samples (list): subset of sample names to
            include
          fastq_subset (int): subset of reads to use for
            FastqScreen
          nthreads (int): number of threads (if not set then
            will be taken from the runner)
          require_tasks (list): list of tasks that the module
            needs to wait for
          verify_runner (JobRunner): runner to use for checks
          compute_runner (JobRunner): runner to use for
            computation
          verbose (bool): enable verbose output
        """
        # Set up fastq_strand.conf file
        setup_fastq_strand_conf = SetupFastqStrandConf(
            "%s: conf file for strandedness (fastq_strand)" %
            project_name,
            project,
            qc_dir=qc_dir,
            organism=organism,
            star_indexes=star_indexes
        )
        p.add_task(setup_fastq_strand_conf,
                   requires=require_tasks)

        # Check outputs for fastq_strand.py
        check_fastq_strand = CheckFastqStrandOutputs(
            "%s: check strandedness outputs (fastq_strand)" %
            project_name,
            project,
            qc_dir,
            setup_fastq_strand_conf.output.fastq_strand_conf,
            fastqs=fastqs,
            read_numbers=read_numbers,
            include_samples=include_samples,
            verbose=verbose
        )
        p.add_task(check_fastq_strand,
                   requires=(setup_fastq_strand_conf,),
                   runner=verify_runner)

        # Run fastq_strand.py
        run_fastq_strand = RunFastqStrand(
            "%s: strandedness (fastq_strand.py)" %
            project_name,
            check_fastq_strand.output.fastq_pairs,
            qc_dir,
            setup_fastq_strand_conf.output.fastq_strand_conf,
            fastq_strand_subset=fastq_subset,
            nthreads=nthreads
        )
        p.add_task(run_fastq_strand,
                   requires=(check_fastq_strand,),
                   runner=compute_runner,
                   envmodules=envmodules)
        return run_fastq_strand

#######################################################################
# Pipeline tasks
#######################################################################

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
    def init(self,project,qc_dir,fastq_strand_conf,fastqs=None,
             read_numbers=None,include_samples=None,verbose=False):
        """
        Initialise the CheckFastqStrandOutputs task.

        Arguments:
          project (AnalysisProject): project to run
            QC for
          qc_dir (str): directory for QC outputs (defaults
            to subdirectory 'qc' of project directory)
          fastq_strand_conf (str): path to the fastq_strand
            config file
          fastqs (list):  explicit list of Fastq files
            to check against (default is to use Fastqs
            from supplied analysis project)
          read_numbers (list): list of read numbers to
            include when checking outputs
          include_samples (list): optional, list of sample
            names to include
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
                      fastqs=self.args.fastqs,
                      read_numbers=self.args.read_numbers)
    def finish(self):
        for result in self.result():
            for fq_pair in result:
                if self.args.include_samples:
                    if self.args.project.fastq_attrs(fq_pair[0]).sample_name \
                       not in self.args.include_samples:
                        continue
                self.output.fastq_pairs.append(fq_pair)
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
             fastq_strand_subset=None,nthreads=None):
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

#######################################################################
# Helper functions
#######################################################################

def check_fastq_strand_outputs(project,qc_dir,fastq_strand_conf,
                               fastqs=None,read_numbers=None):
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
      fastqs (list): optional list of Fastqs to check
        against (defaults to Fastqs from the project)
      read_numbers (list): read numbers to predict
        outputs for

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
    if not fastqs:
        fastqs_in = project.fastqs
    else:
        fastqs_in = fastqs
    fastq_pairs = set()
    for fq_group in group_fastqs_by_name(
            remove_index_fastqs(fastqs_in,
                                project.fastq_attrs),
            fastq_attrs=project.fastq_attrs):
        # Assemble Fastq pairs based on read numbers
        if read_numbers:
            fq_pair = tuple([fq_group[r-1] for r in read_numbers])
        else:
            fq_pair = tuple([fq_group[0]])
        # Strand stats output
        output = os.path.join(qc_dir,
                              fastq_strand_output(fq_pair[0]))
        if not os.path.exists(output):
            fastq_pairs.add(fq_pair)
    return sorted(list(fastq_pairs))
