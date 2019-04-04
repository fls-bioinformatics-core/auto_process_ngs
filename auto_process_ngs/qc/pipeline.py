#!/usr/bin/env python
#
#     qc.pipeline.py: pipelines for running QC
#     Copyright (C) University of Manchester 2019 Peter Briggs
#

"""
Pipeline components for running the QC pipeline.

Pipeline classes:

- QCPipeline

Pipeline task classes:

- SetupQCDirs
- CheckIlluminaQCOutputs
- RunIlluminaQC
- SetupFastqStrandConf
- CheckFastqStrandOutputs
- RunFastqStrand
- ReportQC
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
from ..analysis import AnalysisProject
from ..analysis import AnalysisFastq
from ..analysis import copy_analysis_project
from ..applications import Command
from ..fastq_utils import pair_fastqs_by_name
from ..fastq_utils import remove_index_fastqs
from ..pipeliner import Pipeline
from ..pipeliner import PipelineTask
from ..pipeliner import PipelineFunctionTask
from ..pipeliner import PipelineCommandWrapper
from ..pipeliner import PipelineParam as Param
from ..pipeliner import PipelineFailure
from ..utils import get_organism_list
from .constants import FASTQ_SCREENS
from .outputs import fastqc_output
from .outputs import fastq_screen_output
from .outputs import fastq_strand_output
from .outputs import check_illumina_qc_outputs
from .outputs import check_fastq_strand_outputs
from .outputs import expected_outputs
from .utils import determine_qc_protocol
from .fastq_strand import build_fastq_strand_conf

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
        self.add_param('nthreads',type=int,value=1)
        self.add_param('fastq_subset',type=int)
        self.add_param('fastq_strand_indexes',type=dict)

        # Define runners
        self.add_runner('verify_runner')
        self.add_runner('qc_runner')
        self.add_runner('report_runner')

    def add_project(self,project,qc_dir=None,organism=None,fastq_dir=None,
                    qc_protocol=None,multiqc=False,
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

        # Report details
        self.report("-- Protocol  : %s" % qc_protocol)
        self.report("-- Directory : %s" % project.dirn)
        self.report("-- Fastqs dir: %s" % project.fastq_dir)
        self.report("-- QC dir    : %s" % qc_dir)

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

        # Check illumina_qc.sh is compatible version
        check_illumina_qc_version = CheckIlluminaQCVersion(
            "%s: check illumina_qc.sh version" %
            project_name)
        self.add_task(check_illumina_qc_version,
                      log_dir=log_dir)

        # Check outputs for illumina_qc.sh
        check_illumina_qc = CheckIlluminaQCOutputs(
            "%s: check basic QC outputs (illumina_qc.sh)" %
            project_name,
            project,
            qc_dir,
            qc_protocol=qc_protocol,
            verbose=self.params.VERBOSE
        )
        self.add_task(check_illumina_qc,
                      requires=(setup_qc_dirs,
                                check_illumina_qc_version,),
                      runner=self.runners['verify_runner'],
                      log_dir=log_dir)

        # Run illumina_qc.sh
        run_illumina_qc = RunIlluminaQC(
            "%s: basic QC (illumina_qc.sh)" % project_name,
            check_illumina_qc.output.fastqs,
            qc_dir,
            fastq_screen_subset=self.params.fastq_subset,
            nthreads=self.params.nthreads,
            qc_protocol=qc_protocol,
            fastq_attrs=project.fastq_attrs
        )
        self.add_task(run_illumina_qc,
                      requires=(check_illumina_qc,),
                      runner=self.runners['qc_runner'],
                      log_dir=log_dir)

        # Set up fastq_strand.conf file
        setup_fastq_strand_conf = SetupFastqStrandConf(
            "%s: conf file for strandedness (fastq_strand)" %
            project_name,
            project,
            qc_dir=qc_dir,
            organism=organism,
            star_indexes=self.params.fastq_strand_indexes
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
                      runner=self.runners['qc_runner'],
                      log_dir=log_dir)

        # Make QC report
        report_qc = ReportQC(
            "%s: make QC report" % project_name,
            project,
            qc_dir,
            multiqc=multiqc
        )
        self.add_task(report_qc,
                      requires=(run_illumina_qc,
                                run_fastq_strand,),
                      runner=self.runners['report_runner'],
                      log_dir=log_dir)

    def run(self,nthreads=None,fastq_strand_indexes=None,
            fastq_subset=None,working_dir=None,log_file=None,
            batch_size=None,max_jobs=1,poll_interval=5,
            runners=None,default_runner=None,verbose=False):
        """
        Run the tasks in the pipeline

        Arguments:
          nthreads (int): number of threads/processors to
            use for QC jobs (defaults to 1)
          fastq_strand_indexes (dict): mapping of organism
            IDs to directories with STAR index
          fastq_subset (int): explicitly specify
            the subset size for subsetting running Fastqs
          working_dir (str): optional path to a working
            directory (defaults to temporary directory in
            the current directory)
          log_dir (str): path of directory where log files
            will be written to
          batch_size (int): if set then run commands in
            each task in batches, with each batch running
            this many commands at a time (default is to run
            one command per job)
          max_jobs (int): optional maximum number of
            concurrent jobs in scheduler (defaults to 1)
          poll_interval (float): optional polling interval
            (seconds) to set in scheduler (defaults to 5s)
          runners (dict): mapping of names to JobRunner
            instances; valid names are 'qc_runner',
            'report_runner','verify_runner','default'
          default_runner (JobRunner): optional default
            job runner to use
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
        log_dir = os.path.join(working_dir,"logs")
        scripts_dir = os.path.join(working_dir,"scripts")

        # Execute the pipeline
        status = Pipeline.run(self,
                              working_dir=working_dir,
                              log_dir=log_dir,
                              scripts_dir=scripts_dir,
                              log_file=log_file,
                              batch_size=batch_size,
                              exit_on_failure=PipelineFailure.DEFERRED,
                              params={
                                  'nthreads': nthreads,
                                  'fastq_subset': fastq_subset,
                                  'fastq_strand_indexes': fastq_strand_indexes,
                              },
                              max_jobs=max_jobs,
                              runners=runners,
                              default_runner=default_runner,
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
                           % (self.project.name,
                              stored_protocol,
                              self.args.protocol))
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
        # Store the QC protocol data
        qc_info['protocol'] = self.args.qc_protocol
        qc_info['fastq_dir'] = self.args.project.fastq_dir
        qc_info.save()

class CheckIlluminaQCVersion(PipelineTask):
    """
    Check the illumina_qc.sh version is compatible with the pipeline
    """
    def init(self,compatible_versions=('1.3.2','1.3.3')):
        """
        Initialise the CheckIlluminaQC task
        """
        pass
    def setup(self):
        version = self.illumina_qc_version()
        if version not in self.args.compatible_versions:
            self.fail(message="QC script version is %s, needs %s" %
                      (version,'/'.join(compatible_versions)))
    def illumina_qc_version(self):
        # Return version of illumina_qc.sh script
        status,qc_script_info = Command(
            'illumina_qc.sh',
            '--version').subprocess_check_output()
        if status == 0:
            return qc_script_info.strip().split()[-1]

class CheckIlluminaQCOutputs(PipelineFunctionTask):
    """
    Check the outputs from the illumina_qc.sh script
    """
    def init(self,project,qc_dir,qc_protocol=None,
             verbose=False):
        """
        Initialise the CheckIlluminaQCOutputs task.

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
            missing outputs from illumina_qc.sh under
            the specified QC protocol
        """
        self.add_output('fastqs',list())
    def setup(self):
        self.add_call("Check illumina_qc.sh outputs for %s"
                      % self.args.project.name,
                      check_illumina_qc_outputs,
                      self.args.project,
                      self.args.qc_dir,
                      self.args.qc_protocol)
    def finish(self):
        for result in self.result():
            self.output.fastqs.extend(result)
        if self.output.fastqs:
            if self.args.verbose:
                print "Fastqs with missing QC outputs from " \
                    "illumina_qc.sh:"
                for fq in self.output.fastqs:
                    print "-- %s" % fq
            else:
                print "%s Fastqs with missing QC outputs from " \
                    "illumina_qc.sh" % len(self.output.fastqs)
        else:
            print "No Fastqs with missing QC outputs from " \
                "illumina_qc.sh"

class RunIlluminaQC(PipelineTask):
    """
    Run the illumina_qc.sh script
    """
    def init(self,fastqs,qc_dir,fastq_screen_subset=None,nthreads=1,
             qc_protocol=None,fastq_attrs=None):
        """
        Initialise the RunIlluminaQC task.

        Arguments:
          fastqs (list): list of paths to Fastq files to
            run the illumina_qc.sh script on (it is
            expected that this list will come from the
            CheckIlluminaQCOutputs task)
          qc_dir (str): directory for QC outputs (defaults
            to subdirectory 'qc' of project directory)
          fastq_screen_subset (int): explicitly specify
            the subset size for running Fastq_screen
          nthreads (int): number of threads/processors to
            use (defaults to 1)
          qc_protocol (str): QC protocol to use
          fastq_attrs (BaseFastqAttrs): class to use for
            extracting data from Fastq names
        """
        pass
    def setup(self):
        if not self.args.fastqs:
            print "Nothing to do"
            return
        # Set up the illumina_qc.sh runs for each Fastq
        for fastq in self.args.fastqs:
            cmd = PipelineCommandWrapper(
                "Run illumina_qc.sh for %s" % os.path.basename(fastq),
                'illumina_qc.sh',fastq,
                '--threads',self.args.nthreads,
                '--qc_dir',os.path.abspath(self.args.qc_dir))
            if self.args.fastq_screen_subset is not None:
                cmd.add_args('--subset',self.args.fastq_screen_subset)
            # No screens for R1 reads for single cell
            if self.args.qc_protocol == 'singlecell' and \
               self.args.fastq_attrs(fastq).read_number == 1:
                cmd.add_args('--no-screens')
            # Add the command
            self.add_cmd(cmd)
    def run_illumina_qc(self,fastq,qc_dir,fastq_screen_subset=None,
                        nthreads=1,qc_protocol=None,fastq_attrs=None):
        # How to extract attributes from Fastq names
        if fastq_attrs is None:
            fastq_attrs = AnalysisFastq
        # Build illumina_qc.sh command
        cmd = Command('illumina_qc.sh',fastq,
                      '--threads',nthreads,
                      '--qc_dir',os.path.abspath(qc_dir))
        if fastq_screen_subset is not None:
            cmd.add_args('--subset',fastq_screen_subset)
        # No screens for for in single cell
        if qc_protocol == 'singlecell' and \
           fastq_attrs(fastq).read_number == 1:
            cmd.add_args('--no-screens')
        # Execute the command
        status = cmd.run_subprocess(working_dir=qc_dir)
        if status != 0:
            raise Exception("Error running illumina_qc.sh for %s: "
                            "exit code: %s" % (fastq,status))
        # Check the Fastqc outputs
        failed = False
        for output in [os.path.join(qc_dir,f)
                       for f in fastqc_output(fastq)]:
            if not os.path.exists(output):
                failed = True
        # Check the Fastq_screen outputs
        if qc_protocol == 'singlecell' and \
           fastq_attrs(fastq).read_number == 1:
            # No screens for R1 for single cell
            pass
        else:
            for screen in FASTQ_SCREENS:
                for output in [os.path.join(qc_dir,f)
                               for f in fastq_screen_output(fastq,
                                                            screen)]:
                    if not os.path.exists(output):
                        failed = True
        if failed:
            raise Exception("illumina_qc.sh failed to produce all "
                            "outputs for %s" % (fastq,))

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
            print "Removing existing conf file: %s" % fastq_strand_conf
            os.remove(fastq_strand_conf)
        # Sort out organism(s)
        if not organism:
            organism = project.info.organism
        if not organism:
            print "No organisms specified"
            return
        print "Organism(s): %s" % organism
        # Fetch matching STAR indexes
        if not star_indexes:
            print "No STAR indexes available"
            return
        print "STAR indexes: %s" % star_indexes
        fastq_strand_indexes = build_fastq_strand_conf(
            get_organism_list(organism),
            star_indexes)
        if not fastq_strand_indexes:
            print "No matching indexes for strandedness determination"
            return
        # Create the conf file
        print "Writing conf file: %s" % fastq_strand_conf
        with open(fastq_strand_conf,'w') as fp:
            fp.write("%s\n" % fastq_strand_indexes)
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
            print "No conf file, nothing to check"
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
                print "Fastq pairs with missing QC outputs from " \
                    "fastq_strand.py:"
                for fq_pair in self.output.fastq_pairs:
                    print "-- %s" % (fq_pair,)
            else:
                print "%s Fastq pairs with missing QC outputs from " \
                    "fastq_strand.py" % len(self.output.fastq_pairs)
        else:
            print "No Fastqs with missing QC outputs from " \
                "fastq_strand.py"

class RunFastqStrand(PipelineFunctionTask):
    """
    Run the fastq_strand.py utility
    """
    def init(self,fastq_pairs,qc_dir,fastq_strand_conf,
             fastq_strand_subset=None,nthreads=1,
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
            use (defaults to 1)
          qc_protocol (str): QC protocol to use
        """
        pass
    def setup(self):
        if not self.args.fastq_pairs:
            print "No Fastqs: nothing to do"
            return
        elif not self.args.fastq_strand_conf or \
             not os.path.exists(self.args.fastq_strand_conf):
            print "No conf file: nothing to do"
            return
        for fastq_pair in self.args.fastq_pairs:
            self.add_call("Run fastq_strand.py for %s..." %
                          os.path.basename(fastq_pair[0]),
                          self.run_fastq_strand,
                          fastq_pair,
                          self.args.qc_dir,
                          fastq_strand_subset=
                          self.args.fastq_strand_subset,
                          fastq_strand_conf=self.args.fastq_strand_conf,
                          nthreads=self.args.nthreads,
                          qc_protocol=self.args.qc_protocol)
    def run_fastq_strand(self,fastq_pair,qc_dir,fastq_strand_conf,
                         fastq_strand_subset=None,nthreads=1,
                         qc_protocol=None):
        # Build fastq_strand.py command
        cmd = Command('fastq_strand.py',
                              '-n',nthreads,
                              '--conf',fastq_strand_conf,
                              '--outdir',os.path.abspath(qc_dir))
        if fastq_strand_subset:
            cmd.add_args('--subset',fastq_strand_subset)
        cmd.add_args(*fastq_pair)
        # Execute the command
        status = cmd.run_subprocess(working_dir=qc_dir)
        if status != 0:
            raise Exception("Error running fastq_strand.py for %s: "
                            "exit code: %s" % (fastq_pair,status))
        # Check the outputs
        output = fastq_strand_output(fastq_pair[0])
        failed = (not os.path.exists(os.path.join(qc_dir,output)))
        if failed:
            raise Exception("fastq_strand.py failed to produce "
                            "output for %s" % (fastq_pair,))

class ReportQC(PipelineFunctionTask):
    """
    Generate the QC report
    """
    def init(self,project,qc_dir,report_html=None,fastq_dir=None,
             multiqc=False,zip_outputs=True):
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
          zip_outputs (bool): if True then also generate
            a ZIP archive of the QC reports
        """
        pass
    def setup(self):
        self.add_call("Generate QC report for %s"
                      % self.args.project.name,
                      self.report_qc,
                      self.args.project,
                      self.args.qc_dir,
                      report_html=self.args.report_html,
                      fastq_dir=self.args.fastq_dir,
                      multiqc=self.args.multiqc,
                      zip_outputs=self.args.zip_outputs)
    def report_qc(self,project,qc_dir,report_html=None,
                  fastq_dir=None,multiqc=False,zip_outputs=True):
        # Build the reportqc.py command
        qc_base = os.path.basename(qc_dir)
        if report_html is None:
            out_file = '%s_report.html' % qc_base
        else:
            out_file = report_html
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
        report_cmd = Command(
            "reportqc.py",
            "--qc_dir",qc_dir,
            "--filename",out_file,
            "--title",title)
        if fastq_dir is not None:
            report_cmd.add_args("--fastq_dir",fastq_dir)
        if multiqc:
            report_cmd.add_args("--multiqc")
        if zip_outputs:
            report_cmd.add_args("--zip")
        report_cmd.add_args(project.dirn)
        # Execute the command
        status = report_cmd.run_subprocess(working_dir=qc_dir)
        if status != 0:
            raise Exception("Error running reportqc.py for %s: "
                            "exit code: %s" % (project.name,status))
