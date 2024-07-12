#!/usr/bin/env python3
#
#     fastqc: implements 'fastq_screen' QC module
#     Copyright (C) University of Manchester 2024 Peter Briggs

"""
Implements the 'fastq_screen' QC module:

* FastqScreen: core QCModule class
* CheckFastqScreenOutputs: pipeline task to check outputs
* RunFastqScreen: pipeline task to run FastqScreen
* check_fastq_screen_outputs: helper function for checking outputs
"""

#######################################################################
# Imports
#######################################################################

import os
import logging
from bcftbx.utils import AttributeDictionary
from . import QCModule
from ..fastq_screen import Fastqscreen as FastqscreenOutput
from ..fastq_screen import fastq_screen_output_files
from ..utils import filter_fastqs
from ...fastq_utils import remove_index_fastqs
from ...pipeliner import PipelineTask
from ...pipeliner import PipelineFunctionTask

# Module specific logger
logger = logging.getLogger(__name__)

# Data
from ..constants import FASTQ_SCREENS

#######################################################################
# Core class
#######################################################################

class FastqScreen(QCModule):
    """
    Class for handling the 'fastq_screen' QC module
    """
    name = "fastq_screen"
    mapped_metrics = True
    runners = ("verify_runner",
               "fastq_screen_runner")
    envmodules = ("fastq_screen",)
    
    def __init__(self):
        QCModule.__init__(self)

    @classmethod
    def collect_qc_outputs(self,qc_dir):
        """
        Collect information on FastqScreen outputs

        Returns an AttributeDictionary with the following
        attributes:

        - name: set to 'fastq_screen'
        - software: dictionary of software and versions
        - screen_names: list of associated panel names
        - fastqs: list of associated Fastq names
        - fastqs_for_screen: dictionary of panel names and lists
            of Fastq names associated with each panel
        - output_files: list of associated output files
        - tags: list of associated output classes

        Arguments:
          qc_dir (QCDir): QC directory to examine
        """
        versions = set()
        output_files = list()
        fastqs = set()
        fastqs_for_screen = dict()
        screen_names = set()
        tags = set()
        # Look for legacy screen files
        legacy_screens = list(filter(lambda f:
                                     f.endswith("_screen.txt") or
                                     f.endswith("_screen.png"),
                                     qc_dir.file_list))
        logger.debug("Screens (legacy): %s" % legacy_screens)
        # Look for new-style screen files
        screens = list(filter(lambda f:
                              "_screen_" in os.path.basename(f) and
                              (f.endswith(".txt") or
                               f.endswith(".png")),
                               qc_dir.file_list))
        logger.debug("Screens: %s" % screens)
        print("\t- %d fastq_screen files" % (len(legacy_screens) +
                                             len(screens)))
        if legacy_screens:
            for screen_name in FASTQ_SCREENS:
                # Explicitly check for legacy screens
                for screen in list(filter(lambda s:
                                          s.endswith("_%s_screen.txt" %
                                                     screen_name),
                                          legacy_screens)):
                    fq = qc_dir.fastq_attrs(screen[:-len("_screen.txt")])
                    fastq_name = str(fq)[:-len("_%s" % screen_name)]
                    tags.add("screens_%s%s" %
                             (('i' if fq.is_index_read else 'r'),
                              (fq.read_number
                               if fq.read_number is not None else '1')))
                    # Store general information
                    fastqs.add(fastq_name)
                    screen_names.add(screen_name)
                    versions.add(FastqscreenOutput(screen).version)
                    # Store Fastq against screen name
                    if screen_name not in fastqs_for_screen:
                        fastqs_for_screen[screen_name] = set()
                    fastqs_for_screen[screen_name].add(fastq_name)
            # Store the legacy screen files
            output_files.extend(legacy_screens)
        if screens:
            # Pull out the Fastq names from the .txt files
            for screen in list(filter(lambda s: s.endswith(".txt"),
                                      screens)):
                # Assume that names are 'FASTQ_screen_SCREENNAME.txt'
                fastq_name,screen_name = os.path.basename(screen)\
                                         [:-len(".txt")].\
                                         split("_screen_")
                fq = qc_dir.fastq_attrs(fastq_name)
                tags.add("screens_%s%s" %
                         (('i' if fq.is_index_read else 'r'),
                          (fq.read_number
                           if fq.read_number is not None else '1')))
                # Store general information
                fastqs.add(fastq_name)
                screen_names.add(screen_name)
                versions.add(FastqscreenOutput(screen).version)
                # Store Fastq against screen name
                if screen_name not in fastqs_for_screen:
                    fastqs_for_screen[screen_name] = set()
                fastqs_for_screen[screen_name].add(fastq_name)
            # Store the screen files
            output_files.extend(screens)
        # Return collected information
        if versions:
            software = { 'fastq_screen': sorted(list(versions)) }
        else:
            software = {}
        return AttributeDictionary(
            name='fastq_screen',
            software=software,
            screen_names=sorted(list(screen_names)),
            fastqs=sorted(list(fastqs)),
            fastqs_for_screen={ s: sorted(list(fastqs_for_screen[s]))
                                for s in fastqs_for_screen },
            output_files=output_files,
            tags=sorted(list(tags))
        )

    @classmethod
    def verify(self,params,qc_outputs):
        """
        Verify 'fastq_screen' QC module against outputs

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
        if not params.seq_data_fastqs or not params.fastq_screens:
            # Nothing to check
            return None
        try:
            # Filter Fastq names
            fastqs = filter_fastqs(params.seq_data_reads,
                                   params.seq_data_fastqs)
            # Check outputs exist for each screen
            for screen in params.fastq_screens:
                if screen not in qc_outputs.screen_names:
                    # No outputs associated with screen
                    return False
                # Check outputs exist for each Fastq
                for fq in fastqs:
                    if fq not in qc_outputs.fastqs_for_screen[screen]:
                        return False
            return True
        except KeyError as ex:
            # No FastqScreen outputs present
            return False

    @classmethod
    def add_to_pipeline(self,p,project_name,project,qc_dir,
                        screens,fastqs,read_numbers,include_samples=None,
                        nthreads=None,fastq_subset=None,
                        legacy=False,requires_tasks=[],
                        verify_runner=None,compute_runner=None,
                        envmodules=None,verbose=False):
        """
        Adds tasks for 'fastq_screen' module to pipeline

        Arguments:
          p (Pipeline): pipeline to extend
          project_name (str): name of project
          project (AnalysisProject): project to run module on
          qc_dir (str): path to QC directory
          screens (list): list of screen names
          fastqs (list): Fastqs to run the module on
          read_numbers (list): read numbers to include
          include_samples (list): subset of sample names to
            include
          nthreads (int): number of threads (if not set then
            will be taken from the runner)
          fastq_subset (int): subset of reads to use for
            FastqScreen
          legacy (bool): whether to use legacy naming for output
            files
          verbose (bool): enable verbose output
          require_tasks (list): list of tasks that the module
            needs to wait for
          verify_runner (JobRunner): runner to use for checks
          compute_runner (JobRunner): runner to use for
            computation
        """
        # Check outputs for FastqScreen
        check_fastq_screen = CheckFastqScreenOutputs(
            "%s: check FastqScreen outputs" %
            project_name,
            project,
            qc_dir,
            screens,
            fastqs=fastqs,
            read_numbers=read_numbers,
            include_samples=include_samples,
            fastq_attrs=project.fastq_attrs,
            legacy=legacy,
            verbose=verbose
        )
        p.add_task(check_fastq_screen,
                   requires=requires_tasks,
                   runner=verify_runner)
        # Run FastqScreen
        run_fastq_screen = RunFastqScreen(
            "%s: FastqScreen" % project_name,
            check_fastq_screen.output.fastqs,
            qc_dir,
            screens,
            subset=fastq_subset,
            nthreads=nthreads,
            read_numbers=read_numbers,
            fastq_attrs=project.fastq_attrs,
            legacy=legacy
        )
        p.add_task(run_fastq_screen,
                   requires=(check_fastq_screen,),
                   runner=compute_runner,
                   envmodules=envmodules)
        return run_fastq_screen

#######################################################################
# Pipeline tasks
#######################################################################

class CheckFastqScreenOutputs(PipelineFunctionTask):
    """
    Check the outputs from FastqScreen
    """
    def init(self,project,qc_dir,screens,fastqs=None,
             read_numbers=None,include_samples=None,
             fastq_attrs=None,legacy=False,verbose=False):
        """
        Initialise the CheckFastqScreenOutputs task.

        Arguments:
          project (AnalysisProject): project to run
            QC for
          qc_dir (str): directory for QC outputs (defaults
            to subdirectory 'qc' of project directory)
          screens (mapping): mapping of screen names to
            FastqScreen conf files
          fastqs (list): explicit list of Fastq files
            to check against (default is to use Fastqs
            from supplied analysis project)
          read_numbers (list): read numbers to include
          include_samples (list): optional, list of sample
            names to include
          fastq_attrs (BaseFastqAttrs): class to use for
            extracting data from Fastq names
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
                          fastqs=self.args.fastqs,
                          read_numbers=self.args.read_numbers,
                          legacy=self.args.legacy)
    def finish(self):
        fastqs = set()
        for result in self.result():
            for fq in result:
                if self.args.include_samples:
                    if self.args.fastq_attrs(fq).sample_name \
                       not in self.args.include_samples:
                        continue
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
             read_numbers=None,fastq_attrs=None,legacy=False):
        """
        Initialise the RunFastqScreen task.

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
          read_numbers (list): list of read numbers to
            include when running Fastq Screen
          fastq_attrs (BaseFastqAttrs): class to use for
            extracting data from Fastq names
          legacy (bool): if True then use 'legacy' naming
            convention for output files (default is to
            use new format)
        """
        self.conda("fastq-screen=0.15.3",
                   "bowtie=1.2.3")
        # Also need to specify tbb=2020.2 for bowtie
        # See https://www.biostars.org/p/494922/
        self.conda("tbb=2020.2")
        # perl-gd seems to be needed explicitly
        self.conda("perl-gd")
    def setup(self):
        if not self.args.fastqs:
            print("Nothing to do")
            return
        # Report if legacy naming convention will be used
        if self.args.legacy:
            print("Using legacy FastqScreen output names")
        # Set up the FastqScreen runs for each Fastq
        for fastq in self.args.fastqs:
            if self.args.read_numbers and \
               self.args.fastq_attrs(fastq).read_number not in \
               self.args.read_numbers:
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

#######################################################################
# Helper functions
#######################################################################

#FIXME copied from qc/outputs.py
def check_fastq_screen_outputs(project,qc_dir,screen,fastqs=None,
                               read_numbers=None,legacy=False):
    """
    Return Fastqs missing QC outputs from FastqScreen

    Returns a list of the Fastqs from a project for which
    one or more associated outputs from FastqScreen
    don't exist in the specified QC directory.

    Arguments:
      project (AnalysisProject): project to check the
        QC outputs for
      qc_dir (str): path to the QC directory (relative
        path is assumed to be a subdirectory of the
        project)
      screen (str): screen name to check
      fastqs (list): optional list of Fastqs to check
        against (defaults to Fastqs from the project)
      read_numbers (list): read numbers to define Fastqs
        to predict outputs for; if not set then all
        non-index reads will be included
      legacy (bool): if True then check for 'legacy'-style
         names (defult: False)

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
            # Ignore non-data reads
            continue
        for output in [os.path.join(qc_dir,f)
                       for f in fastq_screen_output_files(fastq,screen,
                                                          legacy=legacy)]:
            if not os.path.exists(output):
                fastqs.add(fastq)
    return sorted(list(fastqs))
