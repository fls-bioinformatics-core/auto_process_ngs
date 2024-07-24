#!/usr/bin/env python3
#
#     rseqc_infer_experiment: implements 'rseqc_infer_experiment' QC module
#     Copyright (C) University of Manchester 2024 Peter Briggs

"""
Implements the 'rseqc_infer_experiment' QC module:

* RseqcInferExperiment: core QCModule class
* RunRSeQCGenebodyCoverage: pipeline task to run 'infer_experiment.py'
"""

#######################################################################
# Imports
#######################################################################

import os
import shutil
import logging
from bcftbx.utils import AttributeDictionary
from . import QCModule
from ..rseqc import InferExperiment as InferExperimentOutput
from ..utils import filter_fastqs
from ..utils import get_bam_basename
from ..utils import read_versions_file
from ...pipeliner import PipelineTask
from ...pipeliner import PipelineParam as Param
from ...utils import normalise_organism_name

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Core class
#######################################################################

class RseqcInferExperiment(QCModule):
    """
    Class for handling the 'rseqc_infer_experiment' QC module
    """
    name = "rseqc_infer_experiment"
    mapped_metrics = True
    require_bam_files = True
    runners = ("rseqc_runner",)
    envmodules = ()
    
    def __init__(self):
        QCModule.__init__(self)

    @classmethod
    def collect_qc_outputs(self,qc_dir):
        """
        Collect information on RSeQC infer_experiment.py outputs

        Returns an AttributeDictionary with the following
        attributes:

        - name: set to 'rseqc_infer_experiment'
        - software: dictionary of software and versions
        - organisms: list of organisms with associated
          outputs
        - bam_files: list of associated BAM file names
        - output_files: list of associated output files
        - tags: list of associated output classes

        Arguments:
          qc_dir (QCDir): QC directory to examine
        """
        software = {}
        bam_files = set()
        output_files = list()
        tags = set()
        # Look for RSeQC infer_experiment.py outputs
        organisms = set()
        infer_experiment_dir = os.path.join(qc_dir.path,
                                            "rseqc_infer_experiment")
        if os.path.isdir(infer_experiment_dir):
            # Look for subdirs with organism names
            for d in filter(
                    lambda dd:
                    os.path.isdir(os.path.join(infer_experiment_dir,dd)),
                    os.listdir(infer_experiment_dir)):
                # Check for outputs
                for f in filter(
                        lambda ff:
                        ff.endswith(".infer_experiment.log"),
                        os.listdir(os.path.join(infer_experiment_dir,d))):
                    name = f[:-len(".infer_experiment.log")]
                    organisms.add(d)
                    bam_files.add(name)
                    output_files.append(os.path.join(infer_experiment_dir,
                                                     d,f))
                # Check for software information
                software = read_versions_file(
                    os.path.join(infer_experiment_dir,
                                 d,"_versions"),
                    software)
        if organisms:
            tags.add("rseqc_infer_experiment")
            if not software:
                software['rseqc:infer_experiment'] = [None]
        # Return collected information
        return AttributeDictionary(
            name='rseqc_infer_experiment',
            software=software,
            bam_files=sorted(list(bam_files)),
            organisms=sorted(list(organisms)),
            output_files=output_files,
            tags=sorted(list(tags))
        )

    @classmethod
    def verify(self,params,qc_outputs):
        """
        Verify 'rseqc_infer_experiment' QC module against outputs

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
        if not params.seq_data_fastqs:
            # Nothing to check
            return None
        if not params.organism:
            # No organism specified
            return None
        if not params.star_index or not params.annotation_bed:
            # No STAR index or annotation
            return None
        if normalise_organism_name(params.organism) not in \
           qc_outputs.organisms:
            return False
        # Filter Fastq names and convert to BAM names
        if params.seq_data_reads:
            fastqs = filter_fastqs(params.seq_data_reads[:1],
                                   params.seq_data_fastqs)
        else:
            # No Fastqs to get BAM names from
            return None
        bams = [get_bam_basename(fq) for fq in fastqs]
        # Check that outputs exist for every BAM
        for bam in bams:
            if bam not in qc_outputs.bam_files:
                return False
        return True

    @classmethod
    def add_to_pipeline(self,p,project_name,qc_dir,bam_files,
                        reference_gene_model,organism_name,
                        required_tasks=[],rseqc_runner=None):
        """
        Adds tasks for 'rseqc_infer_experiment' module to pipeline

        Arguments:
          p (Pipeline): pipeline to extend
          project_name (str): name of project
          qc_dir (str): path to QC directory
          bam_files (list): BAM files to run the module on
          reference_gene_model (str): path to reference gene
            model BED file
          organism_name (str): normalised name for organism
            that BAMs are aligned to
          required_tasks (list): list of tasks that the module
            needs to wait for
          rseqc_runner (JobRunner): runner to use for RSeQC
        """
        out_dir = os.path.join(qc_dir,
                               'rseqc_infer_experiment',
                               organism_name)
        rseqc_infer_experiment = RunRSeQCInferExperiment(
            "%s: infer experiment from BAM files (RSeQC)" %
            project_name,
            bam_files,
            reference_gene_model,
            out_dir)
        p.add_task(rseqc_infer_experiment,
                   require_tasks=required_tasks)
        return rseqc_infer_experiment

#######################################################################
# Pipeline tasks
#######################################################################

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
        get_version = False
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
                get_version = True
        # Get version of RSeQC
        if get_version:
            self.add_cmd("Get RSeQC infer_experiment.py version",
                         """
                         infer_experiment.py --version >_versions 2>&1
                         """)
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
            infer_expt = InferExperimentOutput(infer_expt_log)
            # Dump data associated with previous BAM
            outputs[bam_file] = {
                'paired_end': infer_expt.paired_end,
                'unstranded': infer_expt.unstranded,
                'forward': infer_expt.forward,
                'reverse': infer_expt.reverse,
            }
        # RSeQC version
        if os.path.exists("_versions"):
            rseqc_infer_experiment_version = None
            with open("_versions",'rt') as fp:
                for line in fp:
                    if line.startswith("infer_experiment.py "):
                        # Example: infer_experiment.py 4.0.0
                        rseqc_infer_experiment_version = \
                            ' '.join(line.strip().split(' ')[1:])
            if rseqc_infer_experiment_version:
                with open("_versions",'wt') as fp:
                    fp.write("rseqc:infer_experiment\t%s\n" %
                             rseqc_infer_experiment_version)
            shutil.copy("_versions",self.args.out_dir)
        # Set output
        self.output.experiments.set(outputs)

#######################################################################
# Helper functions
#######################################################################

#TBA
