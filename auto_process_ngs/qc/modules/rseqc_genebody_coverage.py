#!/usr/bin/env python3
#
#     rseqc_genebody_coverage: implements 'rsecq_genebody_coverage' QC module
#     Copyright (C) University of Manchester 2024 Peter Briggs

"""
Implements the 'rseqc_genebody_coverage' QC module:

* RSeqcGenebodyCoverage: core QCModule class
* RunRSeQCGenebodyCoverage: pipeline task to run 'genebody_coverage.py'
"""

#######################################################################
# Imports
#######################################################################

import os
import shutil
import logging
from bcftbx.utils import AttributeDictionary
from . import QCModule
from ..apps.rseqc import rseqc_genebody_coverage_output
from ..utils import read_versions_file
from ...pipeliner import PipelineTask
from ...utils import normalise_organism_name

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Core class
#######################################################################

class RseqcGenebodyCoverage(QCModule):
    """
    Class for handling the 'rseqc_genebody_coverage' QC module
    """
    name = "rseqc_genebody_coverage"
    mapped_metrics = True
    require_bam_files = True
    runners = ("rseqc_runner",)
    envmodules = ()
    
    def __init__(self):
        QCModule.__init__(self)

    @classmethod
    def collect_qc_outputs(self,qc_dir):
        """
        Collect information on RSeQC geneBody_coverage.py outputs

        Returns an AttributeDictionary with the following
        attributes:

        - name: set to 'rseqc_genebody_coverage'
        - software: dictionary of software and versions
        - organisms: list of organisms with associated
          outputs
        - output_files: list of associated output files
        - tags: list of associated output classes

        Arguments:
          qc_dir (QCDir): QC directory to examine
        """
        software = {}
        output_files = list()
        tags = set()
        # Look for RSeQC geneBody_coverage.py outputs
        organisms = set()
        genebody_cov_dir = os.path.join(qc_dir.path,
                                        "rseqc_genebody_coverage")
        if os.path.isdir(genebody_cov_dir):
            # Look for subdirs with organism names
            for d in filter(
                    lambda dd:
                    os.path.isdir(os.path.join(genebody_cov_dir,dd)),
                    os.listdir(genebody_cov_dir)):
                # Check for outputs
                for f in filter(
                        lambda ff:
                        ff.endswith(".geneBodyCoverage.txt"),
                        os.listdir(os.path.join(genebody_cov_dir,d))):
                    name = f[:-len(".geneBodyCoverage.txt")]
                    outputs = rseqc_genebody_coverage_output(
                        name,
                        prefix=os.path.join(genebody_cov_dir,d))
                    if all([os.path.exists(f) for f in outputs]):
                        # All outputs present
                        organisms.add(d)
                        output_files.extend(outputs)
                # Check for software information
                software = read_versions_file(
                    os.path.join(genebody_cov_dir,
                                 d,"_versions"),
                    software)
        if organisms:
            tags.add("rseqc_genebody_coverage")
            if not software:
                software['rseqc:genebody_coverage'] = [None]
        # Return collected information
        return AttributeDictionary(
            name='rseqc_genebody_coverage',
            software=software,
            organisms=sorted(list(organisms)),
            output_files=output_files,
            tags=sorted(list(tags))
        )

    @classmethod
    def verify(self,params,qc_outputs):
        """
        Verify 'rseqc_genebody_coverage' QC module against outputs

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
        if normalise_organism_name(params.organism) \
           not in qc_outputs.organisms:
            return False
        return True

    @classmethod
    def add_to_pipeline(self,p,project_name,qc_dir,bam_files,
                        reference_gene_model,organism_name,
                        required_tasks=[],rseqc_runner=None):
        """
        Adds tasks for 'rseqc_genebody_coverage' module to pipeline

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
        # Run RSeQC gene body coverage
        out_dir = os.path.join(qc_dir,
                               'rseqc_genebody_coverage',
                               organism_name)
        rseqc_gene_body_coverage = RunRSeQCGenebodyCoverage(
            "%s: calculate gene body coverage (RSeQC)" % project_name,
            bam_files,
            reference_gene_model,
            out_dir,
            name=project_name)
        p.add_task(rseqc_gene_body_coverage,
                   requires=required_tasks,
                   runner=rseqc_runner)
        return rseqc_gene_body_coverage

#######################################################################
# Pipeline tasks
#######################################################################

class RunRSeQCGenebodyCoverage(PipelineTask):
    """
    Run RSeQC's 'genebody_coverage.py' on BAM files

    Given a collection of BAM files, runs the RSeQC
    'genebody_coverage.py' utility
    (http://rseqc.sourceforge.net/#genebody-coverage-py).
    """
    def init(self,bam_files,reference_gene_model,out_dir,name="rseqc"):
        """
        Initialise the RunRSeQCGenebodyCoverage task

        Arguments:
          bam_files (list): list of paths to BAM files
            to run genebody_coverage.py on
          reference_gene_model (str): path to BED file
            with the reference gene model data
          out_dir (str): path to a directory where the
            output files will be written
          name (str): optional basename for the output
            files (defaults to 'rseqc')
        """
        # Conda dependencies
        self.conda("rseqc=4.0.0",
                   "r-base=4")
    def setup(self):
        # Check we have BAM files
        if len(self.args.bam_files) < 1:
            print("No BAM files, cannot run RSeQC genebody_coverage.py")
            return
        # Check if outputs already exist
        outputs_exist = True
        for f in rseqc_genebody_coverage_output(self.args.name,
                                                self.args.out_dir):
            outputs_exist = (outputs_exist and os.path.exists(f))
        if outputs_exist:
            print("All outputs exist already, nothing to do")
            return
        # Check for reference gene model
        if self.args.reference_gene_model:
            print("Reference gene model: %s" %
                  self.args.reference_gene_model)
        else:
            print("Reference gene model is not set, cannot run RSeQC "
                  "geneBody_coverage.py")
            return
        # Set up command to run genebody_coverage.py
        self.add_cmd("Run RSeQC geneBody_coverage.py",
                     """
                     # Get version
                     geneBody_coverage.py --version >_versions 2>&1
                     # Run geneBody_coverage
                     geneBody_coverage.py \\
                         -r {reference_gene_model} \\
                         -i {bam_files} \\
                         -f png \\
                         -o {basename}
                     """.format(
                         reference_gene_model=self.args.reference_gene_model,
                         bam_files=','.join(self.args.bam_files),
                         basename=self.args.name))
    def finish(self):
        if not self.args.reference_gene_model or \
           len(self.args.bam_files) < 1:
            return
        # Copy outputs to final location
        if not os.path.exists(self.args.out_dir):
            os.makedirs(self.args.out_dir,exist_ok=True)
        for f in rseqc_genebody_coverage_output(self.args.name,
                                                self.args.out_dir):
            if not os.path.exists(f):
                # Copy new version to ouput location
                shutil.copy(os.path.basename(f),self.args.out_dir)
        # RSeQC version version
        if os.path.exists("_versions"):
            rseqc_genebody_coverage_version = None
            with open("_versions",'rt') as fp:
                for line in fp:
                    if line.startswith("geneBody_coverage.py "):
                        # Example: geneBody_coverage.py 4.0.0
                        rseqc_genebody_coverage_version = \
                            ' '.join(line.strip().split(' ')[1:])
            if rseqc_genebody_coverage_version:
                with open("_versions",'wt') as fp:
                    fp.write("rseqc:genebody_coverage\t%s\n" %
                             rseqc_genebody_coverage_version)
            shutil.copy("_versions",self.args.out_dir)
