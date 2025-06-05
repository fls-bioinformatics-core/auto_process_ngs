#!/usr/bin/env python3
#
#     qualimap_rnaseq: implements 'qualimap_rnaseq' QC module
#     Copyright (C) University of Manchester 2024 Peter Briggs

"""
Implements the 'qualimap_rnaseq' QC module:

* QualimapRnaseq: core QCModule class
* RunQualimapRnaseq: pipeline task to run Qualimap 'rnaseq'
"""

#######################################################################
# Imports
#######################################################################

import os
import shutil
import logging
from bcftbx.utils import AttributeDictionary
from . import QCModule
from ..apps.qualimap import qualimap_rnaseq_output
from ..utils import filter_fastqs
from ..utils import read_versions_file
from ..utils import get_bam_basename
from ...pipeliner import PipelineTask
from ...utils import normalise_organism_name

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Core class
#######################################################################

class QualimapRnaseq(QCModule):
    """
    Class for handling the 'qualimap_rnaseq' QC module
    """
    name = "qualimap_rnaseq"
    mapped_metrics = True
    require_bam_files = True
    runners = ("qualimap_runner",)
    envmodules = ()
    
    def __init__(self):
        QCModule.__init__(self)

    @classmethod
    def collect_qc_outputs(self,qc_dir):
        """
        Collect information on Qualimap 'rnaseq' outputs

        Returns an AttributeDictionary with the following
        attributes:

        - name: set to 'qualimap_rnaseq'
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
        # Look for Qualimap 'rnaseq' outputs
        organisms = set()
        qualimap_dir = os.path.join(qc_dir.path,"qualimap-rnaseq")
        if os.path.isdir(qualimap_dir):
            # Look for subdirs with organism names
            for d in filter(
                    lambda dd:
                    os.path.isdir(os.path.join(qualimap_dir,dd)),
                    os.listdir(qualimap_dir)):
                # Look for subdirs with BAM file names
                organism_dir = os.path.join(qualimap_dir,d)
                for bam in filter(
                        lambda dd:
                        os.path.isdir(os.path.join(organism_dir,dd)),
                        os.listdir(organism_dir)):
                    # Check for Qualimap rnaseq outputs
                    for f in filter(
                            lambda ff:
                            ff == "qualimapReport.html",
                            os.listdir(os.path.join(organism_dir,bam))):
                        outputs = [os.path.join(organism_dir,bam,ff)
                                   for ff in ('qualimapReport.html',
                                              'rnaseq_qc_results.txt')]
                        if all([os.path.exists(ff) for ff in outputs]):
                            # All outputs present
                            organisms.add(d)
                            bam_files.add(bam)
                            output_files.extend(outputs)
                            # Add additional outputs (CSS, images etc)
                            for subdir in ('css',
                                           'images_qualimapReport',
                                           'raw_data_qualimapReport',):
                                dd = os.path.join(organism_dir,bam,subdir)
                                if os.path.exists(dd):
                                    extra_files = [os.path.join(dd,ff)
                                                   for ff in os.listdir(dd)]
                                    output_files.extend(extra_files)
                # Check for software information
                software = read_versions_file(
                    os.path.join(organism_dir,"_versions"),
                    software)
        if organisms:
            tags.add("qualimap_rnaseq")
            if not software:
                software['qualimap'] = [None]
        # Return collected information
        return AttributeDictionary(
            name='qualimap_rnaseq',
            software=software,
            bam_files=sorted(list(bam_files)),
            organisms=sorted(list(organisms)),
            output_files=output_files,
            tags=sorted(list(tags))
        )

    @classmethod
    def verify(self,params,qc_outputs):
        """
        Verify 'qualimap_rnaseq' QC module against outputs

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
        if not params.star_index or not params.annotation_gtf:
            # No STAR index or annotation
            return None
        if normalise_organism_name(params.organism) not in \
           qc_outputs.organisms:
            return False
        # Filter Fastq names and convert to BAM names
        bams = [get_bam_basename(fq)
                for fq in filter_fastqs(params.seq_data_reads[:1],
                                        params.seq_data_fastqs)]
        # Check that outputs exist for every BAM
        for bam in bams:
            if bam not in qc_outputs.bam_files:
                return False
        return True

    @classmethod
    def add_to_pipeline(self,p,project_name,qc_dir,bam_files,
                        gtf_annotation,organism_name,
                        rseqc_infer_experiment_outputs,
                        required_tasks=[],qualimap_runner=None):
        """
        Adds tasks for 'qualimap_rnaseq' module to pipeline

        Arguments:
          p (Pipeline): pipeline to extend
          project_name (str): name of project
          qc_dir (str): path to QC directory
          bam_files (list): BAM files to run the module on
          gtf_annotation (str): path to reference gene
            model GTF file
          organism_name (str): normalised name for organism
            that BAMs are aligned to
          rseqc_infer_experiment_outputs (dict): mapping of
            BAM names to dictionaries of associated RSeCQ
            'infer_experiment.py' outputs
          required_tasks (list): list of tasks that the module
            needs to wait for
          qualimap_runner (JobRunner): runner to use for Qualimap
        """
        # Run Qualimap RNA-seq analysis
        out_dir = os.path.join(qc_dir,
                               'qualimap-rnaseq',
                               organism_name)
        qualimap_rnaseq = RunQualimapRnaseq(
            "%s: Qualimap rnaseq: RNA-seq metrics" % project_name,
            bam_files,
            gtf_annotation,
            out_dir,
            bam_properties=rseqc_infer_experiment_outputs)
        p.add_task(qualimap_rnaseq,
                   requires=required_tasks,
                   runner=qualimap_runner)
        return qualimap_rnaseq

#######################################################################
# Pipeline tasks
#######################################################################

class RunQualimapRnaseq(PipelineTask):
    """
    Run Qualimap's 'rnaseq' module on BAM files

    Given a list of BAM files, for each file runs the
    Qualimap 'rnaseq' module
    (http://qualimap.conesalab.org/doc_html/command_line.html#rna-seq-qc)
    """
    def init(self,bam_files,feature_file,out_dir,
             bam_properties):
        """
        Initialise the RunQualimapRnaseq task

        Arguments:
          bam_files (list): list of paths to BAM files
            to run Qualimap rnaseq on
          feature_file (str): path to GTF file with the
            reference annotation data
          out_dir (str): path to a directory where the
            output files will be written
          bam_properties (mapping): properties for each
            BAM file from RSeQC 'infer_experiment.py'
            (used to determine if BAM is paired and
            what the strand-specificity is)
        """
        self.conda("qualimap=2.2")
        self.java_gc_threads = 1
        self.java_mem_size = '8G'
    def setup(self):
        # Check for feature file
        if self.args.feature_file:
            print("Feature file: %s" % self.args.feature_file)
        else:
            print("Feature file is not set, cannot run Qualimap rnaseq")
            return
        # Check for BAM file properties
        if self.args.bam_properties:
            self.bam_files = self.args.bam_files
        else:
            print("No properties for BAM files, cannot run Qualimap "
                  "rnaseq")
            return
        # Set up Qualimap rnaseq for each BAM file
        get_version = False
        for bam in self.args.bam_files:
            # Output directory for individual BAM file
            bam_name = os.path.basename(bam)[:-4]
            out_dir = os.path.join(self.args.out_dir,bam_name)
            # Check for existing outputs
            outputs_exist = True
            for f in qualimap_rnaseq_output(out_dir):
                outputs_exist = (outputs_exist and os.path.exists(f))
            if outputs_exist:
                # Skip running Qualimap for this BAM
                continue
            # Get properties for BAM file
            paired_end = self.args.bam_properties[bam]['paired_end']
            unstranded = self.args.bam_properties[bam]['unstranded']
            forward = self.args.bam_properties[bam]['forward']
            reverse = self.args.bam_properties[bam]['reverse']
            # Check data are valid
            if paired_end is None or \
               unstranded is None or \
               forward is None or \
               reverse is None:
                print("Bad properties data supplied for %s" %
                      os.path.basename(bam))
                for item in self.args.bam_properties[bam]:
                    print("-- %s: %s" % (item,
                                         self.args.bam_properties[bam][item]))
                print("Skipping Qualimap for this BAM file")
                continue
            # Set sequencing protocol (aka strand specificity)
            # Qualimap sequencing protocol can be
            # 'strand-specific-forward', 'strand-specific-reverse',
            # or 'non-strand-specific'
            if unstranded > forward and unstranded > reverse:
                seq_protocol = "non-strand-specific"
            else:
                ratio = forward/(reverse + 0.001)
                if ratio < 0.2:
                    seq_protocol = "strand-specific-reverse"
                elif ratio > 5:
                    seq_protocol = "strand-specific-forward"
                else:
                    seq_protocol = "non-strand-specific"
            print(f"-- {bam_name}: {seq_protocol} "
                  f"(F{forward:.1f}|R{reverse:.1f}|U{unstranded:.1f})")
            # Run Qualimap
            self.add_cmd("Run qualimap rnaseq on %s" %
                         os.path.basename(bam),
                         """
                         export _JAVA_OPTIONS="-XX:ParallelGCThreads={java_gc_threads} -Xmx{java_mem_size}"
                         if [ ! -z "$DISPLAY" ] ; then
                             echo "DISPLAY set to $DISPLAY"
                             echo "Unsetting to disable interactive window"
                             export DISPLAY=
                         fi
                         qualimap rnaseq \\
                             -bam {bam} \\
                             -gtf {feature_file} \\
                             -p {sequencing_protocol} \\
                             {paired} \\
                             -outdir {out_dir} \\
                             -outformat HTML \\
                             --java-mem-size={java_mem_size}
                         if [ $? -ne 0 ] ; then
                           echo "{bam}: qualimap rnaseq failed"
                           exit 1
                         fi
                         # Check expected outputs
                         for f in {outputs} ; do
                           if [ ! -e {out_dir}/$f ] ; then
                             echo "{bam}: missing output file $f"
                             exit 1
                           fi
                         done
                         # Copy outputs to final location
                         mkdir -p {final_dir}
                         cp -r {out_dir}/* {final_dir}
                         """.format(
                             bam=bam,
                             feature_file=self.args.feature_file,
                             final_dir=os.path.join(self.args.out_dir,bam_name),
                             sequencing_protocol=seq_protocol,
                             paired=('-pe' if paired_end else ''),
                             out_dir=bam_name,
                             outputs=' '.join(qualimap_rnaseq_output()),
                             nthreads=self.runner_nslots,
                             java_gc_threads=self.java_gc_threads,
                             java_mem_size=self.java_mem_size))
            get_version = True
        # Get version of Qualimap
        if get_version:
            self.add_cmd("Get Qualimap version",
                         """
                         qualimap --help >_versions 2>&1
                         """)
    def finish(self):
        missing_outputs = False
        if not self.args.feature_file:
            return
        if not self.args.bam_properties:
            return
        # Qualimap version
        if os.path.exists("_versions"):
            qualimap_version = None
            with open("_versions",'rt') as fp:
                for line in fp:
                    if line.startswith("QualiMap "):
                        # Example: QualiMap v.2.2.2-dev
                        qualimap_version = ' '.join(line.strip().split(' ')[1:])
                        break
            if qualimap_version:
                with open("_versions",'wt') as fp:
                    fp.write("qualimap\t%s\n" % qualimap_version)
            shutil.copy("_versions",self.args.out_dir)
        # Raise failure if errors were encountered
        if missing_outputs:
            self.fail(message="Some outputs are missing")
