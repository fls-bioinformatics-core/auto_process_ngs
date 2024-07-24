#!/usr/bin/env python3
#
#     picard_insert_size_metrics: implements 'picard_insert_size_metrics' QC module
#     Copyright (C) University of Manchester 2024 Peter Briggs

"""
Implements the 'picard_insert_size_metrics' QC module:

* PicardInsertSizeMetrics: core QCModule class
* RunPicardCollectInsertSizeMetrics: pipeline task to run Picard
  'CollectInsertSizeMetrics'
* CollateInsertSizes: pipeline task to collate insert sizes
"""

#######################################################################
# Imports
#######################################################################

import os
import shutil
import logging
from bcftbx.TabFile import TabFile
from bcftbx.utils import AttributeDictionary
from . import QCModule
from ..picard import CollectInsertSizeMetrics
from ..picard import picard_collect_insert_size_metrics_output
from ..utils import filter_fastqs
from ..utils import get_bam_basename
from ..utils import read_versions_file
from ...utils import normalise_organism_name
from ...pipeliner import PipelineTask

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Core class
#######################################################################

class PicardInsertSizeMetrics(QCModule):
    """
    Class for handling the 'picard_insert_size_metrics' QC module
    """
    name = "picard_insert_size_metrics"
    mapped_metrics = True
    require_bam_files = True
    runners = ("picard_runner",)
    envmodules = ()
    
    def __init__(self):
        QCModule.__init__(self)

    @classmethod
    def collect_qc_outputs(self,qc_dir):
        """
        Collect information on picard_insert_size_metrics outputs

        Returns an AttributeDictionary with the following
        attributes:

        - name: set to 'picard_collect_insert_size_metrics'
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
        # Look for Picard CollectInsertSizeMetrics outputs
        organisms = set()
        picard_dir = os.path.join(qc_dir.path,"picard")
        if os.path.isdir(picard_dir):
            # Look for subdirs with organism names
            for d in filter(
                    lambda dd:
                    os.path.isdir(os.path.join(picard_dir,dd)),
                    os.listdir(picard_dir)):
                # Check for outputs
                for f in filter(
                        lambda ff:
                        ff.endswith(".insert_size_metrics.txt"),
                        os.listdir(os.path.join(picard_dir,d))):
                    name = f[:-len(".insert_size_metrics.txt")]
                    outputs = picard_collect_insert_size_metrics_output(
                        name,
                        prefix=os.path.join(picard_dir,d))
                    if all([os.path.exists(f) for f in outputs]):
                        # All outputs present
                        organisms.add(d)
                        bam_files.add(name)
                        output_files.extend(outputs)
                # Check for software information
                software = read_versions_file(
                    os.path.join(picard_dir,d,"_versions"),
                    software)
        if organisms:
            tags.add("picard_insert_size_metrics")
            if not software:
                software['picard'] = [None]
        # Look for collated insert sizes files
        for f in filter(
                lambda ff:
                os.path.basename(ff).startswith("insert_sizes.") and
                os.path.basename(ff).endswith(".tsv"),
                qc_dir.file_list):
            tags.add("collated_insert_sizes")
            output_files.append(f)
            organisms.add(os.path.basename(f)\
                          [len("insert_sizes."):-len(".tsv")])
        # Return collected information
        return AttributeDictionary(
            name=self.name,
            software=software,
            bam_files=sorted(list(bam_files)),
            organisms=sorted(list(organisms)),
            output_files=output_files,
            tags=sorted(list(tags))
        )

    @classmethod
    def verify(self,params,qc_outputs):
        """
        Verify 'picard_insert_size_metrics' QC module against outputs

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
        if not params.star_index:
            # No STAR index
            return None
        if normalise_organism_name(params.organism) not in qc_outputs.organisms:
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
    def add_to_pipeline(self,p,project_name,project,qc_dir,
                        bam_files,organism_name,
                        required_tasks=[],compute_runner=None):
        """
        Adds tasks for 'picard_insert_size_metrics' module to pipeline

        Arguments:
          p (Pipeline): pipeline to extend
          project_name (str): name of project
          project (AnalysisProject): project to run module on
          qc_dir (str): path to QC directory
          bam_files (list): BAM files to run the module on
          organism_name (str): normalised name for organism
            that BAMs are aligned to
          required_tasks (list): list of tasks that the module
            needs to wait for
          compute_runner (JobRunner): runner to use for
            computation
        """
        # Run Picard's CollectInsertSizeMetrics

        out_dir = os.path.join(qc_dir,'picard',organism_name)

        insert_size_metrics = RunPicardCollectInsertSizeMetrics(
            "%s: Picard: collect insert size metrics" %
            project_name,
            bam_files,
            out_dir)
        p.add_task(insert_size_metrics,
                   requires=required_tasks,
                   runner=compute_runner)

        collate_insert_sizes = CollateInsertSizes(
            "%s: collate insert size data" % project_name,
            bam_files,
            out_dir,
            os.path.join(qc_dir,
                         'insert_sizes.%s.tsv' % organism_name))
        p.add_task(collate_insert_sizes,
                   requires=(insert_size_metrics,))
        return collate_insert_sizes

#######################################################################
# Pipeline tasks
#######################################################################

class RunPicardCollectInsertSizeMetrics(PipelineTask):
    """
    Run Picard 'CollectInsertSizeMetrics' on BAM files

    Given a list of BAM files, for each file first runs
    the Picard 'CleanSam' utility (to remove alignments
    that would otherwise cause problems for the insert
    size calculations) and then 'CollectInsertSizeMetrics'
    to generate the insert size metrics.

    Note that this task should only be run on BAM files
    with paired-end data.
    """
    def init(self,bam_files,out_dir):
        """
        Initialise the RunPicardCollectInsertSizeMetrics
        task

        Arguments:
          bam_files (list): list of paths to BAM files
            to run CollectInsertSizeMetrics on
          out_dir (str): path to a directory where the
            output files will be written
        """
        # Conda dependencies
        self.conda("picard=2.27.1",
                   "r-base=4")
        self.java_gc_threads = 1
        self.java_mem_size = '4G'
    def setup(self):
        # Set up commands to run CleanSam and
        # CollectInsertSizeMetrics for each BAM file
        get_version = False
        for bam in self.args.bam_files:
            # Check if outputs already exist
            outputs_exist = True
            for f in picard_collect_insert_size_metrics_output(
                    bam,
                    self.args.out_dir):
                outputs_exist = (outputs_exist and os.path.exists(f))
            if outputs_exist:
                # Skip this BAM
                continue
            # Add command to get insert sizes
            self.add_cmd("%s: collect insert size metrics" %
                         os.path.basename(bam),
                         """
                         export _JAVA_OPTIONS="-XX:ParallelGCThreads={java_gc_threads} -Xmx{java_mem_size}"
                         tmpdir=$(mktemp -d)
                         picard CleanSam \\
                             -I {bam} \\
                             -O $tmpdir/{basename}.bam \\
                             -XX:ActiveProcessorCount={nslots} \\
                         && \\
                         picard CollectInsertSizeMetrics \\
                             -I $tmpdir/{basename}.bam \\
                             -O {basename}.insert_size_metrics.txt \\
                             -H {basename}.insert_size_histogram.pdf \\
                             -XX:ActiveProcessorCount={nslots}
                         """.format(bam=bam,
                                    basename=os.path.basename(bam)[:-4],
                                    nslots=self.runner_nslots,
                                    java_gc_threads=self.java_gc_threads,
                                    java_mem_size=self.java_mem_size))
            get_version = True
        # Get version of Picard
        if get_version:
            self.add_cmd("Get Picard version",
                         """
                         # Get version of CollectInsertSizeMetrics
                         picard CollectInsertSizeMetrics --version >_versions 2>&1
                         # Force zero exit code
                         exit 0
                         """)
    def finish(self):
        # Check if any BAM files were processed
        if not self.args.bam_files:
            return
        # Copy outputs to final location
        if not os.path.exists(self.args.out_dir):
            print("Creating output dir '%s'" % self.args.out_dir)
            os.makedirs(self.args.out_dir)
        for bam in self.args.bam_files:
            for f in picard_collect_insert_size_metrics_output(
                    bam,
                    self.args.out_dir):
                if not os.path.exists(f):
                    # Copy new version to ouput location
                    shutil.copy(os.path.basename(f),self.args.out_dir)
        # Picard version
        if os.path.exists("_versions"):
            picard_version = None
            with open("_versions",'rt') as fp:
                for line in fp:
                    if line.startswith("Version:"):
                        # Example: Version:2.27.1
                        picard_version = ':'.join(line.strip().split(':')[1:])
                        break
            if picard_version:
                with open("_versions",'wt') as fp:
                    fp.write("picard\t%s\n" % picard_version)
            shutil.copy("_versions",self.args.out_dir)

class CollateInsertSizes(PipelineTask):
    """
    Collate insert size metrics data from multiple BAMs

    Gathers together the Picard insert size data from a
    set of BAM files and puts them into a single TSV
    file.
    """
    def init(self,bam_files,picard_out_dir,out_file,delimiter='\t'):
        """
        Initialise the CollateInsertSizes task

        Arguments:
          bam_files (list): list of paths to BAM files
            to get associated insert size data for
          picard_out_dir (str): path to the directory
            containing the Picard CollectInsertSizeMetrics
            output files
          out_file (str): path to the output TSV file
          delimiter (str): specify the delimiter to use
            in the output file
        """
        pass
    def setup(self):
        # Set up a TabFile instance for the collated data
        tf = TabFile(column_names=("Bam file",
                                   "Mean insert size",
                                   "Standard deviation",
                                   "Median insert size",
                                   "Median absolute deviation"))
        metrics_files = {}
        for bam in self.args.bam_files:
            # Get metrics file associated with this BAM file
            outputs = list(filter(lambda f:
                                  f.endswith('.txt') and
                                  os.path.exists(f),
                                  picard_collect_insert_size_metrics_output(
                                      os.path.basename(bam)[:-4],
                                      prefix=self.args.picard_out_dir)))
            if not outputs:
                # No metrics located
                print("%s: no associated Picard insert size "
                      "metrics file found in %s" %
                      (bam,
                       self.args.picard_out_dir))
                continue
            metrics_files[os.path.basename(bam)] = outputs[0]
        # Check there is data to collate
        if not metrics_files:
            print("no insert size metrics files recovered")
            return
        # Set up a TabFile instance for the collated data
        tf = TabFile(column_names=("Bam file",
                                   "Mean insert size",
                                   "Standard deviation",
                                   "Median insert size",
                                   "Median absolute deviation"))
        for bam in sorted(list(metrics_files.keys())):
            # Get mean and median insert sizes
            metrics_file = metrics_files[bam]
            insert_size_metrics = CollectInsertSizeMetrics(metrics_file)
            tf.append(data=(
                bam,
                insert_size_metrics.metrics['MEAN_INSERT_SIZE'],
                insert_size_metrics.metrics['STANDARD_DEVIATION'],
                insert_size_metrics.metrics['MEDIAN_INSERT_SIZE'],
                insert_size_metrics.metrics['MEDIAN_ABSOLUTE_DEVIATION']
            ))
        # Output to file
        print("Writing to %s" % self.args.out_file)
        tf.write(self.args.out_file,
                 include_header=True,
                 delimiter=self.args.delimiter)
