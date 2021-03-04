#!/usr/bin/env python
#
#     barcodes.pipeline.py: pipelines for analysing barcodes
#     Copyright (C) University of Manchester 2019-2020 Peter Briggs
#

"""
Pipeline components for analysing barcodes.

Pipeline classes:

- AnalyseBarcodes

Pipeline task classes:

- SetupBarcodeAnalysisDirs
- CountBarcodes
- ListBarcodeCountFiles
- DetermineMismatches
- ReportBarcodeAnalysis
"""

######################################################################
# Imports
######################################################################

import os
import shutil
import tempfile
from ..analysis import AnalysisFastq
from ..applications import Command
from ..bcl2fastq.utils import get_nmismatches
from ..bcl2fastq.utils import bases_mask_is_valid
from ..bcl2fastq.utils import check_barcode_collisions
from ..pipeliner import Pipeline
from ..pipeliner import PipelineTask
from ..pipeliner import PipelineCommandWrapper
from ..pipeliner import PipelineParam
from ..pipeliner import PipelineFailure
from bcftbx.utils import AttributeDictionary
from bcftbx.FASTQFile import FastqIterator
from bcftbx.IlluminaData import IlluminaData
from bcftbx.IlluminaData import SampleSheet

######################################################################
# Pipeline classes
######################################################################

class AnalyseBarcodes(Pipeline):
    """
    Analyse the barcodes for Fastqs in a sequencing run

    Pipeline to perform barcode analysis on the Fastqs
    produced by bcl2fastq from a sequencing run.
    """
    def __init__(self,bcl2fastq_dir=None,sample_sheet=None):
        """
        Create a new AnalyseBarcodes pipeline instance

        At least one of the bcl2fastq output directory
        or sample sheet must be supplied when the
        pipeline is instantiated.

        If the bcl2fastq output directory is supplied
        on initialisation then it must exist and
        already contain output Fastq files.

        It is possible to set the pipeline up before the
        bcl2fastq outputs have been generated, as long
        as the sample sheet is supplied. The bcl2fastq
        output directory must then be supplied as an
        input when the pipeline is executed via the
        'run' method.

        Arguments:
          bcl2fastq_dir (str): path to the directory
            with outputs from bcl2fastq
          sample_sheet (str): path to the sample sheet
            file
        """
        # Initialise the pipeline superclass
        Pipeline.__init__(self,name="Analyse Barcodes")

        # Internal parameters
        self._bcl2fastq_dir = bcl2fastq_dir
        self._sample_sheet = sample_sheet

        # Define parameters
        self.add_param('bcl2fastq_dir',value=self._bcl2fastq_dir,type=str)
        self.add_param('sample_sheet',value=self._sample_sheet,type=str)
        self.add_param('barcode_analysis_dir',type=str)
        self.add_param('counts_dir',type=str)
        self.add_param('title',type=str)
        self.add_param('lanes',type=list)
        self.add_param('bases_mask',type=str)
        self.add_param('mismatches',type=int)
        self.add_param('cutoff',type=float,value=0.001)
        self.add_param('force',type=bool,value=False)

        # Get a list of projects
        if self._bcl2fastq_dir is not None:
            # Load data from bcl2fastq output
            try:
                analysis_dir = os.path.abspath(
                    os.path.dirname(self._bcl2fastq_dir))
                bcl2fastq_dir = os.path.basename(self._bcl2fastq_dir)
                illumina_data = IlluminaData(analysis_dir,
                                             unaligned_dir=bcl2fastq_dir)
            except Exception as ex:
                raise Exception("Unaligned dir '%s' supplied but can't "
                                "load data" % self._bcl2fastq_dir)
            # Get a list of projects
            projects = [p.name for p in illumina_data.projects]
        elif self._sample_sheet is not None:
            # Load data from sample sheet
            try:
                s = SampleSheet(self._sample_sheet)
                # List of unique project names
                projects = list(set(
                    [d[s.sample_project_column]
                     if d[s.sample_project_column]
                     else d[s.sample_id_column]
                     for d in s]))
            except Exception as ex:
                raise Exception("Sample sheet '%s' supplied but can't "
                                "get a list of project names")
        else:
            raise Exception("Need to supply either unaligned (bcl2fastq "
                            "output) dir or sample sheet")

        self.report("Expecting projects:")
        for p in projects:
            self.report("- %s" % p)

        ####################
        # Build the pipeline
        ####################

        # Setup barcode analysis and counts directories
        setup_barcode_analysis_dir = SetupBarcodeAnalysisDirs(
            "Setup barcode analysis directory",
            self.params.barcode_analysis_dir,
            self.params.counts_dir,
            force=self.params.force)
        self.add_task(setup_barcode_analysis_dir)

        # Load the data from the unaligned/bcl2fastq output dir
        load_illumina_data = LoadIlluminaData(
            "Load Fastq data for barcode analysis",
            self.params.bcl2fastq_dir)
        self.add_task(load_illumina_data)

        # Generate counts for each project
        count_tasks = []
        for project in projects:
            count_barcodes = CountBarcodes(
                "Count barcodes in '%s'" % project,
                load_illumina_data.output.illumina_data,
                project,
                self.params.counts_dir,
                lanes=self.params.lanes)
            self.add_task(count_barcodes,
                          requires=(setup_barcode_analysis_dir,
                                    load_illumina_data))
            count_tasks.append(count_barcodes)

        # Get counts for 'undetermined'
        count_barcodes = CountBarcodes(
            "Count barcodes in 'undetermined'",
            load_illumina_data.output.illumina_data,
            "__undetermined__",
            self.params.counts_dir,
            lanes=self.params.lanes,
            use_project_name="undetermined")
        self.add_task(count_barcodes,
                      requires=(setup_barcode_analysis_dir,
                                load_illumina_data))
        count_tasks.append(count_barcodes)

        # List the counts files
        list_counts_files = ListBarcodeCountFiles(
            "Fetch the barcode counts files",
            self.params.counts_dir)
        self.add_task(list_counts_files,
                      requires=count_tasks)

        # Analyse counts and report the results
        report_barcodes = ReportBarcodeAnalysis(
            "Report barcode analysis",
            list_counts_files.output.counts_files,
            self.params.barcode_analysis_dir,
            sample_sheet=self.params.sample_sheet,
            lanes=self.params.lanes,
            mismatches=self.params.mismatches,
            cutoff=self.params.cutoff,
            title=self.params.title
        )
        self.add_task(report_barcodes,
                      requires=(list_counts_files,))

        # Add final outputs to the pipeline
        self.add_output('report_file',report_barcodes.output.report_file)
        self.add_output('xls_file',report_barcodes.output.xls_file)
        self.add_output('html_file',report_barcodes.output.html_file)

    def run(self,barcode_analysis_dir,bcl2fastq_dir=None,title=None,
            lanes=None,mismatches=None,bases_mask=None,cutoff=None,
            sample_sheet=None,force=False,working_dir=None,log_file=None,
            batch_size=None,max_jobs=1,poll_interval=5,runner=None,
            verbose=False):
        """
        Run the tasks in the pipeline

        Arguments:
          barcode_analysis_dir (str): path to the directory
            to write the analysis results to
          bcl2fastq_dir (str): path to the bcl2fastq outputs
            (must be supplied here if a bcl2fastq output
            directory was not supplied on pipeline creation)
          title (str): optional, title for output reports
          lanes (list): optional, list of lanes to restrict
            the analysis to
          mismatches (int): optional, explicitly specify the
            number of mismatches to allow (default: determine
            number of mismatches automatically)
          bases_mask (str): optional, bases mask used for
            Fastq generation and demultiplexing
          cutoff (float): optional, don't report barcodes with
            a fraction of associated reads below this value
            (e.g. '0.001' excludes barcodes with < 0.1% of
            reads) (default: don't apply a cutoff)
          sample_sheet (str): optional, sample sheet to check
            barcode sequences against
          force (bool): if True then force regeneration of
            counts (default: re-use existing counts)
          working_dir (str): optional path to a working
            directory (defaults: temporary directory in the
            current directory)
          log_file (str): path to write log file to (default:
            don't write a log file)
          batch_size (int): if set then run commands in each
            task in batches, with each batch running this many
            commands at a time (default: run one command per
            job)
          max_jobs (int): optional maximum number of
            concurrent jobs in scheduler (default: 1)
          poll_interval (float): optional polling interval
            (seconds) to set in scheduler (default: 5s)
          runner (JobRunner): JobRunner instance to use to
             run jobs
          verbose (bool): if True then report additional
            information for diagnostics
        """
        # Working directory
        clean_up_on_completion = False
        if working_dir is None:
            working_dir = tempfile.mkdtemp(prefix="__barcode_analysis.",
                                           suffix=".tmp",
                                           dir=os.getcwd())
            clean_up_on_completion = True
        working_dir = os.path.abspath(working_dir)
        if not os.path.exists(working_dir):
            mkdir(working_dir)

        # Log and script directories
        log_dir = os.path.join(working_dir,"logs")
        scripts_dir = os.path.join(working_dir,"scripts")

        # Input bcl2fastq directory
        if self._bcl2fastq_dir is not None:
            if bcl2fastq_dir is None:
                bcl2fastq_dir = self._bcl2fastq_dir
            else:
                raise Exception("Bcl2fastq directory already set")
        else:
            if bcl2fastq_dir is None:
                raise Exception("No bcl2fastq output directory specified")
        bcl2fastq_dir = os.path.abspath(bcl2fastq_dir)

        # Input sample sheet
        if self._sample_sheet is not None:
            if sample_sheet is None:
                sample_sheet = self._sample_sheet
            else:
                raise Exception("Samplesheet already set")
        else:
            if sample_sheet is not None:
                sample_sheet = os.path.abspath(sample_sheet)

        # Barcode analysis and counts directories
        barcode_analysis_dir = os.path.abspath(barcode_analysis_dir)
        counts_dir = os.path.join(barcode_analysis_dir,"counts")

        # Execute the pipeline
        status = Pipeline.run(self,
                              working_dir=working_dir,
                              log_dir=log_dir,
                              scripts_dir=scripts_dir,
                              log_file=log_file,
                              batch_size=batch_size,
                              exit_on_failure=PipelineFailure.IMMEDIATE,
                              params={
                                  'bcl2fastq_dir': bcl2fastq_dir,
                                  'barcode_analysis_dir': barcode_analysis_dir,
                                  'counts_dir': counts_dir,
                                  'title': title,
                                  'lanes': lanes,
                                  'bases_mask': bases_mask,
                                  'mismatches': mismatches,
                                  'cutoff': cutoff,
                                  'sample_sheet': sample_sheet,
                                  'force': force,
                              },
                              max_jobs=max_jobs,
                              default_runner=runner,
                              finalize_outputs=True,
                              verbose=verbose)

        # Clean up working dir
        if status == 0 and clean_up_on_completion:
            shutil.rmtree(working_dir)

        # Return pipeline status
        return status

######################################################################
# Pipeline task classes
######################################################################

class SetupBarcodeAnalysisDirs(PipelineTask):
    """
    Set up the output directories
    """
    def init(self,barcode_analysis_dir,counts_dir,force=False):
        """
        Initialise the SetupBarcodeAnalysisDirs task

        Arguments:
          barcode_analysis_dir (str): final output directory
          counts_dir (str): directory to write counts files to
          force (bool): if True then remove existing counts
            files
        """
        pass
    def setup(self):
        # Create barcode analysis dir
        if not os.path.exists(self.args.barcode_analysis_dir):
            os.mkdir(self.args.barcode_analysis_dir)
        # Clean up existing counts for 'force'
        if os.path.exists(self.args.counts_dir):
            if self.args.force:
                self.report("Removing existing counts data")
                shutil.rmtree(self.args.counts_dir)
        # Create counts dir
        if not os.path.exists(self.args.counts_dir):
            os.mkdir(self.args.counts_dir)

class LoadIlluminaData(PipelineTask):
    """
    Load up an IlluminaData object from the bcl2fastq outputs
    """
    def init(self,bcl2fastq_dir):
        """
        """
        self.add_output('illumina_data',PipelineParam())
    def setup(self):
        # Load the data from bcl2fastq
        print("Loading data from %s" % self.args.bcl2fastq_dir)
        illumina_data = IlluminaData(os.path.dirname(self.args.bcl2fastq_dir),
                                     os.path.basename(self.args.bcl2fastq_dir))
        self.output.illumina_data.set(illumina_data)

class CountBarcodes(PipelineTask):
    """
    Generate barcode counts for a project
    """
    def init(self,illumina_data,project,counts_dir,lanes=None,
             use_project_name=None):
        """
        Initialise the CountBarcodes task

        Arguments:
          illumina_data (IlluminaData): IlluminaData object
            for the Fastq data that should be examined
          project (str): name of project with Fastqs
            to get barcode codes from
          counts_dir (str): directory to write counts
            files to
          lanes (list): optional list of lanes to restrict
            counts generation to
          use_project_name (str): optional alternative name
            for project to use in counts file names (default
            is to use the name from the supplied project)
        """
        pass
    def setup(self):
        # Project name to use
        project_name = self.args.project
        # Acquire the matching IlluminaProject
        project = None
        if project_name != "__undetermined__":
            # Standard projects
            for p in self.args.illumina_data.projects:
                if p.name == project_name:
                    project = p
                    break
        else:
            # Undetermined reads
            project = self.args.illumina_data.undetermined
        # Check project was located
        # If not then exit without an error
        if not project:
            print("Unable to locate project '%s': nothing to do" %
                  project_name)
            return
        # Collect Fastq files
        fastqs = []
        for sample in project.samples:
            for fastq in sample.fastq:
                fq = AnalysisFastq(fastq)
                if fq.is_index_read or fq.read_number != 1:
                    # Only collect R1 Fastqs
                    continue
                if self.args.lanes:
                    if fq.lane_number is not None and \
                       fq.lane_number not in self.args.lanes:
                        # Only collect Fastqs explicitly belonging
                        # to the specified lanes
                        continue
                # Include this Fastq
                fastqs.append(os.path.join(sample.dirn,fastq))
        # Set up counting for each Fastq
        for fastq in fastqs:
            fq = os.path.basename(fastq)
            counts_file = os.path.join(self.args.counts_dir,
                                       "%s.%s.counts" %
                                       (project_name,fq))
            if os.path.exists(counts_file):
                if (os.path.getmtime(counts_file) >
                    os.path.getmtime(fastq)):
                    # Counts file already exists and is newer
                    # than Fastq
                    continue
            # Build count command
            cmd = PipelineCommandWrapper(
                "Run analyse_barcodes.py -c for %s" % fq,
                'analyse_barcodes.py',
                '-o',counts_file,
                '--no-report',
                fastq)
            self.add_cmd(cmd)

class ListBarcodeCountFiles(PipelineTask):
    """
    Collect counts files from a directory
    """
    def init(self,counts_dir):
        """
        Initialise the ListBarcodeCountFiles task

        Arguments:
          counts_dir (str): directory holding the
            counts files

        Outputs:
          counts_files (list): list of counts
            files
        """
        self.add_output('counts_files',list())
    def setup(self):
        for f in os.listdir(self.args.counts_dir):
            if f.endswith(".counts"):
                self.output.counts_files.append(
                    os.path.join(self.args.counts_dir,f))

class ReportBarcodeAnalysis(PipelineTask):
    """
    Perform analysis and reporting of barcode counts
    """
    def init(self,counts_files,barcode_analysis_dir,lanes=None,
             mismatches=None,cutoff=0.001,sample_sheet=None,
             title=None):
        """
        Initialise the ReportBarcodeAnalysis task

        Arguments:
          counts_files (list): counts files to include
            in the analysis
          barcode_analysis_dir (str): path to the directory
            to write the analysis results to
          lanes (list): optional list of lanes to restrict
            the analysis to
          mismatches (int): optional number of mismatches
            to allow when comparing barcodes
          cutoff (float): optional fraction of total barcodes
            below which indexes won't be reported
          sample_sheet (str): optional, sample sheet to check
            barcode sequences against
          title (str): optional, title for output reports

        Outputs:
          report_file (str): path to the report file
          xls_file (str): path to the XLS report
          html_file (str): path to the HTML report
        """
        self.add_output('report_file',PipelineParam(type=str))
        self.add_output('xls_file',PipelineParam(type=str))
        self.add_output('html_file',PipelineParam(type=str))
    def setup(self):
        # Make output filenames
        report_file = os.path.join(self.args.barcode_analysis_dir,
                                   'barcodes.report')
        xls_file = os.path.join(self.args.barcode_analysis_dir,
                                'barcodes.xls')
        html_file = os.path.join(self.args.barcode_analysis_dir,
                                 'barcodes.html')
        # Remove existing copies, if found
        for filen in (report_file,xls_file,html_file):
            if os.path.exists(filen):
                os.remove(filen)
        # Build command to run the barcode analysis
        cmd = PipelineCommandWrapper(
            "Run analyse_barcodes.py to report barcodes",
            'analyse_barcodes.py',
            '--report',report_file,
            '--xls',xls_file,
            '--html',html_file)
        if self.args.sample_sheet:
            cmd.add_args('--sample-sheet',self.args.sample_sheet)
        if self.args.lanes:
            lanes = self.args.lanes
        elif self.args.sample_sheet:
            # Implicitly get lanes from sample sheet
            try:
                lanes = sorted(
                    set([line['Lane']
                         for line in SampleSheet(self.args.sample_sheet)]))
            except KeyError:
                # No lanes
                lanes = None
        else:
            lanes = None
        if lanes:
            cmd.add_args('--lanes',
                         ','.join([str(l) for l in lanes]))
        if self.args.cutoff:
            cmd.add_args('--cutoff',self.args.cutoff)
        if self.args.mismatches:
            cmd.add_args('--mismatches',self.args.mismatches)
        if self.args.title:
            cmd.add_args('--title',self.args.title)
        cmd.add_args('-c')
        cmd.add_args(*self.args.counts_files)
        self.add_cmd(cmd)
        # Update the output parameters
        self.output.report_file.set(report_file)
        self.output.xls_file.set(xls_file)
        self.output.html_file.set(html_file)
    def finish(self):
        # Check the output files exist
        report_file = self.output.report_file.value
        if not os.path.exists(report_file):
            self.fail(message="Missing report file: %s" % report_file)
        xls_file = self.output.xls_file.value
        if not os.path.exists(xls_file):
            self.fail(message="Missing XLS report: %s" % xls_file)
        html_file = self.output.html_file.value
        if not os.path.exists(html_file):
            self.fail(message="Missing HTML report: %s" % html_file)
