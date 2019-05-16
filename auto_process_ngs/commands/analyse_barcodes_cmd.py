#!/usr/bin/env python
#
#     analyse_barcodes_cmd.py: implement analyse_barcodes command
#     Copyright (C) University of Manchester 2019 Peter Briggs
#
#########################################################################

#######################################################################
# Imports
#######################################################################

import os
import shutil
import logging
import tempfile
from ..analysis import AnalysisFastq
from ..applications import Command
from ..bcl2fastq_utils import get_nmismatches
from ..bcl2fastq_utils import bases_mask_is_valid
from ..pipeliner import Pipeline
from ..pipeliner import PipelineTask
from ..pipeliner import PipelineCommandWrapper
from ..pipeliner import PipelineParam
from ..pipeliner import PipelineFailure
from ..simple_scheduler import SimpleScheduler
from ..utils import write_script_file
from bcftbx.utils import AttributeDictionary
from bcftbx.FASTQFile import FastqIterator
from bcftbx import IlluminaData

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Classes
#######################################################################

class BarcodeAnalysis(Pipeline):
    """
    """
    def __init__(self,unaligned_dir=None):

        # Initialise the pipeline superclass
        Pipeline.__init__(self)

        # Define parameters
        self.add_param('barcode_analysis_dir',type=str)
        self.add_param('counts_dir',type=str)
        self.add_param('lanes',type=list)
        self.add_param('sample_sheet',type=str)
        self.add_param('bases_mask',type=str)
        self.add_param('mismatches',type=int)
        self.add_param('cutoff',type=float)
        self.add_param('force',type=bool,value=False)

        # Output
        self._output = AttributeDictionary()

        # Load data from bcl2fastq output
        if not os.path.exists(unaligned_dir):
            raise OSError("'%s': not found" % unaligned_dir)
        analysis_dir = os.path.abspath(os.path.dirname(unaligned_dir))
        unaligned_dir = os.path.basename(unaligned_dir)
        illumina_data = IlluminaData.IlluminaData(analysis_dir,
                                                  unaligned_dir=unaligned_dir)

        # Example Fastq file used for determining mismatches in
        # absence of bases mask
        example_fastq = illumina_data.projects[0].samples[0].fastq_subset(
            read_number=1,full_path=True)[0]

        ####################
        # Build the pipeline
        ####################

        # Setup barcode analysis and counts directories
        setup_barcode_analysis_dir = SetupBarcodeAnalysisDir(
            "Setup barcode analysis directory",
            self.params.barcode_analysis_dir,
            self.params.counts_dir,
            force=self.params.force)
        self.add_task(setup_barcode_analysis_dir)

        # Generate counts for Fastqs in each project
        count_tasks = []
        for project in illumina_data.projects:
            count_barcodes = CountBarcodes(
                "Count barcodes in '%s'" % project.name,
                project,
                self.params.counts_dir,
                lanes=self.params.lanes)
            self.add_task(count_barcodes,
                          requires=(setup_barcode_analysis_dir,))
            count_tasks.append(count_barcodes)

        # Generate counts for undetermined Fastqs
        if illumina_data.undetermined is not None:
            count_barcodes = CountBarcodes(
                "Count barcodes in 'undetermined'",
                illumina_data.undetermined,
                self.params.counts_dir,
                lanes=self.params.lanes,
                use_project_name="undetermined")
            self.add_task(count_barcodes,
                          requires=(setup_barcode_analysis_dir,))
            count_tasks.append(count_barcodes)

        # List the counts files
        list_counts_files = ListBarcodeCountFiles(
            "Fetch the barcode counts files",
            self.params.counts_dir,
            lanes=self.params.lanes)
        self.add_task(list_counts_files,
                      requires=count_tasks)

        # Determine mismatches
        get_n_mismatches = DetermineMismatches(
            "Determine mismatches",
            self.params.mismatches,
            self.params.bases_mask,
            example_fastq)
        self.add_task(get_n_mismatches)

        # Analyse counts and report the results
        report_barcodes = ReportBarcodeAnalysis(
            "Reporting barcode analysis",
            list_counts_files.output.counts_files,
            self.params.barcode_analysis_dir,
            sample_sheet=self.params.sample_sheet,
            lanes=self.params.lanes,
            mismatches=get_n_mismatches.output.mismatches,
            cutoff=self.params.cutoff,
        )
        self.add_task(report_barcodes,
                      requires=(list_counts_files,
                                get_n_mismatches))

        # Add final outputs to the pipeline
        self.add_output('report_file',report_barcodes.output.report_file)
        self.add_output('xls_file',report_barcodes.output.xls_file)
        self.add_output('html_file',report_barcodes.output.html_file)

    def add_output(self,name,value):
        """
        Add an output to the pipeline

        Arguments:
          name (str): name for the output
          value (object): associated object
        """
        self._output[name] = value

    @property
    def output(self):
        """
        Return the output object
        """
        return self._output

    def run(self,barcode_analysis_dir,lanes=None,mismatches=None,
            bases_mask=None,cutoff=None,sample_sheet=None,force=False,
            working_dir=None,log_file=None,batch_size=None,max_jobs=1,
            poll_interval=5,runner=None,verbose=False):
        """
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
                                  'barcode_analysis_dir': barcode_analysis_dir,
                                  'counts_dir': counts_dir,
                                  'lanes': lanes,
                                  'bases_mask': bases_mask,
                                  'mismatches': mismatches,
                                  'cutoff': cutoff,
                                  'sample_sheet': sample_sheet,
                                  'force': force,
                              },
                              max_jobs=max_jobs,
                              default_runner=runner,
                              verbose=verbose)

        # Clean up working dir
        if status == 0 and clean_up_on_completion:
            shutil.rmtree(working_dir)

        # Update the outputs
        for name in self._output:
            try:
                self._output[name] = self._output[name].value
            except AttributeError:
                pass

        # Return pipeline status
        return status

class SetupBarcodeAnalysisDir(PipelineTask):
    def init(self,barcode_analysis_dir,counts_dir,force=False):
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

class CountBarcodes(PipelineTask):
    def init(self,project,counts_dir,lanes=None,use_project_name=None):
        pass
    def setup(self):
        # Project name to use
        project_name = self.args.use_project_name
        if project_name is None:
            project_name = self.args.project.name
        # Collect Fastq files
        fastqs = []
        for sample in self.args.project.samples:
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
                # Already exists
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
    def init(self,counts_dir,lanes=None):
        self.add_output('counts_files',list())
    def setup(self):
        for f in os.listdir(self.args.counts_dir):
            if f.endswith(".counts"):
                self.output.counts_files.append(
                    os.path.join(self.args.counts_dir,f))

class DetermineMismatches(PipelineTask):
    def init(self,mismatches,bases_mask,example_fastq):
        self.add_output('mismatches',PipelineParam(type=int))
    def setup(self):
        if self.args.mismatches is not None:
            # Use the supplied value
            mismatches = self.args.mismatches
        else:
            # Try to determine number of mismatches from bases mask
            bases_mask = self.args.bases_mask
            if bases_mask_is_valid(bases_mask):
                mismatches = get_nmismatches(bases_mask)
            else:
                # Not a valid bases mask
                # Extract from example Fastq file
                for r in FastqIterator(self.args.example_fastq):
                    seq_id = r.seqid
                    break
                if len(seq_id.index_sequence) >= 6:
                    mismatches = 1
                else:
                    mismatches = 0
        # Set the mismatches output value
        self.output.mismatches.set(self.args.mismatches)

class ReportBarcodeAnalysis(PipelineTask):
    def init(self,counts_files,barcode_analysis_dir,lanes=None,
             mismatches=None,cutoff=None,sample_sheet=None):
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
                         for line in
                         IlluminaData.SampleSheet(
                             self.args.sample_sheet).data]))
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

#######################################################################
# Command functions
#######################################################################

def analyse_barcodes(ap,unaligned_dir=None,lanes=None,
                     mismatches=None,cutoff=None,
                     barcode_analysis_dir=None,
                     sample_sheet=None,runner=None,
                     force=False):
    """Analyse the barcode sequences for Fastqs for each specified lane

    Run 'analyse_barcodes.py' for one or more lanes, to analyse the
    barcode index sequences in each lane.
    
    Arguments:
      ap (AutoProcessor): autoprocessor pointing to the analysis
        directory to create Fastqs for
      unaligned_dir (str): if set then use this as the output directory
        for bcl-to-fastq conversion. Default is 'bcl2fastq' (unless
        an alternative is already specified in the config file)
      lanes (list): (optional) specify a list of lane numbers to
        use in the processing; lanes not in the list will be excluded
        (default is to include all lanes)
      mismatches (int): (optional) maximum number of mismatches to
        consider when grouping similar barcodes; default is to determine
        it automatically from the bases mask
      cutoff (float): (optional) exclude barcodes with a smaller fraction
        of associated reads than specified cutoff from reporting (e.g.
        '0.001' excludes barcodes with < 0.1% of reads); default is to
        include all barcodes
      sample_sheet (str): if set then use this as the input samplesheet
        to check barcode sequences against (by default will use the
        sample sheet defined in the parameter file for the run)
      barcode_analysis_dir (str): (optional) explicitly specify the
        subdirectory to use for barcode analysis. Counts will be
        written to and read from the 'counts' subdirectory of this
        directory (defaults to 'barcode_analysis')
      runner (JobRunner): (optional) specify a non-default job runner
        to use for barcode analysis
      force (bool): if True then forces regeneration of any existing
        counts (default is to reuse existing counts)
    """
    # Sort out parameters
    if unaligned_dir is not None:
        ap.params['unaligned_dir'] = unaligned_dir
    elif ap.params['unaligned_dir'] is None:
        ap.params['unaligned_dir'] = 'bcl2fastq'
    bases_mask = ap.params['bases_mask']
    if sample_sheet is None:
        sample_sheet = ap.params.sample_sheet
    # Handle barcode analysis subdirectory
    if barcode_analysis_dir is not None:
        # Create a subdirectory for barcode analysis
        ap.params['barcode_analysis_dir'] = barcode_analysis_dir
    elif ap.params['barcode_analysis_dir'] is None:
        ap.params['barcode_analysis_dir'] = 'barcode_analysis'
    barcode_analysis_dir = ap.params['barcode_analysis_dir']
    if not os.path.isabs(barcode_analysis_dir):
        barcode_analysis_dir = os.path.join(ap.params['analysis_dir'],
                                            barcode_analysis_dir)
    # Create a pipeline for barcode analysis
    barcode_analysis = BarcodeAnalysis(os.path.join(
        ap.params['analysis_dir'],ap.params['unaligned_dir']))
    # Log dir and log file
    ap.set_log_dir(ap.get_log_subdir('analyse_barcodes'))
    log_file = os.path.join(ap.log_dir,"analyse_barcodes.log")
    # Set up runner
    if runner is None:
        runner = ap.settings.general.default_runner
    runner.set_log_dir(ap.log_dir)
    # Run the pipeline
    status = barcode_analysis.run(
        barcode_analysis_dir,
        lanes=lanes,
        mismatches=mismatches,
        cutoff=cutoff,
        sample_sheet=sample_sheet,
        force=force,
        log_file=log_file,
        runner=runner,
        verbose=False)
    # Finish
    if status == 0:
        print "Report written to %s" % barcode_analysis.output.report_file
        print "XLS written to %s" % barcode_analysis.output.xls_file
        print "HTML written to %s" % barcode_analysis.output.html_file
    return status
