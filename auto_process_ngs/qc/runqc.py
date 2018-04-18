#!/usr/bin/env python
#
#     runqc: run and report QC from analysis projects
#     Copyright (C) University of Manchester 2018 Peter Briggs
#

#######################################################################
# Imports
#######################################################################

import os
import logging
import auto_process_ngs.utils as utils
import auto_process_ngs.fileops as fileops
from auto_process_ngs.applications import Command
from auto_process_ngs.settings import Settings
from auto_process_ngs.simple_scheduler import SimpleScheduler

# Module-specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Classes
#######################################################################

class RunQC(object):
    """
    Class for running QC across multiple projects

    Usage example:

    >>> # Set up runner
    >>> runqc = RunQC()
    >>> # Add projects
    >>> for project in project_list:
    >>>    runqc.add_project(project)
    >>> # Execute and get status
    >>> status = runqc.run()
    """
    def __init__(self):
        """
        Create a new RunQC instance
        """
        self._settings = Settings()
        self._projects = []
        self._sched = None

    def add_project(self,project,fastq_dir=None,qc_dir=None,
                    sample_pattern=None,ungzip_fastqs=False):
        """
        Add a project to run the QC for

        Arguments:
          project (AnalysisProject): analysis project
            to run the QC for
          fastq_dir (str): optional, specify the subdir
            with the Fastqs to be be used
          sample_pattern (str): optional, specify a
            glob-style pattern to use to select a subset
            of samples
          qc_dir (str): optional, specify the subdir to
            write the QC outputs to
          ungzip_fastqs (bool): if True then uncompress
            source Fastqs (default: False i.e. don't
            uncompress the Fastqs)
        """
        self._projects.append(ProjectQC(project,
                                        fastq_dir=fastq_dir,
                                        sample_pattern=sample_pattern,
                                        qc_dir=qc_dir,
                                        ungzip_fastqs=ungzip_fastqs))

    def run(self,nthreads=1,fastq_screen_subset=10000,
            report_html=None,multiqc=False,
            qc_runner=None,verify_runner=None,
            report_runner=None,max_jobs=None):
        """
        Schedule and execute QC jobs

        Arguments:
          nthreads (int): the maxiumum number of
            threads/cores to use per job (default: 1)
          fastq_screen_subset (int): the size of
            subset to use with FastQScreen (default:
            10000)
          report_html (str): optional, path to the name of
            the QC report
          multiqc (bool): if True then also run MultiQC
            when generating reports
          qc_runner (JobRunner): job runner to use for
            executing QC
          verify_runner (JobRunner): job runner to use
            for QC verification
          report_runner (JobRunner): job runner to use
            for QC reporting
          max_jobs (int): optional, specify maximum
            number of jobs to run concurrently

        Returns:
          Integer: returns 0 if QC ran to completion
            without problems, non-zero if there was
            an error.
        """
        # Sort out runners
        if qc_runner is None:
            qc_runner = self._settings.runners.qc
        if verify_runner is None:
            verify_runner = \
                self._settings.general.default_runner
        if report_runner is None:
            report_runner = verify_runner
        # Set up and start scheduler
        self._sched = SimpleScheduler(max_concurrent=max_jobs)
        self._sched.start()
        # Initial QC check for each project
        for project in self._projects:
            project.check_qc(self._sched,
                             name="pre_qc_check",
                             runner=verify_runner)
        self._sched.wait()
        # Run QC for each project
        for project in self._projects:
            project.setup_qc(self._sched,
                             nthreads,
                             fastq_screen_subset,
                             qc_runner=qc_runner,
                             verify_runner=verify_runner)
        self._sched.wait()
        # Verify the outputs and generate QC reports
        failed_projects = []
        for project in self._projects:
            if not project.verify():
                failed_projects.append(project)
            else:
                project.report_qc(self._sched,
                                  report_html=report_html,
                                  multiqc=multiqc,
                                  runner=report_runner)
        self._sched.wait()
        self._sched.stop()
        # Check reporting
        for project in self._projects:
            if project not in failed_projects:
                if project.reporting_status == 0:
                    print "Generated QC report for '%s'" % project.name
                else:
                    failed_projects.append(project)
        # Report failed projects
        if failed_projects:
            logger.error("QC failed for one or more samples in the "
                         "following projects:")
            for project in failed_projects:
                logger.error("- %s" % project.name)
            return 1
        # Finish
        return 0

class ProjectQC(object):
    """
    Class for setting up QC jobs for a project
    """
    def __init__(self,project,fastq_dir=None,sample_pattern=None,
                 qc_dir=None,ungzip_fastqs=False):
        """
        Create a new ProjectQC instance

        Arguments:
          project (AnalysisProject): analysis project
            to run the QC for
          fastq_dir (str): optional, specify the subdir
            with the Fastqs to be be used
          sample_pattern (str): optional, specify a
            glob-style pattern to use to select a subset
            of samples
          qc_dir (str): optional, specify the subdir to
            write the QC outputs to
          ungzip_fastqs (bool): if True then uncompress
            source Fastqs (default: False i.e. don't
            uncompress the Fastqs)
        """
        # Clone the supplied project
        self.project = utils.AnalysisProject(project.name,
                                             project.dirn)
        # Unpick the supplied subdirectories
        if qc_dir is None:
            qc_dir = 'qc'
        self.qc_dir = project.setup_qc_dir(qc_dir,
                                           fastq_dir=fastq_dir)
        project.use_qc_dir(self.qc_dir)
        self.fastq_dir = project.qc_info(project.qc_dir).fastq_dir
        project.use_fastq_dir(self.fastq_dir)
        # Log directory
        self.log_dir = os.path.join(project.qc_dir,'logs')
        if not os.path.exists(self.log_dir):
            print "Making QC logs directory: %s" % self.log_dir
            fileops.mkdir(self.log_dir,recursive=True)
        # Other parameters
        self.sample_pattern = sample_pattern
        if self.sample_pattern is None:
            self.sample_pattern = '*'
        self.ungzip_fastqs = ungzip_fastqs
        # Verification
        self.verification_status = None
        self.fastqs_missing_qc = None
        # Reporting
        self.reporting_status = None

    @property
    def name(self):
        """
        Return project name
        """
        return self.project.name

    def check_qc(self,sched,name,wait_for=None,runner=None):
        """
        Check for Fastqs with missing/failed QC outputs

        Schedules a job to run the `reportqc` utility,
        and invokes the `_extract_fastqs` method to
        parse the output and get a list of Fastqs
        which have missing or failed QC outputs.

        The Fastq paths are stored in the
        `fastqs_missing_qc` instance variable.

        The status of the verification can be checked
        via the `verification_status` instance variable:
        None indicates verification isn't completed,
        otherwise zero indicates that verification passed,
        non-zero that it failed.

        Arguments:
          sched (SimpleScheduler): scheduler instance
            to use to run the jobs
          name (str): basename for the job
          wait_for (list): list of jobs or job groups
            to wait for before executing the check
          runner (JobRunner): job runner to use for
            QC verification
        """
        project = self.project
        print "=== Checking QC for '%s' ===" % project.name
        name = "%s.%s" % (name,project.name)
        self.verification_status = None
        self.fastqs_missing_qc = None
        collect_cmd = Command(
            "reportqc.py",
            "--fastq_dir",self.fastq_dir,
            "--qc_dir",self.qc_dir,
            "--verify",
            "--list-unverified",
            project.dirn)
        job = sched.submit(collect_cmd,
                           name="%s" % name,
                           wd=project.dirn,
                           log_dir=self.log_dir,
                           callbacks=(self._extract_fastqs,),
                           wait_for=wait_for,
                           runner=runner)

    def _extract_fastqs(self,name,jobs,sched):
        """
        Internal: callback to get Fastqs from `reportqc`

        Invoked to handle the completion of a job running
        `reportqc ... --verify --list-unverified` and parse
        the output to get a list of Fastqs with missing or
        failed QC.

        The list of Fastqs can be accessed via the
        `fastqs_missing_qc` instance attribute.
        """
        print "Extracting Fastqs with missing/failed QC"
        check_qc = jobs[0]
        print "Exit code: %s" % check_qc.exit_code
        print "Log file: %s" % check_qc.log
        print "Log file contents:"
        with open(check_qc.log,'r') as log:
            print log.read()
        self.verification_status = check_qc.exit_code
        self.fastqs_missing_qc = list()
        with open(check_qc.log,'r') as log:
            for line in log:
                if line.startswith(self.project.dirn):
                    self.fastqs_missing_qc.append(line.rstrip())
                    self.verification_status += 1
        print "Fastqs with missing QC outputs:"
        for fq in self.fastqs_missing_qc:
            print fq

    def setup_qc(self,sched,nthreads,fastq_screen_subset,
                 qc_runner=None,verify_runner=None):
        """
        Set up the QC for the project

        Arguments:
          sched (SimpleScheduler): scheduler instance
            to use to run the jobs
          nthreads (int): the maxiumum number of
            threads/cores to use per job (default: 1)
          fastq_screen_subset (int): the size of
            subset to use with FastQScreen (default:
            10000)
          qc_runner (JobRunner): job runner to use for
            executing QC
          verify_runner (JobRunner): job runner to use
            for QC verification
        """
        project = self.project
        print "=== Setting up QC for '%s' ===" % project.name
        print "Using Fastqs from %s" % project.fastq_dir
        print "Using QC directory %s" % project.qc_dir
        # Loop over samples and queue up those where the QC
        # is missing
        samples = project.get_samples(self.sample_pattern)
        if len(samples) == 0:
            logger.warning("No samples found for QC analysis in "
                           "project '%s'" % project.name)
            return
        print "%d samples matched:" % len(samples)
        for sample in samples:
            print "-- %s" % sample.name
        groups = []
        for sample in samples:
            group = None
            print "Examining files in sample %s" % sample.name
            for fq in sample.fastq:
                print "%s" % fq
                if utils.AnalysisFastq(fq).is_index_read:
                    # Reject index read Fastqs
                    logger.warning("Ignoring index read: %s" %
                                   os.path.basename(fq))
                    continue
                # Check if Fastq is in list of those with
                # missing QC outputs
                if fq not in self.fastqs_missing_qc:
                    logger.debug("\t%s: QC verified" % fq)
                else:
                    print "\t%s: setting up QC run" % os.path.basename(fq)
                    # Create a group if none exists for this sample
                    if group is None:
                        group = sched.group("%s.%s" % (project.name,
                                                       sample.name),
                                            log_dir=self.log_dir)
                    # Create and submit a QC job
                    fastq = os.path.join(project.dirn,'fastqs',fq)
                    label = "illumina_qc.%s.%s" % \
                            (project.name,str(utils.AnalysisFastq(fq)))
                    qc_cmd = Command('illumina_qc.sh',fastq)
                    if self.ungzip_fastqs:
                        qc_cmd.add_args('--ungzip-fastqs')
                    if fastq_screen_subset is None:
                        fastq_screen_subset = 0
                    qc_cmd.add_args(
                        '--threads',nthreads,
                        '--subset',fastq_screen_subset,
                        '--qc_dir',project.qc_dir)
                    job = group.add(qc_cmd,
                                    name=label,
                                    wd=project.dirn,
                                    runner=qc_runner)
                    print "Job: %s" %  job
            # Indicate no more jobs to add
            if group:
                group.close()
                groups.append(group.name)
        # Add verification job
        verify_job = self.check_qc(sched,
                                   "verify_qc",
                                   wait_for=groups,
                                   runner=verify_runner)

    def report_qc(self,sched,report_html=None,zip_outputs=True,
                  multiqc=False,runner=None):
        """
        Generate QC report

        Arguments:
          sched (SimpleScheduler):
          report_html (str): optional, path to the name of
            the QC report
          zip_outputs (bool): if True then also generate ZIP
            archive with the report and QC outputs
          multiqc (bool): if True then also generate MultiQC
            report
          runner (JobRunner): job runner to use QC reporting
        """
        print "=== Reporting QC for '%s' ===" % self.name
        project = self.project
        qc_base = os.path.basename(project.qc_dir)
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
        if self.fastq_dir is not None:
            title = "%s (%s)" % (title,self.fastq_dir)
        title = "%s: QC report" % title
        report_cmd = Command(
            "reportqc.py",
            "--fastq_dir",self.fastq_dir,
            "--qc_dir",self.qc_dir,
            "--filename",out_file,
            "--title",title)
        if zip_outputs:
            report_cmd.add_args("--zip")
        if multiqc:
            report_cmd.add_args("--multiqc")
        report_cmd.add_args(project.dirn)
        job = sched.submit(report_cmd,
                           name="report_qc.%s" % self.name,
                           wd=project.dirn,
                           log_dir=self.log_dir,
                           callbacks=(self._check_report,),
                           runner=runner)

    def _check_report(self,name,jobs,sched):
        """
        Internal: callback to check reporting status

        The status of the reporting can be checked
        via the `reporting_status` instance variable:
        None indicates reporting isn't completed,
        otherwise zero indicates that reporting finished
        okay, and non-zero that it failed.
        """
        print "Checking exit status from reporting"
        report_qc = jobs[0]
        print "Exit code: %s" % report_qc.exit_code
        print "Log file: %s" % report_qc.log
        print "Log file contents:"
        with open(report_qc.log,'r') as log:
            print log.read()
        self.reporting_status = report_qc.exit_code

    def verify(self):
        """
        Verify if the QC completed successfully

        Returns:
          Boolean: True if all QC outputs are verified, False
            if there were problems.
        """
        return (self.verification_status == 0)
