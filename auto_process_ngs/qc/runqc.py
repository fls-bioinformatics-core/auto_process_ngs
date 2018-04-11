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
import bcftbx.utils as bcf_utils
import auto_process_ngs.applications as applications
import auto_process_ngs.utils as utils
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
    def __init__(self,runner=None,max_jobs=None):
        """
        Create a new RunQC instance

        Arguments:
          runner (JobRunner): optional, job runner to
            use
          max_jobs (int): optional, specify maximum
            number of jobs to run concurrently
        """
        self._settings = Settings()
        if max_jobs is None:
            max_jobs = self._settings.general.max_concurrent_jobs
        if runner is None:
            runner = self._settings.runners.qc
        self._projects = []
        self._sched = SimpleScheduler(runner=runner,
                                      max_concurrent=max_jobs)

    def add_project(self,project,fastq_dir=None,qc_dir=None,
                    sample_pattern=None,ungzip_fastqs=False,
                    run_multiqc=True):
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
          run_multiqc (bool): if True then run MultiQC
            at the end of the QC pipeline
        """
        self._projects.append(ProjectQC(project,
                                        fastq_dir=fastq_dir,
                                        sample_pattern=sample_pattern,
                                        qc_dir=qc_dir,
                                        ungzip_fastqs=ungzip_fastqs,
                                        run_multiqc=run_multiqc))

    def run(self,nthreads=1,fastq_screen_subset=10000,
            report_html=None):
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

        Returns:
          Integer: returns 0 if QC ran to completion
            without problems, non-zero if there was
            an error.
        """
        # Start scheduler
        self._sched.start()
        # Add each project
        for project in self._projects:
            project.setup_qc(self._sched,
                             nthreads,
                             fastq_screen_subset)
        # Wait for the scheduler to run all jobs
        self._sched.wait()
        self._sched.stop()
        # Verify the outputs and generate QC reports
        failed_projects = []
        for project in self._projects:
            if not project.verify():
                failed_projects.append(project)
            else:
                project.report(report_html=report_html)
        # Report failed projects
        if failed_projects:
            logger.error("QC failed for one or more samples in the "
                         "following projects:")
            for project in failed_projects:
                logger.error("- %s" % project.project.name)
            return 1
        # Finish
        return 0

class ProjectQC(object):
    """
    Class for setting up QC jobs for a project
    """
    def __init__(self,project,fastq_dir=None,sample_pattern=None,
                 qc_dir=None,ungzip_fastqs=False,
                 run_multiqc=False):
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
          run_multiqc (bool): if True then run MultiQC
            at the end of the QC pipeline
        """
        self.project = project
        self.fastq_dir = fastq_dir
        self.sample_pattern = sample_pattern
        if self.sample_pattern is None:
            self.sample_pattern = '*'
        self.qc_dir = qc_dir
        self.ungzip_fastqs = ungzip_fastqs
        self.run_multiqc = run_multiqc
        self.multiqc_out = None

    def setup_qc(self,sched,nthreads,fastq_screen_subset):
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
        """
        project = self.project
        print "=== Setting up QC for %s ===" % project.name
        print "Subdirectories with Fastqs:"
        for fastq_dir in project.fastq_dirs:
            print "- %s" % fastq_dir
        # Set up qc directory
        project_qc_dir = project.setup_qc_dir(self.qc_dir,
                                              fastq_dir=self.fastq_dir)
        project.use_qc_dir(project_qc_dir)
        print "Using QC directory %s" % project.qc_dir
        # Sort out fastq source
        fastq_dir = project.qc_info(project_qc_dir).fastq_dir
        project.use_fastq_dir(self.fastq_dir)
        print "Gathering Fastqs from %s" % project.fastq_dir
        # Set up the logs directory
        log_dir = os.path.join(project_qc_dir,'logs')
        if not os.path.exists(log_dir):
            print "Making QC logs directory: %s" % log_dir
            bcf_utils.mkdir(log_dir,mode=0775)
        # Loop over samples and queue up those where the QC
        # isn't validated
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
                if utils.AnalysisFastq(fq).is_index_read:
                    # Reject index read Fastqs
                    logger.warning("Ignoring index read: %s" %
                                   os.path.basename(fq))
                    continue
                # TODO: - switch to a different way of checking
                # TODO:   per-sample QC?
                if sample.verify_qc(project_qc_dir,fq):
                    logger.debug("\t%s: QC verified" % fq)
                else:
                    print "\t%s: setting up QC run" % os.path.basename(fq)
                    # Create a group if none exists for this sample
                    if group is None:
                        group = sched.group("%s.%s" % (project.name,
                                                       sample.name),
                                            log_dir=log_dir)
                    # Create and submit a QC job
                    fastq = os.path.join(project.dirn,'fastqs',fq)
                    label = "illumina_qc.%s.%s" % \
                            (project.name,str(utils.AnalysisFastq(fq)))
                    qc_cmd = applications.Command('illumina_qc.sh',fastq)
                    if self.ungzip_fastqs:
                        qc_cmd.add_args('--ungzip-fastqs')
                    if fastq_screen_subset is None:
                        fastq_screen_subset = 0
                    qc_cmd.add_args(
                        '--threads',nthreads,
                        '--subset',fastq_screen_subset,
                        '--qc_dir',project_qc_dir)
                    job = group.add(qc_cmd,name=label,wd=project.dirn)
                    print "Job: %s" %  job
            # Indicate no more jobs to add
            if group:
                group.close()
                groups.append(group.name)
        # TODO: - add job to verify/do reporting (see ICELL8
        # TODO:   pipeline?)
        # TODO: - set a flag to capture the verification status
        # TODO:   (e.g. via callback?)
        # Add MultiQC job (if requested)
        if self.run_multiqc:
            self.multiqc_out = "multi%s_report.html" % \
                               os.path.basename(project_qc_dir)
            # TODO: - check that multiqc report path is correct
            # TODO:   here (I don't think it is)
            if (not os.path.exists(self.multiqc_out)) or groups:
                multiqc_cmd = applications.Command(
                    'multiqc',
                    '--title','%s/%s' % (project.info.run,
                                         project.name),
                    '--filename','./%s' % self.multiqc_out,
                    '--force',
                    project_qc_dir)
                print "Running %s" % multiqc_cmd
                label = "multiqc.%s" % project.name
                job = sched.submit(multiqc_cmd,
                                   name=label,
                                   wd=project.dirn,
                                   log_dir=log_dir,
                                   wait_for=groups)
            else:
                print "MultiQC report '%s': already exists" % \
                    self.multiqc_out

    def verify(self):
        """
        Verify if the QC completed successfully

        Returns:
          Boolean: True if all QC outputs are verified, False
            if there were problems.
        """
        verified = self.project.verify_qc(qc_dir=self.qc_dir)
        if self.run_multiqc:
            multiqc_report = os.path.join(self.project.dirn,
                                          self.multiqc_out)
            if not os.path.exists(multiqc_report):
                logger.warning("Missing MultiQC report for %s" %
                               self.project.name)
                verified = False
        return verified

    def report(self,report_html=None):
        """
        Generate QC report

        Arguments:
          report_html (str): optional, path to the name of
            the QC report.
        """
        qc_base = os.path.basename(self.project.qc_dir)
        if report_html is None:
            out_file = '%s_report.html' % qc_base
        else:
            out_file = report_html
        if not os.path.isabs(out_file):
            out_file = os.path.join(self.project.dirn,out_file)
        if self.project.info.run is None:
            title = "%s" % self.project.name
        else:
            title = "%s/%s" % (self.project.info.run,
                               self.project.name)
        if self.fastq_dir is not None:
            title = "%s (%s)" % (title,self.fastq_dir)
        title = "%s: QC report" % title
        print "QC okay, generating report for %s" % self.project.name
        self.project.qc_report(qc_dir=self.qc_dir,
                               title=title,
                               report_html=out_file)
