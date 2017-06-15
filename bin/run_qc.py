#!/usr/bin/env python
#
#     run_qc.py: run QC pipeline on arbitrary fastq files
#     Copyright (C) University of Manchester 2017 Peter Briggs
#
#########################################################################
#
# run_qc.py
#
#########################################################################

"""
Runs the QC pipeline on an arbitrary set of Fastq files
"""

#######################################################################
# Imports
#######################################################################

import sys
import os
import argparse
from bcftbx.utils import mkdir
from bcftbx.JobRunner import fetch_runner
from auto_process_ngs.utils import AnalysisProject
from auto_process_ngs.applications import Command
from auto_process_ngs.simple_scheduler import SimpleScheduler
from auto_process_ngs.qc.illumina_qc import check_qc_outputs
from auto_process_ngs import envmod

import logging
logging.basicConfig(format='%(levelname) 8s: %(message)s')

import auto_process_ngs.settings
__settings = auto_process_ngs.settings.Settings()

try:
    __modulefiles = __settings.modulefiles['run_qc']
except KeyError:
    # No environment modules specified
    __modulefiles = None

#######################################################################
# Classes
#######################################################################

#######################################################################
# Functions
#######################################################################

def announce(title):
    """
    """
    title = str(title)
    len_title = len(title)
    print "="*len_title
    print title
    print "="*len_title

#######################################################################
# Main program
#######################################################################

if __name__ == "__main__":
    # Make a command line parser
    p = argparse.ArgumentParser()
    p.add_argument("project_dir",metavar="DIR",
                   help="directory with Fastq files to run the "
                   "QC on")
    p.add_argument('--samples',metavar='PATTERN',
                   action='store',dest='sample_pattern',default=None,
                   help="simple wildcard-based pattern specifying a "
                   "subset of samples to run the QC on. If specified "
                   "then only FASTQs with sample names matching "
                   "PATTERN will be examined.")
    p.add_argument('-t','--threads',
                   action='store',dest="nthreads",default=1,
                   help="number of threads to use for QC script "
                   "(default: 1)")
    p.add_argument('-r','--runner',metavar='RUNNER',action='store',
                   dest="runner",default=__settings.runners.qc,
                   help="explicitly specify runner definition for "
                   "running QC script. RUNNER must be a valid job "
                   "runner specification e.g. 'GEJobRunner(-j y)' "
                   "(default: '%s')" % __settings.runners.qc)
    p.add_argument('--modulefiles',action='store',
                   dest='modulefiles',default=None,
                   help="comma-separated list of environment "
                   "modules to load before executing commands "
                   "(overrides any modules specified in the global "
                   "settings)")

    # Parse the command line
    args = p.parse_args()

    # Set up environment
    if args.modulefiles is None:
        modulefiles = __modulefiles
    else:
        modulefiles = args.modulefiles
    if modulefiles is not None:
        announce("Setting up environment")
        for modulefile in modulefiles.split(','):
            envmod.load(modulefile)

    # Job runner
    qc_runner = fetch_runner(args.runner)

    # Load the project
    announce("Loading project data")
    project_dir = os.path.abspath(args.project_dir)
    project_name = os.path.basename(project_dir)
    project = AnalysisProject(project_name,project_dir)

    # Get list of samples
    announce("Acquiring samples")
    project = AnalysisProject(project_name,project_dir)
    if args.sample_pattern is not None:
        samples = project.get_samples(args.sample_pattern)
    else:
        samples = project.samples
    if not samples:
        logging.warning("No samples specified for QC, quitting")
        sys.exit()
    print "%d samples matched" % len(samples)
    for sample in samples:
        print "-- %s" % sample.name

    # Sort out QC dir
    qc_dir = os.path.join(project.dirn,'qc')
    log_dir = os.path.join(qc_dir,'logs')
    mkdir(qc_dir)
    mkdir(log_dir)

    # Run the QC
    announce("Running QC")
    max_jobs = __settings.general.max_concurrent_jobs
    sched = SimpleScheduler(runner=qc_runner,
                            max_concurrent=max_jobs)
    sched.start()
    for sample in samples:
        print "Checking/setting up for sample '%s'" % sample.name
        for fq in sample.fastq:
            if sample.verify_qc(qc_dir,fq):
                print "-- %s: QC ok" % fq
            else:
                print "-- %s: setting up QC" % fq
                qc_cmd = Command('illumina_qc.sh',fq)
                if args.nthreads > 1:
                    qc_cmd.add_args('--threads',args.nthreads)
                job = sched.submit(qc_cmd,
                                   wd=project.dirn,
                                   name="qc.%s" % os.path.basename(fq),
                                   log_dir=log_dir)
                print "Job: %s" % job
    # Wait for the scheduler to run all jobs
    sched.wait()
    sched.stop()

    # Verify the QC
    announce("Verifying QC")
    qc_ok = True
    for sample in samples:
        for fq in sample.fastq:
            if not sample.verify_qc(qc_dir,fq):
                qc_ok = False
                logging.warning("-- %s: QC failed" %
                                os.path.basename(fq))
                present,missing = check_qc_outputs(fq,qc_dir)
                for output in missing:
                    logging.warning("   %s: missing" % output)
    if not qc_ok:
        logging.error("QC failed (see warnings above)")
    else:
        print "QC ok: generating report"
        project.qc_report(force=True)

