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
import logging
from bcftbx.utils import mkdir
from bcftbx.JobRunner import fetch_runner
from auto_process_ngs.utils import AnalysisProject
from auto_process_ngs.applications import Command
from auto_process_ngs.simple_scheduler import SimpleScheduler
from auto_process_ngs.qc.illumina_qc import check_qc_outputs
import auto_process_ngs
import auto_process_ngs.settings
import auto_process_ngs.envmod as envmod

# Module-specific logger
logger = logging.getLogger(__name__)

# Versions and settings
__version__ = auto_process_ngs.get_version()
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
    Print arbitrary string as a title

    Prints the supplied string as a title, e.g.

    >>> announce("Hello!")
    ... ======
    ... Hello!
    .... ======

    Arguments:
      title (str): string to print

    Returns:
      None
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
    p = argparse.ArgumentParser(
        description="Run the QC pipeline standalone on an arbitrary "
        "set of Fastq files.")
    p.add_argument('--version', action='version',
                   version=("%%(prog)s %s" % __version__))
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
                   dest="runner",default=str(__settings.runners.qc),
                   help="explicitly specify runner definition for "
                   "running QC script. RUNNER must be a valid job "
                   "runner specification e.g. 'GEJobRunner(-j y)' "
                   "(default: '%s')" % __settings.runners.qc)
    p.add_argument('--qc_dir',metavar='QC_DIR',
                   action='store',dest='qc_dir',default=None,
                   help="explicitly specify QC output directory. "
                   "NB if a relative path is supplied then it's assumed "
                   "to be a subdirectory of DIR (default: <DIR>/qc)")
    p.add_argument('-f','--filename',metavar='NAME',action='store',
                   dest='filename',default=None,
                   help="file name for output QC report (default: "
                   "<DIR>/<QC_DIR>_report.html)")
    p.add_argument('--fastq_dir',metavar='SUBDIR',
                   action='store',dest='fastq_dir',default=None,
                   help="explicitly specify subdirectory of DIR with "
                   "Fastq files to run the QC on.")
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
    project = AnalysisProject(project_name,project_dir,
                              fastq_dir=args.fastq_dir)
    print "Subdirectories with Fastqs:"
    for fastq_dir in project.fastq_dirs:
        print "- %s" % fastq_dir
    print "Gathering Fastqs from %s" % project.fastq_dir
    if args.sample_pattern is not None:
        samples = project.get_samples(args.sample_pattern)
    else:
        samples = project.samples
    if not samples:
        logger.warning("No samples specified for QC, quitting")
        sys.exit()
    print "%d samples matched" % len(samples)
    for sample in samples:
        print "-- %s" % sample.name

    # Sort out QC dir
    if args.qc_dir is None:
        qc_dir = 'qc'
    else:
        qc_dir = args.qc_dir
    if not os.path.isabs(qc_dir):
        qc_dir = os.path.join(project.dirn,qc_dir)
    print "QC output dir: %s" % qc_dir
    log_dir = os.path.join(qc_dir,'logs')
    mkdir(qc_dir)
    mkdir(log_dir)
    qc_base = os.path.basename(qc_dir)

    # Output file name
    if args.filename is None:
        out_file = '%s_report.html' % qc_base
    else:
        out_file = args.filename
    if not os.path.isabs(out_file):
        out_file = os.path.join(project.dirn,out_file)
    print "QC report: %s" % out_file

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
                qc_cmd.add_args('--qc_dir',qc_dir)
                job = sched.submit(qc_cmd,
                                   wd=project.dirn,
                                   name="%s.%s" % (qc_base,
                                                   os.path.basename(fq)),
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
                logger.warning("-- %s: QC failed" %
                               os.path.basename(fq))
                present,missing = check_qc_outputs(fq,qc_dir)
                for output in missing:
                    logger.warning("   %s: missing" % output)
    if not qc_ok:
        logger.error("QC failed (see warnings above)")
    else:
        print "QC ok: generating report"
        project.qc_report(report_html=out_file,
                          qc_dir=qc_dir,
                          force=True)

