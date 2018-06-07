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
from auto_process_ngs.analysis import AnalysisProject
from auto_process_ngs.qc.illumina_qc import IlluminaQC
from auto_process_ngs.qc.runqc import RunQC
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
    ... ======

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
    p.add_argument('--fastq_screen_subset',metavar='SUBSET',
                   action='store',dest='fastq_screen_subset',
                   default=__settings.qc.fastq_screen_subset,type=int,
                   help="specify size of subset of total reads to use "
                   "for fastq_screen (i.e. --subset option); (default "
                   "%d, set to 0 to use all reads)" %
                   __settings.qc.fastq_screen_subset)
    p.add_argument('--multiqc',action='store_true',dest='run_multiqc',
                   default=False,
                   help="also generate MultiQC report")
    p.add_argument('-t','--threads',action='store',dest="nthreads",
                   type=int,default=__settings.qc.nprocessors,
                   help="number of threads to use for QC script "
                   "(default: %d)" % __settings.qc.nprocessors)
    p.add_argument('-r','--runner',metavar='RUNNER',action='store',
                   dest="runner",default=str(__settings.runners.qc),
                   help="explicitly specify runner definition for "
                   "running QC script. RUNNER must be a valid job "
                   "runner specification e.g. 'GEJobRunner(-j y)' "
                   "(default: '%s')" % __settings.runners.qc)
    p.add_argument('-m','--max-jobs',metavar='N',action='store',
                   dest='max_jobs',type=int,
                   default=__settings.general.max_concurrent_jobs,
                   help="explicitly specify maximum number of "
                   "concurrent QC jobs to run (default %d, change "
                   "in settings file)"
                   % __settings.general.max_concurrent_jobs)
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
    p.add_argument('-b','--batch',metavar='N',action='store',
                   dest='batch_size',type=int, default=None,
                   help="batch QC commands with N commands per job "
                   "(default: no batching)")
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

    # Job runners
    default_runner = __settings.general.default_runner
    qc_runner = fetch_runner(args.runner)

    # Load the project
    announce("Loading project data")
    project_dir = os.path.abspath(args.project_dir)
    project_name = os.path.basename(project_dir)
    project = AnalysisProject(project_name,project_dir)

    # Output file name
    if args.filename is None:
        out_file = None
    else:
        out_file = args.filename
        if not os.path.isabs(out_file):
            out_file = os.path.join(project.dirn,out_file)

    # Handle strand stats
    fastq_strand_conf = None
    if project.info.organism:
        print "Organisms: %s" % project.info.organism
        fastq_strand_indexes = build_fastq_strand_conf(
            project.info.organism.lower().split(','),
            __settings.fastq_strand_indexes)
        if fastq_strand_indexes:
            print "Setting up conf file for strandedness determination"
            fastq_strand_conf = os.path.join(project.dirn,
                                             "fastq_strand.conf")
            with open(fastq_strand_conf,'w') as fp:
                fp.write("%s\n" % fastq_strand_indexes)
        else:
            print "No matching indexes for strandedness determination"
    else:
        print "No organisms specified"

    # Set up QC script
    illumina_qc = IlluminaQC(
        nthreads=args.nthreads,
        fastq_screen_subset=args.fastq_screen_subset,
        fastq_strand_conf=fastq_strand_conf,
        ungzip_fastqs=False)

    # Run the QC
    announce("Running QC")
    runqc = RunQC()
    runqc.add_project(project,
                      fastq_dir=args.fastq_dir,
                      sample_pattern=args.sample_pattern,
                      qc_dir=args.qc_dir,
                      illumina_qc=illumina_qc)
    status = runqc.run(report_html=out_file,
                       multiqc=args.run_multiqc,
                       qc_runner=qc_runner,
                       verify_runner=default_runner,
                       report_runner=default_runner,
                       max_jobs=args.max_jobs,
                       batch_size=args.batch_size)
    if status:
        logger.critical("QC failed (see warnings above)")
    sys.exit(status)

