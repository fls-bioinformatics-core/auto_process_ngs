#!/usr/bin/env python
#
#     run_qc.py: run QC pipeline on arbitrary fastq files
#     Copyright (C) University of Manchester 2017-2020 Peter Briggs
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
import psutil
import math
import logging
from bcftbx.JobRunner import fetch_runner
from bcftbx.JobRunner import SimpleJobRunner
from auto_process_ngs.analysis import AnalysisProject
import auto_process_ngs
import auto_process_ngs.settings
import auto_process_ngs.envmod as envmod
from auto_process_ngs.qc.pipeline import QCPipeline
from auto_process_ngs.tenx_genomics_utils import CELLRANGER_ASSAY_CONFIGS

# QC protocols
from auto_process_ngs.qc.constants import PROTOCOLS

# Module-specific logger
logger = logging.getLogger(__name__)

# Versions and settings
__version__ = auto_process_ngs.get_version()
__settings = auto_process_ngs.settings.Settings()
try:
    __modulefiles = __settings.modulefiles
except KeyError:
    # No environment modules specified
    __modulefiles = None

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
    print("="*len_title)
    print(title)
    print("="*len_title)

#######################################################################
# Main program
#######################################################################

if __name__ == "__main__":
    # Make a command line parser
    p = argparse.ArgumentParser(
        description="Run the QC pipeline standalone on an arbitrary "
        "set of Fastq files.")
    # Defaults
    default_nthreads = __settings.qc.nprocessors
    # Build parser
    p.add_argument('--version', action='version',
                   version=("%%(prog)s %s" % __version__))
    p.add_argument("project_dir",metavar="DIR",
                   help="directory with Fastq files to run the "
                   "QC on")
    p.add_argument('-p','--protocol',metavar='PROTOCOL',
                   action='store',dest='qc_protocol',default=None,
                   choices=PROTOCOLS,
                   help="explicitly specify the QC protocol to use; "
                   "can be one of %s. If not set then protocol will "
                   "be determined automatically based on directory "
                   "contents." %
                   ", ".join(["'%s'" % x for x in PROTOCOLS]))
    p.add_argument('--organism',metavar='ORGANISM',
                   action='store',dest='organism',default=None,
                   help="explicitly specify organism (e.g. 'human', "
                   "'mouse'). Multiple organisms should be separated "
                   "by commas (e.g. 'human,mouse')")
    p.add_argument('--fastq_screen_subset',metavar='SUBSET',
                   action='store',dest='fastq_screen_subset',
                   default=__settings.qc.fastq_screen_subset,type=int,
                   help="specify size of subset of total reads to use "
                   "for fastq_screen (i.e. --subset option); (default "
                   "%d, set to 0 to use all reads)" %
                   __settings.qc.fastq_screen_subset)
    p.add_argument('-t','--threads',action='store',dest="nthreads",
                   type=int,default=None,
                   help="number of threads to use for QC script "
                   "(default: %s)" % ('taken from job runner'
                                      if not default_nthreads
                                      else default_nthreads,))
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
    # Reporting options
    reporting = p.add_argument_group('Output and reporting')
    reporting.add_argument('--qc_dir',metavar='QC_DIR',
                           action='store',dest='qc_dir',default=None,
                           help="explicitly specify QC output directory. "
                           "NB if a relative path is supplied then it's "
                           "assumed to be a subdirectory of DIR (default: "
                           "<DIR>/qc)")
    reporting.add_argument('-f','--filename',metavar='NAME',action='store',
                           dest='filename',default=None,
                           help="file name for output QC report (default: "
                           "<DIR>/<QC_DIR>_report.html)")
    reporting.add_argument('--multiqc',action='store_true',
                           dest='run_multiqc', default=False,
                           help="also generate MultiQC report")
    # Cellranger options
    cellranger = p.add_argument_group('Cellranger/10xGenomics options')
    cellranger.add_argument('--10x_transcriptome',action='append',
                            metavar='ORGANISM=REFERENCE',
                            dest='cellranger_transcriptomes',
                            help="specify cellranger transcriptome "
                            "reference datasets to associate with "
                            "organisms (overrides references defined "
                            "in config file)")
    cellranger.add_argument('--10x_premrna_reference',action='append',
                            metavar='ORGANISM=REFERENCE',
                            dest='cellranger_premrna_references',
                            help="specify cellranger pre-mRNA reference "
                            "datasets to associate with organisms "
                            "(overrides references defined in config "
                            "file)")
    cellranger.add_argument("--10x_chemistry",
                            choices=sorted(CELLRANGER_ASSAY_CONFIGS.keys()),
                            dest="cellranger_chemistry",default="auto",
                            help="assay configuration for 10xGenomics "
                            "scRNA-seq; if set to 'auto' (the default) then "
                            "cellranger will attempt to determine this "
                            "automatically")
    # Pipeline/job options
    pipeline = p.add_argument_group('Job control options')
    pipeline.add_argument('--local',action='store_true',
                          dest='local',default=False,
                          help="run the QC on the local system (overrides "
                          "any runners defined in the configuration or on "
                          "the command line)")
    pipeline.add_argument('-c','--maxcores',metavar='N',action='store',
                          dest='max_cores',type=int,default=None,
                          help="maximum number of cores available for QC "
                          "jobs when using --local (default: unlimited)")
    pipeline.add_argument('-m','--maxmem',metavar='M',action='store',
                          dest='max_mem',type=int,default=None,
                          help="maximum total memory jobs can request at "
                          "once when using --local (in Gbs; default: "
                          "unlimited)")
    pipeline.add_argument('-j','--maxjobs',metavar='N',action='store',
                          dest='max_jobs',type=int,
                          default=__settings.general.max_concurrent_jobs,
                          help="explicitly specify maximum number of "
                          "concurrent QC jobs to run (default %d, change "
                          "in settings file; ignored when using --local)"
                          % __settings.general.max_concurrent_jobs)
    # Advanced options
    advanced = p.add_argument_group('Advanced/debugging options')
    advanced.add_argument('--verbose',action="store_true",
                          dest="verbose",default=False,
                          help="run pipeline in 'verbose' mode")
    advanced.add_argument('-r','--runner',metavar='RUNNER',action='store',
                          dest="runner",default=None,
                          help="explicitly specify runner definition for "
                          "running QC script. RUNNER must be a valid job "
                          "runner specification e.g. 'GEJobRunner(-j y)' "
                          "(default: '%s')" % __settings.runners.qc)
    advanced.add_argument('-b','--batch',metavar='N',action='store',
                          dest='batch_size',type=int, default=None,
                          help="batch QC commands with N commands per job "
                          "(default: no batching)")
    advanced.add_argument('--work-dir',action="store",
                          dest="working_dir",default=None,
                          help="specify the working directory for the "
                          "pipeline operations")

    # Parse the command line
    args = p.parse_args()

    # Set up environment
    envmodules = dict()
    if args.modulefiles is None:
        modulefiles = __modulefiles['run_qc']
    else:
        modulefiles = args.modulefiles
    if modulefiles is not None:
        announce("Setting up environment")
        for modulefile in modulefiles.split(','):
            envmod.load(modulefile)

    # Per task environment modules
    for name in ('illumina_qc',
                 'fastq_strand',
                 'cellranger',
                 'report_qc',):
        try:
            envmodules[name] = __modulefiles[name]
        except KeyError:
            envmodules[name] = None

    # Maximum number of jobs and cores
    max_jobs = args.max_jobs
    max_cores = args.max_cores

    # Cellranger data
    cellranger_settings = __settings['10xgenomics']
    cellranger_transcriptomes = dict()
    if __settings['10xgenomics_transcriptomes']:
        for organism in __settings['10xgenomics_transcriptomes']:
            if organism not in cellranger_transcriptomes:
                cellranger_transcriptomes[organism] = \
                    __settings['10xgenomics_transcriptomes'][organism]
    if args.cellranger_transcriptomes:
        for transcriptome in args.cellranger_transcriptomes:
            organism,reference =  transcriptome.split('=')
            cellranger_transcriptomes[organism] = reference
    cellranger_premrna_references = dict()
    if __settings['10xgenomics_premrna_references']:
        for organism in __settings['10xgenomics_premrna_references']:
            if organism not in cellranger_premrna_references:
                cellranger_premrna_references[organism] = \
                    __settings['10xgenomics_premrna_references'][organism]
    if args.cellranger_premrna_references:
        for premrna_reference in args.cellranger_premrna_references:
            organism,reference =  premrna_reference.split('=')
            cellranger_premrna_references[organism] = reference
    cellranger_atac_references = __settings['10xgenomics_atac_genome_references']

    # Job runners
    announce("Configuring pipeline parameters")
    if args.local:
        # Running in 'local' mode
        # Used local runners and set defaults according to
        # resources available on local system
        print("Running locally: overriding settings in configuration")
        # Set maximum number of slots
        if not max_cores:
            try:
                # If NSLOTS is set in the environment
                # then assume we're running on an SGE
                # node and this sets the maximum number
                # of available cores
                max_cores = int(os.environ['NSLOTS'])
            except KeyError:
                # Set limit from local machine
                max_cores = psutil.cpu_count()
        print("-- Maximum cores: %s" % max_cores)
        # Set maximum memory
        if args.max_mem:
            max_mem = args.max_mem
        else:
            # Maximum memory is scaled by the proportion
            # of the total cores being used
            # Needs to be converted from bytes to Gbs
            max_mem = math.floor(
                float(psutil.virtual_memory().total)/(1024.0**3)
                *float(max_cores)/float(psutil.cpu_count()))
        print("-- Maximum memory: %s Gbs" % max_mem)
        # Set number of threads for QC jobs
        if args.nthreads:
            nthreads = args.nthreads
        else:
            nthreads = min(max_cores,8)
        print("-- Threads for QC: %s" % nthreads)
        # Remove limit on number of jobs
        print("-- Set maximum no of jobs to 'unlimited'")
        max_jobs = None
        # (Re)set cellranger parameters for --local
        print("-- Cellranger will run in jobmode 'local'")
        cellranger_jobmode = "local"
        cellranger_mempercore = None
        cellranger_jobinterval = None
        cellranger_localcores = max_cores
        cellranger_localmem = max_mem
        # Set up local runners
        default_runner = SimpleJobRunner()
        runners = {
            'cellranger_runner': SimpleJobRunner(nslots=cellranger_localcores),
            'qc_runner': SimpleJobRunner(nslots=nthreads),
            'verify_runner': default_runner,
            'report_runner': default_runner,
        }
    else:
        # Set up according to the configuration and
        # command line options
        # Set number of threads for QC jobs
        if args.nthreads:
            nthreads = args.nthreads
        else:
            nthreads = __settings.qc.nprocessors
        # Cellranger settings
        cellranger_jobmode = cellranger_settings.cellranger_jobmode
        cellranger_mempercore = cellranger_settings.cellranger_mempercore
        cellranger_jobinterval = cellranger_settings.cellranger_jobinterval
        cellranger_localcores = cellranger_settings.cellranger_localcores
        cellranger_localmem = cellranger_settings.cellranger_localmem
        # Set up runners
        if args.runner is not None:
            # Runner explicitly supplied on the command line
            print("Setting up runners supplied on command line")
            default_runner = fetch_runner(args.runner)
            runners = {
                'cellranger_runner': default_runner,
                'qc_runner': default_runner,
                'verify_runner': default_runner,
                'report_runner': default_runner,
            }
        else:
            # Runners from configuration
            print("Setting up runners from configuration")
            default_runner = __settings.general.default_runner
            runners = {
                'cellranger_runner': __settings.runners.cellranger,
                'qc_runner': __settings.runners.qc,
                'verify_runner': default_runner,
                'report_runner': default_runner,
            }

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

    # Set up and run the QC pipeline
    announce("Running QC pipeline")
    runqc = QCPipeline()
    runqc.add_project(project,
                      qc_dir=args.qc_dir,
                      fastq_dir=args.fastq_dir,
                      organism=args.organism,
                      qc_protocol=args.qc_protocol,
                      multiqc=args.run_multiqc)
    status = runqc.run(nthreads=nthreads,
                       fastq_subset=args.fastq_screen_subset,
                       fastq_strand_indexes=
                       __settings.fastq_strand_indexes,
                       cellranger_chemistry=\
                       args.cellranger_chemistry,
                       cellranger_transcriptomes=cellranger_transcriptomes,
                       cellranger_premrna_references=\
                       cellranger_premrna_references,
                       cellranger_atac_references=cellranger_atac_references,
                       cellranger_jobmode=cellranger_jobmode,
                       cellranger_maxjobs=max_jobs,
                       cellranger_mempercore=cellranger_mempercore,
                       cellranger_jobinterval=cellranger_jobinterval,
                       cellranger_localcores=cellranger_localcores,
                       cellranger_localmem=cellranger_localmem,
                       max_jobs=max_jobs,
                       max_slots=max_cores,
                       batch_size=args.batch_size,
                       runners=runners,
                       default_runner=default_runner,
                       envmodules=envmodules,
                       working_dir=args.working_dir,
                       verbose=args.verbose)
    if status:
        logger.critical("QC failed (see warnings above)")
    sys.exit(status)
