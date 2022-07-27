#!/usr/bin/env python
#
#     run_qc.py: run QC pipeline on arbitrary fastq files
#     Copyright (C) University of Manchester 2017-2022 Peter Briggs
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
import glob
import tempfile
import shutil
import atexit
import logging
from bcftbx.JobRunner import fetch_runner
from bcftbx.JobRunner import SimpleJobRunner
from auto_process_ngs.analysis import AnalysisProject
from auto_process_ngs.analysis import locate_project_info_file
from auto_process_ngs.metadata import AnalysisProjectInfo
from auto_process_ngs.metadata import AnalysisProjectQCDirInfo
from auto_process_ngs.fastq_utils import IlluminaFastqAttrs
from auto_process_ngs.fastq_utils import group_fastqs_by_name
import auto_process_ngs
import auto_process_ngs.settings
import auto_process_ngs.envmod as envmod
from auto_process_ngs.qc.pipeline import QCPipeline
from auto_process_ngs.qc.utils import report_qc
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

def cleanup_atexit(tmp_project_dir):
    """
    Perform clean up actions on exit

    Removes the temporary project directory
    created for running the QC
    """
    if os.path.isdir(tmp_project_dir):
        print("Removing temporary project directory: %s"
              % tmp_project_dir)
        shutil.rmtree(tmp_project_dir)

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
    p.add_argument("inputs",metavar="DIR | FASTQ [ FASTQ ... ]",
                   nargs="+",
                   help="directory or list of Fastq files to run the "
                   "QC on")
    # Reporting options
    reporting = p.add_argument_group('Output and reporting')
    reporting.add_argument('-n','--name',action='store',
                           help="name for the project (used in report "
                           "title)")
    reporting.add_argument('-o','--out_dir',action='store',
                           help="top-level directory for reports and "
                           "QC output subdirectory (default: current "
                           "working directory)")
    reporting.add_argument('--qc_dir',
                           help="explicitly specify QC output directory. "
                           "NB if a relative path is supplied then it's "
                           "assumed to be a subdirectory of OUT_DIR "
                           "(default: <OUT_DIR>/qc)")
    reporting.add_argument('-f','--filename',action='store',
                           help="file name for output QC report (default: "
                           "<OUT_DIR>/<QC_DIR_NAME>_report.html)")
    reporting.add_argument('-u','--update',action='store_true',
                           help="force QC pipeline to run even if output "
                           "QC directory already exists in <OUT_DIR> "
                           "(default: stop if output QC directory already "
                           "exists)")
    reporting.add_argument('--force',action='store_true',
                           help="force generation of HTML report even if "
                           "pipeline fails (default: don't generate HTML "
                           "report on pipeline failure")
    # Project metadata
    metadata = p.add_argument_group('Metadata')
    metadata.add_argument('--organism',metavar='ORGANISM',
                          action='store',dest='organism',default=None,
                          help="explicitly specify organism (e.g. "
                          "'human', 'mouse'). Multiple organisms "
                          "should be separated by commas (e.g. "
                          "'human,mouse')")
    metadata.add_argument('--library-type',metavar='LIBRARY',
                          action='store',dest='library_type',default=None,
                          help="explicitly specify library type (e.g. "
                          "'RNA-seq', 'ChIP-seq')")
    metadata.add_argument('--single-cell-platform',metavar='PLATFORM',
                          action='store',dest='single_cell_platform',
                          default=None,
                          help="explicitly specify the single cell "
                          "platform (e.g. '10xGenomics Chromium 3'v3')")
    # QC pipeline options
    qc_options = p.add_argument_group('QC options')
    qc_options.add_argument('-p','--protocol',metavar='PROTOCOL',
                            action='store',dest='qc_protocol',
                            default=None,choices=PROTOCOLS,
                            help="explicitly specify the QC protocol to "
                            "use; can be one of %s. If not set then "
                            "protocol will be determined automatically "
                            "based on directory contents and metadata." %
                            ", ".join(["'%s'" % x for x in PROTOCOLS]))
    qc_options.add_argument('--fastq_screen_subset',metavar='SUBSET',
                            action='store',dest='fastq_screen_subset',
                            default=__settings.qc.fastq_screen_subset,
                            type=int,
                            help="specify size of subset of total reads "
                            "to use for fastq_screen (i.e. --subset "
                            "option); (default %d, set to 0 to use all "
                            "reads)" %
                            __settings.qc.fastq_screen_subset)
    qc_options.add_argument('-t','--threads',action='store',
                            dest="nthreads",type=int,default=None,
                            help="number of threads to use for QC script "
                            "(default: %s)" % ('taken from job runner'
                                               if not default_nthreads
                                               else default_nthreads,))
    # Cellranger options
    cellranger = p.add_argument_group('Cellranger/10xGenomics options')
    cellranger.add_argument('--cellranger',action='store',
                            metavar='CELLRANGER_EXE',
                            dest='cellranger_exe',
                            help="explicitly specify path to Cellranger "
                            "executable to use for single library "
                            "analysis")
    cellranger.add_argument('--cellranger-reference',action='store',
                            metavar='REFERENCE',
                            dest='cellranger_reference_dataset',
                            help="specify the path to the reference "
                            "dataset to use when running single libary "
                            "analysis (overrides the organism-specific "
                            "references defined in the config file)")
    cellranger.add_argument("--10x_chemistry",
                            choices=sorted(CELLRANGER_ASSAY_CONFIGS.keys()),
                            dest="cellranger_chemistry",default="auto",
                            help="assay configuration for 10xGenomics "
                            "scRNA-seq; if set to 'auto' (the default) then "
                            "cellranger will attempt to determine this "
                            "automatically")
    cellranger.add_argument("--10x_force_cells",action='store',
                            metavar="N_CELLS",
                            dest="cellranger_force_cells",
                            help="force number of cells for 10xGenomics "
                            "scRNA-seq and scATAC-seq, overriding automatic "
                            "cell detection algorithms (default is to use "
                            "built-in cell detection)")
    # Conda options
    conda = p.add_argument_group("Conda dependency resolution")
    conda.add_argument('--enable-conda',choices=["yes","no"],
                       dest="enable_conda",default=None,
                       help="use conda to resolve task dependencies; can "
                       "be 'yes' or 'no' (default: %s)" %
                       ("yes" if __settings.conda.enable_conda else "no"))
    conda.add_argument('--conda-env-dir',action='store',
                       dest="conda_env_dir",default=__settings.conda.env_dir,
                       help="specify directory for conda enviroments "
                       "(default: %s)" % ('temporary directory'
                                          if not __settings.conda.env_dir else
                                          __settings.conda.env_dir))
    # Pipeline/job options
    pipeline = p.add_argument_group('Job control options')
    pipeline.add_argument('--local',action='store_true',
                          dest='local',default=False,
                          help="run the QC on the local system (overrides "
                          "any runners defined in the configuration or on "
                          "the command line)")
    pipeline.add_argument('-c','--maxcores',metavar='N',action='store',
                          dest='max_cores',type=int,
                          default=__settings.general.max_cores,
                          help="maximum number of cores available for QC "
                          "jobs when using --local (default %s, change in "
                          "in settings file)" %
                          (__settings.general.max_cores
                           if __settings.general.max_cores else 'no limit'))
    pipeline.add_argument('-m','--maxmem',metavar='M',action='store',
                          dest='max_mem',type=int,default=None,
                          help="maximum total memory jobs can request at "
                          "once when using --local (in Gbs; default: "
                          "unlimited)")
    pipeline.add_argument('-j','--maxjobs',metavar='N',action='store',
                          dest='max_jobs',type=int,
                          default=__settings.general.max_concurrent_jobs,
                          help="explicitly specify maximum number of "
                          "concurrent QC jobs to run (default %s, change "
                          "in settings file; ignored when using --local)"
                          % (__settings.general.max_concurrent_jobs
                             if __settings.general.max_concurrent_jobs
                             else 'no limit'))
    pipeline.add_argument('-b','--maxbatches',type=int,action='store',
                          dest='max_batches',metavar='NBATCHES',
                          default=__settings.general.max_batches,
                          help="enable dynamic batching of pipeline "
                          "jobs with maximum number of batches set to "
                          "NBATCHES (default: %s)"
                          % (__settings.general.max_batches
                             if __settings.general.max_batches
                             else 'no batching'))
    # Advanced options
    advanced = p.add_argument_group('Advanced options')
    advanced.add_argument('-r','--runner',metavar='RUNNER',action='store',
                          dest="runner",default=None,
                          help="explicitly specify runner definition for "
                          "running QC components. RUNNER must be a valid job "
                          "runner specification e.g. 'GEJobRunner(-j y)' "
                          "(default: use runners set in configuration)")
    advanced.add_argument('-s','--batch_size',metavar='N',action='store',
                          dest='batch_size',type=int, default=None,
                          help="batch QC commands with N commands per job "
                          "(default: no batching)")
    advanced.add_argument('--ignore-metadata',action="store_true",
                          dest="ignore_metadata",default=False,
                          help="ignore information from project metadata "
                          "file even if one is located (default is to use "
                          "project metadata)")
    advanced.add_argument('--use-legacy-screen-names',choices=['yes','no'],
                          dest="use_legacy_screen_names",default=None,
                          help="use 'legacy' naming convention for "
                          "FastqScreen output files; can be 'yes' or 'no' "
                          "(default: %s)" %
                          ("yes" if __settings.qc.use_legacy_screen_names
                           else "no"))
    advanced.add_argument('--no-multiqc',action="store_true",
                          dest="no_multiqc",default=False,
                          help="turn off generation of MultiQC report")
    # Debugging options
    debugging = p.add_argument_group('Debugging options')
    debugging.add_argument('--verbose',action="store_true",
                           dest="verbose",default=False,
                           help="run pipeline in 'verbose' mode")
    debugging.add_argument('--work-dir',action="store",
                          dest="working_dir",default=None,
                          help="specify the working directory for the "
                          "pipeline operations")
    debugging.add_argument('--no-cleanup',action='store_true',
                           dest="no_cleanup",default=False,
                           help="don't remove the temporary project "
                           "directory on completion (by default the "
                           "temporary directory is deleted)")
    # Deprecated options
    deprecated = p.add_argument_group('Deprecated/redundant options')
    deprecated.add_argument('--multiqc',action='store_true',
                            dest='run_multiqc', default=False,
                            help="redundant: MultiQC report is generated "
                            "by default (use --no-multiqc to disable)")
    deprecated.add_argument('--modulefiles',action='store',
                            dest='modulefiles',default=None,
                            help="deprecated: environment modules should "
                            "be set on per-task basis")
    deprecated.add_argument('--fastq_dir',metavar='SUBDIR',
                            action='store',dest='fastq_dir',default=None,
                            help="unsupported: point DIR directly to Fastq "
                            "directory")
    deprecated.add_argument('--10x_transcriptome',action='append',
                            metavar='ORGANISM=REFERENCE',
                            dest='cellranger_transcriptomes',
                            help="deprecated: use --cellranger-reference "
                            "to specify reference dataset for single "
                            "library analysis")
    deprecated.add_argument('--10x_premrna_reference',action='append',
                            metavar='ORGANISM=REFERENCE',
                            dest='cellranger_premrna_references',
                            help="deprecated: use --cellranger-reference "
                            "to specify reference dataset for single "
                            "library analysis")

    # Parse the command line
    args = p.parse_args()

    # Check for deprecated and unsupported options
    if args.run_multiqc:
        logger.warning("'--multiqc' option is redundant; MultiQC is "
                       "run by default")
        logger.warning("Use '--no-multiqc' to turn off MultiQC report "
                       "generation")
    if args.modulefiles:
        logger.warning("'--modulefiles' option is deprecated and may not "
                       "work")
        logger.warning("Environment modules should be set on a per-task "
                       "basis in the config file")
    if args.cellranger_transcriptomes:
        logger.warning("'--10x_transcriptome' option is deprecated")
        logger.warning("Use '--cellranger-reference' instead")
    if args.cellranger_premrna_references:
        logger.warning("'--10x_premrna_reference' option is deprecated")
        logger.warning("Use '--cellranger-reference' instead")
    if args.fastq_dir:
        logger.fatal("'--fastq_dir' option is no longer supported")
        logger.fatal("Point 'DIR' directly to the directory with the "
                     "Fastq files instead")
        sys.exit(1)

    # Initialise
    project_metadata = AnalysisProjectInfo()
    dir_path = os.getcwd()
    out_dir = args.out_dir
    qc_dir = args.qc_dir
    master_fastq_dir = None
    fastq_attrs = IlluminaFastqAttrs
    extra_files = set()

    # Deal with inputs
    #
    # Possibilities are:
    # - subdirectory in a project
    # - project directory
    # - non-project directory with Fastqs
    # - list of Fastqs
    announce("Locating inputs")
    inputs = []
    for f in args.inputs:
        for ff in glob.glob(os.path.abspath(f)):
            if not os.path.exists(ff):
                # Input not found
                logger.fatal("%s: input not found" % ff)
                sys.exit(1)
            elif os.path.isdir(ff) and len(args.inputs) > 1:
                # Can only be a single directory
                logger.fatal("Input must be a single directory, or a list of "
                             "Fastqs")
                sys.exit(1)
            else:
                inputs.append(ff)
    # Get list of Fastqs from directory
    dir_path = None
    if len(inputs) == 1 and os.path.isdir(inputs[0]):
        dir_path = inputs[0]
        if not os.path.isdir(dir_path):
            logger.fatal("%s: directory not found" % dir_path)
            sys.exit(1)
        # See if directory contains Fastqs
        inputs = [os.path.join(dir_path,f)
                  for f in os.listdir(inputs[0])
                  if (f.endswith('.fastq') or
                      f.endswith('.fq') or
                      f.endswith('.fastq.gz'))]
        if not inputs:
            # No Fastqs, try loading as a project
            inputs = list(AnalysisProject(dir_path).fastqs)
            master_fastq_dir = AnalysisProject(dir_path).fastq_dir
        else:
            # Store the source directory for Fastqs
            master_fastq_dir = dir_path
        # Check we have some Fastqs
        if not inputs:
            logger.fatal("%s: no Fastqs found" % dir_path)
            sys.exit(1)
        # Look for project metadata
        info_file = locate_project_info_file(dir_path)
        # Look for extra files
        for f in ('10x_multiome_libraries.info',
                  '10x_multi_config.csv',):
            ff = os.path.join(dir_path,f)
            if os.path.isfile(ff):
                extra_files.add(ff)
    else:
        # Look for a metadata file based on Fastqs
        for fq in inputs:
            info_file = locate_project_info_file(os.path.dirname(fq))
            if info_file:
                break
    # Filter out index reads
    inputs = [fq for fq in inputs
              if not fastq_attrs(fq).is_index_read]
    if not inputs:
        logger.fatal("No Fastqs found")
        sys.exit(1)

    # Report what was found
    for fqs in group_fastqs_by_name(inputs,fastq_attrs=fastq_attrs):
        print("%s:" % fastq_attrs(fqs[0]).sample_name)
        for fq in fqs:
            print("  %s" % fq)
    print("Located %s Fastq%s" % (len(inputs),
                                  's' if len(inputs) != 1 else ''))
    if info_file:
        print("Located project metadata in %s%s" % (info_file,
                                                    " (will be ignored)"
                                                    if args.ignore_metadata
                                                    else ''))
    else:
        print("Unable to locate project metadata")

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
    for name in ('fastqc',
                 'fastq_screen',
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

    # Conda dependency resolution
    enable_conda = args.enable_conda
    if enable_conda is None:
        enable_conda = __settings.conda.enable_conda
    else:
        enable_conda = (args.enable_conda == "yes")

    # Fastq screens
    if __settings.qc.fastq_screens:
        fastq_screens = dict()
        for screen in __settings.qc.fastq_screens.split(','):
            fastq_screens[screen] = __settings.screens[screen].conf_file
    else:
        fastq_screens = None
    use_legacy_screen_names = args.use_legacy_screen_names
    if use_legacy_screen_names is None:
        use_legacy_screen_names = __settings.qc.use_legacy_screen_names
    else:
        use_legacy_screen_names = (use_legacy_screen_names == "yes")

    # STAR indexes
    star_indexes = dict()
    for organism in __settings.organisms:
        star_index = __settings.organisms[organism].star_index
        if star_index:
            star_indexes[organism] = star_index

    # Cellranger settings
    cellranger_settings = __settings['10xgenomics']

    # Cellranger reference datasets
    cellranger_transcriptomes = dict()
    for organism in __settings.organisms:
        if __settings.organisms[organism].cellranger_reference:
            cellranger_transcriptomes[organism] = \
                __settings.organisms[organism].cellranger_reference
    cellranger_premrna_references = dict()
    for organism in __settings.organisms:
        if __settings.organisms[organism].cellranger_premrna_reference:
            cellranger_premrna_references[organism] = \
                __settings.organisms[organism].cellranger_premrna_reference
    cellranger_atac_references = dict()
    for organism in __settings.organisms:
        if __settings.organisms[organism].cellranger_atac_reference:
            cellranger_atac_references[organism] = \
                __settings.organisms[organism].cellranger_atac_reference
    cellranger_multiome_references = dict()
    for organism in __settings.organisms:
        if __settings.organisms[organism].cellranger_arc_reference:
            cellranger_multiome_references[organism] = \
                __settings.organisms[organism].cellranger_arc_reference

    # Add in organism-specific references supplied on the command line
    if args.cellranger_transcriptomes:
        for transcriptome in args.cellranger_transcriptomes:
            organism,reference =  transcriptome.split('=')
            cellranger_transcriptomes[organism] = reference
    if args.cellranger_premrna_references:
        for premrna_reference in args.cellranger_premrna_references:
            organism,reference =  premrna_reference.split('=')
            cellranger_premrna_references[organism] = reference

    # Single reference supplied on command line
    cellranger_reference_dataset = args.cellranger_reference_dataset
    if cellranger_reference_dataset:
        cellranger_reference_dataset = os.path.abspath(
            cellranger_reference_dataset)

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
        # Set number of threads for STAR jobs
        if args.nthreads:
            nthreads_star = args.nthreads
        else:
            mempercore = max_mem/float(max_cores)
            nthreads_star = int(math.ceil(32.0/mempercore))
        print("-- Threads for STAR: %s" % nthreads_star)
        # Remove limit on number of jobs
        print("-- Set maximum no of jobs to 'unlimited'")
        max_jobs = None
        # (Re)set cellranger parameters for --local
        print("-- Cellranger will run in jobmode 'local'")
        cellranger_jobmode = "local"
        cellranger_mempercore = None
        cellranger_jobinterval = None
        cellranger_localcores = min(max_cores,16)
        cellranger_localmem = max_mem
        print("-- Cellranger localcores: %s" % cellranger_localcores)
        print("-- Cellranger localmem  : %s" % cellranger_localmem)
        # Set up local runners
        default_runner = SimpleJobRunner()
        runners = {
            'cellranger_runner': SimpleJobRunner(nslots=cellranger_localcores),
            'fastqc_runner': SimpleJobRunner(nslots=nthreads),
            'fastq_screen_runner': SimpleJobRunner(nslots=nthreads),
            'star_runner': SimpleJobRunner(nslots=nthreads_star),
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
                'fastqc_runner': default_runner,
                'fastq_screen_runner': default_runner,
                'star_runner': default_runner,
                'verify_runner': default_runner,
                'report_runner': default_runner,
            }
        else:
            # Runners from configuration
            print("Setting up runners from configuration")
            default_runner = __settings.general.default_runner
            runners = {
                'cellranger_runner': __settings.runners.cellranger,
                'fastqc_runner': __settings.runners.fastqc,
                'fastq_screen_runner': __settings.runners.fastq_screen,
                'star_runner': __settings.runners.star,
                'verify_runner': default_runner,
                'report_runner': default_runner,
            }

    # Output directory
    announce("Setting up output destinations")
    if not out_dir:
        if dir_path:
            out_dir = dir_path
        else:
            out_dir = os.getcwd()
    out_dir = os.path.abspath(out_dir)
    print("Output directory: %s" % out_dir)
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    # QC directory
    if not qc_dir:
        qc_dir = 'qc'
    qc_dir = os.path.join(out_dir,qc_dir)
    print("QC directory    : %s" % qc_dir)

    # Check if QC directory already exists
    if os.path.exists(qc_dir):
        if not args.update:
            logger.fatal("QC directory already exists (use --update to "
                         "run QC anyway)")
            sys.exit(1)
        print("Output QC directory already exists, updating")

    # Output file name
    if args.filename is None:
        out_file = "%s_report.html" % os.path.basename(qc_dir)
    else:
        out_file = args.filename
    if not os.path.isabs(out_file):
        out_file = os.path.join(out_dir,out_file)
    print("Output report: %s" % out_file)

    # Build and populate a temporary project directory
    announce("Building temporary project directory")
    project_dir = tempfile.mkdtemp(suffix=".run_qc",dir=os.getcwd())
    print("Building temporary project directory '%s'" % project_dir)
    fastq_dir = os.path.join(project_dir,"fastqs")
    os.mkdir(fastq_dir)
    print("Populating %s" % fastq_dir)
    for fq in inputs:
        # Make symlinks to the Fastq files
        os.symlink(fq,os.path.join(fastq_dir,os.path.basename(fq)))

    # Set up metadata
    if info_file and not args.ignore_metadata:
        project_metadata.load(info_file)
    if args.name:
        # Set project name to user-supplied value
        project_metadata['name'] = args.name
    elif project_metadata.name is None:
        # Set to output directory name if not already set
        project_metadata['name'] = os.path.basename(out_dir)
    if args.organism:
        project_metadata['organism'] = args.organism
    if args.library_type:
        project_metadata['library_type'] = args.library_type
    if args.single_cell_platform:
        project_metadata['single_cell_platform'] = args.single_cell_platform

    # Import extra files specified by the user
    for f in list(extra_files):
        print("Importing %s" % f)
        os.symlink(f,os.path.join(project_dir,os.path.basename(f)))

    # Save out metadata to temporary project dir
    info_file = os.path.join(project_dir,"README.info")
    print("Writing metadata to %s" % info_file)
    project_metadata.save(info_file)

    # Remove the temporary directory on exit
    if not args.no_cleanup:
        print("Registering temporary project directory for "
              "deletion on pipeline completion")
        atexit.register(cleanup_atexit,project_dir)

    # Load the project
    project = AnalysisProject(project_dir)
    print("Loaded project '%s'" % project.name)

    # Set working directory for pipeline
    working_dir = args.working_dir
    if not working_dir:
        working_dir = os.path.join(project_dir,'__run_qc')

    # Set up and run the QC pipeline
    announce("Running QC pipeline")
    runqc = QCPipeline()
    runqc.add_project(project,
                      qc_dir=qc_dir,
                      qc_protocol=args.qc_protocol,
                      report_html=out_file,
                      multiqc=(not args.no_multiqc))
    status = runqc.run(nthreads=nthreads,
                       fastq_screens=fastq_screens,
                       fastq_subset=args.fastq_screen_subset,
                       star_indexes=star_indexes,
                       cellranger_chemistry=\
                       args.cellranger_chemistry,
                       cellranger_force_cells=\
                       args.cellranger_force_cells,
                       cellranger_transcriptomes=cellranger_transcriptomes,
                       cellranger_premrna_references=\
                       cellranger_premrna_references,
                       cellranger_atac_references=cellranger_atac_references,
                       cellranger_arc_references=cellranger_multiome_references,
                       cellranger_jobmode=cellranger_jobmode,
                       cellranger_maxjobs=max_jobs,
                       cellranger_mempercore=cellranger_mempercore,
                       cellranger_jobinterval=cellranger_jobinterval,
                       cellranger_localcores=cellranger_localcores,
                       cellranger_localmem=cellranger_localmem,
                       cellranger_exe=args.cellranger_exe,
                       cellranger_reference_dataset=\
                       cellranger_reference_dataset,
                       cellranger_out_dir=out_dir,
                       max_jobs=max_jobs,
                       max_slots=max_cores,
                       batch_size=args.batch_size,
                       batch_limit=args.max_batches,
                       runners=runners,
                       default_runner=default_runner,
                       envmodules=envmodules,
                       enable_conda=enable_conda,
                       conda_env_dir=args.conda_env_dir,
                       working_dir=working_dir,
                       legacy_screens=use_legacy_screen_names,
                       verbose=args.verbose)
    if status:
        logger.critical("QC failed (see warnings above)")
        if args.force:
            logger.warning("Forcing generation of HTML report")
            report_qc(project,
                      qc_dir=qc_dir,
                      qc_protocol=args.qc_protocol,
                      report_html=out_file,
                      zip_outputs=True,
                      multiqc=(not args.no_multiqc),
                      force=True,
                      runner=runners['report_runner'])

    # Update the QC metadata
    announce("Updating QC metadata")
    qc_info = AnalysisProjectQCDirInfo(filen=os.path.join(qc_dir,
                                                          "qc.info"))
    if qc_info.fastq_dir:
        if master_fastq_dir:
            print("Updating stored Fastq directory for QC: %s" %
                  master_fastq_dir)
        else:
            print("Unsetting stored Fastq directory for QC")
        qc_info['fastq_dir'] = master_fastq_dir
        qc_info.save()
        print("Updated Fastq directory: %s" % (qc_info.fastq_dir
                                               if qc_info.fastq_dir
                                               else '<not set>'))

    # Report locations of final outputs
    announce("QC pipeline completed")
    print("Output directory   : %s" % out_dir)
    print("HTML report        : %s" % out_file)
    print("QC output directory: %s" % qc_dir)
    print("Exit status        : %s" % status)
    # Finish and return exit code from pipeline
    sys.exit(status)
