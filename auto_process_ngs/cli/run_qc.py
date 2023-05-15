#!/usr/bin/env python
#
#     cli/run_qc.py: command line interface for standalone QC pipeline
#     Copyright (C) University of Manchester 2017-2023 Peter Briggs
#
#########################################################################
#
# run_qc.py
#
#########################################################################

"""
Runs the QC pipeline standalone on an arbitrary set of Fastq files
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
from bcftbx.utils import AttributeDictionary
from .. import get_version
from .. import tenx
from .. import icell8
from ..analysis import AnalysisProject
from ..analysis import AnalysisFastq
from ..analysis import locate_project_info_file
from ..metadata import AnalysisProjectInfo
from ..metadata import AnalysisProjectQCDirInfo
from ..fastq_utils import group_fastqs_by_name
from ..settings import Settings
from ..settings import fetch_reference_data
from ..qc.pipeline import QCPipeline
from ..qc.protocols import fetch_protocol_definition
from ..qc.utils import report_qc

# QC protocols
from ..qc.constants import PROTOCOLS

# 10x Genomics assays
from ..tenx import CELLRANGER_ASSAY_CONFIGS

# Module-specific logger
logger = logging.getLogger("run_qc")

#######################################################################
# Classes
#######################################################################

class InfoAction(argparse.Action):
    """
    Custom parser action for the --info option

    Example usage:

    >>> p.add_argument('--info',action=InfoAction,settings=settings)

    where 'settings' should be a populated 'Settings' instance.

    When invoked the action will display information on protocols,
    organisms and other configuration settings, and then exit.
    """
    def __init__(self,option_strings,settings,nargs=None,*args,**kws):
        self.settings = settings
        if nargs is not None:
            raise ValueError("nargs not allowed")
        super(InfoAction,self).__init__(option_strings=option_strings,
                                        nargs=0,*args,**kws)
    def __call__(self,parser,namespace,values,option_string=None):
        display_info(self.settings)
        sys.exit()

#######################################################################
# Functions
#######################################################################

def add_reporting_options(p):
    """
    Reporting options
    """
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

def add_metadata_options(p):
    """
    Metadara options
    """
    metadata = p.add_argument_group('Metadata')
    metadata.add_argument('--organism',metavar='ORGANISM',
                          action='store',dest='organism',default=None,
                          help="explicitly specify organism (e.g. "
                          "'human', 'mouse'). Multiple organisms "
                          "should be separated by commas (e.g. "
                          "'human,mouse'). HINT use the --info option "
                          "to list the defined organisms")
    metadata.add_argument('--library-type',metavar='LIBRARY',
                          action='store',dest='library_type',default=None,
                          help="explicitly specify library type (e.g. "
                          "'RNA-seq', 'ChIP-seq')")
    metadata.add_argument('--single-cell-platform',metavar='PLATFORM',
                          action='store',dest='single_cell_platform',
                          default=None,
                          help="explicitly specify the single cell "
                          "platform (e.g. '10xGenomics Chromium 3'v3')")

def add_pipeline_options(p,fastq_subset_size,default_nthreads):
    """
    QC pipeline options
    """
    qc_options = p.add_argument_group('QC options')
    qc_options.add_argument('-p','--protocol',metavar='PROTOCOL',
                            action='store',dest='qc_protocol',
                            default=None,choices=PROTOCOLS,
                            help="explicitly specify the QC protocol to "
                            "use; can be one of %s. If not set then "
                            "protocol will be determined automatically "
                            "based on directory contents and metadata." %
                            ", ".join(["'%s'" % x for x in PROTOCOLS]))
    qc_options.add_argument('--fastq_subset',metavar='SUBSET',
                            action='store',dest='fastq_subset',
                            default=fastq_subset_size,
                            type=int,
                            help="specify size of subset of reads "
                            "to use for FastQScreen, strandedness, "
                            "coverage etc option); (default %d, set "
                            "to 0 to use all reads)" % fastq_subset_size)
    qc_options.add_argument('-t','--threads',action='store',
                            dest="nthreads",type=int,default=None,
                            help="number of threads to use for QC script "
                            "(default: %s)" % ('taken from job runner'
                                               if not default_nthreads
                                               else default_nthreads,))

def add_reference_data_options(p):
    """
    Reference data options
    """
    refdata = p.add_argument_group('Reference data')
    refdata.add_argument('--star-index',action='store',
                         metavar="INDEX",
                         dest='star_index',
                         help="specify the path to the STAR genome index "
                         "to use when mapping reads for metrics such as "
                         "strandedness etc (overrides the "
                         "organism-specific indexes defined in the config "
                         "file)")
    refdata.add_argument('--gtf',action='store',
                         metavar="GTF",
                         dest='gtf_annotation',
                         help="specify the path to the GTF annotation "
                         "file to use for metrics such as 'qualimap "
                         "rnaseq' (overrides the organism-specific GTF "
                         "files defined in the config file)")

def add_10x_options(p):
    """
    Cellranger/10x Genomics options
    """
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

def add_conda_options(p,enable_conda,conda_env_dir):
    """
    Conda options
    """
    conda = p.add_argument_group("Conda dependency resolution")
    conda.add_argument('--enable-conda',choices=["yes","no"],
                       dest="enable_conda",default=None,
                       help="use conda to resolve task dependencies; can "
                       "be 'yes' or 'no' (default: %s)" %
                       ("yes" if enable_conda else "no"))
    conda.add_argument('--conda-env-dir',action='store',
                       dest="conda_env_dir",default=conda_env_dir,
                       help="specify directory for conda enviroments "
                       "(default: %s)" % ('temporary directory'
                                          if not conda_env_dir else
                                          conda_env_dir))

def add_job_control_options(p,max_cores,max_jobs,max_batches):
    """
    Job control options
    """
    pipeline = p.add_argument_group('Job control options')
    pipeline.add_argument('--local',action='store_true',
                          dest='local',default=False,
                          help="run the QC on the local system (overrides "
                          "any runners defined in the configuration or on "
                          "the command line)")
    pipeline.add_argument('-c','--maxcores',metavar='N',action='store',
                          dest='max_cores',type=int,default=max_cores,
                          help="maximum number of cores available for QC "
                          "jobs when using --local (default %s, change in "
                          "in settings file)" %
                          (max_cores if max_cores else 'no limit'))
    pipeline.add_argument('-m','--maxmem',metavar='M',action='store',
                          dest='max_mem',type=int,default=None,
                          help="maximum total memory jobs can request at "
                          "once when using --local (in Gbs; default: "
                          "unlimited)")
    pipeline.add_argument('-j','--maxjobs',metavar='N',action='store',
                          dest='max_jobs',type=int,default=max_jobs,
                          help="explicitly specify maximum number of "
                          "concurrent QC jobs to run (default %s, change "
                          "in settings file; ignored when using --local)"
                          % (max_jobs if max_jobs else 'no limit'))
    pipeline.add_argument('-b','--maxbatches',type=int,action='store',
                          dest='max_batches',metavar='NBATCHES',
                          default=max_batches,
                          help="enable dynamic batching of pipeline "
                          "jobs with maximum number of batches set to "
                          "NBATCHES (default: %s)"
                          % (max_batches if max_batches else 'no batching'))

def add_advanced_options(p,use_legacy_screen_names):
    """
    Advanced options
    """
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
                          ("yes" if use_legacy_screen_names else "no"))
    advanced.add_argument('--no-multiqc',action="store_true",
                          dest="no_multiqc",default=False,
                          help="turn off generation of MultiQC report")

def add_debug_options(p):
    """
    Debugging options
    """
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

def add_deprecated_options(p):
    """
    Deprecated options
    """
    deprecated = p.add_argument_group('Deprecated/redundant options')
    deprecated.add_argument('--fastq_screen_subset',metavar='SUBSET',
                            action='store',dest='fastq_screen_subset',
                            help="redundant: use the --fastq_subset "
                            "option instead")
    deprecated.add_argument('--force',action='store_true',
                            help="redundant: HTML report generation will "
                            "always be attempted (even when pipeline "
                            "fails)")
    deprecated.add_argument('--multiqc',action='store_true',
                            dest='run_multiqc', default=False,
                            help="redundant: MultiQC report is generated "
                            "by default (use --no-multiqc to disable)")

def display_info(s):
    """
    Displays information about the current configuration

    The information includes the available QC protocols,
    organisms and FastqScreen conf files.

    Arguments:
      s (Settings): populated Settings instance
    """
    # Location of config file
    print("\nConfig file: {config_file}".format(
        config_file=s.settings_file))
    # QC protocols
    print("\nAvailable QC protocols:")
    if PROTOCOLS:
        for name in PROTOCOLS:
            protocol = fetch_protocol_definition(name)
            print("\t{name:20s} {descr}".format(name=protocol.name,
                                                descr=protocol.description))
    else:
        print("\tNo QC protocols defined")
    # Single cell platforms
    sc_platforms = [p for p in tenx.PLATFORMS]
    sc_platforms.append('Parse Evercode')
    sc_platforms.extend([p for p in icell8.PLATFORMS])
    print("\nSingle cell platforms:")
    if sc_platforms:
        for name in sc_platforms:
            print("\t{platform}".format(platform=name))
    else:
        print("\tNo single cell platforms defined")
    # Organisms
    print("\nOrganisms:")
    if s.organisms:
        for name in s.organisms:
            organism = s.organisms[name]
            print("\t{name}".format(name=name))
    else:
        print("\tNo organisms defined")
    # Reference data
    indexes = {
        'star_index': 'STAR index',
        'annotation_gtf': 'GTF file',
        'annotation_bed': 'BED file',
        'cellranger_reference': 'CellRanger',
        'cellranger_premrna_reference': 'CellRanger (pre-mRNA)',
        'cellranger_atac_reference': 'CellRanger-ATAC',
        'cellranger_arc_reference': 'CellRanger-ARC',
        'cellranger_probe_set': 'Flex probeset',
    }
    if s.organisms:
        print("\nReference data\n==============")
        for name in s.organisms:
            print("\n{name}\n{underline}".format(
                name=name.title().replace('_',' '),
                underline='-'*len(name)))
            no_reference_data = True
            for index in indexes:
                if s.organisms[name][index]:
                    no_reference_data = False
                    print("- {index:15s}: {value}".format(
                        index=indexes[index],
                        value=s.organisms[name][index]))
            if no_reference_data:
                print("- no reference data")
    # Fastq screens
    print("\nFastQScreen\n===========")
    if s.qc.fastq_screens:
        for name in s.qc.fastq_screens.split(','):
            conf_file = s.screens[name].conf_file
            print("- {name:15s}: {conf_file}".format(name=name,
                                                    conf_file=conf_file))
    else:
        print("No screens defined")

def process_inputs(input_list):
    """
    Process the inputs and return Fastqs etc

    The inputs can be one of:

    - a subdirectory in a project
    - a project directory
    - a non-project directory with Fastqs
    - a 'raw' list of Fastqs

    The function attempts to determine which type the inputs
    are, generate a list of Fastq files, and locate any related
    filesystem objects (for example a "parent" project directory).

    It returns a dictionary-like object with the following
    elements:

    - 'fastqs': a list of Fastq files
    - 'dir_path': the directory supplied as an input, if any
    - 'info_file': path to an AnalysisProject metadata file
    - 'extra_files': any additional QC-related config files

    Returns:
      AttributeDictionary: elements are 'fastqs', 'dir_path',
        'info_file' and 'extra_files'
    """
    # Initialise
    inputs = []
    dir_path = None
    info_file = None
    extra_files = set()
    # Process list of inputs
    for f in input_list:
        for ff in glob.glob(os.path.abspath(f)):
            if not os.path.exists(ff):
                # Input not found
                logger.fatal("%s: input not found" % ff)
                sys.exit(1)
            elif os.path.isdir(ff) and len(input_list) > 1:
                # Can only be a single directory
                logger.fatal("Input must be a single directory, or a list of "
                             "Fastqs")
                sys.exit(1)
            else:
                inputs.append(ff)
    # Get list of Fastqs from directory
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
    # Return processed inputs
    return AttributeDictionary(
        fastqs=inputs,
        dir_path=dir_path,
        info_file=info_file,
        extra_files=extra_files)

def get_execution_environment():
    """
    Fetch information on the local execution environment

    Interrogates the local system to get information
    on number of cores, memory etc.

    It returns a dictionary-like object with the following
    elements:

    - 'cpu_count': total number of CPUs
    - 'total_mem': total amount of memory (Gb)
    - 'nslots': value of the 'NSLOTS' env variable
    - 'max_cores': maximum available cores
    - 'max_mem': maximum available memory (Gb)
    - 'mem_per_core': memory per core (Gb)

    Available cores is the number of CPUs, or the
    value of 'NSLOTS' if set. Available memory is
    the proportion of total memory scaled by the
    number of available cores. Memory per core
    is the total memory divided by the total number
    of CPUs.

    Returns:
      AttributeDictionary: elements are 'cpu_count',
        'total_mem', 'nslots', 'max_cores', 'max_mem'
        and 'mem_per_core'
    """
    # Get numbers of cores
    cpu_count = psutil.cpu_count()
    try:
        # If NSLOTS is set in the environment
        # then assume we're running on an SGE
        # node and this sets the maximum number
        # of available cores
        max_cores = int(os.environ['NSLOTS'])
        nslots = max_cores
    except KeyError:
        # Set limit from local machine
        max_cores = cpu_count
        nslots = None
    # Get total memory (in Gbs)
    total_mem = float(psutil.virtual_memory().total)/(1024.0**3)
    # Memory per core
    mem_per_core = total_mem/float(cpu_count)
    # Maximum available memory
    # (memory per core times the number of cores
    # being used)
    max_mem = mem_per_core*float(max_cores)
    # Return the values
    return AttributeDictionary(
        cpu_count=cpu_count,
        total_mem=total_mem,
        nslots=nslots,
        max_cores=max_cores,
        max_mem=max_mem,
        mem_per_core=mem_per_core)

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

# Main program

def main():
    # Get configuration settings
    settings = Settings()

    # Make a command line parser
    p = argparse.ArgumentParser(
        description="Run the QC pipeline standalone on an arbitrary "
        "set of Fastq files.")
    # Build parser
    p.add_argument('--version', action='version',
                   version=("%%(prog)s %s" % get_version()))
    p.add_argument('--info',action=InfoAction,settings=settings,
                   help="display information on protocols, organisms "
                   "and other settings (then exit)")
    p.add_argument("inputs",metavar="DIR | FASTQ [ FASTQ ... ]",
                   nargs="+",
                   help="directory or list of Fastq files to run the "
                   "QC on")
    # Add option groups
    add_reporting_options(p)
    add_metadata_options(p)
    add_pipeline_options(p,
                         fastq_subset_size=settings.qc.fastq_subset_size,
                         default_nthreads=settings.qc.nprocessors)
    add_reference_data_options(p)
    add_10x_options(p)
    add_conda_options(p,
                      enable_conda=settings.conda.enable_conda,
                      conda_env_dir=settings.conda.env_dir)
    add_job_control_options(p,
                            max_cores=settings.general.max_cores,
                            max_jobs=settings.general.max_concurrent_jobs,
                            max_batches=settings.general.max_batches)
    add_advanced_options(p,use_legacy_screen_names=
                         settings.qc.use_legacy_screen_names)
    add_debug_options(p)
    add_deprecated_options(p)

    # Parse the command line
    args = p.parse_args()

    # Check for deprecated and unsupported options
    if args.fastq_screen_subset:
        logger.fatal("'--fastq_screen_subset' is redundant; use "
                     "'--fastq_subset' instead")
        sys.exit(1)
    if args.force:
        logger.warning("'--force' option is redundant; HTML report "
                       "generation will always be attempted (even "
                       "if pipeline fails)")
    if args.run_multiqc:
        logger.warning("'--multiqc' option is redundant; MultiQC is "
                       "run by default")
        logger.warning("Use '--no-multiqc' to turn off MultiQC report "
                       "generation")

    # Initialise
    project_metadata = AnalysisProjectInfo()
    out_dir = args.out_dir
    qc_dir = args.qc_dir
    master_fastq_dir = None
    fastq_attrs = AnalysisFastq

    # Deal with inputs
    announce("Locating inputs")
    inputs = process_inputs(args.inputs)
    # Filter out index reads
    inputs.fastqs = [fq for fq in inputs.fastqs
                     if not fastq_attrs(fq).is_index_read]
    if not inputs.fastqs:
        logger.fatal("No Fastqs found")
        sys.exit(1)

    # Report what was found
    for fqs in group_fastqs_by_name(inputs.fastqs,fastq_attrs=fastq_attrs):
        print("%s:" % fastq_attrs(fqs[0]).sample_name)
        for fq in fqs:
            print("  %s" % fq)
    print("Located %s Fastq%s" % (len(inputs.fastqs),
                                  's' if len(inputs.fastqs) != 1 else ''))
    if inputs.info_file:
        print("Located project metadata in %s%s" % (inputs.info_file,
                                                    " (will be ignored)"
                                                    if args.ignore_metadata
                                                    else ''))
    else:
        print("Unable to locate project metadata")

    # Per task environment modules
    envmodules = dict()
    for name in ('fastqc',
                 'fastq_screen',
                 'fastq_strand',
                 'cellranger',
                 'report_qc',):
        try:
            envmodules[name] = settings.modulefiles[name]
        except KeyError:
            envmodules[name] = None

    # Maximum number of jobs and cores
    max_jobs = args.max_jobs
    max_cores = args.max_cores

    # Conda dependency resolution
    enable_conda = args.enable_conda
    if enable_conda is None:
        enable_conda = settings.conda.enable_conda
    else:
        enable_conda = (args.enable_conda == "yes")

    # Fastq screens
    if settings.qc.fastq_screens:
        fastq_screens = dict()
        for screen in settings.qc.fastq_screens.split(','):
            fastq_screens[screen] = settings.screens[screen].conf_file
    else:
        fastq_screens = None
    use_legacy_screen_names = args.use_legacy_screen_names
    if use_legacy_screen_names is None:
        use_legacy_screen_names = settings.qc.use_legacy_screen_names
    else:
        use_legacy_screen_names = (use_legacy_screen_names == "yes")

    # STAR indexes
    star_indexes = fetch_reference_data(settings,'star_index')
    force_star_index = args.star_index
    if force_star_index:
        force_star_index = os.path.abspath(force_star_index)

    # Annotation files
    annotation_bed_files = fetch_reference_data(settings,
                                                'annotation_bed')
    annotation_gtf_files = fetch_reference_data(settings,
                                                'annotation_gtf')
    force_gtf_annotation = args.gtf_annotation
    if force_gtf_annotation:
        force_gtf_annotation = os.path.abspath(force_gtf_annotation)

    # Cellranger settings
    cellranger_settings = settings['10xgenomics']

    # Cellranger reference datasets
    cellranger_transcriptomes = fetch_reference_data(
        settings,
        'cellranger_reference')
    cellranger_premrna_references = fetch_reference_data(
        settings,
        'cellranger_premrna_reference')
    cellranger_atac_references = fetch_reference_data(
        settings,
        'cellranger_atac_reference')
    cellranger_multiome_references = fetch_reference_data(
        settings,
        'cellranger_arc_reference')

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
        local_env = get_execution_environment()
        if not max_cores:
            max_cores = local_env.max_cores
        if args.max_mem:
            max_mem = args.max_mem
        else:
            max_mem = local_env.max_mem
        mempercore = local_env.mem_per_core
        print("-- Maximum cores: %d" % max_cores)
        print("-- Maximum memory: %.1f Gbs" % max_mem)
        print("-- Mem per core: %.1f Gbs" % mempercore)
        # Set the default threads for different jobs
        nthreads = min(max_cores,8)
        nthreads_star = min(max_cores,
                                int(math.ceil(32.0/mempercore)))
        ncores_picard = min(max_cores,
                              int(math.ceil(4.0/mempercore)*2))
        ncores_qualimap = min(max_cores,
                              int(math.ceil(4.0/mempercore)*2))
        # Override if nthreads was explicitly set
        # on the command line
        if args.nthreads:
            nthreads = args.nthreads
            nthreads_star = args.nthreads
            ncores_picard = args.nthreads
            ncores_qualimap = args.nthreads
        print("-- Threads for QC: %s" % nthreads)
        print("-- Threads for STAR: %s" % nthreads_star)
        if nthreads_star*mempercore < 32.0:
            logger.warning("Insufficient memory for STAR?")
        print("-- Cores for Qualimap: %s" % ncores_qualimap)
        if ncores_qualimap*mempercore < 8.0:
            logger.warning("Insufficient memory for Qualimap?")
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
        print("-- Cellranger localcores: %d" % cellranger_localcores)
        print("-- Cellranger localmem  : %.1f Gbs" % cellranger_localmem)
        # Set up local runners
        default_runner = SimpleJobRunner()
        runners = {
            'cellranger_count_runner':
            SimpleJobRunner(nslots=cellranger_localcores),
            'cellranger_multi_runner':
            SimpleJobRunner(nslots=cellranger_localcores),
            'fastqc_runner': SimpleJobRunner(nslots=nthreads),
            'fastq_screen_runner': SimpleJobRunner(nslots=nthreads),
            'picard_runner': SimpleJobRunner(nslots=ncores_picard),
            'qualimap_runner': SimpleJobRunner(nslots=ncores_qualimap),
            'rseqc_runner': SimpleJobRunner(),
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
            nthreads = settings.qc.nprocessors
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
                'cellranger_count_runner': default_runner,
                'cellranger_multi_runner': default_runner,
                'fastqc_runner': default_runner,
                'fastq_screen_runner': default_runner,
                'picard_runner': default_runner,
                'qualimap_runner': default_runner,
                'rseqc_runner': default_runner,
                'star_runner': default_runner,
                'verify_runner': default_runner,
                'report_runner': default_runner,
            }
        else:
            # Runners from configuration
            print("Setting up runners from configuration")
            default_runner = settings.general.default_runner
            runners = {
                'cellranger_count_runner': settings.runners.cellranger_count,
                'cellranger_multi_runner': settings.runners.cellranger_multi,
                'fastqc_runner': settings.runners.fastqc,
                'fastq_screen_runner': settings.runners.fastq_screen,
                'picard_runner': settings.runners.picard,
                'qualimap_runner': settings.runners.qualimap,
                'rseqc_runner': settings.runners.rseqc,
                'star_runner': settings.runners.star,
                'verify_runner': default_runner,
                'report_runner': default_runner,
            }

    # Output directory
    announce("Setting up output destinations")
    if not out_dir:
        if inputs.dir_path:
            out_dir = inputs.dir_path
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
    for fq in inputs.fastqs:
        # Make symlinks to the Fastq files
        os.symlink(fq,os.path.join(fastq_dir,os.path.basename(fq)))

    # Set up metadata
    if inputs.info_file and not args.ignore_metadata:
        project_metadata.load(inputs.info_file)
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
    for f in list(inputs.extra_files):
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
    project = AnalysisProject(project_dir,fastq_attrs=fastq_attrs)
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
                       fastq_subset=args.fastq_subset,
                       star_indexes=star_indexes,
                       annotation_bed_files=annotation_bed_files,
                       annotation_gtf_files=annotation_gtf_files,
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
                       force_star_index=force_star_index,
                       force_gtf_annotation=force_gtf_annotation,
                       legacy_screens=use_legacy_screen_names,
                       verbose=args.verbose)

    # Report if QC pipeline failed
    if status:
        logger.critical("QC failed (see warnings above)")
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
