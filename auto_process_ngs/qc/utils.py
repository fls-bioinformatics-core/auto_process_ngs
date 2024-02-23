#!/usr/bin/env python
#
#     utils: utility classes and functions for QC
#     Copyright (C) University of Manchester 2018-2024 Peter Briggs
#
"""
Provides utility classes and functions for analysis project QC.

Provides the following functions:

- verify_qc: verify the QC run for a project
- report_qc: generate report for the QC run for a project
- get_bam_basename: return the BAM file basename from a Fastq filename
- get_seq_data_samples: identify samples with biological (sequencing)
  data
- set_cell_count_for_project: sets total number of cells for a project
"""

#######################################################################
# Imports
#######################################################################

import os
import logging
from ..analysis import AnalysisFastq
from ..analysis import AnalysisProject
from ..command import Command
from ..metadata import AnalysisProjectQCDirInfo
from ..conda import CondaWrapper
from ..conda import CondaWrapperError
from ..conda import make_conda_env_name
from ..settings import Settings
from ..simple_scheduler import SchedulerJob
from ..tenx.cellplex import CellrangerMultiConfigCsv
from .cellranger import CellrangerCount
from .cellranger import CellrangerMulti

# Module-specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Functions
#######################################################################

def verify_qc(project,qc_dir=None,fastq_dir=None,qc_protocol=None,
              runner=None,log_dir=None):
    """
    Verify the QC run for a project

    Arguments:
      project (AnalysisProject): analysis project
        to verify the QC for
      qc_dir (str): optional, specify the subdir with
        the QC outputs being verified
      fastq_dir (str): optional, specify a non-default
        directory with Fastq files being verified
      qc_protocol (str): optional, QC protocol to
        verify against
      runner (JobRunner): optional, job runner to use
        for running the verification
      log_dir (str): optional, specify a directory to
        write logs to

    Returns:
      Boolean: True if QC passes verification, otherwise
        False.
    """
    # Sort out runners
    if runner is None:
        runner = Settings().general.default_runner
    # Construct command for QC verification
    verify_cmd = Command(
        "reportqc.py",
        "--verify")
    if qc_protocol is not None:
        verify_cmd.add_args("--protocol",qc_protocol)
    if qc_dir is not None:
        verify_cmd.add_args("--qc_dir",qc_dir)
    if fastq_dir is not None:
        verify_cmd.add_args("--fastq_dir",fastq_dir)
    verify_cmd.add_args(project.dirn)
    # Run the command
    verify = SchedulerJob(runner,
                          verify_cmd.command_line,
                          name="verify_qc.%s" % project.name,
                          working_dir=project.dirn,
                          log_dir=log_dir)
    verify.start()
    try:
        verify.wait()
    except KeyboardInterrupt as ex:
        logger.warning("Keyboard interrupt, terminating QC verification")
        verify.terminate()
        raise ex
    # Return boolean based on the exit code
    return (verify.exit_code == 0)

def report_qc(project,qc_dir=None,fastq_dir=None,qc_protocol=None,
              report_html=None,zip_outputs=True,multiqc=False,
              force=False,runner=None,log_dir=None,
              suppress_warning=False):
    """
    Generate report for the QC run for a project

    Arguments:
      project (AnalysisProject): analysis project
        to report the QC for
      qc_dir (str): optional, specify the subdir with
        the QC outputs being reported
      fastq_dir (str): optional, specify a non-default
        directory with Fastq files being verified
      qc_protocol (str): optional, QC protocol to
        verify against
      report_html (str): optional, path to the name of
        the output QC report
      zip_outputs (bool): if True then also generate ZIP
        archive with the report and QC outputs
      multiqc (bool): if True then also generate MultiQC
        report
      force (bool): if True then force generation of
        QC report even if verification fails
      runner (JobRunner): optional, job runner to use
        for running the reporting
      log_dir (str): optional, specify a directory to
        write logs to
      suppress_warning (bool): if True then don't show the
        warning message even when there are missing metrics
        (default: show the warning if there are missing
        metrics)

    Returns:
      Integer: exit code from reporting job (zero indicates
        success, non-zero indicates a problem).
    """
    # Sort out runners
    if runner is None:
        runner = Settings().general.default_runner
    # Basename for the outputs
    if qc_dir is None:
        qc_base = os.path.basename(project.qc_dir)
    else:
        qc_base = os.path.basename(qc_dir)
    # Report HTML file name
    if report_html is None:
        out_file = '%s_report.html' % qc_base
    else:
        out_file = report_html
    if not os.path.isabs(out_file):
        out_file = os.path.join(project.dirn,out_file)
    # Report title
    if project.info.run is None:
        title = "%s" % project.name
    else:
        title = "%s/%s" % (project.info.run,
                           project.name)
    if fastq_dir is not None:
        title = "%s (%s)" % (title,fastq_dir)
    title = "%s: QC report" % title
    # Construct command for reporting
    report_cmd = Command(
        "reportqc.py",
        "--filename",out_file,
        "--title",title)
    if qc_protocol is not None:
        report_cmd.add_args("--protocol",qc_protocol)
    if qc_dir is not None:
        report_cmd.add_args("--qc_dir",qc_dir)
    if fastq_dir is not None:
        report_cmd.add_args("--fastq_dir",fastq_dir)
    if multiqc:
        report_cmd.add_args("--multiqc")
    if zip_outputs:
        report_cmd.add_args("--zip")
    if force:
        report_cmd.add_args("--force")
    if suppress_warning:
        report_cmd.add_args("--suppress-warning")
    report_cmd.add_args(project.dirn)
    # Check if environment modules are defined
    module_load_cmds = None
    if Settings().modulefiles['report_qc']:
        print("Attempting to acquire environment modules for reporting")
        module_load_cmds = []
        try:
            modulepath = os.environ['MODULEPATH']
            if modulepath:
                module_load_cmds.append("export MODULEPATH=%s" % modulepath)
        except KeyError:
            pass
        try:
            envmodules = Settings().modulefiles['report_qc'].split(',')
            for envmodule in envmodules:
                module_load_cmds.append("module load %s" % envmodule)
        except Exception as ex:
            logger.warning("couldn't acquire env modules?: %s" % ex)
        module_load_cmds = '\n'.join(module_load_cmds)
    # Check if conda environments are enabled
    conda_activate_cmd = None
    if Settings().conda.enable_conda:
        print("Attempting to acquire conda environment for reporting")
        # Get location for conda environments
        conda_env_dir = Settings().conda.env_dir
        # Set up conda wrapper
        conda = CondaWrapper(env_dir=conda_env_dir)
        # Get environment for QC reporting
        report_qc_conda_pkgs = ("multiqc=1.8",
                                "pillow",
                                "python=3.8",
                                "numpy=1.20.3")
        env_name = make_conda_env_name(*report_qc_conda_pkgs)
        try:
            conda.create_env(env_name,*report_qc_conda_pkgs)
            conda_env = os.path.join(conda_env_dir,env_name)
            # Script fragment to activate the environment
            conda_activate_cmd = conda.activate_env_cmd(conda_env)
        except CondaWrapperError as ex:
            # Failed to acquire the environment
            logger.warning("failed to acquire conda environment '%s': %s" %
                           (env_name,ex))
    # Wrap the command in a script
    scripts_dir = os.path.join(project.dirn,"ScriptCode")
    if not os.path.isdir(scripts_dir):
        logger.warning("no ScriptCode directory found in '%s'" %
                       project.name)
        scripts_dir = project.dirn
    report_script = os.path.join(scripts_dir,
                                 "report_qc.%s.sh" % project.name)
    prologue = []
    if module_load_cmds:
        prologue.append(module_load_cmds)
    if conda_activate_cmd:
        prologue.append(str(conda_activate_cmd))
    if prologue:
        prologue = '\n'.join(prologue)
    else:
        prologue = None
    report_cmd.make_wrapper_script(filen=report_script,
                                   prologue=prologue,
                                   quote_spaces=True)
    # Locate log dir
    if log_dir is None:
        log_dir = os.path.join(project.dirn,"logs")
        if not os.path.isdir(log_dir):
            log_dir = None
    # Run the command
    report = SchedulerJob(runner,
                          Command('/bin/bash','-l',report_script).command_line,
                          name="report_qc.%s" % project.name,
                          working_dir=project.dirn,
                          log_dir=log_dir)
    report.start()
    try:
        report.wait()
    except KeyboardInterrupt as ex:
        logger.warning("Keyboard interrupt, terminating QC reporting")
        report.terminate()
        raise ex
    # Return the exit code
    return report.exit_code

def get_bam_basename(fastq,fastq_attrs=None):
    """
    Return basename for BAM file from Fastq filename

    Typically this will be the Fastq basename with the
    read ID removed, for example the Fastq filename
    'SM1_S1_L001_R1_001.fastq.gz' will result in the
    BAM basename of 'SM1_S1_L001_001'.

    Arguments:
      fastq (str): Fastq filename; can include leading
        path and extensions (both will be ignored)
      fastq_attrs (BaseFastqAttrs): class for extracting
        data from Fastq names (defaults to 'AnalysisFastq')

    Returns:
      String: basename for BAM file.
    """
    if fastq_attrs is None:
        fastq_attrs = AnalysisFastq
    bam_basename = fastq_attrs(fastq)
    bam_basename.read_number = None
    return str(bam_basename)

def get_seq_data_samples(project_dir,fastq_attrs=None):
    """
    Identify samples with biological (sequencing) data

    Arguments:
      project_dir (str): path to the project directory
      fastq_attrs (BaseFastqAttrs): class for extracting
        data from Fastq names (defaults to 'AnalysisFastq')

    Returns:
      List: list with subset of samples with biological
        data
    """
    # Set up
    if fastq_attrs is None:
        fastq_attrs = AnalysisFastq
    project = AnalysisProject(project_dir,
                              fastq_attrs=fastq_attrs)
    # Initial sample list
    samples = sorted([s.name for s in project.samples])
    # If biological samples explicitly defined in
    # project metadata then use those
    if project.info.biological_samples:
        bio_samples = []
        for s in [str(s).strip()
                  for s in project.info.biological_samples.split(',')]:
            if s not in samples:
                logger.warning("Sample '%s' defined as biological data "
                               "but no sample found with that name?" % s)
            else:
                bio_samples.append(s)
        return bio_samples
    # 10x Genomics CellPlex
    single_cell_platform = project.info.single_cell_platform
    if single_cell_platform:
        if single_cell_platform.startswith("10xGenomics Chromium") and \
           project.info.library_type in ("CellPlex",
                                         "Flex"):
            # CellPlex/Flex
            # Check for a single config file
            config_file = os.path.join(project.dirn,
                                       "10x_multi_config.csv")
            if os.path.exists(config_file):
                config_csv = CellrangerMultiConfigCsv(config_file)
                samples = sorted([s for s in config_csv.gex_libraries
                                  if s in samples])
        elif single_cell_platform.startswith("10xGenomics Chromium") and \
             project.info.library_type == "Single Cell Immune Profiling":
            # Single Cell Immune Profiling
            # Check for multiple config files
            config_files = [os.path.join(project.dirn,f)
                            for f in os.listdir(project.dirn)
                            if (f.startswith("10x_multi_config.") and
                                f.endswith(".csv"))]
            samples_ = []
            for config_file in config_files:
                config_csv = CellrangerMultiConfigCsv(config_file)
                samples_.extend([s for s in config_csv.gex_libraries
                                 if s in samples])
            samples = sorted(samples_)
    return samples

def set_cell_count_for_project(project_dir,qc_dir=None,
                               source="count"):
    """
    Set the total number of cells for a project

    Depending on the specified 'source', sums the number
    of cells for each sample in a project as determined
    from either 'cellranger* count' or 'cellranger multi'.

    Depending the 10x Genomics package and analysis type
    the cell count for individual samples is extracted
    from the 'metrics_summary.csv' file for scRNA-seq
    (i.e. 'cellranger count' or 'cellranger multi'), or
    from the 'summary.csv' file for scATAC (ie.
    'cellranger-atac count').

    The final count is written to the 'number_of_cells'
    metadata item for the project.

    Arguments:
      project_dir (str): path to the project directory
      qc_dir (str): path to QC directory (if not the default
        QC directory for the project)
      source (str): either 'count' or 'multi' (default is
        'count')

    Returns:
      Integer: exit code, non-zero values indicate problems
        were encountered.
    """
    # Set up basic info
    project = AnalysisProject(project_dir)
    if qc_dir is None:
        qc_dir = project.qc_dir
    qc_dir = os.path.abspath(qc_dir)
    print("QC dir: %s" % qc_dir)
    number_of_cells = None
    # Determine which 10x pipeline was used
    pipeline = None
    single_cell_platform = project.info.single_cell_platform
    if single_cell_platform:
        if single_cell_platform.startswith("10xGenomics Chromium 3'"):
            pipeline = "cellranger"
        elif single_cell_platform == "10xGenomics Single Cell ATAC":
            pipeline = "cellranger-atac"
        elif single_cell_platform == "10xGenomics Single Cell Multiome":
            pipeline = "cellranger-arc"
    if not pipeline:
        raise NotImplementedError("Not implemented for platform '%s'"
                                  % single_cell_platform)
    # Fetch information on version and reference data
    cellranger_refdata = None
    cellranger_version = None
    qc_info_file = os.path.join(qc_dir,"qc.info")
    if os.path.exists(qc_info_file):
        qc_info = AnalysisProjectQCDirInfo(filen=qc_info_file)
        try:
            cellranger_refdata = qc_info['cellranger_refdata']
        except KeyError:
            pass
        try:
            cellranger_version = qc_info['cellranger_version']
        except KeyError:
            pass
    else:
        print("%s: not found" % qc_info_file)
    # Determine whether we're handling output from 'multi'
    # or from 'count'
    if source == "multi":
        print("Looking for '%s multi' outputs" % pipeline)
        if not os.path.exists(os.path.join(qc_dir,"cellranger_multi")):
            logger.warning("Unable to set cell count: no data found")
            return
        # Handle outputs from 'multi'
        number_of_cells = 0
        try:
            multi_outs = CellrangerMulti(
                os.path.join(qc_dir,
                             "cellranger_multi",
                             cellranger_version,
                             os.path.basename(
                                 cellranger_refdata)),
                cellranger_exe=pipeline)
            if multi_outs.sample_names:
                for sample in multi_outs.sample_names:
                    print("- %s" % sample)
                    try:
                        ncells = multi_outs.metrics(sample).cells
                        print("  %d cells" % ncells)
                        number_of_cells += ncells
                    except Exception as ex:
                        raise Exception("Failed to add cell count for sample "
                                        "'%s': %s" % (sample,ex))
            else:
                raise Exception("No samples found under %s" %
                                os.path.join(qc_dir,"cellranger_multi"))
        except Exception as ex:
            number_of_cells = None
            logger.warning("Unable to set cell count from data in "
                           "%s: %s" %
                           (os.path.join(qc_dir,"cellranger_multi"),ex))
    elif source == "count":
        print("Looking for '%s count' outputs" % pipeline)
        if not os.path.exists(os.path.join(qc_dir,"cellranger_count")):
            logger.warning("Unable to set cell count: no data found")
            return
        # Handle outputs from 'count'
        # Determine possible locations for outputs
        count_dirs = []
        # New-style with 'version' and 'reference' subdirectories
        if cellranger_version and cellranger_refdata:
            count_dirs.append(os.path.join(qc_dir,
                                           "cellranger_count",
                                           cellranger_version,
                                           os.path.basename(
                                               cellranger_refdata)))
        # Old-style without additional subdirectories
        count_dirs.append(os.path.join(qc_dir,
                                       "cellranger_count"))
        # Check each putative output location in turn
        for count_dir in count_dirs:
            print("Examining %s" % count_dir)
            # Check that the directory exists
            if os.path.exists(count_dir):
                number_of_cells = 0
                try:
                    # Loop over samples and collect cell numbers for
                    # each sample
                    for sample in project.samples:
                        print("- %s" % sample)
                        ncells = None
                        sample_dir = os.path.join(count_dir,
                                                  sample.name)
                        sample_outs = CellrangerCount(sample_dir,
                                                      cellranger_exe=pipeline)
                        for metric in ('Estimated Number of Cells',
                                       'Estimated number of cells',
                                       'annotated_cells',):
                            # Try to fetch metric for cell count
                            try:
                                ncells = sample_outs.metrics.fetch(metric)
                                break
                            except KeyError:
                                pass
                        if ncells is None:
                            number_of_cells = None
                            raise Exception("Failed to add cell count "
                                            "for sample '%s'" % sample.name)
                        else:
                            print("  %d cells" % ncells)
                            number_of_cells += ncells
                    # Extracted cell numbers so break out
                    break
                except Exception as ex:
                    logger.warning("Unable to get cell counts from '%s': %s"
                                   % (count_dir,ex))
    else:
        # No known outputs to get cell counts from
        raise Exception("Unknown source type: '%s'" % source)
    if number_of_cells is not None:
        # Report
        print("Total number of cells: %d" % number_of_cells)
        # Store in the project metadata
        project.info['number_of_cells'] = number_of_cells
        project.info.save()
        return 0
    else:
        # Cell count wasn't set
        return 1
