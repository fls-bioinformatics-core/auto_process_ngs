#!/usr/bin/env python
#
#     update_fastq_stats_cmd.py: implement update_fastq_stats command
#     Copyright (C) University of Manchester 2019-2020 Peter Briggs
#
#########################################################################

#######################################################################
# Imports
#######################################################################

import os
from bcftbx import IlluminaData
from ..applications import Command
from ..bcl2fastq.reporting import ProcessingQCReport
from ..simple_scheduler import SchedulerJob
import logging

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Command functions
#######################################################################

def update_fastq_stats(ap,sample_sheet=None,name=None,stats_file=None,
                       per_lane_stats_file=None,unaligned_dir=None,
                       add_data=False,force=False,nprocessors=None,
                       runner=None):
    """Update statistics for Fastq files

    Updates the statistics for all Fastq files found in the
    'unaligned' directory, by running the 'fastq_statistics.py'
    program.

    Arguments
      ap (AutoProcessor): autoprocessor pointing to the analysis
        directory to create Fastqs for
      sample_sheet (str): path to sample sheet file used in
        bcl-to-fastq conversion (defaults to the sample sheet
        stored in the analysis directory parameters)
      name (str): identifier to use for output stats files
      stats_file (str): path of a non-default file to write the
        statistics to (defaults to 'statistics.info' unless
        over-ridden by local settings)
      per_lane_stats_file (str): path for per-lane statistics
        output file (defaults to 'per_lane_statistics.info'
        unless over-ridden by local settings)
      unaligned_dir (str): output directory for bcl-to-fastq
        conversion
      add_data (bool): if True then add stats to the existing
        stats files (default is to overwrite existing stats
        files)
      force (bool): if True then force update of the stats
        files even if they are newer than the Fastq files
        (by default stats are only updated if they are older
        than the Fastqs)
      nprocessors (int): number of cores to use when running
        'fastq_statistics.py'
      runner (JobRunner): (optional) specify a non-default job
        runner to use for running 'fastq_statistics.py'
    """
    ap.set_log_dir(ap.get_log_subdir("update_fastq_stats"))
    fastq_statistics(ap,
                     unaligned_dir=unaligned_dir,
                     sample_sheet=sample_sheet,
                     name=name,
                     stats_file=stats_file,
                     per_lane_stats_file=per_lane_stats_file,
                     add_data=add_data,
                     force=force,
                     nprocessors=nprocessors,
                     runner=runner)

def fastq_statistics(ap,stats_file=None,per_lane_stats_file=None,
                     unaligned_dir=None,sample_sheet=None,
                     name=None,add_data=False,force=False,
                     nprocessors=None,runner=None):
    """Generate statistics for Fastq files

    Generates statistics for all Fastq files found in the
    'unaligned' directory, by running the 'fastq_statistics.py'
    program.

    Arguments
      ap (AutoProcessor): autoprocessor pointing to the analysis
        directory to create Fastqs for
      stats_file (str): path of a non-default file to write the
        statistics to (defaults to 'statistics.info' unless
        over-ridden by local settings)
      per_lane_stats_file (str): path for per-lane statistics
        output file (defaults to 'per_lane_statistics.info'
        unless over-ridden by local settings)
      unaligned_dir (str): output directory for bcl-to-fastq
        conversion
      sample_sheet (str): path to sample sheet file used in
        bcl-to-fastq conversion
      name (str): identifier to use for output stats files
      add_data (bool): if True then add stats to the existing
        stats files (default is to overwrite existing stats
        files)
      force (bool): if True then force update of the stats
        files even if they are newer than the Fastq files
        (by default stats are only updated if they are older
        than the Fastqs)
      nprocessors (int): number of cores to use when running
        'fastq_statistics.py'
      runner (JobRunner): (optional) specify a non-default job
        runner to use for running 'fastq_statistics.py'
    """
    # Get file names for output files
    if stats_file is None:
        if not name:
            if ap.params['stats_file'] is not None:
                stats_file = ap.params['stats_file']
            else:
                stats_file = 'statistics.info'
        else:
            stats_file = 'statistics.%s.info' % name
    if per_lane_stats_file is None:
        if not name:
            if ap.params['per_lane_stats_file'] is not None:
                per_lane_stats_file = ap.params['per_lane_stats_file']
            else:
                per_lane_stats_file = 'per_lane_statistics.info'
        else:
            per_lane_stats_file = 'per_lane_statistics.%s.info' % name
    if not name:
        per_lane_sample_stats_file = 'per_lane_sample_stats.info'
        full_stats_file = 'statistics_full.info'
    else:
        per_lane_sample_stats_file = 'per_lane_sample_stats.%s.info' % name
        full_stats_file = 'statistics_full.%s.info' % name
    # Sort out unaligned_dir
    if unaligned_dir is None:
        if ap.params.unaligned_dir is None:
            ap.params['unaligned_dir'] = 'bcl2fastq'
        unaligned_dir = ap.params.unaligned_dir
    if not os.path.exists(os.path.join(ap.params.analysis_dir,unaligned_dir)):
        logger.error("Unaligned dir '%s' not found" % unaligned_dir)
    # Check for sample sheet
    if sample_sheet is None:
        sample_sheet = ap.params['sample_sheet']
    if not os.path.exists(sample_sheet):
        logger.error("Sample sheet '%s' not found" % sample_sheet)
        raise Exception("Missing sample sheet")
    # Check if any Fastqs are newer than stats files
    newest_mtime = 0
    for f in (stats_file,per_lane_stats_file,):
        try:
            newest_mtime = max(newest_mtime,
                               os.path.getmtime(f))
        except OSError:
            # Missing file
            newest_mtime = 0
            break
    illumina_data = IlluminaData.IlluminaData(ap.params.analysis_dir,
                                              unaligned_dir)
    if newest_mtime > 0:
        regenerate_stats = False
        for project in illumina_data.projects:
            for sample in project.samples:
                for fq in sample.fastq:
                    if (os.path.getmtime(os.path.join(sample.dirn,fq)) >
                        newest_mtime):
                        regenerate_stats = True
                        break
        if regenerate_stats:
            logger.warning("Fastqs are newer than stats files")
        else:
            logger.warning("Stats files are newer than Fastqs")
            if not force:
                # Don't rerun the stats, just regenerate the report
                processing_qc_html = os.path.join(ap.analysis_dir,
                                                  "processing_qc.html")
                report_processing_qc(ap,processing_qc_html)
                return
    # Set up runner
    if runner is None:
        runner = ap.settings.runners.stats
    runner.set_log_dir(ap.log_dir)
    # Number of cores
    if nprocessors is None:
        if ap.settings.fastq_stats.nprocessors:
            nprocessors = ap.settings.fastq_stats.nprocessors
        else:
            nprocessors = runner.nslots
    # Generate statistics
    fastq_statistics_cmd = Command(
        'fastq_statistics.py',
        '--unaligned',unaligned_dir,
        '--sample-sheet',sample_sheet,
        '--output',os.path.join(ap.params.analysis_dir,stats_file),
        '--per-lane-stats',os.path.join(ap.params.analysis_dir,
                                        per_lane_stats_file),
        '--per-lane-sample-stats',os.path.join(ap.params.analysis_dir,
                                               per_lane_sample_stats_file),
        '--full-stats',os.path.join(ap.params.analysis_dir,
                                    full_stats_file),
        ap.params.analysis_dir,
        '--nprocessors',nprocessors
    )
    if add_data:
        fastq_statistics_cmd.add_args('--update')
    print("Generating statistics: running %s" % fastq_statistics_cmd)
    fastq_statistics_job = SchedulerJob(
        runner,
        fastq_statistics_cmd.command_line,
        name='fastq_statistics',
        working_dir=ap.analysis_dir
    )
    fastq_statistics_job.start()
    try:
        fastq_statistics_job.wait(
            poll_interval=ap.settings.general.poll_interval
        )
    except KeyboardInterrupt as ex:
        logger.warning("Keyboard interrupt, terminating fastq_statistics")
        fastq_statistics_job.terminate()
        raise ex
    exit_code = fastq_statistics_job.exit_code
    print("fastq_statistics completed: exit code %s" % exit_code)
    if exit_code != 0:
        raise Exception("fastq_statistics exited with an error")
    ap.params['stats_file'] = stats_file
    ap.params['per_lane_stats_file'] = per_lane_stats_file
    print("Statistics generation completed: %s" % ap.params.stats_file)
    print("Generating processing QC report")
    processing_qc_html = os.path.join(ap.analysis_dir,
                                      "processing_qc.html")
    report_processing_qc(ap,processing_qc_html)

def report_processing_qc(ap,html_file):
    """
    Generate HTML report for processing statistics

    Arguments:
      ap (AutoProcess): AutoProcess instance to report the
        processing from
      html_file (str): destination path and file name for
        HTML report
    """
    # Per-lane statistics
    per_lane_stats_file = ap.params.per_lane_stats_file
    if per_lane_stats_file is None:
        per_lane_stats_file = "per_lane_statistics.info"
    per_lane_stats_file = get_absolute_file_path(per_lane_stats_file,
                                                 base=ap.analysis_dir)
    # Per lane by sample statistics
    per_lane_sample_stats_file = get_absolute_file_path(
        "per_lane_sample_stats.info",
        base=ap.analysis_dir)
    # Per fastq statistics
    stats_file = get_absolute_file_path("statistics_full.info",
                                        base=ap.analysis_dir)
    if not os.path.exists(stats_file):
        if ap.params.stats_file is not None:
            stats_file = ap.params.stats_file
        else:
            stats_file = "statistics.info"
    stats_file = get_absolute_file_path(stats_file,
                                        base=ap.analysis_dir)
    # Generate the report
    ProcessingQCReport(ap.analysis_dir,
                       stats_file,
                       per_lane_stats_file,
                       per_lane_sample_stats_file).write(html_file)

def get_absolute_file_path(p,base=None):
    """
    Get absolute path for supplied path

    Arguments:
      p (str): path
      base (str): optional, base directory to
        use if p is relative

    Returns:
      String: absolute path for p.
    """
    if not os.path.isabs(p):
        if base is not None:
            p = os.path.join(os.path.abspath(base),p)
        else:
            p = os.path.abspath(p)
    return p
