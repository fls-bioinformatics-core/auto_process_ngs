#!/usr/bin/env python
#
#     analyse_barcodes_cmd.py: implement analyse_barcodes command
#     Copyright (C) University of Manchester 2019,2021 Peter Briggs
#
#########################################################################

#######################################################################
# Imports
#######################################################################

import os
import logging
from ..barcodes.pipeline import AnalyseBarcodes

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Command functions
#######################################################################

def analyse_barcodes(ap,unaligned_dir=None,lanes=None,
                     mismatches=None,cutoff=None,
                     barcode_analysis_dir=None,
                     sample_sheet=None,name=None,
                     runner=None,force=False):
    """Analyse the barcode sequences for Fastqs for each specified lane

    Run 'analyse_barcodes.py' for one or more lanes, to analyse the
    barcode index sequences in each lane.
    
    Arguments:
      ap (AutoProcessor): autoprocessor pointing to the analysis
        directory to create Fastqs for
      unaligned_dir (str): if set then use this as the output directory
        for bcl-to-fastq conversion. Default is 'bcl2fastq' (unless
        an alternative is already specified in the config file)
      lanes (list): (optional) specify a list of lane numbers to
        use in the processing; lanes not in the list will be excluded
        (default is to include all lanes)
      mismatches (int): (optional) maximum number of mismatches to
        consider when grouping similar barcodes; default is to determine
        it automatically
      cutoff (float): (optional) exclude barcodes with a smaller fraction
        of associated reads than specified cutoff from reporting (e.g.
        '0.001' excludes barcodes with < 0.1% of reads); default is to
        include all barcodes
      barcode_analysis_dir (str): (optional) explicitly specify the
        subdirectory to use for barcode analysis. Counts will be
        written to and read from the 'counts' subdirectory of this
        directory (defaults to 'barcode_analysis')
      sample_sheet (str): if set then use this as the input samplesheet
        to check barcode sequences against (by default will use the
        sample sheet defined in the parameter file for the run)
      name (str): (optional) identifier for output directory (if
        'barcode_analysis_dir' not explicitly set) and report title
      runner (JobRunner): (optional) specify a non-default job runner
        to use for barcode analysis
      force (bool): if True then forces regeneration of any existing
        counts (default is to reuse existing counts)
    """
    # Sort out parameters
    if unaligned_dir is not None:
        ap.params['unaligned_dir'] = unaligned_dir
    elif ap.params['unaligned_dir'] is None:
        ap.params['unaligned_dir'] = 'bcl2fastq'
    bases_mask = ap.params['bases_mask']
    if sample_sheet is None:
        sample_sheet = ap.params.sample_sheet
    if not os.path.isabs(sample_sheet):
        sample_sheet = os.path.join(ap.params['analysis_dir'],
                                    sample_sheet)
    # Handle barcode analysis subdirectory
    if barcode_analysis_dir is not None:
        # Create a subdirectory for barcode analysis
        ap.params['barcode_analysis_dir'] = barcode_analysis_dir
    elif name is not None:
        ap.params['barcode_analysis_dir'] = 'barcode_analysis_%s' % name
    elif ap.params['barcode_analysis_dir'] is None:
        ap.params['barcode_analysis_dir'] = 'barcode_analysis'
    barcode_analysis_dir = ap.params['barcode_analysis_dir']
    if not os.path.isabs(barcode_analysis_dir):
        barcode_analysis_dir = os.path.join(ap.params['analysis_dir'],
                                            barcode_analysis_dir)
    # Report title
    title = "Barcode analysis for %s" % ap.metadata.run_name
    if name:
        title = title + " (%s)" % name
    # Create a pipeline for barcode analysis
    barcode_analysis = AnalyseBarcodes(os.path.join(
        ap.params['analysis_dir'],ap.params['unaligned_dir']))
    # Log dir and log file
    ap.set_log_dir(ap.get_log_subdir('analyse_barcodes'))
    log_file = os.path.join(ap.log_dir,"analyse_barcodes.log")
    # Set up runner
    if runner is None:
        runner = ap.settings.general.default_runner
    runner.set_log_dir(ap.log_dir)
    # Get scheduler parameters
    max_jobs = ap.settings.general.max_concurrent_jobs
    poll_interval = ap.settings.general.poll_interval
    # Run the pipeline
    status = barcode_analysis.run(
        barcode_analysis_dir,
        title=title,
        lanes=lanes,
        mismatches=mismatches,
        cutoff=cutoff,
        sample_sheet=sample_sheet,
        force=force,
        log_file=log_file,
        runner=runner,
        max_jobs=max_jobs,
        poll_interval=poll_interval,
        verbose=False)
    # Finish
    if status == 0:
        print("Report written to %s" % barcode_analysis.output.report_file)
        print("XLS written to %s" % barcode_analysis.output.xls_file)
        print("HTML written to %s" % barcode_analysis.output.html_file)
    return status
