#!/usr/bin/env python
#
#     analyse_barcodes_cmd.py: implement auto process analyse_barcodes command
#     Copyright (C) University of Manchester 2019 Peter Briggs
#
#########################################################################

#######################################################################
# Imports
#######################################################################

import os
import shutil
import logging
from ..analysis import AnalysisFastq
from ..applications import Command
from ..bcl2fastq_utils import get_nmismatches
from ..bcl2fastq_utils import bases_mask_is_valid
from ..simple_scheduler import SimpleScheduler
from ..utils import write_script_file
from bcftbx.FASTQFile import FastqIterator
from bcftbx import IlluminaData

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Command functions
#######################################################################

def analyse_barcodes(ap,unaligned_dir=None,lanes=None,
                     mismatches=None,cutoff=None,
                     barcode_analysis_dir=None,
                     sample_sheet=None,runner=None,
                     force=False):
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
        it automatically from the bases mask
      cutoff (float): (optional) exclude barcodes with a smaller fraction
        of associated reads than specified cutoff from reporting (e.g.
        '0.001' excludes barcodes with < 0.1% of reads); default is to
        include all barcodes
      sample_sheet (str): if set then use this as the input samplesheet
        to check barcode sequences against (by default will use the
        sample sheet defined in the parameter file for the run)
      barcode_analysis_dir (str): (optional) explicitly specify the
        subdirectory to use for barcode analysis. Counts will be
        written to and read from the 'counts' subdirectory of this
        directory (defaults to 'barcode_analysis')
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
    # Load data
    illumina_data = ap.load_illumina_data(unaligned_dir=unaligned_dir)
    # Handle barcode analysis subdirectories
    if barcode_analysis_dir is not None:
        # Create a subdirectory for barcode analysis
        ap.params['barcode_analysis_dir'] = barcode_analysis_dir
    elif ap.params['barcode_analysis_dir'] is None:
        ap.params['barcode_analysis_dir'] = 'barcode_analysis'
    barcode_analysis_dir = ap.params['barcode_analysis_dir']
    # Create barcode and count file subdirectories
    barcode_dir = ap.add_directory(barcode_analysis_dir)
    counts_dir = os.path.join(barcode_analysis_dir,'counts')
    if os.path.exists(counts_dir) and force:
        print "Removing existing counts data"
        shutil.rmtree(counts_dir)
    ap.add_directory(counts_dir)
    # Map fastq files to counts files
    counts_files = {}
    for project in illumina_data.projects:
        for sample in project.samples:
            for fq in sample.fastq_subset(read_number=1,
                                          full_path=True):
                counts_files[fq] = os.path.join(
                    counts_dir,
                    "%s.%s.counts" % (project.name,
                                      os.path.basename(fq)))
    if illumina_data.undetermined is not None:
        for sample in illumina_data.undetermined.samples:
            for fq in sample.fastq_subset(read_number=1,
                                          full_path=True):
                counts_files[fq] = os.path.join(
                    counts_dir,
                    "undetermined.%s.counts" % os.path.basename(fq))
    # Subset of fastq files with no corresponding counts
    missing_counts = filter(lambda fq: not os.path.exists(counts_files[fq]),
                            counts_files.keys())
    # Deal with lanes
    lane_numbers = illumina_data.lanes
    if lanes:
        lanes = [int(lane) for lane in lanes]
        print "Requested analysis for lanes: %s" % \
            ', '.join([str(l) for l in lanes])
    else:
        print "No lanes explicitly requested"
    if len(lane_numbers) == 1 and lane_numbers[0] is None:
        lane_numbers = None
        print "No lanes explicitly defined in fastq file names"
    else:
        print "Lanes explicitly defined in file names: %s" % \
            ', '.join([str(l) for l in lane_numbers])
    if lanes is None or not lane_numbers:
        # Need counts from all files
        req_counts = counts_files.keys()
    else:
        # Get subset of files explictly belonging to this lane
        req_counts = filter(lambda fq: utils.AnalysisFastq(fq).lane_number
                            in lanes,
                            counts_files.keys())
    # Check there are files to examine
    if not req_counts:
        logger.warning("No matching files: nothing to do")
        return
    # Log dir
    ap.set_log_dir(ap.get_log_subdir('analyse_barcodes'))
    # Set up runner
    if runner is not None:
        runner = fetch_runner(runner)
    else:
        runner = ap.settings.general.default_runner
    runner.set_log_dir(ap.log_dir)
    # Schedule the jobs needed to do counting
    sched = SimpleScheduler(
        runner=runner,
        max_concurrent=ap.settings.general.max_concurrent_jobs,
        poll_interval=ap.settings.general.poll_interval)
    sched.start()
    # Do counting
    print "Getting counts from fastq files"
    group = sched.group("get_barcode_counts")
    for fq in req_counts:
        if fq in missing_counts:
            # Get counts for this file
            barcode_count_cmd = Command(
                'analyse_barcodes.py',
                '-o',os.path.join(ap.analysis_dir,counts_files[fq]),
                '--no-report',fq)
            print "Running %s" % barcode_count_cmd
            group.add(barcode_count_cmd,
                      name='analyse_barcodes.count.%s.%s' %
                      (os.path.basename(counts_files[fq]).split('.')[0],
                       os.path.basename(counts_files[fq]).split('.')[1]))
    group.close()
    # Do reporting
    report_file = os.path.join(barcode_dir,'barcodes.report')
    xls_file = os.path.join(barcode_dir,'barcodes.xls')
    html_file = os.path.join(barcode_dir,'barcodes.html')
    for filen in (report_file,xls_file,html_file):
        if os.path.exists(filen):
            print "Removing existing file: %s" % filen
            os.remove(filen)
    barcode_report_cmd = Command(
        'analyse_barcodes.py',
        '--report',report_file,
        '--xls',xls_file,
        '--html',html_file)
    # Sample sheet
    if sample_sheet is None:
        sample_sheet = ap.params.sample_sheet
    barcode_report_cmd.add_args('--sample-sheet',sample_sheet)
    # Implicitly set per-lane analysis if none were explicitly
    # requested but some are defined in sample sheet
    if lanes is None:
        try:
            lanes = sorted(
                set([str(line['Lane'])
                     for line in IlluminaData.SampleSheet(sample_sheet).data]))
        except KeyError:
            pass
    if lanes:
        barcode_report_cmd.add_args('--lanes',
                                    ','.join([str(l) for l in lanes]))
    # Cutoff
    if cutoff is not None:
        barcode_report_cmd.add_args('--cutoff',cutoff)
    # Mismatches
    if mismatches is None:
        # Try to determine number of mismatches from
        # stored bases mask
        bases_mask = ap.params.bases_mask
        if bases_mask_is_valid(bases_mask):
            mismatches = get_nmismatches(bases_mask)
        else:
            # No valid stored bases mask - try to extract barcode
            # from Fastq header
            logger.warning("Invalid bases mask: '%s'" %
                           bases_mask)
            fq = illumina_data.\
                 projects[0].\
                 samples[0].\
                 fastq_subset(read_number=1,full_path=True)[0]
            print "Extracting index sequence from %s" % fq
            for r in FastqIterator(fq):
                seq_id = r.seqid
                break
            if len(seq_id.index_sequence) >= 6:
                mismatches = 1
            else:
                mismatches = 0
    barcode_report_cmd.add_args('--mismatches',mismatches)
    # Add the list of count files to process
    barcode_report_cmd.add_args('-c')
    for counts_file in [counts_files[f] for f in req_counts]:
        barcode_report_cmd.add_args(os.path.join(ap.analysis_dir,
                                                 counts_file))
    # Write a script file
    script_file = os.path.join(ap.log_dir,'report_barcodes.sh')
    write_script_file(script_file,barcode_report_cmd,
                      shell='/bin/sh')
    # Submit command
    print "Running %s" % script_file
    sched.submit(Command('sh',script_file),
                 name='report_barcodes',
                 wait_for=('get_barcode_counts',))
    # Wait for the scheduler to run all jobs
    sched.wait()
    sched.stop()
    # Finish
    if os.path.exists(report_file):
        print "Report written to %s" % report_file
    else:
        logger.error("Missing barcode analysis report: %s" %
                     report_file)
    if os.path.exists(xls_file):
        print "XLS written to %s" % xls_file
    else:
        logger.error("Missing barcode analysis XLS report: %s" %
                     xls_file)
    if os.path.exists(html_file):
        print "HTML written to %s" % html_file
    else:
        logger.error("Missing barcode analysis HTML report: %s" %
                     html_file)
