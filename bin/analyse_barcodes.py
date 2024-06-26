#!/usr/bin/env python
#
#     analyse_barcodes.py: analyse index sequences from Illumina FASTQs
#     Copyright (C) University of Manchester 2016-2023 Peter Briggs
#
"""
analyse_barcodes.py

Counts the index sequences (i.e. barcodes) for all reads in one or more
Fastq files, and reports the most numerous. Will also check against the
sequences supplied in a SampleSheet file and highlight those that are
missing or which have very low counts.

"""

import argparse
import sys
import os
import tempfile
import logging
from bcftbx.IlluminaData import IlluminaData
from bcftbx.IlluminaData import IlluminaDataError
from bcftbx.IlluminaData import SampleSheet
from bcftbx.IlluminaData import samplesheet_index_sequence
from bcftbx.FASTQFile import FastqIterator
from bcftbx.utils import parse_lanes
from auto_process_ngs.barcodes.analysis import BarcodeCounter
from auto_process_ngs.barcodes.analysis import Reporter
from auto_process_ngs.barcodes.analysis import report_barcodes
from auto_process_ngs.bcl2fastq.utils import make_custom_sample_sheet
from auto_process_ngs.bcl2fastq.utils import check_barcode_collisions
from auto_process_ngs.fastq_utils import group_fastqs_by_name
from auto_process_ngs.tenx.utils import has_10x_indices

from auto_process_ngs import get_version
__version__ = get_version()

#######################################################################
# Functions
#######################################################################

def count_barcodes_bcl2fastq(dirn):
    """
    Count the barcodes from bcl2fastq output

    """
    try:
        unaligned = os.path.basename(dirn.rstrip(os.sep))
        dirn = os.path.dirname(os.path.abspath(dirn.rstrip(os.sep)))
        illumina_data = IlluminaData(dirn,unaligned_dir=unaligned)
    except IlluminaDataError as ex:
        print("%s: not an Illumina output directory?" % dirn)
        raise ex
    fqs = []
    for p in illumina_data.projects:
        for s in p.samples:
            for fq in s.fastq_subset(read_number=1,full_path=True):
                fqs.append(fq)
    if illumina_data.undetermined:
        for s in illumina_data.undetermined.samples:
            for fq in s.fastq_subset(read_number=1,full_path=True):
                fqs.append(fq)
    return count_barcodes(fqs)

def count_barcodes(fastqs):
    """
    Count the barcodes from multiple fastqs

    """
    print("Reading in %s fastq%s" % (len(fastqs),
                                     ('' if len(fastqs) == 1
                                      else 's')))
    counts = BarcodeCounter()
    for fq in fastqs:
        print("%s" % os.path.basename(fq))
        for r in FastqIterator(fq):
            seq = r.seqid.index_sequence
            lane = int(r.seqid.flowcell_lane)
            counts.count_barcode(seq,lane)
    return counts

def count_sequences(fastqs,start=None,end=None):
    """
    Count sequences (instead of barcodes) from multiple fastqs
    """
    print("Reading in %s fastq%s" % (len(fastqs),
                                     ('' if len(fastqs) == 1
                                      else 's')))
    print("Counting sequences instead of barcodes")
    if start:
        print("Start position: %d" % start)
        # Correct for Python's zero-based counting
        # i.e. first base is actually zero
        start -= 1
    if end:
        print("End position: %d" % end)
    counts = BarcodeCounter()
    for fq in fastqs:
        print("%s" % os.path.basename(fq))
        for r in FastqIterator(fq):
            seq = r.sequence[start:end]
            lane = int(r.seqid.flowcell_lane)
            counts.count_barcode(seq,lane)
    return counts

# Main program
if __name__ == '__main__':
    p = argparse.ArgumentParser(usage=
                                "\n\t%(prog)s FASTQ [FASTQ...]\n"
                                "\t%(prog)s DIR\n"
                                "\t%(prog)s -c COUNTS_FILE [COUNTS_FILE...]",
                                description="Collate and report counts and "
                                "statistics for Fastq index sequences (aka "
                                "barcodes). If multiple Fastq files are "
                                "supplied then sequences will be pooled "
                                "before being analysed. If a single "
                                "directory is supplied then this will be "
                                "assumed to be an output directory from "
                                "bcl2fastq and files will be processed on a "
                                "per-lane basis. If the -c option is "
                                "supplied then the input must be one or more "
                                "file of barcode counts generated previously "
                                "using the -o option.")
    p.add_argument('-v','--version',action='version',
                   version="%%(prog)s %s" % __version__)
    input_and_output = p.add_argument_group("Input and output options")
    input_and_output.add_argument(
        '-c','--counts',
        action='store_true',dest='use_counts',default=False,
        help="input is one or more counts files generated by "
        "previous runs using the '-o/--output' option")
    input_and_output.add_argument(
        '-o','--output',
        action='store',dest='counts_file_out',default=None,
        help="output all counts to tab-delimited file "
        "COUNTS_FILE_OUT. This can be used again in another "
        "run by specifying the '-c' option")
    reporting = p.add_argument_group("Reporting options")
    reporting.add_argument(
        '-l','--lanes',action='store',dest='lanes',default=None,
        help="restrict analysis to the specified lane numbers "
        "(default is to process all lanes). Multiple lanes "
        "can be specified using ranges (e.g. '2-3'), comma-"
        "separated list ('5,7') or a mixture ('2-3,5,7')")
    reporting.add_argument(
        '-m','--mismatches',action='store',dest='mismatches',
        default=None,type=int,
        help="maximum number of mismatches to use when "
        "grouping similar barcodes (will be determined "
        "automatically if samplesheet is supplied, otherwise "
        "defaults to 0)")
    reporting.add_argument(
        '--cutoff',action='store',dest='cutoff',
        default=0.001,type=float,
        help="exclude barcodes/barcode groups from reporting "
        "with a smaller fraction of associated reads than "
        "CUTOFF, e.g. '0.01' excludes barcodes with < 1.0%% "
        "of reads (default: 0.001)")
    reporting.add_argument(
        '-s','--sample-sheet',
        action='store',dest='sample_sheet',default=None,
        help="report best matches against barcodes in "
        "SAMPLE_SHEET")
    reporting.add_argument(
        '-r','--report',
        action='store',dest='report_file',default=None,
        help="write report to REPORT_FILE (otherwise write to "
        "stdout)")
    reporting.add_argument(
        '-x','--xls',
        action='store',dest='xls_file',default=None,
        help="write XLS version of report to XLS_FILE")
    reporting.add_argument(
        '-f','--html',
        action='store',dest='html_file',default=None,
        help="write HTML version of report to HTML_FILE")
    reporting.add_argument(
        '-t','--title',
        action='store',dest='title',default=None,
        help="title for HTML report (default: 'Barcodes "
        "Report')")
    reporting.add_argument(
        '-n','--no-report',
        action='store_true',dest='no_report',default=None,
        help="suppress reporting (overrides --report)")
    advanced = p.add_argument_group("Advanced options")
    advanced.add_argument(
        '--sequences',
        action='store_true',dest='count_seqs',default=False,
        help="count sequences instead of barcodes")
    advanced.add_argument(
        '--seq-start',
        action='store',dest='seq_start',type=int,default=None,
        help="specify first base of sequence to analyse (for "
        "--sequences option; default is start at the first "
        "base position)")
    advanced.add_argument(
        '--seq-end',
        action='store',dest='seq_end',type=int,default=None,
        help="specify last base of sequence to analyse (for "
        "--sequences option; default is start at the first "
        "base position)")
    advanced.add_argument(
        '--minimum_read_fraction',
        action='store',dest='minimum_read_fraction',
        metavar='FRACTION',default=0.000001,type=float,
        help="weed out individual barcodes from initial "
        "analysis which have a smaller fraction of reads "
        "than FRACTION, e.g. '0.001' removes barcodes "
        "with < 0.1%% of reads; speeds up analysis at the "
        "expense of accuracy as reported counts will be "
        "approximate (default: 1.0e-6)")
    # Process command line
    args,extra_args = p.parse_known_args()
    # Report name and version
    print("%s %s" % (os.path.basename(sys.argv[0]),__version__))
    # Anything to do?
    if len(extra_args) == 0:
        if args.use_counts:
            p.error("-c: needs at least one barcode counts file")
        else:
            p.error("needs at least one FASTQ file, a bcl2fastq directory")
    # Set default return value
    retval = 0
    # Determine mode
    if args.use_counts:
        # Read counts from counts file(s)
        counts = BarcodeCounter(*extra_args)
    elif len(extra_args) == 1 and os.path.isdir(extra_args[0]):
        # Generate counts from bcl2fastq output
        counts = count_barcodes_bcl2fastq(extra_args[0])
    elif args.count_seqs:
        # Count sequences from fastq files
        counts = count_sequences(extra_args,
                                 start=args.seq_start,
                                 end=args.seq_end)
    else:
        # Generate counts from fastq files
        print("Fastqs supplied on command line")
        fq_groups = group_fastqs_by_name(extra_args)
        fastqs = []
        for grp in fq_groups:
            # Select one Fastq from each group (to
            # avoid multiple counts)
            fq = grp[0]
            if len(grp) > 1:
                print("Keeping %s from group of %d" % (fq,len(grp)))
            fastqs.append(fq)
        counts = count_barcodes(fastqs)
    # Determine subset of lanes to examine
    if args.lanes is not None:
        print("Lanes supplied on command line: %s" % args.lanes)
        lanes = parse_lanes(args.lanes)
    else:
        print("Taking lanes from counts")
        lanes = counts.lanes
    lanes = sorted(list(set(lanes)))
    # Deal with cutoff
    if args.cutoff == 0.0:
        cutoff = None
    else:
        cutoff = args.cutoff
    # Deal with samplesheet
    if args.sample_sheet:
        sample_sheet = os.path.abspath(args.sample_sheet)
    else:
        sample_sheet = None
    # Report settings
    print("Sample sheet: %s" % sample_sheet)
    print("Lanes       : %s" % lanes)
    print("Cutoff      : %s" % cutoff)
    print("Mismatches  : %s" % args.mismatches)
    # Report the counts
    if not args.no_report:
        reporter = Reporter()
        for lane in lanes:
            # Report for each lane
            print("Reporting barcodes for lane %d" % lane)
            if lane not in counts.lanes:
                logging.error("Requested analysis for lane %d but "
                              "only have counts for lanes %s" %
                              (lane,
                               ','.join([str(l) for l in counts.lanes])))
                retval = 1
                continue
            mismatches = args.mismatches
            # Deal with sample sheet if supplied
            use_sample_sheet = sample_sheet
            if sample_sheet:
                with tempfile.NamedTemporaryFile() as fp:
                    # Make a temporary sample sheet with just the
                    # requested lane
                    s = SampleSheet(sample_sheet)
                    if s.has_lanes:
                        use_lanes = (lane,)
                    else:
                        use_lanes = None
                    s = make_custom_sample_sheet(sample_sheet,
                                                 fp.name,
                                                 lanes=use_lanes)
                    if has_10x_indices(fp.name):
                        # Don't pass the sample sheet to the reporter
                        # for lanes with 10x indices
                        use_sample_sheet = None
                        logging.warning("Lane %s has 10xGenomics-style "
                                        "indices in sample sheet; not "
                                        "matching against samplesheet for "
                                        "this lane" % lane)
                    else:
                        # If mismatches not set then determine from
                        # the barcode lengths in the temporary
                        # samplesheet
                        if mismatches is None:
                            barcode_length = None
                            for line in s:
                                index_sequence = \
                                    samplesheet_index_sequence(line)
                                if index_sequence is None:
                                    # Empty barcode sequence in
                                    # samplesheet
                                    length = 0
                                else:
                                    length = len(index_sequence)
                                if barcode_length is None:
                                    barcode_length = length
                                elif length != barcode_length:
                                    logging.error("Lane %s has a mixture "
                                                  "of barcode lengths" %
                                                  lane)
                                    barcode_length = min(barcode_length,
                                                         length)
                                if barcode_length >= 6:
                                    mismatches = 1
                                else:
                                    mismatches = 0
                                # Check for collisions
                                while mismatches and \
                                      check_barcode_collisions(
                                          fp.name,mismatches):
                                    mismatches = mismatches - 1
            # Check mismatches
            if mismatches is None:
                # Set according to barcode lengths found in counts
                barcode_length = min(counts.barcode_lengths(lane=lane))
                if barcode_length >= 6:
                    mismatches = 1
                else:
                    mismatches = 0
            # Report the analysis
            report_barcodes(counts,
                            lane=lane,
                            cutoff=cutoff,
                            sample_sheet=use_sample_sheet,
                            mismatches=mismatches,
                            reporter=reporter,
                            minimum_read_fraction=
                            args.minimum_read_fraction)
        if args.report_file is not None:
            print("Writing report to %s" % args.report_file)
            reporter.write(filen=args.report_file,title=args.title)
        else:
            reporter.write(title=args.title)
        if args.xls_file is not None:
            print("Writing XLS to %s" % args.xls_file)
            reporter.write_xls(args.xls_file,title=args.title)
        if args.html_file is not None:
            print("Writing HTML to %s" % args.html_file)
            reporter.write_html(args.html_file,title=args.title)
    # Output counts if requested
    if args.counts_file_out is not None:
        counts.write(args.counts_file_out)
    # Finish
    sys.exit(retval)
