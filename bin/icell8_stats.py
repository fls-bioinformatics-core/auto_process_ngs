#!/usr/bin/env python
#
#     icell8_stats.py: collects stats from fastqs from Wafergen iCell8
#     Copyright (C) University of Manchester 2017 Peter Briggs
#
"""
icell8_stats.py

Utility to collect statistics across one or more FASTQ pairs from
Wafergen iCell8.
"""

######################################################################
# Imports
######################################################################

import os
import sys
import argparse
from bcftbx.TabFile import TabFile
from auto_process_ngs.icell8_utils import ICell8WellList
from auto_process_ngs.icell8_utils import ICell8Stats

######################################################################
# Main
######################################################################

if __name__ == "__main__":
    # Handle the command line
    p = argparse.ArgumentParser()
    p.add_argument("FQ_R1",help="R1 FASTQ file")
    p.add_argument("FQ_R2",help="Matching R2 FASTQ file")
    p.add_argument("FQ",nargs='*',help="Additional FASTQ file pairs")
    p.add_argument("-w","--well-list",
                   dest="well_list_file",default=None,
                   help="iCell8 'well list' file")
    p.add_argument("-f","--stats-file",
                   dest="stats_file",
                   help="output statistics file")
    p.add_argument("-a","--append",action='store_true',
                   help="append to statistics file")
    p.add_argument("-s","--suffix",
                   dest="suffix",
                   help="suffix to attach to column names")
    args = p.parse_args()

    # Input Fastqs
    unpaired_fastqs = [args.FQ_R1,args.FQ_R2]
    for fq in args.FQ:
        unpaired_fastqs.append(fq)

    # Well list file
    if args.well_list_file is not None:
        well_list = ICell8WellList(args.well_list_file)
    else:
        well_list = None

    # Gather and report stats
    stats = ICell8Stats(*unpaired_fastqs)

    # Write or append to a file
    if args.stats_file is not None:
        stats_file = os.path.abspath(args.stats_file)
        nreads_col = "Nreads%s" % ('' if args.suffix is None
                                   else args.suffix)
        umis_col = "Unique_UMIs%s" % ('' if args.suffix
                                      is None else args.suffix)
        if not (os.path.isfile(stats_file) and args.append):
            # Create new stats file
            stats_data = TabFile(column_names=('Barcode',
                                               'Sample'))
            for barcode in well_list.barcodes():
                stats_data.append(data=(barcode,
                                        well_list.sample(barcode)))
        else:
            # Append to existing file
            stats_data = TabFile(filen=stats_file,
                                 first_line_is_header=True)
        # Add new columns of data
        stats_data.appendColumn(nreads_col)
        stats_data.appendColumn(umis_col)
        # Populate columns
        for data_line in stats_data:
            barcode = data_line['Barcode']
            try:
                data_line[nreads_col] = stats.nreads(barcode)
                data_line[umis_col] = len(stats.unique_umis(barcode))
            except KeyError:
                pass
        # Write to file
        stats_data.write(filen=stats_file,include_header=True)

    # Report summary
    print "#barcodes   : %s" % len(stats.barcodes())
    print "#reads      : %s" % stats.nreads()
    print "#unique umis: %s" % len(stats.unique_umis())
    print "Finished: %s" % time.ctime()
