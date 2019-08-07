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

import sys
import os
import argparse
import logging
import tempfile
import shutil
import time
from bcftbx.TabFile import TabFile
from auto_process_ngs.fastq_utils import pair_fastqs
from auto_process_ngs.icell8.utils import ICell8WellList
from auto_process_ngs.icell8.utils import ICell8Stats
from auto_process_ngs.icell8.utils import get_batch_size
from auto_process_ngs.icell8.utils import batch_fastqs
from auto_process_ngs.icell8.constants import MAXIMUM_BATCH_SIZE

# Module specific logger
logger = logging.getLogger("icell8_stats")

######################################################################
# Main
######################################################################

if __name__ == "__main__":
    print("[%s] ICell8 stats started" % time.strftime("%Y/%m/%d-%H:%M:%S"))
    # Handle the command line
    p = argparse.ArgumentParser()
    p.add_argument("fastqs",nargs='*',metavar="FASTQ_R1 FASTQ_R2",
                   help="FASTQ file pairs")
    p.add_argument("-w","--well-list",
                   dest="well_list_file",default=None,
                   help="iCell8 'well list' file")
    p.add_argument("-u","--unassigned",action='store_true',
                   help="include 'unassigned' reads")
    p.add_argument("-f","--stats-file",
                   dest="stats_file",
                   help="output statistics file")
    p.add_argument("-a","--append",action='store_true',
                   help="append to statistics file")
    p.add_argument("-s","--suffix",
                   dest="suffix",
                   help="suffix to attach to column names")
    p.add_argument("-n","--nprocessors",
                   type=int,default=1,
                   help="number of processors/cores available for "
                   "statistics generation (default: 1)")
    p.add_argument("-m","--max-batch-size",
                   type=int,default=MAXIMUM_BATCH_SIZE,
                   help="maximum number of reads per batch "
                   "when dividing Fastqs (multicore only; "
                   "default: %d)" % MAXIMUM_BATCH_SIZE)
    p.add_argument("-T","--temporary-directory",
                   action="store",default=None,metavar="DIR",
                   help="use DIR for temporaries, not $TMPDIR "
                   "or /tmp")
    args = p.parse_args()

    # Input Fastqs
    fastqs = args.fastqs

    # Well list file
    if args.well_list_file is not None:
        well_list = ICell8WellList(args.well_list_file)
    else:
        well_list = None

    # Number of cores
    nprocs = args.nprocessors
    print("%d processor%s will be used" % (nprocs,
                                           ('s' if nprocs != 1
                                            else '')))

    # Pair up Fastq files
    fastqs,unpaired = pair_fastqs(fastqs)
    if unpaired:
        print("Unpaired Fastqs specified:")
        for fq in unpaired:
            print("- %s" % fq)
        logging.fatal("Unpaired Fastqs specified")
        sys.exit(1)

    # Only need R1 Fastqs
    fastqs = [pair[0] for pair in fastqs]

    # Set up a working directory
    if args.temporary_directory is not None:
        tmpdir = os.path.abspath(args.temporary_directory)
    else:
        try:
            tmpdir = os.path.abspath(os.environ["TMPDIR"])
        except KeyError:
            tmpdir = None
    working_dir = tempfile.mkdtemp(suffix="icell8_stats",
                                   dir=tmpdir)
    print("Using working dir %s" % working_dir)

    # Split into batches for multiprocessing
    if nprocs > 1:
        try:
            batch_size,nbatches = get_batch_size(
                fastqs,
                max_batch_size=args.max_batch_size,
                min_batches=nprocs)
            batched_fastqs = batch_fastqs(
                fastqs,batch_size,
                basename="icell8_stats",
                out_dir=working_dir)
        except Exception as ex:
            logging.critical("Failed to split Fastqs into batches: "
                             "%s" % ex)
            sys.exit(1)
    else:
        batched_fastqs = fastqs

    # Collect statistics
    stats = ICell8Stats(*batched_fastqs,
                        nprocs=nprocs,
                        verbose=True)

    # Remove the working directory
    shutil.rmtree(working_dir)

    # Report the stats
    if args.stats_file is not None:
        # Output column names
        stats_file = os.path.abspath(args.stats_file)
        nreads_col = "Nreads%s" % (''
                                   if args.suffix is None
                                   else args.suffix)
        umis_col = "Distinct_UMIs%s" % (''
                                        if args.suffix
                                        is None else args.suffix)
        if not (os.path.isfile(stats_file) and args.append):
            # Create new stats file
            if well_list is not None:
                # Initialise barcode and sample names from well list
                stats_data = TabFile(column_names=('Barcode',
                                                   'Sample'))
                for barcode in well_list.barcodes():
                    stats_data.append(data=(barcode,
                                            well_list.sample(barcode)))
            else:
                # Barcodes from collected data
                stats_data = TabFile(column_names=('Barcode',))
                for barcode in stats.barcodes():
                    stats_data.append(data=(barcode,))
        else:
            # Append to an existing file
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
                data_line[umis_col] = len(stats.distinct_umis(barcode))
            except KeyError:
                data_line[nreads_col] = 0
                data_line[umis_col] = 0
        # Deal with 'unassigned' reads
        if args.unassigned:
            # Count reads for barcodes not in list
            unassigned_reads = 0
            unassigned_umis = set()
            if well_list is not None:
                expected_barcodes = well_list.barcodes()
            else:
                expected_barcodes = [l['Barcode'] for l in stats_data]
            for barcode in stats.barcodes():
                if barcode not in expected_barcodes:
                    unassigned_reads += stats.nreads(barcode=barcode)
                    unassigned_umis.update(
                        stats.distinct_umis(barcode=barcode))
            # Check if 'unassigned' is already in stats file
            unassigned = stats_data.lookup('Barcode','Unassigned')
            try:
                data_line = unassigned[0]
            except IndexError:
                # Append the line
                data_line = stats_data.append()
                data_line['Barcode'] = 'Unassigned'
            data_line[nreads_col] = unassigned_reads
            data_line[umis_col] = len(unassigned_umis)
        # Write to file
        stats_data.write(filen=stats_file,include_header=True)

    # Report summary
    print("#barcodes     : %s" % len(stats.barcodes()))
    print("#reads        : %s" % stats.nreads())
    print("[%s] ICell8 stats completed" % time.strftime("%Y/%m/%d-%H:%M:%S"))
