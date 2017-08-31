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
from bcftbx import FASTQFile
from bcftbx.TabFile import TabFile
from auto_process_ngs.fastq_utils import pair_fastqs
from auto_process_ngs.fastq_utils import get_read_number
from auto_process_ngs.applications import Command
from auto_process_ngs.icell8_utils import ICell8WellList
from auto_process_ngs.icell8_utils import ICell8Stats

######################################################################
# Functions
######################################################################

def batch_fastqs(fastqs,nbatches,basename="batched",
                 out_dir=None):
    """
    Splits reads from one or more Fastqs into batches

    Concatenates input Fastq files and then splits
    reads into smaller Fastqs using the external 'batch'
    utility.

    Arguments:
      fastqs (list): list of paths to one or more Fastq
        files to take reads from
      nbatches (int): number of batches to output reads
        into
      basename (str): optional basename to use for the
        output Fastq files (default: 'batched')
      out_dir (str): optional path to a directory where
        the batched Fastqs will be written
    """
    # Count the total number of reads
    print "Fetching read counts:"
    nreads = 0
    for fq in fastqs:
        n = FASTQFile.nreads(fq)
        print "%s:\t%d" % (os.path.basename(fq),n)
        nreads += n
    print "Total reads: %d" % nreads

    # Determine batch size
    batch_size = nreads/nbatches
    if nreads%batch_size:
        # Round up batch size
        batch_size += 1
    assert(batch_size*nbatches >= nreads)
    print "Creating batches of %d reads" % batch_size

    # Check if fastqs are compressed
    gzipped = fastqs[0].endswith('.gz')
    if gzipped:
        batch_cmd = Command('zcat')
    else:
        batch_cmd = Command('cat')

    # Get the read number
    read_number = get_read_number(fastqs[0])
    suffix = ".r%s.fastq" % read_number

    # Build and run the batching command
    batch_cmd.add_args(*fastqs)
    batch_cmd.add_args('|',
                       'split',
                       '-l',batch_size*4,
                       '-d',
                       '-a',3,
                       '--additional-suffix=%s' % suffix,
                       '-',
                       os.path.join(out_dir,"%s.B" % basename))
    batch_script = os.path.join(out_dir,"batch.sh")
    batch_cmd.make_wrapper_script("/bin/bash",
                                  batch_script)

    # Check for successful exit code
    retcode = Command("/bin/bash",
                      batch_script).run_subprocess(
                          working_dir=out_dir)
    if retcode != 0:
        raise Exception("Batching failed: exit code %s" % retcode)

    # Collect and return the batched Fastq names
    batched_fastqs = [os.path.join(out_dir,
                                   "%s.B%03d%s"
                                   % (basename,i,suffix))
                      for i in xrange(0,nbatches)]
    return batched_fastqs

######################################################################
# Main
######################################################################

if __name__ == "__main__":
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
    print "%d processors will be used" % nprocs

    # Pair up Fastq files
    fastqs,unpaired = pair_fastqs(fastqs)
    if unpaired:
        print "Unpaired Fastqs specified:"
        for fq in unpaired:
            print "- %s" % fq
        logging.fatal("Unpaired Fastqs specified")
        sys.exit(1)

    # Only need R1 Fastqs
    fastqs = [pair[0] for pair in fastqs]

    # Set up a working directory
    working_dir = tempfile.mkdtemp(suffix="icell8_stats")
    print "Using working dir %s" % working_dir

    # Split into batches
    try:
        batched_fastqs = batch_fastqs(fastqs,nprocs,
                                      basename="icell8_stats",
                                      out_dir=working_dir)
    except Exception as ex:
        logging.critical("Failed to split Fastqs into batches: "
                         "%s" % ex)
        sys.exit(1)

    # Collect statistics
    stats = ICell8Stats(*batched_fastqs,nprocs=nprocs)

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
    print "#barcodes     : %s" % len(stats.barcodes())
    print "#reads        : %s" % stats.nreads()
