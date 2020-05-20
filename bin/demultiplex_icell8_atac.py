#!/usr/bin/env python
#
#     demultiplex_icell8_atac.py: demultiplex reads from ICELL8 scATAC-seq
#     Copyright (C) University of Manchester 2019-2020 Peter Briggs
#
"""
demultiplex_icell8_atac.py

Utility to demultiplex reads from single-cell ATAC-seq Fastqs produced
by Takara's ICELL8 platform.

Demultiplexing can be specified at the level of samples or barcodes.
"""

######################################################################
# Imports
######################################################################

import sys
import os
import tempfile
import shutil
import json
from argparse import ArgumentParser
try:
    # Python 2
    from itertools import izip as zip
except ImportError:
    pass
from multiprocessing import Pool
from bcftbx.simple_xls import XLSWorkBook
from auto_process_ngs import get_version
from auto_process_ngs.applications import Command
from auto_process_ngs.analysis import AnalysisFastq
from auto_process_ngs.icell8.utils import ICell8WellList
from auto_process_ngs.icell8.atac import split_fastq
from auto_process_ngs.icell8.atac import assign_reads
from auto_process_ngs.icell8.atac import concat_fastqs
from auto_process_ngs.icell8.atac import reverse_complement
from auto_process_ngs.icell8.atac import report

######################################################################
# Constants
######################################################################

BATCH_SIZE = 50000000
BUF_SIZE = 1024*16
UNASSIGNED = "Undetermined"

__version__ = get_version()

######################################################################
# Functions
######################################################################

def multiprocessing_map(f,inputs,n=1):
    """
    Wrap the 'Pool.map' functionality from 'multiprocessing'

    Arguments:
      f (function): function to apply to inputs
      inputs (iterable): set of inputs to invoke
        function 'f' with
      n (int): number of concurrent processes to
        use (default: 1)
    """
    pool = Pool(n)
    results = pool.map(f,inputs)
    pool.close()
    pool.join()
    return results

######################################################################
# Main
######################################################################

if __name__ == "__main__":

    # Set up parser
    p = ArgumentParser(description="Assign reads from ICELL8 ATAC "
                       "R1/R2/I1/I2 Fastq set to barcodes and samples "
                       "in a well list file")
    p.add_argument('-v','--version',action='version',
                   version="%(prog)s "+__version__)
    p.add_argument("well_list",metavar="WELL_LIST",help="Well list file")
    p.add_argument("fastq_r1",metavar="FASTQ_R1",help="FASTQ R1")
    p.add_argument("fastq_r2",metavar="FASTQ_R2",help="FASTQ R2")
    p.add_argument("fastq_i1",metavar="FASTQ_I1",help="FASTQ I1")
    p.add_argument("fastq_i2",metavar="FASTQ_I2",help="FASTQ I2")
    p.add_argument("-o","--output-dir",metavar="OUTDIR",
                   dest="output_dir",default="icell8_atac",
                   help="path to demultiplexed output")
    p.add_argument("-b","--batch_size",metavar="N",
                   dest="batch_size",type=int,default=BATCH_SIZE,
                   help="batch size for splitting index read Fastqs "
                   "(default: %d)" % BATCH_SIZE)
    p.add_argument("-n","--nprocessors",metavar="N",
                   dest="nprocs",type=int,default=1,
                   help="number of processors to use")
    p.add_argument("-m","--mode",action="store",
                   choices=['samples','barcodes'],
                   dest="mode",default="samples",
                   help="demultiplex reads by sample (default) "
                   "or by barcode")
    p.add_argument("--swap-i1-i2",action='store_true',
                   dest="swap_i1_and_i2",
                   help="swap supplied I1 and I2 Fastqs")
    p.add_argument("--reverse-complement",action="store",
                   choices=['i1','i2','both'],
                   dest="reverse_complement",
                   help="reverse complement one or both of the "
                   "indices from the well list file")
    p.add_argument("-u","--update-read-headers",action='store_true',
                   dest="update_read_headers",default=False,
                   help="update read headers in the output Fastqs "
                   "to include the matching index sequence (i.e. "
                   "barcode) from the well list file")
    p.add_argument("--no-demultiplexing",action='store_true',
                   dest="no_demultiplexing",
                   help="don't generate demultiplexed Fastqs "
                   "(only the stats)")
    args = p.parse_args()

    # Report name and version
    print("%s version %s" % (os.path.basename(sys.argv[0]),__version__))

    # Sort out the supplied Fastqs
    fastqs = (args.fastq_r1,
              args.fastq_r2,
              args.fastq_i1,
              args.fastq_i2)
    fastq_r1 = None
    fastq_r2 = None
    fastq_i1 = None
    fastq_i2 = None
    for fastq in fastqs:
        fastq = os.path.abspath(fastq)
        fq = AnalysisFastq(fastq)
        if fq.is_index_read:
            if fq.read_number == 1:
                fastq_i1 = fastq
            else:
                fastq_i2 = fastq
        else:
            if fq.read_number == 1:
                fastq_r1 = fastq
            else:
                fastq_r2 = fastq
    # Report assignments
    report("Assigned Fastqs:")
    report("-- R1: %s" % fastq_r1)
    report("-- R2: %s" % fastq_r2)
    report("-- I1: %s" % fastq_i1)
    report("-- I2: %s" % fastq_i2)

    # Well list file
    well_list_file = os.path.abspath(args.well_list)
    report("Well list: %s" % well_list_file)

    # Report settings
    report("Using %d processors" % args.nprocs)
    if args.reverse_complement:
        report("%s barcode%s from well list will be reverse complemented when "
               "matching to Fastqs" % (("Both I1 and I2"
                                        if args.reverse_complement == 'both'
                                        else args.reverse_complement.upper()),
                                       ('s'
                                        if args.reverse_complement == 'both'
                                        else '')))
    if args.swap_i1_and_i2:
        report("I1 and I2 barcode components will be swapped when matching to "
               "Fastqs")

    # Set up output directory
    output_dir = os.path.abspath(args.output_dir)
    if os.path.exists(output_dir):
        report("Output directory '%s' already exists" % output_dir,
               fp=sys.stderr)
        sys.exit(1)
    os.mkdir(output_dir)
    report("Created output directory: %s" % output_dir)

    # Set up temporary directory
    tmp_dir = tempfile.mkdtemp(prefix="__demultiplex.",dir=os.getcwd())
    report("Created temporary working directory: %s" % tmp_dir)

    # Load the data from the well list
    report("Reading data from %s" % well_list_file)
    well_list = ICell8WellList(well_list_file)
    report("Number of barcodes (cells): %d" % len(well_list.barcodes()))
    report("Number of samples         : %d" % len(well_list.samples()))

    # Split fastqs into batches
    batched_fastqs_dir = os.path.join(tmp_dir,"batched_fastqs")
    os.mkdir(batched_fastqs_dir)
    fastqs = []
    inputs = list()
    for fastq in (fastq_r1,fastq_r2,fastq_i1,fastq_i2):
        inputs.append((fastq,
                       args.batch_size,
                       batched_fastqs_dir,))
    results = multiprocessing_map(split_fastq,inputs,
                                  n=args.nprocs)
    for result in results:
        fastqs.append(result)

    # Build barcode list
    report("Assigning reads to barcodes and samples")
    inputs = list()
    for fastq_r1,fastq_r2,fastq_i1,fastq_i2 in zip(fastqs[0],
                                                   fastqs[1],
                                                   fastqs[2],
                                                   fastqs[3]):
        inputs.append((fastq_r1,fastq_r2,
                       fastq_i1,fastq_i2,
                       well_list_file,
                       args.mode,
                       args.swap_i1_and_i2,
                       args.reverse_complement,
                       args.update_read_headers,
                       tmp_dir,
                       unassigned,))
    results = multiprocessing_map(assign_reads,inputs,
                                  n=args.nprocs)
    report("Collecting outputs from batches")
    batches = list()
    barcode_counts = {}
    sample_counts = {}
    undetermined_barcode_files = list()
    for batch_id,counts,undetermined in results:
        report("Handling batch %s" % batch_id)
        batches.append(batch_id)
        for barcode in counts:
            try:
                barcode_counts[barcode] += counts[barcode]
            except KeyError:
                barcode_counts[barcode] = counts[barcode]
            if barcode == UNASSIGNED:
                sample = UNASSIGNED
            else:
                sample = well_list.sample(barcode)
            try:
                sample_counts[sample] += counts[barcode]
            except KeyError:
                sample_counts[sample] = counts[barcode]
        undetermined_barcode_files.append(undetermined)
    batches = sorted(batches)

    # Report number of reads assigned to each (known) barcode
    barcode_counts_file = os.path.join(tmp_dir,"barcode_counts.txt")
    with open(barcode_counts_file,'w') as fp:
        fp.write("#Sample\tWell list barcode\tFastq barcode\tNreads\n")
        for barcode in well_list.barcodes():
            lineout = list()
            # Sample
            sample = well_list.sample(barcode)
            lineout.append(sample)
            # Original (well list) barcode
            lineout.append(barcode)
            # Equivalent 'raw' barcode from the Fastqs
            i1,i2 = barcode.split('+')
            if args.reverse_complement:
                if args.reverse_complement in ('i1','both'):
                    i1 = reverse_complement(i1)
                if args.reverse_complement in ('i2','both'):
                    i2 = reverse_complement(i2)
            if args.swap_i1_and_i2:
                i1,i2 = i2,i1
            fastq_barcode = "%s+%s" % (i1,i2)
            lineout.append(fastq_barcode)
            # Count
            try:
                count = barcode_counts[barcode]
            except KeyError:
                count = 0
            lineout.append(count)
            # Output to file
            fp.write("%s\n" % '\t'.join([str(x) for x in lineout]))
    report("Counts for assigned barcodes written to %s" % barcode_counts_file)

    # Report number of reads assigned to each sample
    samples = well_list.samples()
    samples.insert(0,UNASSIGNED)
    sample_counts_file = os.path.join(tmp_dir,"sample_counts.txt")
    report("Number of reads assigned to each sample:")
    with open(sample_counts_file,'w') as fp:
        fp.write("#Sample\tNreads\n")
        for sample in samples:
            count = sample_counts[sample]
            report("-- %s: %d" % (sample,count))
            fp.write("%s\t%s\n" % (sample,count))
    report("Sample counts written to %s" % sample_counts_file)

    # Report unassigned barcodes
    report("Sorting and counting unassigned barcodes")
    cmd = Command('sort')
    cmd.add_args(*undetermined_barcode_files)
    cmd.add_args('|',
                 'uniq','-c',
                 '|',
                 'sort','-k',1,'-n','-r')
    # Run the command via a script
    script_file = os.path.join(tmp_dir,"count_undetermined.sh")
    cmd.make_wrapper_script(shell='/bin/bash',
                            filen=script_file,
                            quote_spaces=True)
    log_file = os.path.join(tmp_dir,"count_undetermined.log")
    status = Command("/bin/bash",script_file).run_subprocess(
        log=log_file,
        working_dir=tmp_dir)
    os.remove(script_file)
    # Check return status
    if status != 0:
        report("Failed to sort and count undetermined barcodes",
               fp=sys.stderr)
    # Rewrite to a file
    undetermined_counts_file = os.path.join(tmp_dir,
                                            "undetermined_counts.txt")
    with open(undetermined_counts_file,'w') as fp:
        fp.write("#Barcode\tNreads\n")
        with open(log_file,'r') as fpp:
            for line in fpp:
                count,barcode = line.strip().split()
                fp.write("%s\t%s\n" % (barcode,count))
    os.remove(log_file)
    report("Undetermined barcode counts written to %s" %
           undetermined_counts_file)

    # Generate a JSON file with key data
    json_data = dict()
    json_data['summary'] = dict()
    json_data['summary']['command_line'] = ' '.join(sys.argv)
    json_data['summary']['command_line_arguments'] = sys.argv[1:]
    json_data['summary']['swap_i1_and_i2'] = args.swap_i1_and_i2
    json_data['summary']['reverse_complement'] = args.reverse_complement
    json_data['summary']['well_list_file'] = well_list_file
    json_data['summary']['number_of_barcodes'] = len(well_list.barcodes())
    with open(barcode_counts_file,'r') as fp:
        number_of_barcodes_with_reads = 0
        for line in fp:
            if line.startswith('#'):
                continue
            sample,barcode,fastq_barcode,count = line.strip().split('\t')
            count = int(count)
            if count > 0:
                number_of_barcodes_with_reads += 1
        json_data['summary']['number_of_barcodes_with_reads'] = \
                                        number_of_barcodes_with_reads
    json_data['summary']['number_of_samples'] = len(well_list.samples())
    json_data['reads_per_sample'] = dict()
    for sample in well_list.samples():
        json_data['reads_per_sample'][sample] = sample_counts[sample]
    json_data['reads_per_sample'][UNASSIGNED] = sample_counts[UNASSIGNED]
    json_data['reads_per_barcode'] = dict()
    with open(barcode_counts_file,'r') as fp:
        number_of_barcodes_with_reads = 0
        for line in fp:
            if line.startswith('#'):
                continue
            sample,barcode,fastq_barcode,count = line.strip().split('\t')
            count = int(count)
            json_data['reads_per_barcode'][barcode] = {
                'sample': sample,
                'barcode': barcode,
                'fastq_barcode': fastq_barcode,
                'assigned_reads': count,
            }
    json_data['undetermined_barcodes'] = dict()
    json_data['undetermined_barcodes']['barcodes'] = dict()
    with open(undetermined_counts_file,'r') as fp:
        n_reported = 0
        for line in fp:
            if line.startswith('#'):
                continue
            barcode,count=line.rstrip('\n').split('\t')
            json_data['undetermined_barcodes']['barcodes'][barcode] = int(count)
            n_reported += 1
            if n_reported == 100:
                break
        json_data['undetermined_barcodes']['number_reported'] = n_reported

    # Write data to JSON file
    json_file = os.path.join(output_dir,"icell8_atac_stats.json")
    with open(json_file,'w') as fp:
        json.dump(json_data,fp,indent=2)
    report("Wrote JSON file: %s" % json_file)

    # Write stats to XLSX file
    xlsx_stats_file = os.path.join(output_dir,"icell8_atac_stats.xlsx")
    wb = XLSWorkBook("ICELL8 scATAC-seq")
    # Summary
    ws = wb.add_work_sheet("summary","Summary")
    data = json_data['summary']
    ws.insert_block_data("Command line: %s" % data['command_line'])
    ws.append_row(["Swap I1/I2 fastqs",
                   "%s" % ('YES' if data['swap_i1_and_i2'] else 'NO',)])
    ws.append_row(["Reverse complement",
                   "%s" % (data['reverse_complement'].upper()
                           if data['reverse_complement'] else '',)])
    ws.append_row(["Number of cells",data['number_of_barcodes']])
    ws.append_row(["Number of cells with reads",
                   data['number_of_barcodes_with_reads']])
    # Reads per sample
    data = json_data['reads_per_sample']
    ws = wb.add_work_sheet("samples","Reads per sample")
    ws.append_row(["Sample","Nreads"])
    for sample in data:
        ws.append_row([sample,data[sample]])
    ws.freeze_panes = 'A2'
    # Reads per barcode
    data = json_data['reads_per_barcode']
    ws = wb.add_work_sheet("barcodes","Reads per barcode")
    ws.append_row(["Sample","Well list barcode","Fastq barcode","Nreads"])
    for barcode in data:
        ws.append_row([data[barcode]['sample'],
                       data[barcode]['barcode'],
                       data[barcode]['fastq_barcode'],
                       data[barcode]['assigned_reads']])
    ws.freeze_panes = 'A2'
    # Undetermined barcodes
    data = json_data['undetermined_barcodes']
    ws = wb.add_work_sheet("undetermined","Undetermined barcodes")
    ws.append_row(["Top %d unassigned barcodes" % data['number_reported'],])
    ws.append_row(["Barcode","Nreads"])
    for barcode in data['barcodes']:
        ws.append_row([barcode,
                       data['barcodes'][barcode]])
    ws.freeze_panes = 'A3'
    wb.save_as_xlsx(xlsx_stats_file)
    report("Wrote XLSX file: %s" % xlsx_stats_file)

    # If no demultiplexing was requested then clean up
    # and finish
    if args.no_demultiplexing:
        report("Removing %s" % tmp_dir)
        shutil.rmtree(tmp_dir)
        report("Finished")
        sys.exit()

    # Concatenate and compress Fastqs for each sample
    report("Concatenating and compressing Fastqs")
    inputs = list()
    if args.mode == 'samples':
        # Concatenate across samples
        for index,sample in enumerate(samples):
            for read in ('R1','R2','I1','I2'):
                inputs.append((sample,
                               index,
                               None,
                               read,
                               batches,
                               tmp_dir,
                               output_dir))
    elif args.mode == 'barcodes':
        # Concatenate across barcodes
        for index,sample in enumerate(samples):
            for barcode in well_list.barcodes():
                if well_list.sample(barcode) == sample:
                    for read in ('R1','R2','I1','I2'):
                        inputs.append((sample,
                                       index,
                                       barcode,
                                       read,
                                       batches,
                                       tmp_dir,
                                       output_dir))
    results = multiprocessing_map(concat_fastqs,inputs,
                                  n=args.nprocs)

    # Done
    report("Removing %s" % tmp_dir)
    shutil.rmtree(tmp_dir)
    report("Finished")
