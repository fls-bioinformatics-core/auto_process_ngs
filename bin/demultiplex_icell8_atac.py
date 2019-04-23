#!/usr/bin/env python2
#
#     demultiplex_icell8_atac.py: demultiplex reads from ICELL8 scATAC-seq
#     Copyright (C) University of Manchester 2019 Peter Briggs
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
import time
import tempfile
import glob
import shutil
from argparse import ArgumentParser
from itertools import izip
from collections import defaultdict
from multiprocessing import Pool
from bcftbx.FASTQFile import FastqIterator
from bcftbx.simple_xls import XLSWorkBook
from bcftbx.ngsutils import getreads
from auto_process_ngs.applications import Command
from auto_process_ngs.analysis import AnalysisFastq
from auto_process_ngs.icell8.utils import ICell8WellList
from auto_process_ngs.utils import ProgressChecker
from auto_process_ngs.utils import BufferedOutputFiles

######################################################################
# Constants
######################################################################

BATCH_SIZE = 50000000
BUF_SIZE = 1024*16
UNASSIGNED = "Undetermined"
REVERSE_COMPLEMENTS =  { 'A':'T',
                         'T':'A',
                         'G':'C',
                         'C':'G',
                         'N':'N' }

######################################################################
# Functions
######################################################################

def reverse_complement(s):
    """
    Return reverse complement of a sequence

    Arguments:
      s (str): sequence to be reverse complemented

    Returns:
      String: reverse complement of input sequence
    """
    s1 = ''
    for c in s[::-1]:
        s1 += REVERSE_COMPLEMENTS[c]
    return s1

def report(msg,fp=sys.stdout):
    """
    Write timestamped message

    Arguments:
      msg (string): text to be reported
      fp (file): stream to report to (defaults
        to stdout)
    """
    fp.write("[%s] %s\n" % (time.strftime("%Y-%m-%d %H:%M:%S"),
                            msg))

def split_fastq(args):
    """
    Split Fastq into batches

    Intended to be invoked via 'map' or similar function

    Arguments are supplied in a single list which should
    contain the following items:

    - Fastq: path to Fastq file to split
    - batch_size: size of each batch
    - working_dir: working directory to write batches to

    Arguments:
      args (list): list containing the arguments supplied to
        the splitter

    Returns:
      List: list of batched Fastqs
    """
    # Unpack arguments
    fastq,batch_size,working_dir = args
    fq = AnalysisFastq(fastq)
    label = fq.basename
    # Build split command
    # Use [z]cat FASTQ | split -l BATCH_SIZE*4 -d -a 3 \
    # --additional-suffix=.fastq - BASENAME_
    if batch_size:
        report("[%s] Splitting %s into batches of %d reads" %
               (label,fastq,batch_size))
        if fq.extension.endswith(".gz"):
            cat = "zcat"
        else:
            cat = "cat"
        basename = "%s_B" % fq.basename
        cmd = Command(cat,fastq,
                      '|',
                      'split',
                      '-l',batch_size*4,
                      '-d',
                      '-a',3,
                      '--additional-suffix=.fastq',
                      '-',
                      os.path.join(working_dir,basename))
        # Run the command via a script
        script_file = os.path.join(working_dir,"batch.%s.sh" % fq.basename)
        cmd.make_wrapper_script(shell='/bin/bash',
                                filen=script_file,
                                quote_spaces=True)
        status = Command("/bin/bash",script_file).run_subprocess(
            working_dir=working_dir)
        os.remove(script_file)
        # Check return status
        if status != 0:
            report("[%s] Failed to split %s" % (label,fastq),
                   fp=sys.stderr)
            sys.exit(status)
        # Gather Fastqs
        fastqs = sorted(
            glob.glob(
                os.path.join(working_dir,"%s*.fastq" % basename)))
        report("[%s] Created %d batches" % (label,len(fastqs)))
    else:
        # No batching - just decompress
        report("[%s] Uncompressing/copying %s" % (fastq,label))
        if fq.extension.endswith(".gz"):
            cat = "zcat"
        else:
            cat = "cat"
        cmd = Command(cat,fastq)
        status = cmd.run_subprocess(
            working_dir=working_dir,
            log=os.path.join(tmp_dir,"%s.fastq" % fq.basename))
        # Check return status
        if status != 0:
            report("Failed to uncompress/copy %s" % fastq,
                   fp=sys.stderr)
            sys.exit(status)
        fastqs = ["%s.fastq" % fq.basename]
    return fastqs

def assign_reads(args):
    """
    Assign reads to samples from batched ICELL8 ATAC Fastqs

    Intended to be invoked via 'map' or similar function

    Arguments are supplied in a single list which should
    contain the following items:

    - R1 Fastq: path to R1 Fastq file
    - R2 Fastq: path to R2 Fastq file
    - I1 Fastq: path to I1 Fastq file
    - I2 Fastq: path to I2 Fastq file
    - well list: path to the well list file
    - mode: either 'samples' or 'barcodes'
    - swap_i1_and_i2: boolean indicating whether I1 and I2
      Fastqs should be swapped for matching
    - reverse_complement: either None, 'i1', 'i2' or both
    - working_dir: working directory to write batches to
    - unassigned: basename for output files

    In 'samples' mode assignment is done to samples only;
    in 'barcodes' mode assignment is done to samples and
    barcodes.

    Arguments:
      args (list): list containing the arguments supplied to
        the read assigner

    Returns:
      Tuple: tuple consisting of (batch id,barcode_counts,
        undetermined_barcodes_file).
    """
    # Unpack arguments
    fastq_r1,fastq_r2,fastq_i1,fastq_i2,well_list_file,mode,swap_i1_and_i2,reverse_complement_index,working_dir,unassigned = args
    # Batch ID is the trailing part of the name
    batch_id = AnalysisFastq(fastq_i1).extras.strip('_')
    # Label is sample name plus batch name
    label = "%s/%s" % (AnalysisFastq(fastq_i1).sample_name,batch_id)
    report("[%s] Assigning reads from R1/R2 Fastq pairs based on I1/I2 Fastqs:"
           % label)
    report("[%s] -- R1: %s" % (label,os.path.basename(fastq_r1)))
    report("[%s] -- R2: %s" % (label,os.path.basename(fastq_r2)))
    report("[%s] -- I1: %s" % (label,os.path.basename(fastq_i1)))
    report("[%s] -- I2: %s" % (label,os.path.basename(fastq_i2)))
    report("[%s] -- Well list: %s" % (label,os.path.basename(well_list_file)))
    report("[%s] Mode is '%s'" % (label,mode))
    if swap_i1_and_i2:
        report("[%s] Swapping I1 and I2 Fastqs for matching to well list" %
               label)
    # Check mode
    if mode not in ("samples","barcodes"):
        report("[%s] Unrecognised mode!" % label,fp=sys.stderr)
    # Working directory
    if working_dir is None:
        working_dir = os.getcwd()
    os.mkdir(os.path.join(working_dir,batch_id))
    # Read well list file to get barcodes and lookups
    well_list = ICell8WellList(well_list_file)
    sample_lookup = defaultdict(lambda: unassigned)
    barcode_lookup = defaultdict(lambda: unassigned)
    for sample in well_list.samples():
        barcode_lookup[sample] = list()
    barcodes = well_list.barcodes()
    for barcode in barcodes:
        sample = well_list.sample(barcode)
        sample_lookup[barcode] = sample
        barcode_lookup[sample].append(barcode)
    # Generate adjusted versions of barcodes for matching
    # against barcodes derived from Fastqs
    fastq_barcode_lookup = defaultdict(lambda: None)
    for barcode in barcodes:
        i1,i2 = barcode.split('+')
        if reverse_complement_index:
            if reverse_complement_index in ('i1','both'):
                # Reverse complement the I1 part of each barcode
                i1 = reverse_complement(i1)
            if reverse_complement_index in ('i2','both'):
                # Reverse complement the I2 part of each barcode
                i2 = reverse_complement(i2)
        if swap_i1_and_i2:
            i2,i1 = i1,i2
        fastq_barcode_lookup["%s+%s" % (i1,i2)] = barcode
    # File to write undetermined barcodes to
    undetermined_barcodes_file = os.path.join(working_dir,
                                              batch_id,
                                              "undetermined_barcodes.txt")
    # Set up output files for samples
    samples = well_list.samples()
    samples.insert(0,unassigned)
    fpp = BufferedOutputFiles()
    for read in ('R1','R2','I1','I2'):
        for index,sample in enumerate(samples):
            if mode == 'samples':
                # Output files will only have sample names
                name = "%s_%s" % (sample,read)
                filen = "%s_S%d_%s_001.fastq" % (sample,index,read)
                fpp.open(name,
                         os.path.join(working_dir,batch_id,filen))
            elif mode == 'barcodes':
                # Output files will have sample name plus barcode
                if sample != unassigned:
                    # Standard samples
                    for barcode in barcode_lookup[sample]:
                        name = "%s_%s_%s" % (sample,barcode,read)
                        filen = "%s_S%d_%s_%s_001.fastq" % \
                                (sample,index,barcode,read)
                        fpp.open(name,
                                 os.path.join(working_dir,batch_id,filen))
                else:
                    # Unassigned reads
                    name = "%s_%s" % (sample,read)
                    filen = "%s_S%d_%s_001.fastq" % (sample,index,read)
                    fpp.open(name,
                             os.path.join(working_dir,batch_id,filen))
    barcode_counts = { unassigned: 0, }
    for barcode in well_list.barcodes():
        barcode_counts[barcode] = 0
    # Examine indices and assign reads
    ii = 0
    progress = ProgressChecker(every=1000000)
    if mode == 'samples':
        # Assigning reads to samples
        with open(undetermined_barcodes_file,"w") as fp:
            for r1,r2,i1,i2 in izip(getreads(fastq_r1),
                                    getreads(fastq_r2),
                                    getreads(fastq_i1),
                                    getreads(fastq_i2)):
                # Get barcodes to match against adjusted
                # versions from well list
                fastq_barcode = "%s+%s" % (i1[1],i2[1])
                # Get "real" barcode
                barcode = fastq_barcode_lookup[fastq_barcode]
                # Add to counts
                try:
                    barcode_counts[barcode] += 1
                except KeyError:
                    barcode_counts[unassigned] += 1
                # Determine sample
                sample = sample_lookup[barcode]
                # Write the reads to the appropriate destinations
                fpp.write("%s_R1" % sample,'\n'.join(r1))
                fpp.write("%s_R2" % sample,'\n'.join(r2))
                fpp.write("%s_I1" % sample,'\n'.join(i1))
                fpp.write("%s_I2" % sample,'\n'.join(i2))
                # Write Fastq version of unassigned barcode to file
                if sample == unassigned:
                    fp.write("%s\n" % fastq_barcode)
                # Report progress
                ii += 1
                if progress.check(ii):
                    report("[%s]...%d reads examined" % (label,ii))
    elif mode == 'barcodes':
        # Assigning reads to barcodes
        with open(undetermined_barcodes_file,"w") as fp:
            for r1,r2,i1,i2 in izip(getreads(fastq_r1),
                                    getreads(fastq_r2),
                                    getreads(fastq_i1),
                                    getreads(fastq_i2)):
                # Get barcodes to match against adjusted
                # versions from well list
                fastq_barcode = "%s+%s" % (i1[1],i2[1])
                # Get "real" barcode
                barcode = fastq_barcode_lookup[fastq_barcode]
                # Add to counts
                try:
                    barcode_counts[barcode] += 1
                except KeyError:
                    barcode_counts[unassigned] += 1
                # Determine sample
                sample = sample_lookup[barcode]
                # Write the reads to the appropriate destinations
                if sample != unassigned:
                    # Assign to sample and barcode
                    fpp.write("%s_%s_R1" % (sample,barcode),'\n'.join(r1))
                    fpp.write("%s_%s_R2" % (sample,barcode),'\n'.join(r2))
                    fpp.write("%s_%s_I1" % (sample,barcode),'\n'.join(i1))
                    fpp.write("%s_%s_I2" % (sample,barcode),'\n'.join(i2))
                else:
                    # Write unassigned barcode to file
                    fpp.write("%s_R1" % sample,'\n'.join(r1))
                    fpp.write("%s_R2" % sample,'\n'.join(r2))
                    fpp.write("%s_I1" % sample,'\n'.join(i1))
                    fpp.write("%s_I2" % sample,'\n'.join(i2))
                    # Write Fastq version of unassigned barcode to file
                    fp.write("%s\n" % fastq_barcode)
                # Report progress
                ii += 1
                if progress.check(ii):
                    report("[%s]...%d reads examined" % (label,ii))
    report("[%s] Finished processing batch %s" % (label,batch_id))
    # Close files
    fpp.close()
    # Remove original files
    for fq in (fastq_r1,fastq_r2,fastq_i1,fastq_i2):
        report("[%s] Removing %s" % (label,fq))
        os.remove(fq)
    # Returns tuple with batch ID, barcode counts and
    # file with list of undetermined barcodes
    return (batch_id,barcode_counts,undetermined_barcodes_file)

def concat_fastqs(args):
    """
    Concatenate Fastqs for a sample across batches

    Intended to be invoked via 'map' or similar function

    Arguments are supplied in a single list which should
    contain the following items:

    - sample: name of sample to concatenate Fastqs for
    - index: integer index to assign to the sample in
      output file name
    - barcode: (optional) barcode to concatenate Fastqs
      for (set to None when concatenating across samples)
    - read: read identifier e.g. 'R1' or 'I2'
    - batches: list of batch IDs to concatenate across
    - working_dir: working directory where batches are
      located
    - final_dir: directory to write concatenated Fastq to

    Arguments:
      args (list): list containing the arguments supplied to
        the read assigner

    Returns:
      String: path of concatenated Fastq.
    """
    # Unpack arguments
    sample,index,barcode,read,batches,working_dir,final_dir = args
    label = "%s:%s" % (sample,read)
    if barcode:
        report("[%s] Concatenating batched %s Fastqs for sample '%s' (%s)" %
               (label,read,sample,barcode))
    else:
        report("[%s] Concatenating batched %s Fastqs for sample '%s'" %
               (label,read,sample))
    # Name of output file
    fastq = os.path.join(final_dir,
                         "%s_S%d_%s%s_001.fastq.gz" % 
                         (sample,
                          index,
                          "%s_" % barcode if barcode else "",
                          read))
    # Collect the Fastqs
    fastqs = list()
    for batch_id in batches:
        fq = os.path.join(working_dir,
                          batch_id,
                          "%s_S%d_%s%s_001.fastq" %
                          (sample,
                           index,
                           "%s_" % barcode if barcode else "",
                           read))
        if os.path.exists(fq):
            fastqs.append(fq)
    if not fastqs:
        # Write an empty output file
        report("[%s] No Fastqs found" % label)
    else:
        # Concatenate and pipe into gzip
        cmd = Command('cat')
        cmd.add_args(*fastqs)
        cmd.add_args('|','gzip','-')
        # Run the command via a script
        script_file = os.path.join(working_dir,
                                   "concat_gzip_%s%s_%s.sh" %
                                   (sample,
                                    "%s_" % barcode if barcode else "",
                                    read))
        cmd.make_wrapper_script(shell='/bin/bash',
                                filen=script_file,
                                quote_spaces=True)
        status = Command("/bin/bash",
                         script_file).run_subprocess(
                             log="%s.part" % fastq,
                             working_dir=working_dir)
        os.remove(script_file)
        # Check return status
        if status != 0:
            report("[%s] Failed to concatenate and compress batches" % label,
                   fp=sys.stderr)
        else:
            # Move Fastq to final location
            os.rename("%s.part" % fastq,fastq)
            report("[%s] Wrote %s" % (label,fastq))
    # Remove original files
    for fq in fastqs:
        report("[%s] Removing %s" % (label,fq))
        os.remove(fq)
    # Return the Fastq file name
    return fastq

######################################################################
# Main
######################################################################

if __name__ == "__main__":

    # Set up parser
    p = ArgumentParser()
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
                   help="batch size for splitting index read Fastqs")
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
    p.add_argument("--no-demultiplexing",action='store_true',
                   dest="no_demultiplexing",
                   help="don't generate demultiplexed Fastqs "
                   "(only the stats)")
    args = p.parse_args()

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
    if args.nprocs > 1:
        pool = Pool(args.nprocs)
        results = pool.map(split_fastq,inputs)
        pool.close()
        pool.join()
    else:
        results = map(split_fastq,inputs)
    for result in results:
        fastqs.append(result)

    # Build barcode list
    report("Assigning reads to barcodes and samples")
    inputs = list()
    for fastq_r1,fastq_r2,fastq_i1,fastq_i2 in izip(fastqs[0],
                                                    fastqs[1],
                                                    fastqs[2],
                                                    fastqs[3]):
        inputs.append((fastq_r1,fastq_r2,
                       fastq_i1,fastq_i2,
                       well_list_file,
                       args.mode,
                       args.swap_i1_and_i2,
                       args.reverse_complement,
                       tmp_dir,
                       UNASSIGNED,))
    if args.nprocs > 1:
        pool = Pool(args.nprocs)
        results = pool.map(assign_reads,inputs)
        pool.close()
        pool.join()
    else:
        results = map(assign_reads,inputs)
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

    # Write stats to XLSX file
    report("Making XLSX file with statistics")
    xlsx_stats_file = os.path.join(output_dir,"icell8_atac_stats.xlsx")
    wb = XLSWorkBook("ICELL8 scATAC-seq")
    # Summary
    ws = wb.add_work_sheet("summary","Summary")
    ws.insert_block_data("Command line: %s" % ' '.join(sys.argv[1:]))
    ws.append_row(["Swap I1/I2 fastqs",
                   "%s" % ('YES' if args.swap_i1_and_i2 else 'NO',)])
    ws.append_row(["Reverse complement",
                   "%s" % (args.reverse_complement.upper()
                           if args.reverse_complement else '',)])
    # Reads per sample
    ws = wb.add_work_sheet("samples","Reads per sample")
    with open(sample_counts_file,'r') as fp:
        ws.insert_block_data(fp.read())
    # Reads per barcode
    ws = wb.add_work_sheet("barcodes","Reads per barcode")
    with open(barcode_counts_file,'r') as fp:
        ws.insert_block_data(fp.read())
    # Undetermined barcodes
    ws = wb.add_work_sheet("undetermined","Undetermined barcodes")
    with open(undetermined_counts_file,'r') as fp:
        ws.append_row(data=("Top 100 undetermined barcodes",))
        for i,line in enumerate(fp):
            ws.append_row(data=line.rstrip('\n').split('\t'))
            if i == 100:
                break
    wb.save_as_xlsx(xlsx_stats_file)
    report("Wrote %s" % xlsx_stats_file)

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
    if args.nprocs > 1:
        pool = Pool(args.nprocs)
        results = pool.map(concat_fastqs,inputs)
        pool.close()
        pool.join()
    else:
        results = map(concat_fastqs,inputs)

    # Done
    report("Removing %s" % tmp_dir)
    shutil.rmtree(tmp_dir)
    report("Finished")
