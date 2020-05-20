#!/usr/bin/env python
#
#     icell8_atac.py: utility functions for handling ICELL8 ATAC-seq data
#     Copyright (C) University of Manchester 2019 Peter Briggs
#
"""
Utility functions for handling single-cell ATAC-seq data from
the ICELL8 platform.

Functions:

- report: write a timestamped message
- reverse_complement: get reverse complement of a sequence
- update_fastq_read_index: rewrite index sequence in Fastq read header
- split_fastq: split Fastq into batches
- assign_reads: assign reads to samples from batched ICELL8 ATAC Fastqs
- concat_fastqs: concatenate Fastqs for a sample across batches

"""

######################################################################
# Imports
######################################################################

import sys
import os
import glob
import time
try:
    # Python 2
    from itertools import izip as zip
except ImportError:
    pass
from collections import defaultdict
from bcftbx.FASTQFile import SequenceIdentifier
from bcftbx.ngsutils import getreads
from ..applications import Command
from ..analysis import AnalysisFastq
from ..utils import BufferedOutputFiles
from ..utils import ProgressChecker
from .utils import ICell8WellList

######################################################################
# Constants
######################################################################

REVERSE_COMPLEMENTS =  { 'A':'T',
                         'T':'A',
                         'G':'C',
                         'C':'G',
                         'N':'N' }

######################################################################
# Functions
######################################################################

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

def update_fastq_read_index(read,index_sequence):
    """
    Update the index sequence (aka barcode) in a Fastq read

    Arguments:
      read (list): Fastq read to be updated, as a list
        of lines (with the first element/line being the
        sequence identifier line)
      index_sequence (str): the index sequence to put
        into the read header

    Returns:
      List: the updated Fastq read, as a list of lines.
    """
    seq_id = SequenceIdentifier(read[0])
    seq_id.index_sequence = index_sequence
    read[0] = str(seq_id)
    return read

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
            ##sys.exit(status)
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
    - rewrite_fastq_headers: boolean indicating whether to
      write the matching ICELL8 barcodes into the Fastq
      read headers on output
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
    fastq_r1,fastq_r2,fastq_i1,fastq_i2,well_list_file,mode,swap_i1_and_i2,reverse_complement_index,rewrite_fastq_headers,working_dir,unassigned = args
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
    if rewrite_fastq_headers:
        report("[%s] Rewriting Fastq read headers to include well list "
               "barcodes" % label)
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
            for r1,r2,i1,i2 in zip(getreads(fastq_r1),
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
                # Rewrite read headers to include well list barcode
                if rewrite_fastq_headers and barcode:
                    r1 = update_fastq_read_index(r1,barcode)
                    r2 = update_fastq_read_index(r2,barcode)
                    i1 = update_fastq_read_index(i1,barcode)
                    i2 = update_fastq_read_index(i2,barcode)
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
            for r1,r2,i1,i2 in zip(getreads(fastq_r1),
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
                # Rewrite read headers to include well list barcode
                if rewrite_fastq_headers and barcode:
                    r1 = update_fastq_read_index(r1,barcode)
                    r2 = update_fastq_read_index(r2,barcode)
                    i1 = update_fastq_read_index(i1,barcode)
                    i2 = update_fastq_read_index(i2,barcode)
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
    - lane: (optional) lane number for output Fastq
      (set to None to stop lane number appearing)
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
    sample,index,barcode,lane,read,batches,working_dir,final_dir = args
    label = "%s:%s" % (sample,read)
    if barcode:
        report("[%s] Concatenating batched %s Fastqs for sample '%s' (%s)" %
               (label,read,sample,barcode))
    else:
        report("[%s] Concatenating batched %s Fastqs for sample '%s'" %
               (label,read,sample))
    # Name of output file
    fastq = os.path.join(final_dir,
                         "%s%s_S%d_%s%s_001.fastq.gz" %
                         (sample,
                          "-%s_" % barcode if barcode else "",
                          index,
                          "L%03d_" % int(lane) if lane else "",
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
