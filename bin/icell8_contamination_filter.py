#!/usr/bin/env python
#
#     icell8_contamination_filter.py: remove 'contaminated' reads from fastqs
#     Copyright (C) University of Manchester 2017,2021 Peter Briggs
#
"""
icell8_contamination_filter.py

Utility to process a set of FASTQ pairs from Wafergen iCell8 and
filter out read pairs where the R1 sequence:

- uniquely maps to one of a list of 'contaminant' genomes
- doesn't map at all to one of a list of 'mammalian' genomes

This utility should be run after the trimming step.

The procedure uses fastq_screen to perform the filtering.
"""

######################################################################
# Imports
######################################################################

import os
import sys
import tempfile
import logging
import argparse
import shutil
from bcftbx.utils import mkdir
from bcftbx.FASTQFile import FastqIterator
from bcftbx.utils import strip_ext
from bcftbx.utils import find_program
from auto_process_ngs.command import Command
from auto_process_ngs.utils import OutputFiles
from auto_process_ngs.icell8.utils import ICell8FastqIterator

import logging
logging.basicConfig(format='%(levelname) 8s: %(message)s')

######################################################################
# Functions
######################################################################

def fastq_screen_tag(conf_file,fastq_in,out_dir,
                     aligner=None,threads=1,
                     tempdir=None):
    """
    Run 'fastq_screen' and output tagged fastq file

    Raises an Exception in the event of an error.

    Arguments:
      conf_file (str): path to the fastq_screen .conf file
      fastq_in (str): path to the FASTQ file to screen
      out_dir (str): path to the output directory to put
        the tagged FASTQ in
      aligner (str): optional, name of the aligner to pass
        to fastq_screen (default: don't specify the aligner)
      threads (int): optional, the number of threads to
        use when running fastq_screen (default: 1)
      tempdir (str): optional, directory to create temporary
        working directories in when running fastq_screen

    Returns:
      String: path to the tagged output FASTQ file
    """
    # Make a temporary working directory
    work_dir = tempfile.mkdtemp(suffix='.fastq_screen',
                                dir=tempdir)
    # Build fastq_screen command
    fastq_screen_cmd = Command(
        'fastq_screen',
        '--subset',0,
        '--threads',threads,
        '--conf',conf_file,
        '--tag',
        '--outdir',work_dir)
    if args.aligner is not None:
        fastq_screen_cmd.add_args('--aligner',args.aligner)
    fastq_screen_cmd.add_args(fastq_in)
    print("Running %s" % fastq_screen_cmd)
    # Run the command
    exit_code = fastq_screen_cmd.run_subprocess(working_dir=work_dir)
    if exit_code != 0:
        err_msg = "Screening %s against %s failed (exit code %d)" % \
                  (fastq_in,conf_file,exit_code)
    else:
        # Handle the outputs
        tagged_fastq = os.path.basename(strip_ext(fastq_in,'.fastq')) \
                       + '.tagged.fastq'
        if not os.path.exists(os.path.join(work_dir,tagged_fastq)):
            err_msg = "Failed to generated tagged fastq file %s" % \
                      tagged_fastq
            exit_code = 1
        else:
            os.rename(os.path.join(work_dir,tagged_fastq),
                      os.path.join(out_dir,tagged_fastq))
    # Clean up working directory
    shutil.rmtree(work_dir)
    # Raise exception if there was a problem
    if exit_code != 0:
        raise Exception(err_msg)
    # Return path to tagged file
    return os.path.join(out_dir,tagged_fastq)

def extract_fastq_screen_tag(read):
    """
    Extract the tag string from a read tagged by fastq_screen

    Arguments:
      read (FastqRead): FASTQ read

    Returns:
      String: extracted tag
    """
    # Tags are appended to the index sequence and are of the
    # the form either "#FQST:Human:Mouse:01" (1st read in file)
    # or "#FQST:22" (subsequent reads)
    tag = read.seqid.index_sequence.split('#')[1]
    if not tag.startswith("FQST:"):
        raise Exception("Bad tag in read: %s" % read.seqid)
    return tag.split(':')[-1]

def nohits_tag_from_conf(conf_file):
    """
    Construct a 'nohits' tag from a fastq_screen conf file

    The 'nohits' tag is a string of '0' characters, with
    the number of zeroes equal to the number of 'DATABASE'
    lines in the conf file.

    For example: if there are 3 genomes in the file then
    the 'nohits' tag will look like '000'.
    """
    ndatabases = 0
    with open(conf_file,'r') as conf:
        for line in conf:
            if line.startswith("DATABASE\t"):
                ndatabases += 1
    return '0' * ndatabases

######################################################################
# Main
######################################################################

if __name__ == "__main__":
    # Handle the command line
    p = argparse.ArgumentParser()
    p.add_argument("FQ_R1",help="R1 FASTQ file")
    p.add_argument("FQ_R2",help="Matching R2 FASTQ file")
    p.add_argument("-o","--outdir",
                   dest="out_dir",default=None,
                   help="directory to write output FASTQ files to "
                   "(default: current directory)")
    p.add_argument("-m","--mammalian",
                   dest="mammalian_conf",
                   help="fastq_screen 'conf' file with the "
                   "mammalian genome indices")
    p.add_argument("-c","--contaminant",
                   dest="contaminants_conf",
                   help="fastq_screen 'conf' file with the "
                   "contaminant genome indices")
    p.add_argument("-a","--aligner",
                   dest="aligner",default=None,
                   choices=["bowtie","bowtie2"],
                   help="aligner to use with fastq_screen (default: "
                   "don't specify an aligner)")
    p.add_argument("-n","--threads",
                   type=int,default=1,
                   help="number of threads to run fastq_screen with "
                   "(default: 1)")
    args = p.parse_args()

    # Input FASTQ pair
    fqr1 = os.path.abspath(args.FQ_R1)
    fqr2 = os.path.abspath(args.FQ_R2)

    # Screen files
    mammalian_conf = args.mammalian_conf
    if mammalian_conf is not None:
        mammalian_conf = os.path.abspath(mammalian_conf)
    contaminants_conf = args.contaminants_conf
    if contaminants_conf is not None:
        contaminants_conf = os.path.abspath(contaminants_conf)

    # Check for underlying programs
    required = ["fastq_screen"]
    if args.aligner is not None:
        required.append(args.aligner)
    else:
        logging.warning("Aligner not specified, cannot check")
    for prog in required:
        if find_program(prog) is None:
            logging.critical("couldn't find '%s'" % prog)
            sys.exit(1)

    # Make output dir
    if args.out_dir is not None:
        out_dir = os.path.abspath(args.out_dir)
        mkdir(out_dir)
    else:
        out_dir = os.getcwd()

    # Screen against 'mammalian' genomes
    tagged_fastq = fastq_screen_tag(mammalian_conf,fqr2,
                                    aligner=args.aligner,
                                    threads=args.threads,
                                    out_dir=out_dir,
                                    tempdir=out_dir)
    mammalian_tagged_fq = strip_ext(tagged_fastq,'.fastq') + '.' + \
                          os.path.basename(
                              strip_ext(mammalian_conf,'.conf')) + \
                          '.fastq'
    os.rename(tagged_fastq,mammalian_tagged_fq)

    # Screen against 'contaminants' genomes
    tagged_fastq = fastq_screen_tag(contaminants_conf,fqr2,
                                    aligner=args.aligner,
                                    threads=args.threads,
                                    out_dir=out_dir,
                                    tempdir=out_dir)
    contaminants_tagged_fq = strip_ext(tagged_fastq,'.fastq') + '.' + \
                             os.path.basename(
                                 strip_ext(contaminants_conf,'.conf')) + \
                             '.fastq'
    os.rename(tagged_fastq,contaminants_tagged_fq)

    # Construct fastq_screen tags to match against
    nohits_mammalian = nohits_tag_from_conf(mammalian_conf)
    nohits_contaminants = nohits_tag_from_conf(contaminants_conf)
    ndatabases_mammalian = len(nohits_mammalian)
    ndatabases_contaminants = len(nohits_contaminants)
    print("'nohits' tags: '%s' and '%s'" % (nohits_mammalian,
                                            nohits_contaminants))
    # Output filtered FASTQ pair
    fqr1_out = os.path.basename(strip_ext(fqr1,'.fastq')) \
               + '.filtered.fastq'
    fqr2_out = os.path.basename(strip_ext(fqr2,'.fastq')) \
               + '.filtered.fastq'
    output_fqs = OutputFiles(base_dir=out_dir)
    output_fqs.open('fqr1',fqr1_out)
    output_fqs.open('fqr2',fqr2_out)

    # Filter the iCell8 read pairs against the tagged reads
    for pair,pref,contam in zip(ICell8FastqIterator(fqr1,fqr2),
                                FastqIterator(mammalian_tagged_fq),
                                FastqIterator(contaminants_tagged_fq)):
        # Get the tags
        pref_tag = extract_fastq_screen_tag(pref)
        contam_tag = extract_fastq_screen_tag(contam)
        # Check number of databases are consistent
        if len(pref_tag) != ndatabases_mammalian:
            logging.critical("Mismatch in mammalian tag: "
                             "len('%s') != len('%s')" %
                             (pref_tag,nohits_mammalian))
            sys.exit(1)
        if len(contam_tag) != ndatabases_contaminants:
            logging.critical("Mismatch in contaminant tag: "
                             "len('%s') != len('%s')" %
                             (contam_tag,nohits_contaminants))
            sys.exit(1)
        # Read passes if:
        # -- there is at least one hit on the 'mammalian' genomes, OR
        # -- there are no hits on the 'contaminants' genomes
        read_ok = ((pref_tag != nohits_mammalian) or
                   (contam_tag == nohits_contaminants))
        if read_ok:
            output_fqs.write('fqr1',pair.r1)
            output_fqs.write('fqr2',pair.r2)

    # Close the output files
    output_fqs.close()

    # Remove the tagged fastqs
    os.remove(mammalian_tagged_fq)
    os.remove(contaminants_tagged_fq)
