#!/usr/bin/env python
#
#     build_index.py: utility to build indexes for aligners
#     Copyright (C) University of Manchester 2022 Peter Briggs
#

#######################################################################
# Imports
#######################################################################

import argparse
from bcftbx.JobRunner import fetch_runner
from auto_process_ngs.indexes import IndexBuilder
from auto_process_ngs.settings import Settings
from auto_process_ngs import get_version

__version__ = get_version()
__settings = Settings()

#######################################################################
# Main program
#######################################################################

def main():
    # Create parser to handle command line
    p = argparse.ArgumentParser(description="Generate indexes for "
                                "aligners")
    p.add_argument('-v','--version',action='version',
                   version="%(prog)s "+__version__)
    p.add_argument('aligner',metavar="ALIGNER",
                   choices=['bowtie','bowtie2','star'],
                   help="aligner to build index for")
    p.add_argument('fasta',metavar="FASTA",
                   help="FASTA file with sequence")
    p.add_argument('annotation',metavar="ANNOTATION",nargs='?',
                   help="annotation file (for use with STAR)")
    p.add_argument('-o',dest='out_dir',action='store',
                   help="output directory for indexes")
    # Advanced options
    advanced = p.add_argument_group('Advanced options')
    advanced.add_argument('-r','--runner',metavar='RUNNER',action='store',
                          dest="runner",default=None,
                          help="explicitly specify runner definition for "
                          "building the index. RUNNER must be a valid job "
                          "runner specification e.g. "
                          "'GEJobRunner(-pe smp.pe 8)' (default: use "
                          "appropriate runner from configuration)")
    p.add_argument_group()
    args = p.parse_args()

    # Acquire runner
    if args.runner:
        # Use runner supplied on command line
        runner = fetch_runner(args.runner)
    else:
        if args.aligner in ("bowtie","bowtie2"):
            runner = __settings.runners.fastq_screen
        elif args.aligner == "star":
            runner = __settings.runners.star
    print("Job runner: %s" % runner)

    # Set up index builder
    builder = IndexBuilder(runner)

    # Build indexes
    if args.aligner == "bowtie":
        builder.bowtie(args.fasta,args.out_dir,
                       bowtie_version="1.0.0")
    elif args.aligner == "bowtie2":
        builder.bowtie2(args.fasta,args.out_dir,
                        bowtie2_version="2.4.1")
    elif args.aligner == "star":
        builder.star(args.fasta,args.annotation,args.out_dir,
                     star_version="2.4.2a")

if __name__ == '__main__':
    main()
