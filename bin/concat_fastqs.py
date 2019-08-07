#!/usr/bin/env python
#
#     concat_fastqs.py: concatenate multiple fastq files
#     Copyright (C) University of Manchester 2016,2019 Peter Briggs
#
"""
concat_fastqs.py

Concatenate a list of fastq files.

"""

import argparse
import os
import sys
import logging
from bcftbx.utils import concatenate_fastq_files
import auto_process_ngs

__version__ = auto_process_ngs.get_version()

if __name__ == "__main__":
    # Handle command line
    p = argparse.ArgumentParser(
        description="Concatenate reads from one or more input Fastq files "
        "into a single new file FASTQ_OUT")
    p.add_argument('--version',
                   action='version',version=__version__)
    p.add_argument('-v','--verbose',
                   action="store_true",dest="verbose",
                   default=False,
                   help="verbose output")
    p.add_argument('fastqs',
                   metavar="FASTQ",nargs='+',
                   help="Input FASTQ to concatenate")
    p.add_argument('fastq_out',
                   metavar="FASTQ_OUT",
                   help="Output FASTQ with concatenated reads")
    args = p.parse_args()
    # Sort out inputs
    if len(args.fastqs) < 2:
        p.error("Need to supply at least 2 input fastqs plus output name")
    # Check inputs exist
    for fq in args.fastqs:
        if not os.path.exists(fq):
            logging.critical("Input file '%s' not found" % fq)
            sys.exit(1)
    # Run the concatenation
    try:
        concatenate_fastq_files(args.fastq_out,args.fastqs,
                                verbose=args.verbose)
    except Exception as ex:
        logging.critical("Failed with exception: %s" % ex)
        sys.exit(1)
