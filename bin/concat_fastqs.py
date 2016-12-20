#!/usr/bin/env python
#
#     concat_fastqs.py: concatenate multiple fastq files
#     Copyright (C) University of Manchester 2016 Peter Briggs
#
"""
concat_fastqs.py

Concatenate a list of fastq files.

"""

import optparse
import os
import sys
import logging
from bcftbx.utils import concatenate_fastq_files
import auto_process_ngs

__version__ = auto_process_ngs.get_version()

if __name__ == "__main__":
    # Handle command line
    p = optparse.OptionParser(
        usage="%prog [OPTIONS] FASTQ [FASTQ...] FASTQ_OUT",
        description="Concatenate reads from one or more input Fastq files "
        "into a single new file FASTQ_OUT.",
        version=__version__)
    p.add_option('-v','--verbose',
                 action="store_true",dest="verbose",
                 default=False,
                 help="verbose output")
    opts,args = p.parse_args()
    # Sort out inputs
    if len(args) < 2:
        p.error("Need to supply at least 2 input fastqs plus output name")
    fastq_out = args[-1]
    fastqs = args[:-1]
    # Check inputs exist
    for fq in fastqs:
        if not os.path.exists(fq):
            logging.critical("Input file '%s' not found" % fq)
            sys.exit(1)
    # Run the concatenation
    try:
        concatenate_fastq_files(fastq_out,fastqs,
                                verbose=opts.verbose)
    except Exception as ex:
        logging.critical("Failed with exception: %s" % ex)
        sys.exit(1)
