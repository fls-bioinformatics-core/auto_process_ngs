#!/usr/bin/env python
#
#     assign_barcodes.py: assign barcodes to read headers in fastqs
#     Copyright (C) University of Manchester 2016,2019 Peter Briggs
#
"""
assign_barcodes.py

Utility to extract arbitrary sequence fragments from reads in FASTQ
format files and assign these to the index sequence in the read
headers.

"""
import argparse
from auto_process_ngs.fastq_utils import assign_barcodes_single_end

if __name__ == '__main__':
    p = argparse.ArgumentParser(
        description="Extract arbitrary sequence fragments from reads "
        "in INPUT.fq FASTQ file and assign these as the index (barcode) "
        "sequences in the read headers in OUTPUT.fq.")
    p.add_argument('-n',
                   action='store',dest='n',type=int,default=5,
                   help="remove first N bases from each read and assign "
                   "these as barcode index sequence (default: 5)")
    p.add_argument('fq_in',metavar="INPUT.fq",
                   help="Input FASTQ file")
    p.add_argument('fq_out',metavar="OUTPUT.fq",
                   help="Output FASTQ file")
    args = p.parse_args()
    assign_barcodes_single_end(args.fq_in,args.fq_out,args.n)
