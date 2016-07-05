#!/usr/bin/env python
#
#     assign_barcodes.py: assign barcodes to read headers in fastqs
#     Copyright (C) University of Manchester 2016 Peter Briggs
#
"""
assign_barcodes.py

Utility to extract arbitrary sequence fragments from reads in FASTQ
format files and assign these to the index sequence in the read
headers.

"""
import optparse
from auto_process_ngs.fastq_utils import assign_barcodes_single_end

if __name__ == '__main__':
    p = optparse.OptionParser(usage="%prog [OPTIONS] INPUT.fq OUTPUT.fq",
                              description="Extract arbitrary sequence "
                              "fragments from reads in INPUT.fq FASTQ "
                              "file and assign these as the index "
                              "(barcode) sequences in the read headers "
                              "in OUTPUT.fq.")
    p.add_option('-n',action='store',dest='n',type='int',default=5,
                 help="remove first N bases from each read and assign "
                 "these as barcode index sequence (default: 5)")
    opts,args = p.parse_args()
    if len(args) != 2:
        p.error('Need to specify input and output files')
    assign_barcodes_single_end(args[0],args[1],opts.n)
