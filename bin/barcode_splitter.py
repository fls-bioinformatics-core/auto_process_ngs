#!/usr/bin/env python
#
#     barcode_splitter.py: split reads into fastq files
#     Copyright (C) University of Manchester 2014-15,2019 Peter Briggs
#
#########################################################################
#
# barcode_splitter.py
#
#########################################################################

"""
barcode_splitter.py

Split reads into fastq files based on matching barcode (index) sequences;
can handle single- and paired-end inputs.

The input fastqs can be supplied as an explicit list of files, or as a
single directory which is assumed to contain the outputs from bclToFastq.
In this case the fastqs will be collected automatically.

The files are processed a read at a time and the barcode sequences in the
read headers are matched against those supplied via one or more sequences
supplied via the -b/--barcode argument.

For each barcode there will be an output file called BARCODE.fastq
"""

__version__ = "0.0.7"

import optparse
import logging
import sys
from auto_process_ngs.barcodes.splitter import BarcodeMatcher
from auto_process_ngs.barcodes.splitter import get_fastqs_from_dir
from auto_process_ngs.barcodes.splitter import split_single_end
from auto_process_ngs.barcodes.splitter import split_paired_end

#######################################################################
# Main program
#######################################################################

if __name__ == "__main__":
    p = optparse.OptionParser(usage="\n\t%prog [OPTIONS] FASTQ [FASTQ...]\n"
                              "\t%prog [OPTIONS] FASTQ_R1,FASTQ_R2 [FASTQ_R1,FASTQ_R2...]\n"
                              "\t%prog [OPTIONS] DIR",
                              description="Split reads from one or more input Fastq files "
                              "into new Fastqs based on matching supplied barcodes.")
    p.add_option('-b','--barcode',action='append',dest='index_seq',
                 help="specify index sequence to filter using")
    p.add_option('-m','--mismatches',action='store',dest='n_mismatches',type='int',default=0,
                 help="maximum number of differing bases to allow for two index sequences "
                 "to count as a match. Default is zero i.e. exact matches only")
    p.add_option('-n','--name',action='store',dest='base_name',default=None,
                 help="basename to use for output files")
    p.add_option('-o','--output-dir',action='store',dest='out_dir',
                 help="specify directory for output split Fastqs")
    p.add_option('-u','--unaligned',action='store',dest='unaligned_dir',default=None,
                 help="specify subdirectory with outputs from bcl-to-fastq")
    p.add_option('-l','--lane',action='store',dest='lane',default=None,type='int',
                 help="specify lane to collect and split Fastqs for")
    p.add_option('-p','--paired-end',action='store_true',dest='paired_end',
                 help="input arguments are pairs of Fastq files **NB** deprecated, pairs "
                 "are detected automatically")

    options,args = p.parse_args()
    if options.index_seq is None:
        p.error("No index sequences specified")
    if len(args) == 0:
        p.error("No input files specified")

    matcher = BarcodeMatcher(options.index_seq,
                             max_dist=options.n_mismatches)

    if len(args) == 1 and os.path.isdir(args[0]):
        if options.lane is None:
            p.error("Must supply a lane (-l option)")
        fastqs = get_fastqs_from_dir(args[0],
                                     lane=options.lane,
                                     unaligned_dir=options.unaligned_dir)
    else:
        fastqs = args

    paired_end = ',' in fastqs[0]

    if not paired_end:
        split_single_end(matcher,fastqs,
                         base_name=options.base_name,
                         output_dir=options.out_dir)
    else:
        fastq_pairs = [x.split(',') for x in fastqs] 
        split_paired_end(matcher,fastq_pairs,
                         base_name=options.base_name,
                         output_dir=options.out_dir)

    
