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
import gzip
from bcftbx.FASTQFile import FastqIterator

def assign_barcodes_single_end(fastq_in,fastq_out,n=5):
    """
    """
    if fastq_out.endswith('.gz'):
        fp = gzip.GzipFile(filename=fastq_out,mode='wb')
    else:
        fp = open(fastq_out,'w')
    print "Processing reads from %s" % fastq_in
    nread = 0
    for read in FastqIterator(fastq_in):
        # Extract new barcode sequence
        barcode = read.sequence[:n]
        # Truncate sequence and quality accordingly
        sequence = read.sequence[n:]
        quality = read.quality[n:]
        # Assign new values and write to output
        read.seqid.index_sequence = barcode
        read.sequence = sequence
        read.quality = quality
        fp.write("%s\n" % read)
        nread += 1
    print "Finished (%d reads processed)" % nread

fastq_r1 = """@MISEQ:34:000000000-A7PHP:1:1101:12552:1774 1:N:0:TAAGGCGA
TTTACAACTAGCTTCTCTTTTTCTT
+
>AA?131@C1FCGGGG1BFFGF1F3
@MISEQ:34:000000000-A7PHP:1:1101:16449:1793 1:N:0:TAAGGCGA
TCCCCAGTCTCAGCCCTACTCCACT
+
11>>>11CDFFFBCGGG1AFE11AB
@MISEQ:34:000000000-A7PHP:1:1101:15171:1796 1:N:0:GCCTTACC
CCACCACGCCTGGCTAATTTTTTTT
+
1>1>>AAAAAFAGGGGBFGGGGGG0
@MISEQ:34:000000000-A7PHP:1:1101:18777:1797 1:N:0:TAAGNNNA
CAGCAATATACACTTCACTCTGCAT
+
111>>1FBFFFFGGGGCBGEG3A3D
@MISEQ:34:000000000-A7PHP:1:1101:18622:1812 1:N:0:TAAGGCGA
GATAAAGACAGAGTCTTAATTAAAC
+
11>1>11DFFCFFDGGGB3BF313A
"""

fastq_r1_out = """@MISEQ:34:000000000-A7PHP:1:1101:12552:1774 1:N:0:TTTAC
AACTAGCTTCTCTTTTTCTT
+
31@C1FCGGGG1BFFGF1F3
@MISEQ:34:000000000-A7PHP:1:1101:16449:1793 1:N:0:TCCCC
AGTCTCAGCCCTACTCCACT
+
11CDFFFBCGGG1AFE11AB
@MISEQ:34:000000000-A7PHP:1:1101:15171:1796 1:N:0:CCACC
ACGCCTGGCTAATTTTTTTT
+
AAAAAFAGGGGBFGGGGGG0
@MISEQ:34:000000000-A7PHP:1:1101:18777:1797 1:N:0:CAGCA
ATATACACTTCACTCTGCAT
+
1FBFFFFGGGGCBGEG3A3D
@MISEQ:34:000000000-A7PHP:1:1101:18622:1812 1:N:0:GATAA
AGACAGAGTCTTAATTAAAC
+
11DFFCFFDGGGB3BF313A
"""

import unittest
import os
import tempfile
import shutil
class TestAssignBarcodesSingleEnd(unittest.TestCase):
    def setUp(self):
        # Temporary working dir
        self.wd = tempfile.mkdtemp(suffix='.test_split_single_end')
        # Test file
        self.fastq_in = os.path.join(self.wd,'test.fq')
        open(self.fastq_in,'w').write(fastq_r1)
        # Output file
        self.fastq_out = os.path.join(self.wd,'out.fq')
    def tearDown(self):
        # Remove temporary working dir
        if os.path.isdir(self.wd):
            shutil.rmtree(self.wd)
    def test_assign_barcodes_single_end(self):
        assign_barcodes_single_end(self.fastq_in,
                                   self.fastq_out)
        self.assertEqual(open(self.fastq_out,'r').read(),
                         fastq_r1_out)

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


        
        
