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

import bcftbx.IlluminaData as IlluminaData
import bcftbx.FASTQFile as FASTQFile
from auto_process_ngs.utils import OutputFiles

#########################################################################
# Classes
#########################################################################

import itertools
class HammingMetrics(object):
    """Calculate Hamming distances between two strings
    
    http://en.wikipedia.org/wiki/Hamming_distance
    
    """
    @classmethod
    def hamming_distance(self,s1,s2):
        """Return the Hamming distance between equal-length sequences

        """
        if len(s1) != len(s2):
            raise ValueError("Undefined for sequences of unequal length")
        return sum(ch1 != ch2 for ch1, ch2 in itertools.izip(s1, s2))
    @classmethod
    def hamming_distance_truncate(self,s1,s2):
        """Return the Hamming distance between non-equal-length sequences

        """
        if len(s1) != len(s2):
            l = min(len(s1),len(s2))
            s1 = s1[:l]
            s2 = s2[:l]
        return sum(ch1 != ch2 for ch1, ch2 in itertools.izip(s1, s2))
    @classmethod
    def hamming_distance_with_N(self,s1,s2):
        """Return the Hamming distance between equal-length sequences

        """
        # modified version where N's don't match anything (not even other N's)
        if len(s1) != len(s2):
            raise ValueError("Undefined for sequences of unequal length")
        return sum((ch1 != ch2 or ch1 == 'N' or ch2 == 'N') \
                   for ch1, ch2 in itertools.izip(s1, s2))

class HammingLookup(object):
    """
    """
    def __init__(self,hamming_func=HammingMetrics.hamming_distance):
        self._hamming = hamming_func
        self._distances = dict()
    def dist(self,s1,s2):
        try:
            return self._distances[s1][s2]
        except KeyError:
            try:
                return self._distances[s2][s1]
            except KeyError:
                pass
        if s1 not in self._distances:
            self._distances[s1] = dict()
        self._distances[s1][s2] = self._hamming(s1,s2)
        return self._distances[s1][s2]

class BarcodeMatcher(object):
    """
    """
    def __init__(self,index_seqs,max_dist=0):
        self._hamming = HammingLookup(hamming_func=HammingMetrics.hamming_distance_truncate)
        self._index_seqs = list(index_seqs)
        self._max_dist = max_dist
        for i,seq1 in enumerate(self._index_seqs):
            for seq2 in self._index_seqs[i+1:]:
                if self._hamming.dist(seq1,seq2) <= max_dist:
                    raise Exception,"Index sequence ambiguity: '%s' and '%s' are too " \
                    "similar (differ by %d bases, must be > %d)" % \
                    (seq1,seq2,self._hamming.dist(seq1,seq2),max_dist)
        self._index_seqs.sort()
    def match(self,seq):
        for index_seq in self._index_seqs:
            logging.debug("%s, %s" % (index_seq,seq))
            if self._hamming.dist(index_seq,seq) <= self._max_dist:
                return index_seq
        return None
    @property
    def sequences(self):
        return self._index_seqs

def get_fastqs_from_dir(dirn,lane,unaligned_dir=None):
    """Automatically collect Fastq files for specified lane

    """
    try:
        illumina_data = IlluminaData.IlluminaData(dirn,
                                                  unaligned_dir=unaligned_dir)
    except Exception,ex:
        sys.stderr.write("Unable to read fastqs from %s: %s\n" % (dirn,ex))
        sys.exit(1)
    paired_end = illumina_data.paired_end
    fastqs_r1 = []
    fastqs_r2 = []
    for project in illumina_data.projects:
        for sample in project.samples:
            for fastq in sample.fastq_subset(read_number=1,full_path=True):
                if IlluminaData.IlluminaFastq(fastq).lane_number == lane:
                    fastqs_r1.append(fastq)
            for fastq in sample.fastq_subset(read_number=2,full_path=True):
                if IlluminaData.IlluminaFastq(fastq).lane_number == lane:
                    fastqs_r2.append(fastq)
    if illumina_data.undetermined:
        for sample in illumina_data.undetermined.samples:
            for fastq in sample.fastq_subset(read_number=1,full_path=True):
                if IlluminaData.IlluminaFastq(fastq).lane_number == lane:
                    fastqs_r1.append(fastq)
            for fastq in sample.fastq_subset(read_number=2,full_path=True):
                if IlluminaData.IlluminaFastq(fastq).lane_number == lane:
                    fastqs_r2.append(fastq)
    if not paired_end:
        return fastqs_r1
    fastqs = []
    fastqs_r1.sort()
    fastqs_r2.sort()
    for fq1,fq2 in zip(fastqs_r1,fastqs_r2):
        fastqs.append("%s,%s" % (fq1,fq2))
    return fastqs

def split_single_end(matcher,fastqs,base_name=None,output_dir=None):
    """Split reads from single ended data

    For each fastq file in 'fastqs', check reads against the index
    sequences in the BarcodeMatcher 'matcher' and write to an
    appropriate file.

    """
    if base_name is None:
        base_name = ''
    else:
        base_name = "%s." % base_name
    fp = OutputFiles(base_dir=output_dir)
    for barcode in matcher.sequences:
        fp.open(barcode,"%s%s.fastq" % (base_name,barcode))
    fp.open('undetermined',"%sundetermined.fastq" % base_name)
    # Filter reads
    nread = 0
    for fastq in fastqs:
        print "Processing reads from %s" % fastq
        for read in FASTQFile.FastqIterator(fastq):
            nread += 1
            seq = read.seqid.index_sequence
            if not seq:
                logging.error("No index sequence for read!")
                sys.exit(1)
            assigned_index = matcher.match(seq)
            # Read not assigned
            if assigned_index is None:
                assigned_index = 'undetermined'
            logging.debug("Assigned read #%d to %s" % (nread,assigned_index))
            fp.write(assigned_index,read)
    print "Finished (%d reads processed)" % nread

def split_paired_end(matcher,fastq_pairs,base_name=None,output_dir=None):
    """Split reads from paired end data

    For each fastq file pair in 'fastqs', check reads against the
    index sequences in the BarcodeMatcher 'matcher' and write to an
    appropriate file.

    """
    if base_name is None:
        base_name = ''
    else:
        base_name = "%s." % base_name
    fp = OutputFiles(base_dir=output_dir)
    for barcode in matcher.sequences:
        fp.open((barcode,'R1'),"%s%s_R1.fastq" % (base_name,barcode))
        fp.open((barcode,'R2'),"%s%s_R2.fastq" % (base_name,barcode))
    fp.open(('undetermined','R1'),"%sundetermined_R1.fastq" % base_name)
    fp.open(('undetermined','R2'),"%sundetermined_R2.fastq" % base_name)
    # Filter reads
    nread = 0
    for fq_r1,fq_r2 in fastq_pairs:
        print "Processing reads from fastq pair %s %s" % (fq_r1,fq_r2)
        for read1,read2 in itertools.izip(FASTQFile.FastqIterator(fq_r1),
                                          FASTQFile.FastqIterator(fq_r2)):
            nread += 1
            seq = read1.seqid.index_sequence
            if not seq:
                logging.error("No index sequence for read!")
                sys.exit(1)
            if seq != read2.seqid.index_sequence:
                raise Exception,"Index sequence mismatch between R1 and R2 reads"
            assigned_index = matcher.match(seq)
            # Read not assigned
            if assigned_index is None:
                assigned_index = 'undetermined'
            logging.debug("Assigned read #%d to %s" % (nread,assigned_index))
            fp.write((assigned_index,'R1'),read1)
            fp.write((assigned_index,'R2'),read2)
    print "Finished (%d read pairs processed)" % nread

#######################################################################
# Unit tests
#######################################################################

import unittest
class TestHammingMetrics(unittest.TestCase):
    def test_hamming_exact(self):
        self.assertEqual(HammingMetrics.hamming_distance('AGGTCTA','AGGTCTA'),0)
    def test_hamming_one_mismatch(self):
        self.assertEqual(HammingMetrics.hamming_distance('AGGTCTA','ACGTCTA'),1)
    def test_hamming_two_mismatches(self):
        self.assertEqual(HammingMetrics.hamming_distance('AGGTCTA','ACCTCTA'),2)
class TestHammingWithN(unittest.TestCase):
    def test_hamming_with_N_exact(self):
        self.assertEqual(HammingMetrics.hamming_distance_with_N('AGGTCTA','AGGTCTA'),0)
    def test_hamming_with_N_one_mismatch(self):
        self.assertEqual(HammingMetrics.hamming_distance_with_N('AGGTCTA','ACGTCTA'),1)
        self.assertEqual(HammingMetrics.hamming_distance_with_N('ACNTCTA','ACGTCTA'),1)
        self.assertEqual(HammingMetrics.hamming_distance_with_N('ACNTCTA','ACNTCTA'),1)
    def test_hamming_with_N_two_mismatches(self):
        self.assertEqual(HammingMetrics.hamming_distance_with_N('AGGTCTA','ACCTCTA'),2)
        self.assertEqual(HammingMetrics.hamming_distance_with_N('ANNTCTA','ACCTCTA'),2)
        self.assertEqual(HammingMetrics.hamming_distance_with_N('ANNTCTA','ANNTCTA'),2)
    def test_hamming_different_lengths(self):
        self.assertRaises(ValueError,HammingMetrics.hamming_distance_with_N,
                          'AGGTCTAA','AGGTCTA')

class TestBarcodeMatcher(unittest.TestCase):
    def test_barcodematcher_exact_match(self):
        b = BarcodeMatcher(('TTGCTA','AGGTCT'))
        self.assertEqual(b.match('AGGTCT'),'AGGTCT')
        self.assertEqual(b.match('TTGCTA'),'TTGCTA')
        self.assertEqual(b.match('CGGTCT'),None)
        self.assertEqual(b.match('CGTTCT'),None)
    def test_barcodematcher_one_mismatch(self):
        b = BarcodeMatcher(('TTGCTA','AGGTCT'),max_dist=1)
        self.assertEqual(b.match('AGGTCT'),'AGGTCT')
        self.assertEqual(b.match('TTGCTA'),'TTGCTA')
        self.assertEqual(b.match('CGGTCT'),'AGGTCT')
        self.assertEqual(b.match('CGTTCT'),None)
    def test_barcodematcher_two_mismatches(self):
        b = BarcodeMatcher(('TTGCTA','AGGTCT'),max_dist=2)
        self.assertEqual(b.match('AGGTCT'),'AGGTCT')
        self.assertEqual(b.match('TTGCTA'),'TTGCTA')
        self.assertEqual(b.match('CGGTCT'),'AGGTCT')
        self.assertEqual(b.match('CGTTCT'),'AGGTCT')
    def test_barcodematcher_list_sequences(self):
        b = BarcodeMatcher(('TTGCTA','AGGTCT','GCCTAT'))
        self.assertEqual(b.sequences,['AGGTCT','GCCTAT','TTGCTA'])
    def test_barcodematcher_ambigiuous_sequences(self):
        self.assertRaises(Exception,BarcodeMatcher,('AGGTCTA','AGGTCTA'))
        self.assertRaises(Exception,BarcodeMatcher,('AGGTCTC','AGGTCTA'),max_dist=1)
        self.assertRaises(Exception,BarcodeMatcher,('AGGTCCC','AGGTCTA'),max_dist=2)


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

fastq_r2 = """@MISEQ:34:000000000-A7PHP:1:1101:12552:1774 2:N:0:TAAGGCGA
CTTCTTCTTTTTTTTTTTGTCTATC
+
1>11>13BBD111AAE00/B2B222
@MISEQ:34:000000000-A7PHP:1:1101:16449:1793 2:N:0:TAAGGCGA
TTTTTTTTCTTTTTTTTTTCCTTTT
+
11>>1>110B3BF1A00/A01B22B
@MISEQ:34:000000000-A7PHP:1:1101:15171:1796 2:N:0:GCCTTACC
TTTCTATTCTCCTTCTCATAAAAAA
+
111113333B311B1B133331110
@MISEQ:34:000000000-A7PHP:1:1101:18777:1797 2:N:0:TAAGNNNA
TTTTTTTTTCTTCTATTCTCAGATG
+
1>1>111100BB3B3B22B222121
@MISEQ:34:000000000-A7PHP:1:1101:18622:1812 2:N:0:TAAGGCGA
TTTTTTTTCCATCTTTTTTAATTTT
+
>1>111>1013B1BF3BFE01BB22
"""

class TestSplitSingleEnd(unittest.TestCase):
    def setUp(self):
        # Temporary working dir
        self.wd = tempfile.mkdtemp(suffix='.test_split_single_end')
        # Test file
        self.fastq = os.path.join(self.wd,'test.fq')
        open(self.fastq,'w').write(fastq_r1)
    def tearDown(self):
        # Remove temporary working dir
        if os.path.isdir(self.wd):
            shutil.rmtree(self.wd)
    def test_split_single_end(self):
        matcher = BarcodeMatcher(('TAAGGCGA','GCCTTACC'))
        split_single_end(matcher,(self.fastq,),output_dir=self.wd)
        self.assertTrue(os.path.exists(os.path.join(self.wd,'TAAGGCGA.fastq')))
        self.assertTrue(os.path.exists(os.path.join(self.wd,'GCCTTACC.fastq')))
        self.assertTrue(os.path.exists(os.path.join(self.wd,'undetermined.fastq')))
        self.assertEqual(open(os.path.join(self.wd,'TAAGGCGA.fastq'),'r').read(),
                         """@MISEQ:34:000000000-A7PHP:1:1101:12552:1774 1:N:0:TAAGGCGA
TTTACAACTAGCTTCTCTTTTTCTT
+
>AA?131@C1FCGGGG1BFFGF1F3
@MISEQ:34:000000000-A7PHP:1:1101:16449:1793 1:N:0:TAAGGCGA
TCCCCAGTCTCAGCCCTACTCCACT
+
11>>>11CDFFFBCGGG1AFE11AB
@MISEQ:34:000000000-A7PHP:1:1101:18622:1812 1:N:0:TAAGGCGA
GATAAAGACAGAGTCTTAATTAAAC
+
11>1>11DFFCFFDGGGB3BF313A
""")
        self.assertEqual(open(os.path.join(self.wd,'GCCTTACC.fastq'),'r').read(),
                         """@MISEQ:34:000000000-A7PHP:1:1101:15171:1796 1:N:0:GCCTTACC
CCACCACGCCTGGCTAATTTTTTTT
+
1>1>>AAAAAFAGGGGBFGGGGGG0
""")
        self.assertEqual(open(os.path.join(self.wd,'undetermined.fastq'),'r').read(),
                         """@MISEQ:34:000000000-A7PHP:1:1101:18777:1797 1:N:0:TAAGNNNA
CAGCAATATACACTTCACTCTGCAT
+
111>>1FBFFFFGGGGCBGEG3A3D
""")

class TestSplitPairedEnd(unittest.TestCase):
    def setUp(self):
        # Temporary working dir
        self.wd = tempfile.mkdtemp(suffix='.test_split_paired_end')
        # Test file
        self.fastq_r1 = os.path.join(self.wd,'test_r1.fq')
        self.fastq_r2 = os.path.join(self.wd,'test_r2.fq')
        open(self.fastq_r1,'w').write(fastq_r1)
        open(self.fastq_r2,'w').write(fastq_r2)
    def tearDown(self):
        # Remove temporary working dir
        if os.path.isdir(self.wd):
            shutil.rmtree(self.wd)
    def test_split_paired_end(self):
        matcher = BarcodeMatcher(('TAAGGCGA','GCCTTACC'))
        split_paired_end(matcher,((self.fastq_r1,self.fastq_r2),),output_dir=self.wd)
        self.assertTrue(os.path.exists(os.path.join(self.wd,'TAAGGCGA_R1.fastq')))
        self.assertTrue(os.path.exists(os.path.join(self.wd,'TAAGGCGA_R2.fastq')))
        self.assertTrue(os.path.exists(os.path.join(self.wd,'GCCTTACC_R1.fastq')))
        self.assertTrue(os.path.exists(os.path.join(self.wd,'GCCTTACC_R2.fastq')))
        self.assertTrue(os.path.exists(os.path.join(self.wd,'undetermined_R1.fastq')))
        self.assertTrue(os.path.exists(os.path.join(self.wd,'undetermined_R2.fastq')))
        self.assertEqual(open(os.path.join(self.wd,'TAAGGCGA_R1.fastq'),'r').read(),
                         """@MISEQ:34:000000000-A7PHP:1:1101:12552:1774 1:N:0:TAAGGCGA
TTTACAACTAGCTTCTCTTTTTCTT
+
>AA?131@C1FCGGGG1BFFGF1F3
@MISEQ:34:000000000-A7PHP:1:1101:16449:1793 1:N:0:TAAGGCGA
TCCCCAGTCTCAGCCCTACTCCACT
+
11>>>11CDFFFBCGGG1AFE11AB
@MISEQ:34:000000000-A7PHP:1:1101:18622:1812 1:N:0:TAAGGCGA
GATAAAGACAGAGTCTTAATTAAAC
+
11>1>11DFFCFFDGGGB3BF313A
""")
        self.assertEqual(open(os.path.join(self.wd,'TAAGGCGA_R2.fastq'),'r').read(),
                         """@MISEQ:34:000000000-A7PHP:1:1101:12552:1774 2:N:0:TAAGGCGA
CTTCTTCTTTTTTTTTTTGTCTATC
+
1>11>13BBD111AAE00/B2B222
@MISEQ:34:000000000-A7PHP:1:1101:16449:1793 2:N:0:TAAGGCGA
TTTTTTTTCTTTTTTTTTTCCTTTT
+
11>>1>110B3BF1A00/A01B22B
@MISEQ:34:000000000-A7PHP:1:1101:18622:1812 2:N:0:TAAGGCGA
TTTTTTTTCCATCTTTTTTAATTTT
+
>1>111>1013B1BF3BFE01BB22
""")
        self.assertEqual(open(os.path.join(self.wd,'GCCTTACC_R1.fastq'),'r').read(),
                         """@MISEQ:34:000000000-A7PHP:1:1101:15171:1796 1:N:0:GCCTTACC
CCACCACGCCTGGCTAATTTTTTTT
+
1>1>>AAAAAFAGGGGBFGGGGGG0
""")
        self.assertEqual(open(os.path.join(self.wd,'GCCTTACC_R2.fastq'),'r').read(),
                         """@MISEQ:34:000000000-A7PHP:1:1101:15171:1796 2:N:0:GCCTTACC
TTTCTATTCTCCTTCTCATAAAAAA
+
111113333B311B1B133331110
""")
        self.assertEqual(open(os.path.join(self.wd,'undetermined_R1.fastq'),'r').read(),
                         """@MISEQ:34:000000000-A7PHP:1:1101:18777:1797 1:N:0:TAAGNNNA
CAGCAATATACACTTCACTCTGCAT
+
111>>1FBFFFFGGGGCBGEG3A3D
""")
        self.assertEqual(open(os.path.join(self.wd,'undetermined_R2.fastq'),'r').read(),
                         """@MISEQ:34:000000000-A7PHP:1:1101:18777:1797 2:N:0:TAAGNNNA
TTTTTTTTTCTTCTATTCTCAGATG
+
1>1>111100BB3B3B22B222121
""")

#######################################################################
# Main program
#######################################################################

import optparse
import logging
import sys

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

    
