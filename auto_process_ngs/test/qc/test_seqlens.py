#######################################################################
# Unit tests for qc/seqlens.py
#######################################################################

import unittest
import shutil
import os
import tempfile

from auto_process_ngs.qc.seqlens import SeqLens
from auto_process_ngs.qc.seqlens import get_sequence_lengths

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

class TestSeqLens(unittest.TestCase):
    """
    Tests for the SeqLens class
    """
    def test_seqlens(self):
        """
        SeqLens: load data from dictionary
        """
        seq_len_data = {
            'fastq': '/data/test.fastq',
            'frac_reads_masked': 10.0,
            'frac_reads_padded': 10.0,
            'max_length': 65,
            'mean_length': 42.3,
            'median_length': 33,
            'min_length': 22,
            'nreads': 10,
            'nreads_masked': 1,
            'nreads_padded': 1,
            'seq_lengths_dist': {22: 2, 29: 1, 33: 2, 41: 1, 48: 1, 65: 3},
            'seq_lengths_masked_dist': {22: 1},
            'seq_lengths_padded_dist': {22: 1}
        }
        s = SeqLens(data=seq_len_data)
        self.assertEqual(s.data,seq_len_data)
        self.assertEqual(s.fastq,"/data/test.fastq")
        self.assertEqual(s.nreads,10)
        self.assertEqual(s.min_length,22)
        self.assertEqual(s.max_length,65)
        self.assertEqual(s.mean,42.3)
        self.assertEqual(s.range,"22-65")
        self.assertEqual(s.nmasked,1)
        self.assertEqual(s.frac_masked,10.0)
        self.assertEqual(s.npadded,1)
        self.assertEqual(s.frac_padded,10.0)
        self.assertEqual(s.dist,{22: 2, 29: 1, 33: 2, 41: 1, 48: 1, 65: 3})
        self.assertEqual(s.masked_dist,{22: 1})
        self.assertTrue(s)

class TestGetSequenceLengths(unittest.TestCase):
    """
    Tests for the get_sequence_lengths function
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestGetSequenceLengths')

    def tearDown(self):
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    def test_get_sequence_lengths_from_fastq(self):
        """
        get_sequence_lengths: test outputs with Fastq file
        """
        # Make a test Fastq file
        fastq = os.path.join(self.wd,"test.fastq")
        with open(fastq,'wt') as fp:
            fp.write("""@MISEQ:1:000000000-A2Y1L:1:1101:19264:2433 1:N:0:NAAGGCGATAGATCGC
AGATAGCCGAAGATAAAGAGNTCATAACCGTAA
+
?????BBB@BBBB?BBFFFF#66EAFHHHCEFE
@MISEQ:1:000000000-A2Y1L:1:1101:18667:2435 1:N:0:NAAGGCGATAGATCGC
ATATATTCATCCGCCATTATNA
+
?????BBBDDDDADDDE@FF#6
@MISEQ:1:000000000-A2Y1L:1:1101:17523:2436 1:N:0:NAAGGCGATAGATCGC
CATCACTACCGCTCAGGAATNTGACGGCA
+
?<,<?BBBBBBBBBBBFFFF#6ACECCEC
@MISEQ:1:000000000-A2Y1L:1:1101:15489:2437 1:N:0:NAAGGCGATAGATCGC
GAGCAGTCGGGCTCAGCGCTNTGCAAATTCTAGTTAGAAACTCACAGT
+
5====>/<@@@@@@>@CCCE#66>ACEEEEGGGGGGGFFFEFDFFFFF
@MISEQ:1:000000000-A2Y1L:1:1101:18851:2442 1:N:0:NAAGGCGATAGATCGC
GGTATCCCCCGGCAGTGAGGATGGAGCCATGGTCTGCATCA
+
??,<?BBBDDDDDDD<FFF@FC;FFFBEFHHHCDDHHGHHH
@MISEQ:1:000000000-A2Y1L:1:1101:15290:2442 1:N:0:NAAGGCGATAGATCGC
AAAATAATCCTAAAAAATAACCTCTATGCCGCC
+
?????BBBDDDDDDDDGGGGGGIIIHHFFHHHH
@MISEQ:1:000000000-A2Y1L:1:1101:18106:2444 1:N:0:NAAGGCGATAGATCGC
GTAGTATTCTCATATCACAAGT
+
55,,5?9BBBBB<<BBFFFFFF
/:*:ACE?0:::A::***00::*/?C888??EEE#############
@MISEQ:1:000000000-A2Y1L:1:1101:15892:2446 1:N:0:NAAGGCGATAGATCGC
CTTCCCCACGGCCCAGACACAAGAGACGACCTCCATAAATCTTTTAGA
+
?????BBBDBDDDDDDFFFFFFHIHIHHHHHHIHIFGGHFHHHHIIFH
@MISEQ:1:000000000-A2Y1L:1:1101:17903:2450 1:N:0:TAAGGCGATAGATCGC
GTGCAGGGGGTGTGGTCAATCCACACTGTTGCTGAG
+
=5===<>+5<5<+5=@CC;8CEEEEE;-8ACFDE.7
@MISEQ:1:000000000-A2Y1L:1:1101:15113:2451 1:N:0:TAAGGCGATAGATCGC
TCTCAGATGAGCATGCAGCAGCCCAGACTCGCCCCACGCAGTTTGCCA
+
=,,<=>>>@@@@@9@@CCEE@EE+++6C8-++CECE+>DCC>@@EFFF
""")
        # Expected output
        expected_stats = {
            'fastq': fastq,
            'frac_reads_masked': 0.0,
            'frac_reads_padded': 0.0,
            'max_length': 65,
            'mean_length': 42.3,
            'median_length': 33,
            'min_length': 22,
            'nreads': 10,
            'nreads_masked': 0,
            'nreads_padded': 0,
            'seq_lengths_dist': {22: 2, 29: 1, 33: 2, 41: 1, 48: 1, 65: 3},
            'seq_lengths_masked_dist': {},
            'seq_lengths_padded_dist': {}
        }
        # Get statistics
        stats = get_sequence_lengths(fastq)
        # Check against expected
        self.assertEqual(stats,expected_stats)

    def test_get_sequence_lengths(self):
        """
        get_sequence_lengths: Fastq with masking and padding
        """
        # Make a test Fastq file
        fastq = os.path.join(self.wd,"test.fastq")
        with open(fastq,'wt') as fp:
            fp.write("""@MISEQ:1:000000000-A2Y1L:1:1101:19264:2433 1:N:0:NAAGGCGATAGATCGC
AGATAGCCGAAGATAAAGAGNTCATAACCGTAA
+
?????BBB@BBBB?BBFFFF#66EAFHHHCEFE
@MISEQ:1:000000000-A2Y1L:1:1101:18667:2435 1:N:0:NAAGGCGATAGATCGC
NNNNNNNNNNNNNNNNNNNNNN
+
?????BBBDDDDADDDE@FF#6
@MISEQ:1:000000000-A2Y1L:1:1101:17523:2436 1:N:0:NAAGGCGATAGATCGC
CATCACTACCGCTCAGGAATNTGACGGCA
+
?<,<?BBBBBBBBBBBFFFF#6ACECCEC
@MISEQ:1:000000000-A2Y1L:1:1101:15489:2437 1:N:0:NAAGGCGATAGATCGC
GAGCAGTCGGGCTCAGCGCTNTGCAAATTCTAGTTAGAAACTCACAGT
+
5====>/<@@@@@@>@CCCE#66>ACEEEEGGGGGGGFFFEFDFFFFF
@MISEQ:1:000000000-A2Y1L:1:1101:18851:2442 1:N:0:NAAGGCGATAGATCGC
GGTATCCCCCGGCAGTGAGGATGGAGCCATGGTCTGCATCA
+
??,<?BBBDDDDDDD<FFF@FC;FFFBEFHHHCDDHHGHHH
@MISEQ:1:000000000-A2Y1L:1:1101:15290:2442 1:N:0:NAAGGCGATAGATCGC
AAAATAATCCTAAAAAATAACCTCTATGCCGCC
+
?????BBBDDDDDDDDGGGGGGIIIHHFFHHHH
@MISEQ:1:000000000-A2Y1L:1:1101:18106:2444 1:N:0:NAAGGCGATAGATCGC
GTAGTATTCTCATATNNNNNNN
+
55,,5?9BBBBB<<BBFFFFFF
/:*:ACE?0:::A::***00::*/?C888??EEE#############
@MISEQ:1:000000000-A2Y1L:1:1101:15892:2446 1:N:0:NAAGGCGATAGATCGC
CTTCCCCACGGCCCAGACACAAGAGACGACCTCCATAAATCTTTTAGA
+
?????BBBDBDDDDDDFFFFFFHIHIHHHHHHIHIFGGHFHHHHIIFH
@MISEQ:1:000000000-A2Y1L:1:1101:17903:2450 1:N:0:TAAGGCGATAGATCGC
GTGCAGGGGGTGTGGTCAATCCACACTGTTGCTGAG
+
=5===<>+5<5<+5=@CC;8CEEEEE;-8ACFDE.7
@MISEQ:1:000000000-A2Y1L:1:1101:15113:2451 1:N:0:TAAGGCGATAGATCGC
TCTCAGATGAGCATGCAGCAGCCCAGACTCGCCCCACGCAGTTTGCCA
+
=,,<=>>>@@@@@9@@CCEE@EE+++6C8-++CECE+>DCC>@@EFFF
""")
        # Expected output
        expected_stats = {
            'fastq': fastq,
            'frac_reads_masked': 10.0,
            'frac_reads_padded': 10.0,
            'max_length': 65,
            'mean_length': 42.3,
            'median_length': 33,
            'min_length': 22,
            'nreads': 10,
            'nreads_masked': 1,
            'nreads_padded': 1,
            'seq_lengths_dist': {22: 2, 29: 1, 33: 2, 41: 1, 48: 1, 65: 3},
            'seq_lengths_masked_dist': {22: 1},
            'seq_lengths_padded_dist': {22: 1}
        }
        # Get statistics
        stats = get_sequence_lengths(fastq)
        print(stats['frac_reads_masked'])
        print(stats['frac_reads_padded'])
        # Check against expected
        self.assertEqual(stats,expected_stats)

    def test_get_sequence_lengths_from_empty_fastq(self):
        """
        get_sequence_lengths: test outputs with empty Fastq file
        """
        # Make a test Fastq file
        fastq = os.path.join(self.wd,"test.fastq")
        with open(fastq,'wt') as fp:
            fp.write("")
        # Expected output
        expected_stats = {
            'fastq': fastq,
            'frac_reads_masked': None,
            'frac_reads_padded': None,
            'max_length': None,
            'mean_length': None,
            'median_length': None,
            'min_length': None,
            'nreads': 0,
            'nreads_masked': 0,
            'nreads_padded': 0,
            'seq_lengths_dist': {},
            'seq_lengths_masked_dist': {},
            'seq_lengths_padded_dist': {}
        }
        # Get statistics
        stats = get_sequence_lengths(fastq)
        # Check against expected
        self.assertEqual(stats,expected_stats)

