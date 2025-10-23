#######################################################################
# Tests for bcl2fastq.apps.py module
#######################################################################

import unittest
from auto_process_ngs.bcl2fastq.apps import configureBclToFastq
from auto_process_ngs.bcl2fastq.apps import bcl2fastq2
from auto_process_ngs.bcl2fastq.apps import bclconvert


class TestConfigureBclToFastq(unittest.TestCase):

    def test_configure_bcl_to_fastq(self):
        """Construct 'configureBclToFastq.pl' command lines
        """
        self.assertEqual(configureBclToFastq(
            'Data/Intensities/Basecalls',
            'SampleSheet.csv').command_line,
                         ['configureBclToFastq.pl',
                          '--input-dir','Data/Intensities/Basecalls',
                          '--output-dir','Unaligned',
                          '--sample-sheet','SampleSheet.csv',
                          '--fastq-cluster-count','-1'])
        self.assertEqual(configureBclToFastq(
            'Data/Intensities/Basecalls',
            'SampleSheet.csv',
            output_dir='run/bcl2fastq').command_line,
                         ['configureBclToFastq.pl',
                          '--input-dir','Data/Intensities/Basecalls',
                          '--output-dir','run/bcl2fastq',
                          '--sample-sheet','SampleSheet.csv',
                          '--fastq-cluster-count','-1'])
        self.assertEqual(configureBclToFastq(
            'Data/Intensities/Basecalls',
            'SampleSheet.csv',
            output_dir='run/bcl2fastq',
            ignore_missing_bcl=True).command_line,
                         ['configureBclToFastq.pl',
                          '--input-dir','Data/Intensities/Basecalls',
                          '--output-dir','run/bcl2fastq',
                          '--sample-sheet','SampleSheet.csv',
                          '--fastq-cluster-count','-1',
                          '--ignore-missing-bcl'])
        self.assertEqual(configureBclToFastq(
            'Data/Intensities/Basecalls',
            'SampleSheet.csv',
            configure_bcl_to_fastq_exe=
            "/opt/bin/configureBclToFastq.pl").command_line,
                         ['/opt/bin/configureBclToFastq.pl',
                          '--input-dir','Data/Intensities/Basecalls',
                          '--output-dir','Unaligned',
                          '--sample-sheet','SampleSheet.csv',
                          '--fastq-cluster-count','-1'])

class TestBcl2Fastq(unittest.TestCase):

    def test_bcl2fastq(self):
        """Construct 'bcl2fastq' command lines for bcl2fastq v2.*
        """
        self.assertEqual(bcl2fastq2(
            '/runs/150107_NB123000_0001_ABCX',
            'SampleSheet.csv').command_line,
                         ['bcl2fastq',
                          '--runfolder-dir','/runs/150107_NB123000_0001_ABCX',
                          '--output-dir','Unaligned',
                          '--sample-sheet','SampleSheet.csv'])
        self.assertEqual(bcl2fastq2(
            '/runs/150107_NB123000_0001_ABCX',
            'SampleSheet.csv',
            output_dir='run/bcl2fastq').command_line,
                         ['bcl2fastq',
                          '--runfolder-dir','/runs/150107_NB123000_0001_ABCX',
                          '--output-dir','run/bcl2fastq',
                          '--sample-sheet','SampleSheet.csv'])
        self.assertEqual(bcl2fastq2(
            '/runs/150107_NB123000_0001_ABCX',
            'SampleSheet.csv',
            output_dir='run/bcl2fastq',
            ignore_missing_bcls=True).command_line,
                         ['bcl2fastq',
                          '--runfolder-dir','/runs/150107_NB123000_0001_ABCX',
                          '--output-dir','run/bcl2fastq',
                          '--sample-sheet','SampleSheet.csv',
                          '--ignore-missing-bcls'])
        self.assertEqual(bcl2fastq2(
            '/runs/150107_NB123000_0001_ABCX',
            'SampleSheet.csv',
            output_dir='run/bcl2fastq',
            mismatches=1,
            no_lane_splitting=True).command_line,
                         ['bcl2fastq',
                          '--runfolder-dir','/runs/150107_NB123000_0001_ABCX',
                          '--output-dir','run/bcl2fastq',
                          '--sample-sheet','SampleSheet.csv',
                          '--barcode-mismatches','1',
                          '--no-lane-splitting'])
        self.assertEqual(bcl2fastq2(
            '/runs/150107_NB123000_0001_ABCX',
            'SampleSheet.csv',
            bcl2fastq_exe='/opt/bin/bcl2fastq').command_line,
                         ['/opt/bin/bcl2fastq',
                          '--runfolder-dir','/runs/150107_NB123000_0001_ABCX',
                          '--output-dir','Unaligned',
                          '--sample-sheet','SampleSheet.csv'])


class TestBclConvert(unittest.TestCase):

    def test_bclconvert(self):
        """Construct 'bcl-convert' command lines for BCL Convert v3.*
        """
        self.assertEqual(bclconvert(
            '/runs/150107_NB123000_0001_ABCX',
            '/output/bclconvert').command_line,
                         ['bcl-convert',
                          '--bcl-input-directory',
                          '/runs/150107_NB123000_0001_ABCX',
                          '--output-dir','/output/bclconvert'])
        self.assertEqual(bclconvert(
            '/runs/150107_NB123000_0001_ABCX',
            '/output/bclconvert',
            sample_sheet='SampleSheet.csv').command_line,
                         ['bcl-convert',
                          '--bcl-input-directory',
                          '/runs/150107_NB123000_0001_ABCX',
                          '--output-dir','/output/bclconvert',
                          '--sample-sheet','SampleSheet.csv'])
        self.assertEqual(bclconvert(
            '/runs/150107_NB123000_0001_ABCX',
            '/output/bclconvert',
            no_lane_splitting=True).command_line,
                         ['bcl-convert',
                          '--bcl-input-directory',
                          '/runs/150107_NB123000_0001_ABCX',
                          '--output-dir','/output/bclconvert',
                          '--no-lane-splitting','true'])
        self.assertEqual(bclconvert(
            '/runs/150107_NB123000_0001_ABCX',
            '/output/bclconvert',
            sampleproject_subdirectories=True).command_line,
                         ['bcl-convert',
                          '--bcl-input-directory',
                          '/runs/150107_NB123000_0001_ABCX',
                          '--output-dir','/output/bclconvert',
                          '--bcl-sampleproject-subdirectories','true'])