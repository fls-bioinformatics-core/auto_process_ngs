#######################################################################
# Tests for applications.py module
#######################################################################
from auto_process_ngs.applications import *
import unittest
try:
    from cStringIO import StringIO
except ImportError:
    # cStringIO not available in Python3
    from io import StringIO

class TestBcl2Fastq(unittest.TestCase):

    def test_configure_bcl_to_fastq(self):
        """Construct 'configureBclToFastq.pl' command lines
        """
        self.assertEqual(bcl2fastq.configureBclToFastq(
            'Data/Intensities/Basecalls',
            'SampleSheet.csv').command_line,
                         ['configureBclToFastq.pl',
                          '--input-dir','Data/Intensities/Basecalls',
                          '--output-dir','Unaligned',
                          '--sample-sheet','SampleSheet.csv',
                          '--fastq-cluster-count','-1'])
        self.assertEqual(bcl2fastq.configureBclToFastq(
            'Data/Intensities/Basecalls',
            'SampleSheet.csv',
            output_dir='run/bcl2fastq').command_line,
                         ['configureBclToFastq.pl',
                          '--input-dir','Data/Intensities/Basecalls',
                          '--output-dir','run/bcl2fastq',
                          '--sample-sheet','SampleSheet.csv',
                          '--fastq-cluster-count','-1'])
        self.assertEqual(bcl2fastq.configureBclToFastq(
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
        self.assertEqual(bcl2fastq.configureBclToFastq(
            'Data/Intensities/Basecalls',
            'SampleSheet.csv',
            configureBclToFastq_exe=
            "/opt/bin/configureBclToFastq.pl").command_line,
                         ['/opt/bin/configureBclToFastq.pl',
                          '--input-dir','Data/Intensities/Basecalls',
                          '--output-dir','Unaligned',
                          '--sample-sheet','SampleSheet.csv',
                          '--fastq-cluster-count','-1'])

    def test_bcl2fastq(self):
        """Construct 'bcl2fastq' command lines for bcl2fastq v2.*
        """
        self.assertEqual(bcl2fastq.bcl2fastq2(
            '/runs/150107_NB123000_0001_ABCX',
            'SampleSheet.csv').command_line,
                         ['bcl2fastq',
                          '--runfolder-dir','/runs/150107_NB123000_0001_ABCX',
                          '--output-dir','Unaligned',
                          '--sample-sheet','SampleSheet.csv'])
        self.assertEqual(bcl2fastq.bcl2fastq2(
            '/runs/150107_NB123000_0001_ABCX',
            'SampleSheet.csv',
            output_dir='run/bcl2fastq').command_line,
                         ['bcl2fastq',
                          '--runfolder-dir','/runs/150107_NB123000_0001_ABCX',
                          '--output-dir','run/bcl2fastq',
                          '--sample-sheet','SampleSheet.csv'])
        self.assertEqual(bcl2fastq.bcl2fastq2(
            '/runs/150107_NB123000_0001_ABCX',
            'SampleSheet.csv',
            output_dir='run/bcl2fastq',
            ignore_missing_bcl=True).command_line,
                         ['bcl2fastq',
                          '--runfolder-dir','/runs/150107_NB123000_0001_ABCX',
                          '--output-dir','run/bcl2fastq',
                          '--sample-sheet','SampleSheet.csv',
                          '--ignore-missing-bcls'])
        self.assertEqual(bcl2fastq.bcl2fastq2(
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
        self.assertEqual(bcl2fastq.bcl2fastq2(
            '/runs/150107_NB123000_0001_ABCX',
            'SampleSheet.csv',
            bcl2fastq_exe='/opt/bin/bcl2fastq').command_line,
                         ['/opt/bin/bcl2fastq',
                          '--runfolder-dir','/runs/150107_NB123000_0001_ABCX',
                          '--output-dir','Unaligned',
                          '--sample-sheet','SampleSheet.csv'])

    def test_bclconvert(self):
        """Construct 'bcl-convert' command lines for BCL Convert v3.*
        """
        self.assertEqual(bcl2fastq.bclconvert(
            '/runs/150107_NB123000_0001_ABCX',
            '/output/bclconvert').command_line,
                         ['bcl-convert',
                          '--bcl-input-directory',
                          '/runs/150107_NB123000_0001_ABCX',
                          '--output-dir','/output/bclconvert'])
        self.assertEqual(bcl2fastq.bclconvert(
            '/runs/150107_NB123000_0001_ABCX',
            '/output/bclconvert',
            sample_sheet='SampleSheet.csv').command_line,
                         ['bcl-convert',
                          '--bcl-input-directory',
                          '/runs/150107_NB123000_0001_ABCX',
                          '--output-dir','/output/bclconvert',
                          '--sample-sheet','SampleSheet.csv'])
        self.assertEqual(bcl2fastq.bclconvert(
            '/runs/150107_NB123000_0001_ABCX',
            '/output/bclconvert',
            no_lane_splitting=True).command_line,
                         ['bcl-convert',
                          '--bcl-input-directory',
                          '/runs/150107_NB123000_0001_ABCX',
                          '--output-dir','/output/bclconvert',
                          '--no-lane-splitting','true'])
        self.assertEqual(bcl2fastq.bclconvert(
            '/runs/150107_NB123000_0001_ABCX',
            '/output/bclconvert',
            sampleproject_subdirectories=True).command_line,
                         ['bcl-convert',
                          '--bcl-input-directory',
                          '/runs/150107_NB123000_0001_ABCX',
                          '--output-dir','/output/bclconvert',
                          '--bcl-sampleproject-subdirectories','true'])

class TestGeneral(unittest.TestCase):

    def test_rsync(self):
        """Construct 'rsync' command lines
        """
        self.assertEqual(general.rsync('from','to').command_line,
                         ['rsync','-av','from','to'])
        self.assertEqual(general.rsync('from','user@remote.com:to').command_line,
                         ['rsync','-av','-e','ssh','from','user@remote.com:to'])

    def test_make(self):
        """Construct 'make' command lines
        """
        self.assertEqual(general.make(makefile='makefile',
                                      working_dir='Unaligned',
                                      nprocessors='12').command_line,
                         ['make',
                          '-C','Unaligned',
                          '-j','12',
                          '-f','makefile'])

    def test_ssh_cmd(self):
        """Construct 'ssh' command lines
        """
        self.assertEqual(general.ssh_command('user','example.com',('ls','-l')).command_line,
                         ['ssh','user@example.com','ls','-l'])

    def test_scp(self):
        """Construct 'scp' command lines
        """
        self.assertEqual(
            general.scp('user','example.com','my_file','remotedir').command_line,
            ['scp','my_file','user@example.com:remotedir'])

    def test_scp_recursive(self):
        """Construct 'scp -r' command lines
        """
        self.assertEqual(
            general.scp('user','example.com','my_dir','remotedir',
                        recursive=True).command_line,
            ['scp','-r','my_dir','user@example.com:remotedir'])

