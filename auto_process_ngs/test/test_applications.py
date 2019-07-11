#######################################################################
# Tests for applications.py module
#######################################################################
from auto_process_ngs.applications import *
import unittest
import cStringIO

class TestCommand(unittest.TestCase):

    def test_command_simple(self):
        """Construct command in single call
        """
        cmd = Command('ls','-l','-d')
        self.assertEqual(cmd.command,'ls')
        self.assertEqual(cmd.args,['-l','-d'])
        self.assertEqual(cmd.command_line,['ls','-l','-d'])
        self.assertEqual(str(cmd),'ls -l -d')

    def test_command_extended(self):
        """Construct command in multiple calls
        """
        cmd = Command('ls')
        cmd.add_args('-l','-d')
        self.assertEqual(cmd.command,'ls')
        self.assertEqual(cmd.args,['-l','-d'])
        self.assertEqual(cmd.command_line,['ls','-l','-d'])
        self.assertEqual(str(cmd),'ls -l -d')

    def test_command_handles_forward_slash(self):
        """Check forward slash is handled correctly
        """
        cmd = Command('find','-name','"*"','-exec','grep','"hello"','{}','\;')
        self.assertEqual(cmd.command,'find')
        self.assertEqual(cmd.args,['-name','"*"','-exec','grep','"hello"','{}','\;'])
        self.assertEqual(cmd.command_line,['find','-name','"*"','-exec',
                                           'grep','"hello"','{}','\;'])
        self.assertEqual(str(cmd),'find -name "*" -exec grep "hello" {} \;')

    def test_run_subprocess(self):
        """Run command using subprocess
        """
        cmd = Command('ls','-l')
        cmd.run_subprocess(log='/dev/null')

    def test_subprocess_check_output(self):
        """Run command and capture output
        """
        cmd = Command('echo','hello')
        rc,out = cmd.subprocess_check_output()
        self.assertEqual(rc,0)
        self.assertEqual(out,"hello\n")

    def test_has_exe(self):
        """Check 'has_exe' property works
        """
        cmd = Command('ls','-l')
        self.assertTrue(cmd.has_exe)
        cmd = Command('cannot_possibly_exist_3')
        self.assertFalse(cmd.has_exe)

    def test_make_wrapper_script(self):
        """Check 'make_wrapper_script' method works
        """
        cmd = Command('echo','hello')
        self.assertEqual(cmd.make_wrapper_script(),
                         "echo hello")
        self.assertEqual(cmd.make_wrapper_script(shell="/bin/bash"),
                         "#!/bin/bash\necho hello")
        self.assertEqual(
            cmd.make_wrapper_script(prologue="# Prints hello"),
            "# Prints hello\necho hello")
        self.assertEqual(
            cmd.make_wrapper_script(epilogue="# End of script"),
            "echo hello\n# End of script")
        self.assertEqual(
            cmd.make_wrapper_script(
                shell="/bin/bash",
                prologue="echo \"# $(hostname)\"",
                epilogue="echo \"# $(date)\""),
            "#!/bin/bash\necho \"# $(hostname)\"\n"
            "echo hello\necho \"# $(date)\"")
        fp = cStringIO.StringIO()
        cmd.make_wrapper_script(fp=fp)
        self.assertEqual(fp.getvalue(),"echo hello")

    def test_make_wrapper_script_handle_args_with_spaces(self):
        """Check 'make_wrapper_script' method handles arguments with spaces
        """
        cmd = Command('echo','hello world')
        self.assertEqual(cmd.make_wrapper_script(),
                         "echo hello world")
        self.assertEqual(cmd.make_wrapper_script(quote_spaces=True),
                         "echo 'hello world'")

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

