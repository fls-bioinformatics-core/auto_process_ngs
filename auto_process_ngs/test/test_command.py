#######################################################################
# Tests for command.py module
#######################################################################
from auto_process_ngs.command import Command
import unittest
try:
    from cStringIO import StringIO
except ImportError:
    # cStringIO not available in Python3
    from io import StringIO

class TestCommand(unittest.TestCase):

    def test_command_simple(self):
        """Command: construct command in single call
        """
        cmd = Command('ls','-l','-d')
        self.assertEqual(cmd.command,'ls')
        self.assertEqual(cmd.args,['-l','-d'])
        self.assertEqual(cmd.command_line,['ls','-l','-d'])
        self.assertEqual(str(cmd),'ls -l -d')

    def test_command_extended(self):
        """Command: construct command in multiple calls
        """
        cmd = Command('ls')
        cmd.add_args('-l','-d')
        self.assertEqual(cmd.command,'ls')
        self.assertEqual(cmd.args,['-l','-d'])
        self.assertEqual(cmd.command_line,['ls','-l','-d'])
        self.assertEqual(str(cmd),'ls -l -d')

    def test_command_handles_forward_slash(self):
        """Command: check forward slash is handled correctly
        """
        cmd = Command('find','-name','"*"','-exec','grep','"hello"','{}','\\;')
        self.assertEqual(cmd.command,'find')
        self.assertEqual(cmd.args,['-name','"*"','-exec','grep','"hello"','{}','\\;'])
        self.assertEqual(cmd.command_line,['find','-name','"*"','-exec',
                                           'grep','"hello"','{}','\\;'])
        self.assertEqual(str(cmd),'find -name "*" -exec grep "hello" {} \\;')

    def test_run_subprocess(self):
        """Command: run command using subprocess
        """
        cmd = Command('ls','-l')
        cmd.run_subprocess(log='/dev/null')

    def test_subprocess_check_output(self):
        """Command: run command and capture output
        """
        cmd = Command('echo','hello')
        rc,out = cmd.subprocess_check_output()
        self.assertEqual(rc,0)
        self.assertEqual(out,"hello\n")

    def test_has_exe(self):
        """Command: check 'has_exe' property works
        """
        cmd = Command('ls','-l')
        self.assertTrue(cmd.has_exe)
        cmd = Command('cannot_possibly_exist_3')
        self.assertFalse(cmd.has_exe)

    def test_make_wrapper_script(self):
        """Command: check 'make_wrapper_script' method works
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
        fp = StringIO()
        cmd.make_wrapper_script(fp=fp)
        self.assertEqual(fp.getvalue(),"echo hello")

    def test_make_wrapper_script_handle_args_with_spaces(self):
        """Command: check 'make_wrapper_script' method handles arguments with spaces
        """
        cmd = Command('echo','hello world')
        self.assertEqual(cmd.make_wrapper_script(),
                         "echo hello world")
        self.assertEqual(cmd.make_wrapper_script(quote_spaces=True),
                         "echo 'hello world'")
