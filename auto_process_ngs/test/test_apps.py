#######################################################################
# Tests for apps.py module
#######################################################################
from auto_process_ngs.apps import general
import unittest

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

    def test_ssh_cmd_no_user(self):
        """Construct 'ssh' command lines with no remote user
        """
        self.assertEqual(general.ssh_command(None,'example.com',('ls','-l')).command_line,
                         ['ssh','example.com','ls','-l'])

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

    def test_scp_no_user(self):
        """Construct 'scp' command lines with no remote user
        """
        self.assertEqual(
            general.scp(None,'example.com','my_file','remotedir').command_line,
            ['scp','my_file','example.com:remotedir'])
