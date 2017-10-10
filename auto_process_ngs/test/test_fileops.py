#######################################################################
# Tests for fileops.py module
#######################################################################

import unittest
import os
import tempfile
import shutil
import pwd
import grp
from bcftbx.JobRunner import SimpleJobRunner
from auto_process_ngs.simple_scheduler import SimpleScheduler
from auto_process_ngs.utils import ZipArchive
from auto_process_ngs.fileops import *

class TestLocation(unittest.TestCase):
    """Tests for the Location class
    """
    def test_location_user_server_path(self):
        """fileops.Location: handle user@host.name:/path/to/somewhere
        """
        location = Location('user@host.name:/path/to/somewhere')
        self.assertEqual(location.user,'user')
        self.assertEqual(location.server,'host.name')
        self.assertEqual(location.path,'/path/to/somewhere')
        self.assertTrue(location.is_remote)
        self.assertEqual(str(location),'user@host.name:/path/to/somewhere')
    def test_location_server_path(self):
        """fileops.Location: handle host.name:/path/to/somewhere
        """
        location = Location('host.name:/path/to/somewhere')
        self.assertEqual(location.user,None)
        self.assertEqual(location.server,'host.name')
        self.assertEqual(location.path,'/path/to/somewhere')
        self.assertTrue(location.is_remote)
        self.assertEqual(str(location),'host.name:/path/to/somewhere')
    def test_location_path(self):
        """fileops.Location: handle /path/to/somewhere
        """
        location = Location('/path/to/somewhere')
        self.assertEqual(location.user,None)
        self.assertEqual(location.server,None)
        self.assertEqual(location.path,'/path/to/somewhere')
        self.assertFalse(location.is_remote)
        self.assertEqual(str(location),'/path/to/somewhere')

class FileopsTestCase(unittest.TestCase):
    """Base class for fileops test cases

    Provides setUp, tearDown and _get_scheduler
    methods.
    """
    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
    def tearDown(self):
        if os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)

class TestMkdirFunction(FileopsTestCase):
    """Tests for the 'mkdir' function
    """
    def test_local_mkdir(self):
        """fileops.mkdir: make a local directory
        """
        new_dir = os.path.join(self.test_dir,'new_dir')
        self.assertFalse(os.path.exists(new_dir))
        status = mkdir(new_dir)
        self.assertEqual(status,0)
        self.assertTrue(os.path.isdir(new_dir))

class TestCopyFunction(FileopsTestCase):
    """Tests for the 'copy' function
    """
    def test_local_copy(self):
        """fileops.copy: copy a local file
        """
        src_file = os.path.join(self.test_dir,'test.txt')
        with open(src_file,'w') as fp:
            fp.write("This is a test file")
        self.assertTrue(os.path.isfile(src_file))
        tgt_file = os.path.join(self.test_dir,'test2.txt')
        self.assertFalse(os.path.exists(tgt_file))
        status = copy(src_file,tgt_file)
        self.assertEqual(status,0)
        self.assertTrue(os.path.isfile(tgt_file))
        self.assertEqual(open(tgt_file).read(),
                         "This is a test file")

class TestSetGroupFunction(FileopsTestCase):
    """Tests for the 'set_group' function
    """
    def test_local_set_group(self):
        """fileops.set_group: set group ownership on a local file
        """
        # Get a list of groups
        current_user = pwd.getpwuid(os.getuid()).pw_name
        groups = [g.gr_gid for g in grp.getgrall()
                  if current_user in g.gr_mem]
        print "Available groups: %s" % groups
        if len(groups) < 2:
            raise unittest.SkipTest("user '%s' must be in at least "
                                    "two groups" % current_user)
        # Make test file and get group
        test_file = os.path.join(self.test_dir,'test.txt')
        with open(test_file,'w') as fp:
            fp.write("This is a test file")
        self.assertTrue(os.path.isfile(test_file))
        gid = os.stat(test_file).st_gid
        print "File: %s GID: %s" % (test_file,gid)
        # Get a second group
        new_gid = None
        for group in groups:
            if group != gid:
                new_gid = group
                break
        print "Selected new GID: %s" % new_gid
        self.assertNotEqual(os.stat(test_file).st_gid,new_gid)
        # Change group by name
        new_group = grp.getgrgid(new_gid).gr_name
        print "New group name: %s" % new_group
        status = set_group(new_group,test_file)
        self.assertEqual(status,0)
        print "File: %s GID: %s"  % (test_file,
                                     os.stat(test_file).st_gid)
        self.assertEqual(os.stat(test_file).st_gid,new_gid)

class TestUnzip(FileopsTestCase):
    """Tests for the 'unzip' function
    """
    def test_local_unzip(self):
        """fileops.unzip: unzip a local file
        """
        test_file = os.path.join(self.test_dir,'test.txt')
        with open(test_file,'w') as fp:
            fp.write("This is a test file")
        zip_file = os.path.join(self.test_dir,'test.zip')
        z = ZipArchive(zip_file,
                       contents=(test_file,),
                       relpath=self.test_dir,
                       prefix='files')
        z.close()
        self.assertTrue(os.path.isfile(zip_file))
        out_dir = os.path.join(self.test_dir,'files')
        self.assertFalse(os.path.exists(out_dir))
        status = unzip(zip_file,self.test_dir)
        self.assertEqual(status,0)
        self.assertTrue(os.path.isdir(out_dir))
        out_file = os.path.join(out_dir,'test.txt')
        self.assertTrue(os.path.isfile(out_file))
        self.assertEqual(open(out_file).read(),
                         "This is a test file")

class TestCopyCommand(unittest.TestCase):
    """Tests for the 'copy_command' function
    """
    def test_copy_command_to_local(self):
        """fileops.copy_command: copy file locally
        """
        copy_cmd = copy_command("/here/myfile.txt",
                                "/there/myfile.txt")
        self.assertEqual(copy_cmd.command_line,
                         ['cp',
                          '/here/myfile.txt',
                          '/there/myfile.txt'])
    def test_copy_command_to_remote(self):
        """fileops.copy_command: copy file to remote
        """
        copy_cmd = copy_command("/here/myfile.txt",
                                "pjx@remote.com:/there/myfile.txt")
        self.assertEqual(copy_cmd.command_line,
                         ['scp',
                          '/here/myfile.txt',
                          'pjx@remote.com:/there/myfile.txt'])

class TestSetGroupCommand(unittest.TestCase):
    """Tests for the 'set_group_command' function
    """
    def test_set_group_command_local(self):
        """fileops.set_group_command: set group on local files
        """
        set_group_cmd = set_group_command("adm",
                                          "/here/files")
        self.assertEqual(set_group_cmd.command_line,
                         ['chgrp',
                          '-R',
                          'adm',
                          '/here/files'])
    def test_set_group_command_remote(self):
        """fileops.set_group_command: set group on remote files
        """
        set_group_cmd = set_group_command("adm",
                                          "pjx@remote.com:/there/files")
        self.assertEqual(set_group_cmd.command_line,
                         ['ssh',
                          'pjx@remote.com',
                          'chgrp',
                          '-R',
                          'adm',
                          '/there/files'])

class TestUnzipCommand(unittest.TestCase):
    """Tests for the 'unzip_command' function
    """
    def test_unzip_command_local(self):
        """fileops.unzip_command: unzip local ZIP archive
        """
        unzip_cmd = unzip_command("/here/myarchive.zip",
                                  "/here/unpacked")
        self.assertEqual(unzip_cmd.command_line,
                         ['unzip',
                          '-q',
                          '-o',
                          '-d',
                          '/here/unpacked',
                          '/here/myarchive.zip'])
    def test_unzip_command_remote(self):
        """fileops.unzip_command: unzip remote ZIP archive
        """
        unzip_cmd = unzip_command("pjx@remote.com:/there/myarchive.zip",
                                  "/there/unpacked")
        self.assertEqual(unzip_cmd.command_line,
                         ['ssh',
                          'pjx@remote.com',
                          'unzip',
                          '-q',
                          '-o',
                          '-d',
                          '/there/unpacked',
                          '/there/myarchive.zip'])
