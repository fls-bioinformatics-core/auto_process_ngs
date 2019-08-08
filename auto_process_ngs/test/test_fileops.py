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
        self.assertFalse(location.is_url)
        self.assertEqual(location.url,None)
        self.assertEqual(location.protocol,None)
        self.assertEqual(str(location),'user@host.name:/path/to/somewhere')
    def test_location_server_path(self):
        """fileops.Location: handle host.name:/path/to/somewhere
        """
        location = Location('host.name:/path/to/somewhere')
        self.assertEqual(location.user,None)
        self.assertEqual(location.server,'host.name')
        self.assertEqual(location.path,'/path/to/somewhere')
        self.assertTrue(location.is_remote)
        self.assertFalse(location.is_url)
        self.assertEqual(location.url,None)
        self.assertEqual(location.protocol,None)
        self.assertEqual(str(location),'host.name:/path/to/somewhere')
    def test_location_path(self):
        """fileops.Location: handle /path/to/somewhere
        """
        location = Location('/path/to/somewhere')
        self.assertEqual(location.user,None)
        self.assertEqual(location.server,None)
        self.assertEqual(location.path,'/path/to/somewhere')
        self.assertFalse(location.is_remote)
        self.assertFalse(location.is_url)
        self.assertEqual(location.url,None)
        self.assertEqual(location.protocol,None)
        self.assertEqual(str(location),'/path/to/somewhere')
    def test_location_http_url(self):
        """fileops.Location: handle http://example.com/path/to/url
        """
        location = Location('http://example.com/path/to/url')
        self.assertEqual(location.user,None)
        self.assertEqual(location.server,'example.com')
        self.assertEqual(location.path,'path/to/url')
        self.assertFalse(location.is_remote)
        self.assertTrue(location.is_url)
        self.assertEqual(location.url,'http://example.com/path/to/url')
        self.assertEqual(location.protocol,'http')
        self.assertEqual(str(location),'http://example.com/path/to/url')
    def test_location_https_url(self):
        """fileops.Location: handle https://example.com/path/to/url
        """
        location = Location('https://example.com/path/to/url')
        self.assertEqual(location.user,None)
        self.assertEqual(location.server,'example.com')
        self.assertEqual(location.path,'path/to/url')
        self.assertFalse(location.is_remote)
        self.assertTrue(location.is_url)
        self.assertEqual(location.protocol,'https')
        self.assertEqual(location.url,'https://example.com/path/to/url')
        self.assertEqual(str(location),'https://example.com/path/to/url')
    def test_location_file_url(self):
        """fileops.Location: handle file:///path/to/file
        """
        location = Location('file:///path/to/file')
        self.assertEqual(location.user,None)
        self.assertEqual(location.server,None)
        self.assertEqual(location.path,'/path/to/file')
        self.assertFalse(location.is_remote)
        self.assertTrue(location.is_url)
        self.assertEqual(location.protocol,'file')
        self.assertEqual(location.url,'file:///path/to/file')
        self.assertEqual(str(location),'file:///path/to/file')

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

    def test_local_mkdir_missing_parent(self):
        """fileops.mkdir: fails to make local dir if parent is missing
        """
        new_dir = os.path.join(self.test_dir,'new_dir','sub_dir')
        self.assertFalse(os.path.exists(new_dir))
        status = mkdir(new_dir)
        self.assertEqual(status,1)
        self.assertFalse(os.path.exists(new_dir))

    def test_local_mkdir_recursive(self):
        """fileops.mkdir: make a local directory recursively
        """
        new_dir = os.path.join(self.test_dir,'new_dir','sub_dir')
        self.assertFalse(os.path.exists(new_dir))
        status = mkdir(new_dir,recursive=True)
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

class TestCopytreeFunction(FileopsTestCase):
    """Tests for the 'copy' function
    """
    def test_local_copytree(self):
        """fileops.copytree: recursively copy a local directory
        """
        # Build a source directory
        src_dir = os.path.join(self.test_dir,'test_src')
        os.mkdir(src_dir)
        os.mkdir(os.path.join(src_dir,'subdir'))
        with open(os.path.join(src_dir,'test1.txt'),'w') as fp:
            fp.write("This is test file 1")
        with open(os.path.join(src_dir,'subdir','test2.txt'),'w') as fp:
            fp.write("This is test file 2")
        # Copy contents recursively
        tgt_dir = os.path.join(self.test_dir,'test_dst')
        self.assertFalse(os.path.exists(tgt_dir))
        status = copytree(src_dir,tgt_dir)
        self.assertEqual(status,0)
        self.assertTrue(os.path.isdir(tgt_dir))
        self.assertTrue(os.path.isdir(os.path.join(tgt_dir,'subdir')))
        self.assertEqual(open(os.path.join(tgt_dir,'test1.txt')).read(),
                         "This is test file 1")
        self.assertEqual(open(os.path.join(tgt_dir,'subdir','test2.txt')).read(),
                         "This is test file 2")

class TestSetGroupFunction(FileopsTestCase):
    """Tests for the 'set_group' function
    """
    def test_local_set_group(self):
        """fileops.set_group: set group ownership on a local directory
        """
        # Get a list of groups
        current_user = pwd.getpwuid(os.getuid()).pw_name
        groups = [g.gr_gid for g in grp.getgrall()
                  if current_user in g.gr_mem]
        print("Available groups: %s" % groups)
        if len(groups) < 2:
            raise unittest.SkipTest("user '%s' must be in at least "
                                    "two groups" % current_user)
        # Make test directory and files
        test_dir = os.path.join(self.test_dir,'test.dir')
        os.mkdir(test_dir)
        test_file1 = os.path.join(test_dir,'test1.txt')
        test_file2 = os.path.join(test_dir,'test2.txt')
        test_file3 = os.path.join(test_dir,'test3.txt')
        for f in (test_file1,test_file2,test_file3):
            with open(f,'w') as fp:
                fp.write("This is a test file")
            self.assertTrue(os.path.isfile(f))
        # Get initial group from a test file
        gid = os.stat(test_file1).st_gid
        print("File: %s GID: %s" % (test_file1,gid))
        # Get a second group
        new_gid = None
        for group in groups:
            if group != gid:
                new_gid = group
                break
        print("Selected new GID: %s" % new_gid)
        for f in (test_file1,test_file2,test_file3):
            self.assertNotEqual(os.stat(f).st_gid,new_gid)
        # Change group by name
        new_group = grp.getgrgid(new_gid).gr_name
        print("New group name: %s" % new_group)
        status = set_group(new_group,test_dir)
        self.assertEqual(status,0)
        for f in (test_file1,test_file2,test_file3):
            print("File: %s GID: %s"  % (f,os.stat(f).st_gid))
            self.assertEqual(os.stat(f).st_gid,new_gid)

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

class TestRename(FileopsTestCase):
    """Tests for the 'rename' function
    """
    def test_local_rename_file(self):
        """fileops.rename: rename a local file
        """
        old_file = os.path.join(self.test_dir,'test.txt')
        new_file = os.path.join(self.test_dir,'test2.txt')
        with open(old_file,'w') as fp:
            fp.write("This is a test file")
        status = rename(old_file,new_file)
        self.assertEqual(status,0)
        self.assertFalse(os.path.exists(old_file))
        self.assertTrue(os.path.exists(new_file))

    def test_local_rename_dir(self):
        """fileops.rename: rename a local directory
        """
        old_dir = os.path.join(self.test_dir,'test_dir')
        new_dir = os.path.join(self.test_dir,'test_dir2')
        os.mkdir(old_dir)
        status = rename(old_dir,new_dir)
        self.assertEqual(status,0)
        self.assertFalse(os.path.exists(old_dir))
        self.assertTrue(os.path.exists(new_dir))

class TestExists(FileopsTestCase):
    """Tests for the 'exists' function
    """
    def test_local_exists(self):
        """fileops.exists: test for a local file/directory
        """
        test_file = os.path.join(self.test_dir,'test.txt')
        with open(test_file,'w') as fp:
            fp.write("This is a test file")
        test_dir = os.path.join(self.test_dir,'test_dir')
        os.mkdir(test_dir)
        missing = os.path.join(self.test_dir,'missing')
        self.assertTrue(exists(test_file))
        self.assertTrue(exists(test_dir))
        self.assertFalse(exists(missing))

class TestRemoveFile(FileopsTestCase):
    """Tests for the 'remove_file' function
    """
    def test_local_remove_file(self):
        """fileops.remove_file: remove a local file
        """
        test_file = os.path.join(self.test_dir,'test.txt')
        with open(test_file,'w') as fp:
            fp.write("This is a test file")
        self.assertTrue(os.path.exists(test_file))
        status = remove_file(test_file)
        self.assertEqual(status,0)
        self.assertFalse(os.path.exists(test_file))
    def test_local_remove_file_fails_on_dir(self):
        """fileops.remove_file: cannot remove a local directory
        """
        test_dir = os.path.join(self.test_dir,'test_dir')
        os.mkdir(test_dir)
        self.assertTrue(os.path.exists(test_dir))
        status = remove_file(test_dir)
        self.assertEqual(status,1)
        self.assertTrue(os.path.exists(test_dir))

class TestDiskUsage(FileopsTestCase):
    """Tests for the 'disk_usage' function
    """
    def test_local_disk_usage(self):
        """fileops.disk_usage: get disk usage for a local directory
        """
        test_dir = os.path.join(self.test_dir,'test_dir')
        os.mkdir(test_dir)
        self.assertTrue(os.path.exists(test_dir))
        usage = disk_usage(test_dir)
        self.assertTrue(usage.total > 0)
        self.assertTrue(usage.used > 0)
        self.assertTrue(usage.free > 0)
        self.assertTrue(usage.percent > 0)
    def test_local_remove_file_fails_on_dir(self):
        """fileops.disk_usage: raise OSError when directory doesn't exist
        """
        test_dir = os.path.join(self.test_dir,'test_dir')
        self.assertFalse(os.path.exists(test_dir))
        self.assertRaises(OSError,
                          disk_usage,
                          test_dir)

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

class TestCopytreeCommand(unittest.TestCase):
    """Tests for the 'copytree_command' function
    """
    def test_copytree_command_to_local(self):
        """fileops.copytree_command: copy directory locally
        """
        copytree_cmd = copytree_command("/here",
                                        "/there")
        self.assertEqual(copytree_cmd.command_line,
                         ['cp',
                          '-a',
                          '-n',
                          '/here',
                          '/there'])
    def test_copytree_command_to_remote(self):
        """fileops.copytree_command: copy directory to remote
        """
        copytree_cmd = copytree_command("/here",
                                        "pjx@remote.com:/there")
        self.assertEqual(copytree_cmd.command_line,
                         ['scp',
                          '-r',
                          '/here',
                          'pjx@remote.com:/there'])

class TestSetGroupCommand(unittest.TestCase):
    """Tests for the 'set_group_command' function
    """
    def test_set_group_command_local(self):
        """fileops.set_group_command: set group on local files
        """
        set_group_cmd = set_group_command("adm",
                                          "/here/files")
        self.assertEqual(set_group_cmd.command_line,
                         ['find',
                          '/here/files',
                          '-exec',
                          'chgrp',
                          'adm',
                          '{}',
                          '+'])
    def test_set_group_command_local_verbose(self):
        """fileops.set_group_command: set group on local files (verbose)
        """
        set_group_cmd = set_group_command("adm",
                                          "/here/files",
                                          verbose=True)
        self.assertEqual(set_group_cmd.command_line,
                         ['find',
                          '/here/files',
                          '-exec',
                          'chgrp',
                          '--verbose',
                          'adm',
                          '{}',
                          '+'])

    def test_set_group_command_local_safe(self):
        """fileops.set_group_command: set group on local files ('safe' mode)
        """
        set_group_cmd = set_group_command("adm",
                                          "/here/files",
                                          safe=True)
        self.assertEqual(set_group_cmd.command_line,
                         ['find',
                          '/here/files',
                          '-user',
                          getpass.getuser(),
                          '-exec',
                          'chgrp',
                          'adm',
                          '{}',
                          '+'])

    def test_set_group_command_remote(self):
        """fileops.set_group_command: set group on remote files
        """
        set_group_cmd = set_group_command("adm",
                                          "pjx@remote.com:/there/files")
        self.assertEqual(set_group_cmd.command_line,
                         ['ssh',
                          'pjx@remote.com',
                          'find',
                          '/there/files',
                          '-exec',
                          'chgrp',
                          'adm',
                          '{}',
                          '+'])

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
