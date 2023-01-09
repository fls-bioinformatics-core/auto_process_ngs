#######################################################################
# Tests for utils.py module
#######################################################################

import unittest
import tempfile
import shutil
import gzip
import zipfile
from bcftbx.JobRunner import SimpleJobRunner,GEJobRunner
from bcftbx.utils import find_program
from auto_process_ngs.command import Command
from auto_process_ngs.utils import *

class TestZipArchive(unittest.TestCase):
    """
    Tests for the ZipArchive class
    """
    def setUp(self):
        self.dirn = tempfile.mkdtemp(suffix='TestZipArchive')

    def tearDown(self):
        # Remove the temporary test directory
        shutil.rmtree(self.dirn)

    def test_create_zip_archive(self):
        """ZipArchive: create a new zip archive
        """
        # Make an example directory to zip up
        src_dir = os.path.join(self.dirn,'source')
        items = ('test1',
                 'test2',
                 'sub1/',
                 'sub1/test3',
                 'sub1/test4',
                 'sub1/sub2/',
                 'sub1/sub2/test5',
                 'sub3/',
                 'sub3/test6')
        os.mkdir(src_dir)
        for item in [os.path.join(src_dir,x) for x in items]:
            if item.endswith('/'):
                try:
                    os.mkdir(item)
                except OSError as ex:
                    raise OSError("Failed to make %s" % item)
            else:
                with open(item,'wt') as fp:
                    fp.write('')
        # Create the zip archive
        zip_filename = os.path.join(self.dirn,'test.zip')
        self.assertFalse(os.path.exists(zip_filename))
        z = ZipArchive(zip_filename,relpath=src_dir)
        for item in [os.path.join(src_dir,x) for x in ('test1','sub1')]:
            z.add(item)
        z.close()
        # Check the zip archive exists
        self.assertTrue(os.path.exists(zip_filename),
                        "zip file not created")
        self.assertTrue(zipfile.is_zipfile(zip_filename),
                        "file exists but is not zip file")
        # Check the contents
        namelist = zipfile.ZipFile(zip_filename).namelist()
        print(namelist)
        for item in ('test1',
                     'sub1/test3',
                     'sub1/test4',
                     'sub1/sub2/test5'):
            self.assertTrue((item.rstrip('/') in namelist),
                            "Item '%s' not in zip file" % item)
        for name in namelist:
            self.assertTrue((name in ('test1',
                                      'sub1/test3',
                                      'sub1/test4',
                                      'sub1/sub2/test5')),
                            "'%s' in zip file but shouldn't be"
                            % name)

    def test_create_zip_archive_fails_for_files_outside_relpath(self):
        """ZipArchive: create a new zip archive fails for item outside relpath
        """
        # Make an example directory to zip up
        src_dir = os.path.join(self.dirn,'source')
        items = ('test1',
                 'sub1/',
                 'sub1/test2',)
        os.mkdir(src_dir)
        for item in [os.path.join(src_dir,x) for x in items]:
            if item.endswith('/'):
                try:
                    os.mkdir(item)
                except OSError as ex:
                    raise OSError("Failed to make %s" % item)
            else:
                with open(item,'wt') as fp:
                    fp.write('')
        # Make a file outside the directory
        extra_file = os.path.join(self.dirn,'extra')
        with open(extra_file,'wt') as fp:
            fp.write('')
        # Create the zip archive
        zip_filename = os.path.join(self.dirn,'test.zip')
        self.assertFalse(os.path.exists(zip_filename))
        z = ZipArchive(zip_filename,relpath=src_dir)
        for item in [os.path.join(src_dir,x) for x in ('test1','sub1')]:
            z.add(item)
        # Add file outside the relpath
        self.assertRaises(Exception,
                          z.add,
                          extra_file)

class TestOutputFiles(unittest.TestCase):
    """
    Tests for the OutputFiles class
    """
    def setUp(self):
        # Temporary working dir
        self.wd = tempfile.mkdtemp(suffix='.test_outputfile')
    def tearDown(self):
        # Remove temporary working dir
        if os.path.isdir(self.wd):
            shutil.rmtree(self.wd)
    def test_outputfiles(self):
        out = OutputFiles()
        out.open('test1',os.path.join(self.wd,'test1.txt'))
        out.open('test2',os.path.join(self.wd,'test2.txt'))
        self.assertEqual(
            out.file_name('test1'),os.path.join(self.wd,'test1.txt'))
        self.assertEqual(
            out.file_name('test2'),os.path.join(self.wd,'test2.txt'))
        out.write('test1','Some test text')
        out.write('test2','Some more\ntest text')
        out.close()
        with open(out.file_name('test1'),'rt') as fp:
            self.assertEqual(fp.read(),"Some test text\n")
        with open(out.file_name('test2'),'rt') as fp:
            self.assertEqual(fp.read(),"Some more\ntest text\n")
    def test_outputfiles_with_basedir(self):
        out = OutputFiles(self.wd)
        out.open('test1','test1.txt')
        out.open('test2','test2.txt')
        out.close()
        self.assertEqual(
            out.file_name('test1'),os.path.join(self.wd,'test1.txt'))
        self.assertEqual(
            out.file_name('test2'),os.path.join(self.wd,'test2.txt'))
    def test_outputfiles_contains(self):
        out = OutputFiles()
        out.open('test1',os.path.join(self.wd,'test1.txt'))
        self.assertTrue('test1' in out)
        self.assertFalse('test2' in out)
        out.close()
        self.assertFalse('test1' in out)
        self.assertEqual(out.file_name('test1'),
                         os.path.join(self.wd,'test1.txt'))
        self.assertFalse('test2' in out)
        self.assertRaises(KeyError,out.file_name,'test2')
    def test_outputfiles_append(self):
        out = OutputFiles()
        with open(os.path.join(self.wd,'test1.txt'),'w') as fp:
            fp.write("test1 already exists\n")
        with open(os.path.join(self.wd,'test2.txt'),'w') as fp:
            fp.write("test2 already exists\n")
        out.open('test1',os.path.join(self.wd,'test1.txt'),append=True)
        out.open('test2',os.path.join(self.wd,'test2.txt'))
        self.assertEqual(
            out.file_name('test1'),os.path.join(self.wd,'test1.txt'))
        self.assertEqual(
            out.file_name('test2'),os.path.join(self.wd,'test2.txt'))
        out.write('test1','Some test text')
        out.write('test2','Some more\ntest text')
        out.close()
        with open(out.file_name('test1'),'rt') as fp:
            self.assertEqual(fp.read(),
                             "test1 already exists\nSome test text\n")
        with open(out.file_name('test2'),'rt') as fp:
            self.assertEqual(fp.read(),"Some more\ntest text\n")
    def test_outputfiles_reopen(self):
        out = OutputFiles(base_dir=self.wd)
        out.open('test1','test1.txt')
        out.open('test2','test2.txt')
        self.assertEqual(out.file_name('test1'),
                         os.path.join(self.wd,'test1.txt'))
        self.assertEqual(out.file_name('test2'),
                         os.path.join(self.wd,'test2.txt'))
        out.write('test1','Some test text')
        out.write('test2','Some more\ntest text')
        out.close()
        with open(out.file_name('test1'),'rt') as fp:
            self.assertEqual(fp.read(),"Some test text\n")
        with open(out.file_name('test2'),'rt') as fp:
            self.assertEqual(fp.read(),"Some more\ntest text\n")
        # Reopen files
        out.open('test1',append=True)
        out.open('test2')
        self.assertEqual(out.file_name('test1'),
                         os.path.join(self.wd,'test1.txt'))
        self.assertEqual(out.file_name('test2'),
                         os.path.join(self.wd,'test2.txt'))
        out.write('test1','Some extra test text')
        out.write('test2','Some different\ntest text')
        out.close()
        with open(out.file_name('test1'),'rt') as fp:
            self.assertEqual(fp.read(),
                             "Some test text\nSome extra test text\n")
        with open(out.file_name('test2'),'rt') as fp:
            self.assertEqual(fp.read(),"Some different\ntest text\n")
    def test_outputfiles_len(self):
        out = OutputFiles(base_dir=self.wd)
        self.assertEqual(len(out),0)
        out.open('test1','test1.txt')
        self.assertEqual(len(out),1)
        out.close('test1')
        self.assertEqual(len(out),0)
        out.open('test2','test2.txt')
        out.open('test3','test3.txt')
        self.assertEqual(len(out),2)
        out.open('test1',append=True)
        self.assertEqual(len(out),3)
        out.close()
        self.assertEqual(len(out),0)
        out.open('test2',append=True)
        self.assertEqual(len(out),1)
        out.close()

class TestBufferedOutputFiles(unittest.TestCase):
    """
    Tests for the BufferedOutputFiles class
    """
    def setUp(self):
        # Temporary working dir
        self.wd = tempfile.mkdtemp(suffix='.test_bufferedoutputfiles')
    def tearDown(self):
        # Remove temporary working dir
        if os.path.isdir(self.wd):
            shutil.rmtree(self.wd)
    def test_bufferedoutputfiles(self):
        out = BufferedOutputFiles()
        out.open('test1',os.path.join(self.wd,'test1.txt'))
        out.open('test2',os.path.join(self.wd,'test2.txt'))
        self.assertEqual(
            out.file_name('test1'),os.path.join(self.wd,'test1.txt'))
        self.assertEqual(
            out.file_name('test2'),os.path.join(self.wd,'test2.txt'))
        out.write('test1','Some test text')
        out.write('test2','Some more\ntest text')
        out.close()
        with open(out.file_name('test1'),'rt') as fp:
            self.assertEqual(fp.read(),"Some test text\n")
        with open(out.file_name('test2'),'rt') as fp:
            self.assertEqual(fp.read(),"Some more\ntest text\n")
    def test_bufferedoutputfiles_with_gzip(self):
        out = BufferedOutputFiles()
        out.open('test1',os.path.join(self.wd,'test1.txt.gz'))
        out.open('test2',os.path.join(self.wd,'test2.txt.gz'))
        self.assertEqual(
            out.file_name('test1'),os.path.join(self.wd,'test1.txt.gz'))
        self.assertEqual(
            out.file_name('test2'),os.path.join(self.wd,'test2.txt.gz'))
        out.write('test1','Some test text')
        out.write('test2','Some more\ntest text')
        out.close()
        with gzip.open(out.file_name('test1'),'rt') as fp:
            self.assertEqual(fp.read(),"Some test text\n")
        with gzip.open(out.file_name('test2'),'rt') as fp:
            self.assertEqual(fp.read(),"Some more\ntest text\n")
    def test_bufferedoutputfiles_with_basedir(self):
        out = BufferedOutputFiles(self.wd)
        out.open('test1','test1.txt')
        out.open('test2','test2.txt')
        self.assertEqual(
            out.file_name('test1'),os.path.join(self.wd,'test1.txt'))
        self.assertEqual(
            out.file_name('test2'),os.path.join(self.wd,'test2.txt'))
    def test_bufferedoutputfiles_append(self):
        out = BufferedOutputFiles()
        with open(os.path.join(self.wd,'test1.txt'),'wt') as fp:
            fp.write("test1 already exists\n")
        with open(os.path.join(self.wd,'test2.txt'),'wt') as fp:
            fp.write("test2 already exists\n")
        out.open('test1',os.path.join(self.wd,'test1.txt'),append=True)
        out.open('test2',os.path.join(self.wd,'test2.txt'))
        self.assertEqual(
            out.file_name('test1'),os.path.join(self.wd,'test1.txt'))
        self.assertEqual(
            out.file_name('test2'),os.path.join(self.wd,'test2.txt'))
        out.write('test1','Some test text')
        out.write('test2','Some more\ntest text')
        out.close()
        with open(out.file_name('test1'),'rt') as fp:
            self.assertEqual(fp.read(),
                             "test1 already exists\nSome test text\n")
        with open(out.file_name('test2'),'rt') as fp:
            self.assertEqual(fp.read(),"Some more\ntest text\n")
    def test_bufferedoutputfiles_reopen(self):
        out = BufferedOutputFiles(base_dir=self.wd)
        out.open('test1','test1.txt')
        out.open('test2','test2.txt')
        self.assertEqual(out.file_name('test1'),
                         os.path.join(self.wd,'test1.txt'))
        self.assertEqual(out.file_name('test2'),
                         os.path.join(self.wd,'test2.txt'))
        out.write('test1','Some test text')
        out.write('test2','Some more\ntest text')
        out.close()
        with open(out.file_name('test1'),'rt') as fp:
            self.assertEqual(fp.read(),"Some test text\n")
        with open(out.file_name('test2'),'rt') as fp:
            self.assertEqual(fp.read(),"Some more\ntest text\n")
        # Reopen files
        out.open('test1',append=True)
        out.open('test2')
        self.assertEqual(out.file_name('test1'),
                         os.path.join(self.wd,'test1.txt'))
        self.assertEqual(out.file_name('test2'),
                         os.path.join(self.wd,'test2.txt'))
        out.write('test1','Some extra test text')
        out.write('test2','Some different\ntest text')
        out.close()
        with open(out.file_name('test1'),'rt') as fp:
            self.assertEqual(fp.read(),
                             "Some test text\nSome extra test text\n")
        with open(out.file_name('test2'),'rt') as fp:
            self.assertEqual(fp.read(),"Some different\ntest text\n")
    def test_bufferedoutputfiles_exceed_maximum_open_files(self):
        # Only allow 2 open files at once internally
        # (and with an artifically small buffer size)
        out = BufferedOutputFiles(base_dir=self.wd,
                                  bufsize=2,
                                  max_open_files=2)
        # Declare 4 output files
        out.open('test1','test1.txt')
        out.open('test2','test2.txt')
        out.open('test3','test3.txt')
        out.open('test4','test4.txt')
        # Write content to each
        out.write('test1','test text')
        out.write('test2','different\ntest text')
        out.write('test3','more test text')
        out.write('test4','yet more\ntest text')
        out.write('test1','and a bit more')
        out.write('test2','more but different')
        out.write('test3','different again')
        out.write('test4','and even more different\ntest text')
        out.write('test1','plus some final\ntext')
        out.write('test2','and a last bit')
        out.write('test3','plus a coda')
        out.write('test4','et\nFIN')
        # Close everything
        out.close()
        # Check the contents for each file
        with open(out.file_name('test1'),'rt') as fp:
            self.assertEqual(fp.read(),
                             "test text\n"
                             "and a bit more\n"
                             "plus some final\ntext\n")
        with open(out.file_name('test2'),'rt') as fp:
            self.assertEqual(fp.read(),
                             "different\ntest text\n"
                             "more but different\n"
                             "and a last bit\n")
        with open(out.file_name('test3'),'rt') as fp:
            self.assertEqual(fp.read(),
                             "more test text\n"
                             "different again\n"
                             "plus a coda\n")
        with open(out.file_name('test4'),'rt') as fp:
            self.assertEqual(fp.read(),
                             "yet more\ntest text\n"
                             "and even more different\ntest text\n"
                             "et\nFIN\n")

class TestShowProgressChecker(unittest.TestCase):
    """
    Tests for the ProgressChecker class
    """
    def test_check_progress_every(self):
        """
        ProgressChecker: check every N items
        """
        progress = ProgressChecker(every=1)
        self.assertTrue(progress.check(0))
        self.assertTrue(progress.check(1))
        progress = ProgressChecker(every=5)
        self.assertTrue(progress.check(0))
        self.assertFalse(progress.check(1))
        self.assertTrue(progress.check(5))
        self.assertTrue(progress.check(20))
        self.assertFalse(progress.check(21))
    def test_check_progress_every_percent(self):
        """
        ProgressChecker: check every percentage of items
        """
        progress = ProgressChecker(percent=5,total=12209)
        self.assertTrue(progress.check(0))
        self.assertFalse(progress.check(1))
        self.assertFalse(progress.check(609))
        self.assertTrue(progress.check(610))
    def test_check_progress_every_percent_zero_total(self):
        """
        ProgressChecker: supplied total is zero for percentage
        """
        progress = ProgressChecker(percent=5,total=0)
        self.assertFalse(progress.check(0))
    def test_check_progress_every_percent_zero_total(self):
        """
        ProgressChecker: supplied total is smaller than percentage
        """
        progress = ProgressChecker(percent=10,total=9)
        self.assertTrue(progress.check(0))
        self.assertTrue(progress.check(1))
    def test_get_percentage(self):
        """
        ProgressChecker: check returned percentage
        """
        progress = ProgressChecker(percent=5,total=12209)
        self.assertEqual(progress.percent(0),float(0))
        self.assertEqual(progress.percent(610),610.0/12209.0*100.0)
        self.assertEqual(progress.percent(12209),float(100))

class TestFileLock(unittest.TestCase):
    """
    Tests for the FileLock class
    """
    def setUp(self):
        self.wd = tempfile.mkdtemp(suffix='TestFileLock')

    def tearDown(self):
        # Remove the temporary test directory
        shutil.rmtree(self.wd)

    def test_file_lock_on_file(self):
        """
        FileLock: lock a file
        """
        # Make a file to lock
        f = os.path.join(self.wd,"test.txt")
        with open(f,'wt') as fp:
            fp.write("Text\n")
        # Set up file lock
        flock = FileLock(f)
        self.assertFalse(flock.has_lock)
        # Acquire lock
        flock.acquire()
        self.assertTrue(flock.has_lock)
        # Check that a second lock can't be
        # acquired
        flock2 = FileLock(f)
        self.assertFalse(flock2.has_lock)
        self.assertRaises(FileLockError,
                          flock2.acquire)
        self.assertFalse(flock2.has_lock)
        # Release first lock
        flock.release()
        self.assertFalse(flock.has_lock)
        # Check second lock can now be
        # acquired
        flock2.acquire()
        self.assertTrue(flock2.has_lock)
        flock2.release()
        self.assertFalse(flock2.has_lock)

    def test_file_lock_on_directory(self):
        """
        FileLock: lock a directory
        """
        # Make a directory to lock
        d = os.path.join(self.wd,"test")
        os.mkdir(d)
        # Set up directory lock
        flock = FileLock(d)
        self.assertFalse(flock.has_lock)
        # Acquire lock
        flock.acquire()
        self.assertTrue(flock.has_lock)
        # Check that a second lock can't be
        # acquired
        flock2 = FileLock(d)
        self.assertFalse(flock2.has_lock)
        self.assertRaises(FileLockError,
                          flock2.acquire)
        self.assertFalse(flock2.has_lock)
        # Release first lock
        flock.release()
        self.assertFalse(flock.has_lock)
        # Check second lock can now be
        # acquired
        flock2.acquire()
        self.assertTrue(flock2.has_lock)
        flock2.release()
        self.assertFalse(flock2.has_lock)

    @unittest.skip("Not implemented")
    def test_file_lock_timeout(self):
        """
        FileLock: use timeout when acquiring lock
        """
        pass

    def test_file_lock_using_with(self):
        """
        FileLock: using as part of 'with' statement
        """
        # Make a file to lock
        f = os.path.join(self.wd,"test.txt")
        with open(f,'wt') as fp:
            fp.write("Text\n")
        # 'Probe' file lock instance
        flock = FileLock(f)
        self.assertFalse(flock.has_lock)
        # Acquire lock using with statement
        with FileLock(f) as fl:
            self.assertTrue(fl.has_lock)
            self.assertRaises(FileLockError,
                              flock.acquire)
        # Exiting context should release the lock
        # and enable another file lock to be acquired
        flock.acquire()
        self.assertTrue(flock.has_lock)
        flock.release()

class TestLocation(unittest.TestCase):
    """
    Tests for the Location class
    """
    def test_location_user_server_path(self):
        """
        Location: handle user@host.name:/path/to/somewhere
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
        """
        Location: handle host.name:/path/to/somewhere
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
        """
        Location: handle /path/to/somewhere
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
        """
        Location: handle http://example.com/path/to/url
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
        """
        Location: handle https://example.com/path/to/url
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
        """
        Location: handle file:///path/to/file
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

class TestBasesMaskIsPairedEnd(unittest.TestCase):
    """Tests for the bases_mask_is_paired_end function

    """

    def test_single_end_bases_mask(self):
        """Check examples of single-ended bases masks
        """
        self.assertFalse(bases_mask_is_paired_end('y36'))
        self.assertFalse(bases_mask_is_paired_end('y50,I6'))
        self.assertFalse(bases_mask_is_paired_end('y50,I6n'))

    def test_paired_end_bases_mask(self):
        """Check examples of single-ended bases masks
        """
        self.assertTrue(bases_mask_is_paired_end('y101,I6n,y101'))
        self.assertTrue(bases_mask_is_paired_end('y101,I6,I6,y101'))

class TestGetOrganismList(unittest.TestCase):
    """Tests for the get_organism_list function
    """

    def test_get_organism_list(self):
        """get_organism_list: handle standard names
        """
        self.assertEqual(get_organism_list("Human"),["human"])
        self.assertEqual(get_organism_list("mouse"),["mouse"])
        self.assertEqual(get_organism_list("Human,Mouse"),
                         ["human","mouse"])

    def test_get_organism_list_handle_names_with_spaces(self):
        """get_organism_list: handle names with spaces
        """
        self.assertEqual(get_organism_list("Xenopus tropicalis"),
                         ["xenopus_tropicalis"])
        self.assertEqual(get_organism_list("Homo sapiens,Mus musculus"),
                         ["homo_sapiens","mus_musculus"])

    def test_get_organism_list_handle_extra_whitespace(self):
        """get_organism_list: handle names with extra whitespace
        """
        self.assertEqual(get_organism_list("Homo sapiens, Mus  musculus "),
                         ["homo_sapiens","mus_musculus"])

class TestNormaliseOrganismName(unittest.TestCase):
    """Tests for the normalise_organism_name function
    """

    def test_normalise_organism_name(self):
        """normalise_organism_name: check normalisation
        """
        self.assertEqual(normalise_organism_name("Mouse"),"mouse")
        self.assertEqual(normalise_organism_name("mouse"),"mouse")
        self.assertEqual(normalise_organism_name("Mus musculus"),
                         "mus_musculus")
        self.assertEqual(normalise_organism_name(" Mus  musculus "),
                         "mus_musculus")

class TestSplitUserHostDir(unittest.TestCase):
    """Tests for the split_user_host_dir function

    """

    def test_user_host_dir(self):
        """Check we can split user@host.name:/path/to/somewhere
        """
        user,host,dirn = split_user_host_dir('user@host.name:/path/to/somewhere')
        self.assertEqual(user,'user')
        self.assertEqual(host,'host.name')
        self.assertEqual(dirn,'/path/to/somewhere')

    def test_host_dir(self):
        """Check we can split host.name:/path/to/somewhere
        """
        user,host,dirn = split_user_host_dir('host.name:/path/to/somewhere')
        self.assertEqual(user,None)
        self.assertEqual(host,'host.name')
        self.assertEqual(dirn,'/path/to/somewhere')

    def test_dir(self):
        """Check we can 'split' /path/to/somewhere
        """
        user,host,dirn = split_user_host_dir('/path/to/somewhere')
        self.assertEqual(user,None)
        self.assertEqual(host,None)
        self.assertEqual(dirn,'/path/to/somewhere')

    def test_extra_whitespace(self):
        """Check we can handle additional leading/trailing whitespace
        """
        user,host,dirn = split_user_host_dir('\tuser@host.name:/path/to/somewhere  \n')
        self.assertEqual(user,'user')
        self.assertEqual(host,'host.name')
        self.assertEqual(dirn,'/path/to/somewhere')

    def test_bad_input(self):
        """Check that empty string or None return sensible values
        """
        user,host,dirn = split_user_host_dir('')
        self.assertEqual(user,None)
        self.assertEqual(host,None)
        self.assertEqual(dirn,None)
        user,host,dirn = split_user_host_dir(None)
        self.assertEqual(user,None)
        self.assertEqual(host,None)
        self.assertEqual(dirn,None)

class TestGetNumberedSubdir(unittest.TestCase):
    """
    Tests for the 'get_numbered_subdir' function
    """
    def setUp(self):
        # Make temporary working dir
        self.wd = tempfile.mkdtemp(suffix="TestGetNumberedSubdir")
        # Store cwd
        self.original_dir = os.getcwd()
    def tearDown(self):
        # Restore origin cwd
        os.chdir(self.original_dir)
        # Remove working dir
        shutil.rmtree(self.wd)
    def test_get_numbered_subdir(self):
        """
        get_numbered_subdir: first subdir in parent
        """
        self.assertEqual(
            get_numbered_subdir("test",parent_dir=self.wd),
            "001_test")
    def test_get_numbered_subdir_next_in_sequence(self):
        """
        get_numbered_subdir: next subdir in sequence
        """
        for d in ("001_test","002_test"):
            os.mkdir(os.path.join(self.wd,d))
        self.assertEqual(
            get_numbered_subdir("test",parent_dir=self.wd),
            "003_test")
    def test_get_numbered_subdir_ignore_gaps(self):
        """
        get_numbered_subdir: ignore gaps in sequence
        """
        for d in ("001_test","003_test"):
            os.mkdir(os.path.join(self.wd,d))
        self.assertEqual(
            get_numbered_subdir("test",parent_dir=self.wd),
            "004_test")
    def test_get_numbered_subdir_default_to_cwd(self):
        """
        get_numbered_subdir: parent defaults to CWD
        """
        test_dir1 = os.path.join(self.wd,"test1")
        os.mkdir(test_dir1)
        for d in ("001_test","002_test"):
            os.mkdir(os.path.join(test_dir1,d))
        os.chdir(test_dir1)
        self.assertEqual(
            get_numbered_subdir("test"),"003_test")
        test_dir2 = os.path.join(self.wd,"test2")
        os.mkdir(test_dir2)
        os.chdir(test_dir2)
        self.assertEqual(
            get_numbered_subdir("test"),"001_test")
    def test_get_numbered_subdir_full_path(self):
        """
        get_numbered_subdir: return full path
        """
        self.assertEqual(
            get_numbered_subdir("test",
                                parent_dir=self.wd,
                                full_path=True),
            os.path.join(self.wd,"001_test"))

class TestFindExecutables(unittest.TestCase):
    """
    Tests for the find_executables function

    """
    def setUp(self):
        # Make temporary working dir
        self.wd = tempfile.mkdtemp(suffix="TestCellrangerInfo")
        # Store the original state of PATH env var
        self.original_path = os.environ['PATH']
        # Info function
        def info_func(p):
            name = os.path.basename(p)
            exe = find_program(p)
            version = ''
            output = Command(exe).subprocess_check_output()[1]
            for line in output.split('\n'):
                if line.startswith(name):
                    try:
                        version = line.split()[1]
                    except IndexError:
                        pass
                    break
            return (exe,name.upper(),version)
        self.info_func = info_func

    def tearDown(self):
        # Restore the PATH env var
        os.environ['PATH'] = self.original_path
        # Remove temp dir
        shutil.rmtree(self.wd)

    def _make_executable(self,name,dirn=None,version=None):
        # Create an executable for testing
        # The executable will be created in the working dir,
        # optionally under the specified subdir
        if dirn is not None:
            path = self.wd
            for d in dirn.split(os.sep):
                path = os.path.join(path,d)
                os.mkdir(path)
            dirn = path
        else:
            dirn = self.wd
        exe = os.path.join(dirn,name)
        with open(exe,'w') as fp:
            fp.write(
                "#!/bin/bash\necho %s %s" %
                (name,
                 (version if version is not None else "")))
        os.chmod(exe,0o775)
        return exe

    def _info_func(self):
        # Create an info_func for testing
        def info_func(p):
            name = os.path.basename(p)
            exe = find_program(p)
            version = ''
            output = Command(exe).subprocess_check_output()[1]
            for line in output.split('\n'):
                if line.startswith(name):
                    try:
                        version = line.split()[1]
                    except IndexError:
                        pass
                    break
            return (exe,name.upper(),version)
        return info_func

    def _prepend_path(self,p):
        # Add to start of PATH
        os.environ['PATH'] = "%s:%s" % (p,os.environ['PATH'])

    def test_no_versions(self):
        """
        find_executables: no matching software present
        """
        os.environ['PATH'] = ''
        self.assertEqual(find_executables(("xyxyz",),self.info_func),[])

    def test_single_version(self):
        """
        find_executable: only single version available
        """
        xyxyz = self._make_executable("xyxyz",version="1.0.0")
        self._prepend_path(os.path.dirname(xyxyz))
        self.assertEqual(find_executables(("xyxyz",),self._info_func()),[xyxyz,])
    def test_multiple_versions(self):
        """
        find_executable: multiple versions available
        """
        exes = []
        for v in ('2.0.0','1.2.1','1.0.0'):
            exe = self._make_executable("xyxyz",
                                        dirn="%s/bin" % v,
                                        version=v)
            exes.insert(0,exe)
            self._prepend_path(os.path.dirname(exe))
        self.assertEqual(find_executables(("xyxyz",),self._info_func()),exes)

    def test_require_specific_version(self):
        """
        find_executable: require specific version
        """
        exes = []
        for v in ('2.0.0','1.2.1','1.0.0'):
            exe = self._make_executable("xyxyz",
                                        dirn="%s/bin" % v,
                                        version=v)
            if v == '1.2.1':
                exes.append(exe)
            self._prepend_path(os.path.dirname(exe))
        self.assertEqual(find_executables(("xyxyz",),
                                          self._info_func(),
                                          reqs='1.2.1'),exes)

    def test_required_version_not_found(self):
        """
        find_executable: required version is not found
        """
        for v in ('2.0.0','1.2.1','1.0.0'):
            exe = self._make_executable("xyxyz",
                                        dirn="%s/bin" % v,
                                        version=v)
            self._prepend_path(os.path.dirname(exe))
        self.assertEqual(find_executables(("xyxyz",),
                                          self._info_func(),
                                          reqs='2.0.1'),[])

    def test_multiple_exes(self):
        """
        find_executable: look for multiple executables
        """
        exes = []
        names = ('xyxyz','zyzyx','pqrst')
        for n in names:
            exe = self._make_executable(n)
            if n != 'pqrst':
                exes.insert(0,exe)
        self._prepend_path(self.wd)
        self.assertEqual(find_executables(("zyzyx","xyxyz",),
                                          self._info_func()),exes)

    def test_specify_path(self):
        """
        find_executable: look for multiple executables
        """
        exes = []
        paths = []
        names = ('xyxyz','zyzyx','pqrst')
        for n in names:
            exe = self._make_executable(n,dirn=n)
            if n == 'zyzyx':
                exes.append(exe)
                paths.append(os.path.dirname(exe))
            self._prepend_path(os.path.dirname(exe))
        print(os.environ['PATH'])
        self.assertEqual(find_executables(("zyzyx","xyxyz",),
                                          self.info_func,
                                          paths=paths),exes)

class TestParseVersion(unittest.TestCase):
    """Tests for the parse_version function
    """
    def test_parse_version(self):
        """parse_version splits version numbers
        """
        self.assertEqual(parse_version("2.17.1.4"),
                         (2,17,1,4))
        self.assertEqual(parse_version("1.8.4"),
                         (1,8,4))
        self.assertEqual(parse_version("1.9.rc1"),
                         (1,9,"rc1"))

    def test_compare_versions(self):
        """parse_version compares versions correctly
        """
        self.assertTrue(
            parse_version("2.17.1.4") > parse_version("1.8.4"))
        self.assertTrue(
            parse_version("2.17.1.4") == parse_version("2.17.1.4"))
        self.assertTrue(
            parse_version("1.8.4") < parse_version("2.17.1.4"))
        self.assertTrue(
            parse_version("10.0") > parse_version("9.1"))
        self.assertTrue(
            parse_version("1.10") > parse_version("1.9.rc1"))

    def test_handle_empty_version(self):
        """parse_version handles empty version
        """
        self.assertEqual(parse_version(""),(-99999,))
        self.assertTrue(
            parse_version("") < parse_version("1.8.4"))

class TestPrettyPrintRows(unittest.TestCase):
    """Tests for the pretty_print_rows function

    """
    def test_pretty_print_rows_empty_inputs(self):
        """pretty_print_rows can handle empty list
        """
        self.assertEqual(pretty_print_rows([]),"")

    def test_pretty_print_rows_single_row(self):
        """pretty_print_rows prints a single row
        """
        rows = [['hello:','A salutation']]
        self.assertEqual(pretty_print_rows(rows),
                         "hello: A salutation")

    def test_pretty_print_rows_multiple_rows(self):
        """pretty_print_rows prints multiple rows
        """
        rows = [['-','hello:','A salutation'],
                ['-','goodbye:','The End']]
        self.assertEqual(pretty_print_rows(rows),
                         "- hello:   A salutation\n- goodbye: The End     ")

    def test_pretty_print_rows_prepend(self):
        """pretty_print_rows prepend right-justifies values
        """
        rows = [['-','hello:','A salutation'],
                ['-','goodbye:','The End']]
        self.assertEqual(pretty_print_rows(rows,prepend=True),
                         "-   hello: A salutation\n- goodbye:      The End")

class TestWriteScriptFile(unittest.TestCase):
    """Tests for the write_script_file function

    """
    def setUp(self):
        self.dirn = tempfile.mkdtemp(suffix='TestWriteScripFile')

    def tearDown(self):
        # Remove the temporary test directory
        shutil.rmtree(self.dirn)

    def test_write_script_file(self):
        """write_script_file creates new script file
        """
        script_file = os.path.join(self.dirn,'test.sh')
        self.assertFalse(os.path.exists(script_file))
        write_script_file(script_file,"sleep 50\n#")
        self.assertTrue(os.path.exists(script_file))
        with open(script_file,'rt') as fp:
            self.assertEqual(fp.read(),"sleep 50\n#\n")

    def test_write_script_file_with_shell(self):
        """write_script_file creates new script file with shell
        """
        script_file = os.path.join(self.dirn,'test.sh')
        self.assertFalse(os.path.exists(script_file))
        write_script_file(script_file,"sleep 50\n#",shell='/bin/sh')
        self.assertTrue(os.path.exists(script_file))
        with open(script_file,'rt') as fp:
            self.assertEqual(fp.read(),"#!/bin/sh\nsleep 50\n#\n")

    def test_write_script_file_append(self):
        """write_script_file appends to an existing file
        """
        script_file = os.path.join(self.dirn,'test.sh')
        with open(script_file,'w') as fp:
            fp.write("#\necho Going to sleep\n")
        write_script_file(script_file,"sleep 50\n#",append=True)
        self.assertTrue(os.path.exists(script_file))
        with open(script_file,'rt') as fp:
            self.assertEqual(fp.read(),
                             "#\necho Going to sleep\nsleep 50\n#\n")
