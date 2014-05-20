#!/bin/env python
#
#     data_manager.py: sequencing data curation utility
#     Copyright (C) University of Manchester 2014 Peter Briggs
#
########################################################################
#
# data_manager.py
#
#########################################################################


"""Sequence data and analysis curation utility

"""

#######################################################################
# Module metadata
#######################################################################

__version__ = "0.0.27"

#######################################################################
# Import modules that this module depends on
#######################################################################

import sys
import os
import re
import optparse
import logging
import subprocess
import bcf_utils
import Md5sum

#######################################################################
# Classes
#######################################################################

class DataDir:
    """
    """
    def __init__(self,data_dir):
        """
        """
        data_dir = os.path.abspath(data_dir)
        if not os.path.isdir(data_dir):
            raise OSError, "'%s': not a directory" % data_dir
        self.__data_dir = data_dir
    @property
    def dir(self):
        """
        """
        return self.__data_dir
    def get_size(self,pattern=None,block_size=1,human_readable=False):
        """Return the total size of some or all files

        By default this returns the total size of the data dir
        in bytes.

        Specifying a regex pattern restricts the total to be
        only those files and directories which match the pattern.

        The format of the output value can be modified by
        adjusting the block_size (e.g. 1024 == Kb, 1024*1024 == Mb
        etc) or by setting the human_readable flag to True, in
        which case the value will be expressed in the most
        appropriate units (and with the units appended).

        Note: values returned by this method seem most closely to
        match the output of the 'du --apparent-size' command,
        although they don't appear to match exactly.

        Arguments:
          pattern: if not None, restrict examined files and
            directories to those matching the supplied regex
            pattern (default is to include all files and dirs)
          block_size: report size as number of blocks of
            'block_size' bytes (default is 1)
          human_readable: report size in appropriate units,
            with unit identifier appended e.g. 1G

        Returns
          Total size for specified files/dirs in request units.

        """
        size = 0
        for f in self.walk(pattern=pattern):
            size += os.lstat(f).st_size
        #return int(float(size)/1024)
        if human_readable:
            # Format into human-readable format
            return bcf_utils.format_file_size(size)
        else:
            # Return number of blocks, rounded
            # up or down to nearest integer
            return int(round(float(size)/float(block_size)))
    def walk(self,include_dirs=True,pattern=None):
        """Traverse the directory, subdirectories and files
        
        Arguments:
          include_dirs: if True then yield directories as well
            as files (default)
          pattern: if not None then specifies a regular expression
            pattern which restricts the set of yielded files and
            directories to a subset of those which match the
            pattern
        
        """
        if pattern is not None:
            matcher = re.compile(pattern)
        if include_dirs:
            if pattern is None or matcher.match(self.dir):
                yield self.dir
        for dirpath,dirnames,filenames in os.walk(self.dir):
            if include_dirs:
                for d in dirnames:
                    d1 = os.path.join(dirpath,d)
                    if pattern is None or matcher.match(d1):
                        yield d1
            for f in filenames:
                f1 = os.path.join(dirpath,f)
                if pattern is None or matcher.match(f1):
                    yield f1
    def copy(self,target,dry_run=False,verbose=False):
        """Create a copy of dir using rsync

        Creates a copy of the data directory under directory 'target',
        using the external 'rsync' program.

        The copy will be owned by the user running the operation.
        Other rsync options are:

        -r: recursive copy
        -l: copy symbolic links as links
        -t: preserve timestamps
        -E: preserve executibility
        -D: preserve devices files and special files
        --chmod=u+rwX,g+rwX,o-w: give user and group read and write access
        -o: preserve owner (nb superuser only)

        -p (preserve permissions) is also needed to make the --chmod
        options stick (see rsync manpage).

        Arguments:
          target : target directory to make the copy under
          dry_run: optional, if True then run rsync with --dry-run option
            so no files are copied.
          verbose: optional, if True then run rsync with -v option

        """
        rsync_options = '-rplotDE'
        if verbose:
            rsync_options += 'v'
            stdout = None
        else:
            # Assign stdout to /dev/null so that even error messages
            # will be suppressed
            # This is useful for nosetests, which uses cIOString to
            # capture stdout and so can break when combined with output
            # stream reassignment in subprocess
            stdout = open('/dev/null','w')
        rsync = ['rsync',rsync_options,'--chmod=u+rwX,g+rwX,o-w',self.dir,target]
        if dry_run:
            rsync.insert(1,'--dry-run')
        logging.debug("Rsync command: %s" % rsync)
        try:
            # Use 'check_call' as we want to see output from rsync
            status = subprocess.check_call(rsync,stdout=stdout,stderr=stdout)
        except subprocess.CalledProcessError,ex:
            # Non-zero exit status
            print "%s" % ex
            status = ex.returncode
        logging.debug("Rsync exit status: %s" % status)
        return status
    def find_files(self,regex):
        """
        """
        matcher = re.compile(regex)
        for filen in self.files:
            if matcher.match(filen):
                yield filen
    @property
    def symlinks(self):
        """
        """
        for filen in self.files:
            if os.path.islink(filen):
                yield filen
    @property
    def files(self):
        """
        """
        for dirpath,dirnames,filenames in os.walk(self.dir):
            for f in filenames:
                yield os.path.join(dirpath,f)
    @property
    def dirs(self):
        """
        """
        for dirpath,dirnames,filenames in os.walk(self.dir):
            for d in dirnames:
                yield os.path.join(dirpath,d)
    @property
    def users(self):
        """Return list of user names for this directory

        Walks the directory structure and returns a list of all
        user names which own files and/or directories within it.

        """
        users = []
        for f in self.walk():
            user = bcf_utils.PathInfo(f).user
            if user not in users:
                users.append(user)
        users.sort()
        return users
    def verify(self,dirn):
        """Verify another data directory using this one as a reference

        Returns
          Counter object.

        """
        # Set up DataDir instance for target directory
        data_dir = DataDir(dirn)
        # Set up result counter
        result = Counter()
        result.add_quantity("ok","OK")
        result.add_quantity("md5_failed","MD5 sum checks failed")
        result.add_quantity("links_differ","Symlinks differ")
        result.add_quantity("types_differ","Copy types differ from reference")
        result.add_quantity("missing","Missing from copy")
        result.add_quantity("unreadable_ref","Unreadable reference")
        # Walk the reference directory and check against copy
        for f in self.walk():
            f = bcf_utils.PathInfo(f).relpath(self.dir)
            ref = bcf_utils.PathInfo(f,basedir=self.dir)
            cpy = bcf_utils.PathInfo(f,basedir=data_dir.dir)
            # Check copy exists
            if not ref.is_readable:
                print "UNREADABLE REFERENCE: %s" % f
                result.incr('unreadable_ref')
            elif not cpy.exists:
                print "MISSING: %s" % f
                result.incr('missing')
            elif ref.is_file:
                # Check copy is the same type
                if not cpy.is_file:
                    print "TYPE DIFFERS (NOT FILE): %s" % ref.relpath(self.dir)
                    result.incr('types_differ')
                else:
                    # Compare MD5sums
                    if Md5sum.Md5Checker.md5cmp_files(ref.path,cpy.path) == \
                       Md5sum.Md5Checker.MD5_OK:
                        print "MD5 OK: %s" % ref.relpath(self.dir)
                        result.incr('ok')
                    else:
                        print "MD5 FAILED: %s" % ref.relpath(self.dir)
                        result.incr('md5_failed')
            elif ref.is_link:
                # Check copy is the same type
                if not cpy.is_link:
                    print "TYPE DIFFERS (NOT SYMLINK): %s" % ref.relpath(self.dir)
                    result.incr('types_differ')
                else:
                    # Compare targets
                    if bcf_utils.Symlink(cpy.path).target == \
                       bcf_utils.Symlink(ref.path).target:
                        print "SYMLINK OK: %s" % ref.relpath(self.dir)
                        result.incr('ok')
                    else:
                        print "SYMLINK DIFFERS: %s" % ref.relpath(self.dir)
                        result.incr('links_differ')
            elif ref.is_dir:
                # Check copy is the same type
                if not cpy.is_dir:
                    print "TYPE DIFFERS (NOT DIR): %s" % ref.relpath(self.dir)
                    result.incr('types_differ')
                else:
                    print "DIR OK: %s" % ref.relpath(self.dir)
                    result.incr('ok')
        # Finished
        return result
    def diff(self,dirn):
        """Check whether files & directories are present in another directory
 
        Uses this directory as a reference and checks whether files,
        directories and links are also present in the copy.

        Use 'verify' method for a more rigorous verification which checks MD5
        sums etc.

        Returns
          Counter object.

        """
        # Set up DataDir instance for target directory
        data_dir = DataDir(dirn)
        # Set up result counter
        result = Counter()
        result.add_quantity("ok","OK")
        result.add_quantity("missing","Missing from copy")
        result.add_quantity("types_differ","Copy types differ from reference")
        result.add_quantity("unreadable_ref","Unreadable reference")
        # Walk the reference directory and check against copy
        for f in self.walk():
            f = bcf_utils.PathInfo(f).relpath(self.dir)
            ref = bcf_utils.PathInfo(f,basedir=self.dir)
            cpy = bcf_utils.PathInfo(f,basedir=data_dir.dir)
            # Check copy exists
            if not ref.is_readable:
                print "UNREADABLE REFERENCE: %s" % f
                result.incr('unreadable_ref')
            elif not cpy.exists:
                print "MISSING: %s" % f
                result.incr('missing')
            elif ref.is_file:
                # Check copy is the same type
                if not cpy.is_file:
                    print "TYPE DIFFERS (NOT FILE): %s" % ref.relpath(self.dir)
                    result.incr('types_differ')
                else:
                    print "OK: %s" % ref.relpath(self.dir)
                    result.incr('ok')
            elif ref.is_link:
                # Check copy is the same type
                if not cpy.is_link:
                    print "TYPE DIFFERS (NOT SYMLINK): %s" % ref.relpath(self.dir)
                    result.incr('types_differ')
                else:
                    print "OK: %s" % ref.relpath(self.dir)
                    result.incr('ok')
            elif ref.is_dir:
                # Check copy is the same type
                if not cpy.is_dir:
                    print "TYPE DIFFERS (NOT DIR): %s" % ref.relpath(self.dir)
                    result.incr('types_differ')
                else:
                    print "OK: %s" % ref.relpath(self.dir)
                    result.incr('ok')
        # Finished
        return result

class Counter(bcf_utils.AttributeDictionary):
    """Generic class for counting things

    The Counter class provides a generic object that can be used
    to count discrete arbitrary quantities.

    Basic usage:

    >>> c = Counter() # Create a new Counter object
    >>> c.add_quantity('cats') # Initialise counting of 'cats'
    >>> c.add_quantity('dogs') # Initialise counting of 'dogs'
    
    Then to record something:

    >>> c.incr('cats') # Increments the count of 'cats' by 1
    >>> c.incr('dogs',2) # Increments the count of 'dogs' by 2

    To get the count for something, do:

    >>> c.cats
    
    or
    
    >>> c['dogs']

    Use the 'report' method to get a summary of everything
    counted so far.

    """
    def __init__(self):
        """Create a new Counter instance

        """
        bcf_utils.AttributeDictionary.__init__(self)
        self.__quantities = bcf_utils.OrderedDictionary()
        self.__total = 0

    def add_quantity(self,name,description=None):
        """Register a quantity to be counted

        """
        if name in self:
            raise KeyError,"'%s' already exists" % name
        self[name] = 0
        if description is None:
            description = name
        self.__quantities[name] = description

    def incr(self,name,count=1):
        """Increment a counter

        """
        self[name] += count
        self.__total += count

    def total(self):
        """Return total counts

        """
        return self.__total

    def report(self):
        """Print a report

        """
        print "Summary:"
        if len(self) > 0:
            for name in self.__quantities:
                print "\t%d\t%s" % (self[name], self.__quantities[name])
        print "Total %d" % self.total()

#######################################################################
# Functions
#######################################################################

def get_size(f):
    if os.path.isfile(f):
        return float(os.lstat(f).st_size)/1024
    else:
        size = 0
        for dirpath,dirnames,filenames in os.walk(f):
            for d in dirnames:
                size += get_size(os.path.join(dirpath,d))
            for f in filenames:
                size += get_size(os.path.join(dirpath,f))
        return size

# http://stackoverflow.com/questions/6798097/find-regex-in-python-or-how-to-find-files-whose-whole-name-path-name
def find(dirn,regex):
    matcher = re.compile(regex)
    for dirpath, dirnames, filenames in os.walk(dirn):
        for f in filenames:
            f = os.path.relpath(os.path.join(dirpath, f), dirn)
            if matcher.match(f):
                yield f

def finger_user(name):
    # Run 'finger' command and try to fetch user's real name
    finger = bcf_utils.find_program('finger')
    logging.debug("'finger' executable: %s" % finger)
    if finger is None:
        logging.debug("'finger' not located on this system")
        return None
    p = subprocess.Popen([finger,name], stdout=subprocess.PIPE)
    finger, err = p.communicate()
    real_name = finger.split('\n')[0].split(':')[-1].strip()
    return real_name

#######################################################################
# Tests
#######################################################################
import unittest
from mock_data import TestUtils,ExampleDirLanguages

class TestDataDirCopy(unittest.TestCase):
    """Tests for DataDir class 'copy' functionality
    """
    def setUp(self):
        # Make a test data directory structure
        self.example_dir = ExampleDirLanguages()
        self.wd = self.example_dir.create_directory()
        self.dest = TestUtils.make_dir()
    def tearDown(self):
        # Make sure file is readable
        os.chmod(self.example_dir.path("hello"),0644)
        # Remove the test data directory (and copy)
        self.example_dir.delete_directory()
        TestUtils.remove_dir(self.dest)
    def test_copy(self):
        """DataDir.copy() correctly copies directory

        """
        # Make a copy
        data_dir = DataDir(self.wd)
        status = data_dir.copy(self.dest)
        # Check return status
        self.assertEqual(status,0)
        # Compare contents one way
        target = os.path.join(self.dest,os.path.basename(self.wd))
        for dirpath,dirnames,filenames in os.walk(self.wd):
            for d in dirnames:
                d2 = os.path.join(target,
                                  os.path.relpath(os.path.join(dirpath,d),self.wd))
                self.assertTrue(os.path.isdir(d2))
            for f in filenames:
                f2 = os.path.join(target,
                                  os.path.relpath(os.path.join(dirpath,f),self.wd))
                self.assertTrue(os.path.isfile(f2))
        # Compare contents the other way
        for dirpath,dirnames,filenames in os.walk(target):
            for d in dirnames:
                d2 = os.path.join(self.wd,
                                  os.path.relpath(os.path.join(dirpath,d),target))
                self.assertTrue(os.path.isdir(d2))
            for f in filenames:
                f2 = os.path.join(self.wd,
                                  os.path.relpath(os.path.join(dirpath,f),target))
                self.assertTrue(os.path.isfile(f2))
    def test_copy_unreadable_file(self):
        """DataDir.copy() reports failure for unreadable source file

        """
        # Make an unreadable file in the reference
        os.chmod(self.example_dir.path("hello"),0000)
        # Make a copy
        data_dir = DataDir(self.wd)
        status = data_dir.copy(self.dest)
        # Check return status
        self.assertEqual(status,23)

class TestDataDirWalk(unittest.TestCase):
    """Tests for DataDir class 'walk' functionality
    """
    def setUp(self):
        # Make a test data directory structure
        self.example_dir = ExampleDirLanguages()
        self.wd = self.example_dir.create_directory()
    def tearDown(self):
        # Remove the test data directory
        self.example_dir.delete_directory()
    def test_walk(self):
        """DataDir.walk() traverses all files and directories

        """
        filelist = self.example_dir.filelist(include_dirs=True)
        filelist.append(self.wd)
        data_dir = DataDir(self.wd)
        for f in data_dir.walk():
            self.assertTrue(f in filelist,"%s not expected" % f)
            filelist.remove(f)
        self.assertEqual(len(filelist),0,"Items not returned: %s" %
                         ','.join(filelist))
    def test_walk_no_directories(self):
        """DataDir.walk() traverses all files and ignores directories

        """
        filelist = self.example_dir.filelist(include_dirs=False)
        data_dir = DataDir(self.wd)
        for f in data_dir.walk(include_dirs=False):
            self.assertTrue(f in filelist,"%s not expected" % f)
            filelist.remove(f)
        self.assertEqual(len(filelist),0,"Items not returned: %s" %
                         ','.join(filelist))
    def test_walk_includes_hidden_files_and_directories(self):
        """DataDir.walk() finds 'hidden' files and directories

        """
        self.example_dir.add_file(".hidden_file")
        self.example_dir.add_dir(".hidden_dir")
        self.example_dir.add_file(".hidden_dir/test")
        filelist = self.example_dir.filelist(include_dirs=True)
        filelist.append(self.wd)
        data_dir = DataDir(self.wd)
        for f in data_dir.walk():
            self.assertTrue(f in filelist,"%s not expected" % f)
            filelist.remove(f)
        self.assertEqual(len(filelist),0,"Items not returned: %s" %
                         ','.join(filelist))

class TestDataDirVerify(unittest.TestCase):
    """Tests for the DataDir class 'verify' functionality
    """
    def setUp(self):
        # Make a test and reference data directory structures
        self.example_dir = ExampleDirLanguages()
        self.reference = ExampleDirLanguages()
        self.wd = self.example_dir.create_directory()
        self.ref = self.reference.create_directory()
    def tearDown(self):
        # Remove the test data directories
        self.example_dir.delete_directory()
        self.reference.delete_directory()
    def test_verify_exact_copy(self):
        """DataDir.verify() confirms exact match for exact copy

        """
        ref_dir = DataDir(self.ref)
        result = ref_dir.verify(self.wd)
        self.assertTrue(result.total() > 0)
        self.assertEqual(result.ok,result.total())
        self.assertEqual(result.md5_failed,0)
        self.assertEqual(result.links_differ,0)
        self.assertEqual(result.types_differ,0)
        self.assertEqual(result.missing,0)
        self.assertEqual(result.unreadable_ref,0)
    def test_verify_missing_file(self):
        """DataDir.verify() confirms missing file in copy

        """
        # Remove a file from the copy
        os.remove(self.example_dir.path("hello"))
        # Do the verification
        ref_dir = DataDir(self.ref)
        result = ref_dir.verify(self.wd)
        self.assertTrue(result.total() > 0)
        self.assertNotEqual(result.ok,result.total())
        self.assertEqual(result.md5_failed,0)
        self.assertEqual(result.links_differ,0)
        self.assertEqual(result.types_differ,0)
        self.assertEqual(result.missing,1)
        self.assertEqual(result.unreadable_ref,0)
        self.assertEqual(result.missing+result.ok,result.total())
    def test_verify_missing_empty_dir(self):
        """DataDir.verify() confirms missing empty directory in copy

        """
        # Add an empty directory to the reference
        self.reference.add_dir("empty")
        # Do the verification
        ref_dir = DataDir(self.ref)
        result = ref_dir.verify(self.wd)
        self.assertTrue(result.total() > 0)
        self.assertNotEqual(result.ok,result.total())
        self.assertEqual(result.md5_failed,0)
        self.assertEqual(result.links_differ,0)
        self.assertEqual(result.types_differ,0)
        self.assertEqual(result.missing,1)
        self.assertEqual(result.unreadable_ref,0)
        self.assertEqual(result.missing+result.ok,result.total())
    def test_verify_missing_dir(self):
        """DataDir.verify() confirms missing directory in copy

        """
        # Remove a directory (and its files) from the copy
        os.remove(self.example_dir.path("countries/spain"))
        os.remove(self.example_dir.path("countries/north_wales"))
        os.remove(self.example_dir.path("countries/south_wales"))
        os.remove(self.example_dir.path("countries/iceland"))
        os.rmdir(self.example_dir.path("countries"))
        # Do the verification
        ref_dir = DataDir(self.ref)
        result = ref_dir.verify(self.wd)
        self.assertTrue(result.total() > 0)
        self.assertNotEqual(result.ok,result.total())
        self.assertEqual(result.md5_failed,0)
        self.assertEqual(result.links_differ,0)
        self.assertEqual(result.types_differ,0)
        self.assertEqual(result.missing,5)
        self.assertEqual(result.unreadable_ref,0)
        self.assertEqual(result.missing+result.ok,result.total())
    def test_verify_missing_symlink(self):
        """DataDir.verify() confirms missing symlink in copy

        """
        # Remove a symlink from the copy
        os.remove(self.example_dir.path("hi"))
        # Do the verification
        ref_dir = DataDir(self.ref)
        result = ref_dir.verify(self.wd)
        self.assertTrue(result.total() > 0)
        self.assertNotEqual(result.ok,result.total())
        self.assertEqual(result.md5_failed,0)
        self.assertEqual(result.links_differ,0)
        self.assertEqual(result.types_differ,0)
        self.assertEqual(result.missing,1)
        self.assertEqual(result.unreadable_ref,0)
        self.assertEqual(result.missing+result.ok,result.total())
    def test_verify_replace_file_with_directory(self):
        """DataDir.verify() confirms file replaced with a directory

        """
        # Make a file in the reference and a directory in the copy
        # with the same name
        self.reference.add_file("bonjour","Hello!")
        self.example_dir.add_dir("bonjour")
        # Do the verification
        ref_dir = DataDir(self.ref)
        result = ref_dir.verify(self.wd)
        self.assertTrue(result.total() > 0)
        self.assertNotEqual(result.ok,result.total())
        self.assertEqual(result.md5_failed,0)
        self.assertEqual(result.links_differ,0)
        self.assertEqual(result.types_differ,1)
        self.assertEqual(result.missing,0)
        self.assertEqual(result.unreadable_ref,0)
        self.assertEqual(result.types_differ+result.ok,result.total())
    def test_verify_replace_file_with_symlink(self):
        """DataDir.verify() confirms file replaced with a symlink

        """
        # Make a file in the reference and a symlink in the copy
        # with the same name
        self.reference.add_file("bonjour","Hello!")
        self.example_dir.add_link("bonjour","hello")
        # Do the verification
        ref_dir = DataDir(self.ref)
        result = ref_dir.verify(self.wd)
        self.assertTrue(result.total() > 0)
        self.assertNotEqual(result.ok,result.total())
        self.assertEqual(result.md5_failed,0)
        self.assertEqual(result.links_differ,0)
        self.assertEqual(result.types_differ,1)
        self.assertEqual(result.missing,0)
        self.assertEqual(result.unreadable_ref,0)
        self.assertEqual(result.types_differ+result.ok,result.total())
    def test_verify_replace_symlink_with_file(self):
        """DataDir.verify() confirms symlink replaced with a file

        """
        # Make a symlink in the reference and a file in the copy
        # with the same name
        self.reference.add_link("bonjour","hello")
        self.example_dir.add_file("bonjour","Hello!")
        # Do the verification
        ref_dir = DataDir(self.ref)
        result = ref_dir.verify(self.wd)
        self.assertTrue(result.total() > 0)
        self.assertNotEqual(result.ok,result.total())
        self.assertEqual(result.md5_failed,0)
        self.assertEqual(result.links_differ,0)
        self.assertEqual(result.types_differ,1)
        self.assertEqual(result.missing,0)
        self.assertEqual(result.unreadable_ref,0)
        self.assertEqual(result.types_differ+result.ok,result.total())
    def test_verify_replace_symlink_with_directory(self):
        """DataDir.verify() confirms symlink replaced with a directory

        """
        # Make a symlink in the reference and a directory in the copy
        # with the same name
        self.reference.add_link("bonjour","hello")
        self.example_dir.add_dir("bonjour")
        # Do the verification
        ref_dir = DataDir(self.ref)
        result = ref_dir.verify(self.wd)
        self.assertTrue(result.total() > 0)
        self.assertNotEqual(result.ok,result.total())
        self.assertEqual(result.md5_failed,0)
        self.assertEqual(result.links_differ,0)
        self.assertEqual(result.types_differ,1)
        self.assertEqual(result.missing,0)
        self.assertEqual(result.unreadable_ref,0)
        self.assertEqual(result.types_differ+result.ok,result.total())
    def test_verify_replace_directory_with_file(self):
        """DataDir.verify() confirms directory replaced with a file

        """
        # Make a directory in the reference and a file in the copy
        # with the same name
        self.reference.add_dir("bonjour")
        self.example_dir.add_file("bonjour","Hello!")
        # Do the verification
        ref_dir = DataDir(self.ref)
        result = ref_dir.verify(self.wd)
        self.assertTrue(result.total() > 0)
        self.assertNotEqual(result.ok,result.total())
        self.assertEqual(result.md5_failed,0)
        self.assertEqual(result.links_differ,0)
        self.assertEqual(result.types_differ,1)
        self.assertEqual(result.missing,0)
        self.assertEqual(result.unreadable_ref,0)
        self.assertEqual(result.types_differ+result.ok,result.total())
    def test_verify_replace_directory_with_symlink(self):
        """DataDir.verify() confirms directory replaced with a symlink

        """
        # Make a directory in the reference and a symlink in the copy
        # with the same name
        self.reference.add_dir("bonjour")
        self.example_dir.add_link("bonjour","hello")
        # Do the verification
        ref_dir = DataDir(self.ref)
        result = ref_dir.verify(self.wd)
        self.assertTrue(result.total() > 0)
        self.assertNotEqual(result.ok,result.total())
        self.assertEqual(result.md5_failed,0)
        self.assertEqual(result.links_differ,0)
        self.assertEqual(result.types_differ,1)
        self.assertEqual(result.missing,0)
        self.assertEqual(result.unreadable_ref,0)
        self.assertEqual(result.types_differ+result.ok,result.total())
    def test_verify_md5sums_differ(self):
        """DataDir.verify() confirms files with different MD5 sums

        """
        # Change content of file in copy
        open(self.example_dir.path("icelandic/takk_fyrir"),'w').write("Thanking you")
        # Do the verification
        ref_dir = DataDir(self.ref)
        result = ref_dir.verify(self.wd)
        self.assertTrue(result.total() > 0)
        self.assertNotEqual(result.ok,result.total())
        self.assertEqual(result.md5_failed,1)
        self.assertEqual(result.links_differ,0)
        self.assertEqual(result.types_differ,0)
        self.assertEqual(result.missing,0)
        self.assertEqual(result.unreadable_ref,0)
        self.assertEqual(result.md5_failed+result.ok,result.total())
    def test_verify_symlink_targets_differ(self):
        """DataDir.verify() confirms symlinks with different targets

        """
        # Change link target in copy
        os.remove(self.example_dir.path("hi"))
        os.symlink(self.example_dir.path("goodbye"),
                   self.example_dir.path("hi"))
        # Do the verification
        ref_dir = DataDir(self.ref)
        result = ref_dir.verify(self.wd)
        self.assertTrue(result.total() > 0)
        self.assertNotEqual(result.ok,result.total())
        self.assertEqual(result.md5_failed,0)
        self.assertEqual(result.links_differ,1)
        self.assertEqual(result.types_differ,0)
        self.assertEqual(result.missing,0)
        self.assertEqual(result.unreadable_ref,0)
        self.assertEqual(result.links_differ+result.ok,result.total())
    def test_verify_broken_symlinks_ok(self):
        """DataDir.verify() confirms matching broken symlinks

        """
        # Make matching broken links in reference and copy
        self.reference.add_link("broken.txt","missing.txt")
        self.example_dir.add_link("broken.txt","missing.txt")
        # Do the verification
        ref_dir = DataDir(self.ref)
        result = ref_dir.verify(self.wd)
        self.assertTrue(result.total() > 0)
        self.assertEqual(result.ok,result.total())
        self.assertEqual(result.md5_failed,0)
        self.assertEqual(result.links_differ,0)
        self.assertEqual(result.types_differ,0)
        self.assertEqual(result.missing,0)
        self.assertEqual(result.unreadable_ref,0)
    def test_verify_unreadable_reference(self):
        """DataDir.verify() confirms unreadable reference file

        """
        # Make an unreadable file in the reference
        os.chmod(self.reference.path("hello"),0000)
        # Do the verification
        ref_dir = DataDir(self.ref)
        result = ref_dir.verify(self.wd)
        os.chmod(self.reference.path("hello"),0644)
        self.assertTrue(result.total() > 0)
        self.assertNotEqual(result.ok,result.total())
        self.assertEqual(result.md5_failed,0)
        self.assertEqual(result.links_differ,0)
        self.assertEqual(result.types_differ,0)
        self.assertEqual(result.missing,0)
        self.assertEqual(result.unreadable_ref,1)
        self.assertEqual(result.unreadable_ref+result.ok,result.total())
    def test_verify_unreadable_copy(self):
        """DataDir.verify() confirms unreadable copy

        """
        # Make an unreadable file in the reference
        os.chmod(self.example_dir.path("hello"),0000)
        # Do the verification
        ref_dir = DataDir(self.ref)
        result = ref_dir.verify(self.wd)
        os.chmod(self.example_dir.path("hello"),0644)
        self.assertTrue(result.total() > 0)
        self.assertNotEqual(result.ok,result.total())
        self.assertEqual(result.md5_failed,1)
        self.assertEqual(result.links_differ,0)
        self.assertEqual(result.types_differ,0)
        self.assertEqual(result.missing,0)
        self.assertEqual(result.unreadable_ref,0)
        self.assertEqual(result.md5_failed+result.ok,result.total())
    def test_verify_hidden_files_directories_and_symlinks(self):
        """DataDir.verify() confirms 'hidden' files, directories and links

        """
        # Make hidden file, directory and symlink
        self.reference.add_dir(".hidden")
        self.reference.add_file(".hidden/file")
        self.reference.add_file(".hiding")
        self.reference.add_link(".hiding2","hello")
        # Do the verification
        ref_dir = DataDir(self.ref)
        result = ref_dir.verify(self.wd)
        self.assertTrue(result.total() > 0)
        self.assertNotEqual(result.ok,result.total())
        self.assertEqual(result.md5_failed,0)
        self.assertEqual(result.links_differ,0)
        self.assertEqual(result.types_differ,0)
        self.assertEqual(result.missing,4)
        self.assertEqual(result.unreadable_ref,0)
        self.assertEqual(result.ok+result.missing,result.total())

class TestDataDirDiff(unittest.TestCase):
    """Tests for the DataDir class 'diff' functionality
    """
    def setUp(self):
        # Make a test and reference data directory structures
        self.example_dir = ExampleDirLanguages()
        self.reference = ExampleDirLanguages()
        self.wd = self.example_dir.create_directory()
        self.ref = self.reference.create_directory()
    def tearDown(self):
        # Remove the test data directories
        self.example_dir.delete_directory()
        self.reference.delete_directory()
    def test_diff_exact_copy(self):
        """DataDir.diff() confirms exact match for exact copy

        """
        ref_dir = DataDir(self.ref)
        result = ref_dir.diff(self.wd)
        self.assertTrue(result.total() > 0)
        self.assertEqual(result.ok,result.total())
        self.assertEqual(result.missing,0)
        self.assertEqual(result.types_differ,0)
        self.assertEqual(result.unreadable_ref,0)
    def test_diff_missing_file(self):
        """DataDir.diff() confirms missing file in copy

        """
        # Remove a file from the copy
        os.remove(self.example_dir.path("hello"))
        # Do the verification
        ref_dir = DataDir(self.ref)
        result = ref_dir.diff(self.wd)
        self.assertTrue(result.total() > 0)
        self.assertNotEqual(result.ok,result.total())
        self.assertEqual(result.missing,1)
        self.assertEqual(result.types_differ,0)
        self.assertEqual(result.unreadable_ref,0)
        self.assertEqual(result.missing+result.ok,result.total())
    def test_diff_missing_empty_dir(self):
        """DataDir.diff() confirms missing empty directory in copy

        """
        # Add an empty directory to the reference
        self.reference.add_dir("empty")
        # Do the verification
        ref_dir = DataDir(self.ref)
        result = ref_dir.diff(self.wd)
        self.assertTrue(result.total() > 0)
        self.assertNotEqual(result.ok,result.total())
        self.assertEqual(result.types_differ,0)
        self.assertEqual(result.missing,1)
        self.assertEqual(result.unreadable_ref,0)
        self.assertEqual(result.missing+result.ok,result.total())
    def test_diff_missing_dir(self):
        """DataDir.diff() confirms missing directory in copy

        """
        # Remove a directory (and its files) from the copy
        os.remove(self.example_dir.path("countries/spain"))
        os.remove(self.example_dir.path("countries/north_wales"))
        os.remove(self.example_dir.path("countries/south_wales"))
        os.remove(self.example_dir.path("countries/iceland"))
        os.rmdir(self.example_dir.path("countries"))
        # Do the verification
        ref_dir = DataDir(self.ref)
        result = ref_dir.diff(self.wd)
        self.assertTrue(result.total() > 0)
        self.assertNotEqual(result.ok,result.total())
        self.assertEqual(result.types_differ,0)
        self.assertEqual(result.missing,5)
        self.assertEqual(result.unreadable_ref,0)
        self.assertEqual(result.missing+result.ok,result.total())
    def test_diff_missing_symlink(self):
        """DataDir.diff() confirms missing symlink in copy

        """
        # Remove a symlink from the copy
        os.remove(self.example_dir.path("hi"))
        # Do the verification
        ref_dir = DataDir(self.ref)
        result = ref_dir.diff(self.wd)
        self.assertTrue(result.total() > 0)
        self.assertNotEqual(result.ok,result.total())
        self.assertEqual(result.types_differ,0)
        self.assertEqual(result.missing,1)
        self.assertEqual(result.unreadable_ref,0)
        self.assertEqual(result.missing+result.ok,result.total())
    def test_diff_replace_file_with_directory(self):
        """DataDir.diff() confirms file replaced with a directory

        """
        # Make a file in the reference and a directory in the copy
        # with the same name
        self.reference.add_file("bonjour","Hello!")
        self.example_dir.add_dir("bonjour")
        # Do the verification
        ref_dir = DataDir(self.ref)
        result = ref_dir.diff(self.wd)
        self.assertTrue(result.total() > 0)
        self.assertNotEqual(result.ok,result.total())
        self.assertEqual(result.types_differ,1)
        self.assertEqual(result.missing,0)
        self.assertEqual(result.unreadable_ref,0)
        self.assertEqual(result.types_differ+result.ok,result.total())
    def test_diff_replace_file_with_symlink(self):
        """DataDir.diff() confirms file replaced with a symlink

        """
        # Make a file in the reference and a symlink in the copy
        # with the same name
        self.reference.add_file("bonjour","Hello!")
        self.example_dir.add_link("bonjour","hello")
        # Do the verification
        ref_dir = DataDir(self.ref)
        result = ref_dir.diff(self.wd)
        self.assertTrue(result.total() > 0)
        self.assertNotEqual(result.ok,result.total())
        self.assertEqual(result.types_differ,1)
        self.assertEqual(result.missing,0)
        self.assertEqual(result.unreadable_ref,0)
        self.assertEqual(result.types_differ+result.ok,result.total())
    def test_diff_replace_symlink_with_file(self):
        """DataDir.diff() confirms symlink replaced with a file

        """
        # Make a symlink in the reference and a file in the copy
        # with the same name
        self.reference.add_link("bonjour","hello")
        self.example_dir.add_file("bonjour","Hello!")
        # Do the verification
        ref_dir = DataDir(self.ref)
        result = ref_dir.diff(self.wd)
        self.assertTrue(result.total() > 0)
        self.assertNotEqual(result.ok,result.total())
        self.assertEqual(result.types_differ,1)
        self.assertEqual(result.missing,0)
        self.assertEqual(result.unreadable_ref,0)
        self.assertEqual(result.types_differ+result.ok,result.total())
    def test_diff_replace_symlink_with_directory(self):
        """DataDir.diff() confirms symlink replaced with a directory

        """
        # Make a symlink in the reference and a directory in the copy
        # with the same name
        self.reference.add_link("bonjour","hello")
        self.example_dir.add_dir("bonjour")
        # Do the verification
        ref_dir = DataDir(self.ref)
        result = ref_dir.diff(self.wd)
        self.assertTrue(result.total() > 0)
        self.assertNotEqual(result.ok,result.total())
        self.assertEqual(result.types_differ,1)
        self.assertEqual(result.missing,0)
        self.assertEqual(result.unreadable_ref,0)
        self.assertEqual(result.types_differ+result.ok,result.total())
    def test_diff_replace_directory_with_file(self):
        """DataDir.diff() confirms directory replaced with a file

        """
        # Make a directory in the reference and a file in the copy
        # with the same name
        self.reference.add_dir("bonjour")
        self.example_dir.add_file("bonjour","Hello!")
        # Do the verification
        ref_dir = DataDir(self.ref)
        result = ref_dir.diff(self.wd)
        self.assertTrue(result.total() > 0)
        self.assertNotEqual(result.ok,result.total())
        self.assertEqual(result.types_differ,1)
        self.assertEqual(result.missing,0)
        self.assertEqual(result.unreadable_ref,0)
        self.assertEqual(result.types_differ+result.ok,result.total())
    def test_diff_replace_directory_with_symlink(self):
        """DataDir.diff() confirms directory replaced with a symlink

        """
        # Make a directory in the reference and a symlink in the copy
        # with the same name
        self.reference.add_dir("bonjour")
        self.example_dir.add_link("bonjour","hello")
        # Do the verification
        ref_dir = DataDir(self.ref)
        result = ref_dir.diff(self.wd)
        self.assertTrue(result.total() > 0)
        self.assertNotEqual(result.ok,result.total())
        self.assertEqual(result.types_differ,1)
        self.assertEqual(result.missing,0)
        self.assertEqual(result.unreadable_ref,0)
        self.assertEqual(result.types_differ+result.ok,result.total())

#######################################################################
# Main program
#######################################################################

if __name__ == "__main__":

    p = optparse.OptionParser(usage="%prog [OPTIONS] DATA_DIR",
                              version="%prog "+__version__,
                              description="Sequence data and analysis curation utility. "
                              "DATA_DIR should be a top-level directory holding "
                              "sequencing data and/or bioinformatic analyses.")
    # Copying and managing
    group = optparse.OptionGroup(p,"Copying and managing",
                                 "Options for copying and managing the data directory")
    p.add_option_group(group)
    group.add_option('--copy-to',action='store',dest='dest_dir',default=None,
                     help='copy DATA_DIR into DEST_DIR using rsync')
    group.add_option('--set-group',action='store',dest='new_group',default=None,
                     help='set the group ownership to NEW_GROUP')
    # Checking and verifying
    group = optparse.OptionGroup(p,"Checking and verifying",
                                 "Options for checking and verifying the data directory "
                                 "contents")
    p.add_option_group(group)
    group.add_option('--verify',action='store',dest='ref_dir',default=None,
                     help='verify DATA_DIR against REF_DIR: report MD5 sums for '
                     'common files, link targets for common links, and any files '
                     'or directories missing from DATA_DIR which are present in '
                     'REF_DIR')
    group.add_option('--diff',action='store',dest='diff_dir',default=None,
                     help='check whether DATA_DIR contains files, directories and '
                     'symlinks in DIFF_DIR. Only lists missing items, not whether '
                     'common items are the same - use --verify for more rigorous '
                     'checking')
    group.add_option('--check-group',action='store',dest='group',default=None,
                     help='check that files are owned by GROUP and have group '
                     'read-write permissions')
    group.add_option('--check-symlinks',action='store_true',dest='check_symlinks',
                     help='check for broken and absolute symlinks')
    group.add_option('--check-temporary',action='store_true',dest='check_temporary',
                     help='check for hidden and temporary data')
    # Acquiring information
    group = optparse.OptionGroup(p,"Acquiring information",
                                 "Options for acquiring information about the data "
                                 "directory and its contents")
    p.add_option_group(group)
    group.add_option('--info',action='store_true',dest='info',
                 help='print general information about DATA_DIR')
    group.add_option('--list-users',action='store_true',dest='list_users',default=None,
                     help='print a list of usernames associated with files')
    group.add_option('--find',action='store',dest='regex_pattern',default=None,
                     help='find files that match REGEX_PATTERN')
    # Debugging
    p.add_option('--debug',action='store_true',dest='debug',
                 help='turn on debugging output')
    
    options,args = p.parse_args()
    if len(args) != 1:
        p.error("Need to supply a data directory")
    data_dir = DataDir(args[0])

    # Debugging output
    if options.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Get information
    if options.info:
        print "Gathering information about %s" % data_dir.dir
        data_dir_size = data_dir.get_size()
        print "Size: %fK (%s)" % (float(data_dir_size)/1024,
                                  bcf_utils.format_file_size(data_dir_size))

    # Find files matching pattern
    if options.regex_pattern is not None:
        print "Looking for files matching '%s' in %s" % (options.regex_pattern,
                                                         data_dir.dir)
        for filen in data_dir.find_files(options.regex_pattern):
            print "%s" % os.path.relpath(filen,data_dir.dir)

    # Rsync to working directory
    if options.dest_dir is not None:
        dest_dir = os.path.abspath(options.dest_dir)
        print "Making a copy of %s under %s" % (data_dir.dir,dest_dir)
        status = data_dir.copy(dest_dir,verbose=True)
        print "Copy completed with status %s" % status

    # Diff against a reference directory
    if options.diff_dir is not None:
        print "Diffing %s against reference directory %s" % (data_dir.dir,
                                                             options.diff_dir)
        result = DataDir(options.diff_dir).diff(data_dir.dir)
        # Print report
        result.report()
        status = 0 if result.total() == result.ok else 1
        print "Diff completed with status %s" % status
        

    # Verify against a reference directory
    if options.ref_dir is not None:
        print "Verifying %s against reference directory %s" % (data_dir.dir,
                                                               options.ref_dir)
        result = DataDir(options.ref_dir).verify(data_dir.dir)
        # Print report
        result.report()
        status = 0 if result.total() == result.ok else 1
        print "Verification completed with status %s" % status

    # Check group and permissions
    if options.group:
        print "Checking group ownership for '%s' in %s" % (data_dir.dir,
                                                           options.group)
        status = 0
        group = options.group
        gid = bcf_utils.get_gid_from_group(group)
        if gid is None:
            logging.error("Unable to locate group '%s' on this system" % group)
            sys.exit(1)
        print "Group '%s' guid = %s" % (group,gid)
        print "** NB links will be ignored **"
        header = "File\tOwner\tGroup\tRW"
        for filen in data_dir.walk():
            f = bcf_utils.PathInfo(filen)
            if f.is_link:
                continue
            wrong_group = (f.gid != gid)
            if wrong_group:
                logging.debug("Wrong group (%s):\t%s" % (f.gid,f.relpath(data_dir.dir)))
            group_read_write = (f.is_group_readable and f.is_group_writable)
            if not group_read_write:
                logging.debug("Not group read/writable:\t%s" %
                              f.relpath(data_dir.dir))
            if wrong_group or not group_read_write:
                status = 1
                group_name = f.group
                owner = f.user
                logging.debug("Owner is %s" % owner)
                if header:
                    print header
                    header = None
                print "%s\t%s\t%s\t%s" % (f.relpath(data_dir.dir),
                                          owner,group_name,
                                          '' if group_read_write else 'No')
        print "Group/permissions check completed with status %s" % status

    # Check for temporary and hidden files/directories
    if options.check_temporary:
        print "Checking for temporary/hidden data in %s" % data_dir.dir
        nfiles = 0
        temporary = []
        for filen in data_dir.walk():
            nfiles += 1
            if os.path.basename(filen).startswith('.'):
                temporary.append(filen)
            elif os.path.basename(filen).find('tmp') > -1:
                temporary.append(filen)
        print "%d files & directories examined" % nfiles
        if temporary:
            print "Found %d temporary/hidden data items:" % len(temporary)
            for filen in temporary:
                print "\t%s (%s)" % (os.path.relpath(filen,data_dir.dir),
                                     bcf_utils.format_file_size(get_size(filen)))
        else:
            print "No temporary or hidden files found"

    # Examine links
    if options.check_symlinks:
        print "Checking symlinks in %s" % data_dir.dir
        nlinks = 0
        broken_links = []
        absolute_links = []
        for link in bcf_utils.links(data_dir.dir):
            nlinks += 1
            symlink = bcf_utils.Symlink(link)
            if symlink.is_broken:
                broken_links.append(link)
            elif symlink.is_absolute:
                absolute_links.append(link)
        print "%d symlinks examined" % nlinks
        if broken_links:
            print "Found %d broken links:" % len(broken_links)
            for link in broken_links:
                print "\t%s -> %s" % (os.path.relpath(link,data_dir.dir),
                                      bcf_utils.Symlink(link).target)
        else:
            print "No broken links"
        if absolute_links:
            print "Found %d absolute links:" % len(absolute_links)
            for link in absolute_links:
                print "\t%s -> %s" % (os.path.relpath(link,data_dir.dir),
                                      bcf_utils.Symlink(link).target)
        else:
            print "No absolute links"

    # Set group
    if options.new_group is not None:
        print "Setting group ownership to '%s' on %s" % (options.new_group,
                                                         data_dir.dir)
        new_group = options.new_group
        gid = bcf_utils.get_gid_from_group(new_group)
        if gid is None:
            logging.error("Unable to locate group '%s' on this system" % new_group)
            sys.exit(1)
        print "Group '%s' guid = %s" % (new_group,gid)
        for filen in data_dir.walk():
            f = bcf_utils.PathInfo(filen)
            if not f.is_link:
                f.chown(group=gid)

    # List users
    if options.list_users:
        print "Collecting list of usernames from %s" % data_dir.dir
        for user in data_dir.users:
            real_name = finger_user(user)
            print "%s\t%s" % (user,'?' if real_name is None else real_name)

