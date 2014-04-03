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

__version__ = "0.0.17"

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
    @property
    def platform(self):
        """
        """
        raise NotImplementedError
    def get_size(self):
        """
        """
        size = 0
        for f in self.walk():
            size += os.lstat(f).st_size
        return int(float(size)/1024)
    def walk(self,include_dirs=True):
        """Traverse the directory, subdirectories and files
        
        Arguments:
          include_dirs: if True then yield directories as well
            as files (default)
        
        """
        if include_dirs:
            yield self.dir
        for dirpath,dirnames,filenames in os.walk(self.dir):
            if include_dirs:
                for d in dirnames:
                    yield os.path.join(dirpath,d)
            for f in filenames:
                yield os.path.join(dirpath,f)
    def copy(self,target,dry_run=False):
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

        """
        rsync = ['rsync','-rplotDE','--chmod=u+rwX,g+rwX,o-w',self.dir,target]
        if dry_run:
            rsync.insert(1,'--dry-run')
        logging.debug("Rsync command: %s" % rsync)
        status = subprocess.call(rsync)
        logging.debug("Rsync exit status: %s" % status)
        if status:
            raise Exception, "rsync exit code %s (failure)" % status
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

def cmp_dirs(dir1,dir2):
    for dirpath,dirnames,filenames in os.walk(dir1):
        for d in dirnames:
            d2 = os.path.join(dir2,
                              os.path.relpath(os.path.join(dirpath,d),dir1))
            if not os.path.isdir(d2):
                return False
        for f in filenames:
            f2 = os.path.join(dir2,
                              os.path.relpath(os.path.join(dirpath,f),dir1))
            if not os.path.isfile(f2):
                return False
    return True

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
        # Remove the test data directory (and copy)
        self.example_dir.delete_directory()
        TestUtils.remove_dir(self.dest)
    def test_copy(self):
        """DataDir.copy() correctly copies directory

        """
        data_dir = DataDir(self.wd)
        data_dir.copy(self.dest)
        target = os.path.join(self.dest,os.path.basename(self.wd))
        print "%s" % self.wd
        print "%s" % self.dest
        self.assertTrue(cmp_dirs(data_dir.dir,target))
        self.assertTrue(cmp_dirs(target,data_dir.dir))

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
        

#######################################################################
# Main program
#######################################################################

if __name__ == "__main__":

    p = optparse.OptionParser(usage="%prog [OPTIONS] DATA_DIR",
                              version="%prog "+__version__,
                              description="Sequence data and analysis curation utility. "
                              "DATA_DIR should be a top-level directory holding "
                              "sequencing data and/or bioinformatic analyses.")
    p.add_option('--info',action='store_true',dest='info',
                 help='print information about DATA_DIR')
    p.add_option('--find',action='store',dest='regex_pattern',default=None,
                 help='find files that match REGEX_PATTERN')
    p.add_option('--copy-to',action='store',dest='dest_dir',default=None,
                 help='copy DATA_DIR into DEST_DIR using rsync')
    p.add_option('--check-group',action='store',dest='group',default=None,
                 help='check that files are owned by GROUP and have group '
                 'read-write permissions')
    p.add_option('--check-symlinks',action='store_true',dest='check_symlinks',
                 help='check for broken and absolute symlinks')
    p.add_option('--check-temporary',action='store_true',dest='check_temporary',
                 help='check for hidden and temporary data')
    p.add_option('--md5diff',action='store',dest='ref_dir',default=None,
                 help='compare DATA_DIR against REF_DIR using MD5 sums')
    p.add_option('--set-group',action='store',dest='new_group',default=None,
                 help='set the group ownership to NEW_GROUP')
    p.add_option('--list-users',action='store_true',dest='list_users',default=None,
                 help='print a list of usernames associated with files')
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
        data_dir.copy(dest_dir)

    # Check group and permissions
    if options.group:
        print "Checking group ownership for '%s' in %s" % (data_dir.dir,
                                                           options.group)
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
                logging.debug("Wrong group (%s):\t%s" % (f.gid,
                                                         os.path.relpath(filen,data_dir.dir)))
            group_read_write = (f.is_group_readable and f.is_group_writable)
            if not group_read_write:
                logging.debug("Not group read/writable:\t%s" %
                              os.path.relpath(filen,data_dir.dir))
            if wrong_group or not group_read_write:
                group_name = f.group
                owner = f.user
                logging.debug("Owner is %s" % owner)
                if header:
                    print header
                    header = None
                print "%s\t%s\t%s\t%s" % (os.path.relpath(filen,data_dir.dir),
                                          owner,group_name,
                                          '' if group_read_write else 'No')

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

    # Do MD5 checksum comparison with a reference directory
    if options.ref_dir is not None:
        ref_dir = DataDir(options.ref_dir)
        print "Comparing %s with reference data directory %s" % (data_dir.dir,ref_dir.dir)
        reporter = Md5sum.Md5CheckReporter(
            Md5sum.Md5Checker.md5cmp_dirs(ref_dir.dir,data_dir.dir,
                                          links=Md5sum.Md5Checker.IGNORE_LINKS),
            verbose=True)
        # Summarise
        reporter.summary()

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
            if not os.path.islink(filen):
                os.chown(filen,-1,gid)

    # List users
    if options.list_users:
        print "Collecting list of usernames from %s" % data_dir.dir
        users = []
        for f in data_dir.walk():
            user = bcf_utils.PathInfo(f).user
            if user not in users:
                users.append(user)
        users.sort()
        for user in users:
            real_name = finger_user(user)
            print "%s\t%s" % (user,'?' if real_name is None else real_name)

