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

__version__ = "0.0.13"

#######################################################################
# Import modules that this module depends on
#######################################################################

import sys
import os
import pwd
import grp
import stat
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
        for f in self.walk:
            size += os.lstat(f).st_size
        return int(float(size)/1024)
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

        -p (preserve permissions) is also needed to make the --chmod
        options stick (see rsync manpage).

        Arguments:
          target : target directory to make the copy under
          dry_run: optional, if True then run rsync with --dry-run option
            so no files are copied.

        """
        rsync = ['rsync','-rpltDE','--chmod=u+rwX,g+rwX,o-w',self.dir,target]
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
    @property
    def walk(self):
        """
        """
        yield self.dir
        for dirpath,dirnames,filenames in os.walk(self.dir):
            for d in dirnames:
                yield os.path.join(dirpath,d)
            for f in filenames:
                yield os.path.join(dirpath,f)

class Md5Checker:
    def __init__(self,reference,target):
        self.reference = reference
        self.target = target
        self._passed = []
        self._bad_links = []
        self._missing_files = []
        self._failed_md5sums = []
        self._n_files = 0

    def compare(self,verbose=True):
        """
        """
        ref_dir = self.reference.dir
        tgt_dir = self.target.dir
        for ref_filen in self.reference.files:
            self._n_files += 1
            filen = os.path.join(tgt_dir,os.path.relpath(ref_filen,ref_dir))
            rel_filen = os.path.relpath(filen,tgt_dir)
            if not os.path.exists(filen):
                if os.path.islink(filen):
                    # Deal with bad links
                    if os.readlink(filen) == os.readlink(ref_filen):
                        self.report("OK_LINK",rel_filen,verbose)
                    else:
                        self.report("BAD_LINK",rel_filen,verbose)
                        self._bad_links.append(rel_filen)
                else:
                    self.report("MISSING",rel_filen,verbose)
                    self._missing_files.append(rel_filen)
            else:
                # Compute and compare checksums
                ref_chksum = Md5sum.md5sum(ref_filen)
                chksum = Md5sum.md5sum(filen)
                if chksum != ref_chksum:
                    self.report("FAILED",rel_filen,verbose)
                    self._failed_md5sums.append(rel_filen)
                else:
                    self.report("OK",rel_filen,verbose)

    def report(self,msg,filen,verbose):
        if verbose: print "%s\t%s" % (msg,filen)

    @property
    def n_examined(self):
        """Number of files examined
        """
        return self._n_files
    @property
    def failed_md5sums(self):
        """List of files that failed MD5 sum checks
        """
        return self._failed_md5sums
    @property
    def missing_files(self):
        """List of files missing from the target directory
        """
        return self._missing_files
    @property
    def bad_links(self):
        """List of bad links
        """
        return self._bad_links
    @property
    def n_passed(self):
        """Number of files passing checks
        """
        return self.n_examined \
            - len(self.failed_md5sums) \
            - len(self.missing_files) \
            - len(self.bad_links)

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
    p = subprocess.Popen(['finger',name], stdout=subprocess.PIPE)
    finger, err = p.communicate()
    real_name = finger.split('\n')[0].split(':')[-1].strip()
    return real_name

#######################################################################
# Tests
#######################################################################

import unittest
import tempfile
import shutil
import filecmp

class WorkingDir:
    """Class for building & populating arbitrary directories for tests
    """
    def __init__(self,parent_dir=None):
        """Create a new TestDir instance
        """
        self.dirname = tempfile.mkdtemp(dir=parent_dir)
        logging.debug("Working dir: %s" % self.dirname)
    def rm(self):
        """Remove the entire working directory
        """
        shutil.rmtree(self.dirname)
    def mk_dir(self,*args):
        """Make a subdirectory or series of subdirectories
        """
        subdir = self.dirname
        for arg in args:
            subdir = os.path.join(subdir,arg)
            logging.debug("Making subdir: %s" % subdir)
            try:
                os.mkdir(subdir)
            except OSError,ex:
                if os.path.isdir(subdir): continue
                raise ex
        return subdir
    def rm_dir(self,*args):
        """Remove a subdirectory
        """
        subdir = os.path.join(self.dirname,*args)
        logging.debug("Removing subdir: %s" % subdir)
        shutil.rmtree(subdir)
    def mk_file(self,filen,text=None):
        """Make a file in the working directory
        """
        filen = os.path.join(self.dirname,filen)
        logging.debug("Making filen: %s" % filen)
        fp = open(filen,'w')
        if text is not None:
            fp.write(text)
        fp.close()
        return filen
    def rm_file(self,filen):
        """Remove a file
        """
        filen = os.path.join(self.dirname,filen)
        os.remove(filen)
    def mk_link(self,source,target):
        """Make a symbolic link
        """
        source = os.path.join(self.dirname,source)
        target = os.path.join(self.dirname,target)
        os.symlink(source,target)
        return target

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

def make_languages_data_dir(wd):
    wd.mk_dir('languages')
    wd.mk_file('languages/hello')
    wd.mk_file('languages/goodbye')
    wd.mk_dir('languages/spanish')
    wd.mk_file('languages/spanish/hola')
    wd.mk_file('languages/spanish/adios',text="This means 'goodbye'")
    wd.mk_dir('languages/welsh')
    wd.mk_dir('languages/welsh/north_wales')
    wd.mk_file('languages/welsh/north_wales/maen_ddrwg_gen_i')
    wd.mk_dir('languages/welsh/south_wales')
    wd.mk_file('languages/welsh/south_wales/maen_flin_da_fi')
    return wd

class TestDataDirCopy(unittest.TestCase):
    """Tests for DataDir class 'copy' functionality
    """
    def setUp(self):
        # Make a test data directory structure
        self.wd = WorkingDir()
        make_languages_data_dir(self.wd)
        # Make a temporary destination
        self.dest = WorkingDir()
    def tearDown(self):
        # Remove the test data directory
        self.wd.rm()
        self.dest.rm()
    def test_copy(self):
        """Check copying functionality
        """
        source_dirname = os.path.join(self.wd.dirname,'languages')
        data_dir = DataDir(source_dirname)
        data_dir.copy(self.dest.dirname)
        target_dirname = os.path.join(self.dest.dirname,'languages')
        self.assertTrue(cmp_dirs(source_dirname,target_dirname))
        self.assertTrue(cmp_dirs(target_dirname,source_dirname))

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
        try:
            gid = grp.getgrnam(group).gr_gid
        except KeyError,ex:
            logging.error("Unable to locate group '%s' on this system" % group)
            sys.exit(1)
        print "Group '%s' guid = %s" % (group,gid)
        print "** NB links will be ignored **"
        header = "File\tOwner\tGroup\tRW"
        for filen in data_dir.walk:
            if os.path.islink(filen):
                continue
            st = os.lstat(filen)
            wrong_group = st.st_gid != gid
            if wrong_group:
                logging.debug("Wrong group (%s):\t%s" % (st.st_gid,
                                                         os.path.relpath(filen,data_dir.dir)))
            group_read_write = (st.st_mode & stat.S_IRGRP) and (st.st_mode & stat.S_IWGRP)
            if not group_read_write:
                logging.debug("Not group read/writable:\t%s" %
                              os.path.relpath(filen,data_dir.dir))
            if wrong_group or not group_read_write:
                try:
                    group_name = grp.getgrgid(st.st_gid).gr_name
                except KeyError:
                    group_name = st.st_gid
                try:
                    owner = pwd.getpwuid(st.st_uid).pw_name
                except KeyError:
                    owner = st.st_uid
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
        for filen in data_dir.walk:
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
        for link in data_dir.symlinks:
            nlinks += 1
            if not os.path.exists(os.path.realpath(link)):
                broken_links.append(link)
            elif os.path.isabs(os.readlink(link)):
                absolute_links.append(link)
        print "%d symlinks examined" % nlinks
        if broken_links:
            print "Found %d broken links:" % len(broken_links)
            for link in broken_links:
                print "\t%s -> %s" % (os.path.relpath(link,data_dir.dir),
                                      os.readlink(link))
        else:
            print "No broken links"
        if absolute_links:
            print "Found %d absolute links:" % len(absolute_links)
            for link in absolute_links:
                print "\t%s -> %s" % (os.path.relpath(link,data_dir.dir),
                                      os.readlink(link))
        else:
            print "No absolute links"

    # Do MD5 checksum comparison with a reference directory
    if options.ref_dir is not None:
        ref_dir = DataDir(options.ref_dir)
        print "Comparing %s with reference data directory %s" % (data_dir.dir,ref_dir.dir)
        md5check = Md5Checker(ref_dir,data_dir)
        md5check.compare()
        print "Summary"
        print "\t%d files examined" % md5check.n_examined
        print "\tFAILED\t%d" % len(md5check.failed_md5sums)
        print "\tMISSING\t%d" % len(md5check.missing_files)
        print "\tBAD_LINKS\t%d" % len(md5check.bad_links)
        print "\tOK\t%d" % md5check.n_passed

    # Set group
    if options.new_group is not None:
        print "Setting group ownership to '%s' on %s" % (options.new_group,
                                                         data_dir.dir)
        new_group = options.new_group
        try:
            gid = grp.getgrnam(new_group).gr_gid
        except KeyError,ex:
            logging.error("Unable to locate group '%s' on this system" % new_group)
            sys.exit(1)
        print "Group '%s' guid = %s" % (new_group,gid)
        for filen in data_dir.walk:
            if not os.path.islink(filen):
                os.chown(filen,-1,gid)

    # List users
    if options.list_users:
        print "Collecting list of usernames from %s" % data_dir.dir
        users = []
        for f in data_dir.walk:
            st = os.lstat(f)
            try:
                user = pwd.getpwuid(st.st_uid).pw_name
            except KeyError:
                user = st.st_uid
            if user not in users:
                users.append(user)
        users.sort()
        for user in users:
            real_name = finger_user(user)
            print "%s\t%s" % (user,real_name)

