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

__version__ = "0.0.1"

#######################################################################
# Import modules that this module depends on
#######################################################################

import sys
import os
import grp
import stat
import re
import optparse
import logging
import subprocess
import bcf_utils

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
            raise OSError, "'%s': not a directory"
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
        """
        """
        rsync = ['rsync','-rltDE','--chmod=u+rwX,g+rwX,o-w',self.dir,target]
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
        for dirpath,dirnames,filenames in os.walk(self.dir):
            for d in dirnames:
                yield os.path.join(dirpath,d)
            for f in filenames:
                yield os.path.join(dirpath,f)

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
    p.add_option('--copy-to',action='store',dest='dest_dir',default=None,
                 help='copy DATA_DIR into DEST_DIR using rsync')
    p.add_option('--check-group',action='store',dest='check_group',default=None,
                 help='check that files are owned by CHECK_GROUP and have group '
                 'read-write permissions')
    p.add_option('--check-symlinks',action='store_true',dest='check_symlinks',
                 help='check for broken symlinks')
    p.add_option('--check-temporary',action='store_true',dest='check_temporary',
                 help='check for hidden and temporary data')
    p.add_option('--set-group',action='store',dest='new_group',default=None,
                 help='set the group ownership to NEW_GROUP')
    p.add_option('--debug',action='store_true',dest='debug',
                 help='turn on debugging output')
    
    options,args = p.parse_args()
    if len(args) != 1:
        p.error("Need to supply a data directory")
    data_dir = DataDir(args[0])

    # Debugging output
    if options.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Get initial size
    if options.info:
        print "Gathering information about %s" % data_dir.dir
        data_dir_size = data_dir.get_size()
        print "Size: %sK (%s)" % (data_dir_size,
                                  bcf_utils.format_file_size(data_dir_size))

    # Rsync to working directory
    if options.dest_dir is not None:
        dest_dir = os.path.abspath(options.dest_dir)
        print "Making a copy of %s under %s" % (data_dir.dir,dest_dir)
        data_dir.copy(dest_dir)

    # Check group and permissions
    if options.check_group:
        print "Checking group ownership for '%s' in %s" % (data_dir.dir,
                                                           options.check_group)
        group = options.check_group
        try:
            gid = grp.getgrnam(group).gr_gid
        except KeyError,ex:
            logging.error("Unable to locate group '%s' on this system" % group)
            sys.exit(1)
        print "Group '%s' guid = %s" % (group,gid)
        for filen in data_dir.walk:
            st = os.stat(filen)
            if st.st_gid != gid:
                print "Wrong group:\t%s" % os.path.relpath(filen,data_dir.dir)
            if not ((st.st_mode & stat.S_IRGRP) and (st.st_mode & stat.S_IWGRP)):
                print "Not group read/writable:\t%s" % os.path.relpath(filen,data_dir.dir)

    # Check for temporary and hidden files/directories
    if options.check_temporary:
        print "Checking for temporary/hidden data in %s" % data_dir.dir
        for filen in data_dir.walk:
            if os.path.basename(filen).startswith('.'):
                print "\t%s (%s)" % (os.path.relpath(filen,data_dir.dir),
                                     bcf_utils.format_file_size(get_size(filen)))
            elif os.path.basename(filen).find('tmp') > -1:
                print "\t%s (%s)" % (os.path.relpath(filen,data_dir.dir),
                                     bcf_utils.format_file_size(get_size(filen)))

    # Examine links
    if options.check_symlinks:
        print "Checking symlinks in %s" % data_dir.dir
        broken_links = []
        absolute_links = []
        for link in data_dir.symlinks:
            filen = os.path.realpath(link)
            if not os.path.exists(filen):
                broken_links.append(link)
            elif os.path.isabs(filen):
                absolute_links.append(link)
        if broken_links:
            print "Broken links:"
            for link in broken_links:
                print "\t%s -> %s" % (os.path.relpath(link,data_dir.dir),
                                      os.path.realpath(link))
        else:
            print "No broken links"
        if absolute_links:
            print "Absolute links:"
            for link in absolute_links:
                print "\t%s -> %s" % (os.path.relpath(link,data_dir.dir),
                                      os.path.realpath(link))
        else:
            print "No absolute links"

    # Set group
    if options.new_group is not None:
        # See http://stackoverflow.com/questions/5994840/how-to-change-the-user-and-group-permissions-for-a-directory-by-name
        raise NotImplementedError
