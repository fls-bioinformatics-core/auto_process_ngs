#!/usr/bin/env python
#
#     fileops: single interface for file ops on local and remote systems
#     Copyright (C) University of Manchester 2017 Peter Briggs
#
########################################################################
#
# fileops.py
#
#########################################################################

"""
fileops

Utility functions providing a single interface for performing
various file system operations (e.g. make a directory, copy
files etc) trasparently on either a local or a remote system.

Classes:

- Location: extracts information from a location specifier

Functions:

- mkdir: create a directory
- copy: copy a file
- set_group: set the group on a file or directory
- unzip: unpack a ZIP archive
"""

########################################################################
# Imports
#########################################################################

import os
import shutil
import logging
import bcftbx.utils as bcftbx_utils
import applications
from utils import split_user_host_dir

# Module specific logger
logger = logging.getLogger(__name__)

########################################################################
# Classes
#########################################################################

class Location(object):
    """
    Class for examining a file-system location specifier

    A location specifier can be a local or a remote file or
    directory. The general form is:

    ``[[user@]server:]path``

    For a local location, only the 'path' component needs to
    be supplied.

    For a remote location, 'server' and 'path' must be
    supplied, while 'user' is optional.

    The following properties are available:

    - user: the user name (or None if not specified)
    - server: the server name (or None if not specified)
    - path: the path component
    - is_remote: True if the location is remote, False if it
      is local
    """
    def __init__(self,location):
        """
        Create a new Location instance

        Arguments:
          location (str): location specifer of the form
            '[[user@]server:]path'
        """
        self._location = location
        self._user,self._server,self._path = split_user_host_dir(
            self._location)
    @property
    def user(self):
        """
        Return 'user' part of '[[user@]server:]path'
        """
        return self._user
    @property
    def server(self):
        """
        Return 'server' part of '[[user@]server:]path'
        """
        return self._server
    @property
    def path(self):
        """
        Return 'path' part of '[[user@]server:]path'
        """
        return self._path
    @property
    def is_remote(self):
        """
        Check if location is on a remote server
        """
        return (self._server is not None)
    def __repr__(self):
        return self._location

########################################################################
# Functions
#########################################################################

def mkdir(newdir):
    """
    Create a directory

    The new directory should be identified using a
    specifier of the form '[[USER@]HOST:]NEWDIR'.

    Arguments:
      newdir (str): location of the new directory (can
        be on local or remote system)
    """
    newdir = Location(newdir)
    if not newdir.is_remote:
        # Local directory
        bcftbx_utils.mkdir(newdir.path)
    else:
        # Remote directory
        try:
            mkdir_cmd = applications.general.ssh_command(
                newdir.user,
                newdir.server,
                ('mkdir',newdir.path))
            print "Running %s" % mkdir_cmd
            mkdir_cmd.run_subprocess()
        except Exception as ex:
            raise Exception(
                "Exception making remote directory %s: %s" %
                (newdir,ex))

def copy(src,dest):
    """
    Copy a file

    Copies a local file to a local or remote location.

    Arguments:
      src (str): local file to copy
      dest (str): destination (file or directory)
        on a local or remote system, identified by
        a specifier of the form '[[USER@]HOST:]DEST'
    """
    dest = Location(dest)
    if not dest.is_remote:
        # Local copy
        shutil.copy(src,dest.path)
    else:
        try:
            # Remote copy
            copy_cmd = applications.general.scp(
                dest.user,
                dest.server,
                src,dest.path)
            print "Running %s" % copy_cmd
            copy_cmd.run_subprocess()
        except Exception as ex:
            raise Exception("Exception copying %s to remote "
                            "destination %s: %s" %
                            (src,dest,ex))

def set_group(group,path):
    """
    Set the group for a file or directory

    'path' can be a file or directory on a local
    or remote system; if it is a directory then it
    will operate recursively i.e. all subdirectories
    and files will also have their group changed.

    'group' is the name of the group to change
    ownership to (must exist on the target system).

    Arguments:
      group (str): name of the new group
      path (str): path to the file or directory
        to change the group of, identified by a
        specifier of the form '[[USER@]HOST:]PATH'
    """
    path = Location(path)
    if not path.is_remote:
        # Set the group for local files
        gid = bcftbx_utils.get_gid_from_group(group)
        if gid is None:
            raise Exception("Failed to get gid for group '%s'" % group)
        for f in bcftbx_utils.walk(path.path,include_dirs=True):
            logger.debug("Updating group for %s" % f)
            os.lchown(f,-1,gid)
    else:
        try:
            # Set group for remote files
            chmod_cmd = applications.general.ssh_command(
                path.user,
                path.server,
                ('chgrp','-R',group,path.path))
            print "Running %s" % chmod_cmd
            chmod_cmd.run_subprocess()
        except Exception as ex:
            raise Exception(
                "Exception changing group to '%s' on remote "
                "destination %s: %s" %
                (group,path,ex))

def unzip(zip_file,dest):
    """
    Unpack ZIP archive file on local or remote system

    Arguments:
      zip_file (str): ZIP archive file identified
        using a specifier of the form
        '[[USER@]HOST:]ZIP_FILE'
      dest (str): path to extract the archive
        contents to (on the same system as the ZIP
        archive)
    """
    zip_file = Location(zip_file)
    unzip_cmd = applications.Command('unzip',
                                     '-q',
                                     '-o',
                                     '-d',dest,
                                     zip_file.path)
    if zip_file.is_remote:
        # Wrap in ssh command for remote zip file
        unzip_cmd = applications.general.ssh_command(
            zip_file.user,
            zip_file.server,
            unzip_cmd.command_line)
    try:
        print "Running %s" % unzip_cmd
        unzip_cmd.run_subprocess()
    except Exception as ex:
        raise Exception("Failed to unzip %s: %s" %
                        (zip_file,ex))
