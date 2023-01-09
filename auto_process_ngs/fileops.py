#!/usr/bin/env python
#
#     fileops: single interface for file ops on local and remote systems
#     Copyright (C) University of Manchester 2017,2023 Peter Briggs
#
########################################################################
#
# fileops.py
#
#########################################################################

"""
Utility functions providing a single interface for performing
various file system operations (e.g. make a directory, copy
files etc) transparently on either a local or a remote system.

The following functions perform specific operations directly:

- mkdir: create a directory
- copy: copy a file
- copytree: recursively copy a directory
- set_group: set the group on a file or directory
- unzip: unpack a ZIP archive
- rename: rename (move) a file or directory
- exists: test if a file or directory exists
- isdir: test if a path is a directory
- remove_file: remove (delete) a file
- remove_dir: remove (delete) a directory
- disk_usage: get info on disk usage for a path

These functions generate commands that can be executed e.g.
via a scheduler, to perform the required operations:

- copy_command: generate command to perform copy operation
- copytree_command: generate command to perform recursive copy
  operation
- set_group_command: generate command to perform group set
  operation
- unzip_command: generate command to unpack a ZIP archive
"""

########################################################################
# Imports
#########################################################################

import os
import shutil
import getpass
import psutil
import collections
import logging
import bcftbx.utils as bcftbx_utils
from . import applications
from .utils import Location

# Module specific logger
logger = logging.getLogger(__name__)

########################################################################
# Command execution functions
#########################################################################

def _run_command(cmd,sched=None):
    """
    Run the command, either locally or via a scheduler

    Arguments:
      cmd (Command): command to run
      sched (SimpleScheduler): optional, a scheduler
        to use to run the command
    """
    print("Running %s" % cmd)
    if sched is None:
        retcode,output = cmd.subprocess_check_output()
    else:
        job = sched.submit(cmd)
        job.wait()
        retcode = job.exit_code
    return retcode

def mkdir(newdir,recursive=False):
    """
    Create a directory

    The new directory should be identified using a
    specifier of the form '[[USER@]HOST:]NEWDIR'.

    The parent directories must already exist, unless
    the 'recursive' argument is set (in which case
    all missing parent directories will also be
    created).

    Arguments:
      newdir (str): location of the new directory (can
        be on local or remote system)
      recursive (bool): if True then also create
        missing parent directories (default is not
        to create missing directories)
    """
    newdir = Location(newdir)
    mkdir_cmd = applications.Command('mkdir')
    if recursive:
        mkdir_cmd.add_args('-p')
    mkdir_cmd.add_args(newdir.path)
    if newdir.is_remote:
        # Remote directory
        mkdir_cmd = applications.general.ssh_command(
            newdir.user,
            newdir.server,
            mkdir_cmd.command_line)
    try:
        return _run_command(mkdir_cmd)
    except Exception as ex:
        raise Exception(
            "Exception making directory %s: %s" %
            (newdir,ex))

def copy(src,dest,link=False):
    """
    Copy a file

    Copies a local file to a local or remote location.

    Arguments:
      src (str): local file to copy
      dest (str): destination (file or directory)
        on a local or remote system, identified by
        a specifier of the form '[[USER@]HOST:]DEST'
      link (bool): hard link files instead of
        copying (ignored for remote copies)

    """
    copy_cmd = copy_command(src,dest,link=link)
    try:
        return _run_command(copy_cmd)
    except Exception as ex:
        raise Exception("Exception copying %s to %s: "
                        "%s" % (src,dest,ex))

def copytree(src,dest):
    """
    Recursively copy a local directory tree

    Recursively copies an entire local directory tree
    rooted at 'src', to a local or remote destination
    directory 'dest'.

    Note that if 'dest' already exists then 'src' will
    be copied into it as a subdirectory.

    Arguments:
      src (str): local directory to copy
      dest (str): destination directory)
        on a local or remote system, identified by
        a specifier of the form '[[USER@]HOST:]DEST'
    """
    copytree_cmd = copytree_command(src,dest)
    try:
        return _run_command(copytree_cmd)
    except Exception as ex:
        raise Exception("Exception copying %s to %s: "
                        "%s" % (src,dest,ex))

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
    chmod_cmd = set_group_command(group,path)
    try:
        return _run_command(chmod_cmd)
    except Exception as ex:
        raise Exception(
            "Exception changing group to '%s' for "
            "destination %s: %s" % (group,path,ex))

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
    unzip_cmd = unzip_command(zip_file,dest)
    try:
        return _run_command(unzip_cmd)
    except Exception as ex:
        raise Exception("Failed to unzip %s: %s" %
                        (zip_file,ex))

def rename(src,dst):
    """
    Rename (move) a file or directory

    Arguments:
      src (str): path to file or directory to
        rename
      dst (str): path to rename 'src' to
    """
    src = Location(src)
    dst = Location(dst)
    # Sanity check: if destination is remote then
    # must be on same server as source
    if dst.is_remote:
        if dst.server != src.server:
            raise Exception("Rename: can't rename on different "
                            "servers")
    # Build generic system command
    rename_cmd = applications.Command('mv',
                                      src.path,
                                      dst.path)
    if src.is_remote:
        # Renaming file on remote system
        rename_cmd = applications.general.ssh_command(
            src.user,
            src.server,
            rename_cmd.command_line)
    # Run command and return
    retval,output = rename_cmd.subprocess_check_output()
    return retval

def exists(path):
    """
    Test if a file or directory exists

    Arguments:
      path (str): path to file or directory to
        check existence of

    Returns:
      Boolean: True if file or directory exists,
        False otherwise.
    """
    path = Location(path)
    test_cmd = applications.Command('test',
                                    '-e',
                                    path.path)
    if path.is_remote:
        # Run test on remote system
        test_cmd = applications.general.ssh_command(
            path.user,
            path.server,
            test_cmd.command_line)
    retval,output = test_cmd.subprocess_check_output()
    return (retval == 0)

def isdir(path):
    """
    Test if a path is a directory

    Arguments:
      path (str): path to check

    Returns:
      Boolean: True if path is a directory,
        False otherwise.
    """
    path = Location(path)
    test_cmd = applications.Command('test',
                                    '-d',
                                    path.path)
    if path.is_remote:
        # Run test on remote system
        test_cmd = applications.general.ssh_command(
            path.user,
            path.server,
            test_cmd.command_line)
    retval,output = test_cmd.subprocess_check_output()
    return (retval == 0)

def remove_file(path):
    """
    Remove (delete) a file

    Arguments:
      path (str): path to file to delete

    Returns:
      Integer: zero on success, non-zero on
        failure.
    """
    path = Location(path)
    rm_cmd = applications.Command('rm',
                                  '-f',
                                  path.path)
    if path.is_remote:
        # Run removal on remote system
        rm_cmd = applications.general.ssh_command(
            path.user,
            path.server,
            rm_cmd.command_line)
    retval,output = rm_cmd.subprocess_check_output()
    return retval

def remove_dir(path):
    """
    Remove (delete) a directory

    Arguments:
      path(str): path to directory to delete

    Returns:
      Integer: zero on success, non-zero on
        failure.
    """
    path = Location(path)
    if not isdir(path.path):
        return 1
    rm_cmd = applications.Command('rm',
                                  '-rf',
                                  path.path)
    if path.is_remote:
        # Run removal on remote system
        rm_cmd = applications.general.ssh_command(
            path.user,
            path.server,
            rm_cmd.command_line)
    retval,output = rm_cmd.subprocess_check_output()
    return retval

def disk_usage(path):
    """
    Get disk usage for a path

    Wraps the psutil 'disk_usage' function for local paths,
    and runs the 'df' command for paths on remote systems.
    Raises OSError if the path doesn't exist.

    Arguments:
      path (str): path to directory

    Returns:
      NamedTuple: NamedTuple with fields 'total', 'used',
        'free' and 'percent'.
    """
    path = Location(path)
    if not path.is_remote:
        return psutil.disk_usage(path.path)
    else:
        # Check path exists on remote system
        if not exists(path.path):
            raise OSError("[Errno 2] No such file or directory: '%s'" %
                          path)
        # Run df command on remote system
        # Use --block-size=1 to get same output as
        df_cmd = applications.Command('df',
                                      '--block-size=1',
                                      '--output=size,used,avail,pcent',
                                      path.path)
        df_cmd = applications.general.ssh_command(
            path.user,
            path.server,
            df_cmd.command_line)
        retval,output = df_cmd.subprocess_check_output()
        # Process the output
        # df output looks like e.g.:
        #   1B-blocks        Used       Avail Use%
        # 78729973760 40163758080 34522906624  54%
        fields = output.split('\n')[1].split()
        sdiskusage = collections.namedtuple("sdiskusage",
                                            "total used free percent")
        return sdiskusage(total=int(fields[0]),
                          used=int(fields[1]),
                          free=int(fields[2]),
                          percent=float(fields[3].strip('%')))

########################################################################
# Command generation functions
#########################################################################

def copy_command(src,dest,link=False):
    """
    Generate command to copy a file

    Creates a command which copies a local file to a
    local or remote location.

    Arguments:
      src (str): local file to copy
      dest (str): destination (file or directory)
        on a local or remote system, identified by
        a specifier of the form '[[USER@]HOST:]DEST'
      link (bool): hard link files instead of
        copying (ignored for remote copies)

    Returns:
      Command: Command instance that can be used to
        perform the copy operation.
    """
    dest = Location(dest)
    if not dest.is_remote:
        # Local-to-local copy
        copy_cmd = applications.Command('cp')
        if link:
            copy_cmd.add_args('-l')
        copy_cmd.add_args(src,
                          dest.path)
    if dest.is_remote:
        # Local-to-remote copy
        copy_cmd = applications.general.scp(
                dest.user,
                dest.server,
                src,dest.path)
    return copy_cmd

def copytree_command(src,dest):
    """
    Generate command to recursively copy a directory tree

    Creates a command which recursively copies an entire
    local directory tree rooted at 'src', to a local or
    remote destination directory 'dest'.

    Note that if 'dest' already exists then 'src' will
    be copied into it as a subdirectory.

    Arguments:
      src (str): local directory to copy
      dest (str): destination directory)
        on a local or remote system, identified by
        a specifier of the form '[[USER@]HOST:]DEST'

    Returns:
      Command: Command instance that can be used to
        perform the copy operation.
    """
    dest = Location(dest)
    if not dest.is_remote:
        # Local-to-local copy
        copytree_cmd = applications.Command('cp',
                                            '-a',
                                            '-n',
                                            src,
                                            dest.path)
    else:
        # Local-to-remote copy
        copytree_cmd = applications.general.scp(
            dest.user,
            dest.server,
            src,dest.path,
            recursive=True)
    return copytree_cmd

def set_group_command(group,path,verbose=False,safe=False):
    """
    Generate command to set group on file or directory

    Creates a command which sets the group on a
    file or directory.

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
      verbose (bool): if True then output a
        diagnostic for every file processed
        (default: don't output diagnostics)
      safe (bool): if True then run in 'safe'
        mode i.e. only files owned by user will
        be have group updated (default: try to
        set on all files and directories)

    Returns:
      Command: Command instance that can be used to
        set the group.
    """
    path = Location(path)
    username = path.user
    if username is None:
        username = getpass.getuser()
    chmod_cmd = applications.Command('find',
                                     path.path)
    if safe:
        chmod_cmd.add_args('-user',username)
    chmod_cmd.add_args('!',
                       '-type',
                       'l',
                       '-exec',
                       'chgrp')
    if verbose:
        chmod_cmd.add_args('--verbose')
    chmod_cmd.add_args(group,
                       '{}',
                       '+')
    if path.is_remote:
        # Set group for remote files
        chmod_cmd = applications.general.ssh_command(
            path.user,
            path.server,
            chmod_cmd.command_line)
    return chmod_cmd

def unzip_command(zip_file,dest):
    """
    Generate command to unpack ZIP archive file

    The ZIP archive can be on a local or a remote
    system.

    Arguments:
      zip_file (str): ZIP archive file identified
        using a specifier of the form
        '[[USER@]HOST:]ZIP_FILE'
      dest (str): path to extract the archive
        contents to (on the same system as the ZIP
        archive)

    Returns:
      Command: Command instance that can be used to
        perform the unzip operation.
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
    return unzip_cmd
