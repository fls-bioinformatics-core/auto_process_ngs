#!/usr/bin/env python
#
#     applications.py: utilities for running command line applications
#     Copyright (C) University of Manchester 2013-2025 Peter Briggs
#
########################################################################
#
# applications.py
#
#########################################################################

"""
Utilities for generating and executing command lines to run various
command line applications.

Provides the following classes with static methods which will return a
``Command`` instance for running a specific command line application:

* ``general``

For example, to create a ``Command`` object representing the command line
for a simple mirroring ``rsync`` job:

>>> rsync = general.rsync('source','target',mirror=True)
>>> rsync
rsync -av --delete-after source target
>>> rsync.command_line
['rsync', '-av', '--delete-after', 'source', 'target']

The resulting command line string or list can be fed to another function
or class, or it can be executed directly via the subprocess module using
the 'run_subprocess' method of the Command object, e.g:

>>> rsync.run_subprocess()
"""

#######################################################################
# Import modules that this module depends on
#######################################################################

import re
from .command import Command

#######################################################################
# Classes
#######################################################################


class general:
    """General command line applications (e.g. rsync, make)
 
    Provides static methods to create Command instances for a class of
    'general' command line applications:

    * rsync
    * make
    * ssh_command
    * scp
    """

    @staticmethod
    def rsync(source,target,dry_run=False,mirror=False,chmod=None,
              prune_empty_dirs=False,escape_spaces=True,
              extra_options=None):
        """Generate Command instance for 'rsync' command

        Create a Command instance to run the 'rsync' command line, to
        recursively copy/sync one directory (the 'source') into another
        (the 'target').

        The target can be a local directory or on a remote system (in which
        case it should be qualified with a user and hostname i.e.
        'user@hostname:target').

        Arguments:
          source: the directory being copied/sync'ed
          target: the directory the source will be copied into
          dry_run: run rsync using --dry-run option i.e. no files will
            be copied/sync'ed, just reported
          mirror: if True then run rsync in 'mirror' mode i.e. with
            --delete-after option (to remove files from the target
            that have also been removed from the source)
          chmod: optional, mode specification to be applied to the copied
            files e.g. chmod='u+rwX,g+rwX,o-w
          prune_empty_dirs: optional, don't include empty target
            directories i.e. -m option
          escape_spaces: optional, if True then add '\' in front of
            spaces in file name paths to escape them (default)
          extra_options: optional, a list of additional rsync options to be
            added to the command (e.g. --include and --exclude filter
            patterns)

        Returns:
          Command object.

        """
        # Collect rsync command options
        rsync_cmd = Command('rsync','-av')
        # Dry run mode
        if dry_run:
            rsync_cmd.add_args('--dry-run')
        # Check for remote source or target (requires ssh protocol)
        if re.compile(r'^([^@]*@)?[^:]*:').match(target) or \
           re.compile(r'^([^@]*@)?[^:]*:').match(source):
            rsync_cmd.add_args('-e','ssh')
        # Mirroring
        if mirror:
            rsync_cmd.add_args('--delete-after')
        # Don't include empty target directories
        if prune_empty_dirs:
            rsync_cmd.add_args('-m')
        # Set mode of target files and directories
        if chmod is not None:
            rsync_cmd.add_args('--chmod=%s' % chmod)
        # Additional options
        if extra_options is not None:
            rsync_cmd.add_args(*extra_options)
        # Escape spaces in source and target, if required
        if escape_spaces:
            if ' ' in source:
                source = source.replace(' ','\ ')
            if ' ' in target:
                target = target.replace(' ','\ ')
        # Make rsync command
        rsync_cmd.add_args(source,target)
        return rsync_cmd

    @staticmethod
    def make(makefile=None,working_dir=None,nprocessors=None):
        """Generate Command instance for 'make' command

        Creates a Command instance to run 'make'.
        
        Arguments:
          makefile: optional, name of input Makefile (-f)
          working_dir: optional, specify the working directory to change to
            (-C)
          nprocessors: optional, specify number of processors to use (-j)

        Returns:
          Command object.

        """
        make_cmd = Command('make')
        if working_dir is not None:
            make_cmd.add_args('-C',working_dir)
        if nprocessors is not None:
            make_cmd.add_args('-j',nprocessors)
        if makefile is not None:
            make_cmd.add_args('-f',makefile)
        return make_cmd

    @staticmethod
    def ssh_command(user,server,cmd):
        """Generate Command instance for 'ssh' to execute a remote command

        Creates a Command instance to run 'ssh ... COMMAND'.

        Arguments:
          user: name of the remote user
          server: name of the server
          cmd: command to execute on the server via ssh

        Returns:
          Command object.

        """
        ssh_command = Command('ssh','%s%s' % ('%s@' % user if user else '',
                                              server))
        ssh_command.add_args(*cmd)
        return ssh_command

    @staticmethod
    def scp(user,server,source,target,recursive=False):
        """Generate Command instance for 'scp'

        Creates a Command instance to run 'scp' to copy to another system.

        Arguments:
          user: name of the remote user
          server: name of the server
          source: source file on local system
          target: target destination on remote system
          recursive: optional, if True then copy source
            recursively (i.e. specify the '-r' option)

        Returns:
          Command object.

        """
        scp_command = Command('scp')
        if recursive:
            scp_command.add_args('-r')
        scp_command.add_args(source,
                             '%s%s:%s' % ('%s@' % user if user else '',
                                          server,target))
        return scp_command
