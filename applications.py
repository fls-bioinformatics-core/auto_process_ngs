#!/bin/env python
#
#     applications.py: utilities for running command line applications
#     Copyright (C) University of Manchester 2013 Peter Briggs
#
########################################################################
#
# applications.py
#
#########################################################################

"""applications.py

Utility classes and functions for generating and executing command lines
to run various command line applications.

Static classes provide methods for building command lines for various
NGS applications, in the form of 'Command' instances.

For example, to create a Command object representing the command line for
a simple mirroring 'rsync' job:

>>> import applications
>>> rsync = applications.general.rsync('source','target',mirror=True)
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
# Module metadata
#######################################################################

__version__ = "0.0.6"

#######################################################################
# Import modules that this module depends on
#######################################################################

import sys
import re
import subprocess
import logging
import bcf_utils

#######################################################################
# Classes
#######################################################################

class Command:
    """Class for creating and executing command lines

    The Command class is intended to help with building and
    executing command lines.

    For example to create a command line to run the Linux
    'ls' command:

    >>> ls = Command('ls')

    To add some arguments:

    >>> ls.add_args('-l','-tr')

    To see the resulting command line:

    >>> str(ls)

    To run using the subprocess module do:

    >>> ls.run_subprocess()

    """
    def __init__(self,command,*args):
        """Create a new Command instance

        Arguments:
          command: the initial command i.e. program name
          args   : optional, one or more additional command
                   line arguments

        """
        self._cmd = str(command)
        if not args:
            self._args = []
        else:
            self._args = [str(x) for x in args]
    def add_args(self,*args):
        for arg in args:
            self._args.append(str(arg))

    @property
    def command(self):
        """Return the command
        
        """
        return self._cmd

    @property
    def args(self):
        """Return the arguments as a list
        
        """
        return self._args

    @property
    def command_line(self):
        """Return the full command line as a list
        
        """
        command_line = [self.command]
        command_line.extend(self.args)
        return command_line

    def __contains__(self,item):
        """Implement container functionality for 'item in x'

        """
        return (item in self.command_line)

    def __getitem__(self,key):
        """Implement container functionality for 'x[key]'

        """
        return self.command_line[key]

    def __len__(self):
        """Implement container functionality for 'len(x)'

        """
        return len(self.command_line)

    def __repr__(self):
        """Return the full command line as a string

        """
        return ' '.join(self.command_line)

    def run_subprocess(self,log=None,err=None,working_dir=None):
        """Run the command using the subprocess module

        This runs the command using the subprocess.popen() function
        and wais for it to finish.

        Arguments:
          log: optional, name of file to write stdout to (defaults
            to sys.stdout)
          err: optional, name of file to write stderr to (defaults
            to same location as log, if that was specified, or else
            to sys.stderr)
          working_dir: optional, working directory to use (defaults
            to current directory

        Returns:
          Return code from subprocess.popen() call.

        """
        # Deal with output destinations
        if log is None:
            fpout = sys.stdout
        else:
            logging.debug("Writing stdout to %s" % log)
            fpout = open(log,'w')
        if err is None:
            if log is not None:
                fperr = subprocess.STDOUT
            else:
                fperr = sys.stderr
        else:
            fperr = open(err,'w')
            logging.debug("Writing stderr to %s" % err)
        # Execute command and wait for finish
        try:
            p = subprocess.Popen(self.command_line,
                                 cwd=working_dir,stdout=fpout,stderr=fperr)
            returncode = p.wait()
        except KeyboardInterrupt,ex:
            # Handle keyboard interrupt while rsync is running
            logging.warning("KeyboardInterrupt: stopping command subprocess")
            p.kill()
            returncode = -1
        return returncode

class bcl2fastq:
    """Bcl to fastq conversion line applications
 
    Provides static methods to create Command instances for command line
    applications used in bcl to fastq conversion:

    configureBclToFastq

    """

    @staticmethod
    def configureBclToFastq(basecalls_dir,sample_sheet,output_dir="Unaligned",
                            mismatches=None,
                            bases_mask=None,
                            force=False,
                            ignore_missing_bcl=False,
                            ignore_missing_stats=False,
                            ignore_missing_control=False):
        """Generate Command instance for 'configureBclToFastq.pl' script

        Creates a Command instance to run the CASAVA 'configureBclToFastq.pl'
        script (which generates a Makefile to perform the bcl to fastq
        conversion).
        
        Arguments:
          basecalls_dir: path to the top-level directory holding the bcl
            files (typically 'Data/Intensities/Basecalls/' subdirectory)
          sample_sheet: path to the sample sheet file to use
          output_dir: optional, path to the output directory. Defaults to
            'Unaligned'. If this directory already exists then the
            conversion will fail unless the force option is set to True
          mismatches: optional, specify maximum number of mismatched bases
            allowed for matching index sequences during multiplexing.
            Recommended values are zero for indexes shorter than 6 base
            pairs, 1 for indexes of 6 or longer
            (If not specified and bases_mask is supplied then mismatches
            will be derived automatically from the bases mask string)
          bases_mask: optional, specify string indicating how to treat
            each cycle within each read e.g. 'y101,I6,y101'
          force: optional, if True then force overwrite of an existing
            output directory (default is False)
          ignore_missing_bcl: optional, if True then interpret missing bcl
            files as no call (default is False)
          ignore_missing_stats: optional, if True then fill in with zeroes
            when *.stats files are missing (default is False)
          ignore_missing_control: optional, if True then interpret missing
            control files as not-set control bits (default is False)

        Returns:
          Command object.

        """
        configure_cmd = Command('configureBclToFastq.pl',
                                '--input-dir',basecalls_dir,
                                '--output-dir',output_dir,
                                '--sample-sheet',sample_sheet,
                                '--fastq-cluster-count','-1')
        if bases_mask is not None:
            configure_cmd.add_args('--use-bases-mask',bases_mask)
        if mismatches is not None:
            configure_cmd.add_args('--mismatches',mismatches)
        else:
            # Nmismatches not supplied, derive from bases mask
            if bases_mask is not None:
                configure_cmd.add_args('--mismatches',get_nmismatches(bases_mask))
        if force:
            configure_cmd.add_args('--force')
        if ignore_missing_bcl:
            configure_cmd.add_args('--ignore-missing-bcl')
        if ignore_missing_stats:
            configure_cmd.add_args('--ignore-missing-stats')
        if ignore_missing_control:
            configure_cmd.add_args('--ignore-missing-control')
        return configure_cmd

class general:
    """General command line applications (e.g. rsync, make)
 
    Provides static methods to create Command instances for a class of
    'general' command line applications:

    rsync
    make

    """

    @staticmethod
    def rsync(source,target,dry_run=False,mirror=False,chmod=None,
              prune_empty_dirs=False,extra_options=None):
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
        ssh_command = Command('ssh','%s@%s' % (user,server))
        ssh_command.add_args(*cmd)
        return ssh_command

    @staticmethod
    def scp(user,server,source,target):
        """Generate Command instance for 'scp'

        Creates a Command instance to run 'scp' to copy to another system.

        Arguments:
          user: name of the remote user
          server: name of the server
          source: source file on local system
          target: target destination on remote system

        Returns:
          Command object.

        """
        scp_command = Command('scp',source,'%s@%s:%s' % (user,server,target))
        return scp_command
