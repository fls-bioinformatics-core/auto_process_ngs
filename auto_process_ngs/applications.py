#!/usr/bin/env python
#
#     applications.py: utilities for running command line applications
#     Copyright (C) University of Manchester 2013-17,2019-2021 Peter Briggs
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

__version__ = "0.0.10"

#######################################################################
# Import modules that this module depends on
#######################################################################

import sys
import re
import subprocess
import logging
import tempfile
from bcl2fastq.utils import get_nmismatches
from bcftbx.utils import find_program

#######################################################################
# Classes
#######################################################################

class Command(object):
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
        """Append arguments to the command

        """
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

    @property
    def has_exe(self):
        """Check if the command executable exists

        """
        return (find_program(self.command) is not None)

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
        """Run the command using subprocess.Popen

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
          Return code from subprocess.Popen() call.

        """
        # Deal with output destinations
        if log is None:
            fpout = None
        else:
            logging.debug("Writing stdout to %s" % log)
            fpout = open(log,'wt')
        if err is None:
            if log is not None:
                fperr = subprocess.STDOUT
            else:
                fperr = None
        else:
            fperr = open(err,'wt')
            logging.debug("Writing stderr to %s" % err)
        # Execute command and wait for finish
        try:
            p = subprocess.Popen(self.command_line,
                                 cwd=working_dir,stdout=fpout,stderr=fperr)
            returncode = p.wait()
        except KeyboardInterrupt as ex:
            # Handle keyboard interrupt while process is running
            logging.warning("KeyboardInterrupt: stopping command subprocess")
            p.kill()
            returncode = -1
        # Close files before finishing
        if log:
            fpout.close()
        if err:
            fperr.close()
        return returncode

    def subprocess_check_output(self,include_err=True,working_dir=None):
        """Run the command and capture the output

        This runs the command using the subprocess.call() and
        captures the output in a string which is returned to the
        calling subprogram (along with the return code from the
        command).

        Nb it would be better if we could use the subprocess.check_output
        function (but this is not available in Python 2.6).

        Arguments:
          include_err: optional, if True then stderr is included in
            the output (default); otherwise stderr is discarded.
          working_dir: optional, working directory to use (defaults
            to current directory

        Returns:
          Tuple of (returncode,output).

        """
        # Include stderr (or not)
        if include_err:
            stderr = subprocess.STDOUT
        else:
            stderr = None
        # Create a temporary file-like object to capture output
        with tempfile.TemporaryFile(mode='w+t') as ftmp:
            status = subprocess.call(self.command_line,stdout=ftmp,
                                     cwd=working_dir,stderr=stderr)
            # Read the output
            ftmp.seek(0)
            output = ftmp.read()
        return (status,output)

    def make_wrapper_script(self,shell=None,filen=None,fp=None,
                            prologue=None,epilogue=None,
                            quote_spaces=False):
        """Wrap the command in a script

        Returns a string which can be injected into a file and
        run as a script.

        Arguments:
          shell (str): optional, if set then will be written
            to the wrapper script shebang (#!)
          filen (str): optional, if set then wrapper script will
            be written to a file with this path
          fp (File): optional, if set then must be a File-like
            object opened for writing, to which the wrapper script
            will be written
          prologue (str): optional, if set then will be written
            into the script before the command
          epilogue (str): optional, if set then will be written
            into the script after the command
          quote_spaces (str): if True then arguments containing
            whitespace will be wrapped in quotes

        Returns:
          String: the wrapper script contents.
        """
        # Handle quoting of spaces
        if quote_spaces:
            args = []
            for arg in self.command_line:
                if (' ' in str(arg)):
                    args.append("'%s'" % arg)
                else:
                    args.append(str(arg))
        else:
            args = [str(arg) for arg in self.command_line]
        # Build script
        script = []
        if shell is not None:
            script.append("#!%s" % shell)
        if prologue is not None:
            script.append("%s" % prologue)
        script.append(' '.join(args))
        if epilogue is not None:
            script.append("%s" % epilogue)
        script = '\n'.join(script)
        # Write to file
        if fp is not None:
            fp.write(script)
        if filen is not None:
            with open(filen,'wt') as fp:
                fp.write(script)
        return script

class bcl2fastq(object):
    """Bcl to fastq conversion line applications
 
    Provides static methods to create Command instances for command line
    applications used in bcl to fastq conversion:

    configureBclToFastq
    bcl2fastq2

    """

    @staticmethod
    def configureBclToFastq(basecalls_dir,sample_sheet,output_dir="Unaligned",
                            mismatches=None,
                            bases_mask=None,
                            force=False,
                            ignore_missing_bcl=False,
                            ignore_missing_stats=False,
                            ignore_missing_control=False,
                            configureBclToFastq_exe=None):
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
          configureBclToFastq_exe: optional, if set then will be taken
            as the name/path for the 'configureBclToFastq.pl' script

        Returns:
          Command object.

        """
        if configureBclToFastq_exe is None:
            configureBclToFastq_exe = 'configureBclToFastq.pl'
        configure_cmd = Command(configureBclToFastq_exe,
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

    @staticmethod
    def bcl2fastq2(run_dir,sample_sheet,output_dir="Unaligned",
                   mismatches=None,
                   bases_mask=None,
                   ignore_missing_bcl=False,
                   no_lane_splitting=False,
                   minimum_trimmed_read_length=None,
                   mask_short_adapter_reads=None,
                   create_fastq_for_index_reads=False,
                   loading_threads=None,
                   demultiplexing_threads=None,
                   processing_threads=None,
                   writing_threads=None,
                   bcl2fastq_exe=None):
        """
        Generate Command instance for 'bcl2fastq' program (v2.*)

        Creates a Command instance to run the Illumina 'bcl2fastq'
        program (for versions 2.*).

        Arguments:
          run: path to the top-level directory for the run
          sample_sheet: path to the sample sheet file to use
          output_dir: optional, path to the output directory. Defaults to
            'Unaligned'
          mismatches: optional, specify maximum number of mismatched bases
            allowed for matching index sequences during multiplexing.
            Recommended values are zero for indexes shorter than 6 base
            pairs, 1 for indexes of 6 or longer
            (If not specified and bases_mask is supplied then mismatches
            will be derived automatically from the bases mask string)
          bases_mask: optional, specify string indicating how to treat
            each cycle within each read e.g. 'y101,I6,y101'
          ignore_missing_bcl: optional, if True then interpret missing bcl
            files as no call (default is False)
          no_lane_splitting: optional, if True then don't split FASTQ
            files by lane (--no-lane-splitting) (default is False)
          minimum_trimmed_read_length: optional, specify minimum length
            for reads after adapter trimming (shorter reads will be padded
            with Ns to make them long enough)
          mask_short_adapter_reads: optional, specify the minimum length
            of ACGT bases that must be present in a read after adapter
            trimming for it not to be masked completely with Ns.
          create_fastq_for_index_reads: optional, if True then also create
            Fastq files for index reads (default, don't create index read
            Fastqs) (--create-fastq-for-index-reads)
          loading_threads: optional, specify number of threads to use
            for loading bcl data (--loading-threads)
          demultiplexing_threads: optional, specify number of threads to
            use for demultiplexing (--demultiplexing-threads)
          processing_threads: optional, specify number of threads to
            use for processing (--processing-threads)
          writing_threads: optional, specify number of threads to
            use for writing FASTQ data (--writing-threads)
          bcl2fastq_exe: optional, if set then specifies the name/path
            of the bcl2fastq executable to use

        Returns:
          Command object.

        """
        if bcl2fastq_exe is None:
            bcl2fastq_exe = 'bcl2fastq'
        bcl2fastq_cmd = Command(bcl2fastq_exe,
                                '--runfolder-dir',run_dir,
                                '--output-dir',output_dir,
                                '--sample-sheet',sample_sheet)
        if bases_mask is not None:
            bcl2fastq_cmd.add_args('--use-bases-mask',bases_mask)
        if mismatches is not None:
            bcl2fastq_cmd.add_args('--barcode-mismatches',mismatches)
        else:
            # Nmismatches not supplied, derive from bases mask
            if bases_mask is not None:
                bcl2fastq_cmd.add_args('--barcode-mismatches',
                                       get_nmismatches(bases_mask))
        if ignore_missing_bcl:
            bcl2fastq_cmd.add_args('--ignore-missing-bcls')
        if no_lane_splitting:
            bcl2fastq_cmd.add_args('--no-lane-splitting')
        if minimum_trimmed_read_length is not None:
            bcl2fastq_cmd.add_args('--minimum-trimmed-read-length',
                                   minimum_trimmed_read_length)
        if mask_short_adapter_reads is not None:
            bcl2fastq_cmd.add_args('--mask-short-adapter-reads',
                                   mask_short_adapter_reads)
        if create_fastq_for_index_reads:
            bcl2fastq_cmd.add_args('--create-fastq-for-index-read')
        if loading_threads is not None:
            bcl2fastq_cmd.add_args('-r',loading_threads)
        if demultiplexing_threads is not None:
            bcl2fastq_cmd.add_args('-d',demultiplexing_threads)
        if processing_threads is not None:
            bcl2fastq_cmd.add_args('-p',processing_threads)
        if  writing_threads is not None:
            bcl2fastq_cmd.add_args('-w',writing_threads)
        return bcl2fastq_cmd

class general(object):
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
                             '%s@%s:%s' % (user,server,target))
        return scp_command
