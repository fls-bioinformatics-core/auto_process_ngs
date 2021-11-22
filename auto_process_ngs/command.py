#!/usr/bin/env python
#
#     command.py: utilities for running command line applications
#     Copyright (C) University of Manchester 2021 Peter Briggs
#
########################################################################
#
# command.py
#
#########################################################################

"""
Provides a single utility class `Command` which can be used
to build command lines to execute applications.
"""

#######################################################################
# Imports
#######################################################################

import subprocess
import logging
from bcftbx.utils import find_program

# Module specific logger
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

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
            logger.debug("Writing stdout to %s" % log)
            fpout = open(log,'wt')
        if err is None:
            if log is not None:
                fperr = subprocess.STDOUT
            else:
                fperr = None
        else:
            fperr = open(err,'wt')
            logger.debug("Writing stderr to %s" % err)
        # Execute command and wait for finish
        try:
            p = subprocess.Popen(self.command_line,
                                 cwd=working_dir,
                                 stdout=fpout,
                                 stderr=fperr)
            returncode = p.wait()
        except KeyboardInterrupt as ex:
            # Handle keyboard interrupt while process is running
            logger.warning("KeyboardInterrupt: stopping command "
                           "subprocess")
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

        This runs the command using ``subprocess.check_output`` and
        returns the output (along with the return code from the
        command).

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
        try:
            output = subprocess.check_output(self.command_line,
                                             cwd=working_dir,
                                             stderr=stderr,
                                             universal_newlines=True)
            status = 0
        except subprocess.CalledProcessError as ex:
            output = ex.output
            status = ex.returncode
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
