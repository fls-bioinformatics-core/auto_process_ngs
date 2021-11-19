#!/usr/bin/env python
#
#     pipeliner.py: utilities for managing conda environments
#     Copyright (C) University of Manchester 2021 Peter Briggs
#
"""
Module providing utility classes and functions to help with managing
``conda`` environments:

- CondaWrapper: wrapper for ``conda``, including environment creation
- CondaWrapperException: base class for exceptions from CondaWrapper
- CondaCreateEnvError: exception for errors when creating environments
- make_conda_env_name: construct consistent names for conda environments

"""

######################################################################
# Imports
######################################################################

import os
import logging
from bcftbx.JobRunner import ResourceLock
from bcftbx.utils import find_program
from .command import Command

# Module specific logger
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

######################################################################
# Module constants
######################################################################

# Default channels for conda dependency resolution
DEFAULT_CONDA_CHANNELS = (
    'conda-forge',
    'bioconda',
    'defaults',
)

######################################################################
# Classes
######################################################################

class CondaWrapper(object):
    """
    Class for installing conda and creating environments

    Example usage:

    >>> conda = Conda('/usr/local/miniconda/bin/conda',
    ...               channels=('bioconda','conda-forge'),
    ...               env_dir='_conda_envs')
    >>> conda.install()
    >>> conda.create_env('multiqc@1.9','multiqc=1.9')
    """
    # Details of Miniconda installer
    _miniconda = {
        'url':
        "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh",
        'sha256':
        "1314b90489f154602fd794accfc90446111514a5a72fe1f71ab83e07de9504a7",
    }
    def __init__(self,conda=None,env_dir=None,channels=None):
        """
        Create a new CondaWrapper instance

        Arguments:
          conda (str): path to conda executable
          env_dir (str): optional, non-default directory
            for conda environments
          channels (list): optional, list of non-default
            channels to use for installing packages
        """
        # Conda executable
        if conda is None:
            conda = find_program("conda")
        self._conda = conda
        if self._conda:
            self._conda = os.path.abspath(self._conda)
            conda_dir = os.sep.join(self._conda.split(os.sep)[:-2])
        else:
            conda_dir = None
        self._conda_dir = conda_dir
        # Default location for environments
        if env_dir:
            env_dir = os.path.abspath(env_dir)
        elif self._conda_dir:
            env_dir = os.path.join(self._conda_dir,'envs')
        self._env_dir = env_dir
        # Channels
        if channels:
            channels = [c for c in channels]
        elif channels is None:
            channels = DEFAULT_CONDA_CHANNELS
        else:
            channels = list()
        self._channels = channels
        # Lock for blocking operations
        self._lock_manager = ResourceLock()

    @property
    def conda(self):
        """
        Path to conda executable
        """
        return self._conda

    @property
    def version(self):
        """
        Return the conda version
        """
        if not self.is_installed:
            return None
        version_cmd = Command(self.conda,'--version')
        output = version_cmd.subprocess_check_output()[1]
        try:
            return output.split()[1].strip()
        except Exception as ex:
            logger.warning("Unable to get conda version")

    @property
    def is_installed(self):
        """
        Check whether conda is installed
        """
        if self.conda:
            return os.path.exists(self.conda)
        return False

    @property
    def env_dir(self):
        """
        Path to the directory for conda environments
        """
        return self._env_dir

    @property
    def list_envs(self):
        """
        Return a list of environments in the 'envs' directory
        """
        if self._env_dir and os.path.exists(self._env_dir):
            return sorted(os.listdir(self._env_dir))
        else:
            return []

    def install(self):
        """
        Install conda

        Downloads and runs the miniconda installer if
        conda is not already installed
        """
        if self.is_installed:
            # Already installed
            return
        # Download miniconda
        tmpdir = tempfile.mkdtemp()
        miniconda_installer = os.path.join(
            tmpdir,
            os.path.basename(self._miniconda['url'])
        )
        try:
            url = urlopen(self._miniconda['url'])
            with open(miniconda_installer,'wb') as fp:
                fp.write(url.read())
        except URLError as ex:
            raise Exception("Failed to download miniconda installer from "
                            "'%s': %s" (self._miniconda['url']))
        # Run the installer in silent mode
        install_cmd = Command('bash',
                              miniconda_installer,
                              '-b',
                              '-p',self._conda_dir)
        install_cmd.run_subprocess()

    def create_env(self,name,*packages,**args):
        """
        Create a new conda environment

        Arguments:
          name (str): name of the new environment
          packages (list): package specifications for
            packages to install (e.g. 'fastqc=0.11.3')
        """
        # Get a lock on create operation
        lock = self._lock_manager.acquire("conda.create_env")
        if name in self.list_envs:
            # Environment already exists
            logger.warning("'%s': environment already exists" % name)
        else:
            # Create new environment
            if not self.is_installed:
                raise PipelineError("Can't create environment: conda not "
                                "installed")
            # Base command
            create_cmd = Command(self.conda,'create','-y',)
            # Channels
            if self._channels:
                create_cmd.add_args('--override-channels')
                for channel in self._channels:
                    create_cmd.add_args('-c',channel)
            # Create environment in specific location
            if self._env_dir:
                create_cmd.add_args('--prefix',
                                    os.path.join(self._env_dir,name))
            else:
                create_cmd.add_args('-n',name)
            # Packages to install
            create_cmd.add_args(*packages)
            # Run the command
            status,output = create_cmd.subprocess_check_output()
            if status != 0:
                # Release the lock
                self._lock_manager.release(lock)
                # Raise exception
                raise CondaCreateEnvError(
                    status=status,
                    env_name=name,
                    cmdline=str(create_cmd),
                    output=output)
        if self._env_dir:
            env_name = os.path.join(self._env_dir,name)
        else:
            env_name = name
        # Release the lock
        self._lock_manager.release(lock)
        # Return name/path to conda environment
        return env_name

    def activate_env_cmd(self,name):
        """
        Fetch command to activate a conda environment

        Arguments:
          name (str): name/path of the environment to
            activate

        Returns:
          Command: Command instance to activate the
            named environment.
        """
        if self._conda_dir:
            return Command(
                'source',
                os.path.join(self._conda_dir,'bin','activate'),
                name)
        else:
            return None

######################################################################
# Custom exceptions
######################################################################

class CondaWrapperError(Exception):
    """
    Base class for conda-specific exceptions
    """

class CondaCreateEnvError(CondaWrapperError):
    """
    Exception raised when CondaWrapper class fails
    to create a new environment

    Arguments:
      message (str): error message
      env_name (str): name of the environment
      status (int): status returned by the conda command
      cmdline (str): command line for the conda command
        that generated the error
      output (str): output from the conda command
    """
    def __init__(self,message=None,env_name=None,status=None,
                 cmdline=None,output=None):
        if message is None:
            message = "Unable to create environment"
        self.message = message
        self.env_name = env_name
        self.status = status
        self.cmdline = cmdline
        self.output = output
        CondaWrapperError.__init__(self,self.message)

######################################################################
# Utility functions
######################################################################

def make_conda_env_name(*pkgs):
    """
    Return a name for a conda env based on package list

    Constructs a name for a conda environment, based on
    the packages specified for inclusion in the
    environment.

    The name consists of components of the form
    "NAME[@VERSION]" for each package specification
    "NAME[=VERSION]" (e.g. 'star=2.4.2a' becomes
    'star@2.4.2a').

    Components are joined by '+' and are sorted by package
    name (regardless of the order they are supplied in).

    Arguments:
      pkgs (list): list of conda package specifiers

    Returns:
      String: environment name constructed from package
        list.
    """
    return '+'.join(sorted([str(p) for p in pkgs])).replace('=','@')
