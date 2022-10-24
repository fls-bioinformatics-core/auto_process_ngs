#!/usr/bin/env python
#
#     indexes.py: utilities for building sequencing indexes
#     Copyright (C) University of Manchester 2022 Peter Briggs
#

"""
Provides a utility class ``IndexBuilder`` with methods for building
indexes for various aligners:

- bowtie
- bowtie2
- STAR
"""

#######################################################################
# Imports
#######################################################################

import os
import shutil
import logging
import tempfile
import atexit
from .command import Command
from .conda import CondaWrapper
from .conda import CondaWrapperError
from .conda import make_conda_env_name
from .settings import Settings
from .simple_scheduler import SchedulerJob

# Module-specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Classes
#######################################################################

class IndexBuilder:
    """
    Utility class for building aligner indexes

    Provides the following methods:

    - ``bowtie``
    - ``bowtie2``
    - ``STAR``

    Arguments:
      runner (JobRunner): JobRunner instance that will be
        used to run the index command
    """
    def __init__(self,runner):
        """
        Create a new IndexBuilder instance
        """
        self._runner = runner

    def _run(self,build_cmd,working_dir,conda_pkgs=None):
        """
        Internal: run an index building command

        Arguments:
          build_cmd (Command): command to build the index
          working_dir (str): path to a directory to use
            as the workspace
          conda_pkgs (list): conda packages required to
            build the index (optional)
        """
        # Conda packages
        conda_activate_cmd = None
        if conda_pkgs:
            # Report dependencies
            print("Require conda packages:")
            for pkg in conda_pkgs:
                print("- %s" % pkg)
            # Check if conda environments are enabled
            if Settings().conda.enable_conda:
                # Get location for conda environments
                conda_env_dir = Settings().conda.env_dir
                # Set up conda wrapper
                conda = CondaWrapper(env_dir=conda_env_dir)
                env_name = make_conda_env_name(*conda_pkgs)
                try:
                    conda.create_env(env_name,*conda_pkgs)
                    conda_env = os.path.join(conda_env_dir,env_name)
                    # Script fragment to activate the environment
                    conda_activate_cmd = conda.activate_env_cmd(conda_env)
                except CondaWrapperError as ex:
                    # Failed to acquire the environment
                    logger.warning("failed to acquire conda environment "
                                   "'%s': %s" % (env_name,ex))
            else:
                print("Conda dependency resolution not enabled")
        # Wrap the command in a script
        script_file = os.path.join(working_dir,
                                   "build_indexes.%s.sh" %
                                   build_cmd.command)
        build_cmd.make_wrapper_script(filen=script_file,
                                      prologue=conda_activate_cmd,
                                      quote_spaces=True)
        # Run the wrapper script
        build_indexes = SchedulerJob(
            self._runner,
            Command('/bin/bash','-l',script_file).command_line,
            name="build_indexes.%s" % build_cmd.command,
            working_dir=working_dir,
            log_dir=working_dir)
        build_indexes.start()
        try:
            build_indexes.wait()
        except KeyboardInterrupt as ex:
            logger.warning("Keyboard interrupt, terminating QC reporting")
            build_indexes.terminate()
            raise ex
        # Report errors
        if build_indexes.exit_code != 0:
            logger.critical("Index generation failed")
            print("STDOUT:")
            with open(build_indexes.log,'rt') as fp:
                print(fp.read())
            try:
                # Older versions of Job object
                # don't have 'err' property?
                print("STDERR:")
                with open(build_indexes.err,'rt') as fp:
                    print(fp.read())
            except AttributeError:
                logger.warning("Stderr not available")
        # Return the exit code
        return build_indexes.exit_code

    def get_working_dir(self,remove_on_exit=True):
        """
        Create and return a temporary working directory

        If 'remove_on_exit' is True then the working
        directory will be deleted on exit.
        """
        working_dir = tempfile.mkdtemp(prefix="__build_index.",
                                       suffix=".tmp",
                                       dir=os.getcwd())
        if remove_on_exit:
            atexit.register(self.remove_working_dir,working_dir)
        return working_dir

    def remove_working_dir(self,working_dir):
        """
        Remove working directory
        """
        if os.path.isdir(working_dir):
            print("Removing %s" % working_dir)
            shutil.rmtree(working_dir)

    def bowtie(self,fasta,out_dir,ebwt_basename=None,
               bowtie_version=None):
        """
        Build index for bowtie

        Arguments:
          fasta (str): path to Fasta file
          out_dir (str): path to output directory
            for index files
          ebwt_basename (str): optional base name for
            the output index (.ebwt) files
          bowtie_version (str): specify the version
            of Bowtie to use (only if conda dependency
            resolution being used)
        """
        # Fasta file
        fasta = os.path.abspath(fasta)
        # ebwt base name
        if not ebwt_basename:
            # Use Fasta file name
            ebwt_basename = os.path.splitext(
                os.path.basename(fasta))[0]
        # Output directory
        out_dir = os.path.abspath(out_dir)
        # Conda packages
        bowtie_pkg = "bowtie"
        if bowtie_version:
            bowtie_pkg += ("=%s" % bowtie_version)
        # Working directory
        working_dir = self.get_working_dir()
        # Command to build index
        build_index_cmd = Command("bowtie-build")
        build_index_cmd.add_args("-f",fasta,
                                 ebwt_basename)
        # Run the command
        print("%s" % build_index_cmd)
        ret_code = self._run(build_index_cmd,
                             working_dir,
                             conda_pkgs=(bowtie_pkg,))
        if ret_code != 0:
            raise Exception("%s: returned non-zero exit code (%s)"
                            % (build_index_cmd,ret_code))
        # Copy .ebwt files to final location
        if not os.path.isdir(out_dir):
            os.makedirs(out_dir)
        for f in filter(lambda x: x.endswith('.ebwt'),
                        os.listdir(working_dir)):
            shutil.copy(os.path.join(working_dir,f),
                        os.path.join(out_dir,f))
        print("Index files in %s" % out_dir)

    def bowtie2(self,fasta,out_dir,bt2_basename=None,
                nthreads=None,bowtie2_version=None):
        """
        Build index for bowtie2

        Arguments:
          fasta (str): path to Fasta file
          out_dir (str): path to output directory
            for index files
          bt2_basename (str): optional base name for
            the output index (.bt2) files
          nthreads (int): specify number of threads to
            use when making index (defaults to the
            number defined in the job runner)
          bowtie2_version (str): specify the version
            of Bowtie2 to use (only if conda dependency
            resolution being used)
        """
        # Fasta file
        fasta = os.path.abspath(fasta)
        # bt2 base name
        if not bt2_basename:
            # Use Fasta file name
            bt2_basename = os.path.splitext(
                os.path.basename(fasta))[0]
        # Output directory
        out_dir = os.path.abspath(out_dir)
        # Number of threads
        if not nthreads:
            nthreads = self._runner.nslots
        # Conda packages
        bowtie2_pkg = "bowtie2"
        if bowtie2_version:
            bowtie2_pkg += ("=%s" % bowtie2_version)
        # Working directory
        working_dir = self.get_working_dir()
        # Command to build index
        build_index_cmd = Command("bowtie2-build")
        if nthreads and nthreads > 1:
            build_index_cmd.add_args("--threads",nthreads)
        build_index_cmd.add_args("-f",fasta,
                                 bt2_basename)
        # Run the command
        print("%s" % build_index_cmd)
        ret_code = self._run(build_index_cmd,
                             working_dir,
                             conda_pkgs=(bowtie2_pkg,))
        if ret_code != 0:
            raise Exception("%s: returned non-zero exit code (%s)"
                            % (build_index_cmd,ret_code))
        # Copy .bt2 files to final location
        if not os.path.isdir(out_dir):
            os.makedirs(out_dir)
        for f in filter(lambda x: x.endswith('.bt2'),
                        os.listdir(working_dir)):
            shutil.copy(os.path.join(working_dir,f),
                        os.path.join(out_dir,f))
        print("Index files in %s" % out_dir)

    def STAR(self,fasta,annotation,out_dir,memory_limit=None,
             overhang=100,nthreads=None,star_version=None):
        """
        Build index for STAR

        Arguments:
          fasta (str): path to Fasta file
          annotation (str): path to annotation file
          out_dir (str): path to output directory
            for index files
          memory_limit (int): specify memory limit
            (optional; default: no limit specified)
          overhang (int): specify overhang (default:
            100)
          nthreads (int): specify number of threads to
            use when making index (defaults to the
            number defined in the job runner)
          star_version (str): specify the version of
            STAR to use (only if conda dependency
            resolution being used)
        """
        # Fasta file
        fasta = os.path.abspath(fasta)
        # Annotation file
        annotation = os.path.abspath(annotation)
        # Output directory
        out_dir = os.path.abspath(out_dir)
        # Number of threads
        if not nthreads:
            nthreads = self._runner.nslots
        # Conda packages
        star_pkg = "star"
        if star_version:
            star_pkg += ("=%s" % star_version)
        # Working directory
        working_dir = self.get_working_dir()
        # STAR index directory
        star_dir = os.path.join(working_dir,"star_index")
        os.makedirs(star_dir)
        # Command to build index
        build_index_cmd = Command("STAR",
                                  "--runMode","genomeGenerate")
        if nthreads and nthreads > 1:
            build_index_cmd.add_args("--runThreadN",nthreads)
        build_index_cmd.add_args("--genomeFastaFiles",fasta,
                                 "--sjdbGTFfile",annotation,
                                 "--sjdbOverhang",overhang,
                                 "--genomeDir",star_dir)
        if memory_limit:
            build_index_cmd.add_args("--limitGenomeGenerateRAM",
                                     memory_limit)
        # Run the command
        print("%s" % build_index_cmd)
        ret_code = self._run(build_index_cmd,
                             working_dir,
                             conda_pkgs=(star_pkg,))
        if ret_code != 0:
            raise Exception("%s: returned non-zero exit code (%s)"
                            % (build_index_cmd,ret_code))
        # Copy STAR index files to final location
        if not os.path.isdir(out_dir):
            os.makedirs(out_dir)
        for f in os.listdir(star_dir):
            shutil.copy(os.path.join(star_dir,f),
                        os.path.join(out_dir,f))
        print("Index files in %s" % out_dir)
