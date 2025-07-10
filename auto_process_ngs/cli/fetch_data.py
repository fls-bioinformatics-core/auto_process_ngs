#!/usr/bin/env python
#
#     cli/fetch_data.py: utility for fetching data files
#     Copyright (C) University of Manchester 2025 Peter Briggs
#
import os
import shutil
import argparse
import tempfile
from .. import applications
from .. import fileops
from .. import get_version
from ..settings import Settings
from ..simple_scheduler import SchedulerJob
from ..utils import Location
from bcftbx.JobRunner import fetch_runner
import bcftbx.utils as bcf_utils

# Logging
import logging
logging.basicConfig()
logger = logging.getLogger(__name__)


def run_command(cmd, runner, working_dir=None):
    """
    Run a command

    Given a Command instance, execute that comamnd using
    the specified job runner.

    Arguments:
      cmd (Command): the command line to be executed
      runner (JobRunner): the job runner to use
      working_dir (str): path to the working directory
        (or None to default to CWD)

    Returns:
      Integer: exit code from the job.
    """
    if working_dir is None:
        working_dir = os.getcwd()
    # Run command using job runner
    job = SchedulerJob(runner,
                       cmd.command_line,
                       name=cmd.command,
                       working_dir=working_dir,
                       log_dir=working_dir)
    job_id = job.start()
    job.wait()
    # Echo stdout and stderr and remove logs
    with open(job.log, "rt") as fp:
        output = fp.read()
        print(f"STDOUT:\n{output}")
        os.remove(job.log)
    if job.err:
        with open(job.err, "rt") as fp:
            output = fp.read()
            print(f"STDERR:\n{output}")
        os.remove(job.err)
    if job.exit_code != 0:
        # Job didn't complete successfully
        logger.error(f"Error running {cmd}: {job.exit_code}")
    return job.exit_code

def run_rsync(src, dst, runner, working_dir):
    """
    Run 'rsync' command to copy files and directories

    Arguments:
      src (str): path of source directory
      dst (str): path of destination directory
      runner (JobRunner): the job runner to use
      working_dir (str): path to the working directory
        (or None to default to CWD)
    """
    rsync = applications.general.rsync(src, dst, escape_spaces=False)
    print(f"Running {rsync.command_line}")
    return run_command(rsync, runner, working_dir=working_dir)

def copy_dir_contents(src, dst, replace_spaces=True, flatten=False,
                      overwrite=False):
    """
    Copy the contents of one directory into another

    Arguments:
      src (str): path of source directory
      dst (str): path of destination directory
      replace_spaces (bool): if True (default) then
        replace spaces in source names with
        underscores in the destination names
      flatten (bool): if True then don't replicate
        the source directory structure (default is
        to retain the directory structure)
      overwrite (bool): if True then replace (i.e.
        overwrite) existing destination files if an
        imported file has the same name (default is
        to skip existing files)
    """
    for f in bcf_utils.walk(src):
        if f == src:
            # Don't try to copy the top-level dir
            continue
        # Make destination name
        if flatten:
            if os.path.isdir(f):
                # Ignore directories when flattening
                continue
            ff = os.path.basename(f)
        else:
            ff = os.path.relpath(f, src)
        if replace_spaces:
            ff = ff.replace(" ","_")
        ff = os.path.join(dst, ff)
        if os.path.exists(ff):
            if os.path.isdir(ff):
                # Always skip existing directories
                continue
            if not overwrite:
                # Skip existing files
                logger.warning(f"{ff}: file with this name already "
                               f"exists, skipping")
                continue
            else:
                # Overwrite existing file
                logger.warning(f"{ff}: overwriting existing file")
        print(f"{os.path.relpath(ff, dst)}")
        if os.path.isdir(f):
            os.mkdir(ff)
        else:
            # Copy the file directly
            shutil.copyfile(f, ff, follow_symlinks=False)

def copy_file(src, dst, replace_spaces=True, overwrite=False):
    """
    Copy a file

    Arguments:
      src (str): path of source file
      dst (str): path of destination (can be an
        existing directory or a new file)
      replace_spaces (bool): if True (default) then
        replace spaces in source name with
        underscores in the destination name
      overwrite (bool): if True then replace (i.e.
        overwrite) existing destination file if the
        imported file has the same name (default is
        to skip existing file)
    """
    if os.path.isdir(dst):
        # Copy into directory preserving file name
        f = os.path.basename(src)
        if replace_spaces:
            f = f.replace(" ","_")
        dst = os.path.join(dst, f)
        if os.path.exists(dst):
            if not overwrite:
                # Skip existing file
                logger.warning(f"{dst}: file with this name already "
                               f"exists, skipping")
                return
            else:
                # Overwrite existing file
                logger.warning(f"{dst}: overwriting existing file")
    # Copy file
    print(f"{os.path.basename(dst)}")
    shutil.copyfile(src, dst)

# Main function

def main():
    """
    Implements the 'fetch_data.py' CLI utility
    """
    # Load configuration
    settings = Settings()

    # Collect defaults
    default_runner = settings.runners.rsync

    # Command line
    p = argparse.ArgumentParser(
        description="Copy files and directories from arbitrary "
        "locations to the local system")
    p.add_argument('--version', action='version',
                   version=("%%(prog)s %s" % get_version()))
    p.add_argument('--flatten', action='store_true',
                   help="copy files without replicating the source "
                   "directory structure")
    p.add_argument('--overwrite', action='store_true',
                   help="overwrite existing files (default is to skip "
                   "existing files)")
    p.add_argument('--runner',action='store',
                   help="specify the job runner to use for executing "
                   "'rsync' operations (defaults to job runner defined "
                   "for copying in config file [%s])" % default_runner)
    p.add_argument("src", metavar="SOURCE",
                   help="source data (file or directory) to copy; can "
                   "be on a local or remote file system. NB if source is "
                   "a directory then the contents are copied (not the "
                   "top-level directory)")
    p.add_argument("dst", metavar="DEST",
                   help="destination directory on local file system")
    args = p.parse_args()

    # Source
    src = args.src
    remote_src = Location(src).is_remote
    src_is_dir = fileops.isdir(src)
    print(f"Source     : {src}")
    if src_is_dir:
        print(f"           : (directory)")
    if remote_src:
        print(f"           : (remote)")

    # Destination
    dst = os.path.abspath(args.dst)
    print(f"Destination: {dst}")
    dst_exists = os.path.exists(dst)

    # Check there's something to transfer
    if not fileops.exists(src):
        if remote_src:
            logger.error(f"Source doesn't exist (or can't be reached)")
        else:
            logger.error(f"Source doesn't exist")
        return 1

    # Check source and destination compatibility
    if dst_exists:
        if src_is_dir and not os.path.isdir(dst):
            logger.error(f"Source is a directory but destination is "
                         f"not")
            return 1

    # Get runner for rsync and copy jobs
    if args.runner:
        runner = fetch_runner(args.runner)
    else:
        runner = default_runner

    # Fetch the data
    if src_is_dir:
        # Make final directory if it doesn't exist
        if not os.path.exists(dst):
            os.mkdir(dst)
            print(f"Made destination directory '{dst}'")
        # Directory copy
        print("Copying directory contents")
        if remote_src:
            # Rsync remote source to local temporary dir
            with tempfile.TemporaryDirectory() as d:
                # Do rsync
                print(f"Using temporary directory: {d}")
                if src.endswith("/"):
                    rsync_src = src
                else:
                    rsync_src = src + "/"
                status = run_rsync(rsync_src, d, runner, d)
                if status != 0:
                    # Rsync job failed
                    logger.error(f"rsyncing directory contents failed")
                    return 1
                # Copy the contents to final location
                print(f"Copying into '{dst}'")
                copy_dir_contents(d, dst, flatten=args.flatten,
                                  overwrite=args.overwrite)
        else:
            # Source is local directory
            # Do direct copy
            print(f"Copying into '{dst}'")
            copy_dir_contents(src, dst, flatten=args.flatten,
                              overwrite=args.overwrite)
    else:
        # File copy
        print("Copying file")
        if remote_src:
            # Rsync remote file to local temporary dir
            with tempfile.TemporaryDirectory() as d:
                print(f"Using temporary directory: {d}")
                # Do rsync
                status = run_rsync(src, d, runner, d)
                if status != 0:
                    # Rsync job failed
                    logger.error(f"rsyncing file failed")
                    return 1
                # Copy file to final location
                src = os.path.join(d, os.path.basename(Location(src).path))
                copy_file(src, dst, overwrite=args.overwrite)
        else:
            # Local copy
            copy_file(src, dst, overwrite=args.overwrite)
