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

def copy_dir_contents(src, dst, replace_spaces=True):
    """
    Copy the contents of one directory into another

    Arguments:
      src (str): path of source directory
      dst (str): path of destination directory
      replace_spaces (bool): if True (default) then
        replace spaces in source names with
        underscores in the destination names
    """
    for f in bcf_utils.walk(src):
        if f == src:
            # Don't try to copy the top-level dir
            continue
        # Make destination name
        ff = os.path.relpath(f, src)
        if replace_spaces:
            ff = ff.replace(" ","_")
        ff = os.path.join(dst, ff)
        if os.path.exists(ff):
            # Skip existing files
            logger.warning(f"{ff}: file with this name already exists, "
                           f"skipping")
        else:
            print(f"{os.path.relpath(ff, dst)}")
            if os.path.isdir(f):
                os.mkdir(ff)
            else:
                # Copy the file directly
                shutil.copyfile(f, ff, follow_symlinks=False)

# Main function

def main():
    """
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
    # FIXME add support for flattening directory structure
    ##p.add_argument('--flatten', action='store_true',
    ##               help="copy files without replicating the source "
    ##               "directory structure")
    # FIXME add support for overwriting files at the destination
    ##p.add_argument('--overwrite', action='store_true',
    ##               help="overwrite existing files (default is to skip "
    ##               "existing files)")
    p.add_argument('--runner',action='store',
                   help="specify the job runner to use for executing "
                   "'rsync' operations (defaults to job runner defined "
                   "for copying in config file [%s])" % default_runner)
    p.add_argument("src", metavar="SOURCE",
                   help="source data (file or directory) to copy; can "
                   "be on a local or remote file system")
    p.add_argument("dst", metavar="DEST",
                   help="destination on local file system")
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
        # Directory copy
        # Rsync contents to a temporary dir
        # FIXME probably only want to do this for remote
        # FIXME sources in future
        print("Copying directory contents")
        with tempfile.TemporaryDirectory() as d:
            # Do rsync
            print(f"Temporary directory: {d}")
            if src.endswith("/"):
                rsync_src = src
            else:
                rsync_src = src + "/"
            rsync = applications.general.rsync(rsync_src, d,
                                               escape_spaces=False)
            print(f"Running {rsync.command_line}")
            # Run rsync command via job runner
            rsync_job = SchedulerJob(runner,
                                     rsync.command_line,
                                     name="rsync_data",
                                     working_dir=d,
                                     log_dir=d)
            job_id = rsync_job.start()
            rsync_job.wait()
            with open(rsync_job.log, "rt") as fp:
                output = fp.read()
            if rsync_job.err:
                with open(rsycn_job.err, "rt") as fp:
                    output += fp.read()
            if rsync_job.exit_code != 0:
                # Rsync didn't complete successfully
                logger.error(f"Error running rsync: {output}")
                return 1
            else:
                print(f"Rsync ok: {output}")
            os.remove(rsync_job.log)
            if rsync_job.err:
                os.remove(rsync_job.err)
            # Make final directory if it doesn't exist
            if not os.path.exists(dst):
                os.mkdir(dst)
                print(f"Made destination directory '{dst}'") 
            # Copy the contents to final location
            print(f"Copying into '{dst}'")
            copy_dir_contents(d, dst)
    else:
        # File copy
        print("Copying file")
        # Rsync file to a temporary dir
        with tempfile.TemporaryDirectory() as d:
            print(f"Temporary directory: {d}")
            # Do rsync
            rsync = applications.general.rsync(src, d,
                                               escape_spaces=False)
            print(f"Running {rsync.command_line}")
            retcode, output = rsync.subprocess_check_output()
            if retcode != 0:
                # Rsync didn't complete successfully
                logger.error(f"Error running rsync: {output}")
                return 1
            else:
                print(f"Rsync ok: {output}")
            # Copy file to final location
            f = os.path.basename(Location(src).path)
            if os.path.isdir(dst):
                print(f"{f}")
                shutil.copyfile(os.path.join(d, f),
                                os.path.join(dst, f))
            else:
                print(f"{os.path.basename(dst)}")
                shutil.copyfile(os.path.join(d, f), dst)

