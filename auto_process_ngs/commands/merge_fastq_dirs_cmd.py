#!/usr/bin/env python
#
#     merge_fastq_dirs_cmd.py: implement auto process merge_fastq_dirs command
#     Copyright (C) University of Manchester 2019 Peter Briggs
#
#########################################################################

#######################################################################
# Imports
#######################################################################

import os
import shutil
import logging
from ..applications import Command
from ..fileops import copytree_command
from ..fileops import copy_command
from ..simple_scheduler import SimpleScheduler
from bcftbx import IlluminaData
from bcftbx.utils import list_dirs
from bcftbx.utils import mkdir

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Command functions
#######################################################################

def merge_fastq_dirs(ap,primary_unaligned_dir,output_dir=None,
                     dry_run=False):
    """
    Combine multiple 'unaligned' output directories into one
    
    This method combines the output from multiple runs of
    CASAVA/bcl2fastq into a single 'unaligned'-equivalent
    directory.

    Currently it operates in an automatic mode and should
    detect additional 'unaligned' dirs on its own.

    Arguments:
      ap (AutoProcessor): autoprocessor pointing to the parent
        analysis directory
      primary_unaligned_dir (str): the 'unaligned' dir that
        data from from all others will be put into (relative
        path), unless overridden by 'output_dir' argument
      output_dir (str): optional, new 'unaligned' dir that
        will be created to hold merged data (relative path,
        defaults to 'primary_unaligned_dir')
      dry_run (boolean): if True then just report operations
        that would have been performed.
    """
    if primary_unaligned_dir is None:
        raise Exception("Primary unaligned dir not defined")
    # Output directory
    if output_dir is None:
        output_dir = primary_unaligned_dir
    print("Fastqs will be merged into '%s'" % output_dir)
    # Collect unaligned dirs
    print("Collecting bcl2fastq directories")
    primary_illumina_data = None
    unaligned_dirs = {}
    for dirn in list_dirs(ap.analysis_dir):
        try:
            illumina_data = IlluminaData.IlluminaData(ap.analysis_dir,
                                                      unaligned_dir=dirn)
            if dirn == primary_unaligned_dir:
                print("* %s (primary dir)" % dirn)
                primary_illumina_data = illumina_data
            elif dirn.endswith(".bak") or dirn.startswith("save."):
                print("Ignoring %s" % dirn)
            else:
                print("* %s" % dirn)
                unaligned_dirs[dirn] = illumina_data
        except Exception as ex:
            logger.debug("Rejecting %s: %s" % (dirn,ex))
    # Check primary unaligned dir
    if primary_illumina_data is None:
        raise Exception("Primary dir '%s' doesn't exist, or doesn't "
                        "contain data?" % primary_unaligned_dir)
    # Is there anything to do?
    if not unaligned_dirs:
        print("No extra bcl2fastq output directories found, nothing to do")
        return 0
    # Make log directory and set up scheduler (if not dry run)
    if not dry_run:
        ap.set_log_dir(ap.get_log_subdir('merge_fastq_dirs'))
        runner = ap.settings.general.default_runner
        runner.set_log_dir(ap.log_dir)
        sched = SimpleScheduler(
            runner=runner,
            max_concurrent=ap.settings.general.max_concurrent_jobs,
            poll_interval=ap.settings.general.poll_interval
        )
        sched.start()
        jobs = []
    # Top-level for undetermined reads
    if primary_illumina_data.undetermined.dirn != \
       primary_illumina_data.unaligned_dir:
        undetermined_dir = os.path.basename(
            primary_illumina_data.undetermined.dirn)
    else:
        undetermined_dir = None
    # Do sanity checks before proceeding
    print("Checking primary data directory")
    fmt = primary_illumina_data.format
    paired_end = primary_illumina_data.paired_end
    no_lane_splitting = (len(primary_illumina_data.lanes) == 1) \
                        and (primary_illumina_data.lanes[0] is None)
    print("* Format: %s" % fmt)
    print("* no-lane-splitting: %s" % ('yes' if no_lane_splitting
                                       else 'no'))
    print("* paired-end: %s" % ('yes' if paired_end else 'no'))
    print("* undetermined dir: %s" % undetermined_dir)
    consistent_data = True
    for unaligned_dir in unaligned_dirs:
        illumina_data = unaligned_dirs[unaligned_dir]
        fmt0 = illumina_data.format
        no_lane_splitting0 = (len(illumina_data.lanes) == 1) \
                             and (primary_illumina_data.lanes[0] is None)
        if (fmt0 != fmt) or (no_lane_splitting0 != no_lane_splitting):
            print("!!! %s: inconsistent format to primary data dir !!!" %
                  unaligned_dir)
            consistent_data = False
    if not consistent_data:
        raise Exception("Data directories not consistent with primary "
                        "dir '%s'" % primary_unaligned_dir)
    # Collect the projects from the extra directories
    projects = []
    undetermined = []
    for unaligned_dir in unaligned_dirs:
        print("Examining projects in %s:" % unaligned_dir)
        illumina_data = unaligned_dirs[unaligned_dir]
        for project in illumina_data.projects:
            if not filter(lambda p: p.name == project.name,
                          projects):
                print("- %s: will be merged in" % project.name)
                projects.append(project)
            else:
                raise Exception("collision: %s already exists" %
                                project.name)
        # Deal with undetermined reads
        if illumina_data.undetermined is not None:
            print("Examining undetermined samples:")
            if no_lane_splitting:
                # No lane info: should merge undetermined fastqs
                for sample in illumina_data.undetermined.samples:
                    print("- %s: reads will be concatenated" % sample.name)
                    undetermined.append(sample)
            else:
                for sample in illumina_data.undetermined.samples:
                    if not filter(lambda s: s.name == sample.name,
                                  undetermined):
                        print("- %s: will be merged in" % sample.name)
                        undetermined.append(sample)
                    else:
                        raise Exception("collision: %s already exists" %
                                        sample.name)
        else:
            print("No undetermined samples")
    # Collect any remaining projects from the primary
    # unaligned directory
    print("Examining projects in primary dir %s:" %
          primary_unaligned_dir)
    for project in primary_illumina_data.projects:
        if not filter(lambda p: p.name == project.name,
                      projects):
            print("- %s: will be merged in" % project.name)
            projects.append(project)
        else:
            print("- %s: already exists, will be discarded" %
                  project.name)
    # Sort out the undetermined reads
    print("Examining undetermined samples:")
    if no_lane_splitting:
        # No lane info: should merge undetermined fastqs
        for sample in primary_illumina_data.undetermined.samples:
            print("- %s: reads will be concatenated" % sample.name)
            undetermined.insert(0,sample)
    else:
        for sample in primary_illumina_data.undetermined.samples:
            if not filter(lambda s: s.name == sample.name,
                          undetermined):
                print("- %s: will be merged in" % sample.name)
                undetermined.insert(0,sample)
            else:
                print("- %s: already exists, will be discarded" %
                      sample.name)
    # Make a new directory for the merging
    merge_dir = os.path.join(ap.analysis_dir,
                             output_dir + ".new")
    if undetermined_dir is not None:
        merge_undetermined_dir = os.path.join(merge_dir,
                                              undetermined_dir)
    else:
        merge_undetermined_dir = merge_dir
    if not dry_run:
        print("Making temporary merge directory %s" % merge_dir)
        mkdir(merge_dir)
        if not os.path.exists(merge_undetermined_dir):
            print("Making directory for undetermined %s" %
                  merge_undetermined_dir)
            mkdir(merge_undetermined_dir)
    # Copy the projects
    print("Importing projects:")
    for project in projects:
        print("- %s" % project.name)
        project_dir = os.path.join(merge_dir,
                                   os.path.basename(project.dirn))
        cmd = copytree_command(project.dirn,project_dir)
        print("- Running %s" % cmd)
        if not dry_run:
            job = sched.submit(cmd,
                               name="copy_project.%s" % project.name,
                               wd=merge_dir)
            print("Job: %s" % job)
            jobs.append(job)
    # Handle the undetermined reads
    print("Dealing with undetermined reads:")
    if no_lane_splitting:
        # No lane info: merge undetermined fastqs
        if len(undetermined) == 1:
            # Only one undetermined sample - copy Fastqs
            for read in (1,2):
                if read == 2 and not paired_end:
                    break
                fastqs = sample.fastq_subset(read_number=read,
                                             full_path=True)
                for fq in fastqs:
                    cmd = copy_command(fq,merge_undetermined_dir)
                    print("- Running %s" % cmd)
                    if not dry_run:
                        job = sched.submit(cmd,
                                           name="copy_undetermined.R%s" % read,
                                           wd=merge_dir)
                        print("Job: %s" % job)
                        jobs.append(job)
        else:
            # Multiple undetermined samples - concat Fastqs
            for read in (1,2):
                if read == 2 and not paired_end:
                    break
                cmd = Command('concat_fastqs.py')
                for sample in undetermined:
                    fastqs = sample.fastq_subset(read_number=read,
                                                 full_path=True)
                    cmd.add_args(*fastqs)
                cmd.add_args(os.path.join(
                    merge_undetermined_dir,
                    "Undetermined_S0_R%s_001.fastq.gz" % read))
                print("- Running %s" % cmd)
                if not dry_run:
                    job = sched.submit(cmd,
                                       name="merge_undetermined.R%s" % read,
                                       wd=merge_dir)
                    print("Job: %s" % job)
                    jobs.append(job)
    else:
        for sample in undetermined:
            print("- %s" % sample.name)
            if fmt == "bcl2fastq2":
                # Hardlink copy fastqs directly
                sample_dir = merge_undetermined_dir
                if not dry_run:
                    for fq in sample.fastq:
                        src_fq = os.path.join(sample.dirn,fq)
                        dst_fq = os.path.join(sample_dir,fq)
                        os.link(src_fq,dst_fq)
            else:
                # Just copy directory tree wholesale
                sample_dir = os.path.join(merge_undetermined_dir,
                                          os.path.basename(sample.dirn))
                cmd = copytree_command(sample.dirn,sample_dir)
                print("- Running %s" % cmd)
                if not dry_run:
                    job = sched.submit(cmd,
                                       name="copy_sample_dir.%s" %
                                       sample.name,
                                       wd=merge_dir)
                    print("Job: %s" % job.name)
                    jobs.append(job)
    # Make expected subdirs for bcl2fastq2
    if not dry_run and fmt == "bcl2fastq2":
        for dirn in ('Reports','Stats'):
            mkdir(os.path.join(merge_dir,dirn))
            # Add a hidden placeholder to preserve these directories
            # on rsync -m (prune empty dirs)
            with open(os.path.join(merge_dir,dirn,'.placeholder'),'w') as fp:
                fp.write("")
    # Wait for scheduler jobs to complete
    if not dry_run:
        sched.wait()
        sched.stop()
        # Check job exit status
        exit_status = 0
        for j in jobs:
            exit_status += j.exit_status
            if j.exit_status != 0:
                logger.warning("Job failed: %s" % j)
        if exit_status:
            logger.critical("One or more jobs failed (non-zero "
                            "exit status)")
            return exit_status
    # Move all the 'old' directories out of the way
    all_unaligned = [u for u in unaligned_dirs]
    all_unaligned.append(primary_unaligned_dir)
    for unaligned_dir in all_unaligned:
        unaligned_backup = os.path.join(ap.analysis_dir,
                                        "save.%s" %
                                        unaligned_dir)
        print("Moving %s to %s" % (unaligned_dir,
                                   unaligned_backup))
        if not dry_run:
            shutil.move(os.path.join(ap.analysis_dir,unaligned_dir),
                        unaligned_backup)
    # Rename the merged directory
    print("Renaming %s to %s" % (merge_dir,output_dir))
    if not dry_run:
        shutil.move(merge_dir,
                    os.path.join(ap.analysis_dir,output_dir))
    # Reset the bcl2fastq dir
    if not dry_run:
        ap.params['unaligned_dir'] = output_dir
    # Make a new 'projects.info' metadata file
    project_metadata_file = os.path.join(ap.analysis_dir,
                                         'projects.info')
    if os.path.exists(project_metadata_file):
        print("Moving existing projects.info file out of the way")
        if not dry_run:
            os.rename(project_metadata_file,
                      os.path.join(ap.analysis_dir,
                                   'save.projects.info'))
    print("Creating new projects.info file")
    if not dry_run:
        ap.make_project_metadata_file()
    return 0
