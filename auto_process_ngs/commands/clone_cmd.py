#!/usr/bin/env python
#
#     clone_cmd.py: implement auto process clone command
#     Copyright (C) University of Manchester 2019 Peter Briggs
#
#########################################################################

#######################################################################
# Imports
#######################################################################

import os
import logging
import shutil
from ..analysis import AnalysisProject
from ..metadata import AnalysisDirParameters
import bcftbx.utils as bcf_utils

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Command functions
#######################################################################

def clone(ap,clone_dir,copy_fastqs=False):
    """
    Make a 'clone' (i.e. copy) of an analysis directory

    Makes a copy of the analysis directory, including metadata
    and parameters but ignoring any project subdirectories.

    By default the 'unaligned' directory in the new directory is
    simply a symlink from the original directory; set the
    'copy_fastqs' to make copies instead.

    Arguments
      ap (AutoProcessor): autoprocessor pointing to the parent
        analysis directory
      clone_dir (str): path to the new directory to create as a
        clone (must not already exist).
      copy_fastqs (boolean): set to True to copy the Fastq files
        (otherwise default behaviour is to make symlinks)
    """
    clone_dir = os.path.abspath(clone_dir)
    if os.path.exists(clone_dir):
        # Directory already exists
        logger.critical("Target directory '%s' already exists" %
                        clone_dir)
        raise Exception("Clone failed: target directory '%s' "
                        "already exists" % clone_dir)
    bcf_utils.mkdir(clone_dir)
    # Copy metadata and parameters
    for f in (ap.metadata_file,ap.parameter_file):
        if os.path.exists(f):
            shutil.copy(f,os.path.join(clone_dir,os.path.basename(f)))
    # Link to or copy fastqs
    unaligned_dir = os.path.join(ap.analysis_dir,ap.params.unaligned_dir)
    if os.path.isdir(unaligned_dir):
        clone_unaligned_dir = os.path.join(clone_dir,
                                           os.path.basename(unaligned_dir))
        if not copy_fastqs:
            # Link to unaligned dir
            print "Symlinking %s" % clone_unaligned_dir
            os.symlink(unaligned_dir,clone_unaligned_dir)
        else:
            # Copy unaligned dir
            print "Copying %s" % clone_unaligned_dir
            shutil.copytree(unaligned_dir,clone_unaligned_dir)
    else:
        print "No 'unaligned' dir found"
    # Duplicate project directories
    projects = ap.get_analysis_projects()
    if projects:
        print "Duplicating project directories:"
        for project in ap.get_analysis_projects():
            print "-- %s" % project.name
            fastqs = project.fastqs
            new_project = AnalysisProject(
                project.name,
                os.path.join(clone_dir,project.name),
                user=project.info.user,
                PI=project.info.PI,
                library_type=project.info.library_type,
                single_cell_platform=project.info.single_cell_platform,
                organism=project.info.organism,
                run=project.info.run,
                comments=project.info.comments,
                platform=project.info.platform)
            new_project.create_directory(fastqs=fastqs,
                                         link_to_fastqs=(not copy_fastqs))
    # Copy additional files, if found
    for f in ("SampleSheet.orig.csv",
              ap.params.sample_sheet,
              ap.params.stats_file,
              ap.params.project_metadata,):
        srcpath = os.path.join(ap.analysis_dir,f)
        if os.path.exists(srcpath):
            shutil.copy(srcpath,clone_dir)
    # Create the basic set of subdirectories
    for subdir in ('logs','ScriptCode',):
        bcf_utils.mkdir(os.path.join(clone_dir,subdir))
    # Update the settings
    parameter_file = os.path.join(clone_dir,
                                  os.path.basename(ap.parameter_file))
    params = AnalysisDirParameters(filen=os.path.join(
        clone_dir,
        os.path.basename(ap.parameter_file)))
    for s in ("sample_sheet",):
        params[s] = os.path.join(clone_dir,
                                 os.path.relpath(params[s],
                                                 ap.analysis_dir))
    params.save()
