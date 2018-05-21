#!/usr/bin/env python
#
#     setup_analysis_dirs_cmd.py: implement 'setup_analysis_dirs' command
#     Copyright (C) University of Manchester 2018 Peter Briggs
#
#########################################################################

#######################################################################
# Imports
#######################################################################

import os
import logging
import bcftbx.IlluminaData as IlluminaData
import auto_process_ngs.analysis as analysis

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Command functions
#######################################################################

def setup_analysis_dirs(ap,
                        unaligned_dir=None,
                        project_metadata_file=None,
                        ignore_missing_metadata=False,
                        short_fastq_names=False,
                        link_to_fastqs=False,
                        projects=None,
                        undetermined_project=None):
    """
    Construct and populate project analysis directories

    Arguments:
      unaligned_dir (str): optional, name of 'unaligned'
        subdirectory (defaults to value stored in parameters)
      project_metadata_file (str): optional, name of the
        'projects.info' metadata file to take project
        information from
      ignore_missing_metadata (bool): if True then make
        project directories for all projects even if metadata
        hasn't been set (default is to stop if metadata isn't
        set)
      short_fastq_names (bool): if True then use 'short' Fastq
        names (default is to use full Fastq names as output
        from bcl2fastq)
      link_to_fastqs (bool): if True then make symbolic links
        to the Fastq files in the source 'unaligned' subdir
        (default is to make hard links)
      projects (list): optional, subset of projects to create
        analysis dirs for (default is to attempt to create
        directories for all projects in the metadata file).
      undetermined_project (str): optional, specify name for
         project directory to create with 'undetermined' Fastqs
        (defaults to 'undetermined')
    """
    # Source location for fastq files
    if unaligned_dir is None:
        unaligned_dir = ap.params.unaligned_dir
    if unaligned_dir is None:
        logger.error("No unaligned directory, cannot build analysis "
                     "directories")
        raise Exception("Cannot build analysis directories")
    illumina_data = ap.load_illumina_data(unaligned_dir=unaligned_dir)
    # Project metadata file
    if project_metadata_file is None:
        project_metadata_file = ap.params.project_metadata
    if project_metadata_file is None:
        project_metadata_file = 'projects.info'
    if not os.path.exists(os.path.join(ap.params.analysis_dir,
                                       project_metadata_file)):
        logger.warning("No project metadata file '%s' found, "
                       "attempting to create" % project_metadata_file)
        ap.make_project_metadata_file(project_metadata_file)
        logger.warning("Update '%s' and rerun" % project_metadata_file)
        return 1
    project_metadata = ap.load_project_metadata(
        project_metadata_file=project_metadata_file,
        check=True)
    # Sanity check that the project data file has been populated
    got_project_data = True
    for line in project_metadata:
        for item in ('User','PI','Organism','Library',):
            if line[item] == '.':
                logger.warning("Missing data from %s for '%s': %s" %
                               (ap.params.project_metadata,
                                line['Project'],item))
                got_project_data = False
    if not got_project_data:
        if ignore_missing_metadata:
            logger.warning("Missing project metadata")
        else:
            logger.error("Missing project metadata")
            raise Exception("Missing project metadata")
    # Create the projects
    n_projects = 0
    for line in project_metadata:
        # Acquire the run name
        run_name = ap.run_name
        # Look up project data
        project_name = line['Project']
        user = line['User']
        PI = line['PI']
        organism = line['Organism']
        library_type = line['Library']
        single_cell_platform = line['SC_Platform']
        comments = line['Comments']
        # Check it's in the list
        if projects and project_name not in projects:
            logger.warning("Skipping '%s'" % project_name)
            continue
        # Create the project
        project = analysis.AnalysisProject(
            project_name,
            os.path.join(ap.analysis_dir,
                         project_name),
            user=user,
            PI=PI,
            organism=organism,
            library_type=library_type,
            single_cell_platform=single_cell_platform,
            run=run_name,
            comments=comments,
            platform=ap.metadata.platform)
        if project.exists:
            logging.warning("Project '%s' already exists, skipping" %
                            project.name)
            continue
        print "Creating project: '%s'" % project_name
        try:
            project.create_directory(
                illumina_data.get_project(project_name),
                short_fastq_names=short_fastq_names,
                link_to_fastqs=link_to_fastqs)
            n_projects += 1
        except IlluminaData.IlluminaDataError as ex:
            logger.warning("Failed to create project '%s': %s" %
                           (project_name,ex))
    # Tell us how many were made
    print "Created %d project%s" % (n_projects,'s' if n_projects != 1 else '')
    # Also set up analysis directory for undetermined reads
    if undetermined_project is None:
        undetermined_project = 'undetermined'
    undetermined = illumina_data.undetermined
    if illumina_data.undetermined is not None:
        undetermined = analysis.AnalysisProject(
            undetermined_project,
            os.path.join(
                ap.analysis_dir,
                undetermined_project),
            run=run_name,
            comments="Analysis of reads "
            "with undetermined indices",
            platform=ap.metadata.platform)
        if not undetermined.exists:
            print "Creating directory '%s' for analysing reads " \
                "with undetermined indices" % undetermined.name
            undetermined.create_directory(illumina_data.undetermined,
                                          link_to_fastqs=link_to_fastqs)
        else:
            logger.warning("'%s' directory already exists, skipping" %
                           undetermined.name)
    return 0
