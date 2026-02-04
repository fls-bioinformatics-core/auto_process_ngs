#!/usr/bin/env python
#
#     setup_analysis_dirs_cmd.py: implement 'setup_analysis_dirs' command
#     Copyright (C) University of Manchester 2018-2026 Peter Briggs
#
#########################################################################

#######################################################################
# Imports
#######################################################################

import os
import ast
import fnmatch
import shutil
import json
import logging
import bcftbx.IlluminaData as IlluminaData
from .. import analysis
from .. import tenx
from ..applications import identify_application
from ..utils import normalise_organism_name

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Command functions
#######################################################################

def setup_analysis_dirs(ap,
                        name=None,
                        unaligned_dir=None,
                        project_metadata_file=None,
                        ignore_missing_metadata=False,
                        short_fastq_names=False,
                        link_to_fastqs=False,
                        projects=None,
                        undetermined_project=None,
                        custom_metadata_items=None):
    """
    Construct and populate project analysis directories

    Arguments:
      ap (AutoProcess): AutoProcess instance
      unaligned_dir (str): optional, name of 'unaligned'
        subdirectory (defaults to value stored in parameters)
      name (str): (optional) identifier to append to output
        project directories
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
      custom_metadata_items (list): optional, list of strings
        defining additional custom metadata items to add to the
        core metadata items for each project (overrides custom
        items specified in configuration file)
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
        project_metadata_file=project_metadata_file)
    # Filter out the commented projects
    project_metadata = [p for p in project_metadata
                        if not p['Project'].startswith('#')]
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
        # Append the identifier, if supplied
        if name:
            new_project_name = project_name + '_' + name
        else:
            new_project_name = project_name
        print("Setting up analysis directory for project '%s'" %
              new_project_name)
        # Look up the application that matches the platform and library
        platform = line["SC_Platform"]
        if platform == ".":
            platform = None
        library_type = line["Library"]
        if library_type == ".":
            library_type = None
        print(f"-- Platform: '{platform}' Library: '{library_type}'")
        application = identify_application(platform, library_type)
        if application:
            print("-- identified application from platform and library")
            try:
                templates = application["setup"]["templates"]
            except KeyError:
                templates = []
            try:
                subdirs = application["setup"]["directories"]
            except KeyError:
                subdirs = []
        else:
            # Unable to identify application
            logger.error("Unknown application for platform/library combination for '%s': "
                         "'%s/%s'" % (line['Project'], platform, library_type))
            raise Exception("Unable to identify matching application")
        # Legacy failure for unidentified 10x applications
        if platform and application["platforms"] == ["*"]:
            logger.error("Unidentified platform for '%s': '%s'" % (line['Project'],
                                                                   platform))
            raise Exception("Unable to identify matching application for specific platform")
        # Custom project metadata items
        if custom_metadata_items is None:
            custom_metadata_items = ap.custom_project_metadata
        if custom_metadata_items:
            print(f"-- including extra metadata items for projects: {', '.join(custom_metadata_items)}")
        # Create the project
        project = analysis.AnalysisProject(
            new_project_name,
            os.path.join(ap.analysis_dir,
                         new_project_name),
            user=user,
            PI=PI,
            organism=organism,
            library_type=library_type,
            single_cell_platform=single_cell_platform,
            run=run_name,
            comments=comments,
            platform=ap.metadata.platform,
            sequencer_model=ap.metadata.sequencer_model,
            custom_metadata_items=custom_metadata_items)
        if project.exists:
            logging.warning("Project '%s' already exists, skipping" %
                            project.name)
            continue
        print("-- creating project: '%s'" % new_project_name)
        try:
            project.create_directory(
                illumina_data.get_project(project_name),
                short_fastq_names=short_fastq_names,
                link_to_fastqs=link_to_fastqs)
            n_projects += 1
        except IlluminaData.IlluminaDataError as ex:
            logger.warning("Failed to create project '%s': %s" %
                           (new_project_name,ex))
            continue
        # Create template files
        for template in templates:
            # Parse template definition
            if "(" in template and template.endswith(")"):
                try:
                    template, params = template.rstrip(")").split("(")
                    params = [p.strip().split("=") for p in params.split(",")]
                except Exception:
                    raise Exception("Could not parse template '%s'" % template)
            else:
                params = []
            # Set parameters
            include_probeset = False
            template_library_type = library_type
            for param in params:
                if param[0] == "include_probeset":
                    include_probeset = bool(param[1].lower() in ("true", "yes"))
                elif param[0] == "library_type":
                    template_library_type = param[1]
            # Set up the templates
            if template == "10x_multiome_libraries":
                # Config file for 10x single cell multiome
                # analyses
                f = "10x_multiome_libraries.info.template"
                print("-- making %s" % f)
                f = os.path.join(project.dirn,f)
                try:
                    with open(f,'wt') as fp:
                        fp.write(
                            "## 10x_multiome_libraries.info\n"
                            "## Link samples with complementary samples\n"
                            "## (can be in other runs and/or projects)\n"
                            "## See https://auto-process-ngs.readthedocs.io/en/latest/using/"
                            "setup_analysis_dirs.html#xgenomics-single-cell-multiome-linked-samples\n")
                        for sample in project.samples:
                            fp.write("#%s\t[RUN:][PROJECT][/SAMPLE]\n" % sample)
                except Exception as ex:
                    logger.warning("Failed to create '%s': %s" % (f,ex))
            elif template == "10x_multi_config":
                # Config file for running 10x cellranger multi
                print("-- setting up cellranger 'multi' config file for '%s/%s'" %
                      (single_cell_platform, library_type))
                # Identify organism
                organism_id = normalise_organism_name(organism)
                print("-- ...organism '%s' (ID '%s')" % (organism, organism_id,))
                # Target cellranger version
                try:
                    processing_software = ast.literal_eval(
                        ap.metadata.processing_software)
                    cellranger_version = processing_software['cellranger'][2]
                except Exception as ex:
                    print("-- unable to determine Cellranger version used "
                        "for processing: %s" % ex)
                    cellranger_version = None
                if not cellranger_version:
                    cellranger_version = tenx.DEFAULT_CELLRANGER_VERSION
                print("-- ...target Cellranger version '%s'" % cellranger_version)
                # Reference transcriptome
                try:
                    reference_dataset = ap.settings.organisms[organism_id]. \
                        cellranger_reference
                    if reference_dataset:
                        print("-- ...found %s" % reference_dataset)
                    else:
                        print("-- ...no matching transcriptome data")
                except Exception as ex:
                    logger.warning("Failed to locate 10xGenomics reference "
                                   "transcriptome for project '%s': %s" %
                                   (project_name, ex))
                    reference_dataset = None
                # Probe set
                if include_probeset:
                    try:
                        probe_set = ap.settings.organisms[organism_id].\
                                    cellranger_probe_set
                        if probe_set:
                            print("-- ...found %s" % probe_set)
                        else:
                            print("-- ...no matching probe set")
                    except Exception as ex:
                        logger.warning("Failed to locate 10xGenomics reference "
                                       "probe set for project '%s': %s" %
                                       (project_name,ex))
                        probe_set = None
                else:
                    probe_set = None
                # Make the template file
                f = "10x_multi_config.csv.template"
                print("-- making %s" % f)
                f = os.path.join(project.dirn,f)
                try:
                    tenx.utils.make_multi_config_template(
                        f,
                        reference=reference_dataset,
                        probe_set=probe_set,
                        fastq_dir=project.fastq_dir,
                        samples=[s.name for s in project.samples],
                        library_type=template_library_type,
                        cellranger_version=cellranger_version)
                except Exception as ex:
                    logger.warning("Error when attempting to create '%s': "
                                   "%s" % (f,ex))
                    if os.path.exists(f):
                        logger.warning("Template file was created but may be "
                                       "incomplete")
                    else:
                        logger.warning("Template file was not created")

        # Create subdirectories
        for subdir in subdirs:
            print(f"-- making '{subdir}' directory")
            d = os.path.join(project.dirn, subdir)
            os.mkdir(d)
    # Report how many projects were made
    print("Created %d project%s" % (n_projects,
                                    's' if n_projects != 1 else ''))
    # Also set up analysis directory for undetermined reads
    if undetermined_project is None:
        undetermined_project = 'undetermined'
        if name:
            undetermined_project = undetermined_project + '_' + name
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
            platform=ap.metadata.platform,
            custom_metadata_items=custom_metadata_items)
        if not undetermined.exists:
            print("Creating directory '%s' for analysing reads "
                  "with undetermined indices" % undetermined.name)
            undetermined.create_directory(illumina_data.undetermined,
                                          link_to_fastqs=link_to_fastqs)
        else:
            logger.warning("'%s' directory already exists, skipping" %
                           undetermined.name)
    return 0
