#!/usr/bin/env python
#
#     galaxy.py: utility classes & functions for interacting with Galaxy
#     Copyright (C) University of Manchester 2016 Peter Briggs
#
"""
Set of utility classes and functions for interacting with a Galaxy
instance:

- LibraryFolder: class representing a Galaxy data library folder
- LibraryDataset: class representing a Galaxy dataset
- get_library_contents: fetches contents of a data library
- build_library_directory: creates a directory structure with
  uncompressed FASTQs
- create_data_library: builds and populates a Galaxy data library

"""

#######################################################################
# Imports
########################################################################

import os
import time
import logging
from bioblend import galaxy
from nebulizer.core import get_galaxy_instance
from nebulizer.core import turn_off_urllib3_warnings
from nebulizer.libraries import library_id_from_name
from nebulizer.libraries import folder_id_from_name
from nebulizer.libraries import create_folder
from nebulizer.libraries import add_library_datasets
from nebulizer.libraries import split_library_path
from bcftbx.FASTQFile import get_fastq_file_handle
from bcftbx.utils import mkdir
from .utils import split_user_host_dir

#######################################################################
# Classes
#######################################################################

# These are classes and functions that should be
# considered for incorporation into nebulizer.libraries

class LibraryFolder:
    # FIXME: rename to Folder
    # FIXME: move to nebulizer.libraries
    def __init__(self,folder_data):
        self.id = folder_data['id']
        self.name = folder_data['name']
        self.type = folder_data['type']
        self.url = folder_data['url']

class LibraryDataset:
    # FIXME: rename to Dataset
    # FIXME: move to nebulizer.libraries
    def __init__(self,dataset_data):
        self.id = dataset_data['id']
        self.name = dataset_data['name']
        self.type = dataset_data['type']
        self.url = dataset_data['url']

#######################################################################
# Functions
#######################################################################

def get_library_contents(gi,path):
    """
    Return the contents of a data library

    Returns a list of LibraryFolder and LibraryDataset
    instances representing the contents of the specified
    Galaxy data library.

    FIXME: should filter on full path (currently
    lists everything)

    Arguments:
      gi (GalaxyInstance): bioblend GalaxyInstance
      path (str): path of the data library to fetch
        the contents of

    Returns:
      list: list of folders and datasets.

    """
    logging.debug("Path '%s'" % path)
    lib_client = galaxy.libraries.LibraryClient(gi)
    library_name,folder_path = split_library_path(path)
    logging.debug("library_name '%s'" % library_name)
    library_id = library_id_from_name(gi,library_name)
    if library_id is None:
        print "No library '%s'" % library_name
        return
    # Get library contents
    contents = []
    for item in lib_client.show_library(library_id,contents=True):
        if item['type'] == 'folder':
            contents.append(LibraryFolder(item))
        else:
            contents.append(LibraryDataset(item))
    return contents

def build_library_directory(analysis_dir,dest,projects=None):
    """
    Build and populate data library directory on server

    Arguments:
      analysis_dir (AnalysisDir): analysis directory to export
        files from
      dest (str): location of top-level data library directory
      projects (list): list of projects to export (default is to
        export all projects)
    
    """
    # Create and populate internal directory structure on server
    user,server,dirn = split_user_host_dir(dest)
    remote = (server is not None)
    path = os.path.join(dirn,analysis_dir.run_name)
    if remote:
        logging.critical("Dealing with remote systems not implemented")
        raise NotImplementedError("Cannot build library directory on remote system")
    run_path = os.path.join(dirn,analysis_dir.run_name)
    print "Creating %s" % run_path
    mkdir(run_path)
    for project in analysis_dir.get_projects(
            include_undetermined=False):
        if projects is not None and project.name not in projects:
            print "Ignoring project '%s'" % project.name
            continue
        project_path = os.path.join(run_path,project.name)
        print "Creating %s" % project_path
        mkdir(project_path)
        print "Populating with uncompressed Fastqs:"
        for sample in project.samples:
            for fq in sample.fastq:
                fqcp = os.path.join(project_path,
                                    os.path.basename(fq))
                if fqcp.endswith('.gz'):
                    fqcp = fqcp[0:-3]
                if os.path.exists(fqcp):
                    print "-- found: %s" % fqcp
                    continue
                print "-- %s" % fqcp
                with get_fastq_file_handle(fq) as fp:
                    with open(fqcp,'wb') as fpcp:
                        while True:
                            data = fp.read(102400)
                            if not data:
                                break
                            fpcp.write(data)

def create_data_library(galaxy_url,library_name,analysis_dir,dest,
                        projects=None,no_verify=False,
                        wait_interval=5.0):
    """
    Create and populate data library in Galaxy

    Arguments:
      galaxy_url (str): URL or alias of Galaxy server
      library_name (str): name of the library on the server
        to export the files to
      analysis_dir (AnalysisDir): analysis directory to
        export the files from
      dest (str): location of top-level data library directory
      projects (list): list of projects to export (default is to
        export all projects)
      no_verify (boolean): True to disable SSL certificate
        checking when connecting to Galaxy server (default is
        to verify certificate)
      wait_interval (float): number of seconds to wait for
        upload of file to complete

    """
    # Split up destination path
    user,server,dirn = split_user_host_dir(dest)
    remote = (server is not None)
    # Turn off SSL certificate verification?
    if no_verify:
        logging.warning("SSL certificate verification disabled")
        turn_off_urllib3_warnings()
    # Create data library structure in Galaxy
    print "Fetching Galaxy instance for %s" % galaxy_url
    gi = get_galaxy_instance(galaxy_url,verify_ssl=(not no_verify))
    if gi is None:
        logging.critical("%s: failed to connect to Galaxy instance" %
                         galaxy_url)
        raise GalaxyUploadException("%s: failed to connect to Galaxy "
                                    "instance" % galaxy_url)
    # Create the data library
    print "Creating folder for run in Galaxy"
    run_path = '/'.join((library_name,analysis_dir.run_name))
    library = library_id_from_name(gi,library_name)
    if library is None:
        logging.critical("%s: library not found" % library_name)
        raise GalaxyUploadException("%s: library not found" %
                                    library_name)
    run_folder = folder_id_from_name(gi,library,
                                     analysis_dir.run_name)
    if run_folder is None:
        description = "Data for %s run #%s datestamped %s" \
                      % (analysis_dir.metadata.platform.upper(),
                         analysis_dir.metadata.run_number,
                         analysis_dir.date_stamp)
        run_folder = create_folder(gi,run_path,description)
        if run_folder is None:
            logging.critical("%s: failed to create folder" % run_path)
            raise GalaxyUploadException("%s: failed to create folder" %
                                        run_path)
        else:
            print "Created run folder: '%s' '%s'" % (run_path,description)
    else:
        logging.warning("%s: run folder already exists"
                        % analysis_dir.run_name)
    for project in analysis_dir.get_projects(
            include_undetermined=False):
        if projects is not None and project.name not in projects:
            print "Ignoring project '%s'" % project.name
            continue
        project_name = "Fastqs (%s: %s)" % (project.name,
                                            project.info.organism.replace('/',','))
        project_path = '/'.join((run_path,project_name))
        if project.info.organism is not None:
            description = "%s: %s" % (project.name,
                                      project.info.organism.replace('/',','))
        else:
            description = "%s" % project.name
        project_folder = folder_id_from_name(gi,library,
                                             os.path.join(analysis_dir.run_name,
                                                          project_name))
        if project_folder is not None:
            print "Found existing subfolder for %s" % project.name
        else:
            print "Creating subfolder for project '%s'" % project.name
            print "-- name       : %s" % project_path
            print "-- description: %s" % description
            project_folder = create_folder(gi,project_path,description)
            if project_folder is None:
                logging.critical("%s: failed to create folder" %
                                 project_path)
                raise GalaxyUploadException("%s: failed to create folder" %
                                            project_path)
        print "Populating project folder:"
        for sample in project.samples:
            fastqs = []
            for fq in sample.fastq:
                fqcp = os.path.join(dirn,
                                    analysis_dir.run_name,
                                    project.name,
                                    os.path.basename(fq))
                if fqcp.endswith('.gz'):
                    fqcp = fqcp[0:-3]
                matches = filter(lambda x: os.path.basename(x.name) == 
                                 os.path.basename(fqcp),
                                 get_library_contents(gi,project_path))
                if len(matches):
                    print "-- found: %s" % fqcp
                else:
                    print "-- %s" % fqcp
                    fastqs.append(fqcp)
            if fastqs:
                for fq in fastqs:
                    # Add one at a time and pause to try
                    # and prevent overloading the server
                    add_library_datasets(gi,
                                         project_path,
                                         [fq,],
                                         from_server=True,
                                         link_only=True,
                                         file_type='fastqsanger')
                    time.sleep(wait_interval)

class GalaxyUploadException(Exception):
    """Used to indicate a problem with uploading data to Galaxy
    """
