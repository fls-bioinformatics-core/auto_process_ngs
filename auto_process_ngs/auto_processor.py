#!/usr/bin/env python
#
#     auto_processor.py: automated processing of Illumina sequence data
#     Copyright (C) University of Manchester 2013-2019 Peter Briggs
#
#########################################################################
#
# auto_processor.py
#
#########################################################################

#######################################################################
# Imports
#######################################################################

import sys
import os
import subprocess
import logging
import shutil
import uuid
import time
import ast
import gzip
import urllib2
import bcftbx.IlluminaData as IlluminaData
import bcftbx.TabFile as TabFile
import bcftbx.utils as bcf_utils
import bcftbx.htmlpagewriter as htmlpagewriter
from bcftbx.JobRunner import fetch_runner
from bcftbx.FASTQFile import FastqIterator
from . import commands
from .analysis import AnalysisProject
from .analysis import run_reference_id
from .metadata import AnalysisDirParameters
from .metadata import AnalysisDirMetadata
from .metadata import ProjectMetadataFile
from .utils import edit_file
from .utils import get_numbered_subdir
from .bcl2fastq_utils import get_sequencer_platform
from .samplesheet_utils import check_and_warn
from .settings import Settings
from .exceptions import MissingParameterFileException
from . import get_version

#######################################################################
# Decorators
#######################################################################

def add_command(name,f):
    """
    Add a method to a class

    Implements an '@add_command' decorator which can be
    used to add a function to a class as a new method
    (aka 'command').

    For example:

    >>> def hello(cls):
    ...    print("Hello %s" % cls.person)
    ...
    >>> @command("greeting",hello)
    ... class Example:
    ...   def __init__(self):
    ...      self.person = "World"
    ...
    >>> Example().greeting()
    Running 'greeting' command
    Hello World
    'greeting': finished

    The function must accept a class instance as the
    first argument.
    """
    def wrapped_func(*args,**kws):
        # Wraps execution of the supplied
        # function to trap exceptions and
        # add additional commentary
        print("[%s] Running '%s' command" % (timestamp(),name))
        try:
            ret = f(*args,**kws)
        except Exception as ex:
            logging.fatal("%s: %s" % (name,ex))
            ret = 1
        else:
            print("[%s] %s: finished" % (timestamp(),name))
        return ret
    def timestamp():
        # Return the current time
        return time.strftime("%Y-%m-%d %H:%M:%S")
    def wrapper(cls):
        # Adds the supplied function to
        # to the class
        setattr(cls,name,wrapped_func)
        return cls
    return wrapper

#######################################################################
# Classes
#######################################################################

@add_command("setup",commands.setup)
@add_command("make_fastqs",commands.make_fastqs)
@add_command("analyse_barcodes",commands.analyse_barcodes)
@add_command("merge_fastq_dirs",commands.merge_fastq_dirs)
@add_command("setup_analysis_dirs",commands.setup_analysis_dirs)
@add_command("run_qc",commands.run_qc)
@add_command("publish_qc",commands.publish_qc)
@add_command("archive",commands.archive)
@add_command("report",commands.report)
@add_command("update_fastq_stats",commands.update_fastq_stats)
@add_command("import_project",commands.import_project)
@add_command("clone",commands.clone)
class AutoProcess(object):
    """
    Class implementing an automatic fastq generation and QC
    processing procedure for Illumina sequencing data

    """
    def __init__(self,analysis_dir=None,settings=None,
                 allow_save_params=True):
        """
        Create a new AutoProcess instance

        Arguments:
          analysis_dir (str): name/path for existing analysis
            directory
          settings (Settings): optional, if supplied then should
            be a Settings instance; otherwise use a default
            instance populated from the installation-specific
            'auto_process.ini' file
          allow_save_params (bool): if True then allow updates
            to parameters to be saved back to the parameter file
            (this is the default)

        """
        # Initialise
        self._master_log_dir = "logs"
        self._log_dir = self._master_log_dir
        # Load configuration settings
        if settings is None:
            settings = Settings()
        self.settings = settings
        # Create empty parameter and metadata set
        self.params = AnalysisDirParameters()
        self.metadata = AnalysisDirMetadata()
        # Set flags to indicate whether it's okay to save parameters
        self._save_params = False
        self._save_metadata = False
        # Set where the analysis directory actually is
        self.analysis_dir = analysis_dir
        if self.analysis_dir is not None:
            # Load parameters
            self.analysis_dir = os.path.abspath(self.analysis_dir)
            try:
                self.load_parameters(allow_save=allow_save_params)
            except MissingParameterFileException as ex:
                # No parameter file
                logging.warning("Failed to load parameters: %s (ignored)" % ex)
                logging.warning("Perhaps this is not an auto_process project?")
                # Attempt to detect existing data directory
                self.params['unaligned_dir'] = self.detect_unaligned_dir()
                if self.params.unaligned_dir is None:
                    logging.warning("Unable to find subdirectory containing data")
            except Exception as ex:
                logging.error("Failed to load parameters: %s" % ex)
                logging.error("Stopping")
                sys.exit(1)
            self.params['analysis_dir'] = self.analysis_dir
            # Load metadata
            try:
                self.load_metadata(allow_save=allow_save_params)
            except MissingParameterFileException as ex:
                # No metadata file
                logging.warning("Failed to load metadata: %s (ignored)" % ex)
                logging.warning("Consider running metadata --update?")
            except Exception as ex:
                # Some other problem
                logging.error("Failed to load metadata: %s" % ex)
                logging.error("Stopping")
                sys.exit(1)

    def add_directory(self,sub_dir):
        # Add a directory to the AutoProcess object
        dirn = os.path.join(self.analysis_dir,sub_dir)
        self.create_directory(dirn)
        return dirn

    def create_directory(self,dirn):
        # Make the specified directory, and any leading directories
        # that don't already exist
        if not os.path.exists(dirn):
            dir_path = os.sep
            for sub_dir in dirn.split(os.sep):
                dir_path = os.path.join(dir_path,sub_dir)
                if not os.path.exists(dir_path):
                    print("Making %s" % dir_path)
                    bcf_utils.mkdir(dir_path)

    def load_parameters(self,allow_save=True):
        """
        Load parameter values from file

        Arguments:
          allow_save (boolean): if True then allow params to be
            saved back to the parameter file (the default);
            otherwise don't allow save.

        """
        # check for parameter file
        if not self.has_parameter_file:
            raise MissingParameterFileException(
                "No parameter file %s" % self.parameter_file)
        # Read contents of parameter file and assign values
        logging.debug("Loading parameters from %s" %
                      self.parameter_file)
        self.params.load(self.parameter_file,strict=False)
        # File exists and can be read so set flag accordingly
        self._save_params = allow_save

    def save_parameters(self,alt_parameter_file=None,force=False):
        """
        Save parameters to file

        Arguments:
          alt_parameter_file (str): optional, path to an
            'alternative' parameter file; otherwise
            parameters are saved to the default file for the
            processing directory.
          force (boolean): if True then force the parameters
            to be saved even if saving was previously
            turned off (default is False i.e. don't force
            save).

        """
        if self._save_params or force:
            if alt_parameter_file is None:
                self.params.save(self.parameter_file)
            else:
                self.params.save(alt_parameter_file)

    def load_metadata(self,allow_save=True):
        """
        Load metadata values from file

        Arguments:
          allow_save (boolean): if True then allow metadata
            items to be saved back to the metadata file (the
            default); otherwise don't allow save.

        """
        # check for metadata file
        if not os.path.exists(self.metadata_file):
            raise MissingParameterFileException(
                "No metadata file %s" % self.metadata_file)
        # Read contents of metadata file and assign values
        logging.debug("Loading metadata from %s" % self.metadata_file)
        self.metadata.load(self.metadata_file)
        # File exists and can be read so set flag accordingly
        self._save_metadata = allow_save

    def save_metadata(self,alt_metadata_file=None,force=False):
        """
        Save metadata to file

        Arguments:
          alt_metadata_file (str): optional, path to an
            'alternative' metadata file; otherwise
            metadata are saved to the default file for the
            processing directory.
          force (boolean): if True then force the metadata
            to be saved even if saving was previously
            turned off (default is False i.e. don't force
            save).

        """
        if self._save_metadata or force:
            if alt_metadata_file is None:
                self.metadata.save(self.metadata_file)
            else:
                self.metadata.save(alt_metadata_file)

    def update_metadata(self):
        """
        Migrates and updates metadata values

        """
        # Migrate missing values from parameter file
        if self.has_parameter_file:
            # Migrate relevant values across
            print("Migrating metadata values from parameter file")
            for param in ('platform','run_number','source','assay'):
                if param not in self.params:
                    continue
                if self.metadata[param] is None:
                    logging.debug("Importing metadata item '%s': set to "
                              "'%s'" % (param,self.params[param]))
                    print("Importing metadata item '%s'" % param)
                    self.metadata[param] = self.params[param]
        # Run name
        if self.metadata.run_name is None:
            print("Attempting to set missing 'run_name' metadata item")
            self.metadata['run_name'] = self.run_name
        # Instrument-related metadata
        if self.metadata.instrument_name is None or \
           self.metadata.instrument_datestamp is None or \
           self.metadata.instrument_run_number is None:
            print("Attempting to set missing instrument metadata items")
            # Extract from run name
            try:
                datestamp,instrument,run_number,\
                    flow_cell_prefix,flow_cell_id = \
                    IlluminaData.split_run_name_full(self.run_name)
                if self.metadata.instrument_name is None:
                    self.metadata['instrument_name'] = instrument
                if self.metadata.instrument_datestamp is None:
                    self.metadata['instrument_datestamp'] = datestamp
                if self.metadata.instrument_run_number is None:
                    self.metadata['instrument_run_number'] = run_number
                if self.metadata.instrument_flow_cell_id is None:
                    self.metadata['instrument_flow_cell_id'] = \
                        flow_cell_prefix + flow_cell_id
            except Exception as ex:
                logging.warning("Unable to extract missing instrument metadata "
                                "from run name")
        # Sequencing platform
        if self.metadata.platform is None:
            # Attempt to look up the instrument name
            platform = get_sequencer_platform(
                self.analysis_dir,
                instrument=self.metadata.instrument_name,
                settings=self.settings)
            if platform:
                print("Setting 'platform' metadata item to %s" %
                      platform)
                self.metadata['platform'] = platform

    def edit_samplesheet(self):
        """
        Bring up SampleSheet in an editor
        """
        # Fetch the sample sheet
        sample_sheet_file = self.params.sample_sheet
        if sample_sheet_file is None:
            logging.error("No sample sheet file to edit")
            return
        edit_file(sample_sheet_file)
        # Check updated sample sheet and issue warnings
        if check_and_warn(sample_sheet_file=sample_sheet_file):
            logging.error("Sample sheet may have problems, see warnings above")

    def init_readme(self):
        """
        Create a new README file
        """
        if self.readme_file is None:
            readme_file = os.path.join(self.analysis_dir,'README')
            print("Initialising %s" % readme_file)
            with open(readme_file,'w') as fp:
                title = "Processing notes for %s" % \
                        os.path.basename(self.analysis_dir)
                fp.write("%s\n%s\n" % (title,'='*len(title)))
        else:
            logging.warning("'%s' already exists" % self.readme_file)

    def edit_readme(self):
        """
        Bring up README in an editor
        """
        if self.readme_file is None:
            logging.error("No README file to edit")
            return
        edit_file(self.readme_file,
                  append="\n[%s]" % time.ctime())

    def load_illumina_data(self,unaligned_dir=None):
        # Load and return an IlluminaData object
        if unaligned_dir is None:
            unaligned_dir = self.params.unaligned_dir
        if unaligned_dir is None:
            logging.error("Unaligned directory not specified, cannot load data")
            return None
        return IlluminaData.IlluminaData(self.analysis_dir,
                                         unaligned_dir=unaligned_dir)

    def load_project_metadata(self,project_metadata_file='projects.info',
                              check=True,update=False):
        # Load data from 'projects.info' metadata file which lists
        # and describes projects
        # check: if True then check existing metadata for consistency with fastq files
        # update: if True then update inconsistent metadata (i.e. add missing projects
        #         and remove ones that are inconsistent); implies 'check=True'
        if project_metadata_file is not None:
            filen = os.path.join(self.params.analysis_dir,project_metadata_file)
        else:
            filen = None
        logging.debug("Project metadata file: %s" % filen)
        try:
            illumina_data = self.load_illumina_data()
        except IlluminaData.IlluminaDataError as ex:
            logging.warning("Failed to load data from bcl2fastq output "
                            "(ignored): %s" % ex)
            illumina_data = None
        projects_from_dirs = self.get_analysis_projects_from_dirs()
        if filen is not None and os.path.exists(filen):
            # Load existing file and check for consistency
            logging.debug("Loading project metadata from existing file")
            project_metadata = ProjectMetadataFile(filen)
        else:
            # First try to populate basic metadata from existing projects
            logging.debug("Metadata file not found, guessing basic data")
            project_metadata = ProjectMetadataFile()
            projects = projects_from_dirs
            if not projects:
                # Get information from fastq files
                logging.warning("No existing project directories detected")
                logging.debug("Use fastq data from 'unaligned' directory")
                if illumina_data is None:
                    # Can't even get fastq files
                    logging.warning("Failed to load fastq data from '%s'" %
                                    self.params.unaligned_dir)
                else:
                    projects = illumina_data.projects
            # Populate the metadata file list of projects
            logging.debug("Project\tSample\tFastq")
            for project in projects:
                project_name = project.name
                sample_names = []
                for sample in project.samples:
                    sample_name = sample.name
                    for fastq in sample.fastq:
                        logging.debug("%s\t%s\t%s" % (project_name,sample_name,fastq))
                    sample_names.append(sample_name)
                project_metadata.add_project(project_name,sample_names)
            # Turn off redundant checking/updating
            check = False
            update = False
        # Perform consistency check or update
        if check or update:
            # Check that each project listed actually exists
            bad_projects = []
            for line in project_metadata:
                pname = line['Project']
                test_project = AnalysisProject(
                    pname,os.path.join(self.analysis_dir,pname))
                if not test_project.is_analysis_dir:
                    # Project doesn't exist
                    logging.warning("Project '%s' listed in metadata file doesn't exist" \
                                    % pname)
                    bad_projects.append(line)
            # Remove bad projects
            if update:
                logging.debug("Removing non-existent projects")
                for bad_project in bad_projects:
                    del(bad_project)
            # Check that all actual projects are listed
            for project in projects_from_dirs:
                if len(project_metadata.lookup('Project',project.name)) == 0:
                    # Project not listed
                    if project.name != 'undetermined':
                        logging.warning("Project '%s' not listed in metadata file" %
                                        project.name)
                    if update:
                        # Add line for unlisted project
                        logging.debug("Adding basic data for project '%s'" % project.name)
                        sample_names = []
                        for sample in project.samples:
                            sample_name = sample.name
                            for fastq in sample.fastq:
                                logging.debug("%s\t%s\t%s" % (project.name,sample_name,fastq))
                            sample_names.append(sample_name)
                        project_metadata.add_project(project.name,sample_names)
        # Return the metadata object
        return project_metadata

    def update_project_metadata_file(self,unaligned_dir=None,
                                     project_metadata_file='projects.info'):
        """
        Update project metadata file from bcl2fastq outputs

        Updates the contents of the project metadata file
        (default: "projects.info") from a bcl-to-fastq output
        directory, by adding new entries for projects in the
        bcl-to-fastq outputs which don't currently appear.

        Arguments:
          unaligned_dir (str): path to the bcl-to-fastq
            output directory relative to the analysis dir.
            Defaults to the unaligned dir stored in the
            analysis directory parameter file.
          project_metatadata_file (str): optional, path to
            the project metadata file to update
        """
        if project_metadata_file is not None:
            self.params['project_metadata'] = project_metadata_file
        print("Project metadata file: %s" % self.params.project_metadata)
        filen = os.path.join(self.analysis_dir,
                             self.params.project_metadata)
        if unaligned_dir is not None:
            self.params['unaligned_dir'] = unaligned_dir
        print("Unaligned_dir: %s" % self.params.unaligned_dir)
        illumina_data = IlluminaData.IlluminaData(
            self.analysis_dir,
            unaligned_dir=self.params.unaligned_dir)
        if os.path.exists(filen):
            # Load data from existing file
            print("Loading project metadata from existing file: %s" % filen)
            project_metadata = ProjectMetadataFile(filen)
        else:
            # New (empty) metadata file
            print("Creating new project metadata file: %s" % filen)
            project_metadata = ProjectMetadataFile()
        # Populate/update
        for project in illumina_data.projects:
            project_name = project.name
            sample_names = [s.name for s in project.samples]
            if project_name not in project_metadata:
                project_metadata.add_project(project_name,sample_names)
            else:
                project_metadata.update_project(project_name,
                                                sample_names=sample_names)
        # Save
        project_metadata.save(filen)

    def detect_unaligned_dir(self):
        # Attempt to detect an existing 'bcl2fastq' or 'Unaligned' directory
        # containing data from bcl2fastq
        for test_unaligned in ('bcl2fastq','Unaligned'):
            if os.path.isdir(os.path.join(self.analysis_dir,test_unaligned)):
                logging.debug("Testing subdirectory '%s' to see if it has sequence data" % 
                              test_unaligned)
                try:
                    IlluminaData.IlluminaData(self.analysis_dir,
                                              unaligned_dir=test_unaligned)
                    print("Setting 'unaligned_dir' parameter to %s" %
                          test_unaligned)
                    return test_unaligned
                except IlluminaData.IlluminaDataError as ex:
                    logging.debug("Unable to load data from %s" % test_unaligned)
        # Unable to detect existing data directory
        return None

    def set_log_dir(self,path):
        """
        (Re)set the path for the log directory

        If supplied ``path`` is relative then make a
        subdirectory in the existing log directory

        Arguments:
          path (str): path for the log directory

        Returns:
          String: Full path for the new log directory.

        """

        # If the new directory doesn't already exist then
        # create it
        if os.path.isabs(path):
            self._log_dir = path
        else:
            self._log_dir = os.path.join(self.analysis_dir,
                                         self._master_log_dir,
                                         path)
        return self.log_dir

    def log_path(self,*args):
        # Return path appended to log directory
        # Use for getting paths of files under the logs directory
        return os.path.join(self.log_dir,*args)

    def get_log_subdir(self,name):
        """
        Return the name for a new log subdirectory

        Subdirectories are named as NNN_<name>  e.g.
        001_setup, 002_make_fastqs etc

        Arguments:
          name (str): name for the subdirectory
            (typically the name of the processing
            stage that will produce logs to be
            written to the subdirs

        Returns:
          String: name for the new log subdirectory
            (nb not the full path).
        """
        return get_numbered_subdir(
            name,
            parent_dir=os.path.join(self.analysis_dir,self._master_log_dir))

    def __del__(self):
        """
        Implement __del__ method

        Peforms clean up operations (e.g. save parameters,
        remove temporary files etc) when the AutoProcess
        object is destroyed.

        """
        if self.analysis_dir is None:
            return
        try:
            if not os.path.exists(self.analysis_dir):
                logging.warning("Analysis dir '%s' not found" %
                                self.analysis_dir)
                return
            tmp_dir = os.path.join(self.analysis_dir,'tmp')
            if os.path.isdir(tmp_dir):
                logging.debug("Removing %s" % tmp_dir)
                import shutil
                shutil.rmtree(tmp_dir)
            logging.debug("Saving parameters to file")
            self.save_parameters()
            logging.debug("Saving metadata to file")
            self.save_metadata()
        except Exception as ex:
            logging.warning("Exception trying to delete "
                            "AutoProcess instance: %s" %
                            ex)

    @property
    def run_name(self):
        # Return run name
        if self.metadata.run_name is not None:
            return self.metadata.run_name
        elif self.params.data_dir is not None:
            return os.path.basename(self.params.data_dir)
        else:
            run_name = os.path.basename(self.params.analysis_dir)
            if run_name.endswith('_analysis'):
                # Strip trailing _analysis
                run_name = run_name[:-len('_analysis')]
            return run_name

    @property
    def log_dir(self):
        # Generate and return full path to log directory
        return self.add_directory(self._log_dir)

    @property
    def tmp_dir(self):
        # Generate and return full path to tmp directory
        return self.add_directory('tmp')

    @property
    def script_code_dir(self):
        # Generate and return full path to ScriptCode directory
        script_code = self.add_directory('ScriptCode')
        # Put a README file in ScriptCode to make sure it's
        # not pruned on subsequent rsync operations
        readme = os.path.join(script_code,'README.txt')
        if not os.path.exists(readme):
            open(readme,'w').write("The ScriptCode directory is a "
                                   "place to put custom scripts and programs\n")

    @property
    def readme_file(self):
        # If the analysis dir contains a README file then
        # return the full path; otherwise return None
        readme_file = None
        for name in ('README','README.txt'):
            readme_file = os.path.join(self.analysis_dir,name)
            if os.path.isfile(readme_file):
                return readme_file
        # No match found
        return None

    @property
    def run_reference_id(self):
        """
        Return a run reference id (e.g. 'HISEQ_140701/242#22')
        """
        return run_reference_id(
            self.run_name,
            platform=self.metadata.platform,
            facility_run_number=self.metadata.run_number)

    @property
    def parameter_file(self):
        """
        Return name of parameter file ('auto_process.info')
        """
        return os.path.join(self.analysis_dir,'auto_process.info')

    @property
    def has_parameter_file(self):
        """
        Indicate if there is a parameter file (typically auto_process.info)
        """
        return os.path.exists(os.path.join(self.parameter_file))

    @property
    def metadata_file(self):
        """
        Return name of metadata file ('metadata.info')

        """
        return os.path.join(self.analysis_dir,'metadata.info')

    @property
    def paired_end(self):
        """
        Check if run is paired end

        The endedness of the run is checked as follows:

        - If there are analysis project directories then the
          ended-ness is determined by checking the contents of
          these directories
        - If there are no project directories then the
          ended-ness is determined from the contents of the
          'unaligned' directory

        Returns:
          Boolean: True if run is paired end, False if single end,
            None if endedness cannot be determined
        """
        projects = self.get_analysis_projects_from_dirs()
        if projects:
            return reduce(lambda x,y: x and y.info.paired_end,
                          projects,True)
        else:
            try:
                return self.load_illumina_data().paired_end
            except IlluminaData.IlluminaDataError:
                return None

    def print_values(self,data):
        """
        Print key/value pairs from a dictionary

        """
        values = bcf_utils.OrderedDictionary()
        values['Run reference'] = self.run_reference_id
        for i in data:
            values[i] = data[i]
        field_width = max([len(i) for i in values])
        for item in values:
            print("%s: %s" % (item+' '*(field_width-len(item)),
                              values[item]))

    def set_param(self,key,value):
        """
        Set an analysis directory parameter

        Arguments:
          key (str): parameter name
          value (object): value to assign to the parameter

        """
        if key in self.params:
            print("Setting parameter '%s' to '%s'" % (key,value))
            self.params[key] = value
        else:
            raise KeyError("Parameter 'key' not found" % key)

    def print_params(self):
        """
        Print the current parameter settings

        """
        if self.has_parameter_file:
            print("Parameters in %s:" % (os.path.basename(
                self.parameter_file)))
        else:
            print("No parameters file found")
        self.print_values(self.params)

    def set_metadata(self,key,value):
        """
        Set an analysis directory metadata item

        Arguments:
          key (str): parameter name
          value (object): value to assign to the parameter

        """
        if key in self.metadata:
            print("Setting metadata item '%s' to '%s'" % (key,value))
            self.metadata[key] = value
        else:
            raise KeyError("Metadata item 'key' not found" % key)

    def print_metadata(self):
        """
        Print the metadata items and associated values

        """
        if os.path.exists(self.metadata_file):
            print("Metadata in %s:" % (os.path.basename(
                self.metadata_file)))
        else:
            print("No metadata file found")
        self.print_values(self.metadata)

    def make_project_metadata_file(self,project_metadata_file='projects.info'):
        # Generate a project metadata file based on the fastq
        # files and directory structure
        project_metadata = self.load_project_metadata(
            project_metadata_file=project_metadata_file,
            update=True)
        # Save to file
        filen = os.path.join(self.params.analysis_dir,project_metadata_file)
        project_metadata.save(filen)
        self.params['project_metadata'] = project_metadata_file
        print("Saving project metadata to %s" % self.params.project_metadata)

    def get_analysis_projects(self,pattern=None):
        """
        Return the analysis projects in a list

        By default returns all projects within the analysis
        directory (including 'undetermined') which are listed
        in the 'projects.info' metadata file.

        If the 'pattern' is not None then it should be a simple
        pattern used to match against available names to select
        a subset of projects (see bcf_utils.name_matches).

        If any project in 'projects.info' doesn't have a
        matching analysis directory then an exception is
        raised.

        Note:

        - If there is no 'projects.info' file then the projects
          are taken from those in the 'unaligned' directory of
          the analysis directory.
        - If there is no 'unaligned' directory then the projects
          are determined from the subdirectories in the analysis
          directory.

        Arguments:
          pattern (str): optional pattern to select a subset
            of projects (default: select all projects)

        Returns:
          List: list of AnalysisProject instances.
        """
        project_metadata = self.load_project_metadata(
            self.params.project_metadata)
        projects = []
        if pattern is None:
            pattern = '*'
        for line in project_metadata:
            name = line['Project']
            if not bcf_utils.name_matches(name,pattern):
                # Name failed to match, ignore
                continue
            logging.debug("Acquiring data for project %s" % name)
            # Look for a matching project directory
            project_dir = None
            dirs = bcf_utils.list_dirs(self.analysis_dir,startswith=name)
            logging.debug("Possible matching directories: %s" % dirs)
            if len(dirs) == 1:
                # Just a single match
                project_dir = dirs[0]
            else:
                # Multiple matches, look for an exact match
                for d in dirs:
                    if d == name:
                        project_dir = name
                    break
            if project_dir is None:
                logging.error("Unable to resolve directory for project "
                              "'%s'" % name)
                logging.error("Possible dirs: %s" % dirs)
                raise Exception("Unable to resolve directory for project "
                                "'%s'" % name)
            # Attempt to load the project data
            project_dir = os.path.join(self.analysis_dir,project_dir)
            projects.append(AnalysisProject(name,project_dir))
        # Add undetermined reads directory
        if bcf_utils.name_matches('undetermined',pattern):
            undetermined_analysis = self.undetermined()
            if undetermined_analysis is not None and \
               'undetermined' not in [p.name for p in projects]:
                projects.append(undetermined_analysis)
        return projects

    def get_analysis_projects_from_dirs(self,pattern=None):
        """
        Return a list of AnalysisProjects in the analysis directory

        Tests each of the subdirectories in the top-level of the
        analysis directory and rejects any that appear to be
        CASVAVA/bcl2fastq outputs or which don't successfully load
        as AnalysisProject instances.

        Unlike the `get_analysis_projects` method, no checking
        against the project metadata (typically in 'projects.info')
        is performed.

        If the 'pattern' is not None then it should be a simple
        pattern used to match against available names to select
        a subset of projects (see bcf_utils.name_matches).

        Arguments:
          pattern (str): optional pattern to select a subset
            of projects (default: select all projects)

        Returns:
          List: list of AnalysisProject instances.
        """
        logging.debug("Testing subdirectories to determine analysis projects")
        projects = []
        if pattern is None:
            pattern = '*'
        # Try loading each subdirectory as a project
        for dirn in bcf_utils.list_dirs(self.analysis_dir):
            # Test for bcl2fastq output
            try:
                IlluminaData.IlluminaData(self.analysis_dir,
                                          unaligned_dir=dirn)
                logging.debug("* %s: rejected" % dirn)
                continue
            except IlluminaData.IlluminaDataError:
                pass
            except Exception as ex:
                logging.debug("Exception when attempting to load "
                              "subdir '%s' as CASAVA/bcl2fastq output "
                              "(ignored): %s" % (dirn,ex))
            # Try loading as a project
            test_project = AnalysisProject(
                dirn,os.path.join(self.analysis_dir,dirn))
            if test_project.is_analysis_dir:
                logging.debug("* %s: analysis directory" % dirn)
                if bcf_utils.name_matches(test_project.name,
                                              pattern):
                    projects.append(test_project)
            else:
                logging.debug("* %s: rejected" % dirn)
        return projects

    def undetermined(self):
        # Return analysis project directory for undetermined indices
        # or None if not found
        dirs = bcf_utils.list_dirs(self.analysis_dir,matches='undetermined')
        if len(dirs) == 0:
            logging.debug("No undetermined analysis directory found")
            return None
        elif len(dirs) > 1:
            raise Exception("Found multiple undetermined analysis "
                            "directories: %s" % ' '.join(dirs))
        # Attempt to load the analysis project data
        undetermined_dir = os.path.join(self.analysis_dir,dirs[0])
        return AnalysisProject(dirs[0],undetermined_dir)

    def log_analysis(self):
        # Add a record of the analysis to the logging file
        raise NotImplementedError

    def check_metadata(self,items):
        """
        Check that metadata items are set

        For metadata items supplied as an iterable in 'items',
        check that each is set to a non-null value. Report
        those that are null.

        Return False if one or more are null; otherwise return
        True.

        """
        # Check metadata
        metadata_ok = True
        for item in items:
            if item in self.metadata.null_items():
                metadata_ok = False
                logging.warning("Metadata item '%s' is not set" % item)
        return metadata_ok
