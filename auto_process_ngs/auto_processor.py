#!/usr/bin/env python
#
#     auto_processor.py: automated processing of Illumina sequence data
#     Copyright (C) University of Manchester 2013-2018 Peter Briggs
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
import bcftbx.IlluminaData as IlluminaData
import bcftbx.TabFile as TabFile
import bcftbx.utils as bcf_utils
import bcftbx.htmlpagewriter as htmlpagewriter
from bcftbx.JobRunner import fetch_runner
import config
import commands
import applications
import analysis
import metadata
import fileops
import utils
import simple_scheduler
import bcl2fastq_utils
import samplesheet_utils
from icell8.utils import get_icell8_bases_mask
import tenx_genomics_utils
from .settings import Settings
from .qc.processing import report_processing_qc
from .exceptions import MissingParameterFileException
from auto_process_ngs import get_version

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
    ...    print "Hello %s" % cls.person
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
        print "Running '%s' command" % name
        try:
            ret = f(*args,**kws)
        except Exception as ex:
            logging.fatal("%s: %s" % (name,ex))
            ret = 1
        else:
            print "%s: finished" % name
        return ret
    def wrapper(cls):
        # Adds the supplied function to
        # to the class
        setattr(cls,name,wrapped_func)
        return cls
    return wrapper

#######################################################################
# Classes
#######################################################################

@add_command("setup_analysis_dirs",commands.setup_analysis_dirs)
@add_command("run_qc",commands.run_qc)
@add_command("publish_qc",commands.publish_qc_cmd.publish_qc)
@add_command("archive",commands.archive_cmd.archive)
@add_command("report",commands.report_cmd.report)
class AutoProcess:
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
            'settings.ini' file
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
        self.params = metadata.AnalysisDirParameters()
        self.metadata = metadata.AnalysisDirMetadata()
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
            except MissingParameterFileException, ex:
                # No parameter file
                logging.warning("Failed to load parameters: %s (ignored)" % ex)
                logging.warning("Perhaps this is not an auto_process project?")
                # Attempt to detect existing data directory
                self.params['unaligned_dir'] = self.detect_unaligned_dir()
                if self.params.unaligned_dir is None:
                    logging.warning("Unable to find subdirectory containing data")
            except Exception, ex:
                logging.error("Failed to load parameters: %s" % ex)
                logging.error("Stopping")
                sys.exit(1)
            self.params['analysis_dir'] = self.analysis_dir
            # Load metadata
            try:
                self.load_metadata(allow_save=allow_save_params)
            except MissingParameterFileException, ex:
                # No metadata file
                logging.warning("Failed to load metadata: %s (ignored)" % ex)
                logging.warning("Consider running metadata --update?")
            except Exception, ex:
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
                    print "Making %s" % dir_path
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
            print "Migrating metadata values from parameter file"
            for param in ('platform','run_number','source','assay'):
                if param not in self.params:
                    continue
                if self.metadata[param] is None:
                    logging.debug("Importing metadata item '%s': set to "
                              "'%s'" % (param,self.params[param]))
                    print "Importing metadata item '%s'" % param
                    self.metadata[param] = self.params[param]
        # Run name
        if self.metadata.run_name is None:
            print "Attempting to set missing 'run_name' metadata item"
            self.metadata['run_name'] = self.run_name
        # Instrument-related metadata
        if self.metadata.instrument_name is None or \
           self.metadata.instrument_datestamp is None or \
           self.metadata.instrument_run_number is None:
            print "Attempting to set missing instrument metadata items"
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
            platform = bcl2fastq_utils.get_sequencer_platform(
                self.analysis_dir,
                instrument=self.metadata.instrument_name,
                settings=self.settings)
            if platform:
                print "Setting 'platform' metadata item to %s" % \
                    platform
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
        utils.edit_file(sample_sheet_file)
        # Check updated sample sheet and issue warnings
        if samplesheet_utils.check_and_warn(sample_sheet_file=sample_sheet_file):
            logging.error("Sample sheet may have problems, see warnings above")

    def init_readme(self):
        """
        Create a new README file
        """
        if self.readme_file is None:
            readme_file = os.path.join(self.analysis_dir,'README')
            print "Initialising %s" % readme_file
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
        utils.edit_file(self.readme_file,
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
        except IlluminaData.IlluminaDataError,ex:
            logging.warning("Failed to load data from bcl2fastq output "
                            "(ignored): %s" % ex)
            illumina_data = None
        projects_from_dirs = self.get_analysis_projects_from_dirs()
        if filen is not None and os.path.exists(filen):
            # Load existing file and check for consistency
            logging.debug("Loading project metadata from existing file")
            project_metadata = metadata.ProjectMetadataFile(filen)
        else:
            # First try to populate basic metadata from existing projects
            logging.debug("Metadata file not found, guessing basic data")
            project_metadata = metadata.ProjectMetadataFile()
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
                test_project = analysis.AnalysisProject(
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
        print "Project metadata file: %s" % self.params.project_metadata
        filen = os.path.join(self.analysis_dir,
                             self.params.project_metadata)
        if unaligned_dir is not None:
            self.params['unaligned_dir'] = unaligned_dir
        print "Unaligned_dir: %s" % self.params.unaligned_dir
        illumina_data = IlluminaData.IlluminaData(
            self.analysis_dir,
            unaligned_dir=self.params.unaligned_dir)
        if os.path.exists(filen):
            # Load data from existing file
            print "Loading project metadata from existing file: %s" % filen
            project_metadata = metadata.ProjectMetadataFile(filen)
        else:
            # New (empty) metadata file
            print "Creating new project metadata file: %s" % filen
            project_metadata = metadata.ProjectMetadataFile()
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
                    print "Setting 'unaligned_dir' parameter to %s" % test_unaligned
                    return test_unaligned
                except IlluminaData.IlluminaDataError, ex:
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
        return utils.get_numbered_subdir(
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
        """Return a run reference id e.g. 'HISEQ_140701/242#22'

        The run reference code is a code that identifies the sequencing
        run, and has the general form:

        PLATFORM_DATESTAMP[/INSTRUMENT_RUN_NUMBER]#FACILITY_RUN_NUMBER

        - PLATFORM is always uppercased e.g. HISEQ, MISEQ, GA2X
        - DATESTAMP is the YYMMDD code e.g. 140701
        - INSTRUMENT_RUN_NUMBER is the run number that forms part of the
          run directory e.g. for '140701_SN0123_0045_000000000-A1BCD'
          it is '45'
        - FACILITY_RUN_NUMBER is the run number that has been assigned
          by the facility

        Note that the instrument run number is only used if it differs
        from the facility run number.

        """
        # Extract information from run name
        run_name = self.run_name
        try:
            datestamp,instrument,run_number = IlluminaData.split_run_name(run_name)
        except Exception, ex:
            logging.warning("Unable to extract information from run name '%s'" \
                            % run_name)
            logging.warning("Exception: %s" % ex)
            instrument = None
            date_stamp = None
            run_number = None
        # Update items from metadata
        if self.metadata.instrument_name:
            instrument = self.metadata.instrument_name
        if self.metadata.instrument_datestamp:
            date_stamp = self.metadata.instrument_datestamp
        if self.metadata.instrument_run_number:
            run_number = self.metadata.instrument_run_number
        if run_number is not None:
            run_number = str(run_number).lstrip('0')
        # Platform
        if self.metadata.platform is not None:
            platform = self.metadata.platform.upper()
        else:
            platform = None
        # Facility run number
        if self.metadata.run_number is not None:
            facility_run_number = str(self.metadata.run_number)
        else:
            facility_run_number = None
        # Construct the reference id
        if platform is not None:
            run_id = platform
            if datestamp is not None:
                run_id += "_%s" % datestamp
            if run_number is not None:
                try:
                    if run_number != facility_run_number:
                        run_id += "/%s" % run_number
                except ValueError:
                    run_id += "/%s" % run_number
            if facility_run_number is not None:
                run_id += "#%s" % facility_run_number
        else:
            run_id = run_name
        return run_id

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

    def setup(self,data_dir,analysis_dir=None,sample_sheet=None):
        """
        Set up the initial analysis directory

        This does all the initialisation of the analysis directory
        and processing parameters

        Arguments:
          data_dir (str): source data directory
          analysis_dir (str): corresponding analysis directory
          sample_sheet (str): name and location of non-default sample sheet
            file

        """
        data_dir = data_dir.rstrip(os.sep)
        if not fileops.exists(data_dir):
            raise Exception("Data directory '%s' not found" %
                            data_dir)
        if not fileops.Location(data_dir).is_remote:
            data_dir = os.path.abspath(data_dir)
        run_name = os.path.basename(data_dir)
        if analysis_dir is None:
            analysis_dir = os.path.join(
                os.getcwd(),run_name)+'_analysis'
        else:
            analysis_dir = os.path.abspath(analysis_dir)
        # Create the analysis directory structure
        if not os.path.exists(analysis_dir):
            # Make a temporary analysis dir
            tmp_analysis_dir = os.path.join(
                os.path.dirname(analysis_dir),
                ".%s.%s" % (os.path.basename(analysis_dir),
                            uuid.uuid4()))
            self.analysis_dir = tmp_analysis_dir
            logging.debug("Creating temp directory '%s'" %
                          self.analysis_dir)
            # Create directory structure
            self.create_directory(self.analysis_dir)
            self.log_dir
            self.script_code_dir
        else:
            # Directory already exists
            logging.warning("Analysis directory '%s' already exists" %
                            analysis_dir)
            self.analysis_dir = analysis_dir
            # check for parameter file
            if self.has_parameter_file:
                self.load_parameters()
            else:
                logging.warning("No parameter file found in %s" %
                                self.analysis_dir)
        # Run datestamp, instrument name and instrument run number
        try:
            datestamp,instrument,run_number,flow_cell_prefix,flow_cell_id = \
                IlluminaData.split_run_name_full(run_name)
            run_number = run_number.lstrip('0')
            flow_cell = flow_cell_prefix + flow_cell_id
        except Exception as ex:
            logging.warning("Unable to extract information from run name '%s'" \
                            % run_name)
            logging.warning("Exception: %s" % ex)
            datestamp = None
            instrument= None
            run_number = None
            flow_cell = None
        # Identify missing data and attempt to acquire
        # Sequencing platform
        platform = self.metadata.platform
        if platform is None:
            platform = bcl2fastq_utils.get_sequencer_platform(
                data_dir,
                instrument=instrument,
                settings=self.settings)
        print "Platform identified as '%s'" % platform
        # Log dir
        self.set_log_dir(self.get_log_subdir('setup'))
        # Custom SampleSheet.csv file
        custom_sample_sheet = self.params.sample_sheet
        if custom_sample_sheet is None:
            print "Acquiring sample sheet..."
            if sample_sheet is None:
                targets = ('Data/Intensities/BaseCalls/SampleSheet.csv',
                           'SampleSheet.csv',)
            else:
                targets = (sample_sheet,)
            # Try each possibility until one sticks
            for target in targets:
                target = fileops.Location(target)
                if target.is_remote:
                    sample_sheet = str(target)
                else:
                    if os.path.isabs(target.path):
                        sample_sheet = target.path
                    else:
                        sample_sheet = os.path.join(data_dir,
                                                    target.path)
                print "Trying '%s'" % sample_sheet
                tmp_sample_sheet = os.path.join(self.tmp_dir,
                                                os.path.basename(target.path))
                rsync = applications.general.rsync(sample_sheet,
                                                   self.tmp_dir)
                print "%s" % rsync
                status = rsync.run_subprocess(log=self.log_path('rsync.sample_sheet.log'))
                if status != 0:
                    logging.warning("Failed to fetch sample sheet '%s'"
                                    % sample_sheet)
                    tmp_sample_sheet = None
                else:
                    break
            # Bail out if no sample sheet was acquired
            if tmp_sample_sheet is None:
                shutil.rmtree(self.analysis_dir)
                self.analysis_dir = None
                raise Exception("Unable to acquire sample sheet")
            # Keep a copy of the original sample sheet
            original_sample_sheet = os.path.join(self.analysis_dir,
                                                 'SampleSheet.orig.csv')
            print "Copying original sample sheet to %s" % original_sample_sheet
            shutil.copyfile(tmp_sample_sheet,original_sample_sheet)
            # Set the permissions for the original SampleSheet
            os.chmod(original_sample_sheet,0664)
            # Process acquired sample sheet
            custom_sample_sheet = os.path.join(self.analysis_dir,
                                               'custom_SampleSheet.csv')
            sample_sheet = bcl2fastq_utils.make_custom_sample_sheet(
                tmp_sample_sheet,custom_sample_sheet)
        else:
            sample_sheet = IlluminaData.SampleSheet(custom_sample_sheet)
            original_sample_sheet = os.path.join(self.analysis_dir,
                                                 'SampleSheet.orig.csv')
        print "Sample sheet '%s'" % custom_sample_sheet
        # Library Prep Kit/Assay
        assay = None
        for item in ('Assay','Library Prep Kit'):
            try:
                assay = IlluminaData.SampleSheet(original_sample_sheet).header[item]
                break
            except KeyError:
                logging.warning("No element '%s' found in sample sheet"
                                % item)
        # Bases mask
        print "Bases mask set to 'auto' (will be determined at run time)"
        bases_mask = "auto"
        # Generate and print predicted outputs and warnings
        print samplesheet_utils.predict_outputs(sample_sheet=sample_sheet)
        samplesheet_utils.check_and_warn(sample_sheet=sample_sheet)
        # Move analysis dir to final location if necessary
        if self.analysis_dir != analysis_dir:
            logging.debug("Moving %s to final directory" % self.analysis_dir)
            os.rename(self.analysis_dir,analysis_dir)
            self.analysis_dir = analysis_dir
            # Update the custom sample sheet path
            custom_sample_sheet = os.path.join(
                analysis_dir,
                os.path.basename(custom_sample_sheet))
            print "Created analysis directory '%s'" % self.analysis_dir
        # Store the parameters
        self.params['data_dir'] = data_dir
        self.params['analysis_dir'] = self.analysis_dir
        self.params['sample_sheet'] = custom_sample_sheet
        self.params['bases_mask'] = bases_mask
        # Store the metadata
        self.metadata['run_name'] = self.run_name
        self.metadata['platform'] = platform
        self.metadata['instrument_name'] = instrument
        self.metadata['instrument_datestamp'] = datestamp
        self.metadata['instrument_run_number'] = run_number
        self.metadata['instrument_flow_cell_id'] = flow_cell
        self.metadata['assay'] = assay
        # Set flags to allow parameters etc to be saved back
        self._save_params = True
        self._save_metadata = True

    def setup_from_fastq_dir(self,analysis_dir,fastq_dir):
        # Do setup for an existing directory containing fastq files
        # with the same structure as that produced by CASAVA and bcl2fastq
        #
        # Assumes that the files are in a subdirectory of the analysis
        # directory specified by the 'fastq_dir' argument, and
        # that within that they are arranged in the structure
        # 'Project_<name>/Sample_<name>/<fastq>'
        self.analysis_dir = os.path.abspath(analysis_dir)
        # Create directory structure
        self.create_directory(self.analysis_dir)
        # Get information
        print "Identifying platform from data directory name"
        platform = bcl2fastq_utils.get_sequencer_platform(analysis_dir)
        # Store the parameters
        self.params['analysis_dir'] = self.analysis_dir
        self.params['unaligned_dir'] = fastq_dir
        # Store the metadata
        self.metadata['platform'] = platform
        # Set flag to allow parameters etc to be saved back
        self._save_params = True
        self._save_metadata = True
        # Generate statistics
        self.generate_stats()
        # Make a 'projects.info' metadata file
        self.make_project_metadata_file()

    def clone(self,clone_dir,copy_fastqs=False):
        """
        Make a 'clone' (i.e. copy) to new directory

        Makes a copy of the processing directory, including
        metadata and parameters but ignoring any project
        subdirectories.

        By default the 'unaligned' directory in the new
        directory is simply a symlink from the original
        directory; set the 'copy_fastqs' to make copies
        instead.

        Arguments
          clone_dir (str): path to the new directory to
            create as a clone (must not already exist).
          copy_fastqs (boolean): set to True to copy the
            fastq files (otherwise default behaviour is
            to make a symlink)

        """
        clone_dir = os.path.abspath(clone_dir)
        if os.path.exists(clone_dir):
            # Directory already exists
            logging.critical("Target directory '%s' already exists" % clone_dir)
            raise Exception("Clone failed, target directory already exists")
        self.create_directory(clone_dir)
        # Copy metadata and parameters
        self.save_metadata(os.path.join(clone_dir,
                                        os.path.basename(
                                            self.metadata_file)),
                           force=True)
        self.save_parameters(os.path.join(clone_dir,
                                          os.path.basename(
                                              self.parameter_file)),
                             force=True)
        # Link to or copy fastqs
        unaligned_dir = os.path.join(self.analysis_dir,self.params.unaligned_dir)
        if os.path.isdir(unaligned_dir):
            clone_unaligned_dir = os.path.join(clone_dir,
                                               os.path.basename(
                                                   self.params.unaligned_dir))
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
        projects = self.get_analysis_projects()
        if projects:
            print "Duplicating project directories:"
            for project in self.get_analysis_projects():
                print "-- %s" % project.name
                fastqs = project.fastqs
                new_project = analysis.AnalysisProject(
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
                  self.params.sample_sheet,
                  self.params.stats_file,
                  self.params.project_metadata):
            srcpath = os.path.join(self.analysis_dir,f)
            if os.path.exists(srcpath):
                shutil.copy(srcpath,clone_dir)
        # Create the basic set of subdirectories
        d = AutoProcess(analysis_dir=clone_dir)
        d.log_dir
        d.script_code_dir
        # Update the settings
        for s in ("sample_sheet",):
            d.params[s] = os.path.join(d.analysis_dir,
                                       os.path.relpath(d.params[s],
                                                       self.analysis_dir))
        d.save_parameters()

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
            print "%s: %s" % (item+' '*(field_width-len(item)),
                              values[item])

    def set_param(self,key,value):
        """
        Set an analysis directory parameter

        Arguments:
          key (str): parameter name
          value (object): value to assign to the parameter

        """
        if key in self.params:
            print "Setting parameter '%s' to '%s'" % (key,value)
            self.params[key] = value
        else:
            raise KeyError("Parameter 'key' not found" % key)

    def print_params(self):
        """
        Print the current parameter settings

        """
        if self.has_parameter_file:
            print "Parameters in %s:" % (os.path.basename(self.parameter_file))
        else:
            print "No parameters file found"
        self.print_values(self.params)

    def set_metadata(self,key,value):
        """
        Set an analysis directory metadata item

        Arguments:
          key (str): parameter name
          value (object): value to assign to the parameter

        """
        if key in self.metadata:
            print "Setting metadata item '%s' to '%s'" % (key,value)
            self.metadata[key] = value
        else:
            raise KeyError("Metadata item 'key' not found" % key)

    def print_metadata(self):
        """
        Print the metadata items and associated values

        """
        if os.path.exists(self.metadata_file):
            print "Metadata in %s:" % (os.path.basename(self.metadata_file))
        else:
            print "No metadata file found"
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
        print "Saving project metadata to %s" % self.params.project_metadata

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
            projects.append(analysis.AnalysisProject(name,project_dir))
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
            test_project = analysis.AnalysisProject(
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
            raise Exception, "Found multiple undetermined analysis directories: %s" \
                % ' '.join(dirs)
        # Attempt to load the analysis project data
        undetermined_dir = os.path.join(self.analysis_dir,dirs[0])
        return analysis.AnalysisProject(dirs[0],undetermined_dir)
        
    def make_fastqs(self,protocol='standard',platform=None,
                    unaligned_dir=None,sample_sheet=None,lanes=None,
                    ignore_missing_bcl=False,ignore_missing_stats=False,
                    skip_rsync=False,remove_primary_data=False,
                    nprocessors=None,require_bcl2fastq_version=None,
                    bases_mask=None,no_lane_splitting=None,
                    minimum_trimmed_read_length=None,
                    mask_short_adapter_reads=None,
                    generate_stats=True,stats_file=None,
                    per_lane_stats_file=None,
                    skip_fastq_generation=False,
                    only_fetch_primary_data=False,
                    create_empty_fastqs=None,runner=None,
                    cellranger_jobmode="local",
                    cellranger_mempercore=None,
                    cellranger_maxjobs=None,
                    cellranger_jobinterval=None,
                    cellranger_localcores=None,
                    cellranger_localmem=None,
                    cellranger_ignore_dual_index=False):
        """Create and summarise FASTQ files

        Wrapper for operations related to FASTQ file generation and analysis.
        The operations are typically:
 
        - get primary data (BCL files)
        - run bcl-to-fastq conversion
        - generate statistics

        If the number of processors and the job runner are not explicitly
        specified then these are taken from the settings for the bcl2fastq
        and the statistics generation steps, which may differ from each other.
        However if either of these values are set explicitly then the same
        values will be used for both steps.

        Arguments:
          protocol            : if set then specifies the protocol to use
                                for fastq generation, otherwise use the
                                'standard' bcl2fastq protocol
          platform            : if set then specifies the sequencing platform
                                (otherwise platform will be determined from the
                                primary data)
          unaligned_dir       : if set then use this as the output directory for
                                bcl-to-fastq conversion. Default is 'bcl2fastq' (unless
                                an alternative is already specified in the config file)
          sample_sheet        : if set then use this as the input samplesheet
          lanes               : (optional) specify a list of lane numbers to
                                use in the processing; lanes not in the list
                                will be excluded (default is to include all
                                lanes)
          nprocessors         : number of processors to run bclToFastq.py with
          ignore_missing_bcl  : if True then run bcl2fastq with --ignore-missing-bcl
          ignore_missing_stats: if True then run bcl2fastq with --ignore-missing-stats
          skip_rsync          : if True then don't rsync primary data at the start of
                                bcl2fastq conversion
          remove_primary_data : if True then remove primary data at the end of bcl2fastq
                                conversion (default is to keep it)
          generate_stats      : if True then (re)generate statistics file for fastqs
          require_bcl2fastq_version: (optional) specify bcl2fastq version to use
                                Should be a string of the form '1.8.4' or
                                '>2.0'. Set to None to automatically determine
                                bcl2fastq version.
          bases_mask          : if set then use this as an alternative bases mask setting
          no_lane_splitting   : if True then run bcl2fastq with --no-lane-splitting
          minimum_trimmed_read_length: if set then specify minimum length
                                for reads after adapter trimming (shorter reads will
                                be padded with Ns to make them long enough)
          mask_short_adapter_reads: if set then specify the minimum length
                                of ACGT bases that must be present in a read after
                                adapter trimming for it not to be masked completely
                                with Ns.
          stats_file          : if set then use this as the name of the output
                                per-fastq stats file.
          per_lane_stats_file : if set then use this as the name of the output
                                per-lane stats file.
          skip_fastq_generation: if True then don't perform fastq generation
          only_fetch_primary_data: if True then fetch primary data, don't do anything else
          create_empty_fastqs : if True then create empty 'placeholder' fastq
                                files for any missing fastqs after bcl2fastq
                                (must have completed with zero exit status)
          runner              : (optional) specify a non-default job runner to use for
                                fastq generation
          cellranger_jobmode  : (optional) job mode to run cellranger in; defaults
                                to "local" (10xGenomics Chromium SC data only)
          cellranger_mempercore: (optional) memory assumed per core (in Gbs)
                                (10xGenomics Chromium SC data only)
          cellranger_maxjobs  : (optional) maxiumum number of concurrent jobs
                                to run (10xGenomics Chromium SC data only)
          cellranger_jobinterval: (optional) how often jobs are submitted (in
                                ms) (10xGenomics Chromium SC data only)
          cellranger_localcores: (optional) maximum number of cores cellranger
                                can request in jobmode 'local' (10xGenomics Chromium
                                SC data only)
          cellranger_localmem : (optional) maximum memory cellranger can request in
                                jobmode 'local' (10xGenomics Chromium SC data
                                only)
          cellranger_ignore_dual_index: (optional) on a dual-indexed flowcell where
                                the second index was not used for the 10x sample,
                                ignore it (10xGenomics Chromium SC data only)
        """
        # Report protocol
        print "Protocol              : %s" % protocol
        if protocol not in ('standard','icell8','10x_chromium_sc'):
            raise Exception("Unknown protocol: '%s'" % protocol)
        # Unaligned dir
        if unaligned_dir is not None:
            self.params['unaligned_dir'] = unaligned_dir
        elif self.params['unaligned_dir'] is None:
            self.params['unaligned_dir'] = 'bcl2fastq'
        print "Output dir            : %s" % self.params.unaligned_dir
        # Sample sheet
        if sample_sheet is None:
            sample_sheet = self.params.sample_sheet
        if not os.path.isabs(sample_sheet):
            sample_sheet = os.path.join(self.analysis_dir,sample_sheet)
        if not os.path.isfile(sample_sheet):
            raise Exception("Missing sample sheet '%s'" % sample_sheet)
        self.params['sample_sheet'] = sample_sheet
        print "Source sample sheet   : %s" % self.params.sample_sheet
        # Check requested lanes are actually present
        print "Lanes                 : %s" % ('all' if lanes is None
                                              else
                                              ','.join([str(l)
                                                        for l in lanes]))
        if lanes is not None:
            s = IlluminaData.SampleSheet(self.params.sample_sheet)
            if not s.has_lanes:
                raise Exception("Requested subset of lanes but "
                                "samplesheet doesn't contain any "
                                "lane information")
            samplesheet_lanes = list(set([l['Lane'] for l in s]))
            for l in lanes:
                if l not in samplesheet_lanes:
                    raise Exception("Requested lane '%d' not present "
                                    "in samplesheet" % l)
        # Make a temporary sample sheet
        if lanes:
            lanes_id = ".L%s" % ''.join([str(l) for l in lanes])
        else:
            lanes_id = ""
        sample_sheet = os.path.join(self.tmp_dir,
                                    "SampleSheet%s.%s.csv" %
                                    (lanes_id,
                                     time.strftime("%Y%m%d%H%M%S")))
        bcl2fastq_utils.make_custom_sample_sheet(self.params.sample_sheet,
                                                 sample_sheet,
                                                 lanes=lanes)
        # Check the temporary sample sheet
        print "Checking temporary sample sheet"
        invalid_barcodes = samplesheet_utils.SampleSheetLinter(
            sample_sheet_file=sample_sheet).has_invalid_barcodes()
        if invalid_barcodes:
            logging.error("Invalid barcodes detected")
            for line in invalid_barcodes:
                logging.critical("%s" % line)
        invalid_characters = samplesheet_utils.SampleSheetLinter(
            sample_sheet_file=sample_sheet).has_invalid_characters()
        if invalid_characters:
            logging.critical("Invalid non-printing/non-ASCII characters "
                             "detected")
        if invalid_barcodes or invalid_characters:
            raise Exception("Errors detected in generated sample sheet")
        # Adjust verification settings for 10xGenomics Chromium SC
        # data if necessary
        verify_include_sample_dir = False
        if tenx_genomics_utils.has_chromium_sc_indices(sample_sheet):
            if protocol == '10x_chromium_sc':
                # Force inclusion of sample-name subdirectories
                # when verifying Chromium SC data
                print "Sample sheet includes Chromium SC indices"
                verify_include_sample_dir = True
            else:
                # Chromium SC indices detected but not using
                # 10x_chromium_sc protocol
                raise Exception("Detected 10xGenomics Chromium SC indices "
                                "in generated sample sheet but protocol "
                                "'%s' has been specified; must use "
                                "'10x_chromium_sc' for these indices" %
                                protocol)
        # Check for pre-existing Fastq outputs
        if self.verify_fastq_generation(
                unaligned_dir=self.params.unaligned_dir,
                lanes=lanes,
                include_sample_dir=verify_include_sample_dir):
            print "Expected Fastq outputs already present"
            skip_rsync = True
            skip_fastq_generation = True
        # Check if there's anything to do
        if (skip_rsync and skip_fastq_generation) and not generate_stats:
            print "Nothing to do"
            return
        # Log dir
        log_dir = 'make_fastqs'
        if protocol != 'standard':
            log_dir += "_%s" % protocol
        if lanes:
            log_dir += "_L%s" % ''.join([str(l) for l in sorted(lanes)])
        self.set_log_dir(self.get_log_subdir(log_dir))
        # Fetch primary data
        if not skip_rsync:
            if self.get_primary_data() != 0:
                logging.error("Failed to acquire primary data")
                raise Exception, "Failed to acquire primary data"
            if only_fetch_primary_data:
                return
        # Deal with platform information
        if not platform:
            platform = self.metadata.platform
        # Do fastq generation using the specified protocol
        if not skip_fastq_generation:
            # Set primary data location and report info
            primary_data_dir = os.path.join(
                self.params.primary_data_dir,
                os.path.basename(self.params.data_dir))
            print "Primary data dir      : %s" % primary_data_dir
            try:
                illumina_run = IlluminaData.IlluminaRun(primary_data_dir,
                                                        platform=platform)
            except IlluminaData.IlluminaDataPlatformError as ex:
                logging.critical("Error loading primary data: %s" % ex)
                if platform is None:
                    logging.critical("Try specifying platform using --platform?")
                else:
                    logging.critical("Check specified platform is valid (or "
                                     "omit --platform")
                raise Exception("Error determining sequencer platform")
            print "Platform              : %s" % illumina_run.platform
            print "Bcl format            : %s" % illumina_run.bcl_extension
            # Set platform in metadata
            self.metadata['platform'] = illumina_run.platform
            # Bases mask
            if bases_mask is None:
                bases_mask = self.params.bases_mask
            print "Bases mask setting    : %s" % bases_mask
            if protocol != '10x_chromium_sc':
                if bases_mask == "auto":
                    print "Determining bases mask from RunInfo.xml"
                    bases_mask = bcl2fastq_utils.get_bases_mask(
                        illumina_run.runinfo_xml,
                        sample_sheet)
                    if not bcl2fastq_utils.bases_mask_is_valid(bases_mask):
                        raise Exception("Invalid bases mask: '%s'" %
                                        bases_mask)
            self.params.bases_mask = bases_mask
            # Do fastq generation according to protocol
            if protocol == 'icell8':
                # ICell8 data
                # Update bcl2fastq settings appropriately
                print "Updating read trimming and masking for ICell8"
                minimum_trimmed_read_length = 21
                mask_short_adapter_reads = 0
                # Reset the default bases mask
                bases_mask = IlluminaData.IlluminaRunInfo(
                    illumina_run.runinfo_xml).bases_mask
                bases_mask = get_icell8_bases_mask(
                    bases_mask,
                    sample_sheet=sample_sheet)
                if not bcl2fastq_utils.bases_mask_is_valid(bases_mask):
                    raise Exception("Invalid bases mask: '%s'" %
                                    bases_mask)
                self.params.bases_mask = bases_mask
                # Switch to standard protocol
                protocol = 'standard'
            if protocol == 'standard':
                # Standard protocol
                try:
                    exit_code = self.bcl_to_fastq(
                        unaligned_dir=self.params.unaligned_dir,
                        sample_sheet=sample_sheet,
                        primary_data_dir=primary_data_dir,
                        require_bcl2fastq=require_bcl2fastq_version,
                        bases_mask=bases_mask,
                        ignore_missing_bcl=ignore_missing_bcl,
                        ignore_missing_stats=ignore_missing_stats,
                        no_lane_splitting=no_lane_splitting,
                        minimum_trimmed_read_length=minimum_trimmed_read_length,
                        mask_short_adapter_reads=mask_short_adapter_reads,
                        nprocessors=nprocessors,
                        runner=runner)
                except Exception,ex:
                    raise Exception("Bcl2fastq stage failed: '%s'" % ex)
            elif protocol == '10x_chromium_sc':
                # 10xGenomics Chromium SC
                if bases_mask == 'auto':
                    bases_mask = None
                try:
                    # Check we have cellranger
                    cellranger = bcf_utils.find_program('cellranger')
                    if not cellranger:
                        raise Exception("No cellranger package found")
                    print "Using cellranger %s: %s" % (
                        tenx_genomics_utils.cellranger_info(cellranger)[-1],
                        cellranger)
                    # Check we have bcl2fastq
                    bcl2fastq = bcf_utils.find_program('bcl2fastq')
                    if not bcl2fastq:
                        raise Exception("No bcl2fastq package found")
                    bcl2fastq = bcl2fastq_utils.available_bcl2fastq_versions(
                        paths=(os.path.dirname(bcl2fastq),),
                        reqs='>=2.17')
                    if not bcl2fastq:
                        raise Exception("No appropriate bcl2fastq software "
                                        "located")
                    bcl2fastq = bcl2fastq[0]
                    bcl2fastq_info = bcl2fastq_utils.bcl_to_fastq_info(
                        bcl2fastq)
                    print "Using bcl2fastq %s: %s" % (bcl2fastq_info[-1],
                                                      bcl2fastq)
                    # Store info on bcl2fastq package
                    self.metadata['bcl2fastq_software'] = bcl2fastq_info
                    # Put a copy of sample sheet in the log directory
                    shutil.copy(sample_sheet,self.log_dir)
                    # Run cellranger mkfastq
                    exit_code = tenx_genomics_utils.run_cellranger_mkfastq(
                        sample_sheet=sample_sheet,
                        primary_data_dir=primary_data_dir,
                        output_dir=self.params.unaligned_dir,
                        lanes=(None if lanes is None
                               else ','.join([str(l) for l in lanes])),
                        bases_mask=bases_mask,
                        ignore_dual_index=cellranger_ignore_dual_index,
                        cellranger_jobmode=cellranger_jobmode,
                        cellranger_maxjobs=cellranger_maxjobs,
                        cellranger_mempercore=cellranger_mempercore,
                        cellranger_jobinterval=cellranger_jobinterval,
                        cellranger_localcores=cellranger_localcores,
                        cellranger_localmem=cellranger_localmem,
                        log_dir=self.log_dir)
                except Exception,ex:
                    raise Exception("'cellranger mkfastq' stage failed: "
                                    "'%s'" % ex)
            else:
                # Unknown protocol
                raise Exception("Unknown protocol '%s'" % protocol)
            # Check the outputs
            if exit_code != 0:
                raise Exception("Fastq generation finished with error: "
                                "exit code %d" % exit_code)
            if not self.verify_fastq_generation(
                    lanes=lanes,
                    include_sample_dir=verify_include_sample_dir):
                # Check failed
                logging.error("Failed to verify output Fastqs against "
                              "sample sheet")
                # Try to load the data from unaligned dir
                try:
                    illumina_data = IlluminaData.IlluminaData(
                        self.analysis_dir,
                        unaligned_dir=self.params.unaligned_dir)
                except IlluminaData.IlluminaDataError as ex:
                    raise Exception("Unable to load data from %s: %s"
                                    % (self.params.unaligned_dir,ex))
                # Generate a list of missing Fastqs
                missing_fastqs = IlluminaData.list_missing_fastqs(
                    illumina_data,
                    sample_sheet,
                    include_sample_dir=verify_include_sample_dir)
                assert(len(missing_fastqs) > 0)
                missing_fastqs_file = os.path.join(self.log_dir,
                                                   "missing_fastqs.log")
                print "Writing list of missing Fastq files to %s" % \
                    missing_fastqs_file
                with open(missing_fastqs_file,'w') as fp:
                    for fq in missing_fastqs:
                        fp.write("%s\n" % fq)
                # Create empty FASTQs
                if create_empty_fastqs is None:
                    try:
                        create_empty_fastqs = \
                            self.settings.platform[self.metadata.platform].\
                            create_empty_fastqs
                    except (KeyError,AttributeError):
                        create_empty_fastqs = \
                            self.settings.bcl2fastq.create_empty_fastqs
                if create_empty_fastqs:
                    logging.warning("Making 'empty' placeholder Fastqs")
                    for fq in missing_fastqs:
                        fastq = os.path.join(self.analysis_dir,
                                             self.params.unaligned_dir,fq)
                        print "-- %s" % fastq
                        if not os.path.exists(os.path.dirname(fastq)):
                            bcf_utils.mkdirs(os.path.dirname(fastq))
                        with gzip.GzipFile(filename=fastq,mode='wb') as fp:
                            fp.write('')
                else:
                    raise Exception("Fastq generation failed to produce "
                                    "expected outputs")
        # Generate statistics
        if generate_stats:
            self.generate_stats(stats_file=stats_file,
                                per_lane_stats_file=per_lane_stats_file,
                                unaligned_dir=self.params.unaligned_dir,
                                nprocessors=nprocessors,
                                runner=runner)
        # Make a 'projects.info' metadata file
        if lanes:
            self.update_project_metadata_file()
        else:
            self.make_project_metadata_file()
        # Remove primary data
        if remove_primary_data:
            self.remove_primary_data()

    def get_primary_data(self,runner=None):
        """Acquire the primary sequencing data (i.e. BCL files)

        Copies the primary sequencing data (bcl files etc) to a local area
        using rsync.

        Arguments:
          runner: (optional) specify a non-default job runner to use
            for primary data rsync

        """
        # Source and target directories
        data_dir = self.params.data_dir
        self.params["primary_data_dir"] = self.add_directory('primary_data')
        # Set up runner
        if runner is not None:
            runner = fetch_runner(runner)
        else:
            runner = self.settings.runners.rsync
        runner.set_log_dir(self.log_dir)
        # Run rsync command
        rsync = applications.general.rsync(data_dir,
                                           self.params.primary_data_dir,
                                           prune_empty_dirs=True,
                                           extra_options=('--copy-links',
                                                          '--include=*/',
                                                          '--include=Data/**',
                                                          '--include=RunInfo.xml',
                                                          '--include=SampleSheet.csv',
                                                          '--include=RTAComplete.txt',
                                                          '--include=runParameters.xml',
                                                          '--include=RunParameters.xml',
                                                          '--exclude=*'))
        print "Running %s" % rsync
        rsync_job = simple_scheduler.SchedulerJob(runner,
                                                  rsync.command_line,
                                                  name='rsync.primary_data',
                                                  working_dir=os.getcwd())
        rsync_job.start()
        try:
            rsync_job.wait()
        except KeyboardInterrupt:
            logging.warning("Keyboard interrupt, terminating primary data "
                            "rsync operation")
            rsync_job.terminate()
            return -1
        exit_code = rsync_job.exit_code
        print "rsync of primary data completed: exit code %s" % exit_code
        if exit_code != 0:
            logging.error("Failed to acquire primary data (non-zero "
                          "exit code returned)")
        return exit_code

    def bcl_to_fastq(self,unaligned_dir,sample_sheet,primary_data_dir,
                     require_bcl2fastq=None,bases_mask=None,
                     ignore_missing_bcl=False,ignore_missing_stats=False,
                     no_lane_splitting=None,minimum_trimmed_read_length=None,
                     mask_short_adapter_reads=None,nprocessors=None,
                     runner=None):
        """Generate FASTQ files from the raw BCL files

        Performs FASTQ generation from raw BCL files produced by an Illumina
        sequencer, by running the external 'bclToFastq.py' program (which
        wraps the 'configureBclToFastq' and 'make' steps).

        Arguments:
          unaligned_dir: output directory for bcl-to-fastq conversion
          sample_sheet: input sample sheet file
          primary_data_dir: path to the top-level directory holding
            the sequencing data
          require_bcl2fastq: if set then should be a string of the form
            '1.8.4' or '>2.0' explicitly specifying the version of
            bcl2fastq to use. (Default to use specifications from the
            settings)
          bases_mask: if set then use this as an alternative bases mask setting
          ignore_missing_bcl: if True then run bcl2fastq with --ignore-missing-bcl
          ignore_missing_stats: if True then run bcl2fastq with --ignore-missing-stats
          no_lane_splitting: if True then run bcl2fastq with --no-lane-splitting
          minimum_trimmed_read_length: if set then supply to bcl2fastq with
            --minimum-trimmed-read-length N
          mask_short_adapter_reads: if set then supply to bcl2fastq with
            --mask-short-adapter-reads N
          nprocessors: number of processors to run bclToFastq.py with
          runner: (optional) specify a non-default job runner to use for fastq
            generation
        """
        # Directories
        analysis_dir = self.params.analysis_dir
        # Bases mask
        if bases_mask is None:
            bases_mask = self.params.bases_mask
        else:
            self.params['bases_mask'] = bases_mask
        # Check for basic information needed to do bcl2fastq conversion
        if self.params.data_dir is None:
            raise Exception("No source data directory")
        if bases_mask is None:
            raise Exception("No bases mask")
        # Number of cores
        if nprocessors is None:
            nprocessors = self.settings.bcl2fastq.nprocessors
        # Whether to use lane splitting
        if no_lane_splitting is None:
            try:
                no_lane_splitting = self.settings.platform[self.metadata.platform].no_lane_splitting
            except (KeyError,AttributeError):
                pass
            if no_lane_splitting is None:
                no_lane_splitting = self.settings.bcl2fastq.no_lane_splitting
        # Determine which bcl2fastq software to use
        if require_bcl2fastq is None:
            try:
                require_bcl2fastq = self.settings.platform[self.metadata.platform].bcl2fastq
            except (KeyError,AttributeError):
                pass
            if require_bcl2fastq is None:
                require_bcl2fastq = self.settings.bcl2fastq.default_version
        if require_bcl2fastq is not None:
            print "Platform '%s' requires bcl2fastq version %s" \
                % (self.metadata.platform,require_bcl2fastq)
        else:
            logging.warning("No bcl2fastq version explicitly specified")
        bcl2fastq = bcl2fastq_utils.available_bcl2fastq_versions(
            require_bcl2fastq)
        if bcl2fastq:
            bcl2fastq_exe = bcl2fastq[0]
            bcl2fastq_info = bcl2fastq_utils.bcl_to_fastq_info(bcl2fastq_exe)
        else:
            raise Exception("No appropriate bcl2fastq software located")
        # Store info on bcl2fastq package
        self.metadata['bcl2fastq_software'] = bcl2fastq_info
        # Generate temporary sample sheet with required format
        fmt = bcl2fastq_utils.get_required_samplesheet_format(bcl2fastq_info[2])
        tmp_sample_sheet = os.path.join(self.tmp_dir,
                                        "SampleSheet.%s.%s.csv" %
                                        (fmt,
                                         time.strftime("%Y%m%d%H%M%S")))
        print "Generating '%s' format sample sheet: %s" % (fmt,tmp_sample_sheet)
        bcl2fastq_utils.make_custom_sample_sheet(sample_sheet,
                                                 tmp_sample_sheet,
                                                 fmt=fmt)
        # Put a copy in the log directory
        shutil.copy(tmp_sample_sheet,self.log_dir)
        # Create bcl2fastq directory
        bcl2fastq_dir = self.add_directory(unaligned_dir)
        # Determine initial number of mismatches
        nmismatches = bcl2fastq_utils.get_nmismatches(bases_mask)
        # Check for barcode collisions
        collisions = bcl2fastq_utils.check_barcode_collisions(tmp_sample_sheet,
                                                              nmismatches)
        if collisions:
            # Report problem barcodes
            logging.warning("Barcode collisions detected using %d mismatches"
                            % nmismatches)
            for collision in collisions:
                logging.warning("Barcode collision for barcodes: %s, %s"
                                % (collision[0],collision[1]))
            # Reduce mismatches to try and address the collisions
            logging.warning("Attempting to address by adjusting #mismatches")
            while nmismatches > 0 and collisions:
                nmismatches -= 1
                collisions = bcl2fastq_utils.check_barcode_collisions(
                    tmp_sample_sheet,
                    nmismatches)
                if not collisions:
                    print "No collisions using %d mismatches" % nmismatches
                    break
            else:
                # Unable to address collisions, bail out
                raise Exception("Barcode collisions with zero mismatches "
                                "(duplicated indexes?): unable to proceed")
        else:
            print "No barcode collisions detected using %d mismatches" % \
                nmismatches
        # Report values and settings
        print "Bcl-to-fastq exe      : %s" % bcl2fastq_exe
        print "Bcl-to-fastq version  : %s %s" % (bcl2fastq_info[1],
                                               bcl2fastq_info[2])
        print "Sample sheet          : %s" % os.path.basename(tmp_sample_sheet)
        print "Bases mask            : %s" % bases_mask
        print "Nmismatches           : %d" % nmismatches
        print "Nprocessors           : %s" % nprocessors
        print "Ignore missing bcl    : %s" % ignore_missing_bcl
        print "Ignore missing stats  : %s" % ignore_missing_stats
        print "No lane splitting     : %s" % no_lane_splitting
        print "Min trimmed read len  : %s" % minimum_trimmed_read_length
        print "Mask short adptr reads: %s" % mask_short_adapter_reads
        # Set up runner
        if runner is not None:
            runner = fetch_runner(runner)
        else:
            runner = self.settings.runners.bcl2fastq
        runner.set_log_dir(self.log_dir)
        # Run bcl2fastq
        bcl2fastq = applications.Command('bclToFastq.py',
                                         '--nprocessors',nprocessors,
                                         '--use-bases-mask',bases_mask,
                                         '--nmismatches',nmismatches,
                                         '--ignore-missing-control')
        if ignore_missing_bcl:
            bcl2fastq.add_args('--ignore-missing-bcl')
        if ignore_missing_stats:
            bcl2fastq.add_args('--ignore-missing-stats')
        if no_lane_splitting:
            bcl2fastq.add_args('--no-lane-splitting')
        if minimum_trimmed_read_length is not None:
            bcl2fastq.add_args('--minimum-trimmed-read-length',
                               minimum_trimmed_read_length)
        if mask_short_adapter_reads is not None:
            bcl2fastq.add_args('--mask-short-adapter-reads',
                               mask_short_adapter_reads)
        bcl2fastq.add_args('--platform',
                           self.metadata.platform,
                           '--bcl2fastq_path',
                           bcl2fastq_exe,
                           primary_data_dir,
                           bcl2fastq_dir,
                           tmp_sample_sheet)

        print "Running %s" % bcl2fastq
        bcl2fastq_job = simple_scheduler.SchedulerJob(runner,
                                                      bcl2fastq.command_line,
                                                      name='bclToFastq',
                                                      working_dir=os.getcwd())
        bcl2fastq_job.start()
        try:
            bcl2fastq_job.wait()
        except KeyboardInterrupt,ex:
            logging.warning("Keyboard interrupt, terminating bcl2fastq")
            bcl2fastq_job.terminate()
            raise ex
        exit_code = bcl2fastq_job.exit_code
        print "bcl2fastq completed: exit code %s" % exit_code
        if exit_code != 0:
            logging.error("bcl2fastq exited with an error")
        return exit_code

    def generate_stats(self,stats_file=None,per_lane_stats_file=None,
                       unaligned_dir=None,add_data=False,nprocessors=None,
                       runner=None):
        """Generate statistics for FASTQ files

        Generates statistics for all FASTQ files found in the
        'unaligned' directory, by running the 'fastq_statistics.py'
        program.

        Arguments
          stats_file: (optional) specify the name and path of
            a non-default file to write the statistics to
            (defaults to 'statistics.info' unless over-ridden by
            local settings)
          per_lane_stats_file: (optional) path for per-lane
            statistics output file (defaults to
            'per_lane_statistics.info' unless over-ridden by
            local settings)
          unaligned_dir: (optional) where to look for Fastq files
            from bcl2fastq
          add_data: (optional) if True then add stats to the existing
            stats files (default is to overwrite existing stats
            files)
          nprocessors: (optional) number of cores to use when
            running 'fastq_statistics.py'
          runner: (optional) specify a non-default job runner to
            use for fastq_statistics.py

        """
        # Get file names for output files
        if stats_file is None:
            if self.params['stats_file'] is not None:
                stats_file = self.params['stats_file']
            else:
                stats_file='statistics.info'
        if per_lane_stats_file is None:
            if self.params['per_lane_stats_file'] is not None:
                per_lane_stats_file = self.params['per_lane_stats_file']
            else:
                per_lane_stats_file='per_lane_statistics.info'
        # Sort out unaligned_dir
        if unaligned_dir is None:
            if self.params.unaligned_dir is None:
                self.params['unaligned_dir'] = 'bcl2fastq'
            unaligned_dir = self.params.unaligned_dir
        if not os.path.exists(os.path.join(self.params.analysis_dir,unaligned_dir)):
            logging.error("Unaligned dir '%s' not found" % unaligned_dir)
        # Set up runner
        if runner is not None:
            runner = fetch_runner(runner)
        else:
            runner = self.settings.runners.stats
        runner.set_log_dir(self.log_dir)
        # Number of cores
        if nprocessors is None:
            nprocessors = self.settings.fastq_stats.nprocessors
        # Generate statistics
        fastq_statistics = applications.Command('fastq_statistics.py',
                                                '--unaligned',unaligned_dir,
                                                '--output',
                                                os.path.join(self.params.analysis_dir,
                                                             stats_file),
                                                '--per-lane-stats',
                                                os.path.join(self.params.analysis_dir,
                                                             per_lane_stats_file),
                                                self.params.analysis_dir,
                                                '--nprocessors',nprocessors)
        if add_data:
            fastq_statistics.add_args('--update')
        print "Generating statistics: running %s" % fastq_statistics
        fastq_statistics_job = simple_scheduler.SchedulerJob(runner,
                                                             fastq_statistics.command_line,
                                                             name='fastq_statistics',
                                                             working_dir=self.analysis_dir)
        fastq_statistics_job.start()
        try:
            fastq_statistics_job.wait()
        except KeyboardInterrupt,ex:
            logging.warning("Keyboard interrupt, terminating fastq_statistics")
            fastq_statistics_job.terminate()
            raise ex
        exit_code = fastq_statistics_job.exit_code
        print "fastq_statistics completed: exit code %s" % exit_code
        if exit_code != 0:
            raise Exception("fastq_statistics exited with an error")
        self.params['stats_file'] = stats_file
        self.params['per_lane_stats_file'] = per_lane_stats_file
        print "Statistics generation completed: %s" % self.params.stats_file
        print "Generating processing QC report"
        processing_qc_html = os.path.join(self.analysis_dir,
                                          "processing_qc.html")
        report_processing_qc(self,processing_qc_html)
        print "Finished"

    def analyse_barcodes(self,unaligned_dir=None,lanes=None,
                         mismatches=None,cutoff=None,
                         barcode_analysis_dir=None,
                         sample_sheet=None,runner=None,
                         force=False):
        """Analyse the barcode sequences for FASTQs for each specified lane

        Run 'analyse_barcodes.py' for one or more lanes, to analyse the
        barcode index sequences in each lane.

        Arguments:
          unaligned_dir: if set then use this as the output directory for
            bcl-to-fastq conversion. Default is 'bcl2fastq' (unless an
            alternative is already specified in the settings)
          lanes: a list of lane numbers (integers) to perform the analysis
            for. Default is to analyse all lanes.
          mismatches: optional, maximum number of mismatches to consider
            when grouping similar barcodes; default is to determine it
            automatically from the bases mask
          cutoff: optional, exclude barcodes with a smaller fraction of
            associated reads than specified cutoff from reporting (e.g.
            '0.001' excludes barcodes with < 0.1% of reads); default is to
            include all barcodes
          sample_sheet: optional, explicitly specify a sample sheet to
            check barcode sequences against (by default will use the
            sample sheet defined in the parameter file for the run)
          barcode_analysis_dir: optional, explicitly specify the
            subdirectory to use for barcode analysis. Counts will be
            written to and read from the 'counts' subdirectory of this
            directory (defaults to 'barcode_analysis')
          runner: set a non-default job runner
          force: if True then forces regeneration of any existing counts
            (default is to reuse existing counts).
        
        """
        # Sort out parameters
        if unaligned_dir is not None:
            self.params['unaligned_dir'] = unaligned_dir
        elif self.params['unaligned_dir'] is None:
            self.params['unaligned_dir'] = 'bcl2fastq'
        # Load data
        illumina_data = self.load_illumina_data(unaligned_dir=unaligned_dir)
        # Handle barcode analysis subdirectories
        if barcode_analysis_dir is not None:
            # Create a subdirectory for barcode analysis
            self.params['barcode_analysis_dir'] = barcode_analysis_dir
        elif self.params['barcode_analysis_dir'] is None:
            self.params['barcode_analysis_dir'] = 'barcode_analysis'
        barcode_analysis_dir = self.params['barcode_analysis_dir']
        # Create barcode and count file subdirectories
        barcode_dir = self.add_directory(barcode_analysis_dir)
        counts_dir = os.path.join(barcode_analysis_dir,'counts')
        if os.path.exists(counts_dir) and force:
            print "Removing existing counts data"
            shutil.rmtree(counts_dir)
        self.add_directory(counts_dir)
        # Map fastq files to counts files
        counts_files = {}
        for project in illumina_data.projects:
            for sample in project.samples:
                for fq in sample.fastq_subset(read_number=1,
                                              full_path=True):
                    counts_files[fq] = os.path.join(
                        counts_dir,
                        "%s.%s.counts" % (project.name,
                                          os.path.basename(fq)))
        if illumina_data.undetermined is not None:
            for sample in illumina_data.undetermined.samples:
                for fq in sample.fastq_subset(read_number=1,
                                              full_path=True):
                    counts_files[fq] = os.path.join(
                        counts_dir,
                        "undetermined.%s.counts" % os.path.basename(fq))
        # Subset of fastq files with no corresponding counts
        missing_counts = filter(lambda fq: not os.path.exists(counts_files[fq]),
                                counts_files.keys())
        # Deal with lanes
        lane_numbers = illumina_data.lanes
        if lanes:
            lanes = [int(lane) for lane in lanes]
            print "Requested analysis for lanes: %s" % \
                ', '.join([str(l) for l in lanes])
        else:
            print "No lanes explicitly requested"
        if len(lane_numbers) == 1 and lane_numbers[0] is None:
            lane_numbers = None
            print "No lanes explicitly defined in fastq file names"
        else:
            print "Lanes explicitly defined in file names: %s" % \
                ', '.join([str(l) for l in lane_numbers])
        if lanes is None or not lane_numbers:
            # Need counts from all files
            req_counts = counts_files.keys()
        else:
            # Get subset of files explictly belonging to this lane
            req_counts = filter(lambda fq: utils.AnalysisFastq(fq).lane_number
                                in lanes,
                                counts_files.keys())
        # Check there are files to examine
        if not req_counts:
            logging.warning("No matching files: nothing to do")
            return
        # Log dir
        self.set_log_dir(self.get_log_subdir('analyse_barcodes'))
        # Set up runner
        if runner is not None:
            runner = fetch_runner(runner)
        else:
            runner = self.settings.general.default_runner
        runner.set_log_dir(self.log_dir)
        # Schedule the jobs needed to do counting
        sched = simple_scheduler.SimpleScheduler(
            runner=runner,
            max_concurrent=self.settings.general.max_concurrent_jobs)
        sched.start()
        # Do counting
        print "Getting counts from fastq files"
        group = sched.group("get_barcode_counts")
        for fq in req_counts:
            if fq in missing_counts:
                # Get counts for this file
                barcode_count_cmd = applications.Command(
                    'analyse_barcodes.py',
                    '-o',os.path.join(self.analysis_dir,counts_files[fq]),
                    '--no-report',fq)
                print "Running %s" % barcode_count_cmd
                group.add(barcode_count_cmd,
                          name='analyse_barcodes.count.%s.%s' %
                          (os.path.basename(counts_files[fq]).split('.')[0],
                           os.path.basename(counts_files[fq]).split('.')[1]))
        group.close()
        # Do reporting
        report_file = os.path.join(barcode_dir,'barcodes.report')
        xls_file = os.path.join(barcode_dir,'barcodes.xls')
        html_file = os.path.join(barcode_dir,'barcodes.html')
        for filen in (report_file,xls_file,html_file):
            if os.path.exists(filen):
                print "Removing existing file: %s" % filen
                os.remove(filen)
        barcode_report_cmd = applications.Command(
            'analyse_barcodes.py',
            '--report',report_file,
            '--xls',xls_file,
            '--html',html_file)
        # Sample sheet
        if sample_sheet is None:
            sample_sheet = self.params.sample_sheet
        barcode_report_cmd.add_args('--sample-sheet',sample_sheet)
        # Implicitly set per-lane analysis if none were explicitly
        # requested but some are defined in sample sheet
        if lanes is None:
            try:
                lanes = sorted(
                    set([str(line['Lane'])
                         for line in IlluminaData.SampleSheet(sample_sheet).data]))
            except KeyError:
                pass
        if lanes:
            barcode_report_cmd.add_args('--lanes',
                                        ','.join([str(l) for l in lanes]))
        # Cutoff
        if cutoff is not None:
            barcode_report_cmd.add_args('--cutoff',cutoff)
        # Mismatches
        if mismatches is None:
            mismatches = bcl2fastq_utils.get_nmismatches(
                self.params.bases_mask)
        barcode_report_cmd.add_args('--mismatches',mismatches)
        # Add the list of count files to process
        barcode_report_cmd.add_args('-c')
        for counts_file in [counts_files[f] for f in req_counts]:
            barcode_report_cmd.add_args(os.path.join(self.analysis_dir,
                                                     counts_file))
        # Write a script file
        script_file = os.path.join(self.log_dir,'report_barcodes.sh')
        utils.write_script_file(script_file,barcode_report_cmd,
                                shell='/bin/sh')
        # Submit command
        print "Running %s" % script_file
        sched.submit(applications.Command('sh',script_file),
                     name='report_barcodes',
                     wait_for=('get_barcode_counts',))
        # Wait for the scheduler to run all jobs
        sched.wait()
        sched.stop()
        # Finish
        if os.path.exists(report_file):
            print "Report written to %s" % report_file
        else:
            logging.error("Missing barcode analysis report: %s" %
                          report_file)
        if os.path.exists(xls_file):
            print "XLS written to %s" % xls_file
        else:
            logging.error("Missing barcode analysis XLS report: %s" %
                          xls_file)
        if os.path.exists(html_file):
            print "HTML written to %s" % html_file
        else:
            logging.error("Missing barcode analysis HTML report: %s" %
                          html_file)

    def remove_primary_data(self):
        """Remove primary data

        """
        primary_data = os.path.join(self.params.primary_data_dir,
                                    os.path.basename(self.params.data_dir))
        if os.path.isdir(primary_data):
            print "Removing copy of primary data in %s" % primary_data
            shutil.rmtree(primary_data)

    def verify_fastq_generation(self,unaligned_dir=None,lanes=None,
                                include_sample_dir=False):
        """Check that generated Fastqs match sample sheet predictions

        Arguments:
          unaligned_dir (str): explicitly specify the bcl2fastq output
            directory to check
          lanes (list): specify a list of lane numbers (integers) to
            check (others will be ignored)
          include_sample_dir (bool): if True then include a
            'sample_name' directory level when checking for
            bcl2fastq2 outputs, even if one shouldn't be present

        Returns:
          True if outputs match sample sheet, False otherwise.
 
        """
        if unaligned_dir is None:
            if self.params.unaligned_dir is not None:
                unaligned_dir = self.params.unaligned_dir
            else:
                logging.debug("Bcl2fastq output directory not defined")
                return False
        else:
            logging.warning("Checking custom bcl2fastq output directory '%s'" %
                            unaligned_dir)
        bcl_to_fastq_dir = os.path.join(self.analysis_dir,unaligned_dir)
        if not os.path.isdir(bcl_to_fastq_dir):
            # Directory doesn't exist
            return False
        # Make a temporary sample sheet to verify against
        tmp_sample_sheet = os.path.join(self.tmp_dir,
                                        "SampleSheet.verify.%s.csv" %
                                        time.strftime("%Y%m%d%H%M%S"))
        bcl2fastq_utils.make_custom_sample_sheet(self.params.sample_sheet,
                                                 tmp_sample_sheet,
                                                 lanes=lanes)
        # Try to create an IlluminaData object
        try:
            illumina_data = IlluminaData.IlluminaData(
                self.analysis_dir,
                unaligned_dir=unaligned_dir)
        except IlluminaData.IlluminaDataError as ex:
            # Failed to initialise
            logging.warning("Failed to get information from %s: %s" %
                            (bcl_to_fastq_dir,ex))
            return False
        # Do check
        return IlluminaData.verify_run_against_sample_sheet(
            illumina_data,
            tmp_sample_sheet,
            include_sample_dir=include_sample_dir)

    def merge_fastq_dirs(self,primary_unaligned_dir,output_dir=None,
                         dry_run=False):
        """
        Combine multiple 'unaligned' output directories into one

        This method combines the output from multiple runs of
        CASAVA/bcl2fastq into a single 'unaligned'-equivalent
        directory.

        Currently it operates in an automatic mode and should
        detect additional 'unaligned' dirs on its own.

        Arguments:
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
        print "Fastqs will be merged into '%s'" % output_dir
        # Collect unaligned dirs
        print "Collecting bcl2fastq directories"
        primary_illumina_data = None
        unaligned_dirs = {}
        for dirn in bcf_utils.list_dirs(self.analysis_dir):
            try:
                illumina_data = IlluminaData.IlluminaData(self.analysis_dir,
                                                          unaligned_dir=dirn)
                if dirn == primary_unaligned_dir:
                    print "* %s (primary dir)" % dirn
                    primary_illumina_data = illumina_data
                elif dirn.endswith(".bak") or dirn.startswith("save."):
                    print "Ignoring %s" % dirn
                else:
                    print "* %s" % dirn
                    unaligned_dirs[dirn] = illumina_data
            except Exception as ex:
                logging.debug("Rejecting %s: %s" % (dirn,ex))
        # Check primary unaligned dir
        if primary_illumina_data is None:
            raise Exception("Primary dir '%s' doesn't exist, or doesn't "
                            "contain data?" % primary_unaligned_dir)
        # Is there anything to do?
        if not unaligned_dirs:
            print "No extra bcl2fastq output directories found, nothing to do"
            return 0
        # Make log directory and set up scheduler (if not dry run)
        if not dry_run:
            self.set_log_dir(self.get_log_subdir('merge_fastq_dirs'))
            runner = self.settings.general.default_runner
            runner.set_log_dir(self.log_dir)
            sched = simple_scheduler.SimpleScheduler(
                runner=runner,
                max_concurrent=self.settings.general.max_concurrent_jobs)
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
        print "Checking primary data directory"
        fmt = primary_illumina_data.format
        paired_end = primary_illumina_data.paired_end
        no_lane_splitting = (len(primary_illumina_data.lanes) == 1) \
                            and (primary_illumina_data.lanes[0] is None)
        print "* Format: %s" % fmt
        print "* no-lane-splitting: %s" % ('yes' if no_lane_splitting
                                           else 'no')
        print "* paired-end: %s" % ('yes' if paired_end else 'no')
        print "* undetermined dir: %s" % undetermined_dir
        consistent_data = True
        for unaligned_dir in unaligned_dirs:
            illumina_data = unaligned_dirs[unaligned_dir]
            fmt0 = illumina_data.format
            no_lane_splitting0 = (len(illumina_data.lanes) == 1) \
                                 and (primary_illumina_data.lanes[0] is None)
            if (fmt0 != fmt) or (no_lane_splitting0 != no_lane_splitting):
                print "!!! %s: inconsistent format to primary data dir !!!" \
                    % unaligned_dir
                consistent_data = False
        if not consistent_data:
            raise Exception("Data directories not consistent with primary "
                            "dir '%s'" % primary_unaligned_dir)
        # Collect the projects from the extra directories
        projects = []
        undetermined = []
        for unaligned_dir in unaligned_dirs:
            print "Examining projects in %s:" % unaligned_dir
            illumina_data = unaligned_dirs[unaligned_dir]
            for project in illumina_data.projects:
                if not filter(lambda p: p.name == project.name,
                              projects):
                    print "- %s: will be merged in" % project.name
                    projects.append(project)
                else:
                    raise Exception("collision: %s already exists" %
                                    project.name)
            # Deal with undetermined reads
            if illumina_data.undetermined is not None:
                print "Examining undetermined samples:"
                if no_lane_splitting:
                    # No lane info: should merge undetermined fastqs
                    for sample in illumina_data.undetermined.samples:
                        print "- %s: reads will be concatenated" % sample.name
                        undetermined.append(sample)
                else:
                    for sample in illumina_data.undetermined.samples:
                        if not filter(lambda s: s.name == sample.name,
                                      undetermined):
                            print "- %s: will be merged in" % sample.name
                            undetermined.append(sample)
                        else:
                            raise Exception("collision: %s already exists" %
                                            sample.name)
            else:
                print "No undetermined samples"
        # Collect any remaining projects from the primary
        # unaligned directory
        print "Examining projects in primary dir %s:" \
            % primary_unaligned_dir
        for project in primary_illumina_data.projects:
            if not filter(lambda p: p.name == project.name,
                          projects):
                print "- %s: will be merged in" % project.name
                projects.append(project)
            else:
                print "- %s: already exists, will be discarded" \
                    % project.name
        # Sort out the undetermined reads
        print "Examining undetermined samples:"
        if no_lane_splitting:
            # No lane info: should merge undetermined fastqs
            for sample in primary_illumina_data.undetermined.samples:
                print "- %s: reads will be concatenated" % sample.name
                undetermined.insert(0,sample)
        else:
            for sample in primary_illumina_data.undetermined.samples:
                if not filter(lambda s: s.name == sample.name,
                              undetermined):
                    print "- %s: will be merged in" % sample.name
                    undetermined.insert(0,sample)
                else:
                    print "- %s: already exists, will be discarded" \
                        % sample.name
        # Make a new directory for the merging
        merge_dir = os.path.join(self.analysis_dir,
                                 output_dir + ".new")
        if undetermined_dir is not None:
            merge_undetermined_dir = os.path.join(merge_dir,
                                                  undetermined_dir)
        else:
            merge_undetermined_dir = merge_dir
        if not dry_run:
            print "Making temporary merge directory %s" % merge_dir
            bcf_utils.mkdir(merge_dir)
            if not os.path.exists(merge_undetermined_dir):
                print "Making directory for undetermined %s" \
                    %  merge_undetermined_dir
                bcf_utils.mkdir(merge_undetermined_dir)
        # Copy the projects
        print "Importing projects:"
        for project in projects:
            print "- %s" % project.name
            project_dir = os.path.join(merge_dir,
                                       os.path.basename(project.dirn))
            cmd = fileops.copytree_command(project.dirn,project_dir)
            print "- Running %s" % cmd
            if not dry_run:
                job = sched.submit(cmd,
                                   name="copy_project.%s" % project.name,
                                   wd=merge_dir)
                print "Job: %s" % job
                jobs.append(job)
        # Handle the undetermined reads
        print "Dealing with undetermined reads:"
        if no_lane_splitting:
            # No lane info: merge undetermined fastqs
            for read in (1,2):
                if read == 2 and not paired_end:
                    break
                cmd = applications.Command('concat_fastqs.py')
                for sample in undetermined:
                    fastqs = sample.fastq_subset(read_number=read,
                                                 full_path=True)
                    cmd.add_args(*fastqs)
                cmd.add_args(os.path.join(
                    merge_undetermined_dir,
                    "Undetermined_S0_R%s_001.fastq.gz" % read))
                print "- Running %s" % cmd
                if not dry_run:
                    job = sched.submit(cmd,
                                       name="merge_undetermined.R%s" % read,
                                       wd=merge_dir)
                    print "Job: %s" % job
                    jobs.append(job)
        else:
            for sample in undetermined:
                print "- %s" % sample.name
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
                    cmd = fileops.copytree_command(sample.dirn,sample_dir)
                    print "- Running %s" % cmd
                    if not dry_run:
                        job = sched.submit(cmd,
                                           name="copy_sample_dir.%s" %
                                           sample.name,
                                           wd=merge_dir)
                        print "Job: %s" % job.name
                        jobs.append(job)
        # Make expected subdirs for bcl2fastq2
        if not dry_run and fmt == "bcl2fastq2":
            for dirn in ('Reports','Stats'):
               bcf_utils.mkdir(os.path.join(merge_dir,dirn))
        # Wait for scheduler jobs to complete
        if not dry_run:
            sched.wait()
            sched.stop()
            # Check job exit status
            exit_status = 0
            for j in jobs:
                exit_status += j.exit_status
                if j.exit_status != 0:
                    logging.warning("Job failed: %s" % j)
            if exit_status:
                logging.critical("One or more jobs failed (non-zero "
                                 "exit status)")
                return exit_status
        # Move all the 'old' directories out of the way
        all_unaligned = [u for u in unaligned_dirs]
        all_unaligned.append(primary_unaligned_dir)
        for unaligned_dir in all_unaligned:
            unaligned_backup = os.path.join(self.analysis_dir,
                                            "save.%s" %
                                            unaligned_dir)
            print "Moving %s to %s" % (unaligned_dir,
                                       unaligned_backup)
            if not dry_run:
                shutil.move(os.path.join(self.analysis_dir,unaligned_dir),
                            unaligned_backup)
        # Rename the merged directory
        print "Renaming %s to %s" % (merge_dir,output_dir)
        if not dry_run:
            shutil.move(merge_dir,
                        os.path.join(self.analysis_dir,output_dir))
        # Reset the bcl2fastq dir
        if not dry_run:
            self.params['unaligned_dir'] = output_dir
        # Make a new 'projects.info' metadata file
        project_metadata_file = os.path.join(self.analysis_dir,
                                             'projects.info')
        if os.path.exists(project_metadata_file):
            print "Moving existing projects.info file out of the way"
            if not dry_run:
                os.rename(project_metadata_file,
                          os.path.join(self.analysis_dir,
                                       'save.projects.info'))
        print "Creating new projects.info file"
        if not dry_run:
            self.make_project_metadata_file()
        return 0

    def log_analysis(self):
        # Add a record of the analysis to the logging file
        raise NotImplementedError

    def import_project(self,project_dir):
        """
        Import a project directory into this analysis directory

        Arguments:
          project_dir (str): path to project directory to be
            imported

        """
        # Check that target directory exists
        project_dir = os.path.abspath(project_dir)
        # Check that project doesn't already exist
        project_name = os.path.basename(project_dir)
        project_metadata = self.load_project_metadata()
        if project_name in [p['Project'] for p in project_metadata] or \
           analysis.AnalysisProject(project_name,
                                    os.path.join(self.analysis_dir,
                                                 project_name)).exists:
            raise Exception("Project called '%s' already exists" %
                            project_name)
        # Load target as a project
        project = analysis.AnalysisProject(project_name,project_dir)
        # Rsync the project directory
        print "Importing project directory contents for '%s'" % project_name
        try:
            excludes = ['--exclude=tmp.*',
                        '--exclude=qc_report.*']
            rsync = applications.general.rsync(project_dir,
                                               self.analysis_dir,
                                               extra_options=excludes)
            print "Running %s" % rsync
            status = rsync.run_subprocess(log=self.log_path('import_project.rsync.log'))
        except Exception as ex:
            logging.error("Exception importing project: %s" % ex)
            raise ex
        if status != 0:
            raise Exception("Failed to import project from %s (status %s)" %
                            (project_dir,status))
        # Update the projects.info metadata file
        print "Updating projects.info file with imported project"
        project_metadata = self.load_project_metadata()
        sample_names = [s.name for s in project.samples]
        project_metadata.add_project(project_name,
                                     sample_names,
                                     user=project.info.user,
                                     library_type=project.info.library_type,
                                     single_cell_platform=project.info.single_cell_platform,
                                     organism=project.info.organism,
                                     PI=project.info.PI,
                                     comments=project.info.comments)
        project_metadata.save()
        # Report
        print "Projects now in metadata file:"
        for p in project_metadata:
            print "- %s" % p['Project']
        # Update the QC report
        try:
            project = self.get_analysis_projects(pattern=project_name)[0]
        except Exception as ex:
            logging.error("Exception when trying to acquire project %s: %s"
                          % (project_name,ex))
            return
        if project.qc is None:
            print "No QC for %s" % project_name
        else:
            if project.qc.verify():
                try:
                    project.qc_report()
                    print "Updated QC report for %s" % project_name
                except Exception, ex:
                    logging.error("import_project: failed to generate QC "
                                  "report for %s" % project_name)

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
