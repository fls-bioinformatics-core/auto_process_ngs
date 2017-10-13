#!/usr/bin/env python
#
#     auto_processor.py: core module for automated processing of Illumina sequence data
#     Copyright (C) University of Manchester 2013-17 Peter Briggs
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
import time
import ast
import gzip
import string
import bcftbx.IlluminaData as IlluminaData
import bcftbx.platforms as platforms
import bcftbx.TabFile as TabFile
import bcftbx.utils as bcf_utils
import bcftbx.htmlpagewriter as htmlpagewriter
from bcftbx.JobRunner import fetch_runner
import config
import applications
import fileops
import utils
import simple_scheduler
import bcl2fastq_utils
import samplesheet_utils
import icell8_utils
import tenx_genomics_utils
import settings
from .qc.processing import report_processing_qc
from .exceptions import MissingParameterFileException
from auto_process_ngs import get_version

#######################################################################
# Classes
#######################################################################

class AutoProcess:
    """
    Class implementing an automatic fastq generation and QC
    processing procedure for Illumina sequencing data

    """
    def __init__(self,analysis_dir=None,allow_save_params=True):
        """
        Create a new AutoProcess instance

        Arguments:
          analysis_dir (str): name/path for existing analysis
            directory
          allow_save_params (bool): if True then allow updates
            to parameters to be saved back to the parameter file
            (this is the default)

        """
        # Initialise
        self._log_dir = 'logs'
        # Load configuration settings
        self.settings = settings.Settings()
        # Create empty parameter and metadata set
        self.params = utils.AnalysisDirParameters()
        self.metadata = utils.AnalysisDirMetadata()
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
        # Make guesses for other metadata
        if self.metadata.platform is None:
            # Attempt to detect sequencing platform
            self.metadata['platform'] = \
                platforms.get_sequencer_platform(self.analysis_dir)
            if self.metadata.platform is None:
                logging.warning("Unable to identify platform from "
                                "directory name")
            else:
                print "Setting 'platform' metadata item to %s" % \
                    self.metadata.platform
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
                datestamp,instrument,run_number = \
                    IlluminaData.split_run_name(self.run_name)
                if self.metadata.instrument_name is None:
                    self.metadata['instrument_name'] = instrument
                if self.metadata.instrument_datestamp is None:
                    self.metadata['instrument_datestamp'] = datestamp
                if self.metadata.instrument_run_number is None:
                    self.metadata['instrument_run_number'] = run_number
            except Exception, ex:
                logging.warning("Unable to extract missing instrument metadata "
                                "from run name")

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
            project_metadata = utils.ProjectMetadataFile(filen)
        else:
            # First try to populate basic metadata from existing projects
            logging.debug("Metadata file not found, guessing basic data")
            project_metadata = utils.ProjectMetadataFile()
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
                test_project = utils.AnalysisProject(
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
            project_metadata = utils.ProjectMetadataFile(filen)
        else:
            # New (empty) metadata file
            print "Creating new project metadata file: %s" % filen
            project_metadata = utils.ProjectMetadataFile()
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

    def get_analysis_projects_from_dirs(self):
        """
        Return a list of AnalysisProjects in the analysis directory

        Tests each of the subdirectories in the top-level of the
        analysis directory and rejects any that appear to be
        CASVAVA/bcl2fastq outputs or which don't successfully load
        as AnalysisProject instances.

        Returns:
          List: list of AnalysisProject instances.

        """
        logging.debug("Testing subdirectories to determine analysis projects")
        projects = []
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
                logging.warning("Exception when attempting to load "
                                "subdir '%s' as CASAVA/bcl2fastq output "
                                "(ignored): %s" % (dirn,ex))
            # Try loading as a project
            test_project = utils.AnalysisProject(
                dirn,os.path.join(self.analysis_dir,dirn))
            if test_project.is_analysis_dir:
                logging.debug("* %s: analysis directory" % dirn)
                projects.append(test_project)
            else:
                logging.debug("* %s: rejected" % dirn)
        return projects

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
            self._log_dir = self.log_path(path)
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
        return utils.get_numbered_subdir(name,
                                         parent_dir=self.log_dir)

    def __del__(self):
        """
        Implement __del__ method

        Peforms clean up operations (e.g. save parameters,
        remove temporary files etc) when the AutoProcess
        object is destroyed.

        """
        if self.analysis_dir is None:
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
        if analysis_dir is None:
            self.analysis_dir = os.path.join(os.getcwd(),
                                             os.path.basename(data_dir))+'_analysis'
        else:
            self.analysis_dir = os.path.abspath(analysis_dir)
        # Create the analysis directory structure
        if not os.path.exists(self.analysis_dir):
            # Create directory structure
            self.create_directory(self.analysis_dir)
            self.log_dir
            self.script_code_dir
        else:
            # Directory already exists
            logging.warning("Analysis directory already exists")
            # check for parameter file
            if self.has_parameter_file:
                self.load_parameters()
            else:
                logging.warning("No parameter file found in %s" %
                                self.analysis_dir)
        # Identify missing data and attempt to acquire
        # Sequencing platform
        platform = self.metadata.platform
        if platform is None:
            print "Identifying platform from data directory name"
            platform = platforms.get_sequencer_platform(data_dir)
        print "Platform identified as '%s'" % platform
        # Run datestamp, instrument name and instrument run number
        try:
            datestamp,instrument,run_number = IlluminaData.split_run_name(
                os.path.basename(self.analysis_dir))
            run_number = run_number.lstrip('0')
        except Exception, ex:
            logging.warning("Unable to extract information from run name '%s'" \
                            % os.path.basename(self.analysis_dir))
            logging.warning("Exception: %s" % ex)
            datestamp = None
            instrument= None
            run_number = None
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
                logging.error("Unable to acquire sample sheet")
                return
            # Keep a copy of the original sample sheet
            original_sample_sheet = os.path.join(self.analysis_dir,
                                                 'SampleSheet.orig.csv')
            print "Copying original sample sheetd to %s" % original_sample_sheet
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
        # Assay type (= kit)
        try:
            assay = IlluminaData.SampleSheet(original_sample_sheet).header['Assay']
        except KeyError:
            logging.warning("No element 'Assay' found in sample sheet")
            assay = None
        # Bases mask
        bases_mask = self.params.bases_mask
        if bases_mask is None:
            print "Acquiring RunInfo.xml to determine bases mask..."
            tmp_run_info = os.path.join(self.tmp_dir,'RunInfo.xml')
            rsync = applications.general.rsync(os.path.join(data_dir,
                                                            'RunInfo.xml'),
                                               self.tmp_dir)
            status = rsync.run_subprocess(log=self.log_path('rsync.run_info.log'))
            bases_mask = bcl2fastq_utils.get_bases_mask(tmp_run_info,
                                                        custom_sample_sheet)
            os.remove(tmp_run_info)
        print "Corrected bases mask: %s" % bases_mask
        # Generate and print predicted outputs and warnings
        print samplesheet_utils.predict_outputs(sample_sheet=sample_sheet)
        samplesheet_utils.check_and_warn(sample_sheet=sample_sheet)
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
        platform = platforms.get_sequencer_platform(analysis_dir)
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
                new_project = utils.AnalysisProject(
                    project.name,
                    os.path.join(clone_dir,project.name),
                    user=project.info.user,
                    PI=project.info.PI,
                    library_type=project.info.library_type,
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
        # Return the analysis projects in a list
        #
        # By default returns all projects within the analysis directory
        # (including 'undetermined').
        #
        # If the 'pattern' is not None then it should be a simple pattern
        # used to match against available names to select a subset of
        # projects (see bcf_utils.name_matches).
        project_metadata = self.load_project_metadata(self.params.project_metadata)
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
                logging.error("Unable to resolve directory for project '%s'" % name)
                logging.error("Possible dirs: %s" % dirs)
                raise Exception("Unable to resolve directory for project '%s'" % name)
            # Attempt to load the project data
            project_dir = os.path.join(self.analysis_dir,project_dir)
            projects.append(utils.AnalysisProject(name,project_dir))
        # Add undetermined reads directory
        if bcf_utils.name_matches('undetermined',pattern):
            undetermined_analysis = self.undetermined()
            if undetermined_analysis is not None and \
               'undetermined' not in [p.name for p in projects]:
                projects.append(undetermined_analysis)
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
        return utils.AnalysisProject(dirs[0],undetermined_dir)
        
    def make_fastqs(self,protocol='standard',
                    unaligned_dir=None,sample_sheet=None,lanes=None,
                    ignore_missing_bcl=False,ignore_missing_stats=False,
                    skip_rsync=False,remove_primary_data=False,
                    nprocessors=None,require_bcl2fastq_version=None,
                    bases_mask=None,no_lane_splitting=None,
                    minimum_trimmed_read_length=None,
                    mask_short_adapter_reads=None,
                    generate_stats=True,stats_file=None,
                    per_lane_stats_file=None,
                    report_barcodes=False,barcodes_file=None,
                    skip_fastq_generation=False,
                    only_fetch_primary_data=False,
                    create_empty_fastqs=None,runner=None,
                    cellranger_jobmode=None,
                    cellranger_mempercore=None,
                    cellranger_maxjobs=None,
                    cellranger_jobinterval=None):
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
          report_barcodes     : if True then analyse barcodes in outputs (default is False
                                i.e. don't do barcode analyses)
          barcodes_file       : if set then use this as the name of the report file for
                                barcode sequences analysis
          skip_fastq_generation: if True then don't perform fastq generation
          only_fetch_primary_data: if True then fetch primary data, don't do anything else
          create_empty_fastqs : if True then create empty 'placeholder' fastq
                                files for any missing fastqs after bcl2fastq
                                (must have completed with zero exit status)
          runner              : (optional) specify a non-default job runner to use for
                                fastq generation
          cellranger_jobmode  : (optional) job mode to run cellranger in
                                (10xGenomics Chromium SC data only)
          cellranger_mempercore: (optional) memory assumed per core (in Gbs)
                                (10xGenomics Chromium SC data only)
          cellranger_maxjobs  : (optional) maxiumum number of concurrent jobs
                                to run (10xGenomics Chromium SC data only)
          cellranger_jobinterval: (optional) how often jobs are submitted (in
                                ms) (10xGenomics Chromium SC data only)

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
        # Adjust verification settings for 10xGenomics Chromium SC
        # data if necessary
        verify_include_sample_dir = False
        if protocol == '10x_chromium_sc':
            if tenx_genomics_utils.has_chromium_sc_indices(sample_sheet):
                # Force inclusion of sample-name subdirectories
                # when verifying Chromium SC data
                print "Sample sheet includes Chromium SC indices"
                verify_include_sample_dir = True
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
        # Do fastq generation using the specified protocol
        if not skip_fastq_generation:
            # Set primary data location and report info
            primary_data_dir = os.path.join(
                self.params.primary_data_dir,
                os.path.basename(self.params.data_dir))
            print "Primary data dir      : %s" % primary_data_dir
            illumina_run = IlluminaData.IlluminaRun(primary_data_dir)
            print "Platform              : %s" % illumina_run.platform
            print "Bcl format            : %s" % illumina_run.bcl_extension
            if protocol == 'icell8':
                # ICell8 data
                # Update bcl2fastq settings appropriately
                print "Updating read trimming and masking for ICell8"
                minimum_trimmed_read_length = 21
                mask_short_adapter_reads = 0
                # Reset the default bases mask
                bases_mask = IlluminaData.IlluminaRunInfo(
                    illumina_run.runinfo_xml).bases_mask
                bases_mask = icell8_utils.get_icell8_bases_mask(bases_mask)
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
                        cellranger_jobmode=cellranger_jobmode,
                        cellranger_maxjobs=cellranger_maxjobs,
                        cellranger_mempercore=cellranger_mempercore,
                        cellranger_jobinterval=cellranger_jobinterval,
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
                    include_sample_dir=include_sample_dir)
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
                    logging.warning("Making 'empty' Fastqs as placeholders")
                    for fq in missing_fastqs:
                        fastq = os.path.join(self.analysis_dir,
                                             self.params.unaligned_dir,fq)
                        print "-- %s" % fastq
                    with gzip.GzipFile(filename=fastq,mode='wb') as fp:
                        fp.write('')
                raise Exception("Fastq generation failed to produce "
                                "expected outputs")
        # Generate statistics
        if generate_stats:
            self.generate_stats(stats_file=stats_file,
                                per_lane_stats_file=per_lane_stats_file,
                                unaligned_dir=self.params.unaligned_dir,
                                nprocessors=nprocessors,
                                runner=runner)
        # Count and report barcode sequences
        if report_barcodes:
            self.report_barcodes(barcodes_file)
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
        bcl2fastq.add_args('--bcl2fastq_path',
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
        processing_qc_html = "processing_qc.html"
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
                    '-o',counts_files[fq],
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
            barcode_report_cmd.add_args(counts_file)
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

    def report_barcodes(self,report_file=None):
        """Count and report barcode sequences in FASTQs for each lane

        Runs the 'count_barcodes.py' program to count unique barcode index
        sequences for all FASTQ files in each lane.

        Arguments:
          report_file: (optional) specify the name and path of
            the report file; otherwise this defaults to a file called
            'index_sequences.report'.

        """
        # Output files
        if report_file is None:
            barcode_report = os.path.join(self.analysis_dir,'index_sequences.report')
        else:
            barcode_report = report_file
        barcode_counts = os.path.join(self.analysis_dir,'index_sequences.counts')
        # Set up runner
        runner = self.settings.runners.stats
        runner.set_log_dir(self.log_dir)
        # Run count_barcodes.py
        count_barcodes = applications.Command('count_barcodes.py',
                                              '-o',barcode_counts,
                                              '-r',barcode_report,
                                                os.path.join(self.analysis_dir,
                                                             self.params.unaligned_dir))
        print "Counting barcode index sequences: running %s" % count_barcodes    
        count_barcodes_job = simple_scheduler.SchedulerJob(runner,
                                                           count_barcodes.command_line,
                                                           name='count_barcodes',
                                                           working_dir=self.analysis_dir)
        count_barcodes_job.start()
        count_barcodes_job.wait()
        print "Barcode counting completed"

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
            sched = simple_scheduler.SimpleScheduler(runner=runner)
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

    def setup_analysis_dirs(self,unaligned_dir=None,
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
            hasn't been set (defaults is to stop if metadata isn't
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
            unaligned_dir = self.params.unaligned_dir
        if unaligned_dir is None:
            logging.error("No unaligned directory, cannot build analysis "
                          "directories")
            raise Exception("Cannot build analysis directories")
        illumina_data = self.load_illumina_data(unaligned_dir=unaligned_dir)
        # Project metadata file
        if project_metadata_file is None:
            project_metadata_file = self.params.project_metadata
        if project_metadata_file is None:
            project_metadata_file = 'projects.info'
        if not os.path.exists(os.path.join(self.params.analysis_dir,
                                           project_metadata_file)):
            logging.warning("No project metadata file '%s' found, "
                            "attempting to create" % project_metadata_file)
            self.make_project_metadata_file(project_metadata_file)
            logging.warning("Update '%s' and rerun" % project_metadata_file)
            return
        project_metadata = self.load_project_metadata(
            project_metadata_file=project_metadata_file,
            check=True)
        # Sanity check that the project data file has been populated
        got_project_data = True
        for line in project_metadata:
            for item in ('User','PI','Organism','Library',):
                if line[item] == '.':
                    logging.warning("Missing data from %s for '%s': %s" %
                                    (self.params.project_metadata,
                                     line['Project'],item))
                    got_project_data = False
        if not got_project_data:
            if ignore_missing_metadata:
                logging.warning("Missing project metadata")
            else:
                logging.error("Missing project metadata")
                raise Exception("Missing project metadata")
        # Create the projects
        n_projects = 0
        for line in project_metadata:
            # Acquire the run name
            run_name = self.run_name
            # Look up project data
            project_name = line['Project']
            user = line['User']
            PI = line['PI']
            organism = line['Organism']
            library_type = line['Library']
            comments = line['Comments']
            # Check it's in the list
            if projects and project_name not in projects:
                logging.warning("Skipping '%s'" % project_name)
                continue
            # Create the project
            project = utils.AnalysisProject(project_name,
                                            os.path.join(self.analysis_dir,
                                                         project_name),
                                            user=user,
                                            PI=PI,
                                            organism=organism,
                                            library_type=library_type,
                                            run=run_name,
                                            comments=comments,
                                            platform=self.metadata.platform)
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
                logging.warning("Failed to create project '%s': %s" %
                                (project_name,ex))
        # Tell us how many were made
        print "Created %d project%s" % (n_projects,'s' if n_projects != 1 else '')
        # Also set up analysis directory for undetermined reads
        if undetermined_project is None:
            undetermined_project = 'undetermined'
        undetermined = illumina_data.undetermined
        if illumina_data.undetermined is not None:
            undetermined = utils.AnalysisProject(undetermined_project,
                                                 os.path.join(
                                                     self.analysis_dir,
                                                     undetermined_project),
                                                 run=run_name,
                                                 comments="Analysis of reads "
                                                 "with undetermined indices",
                                                 platform=self.metadata.platform)
            if not undetermined.exists:
                print "Creating directory '%s' for analysing reads " \
                    "with undetermined indices" % undetermined.name
                undetermined.create_directory(illumina_data.undetermined,
                                              link_to_fastqs=link_to_fastqs)
            else:
                logging.warning("'%s' directory already exists, skipping" %
                                undetermined.name)

    def run_qc(self,projects=None,max_jobs=4,ungzip_fastqs=False,
               fastq_screen_subset=100000,nthreads=1,
               runner=None,fastq_dir=None,qc_dir=None,
               report_html=None,run_multiqc=True):
        """Run QC pipeline script for projects

        Run the illumina_qc.sh script to perform QC on projects.

        Note that if all QC outputs already exist for a project then
        the QC will *not* be run for that project.

        A subset of projects can be selected for QC by setting the
        'projects' argument to a name or pattern, only matching
        projects will be examined.

        Arguments:
          projects: specify a pattern to match one or more projects to
                    run the QC for (default is to run QC for all
                    projects)
          max_jobs: maximum number of jobs that will be scheduled
                    to run at one time (passed to the scheduler;
                    default is 4, set to zero to remove the limit)
          ungzip_fastqs: if True then run the QC script with the
                    '--ungzip-fastqs' option to create decompressed
                    copies of any fastq.gz inputs (default is False,
                    don't decompress the input files)
          fastq_screen_subset: subset of reads to use in fastq_screen
                    (default is 100000, set to zero or None to use
                    all reads)
          nthreads: (optional) specify number of threads to run the
                    QC pipeline with (default is 1)
          runner:   (optional) specify a non-default job runner to
                    use for the QC.
          fastq_dir: (optional) specify the subdirectory to take the
                    Fastq files from; will be used for all projects
                    that are processed (default is 'fastqs')
          qc_dir:   (optional) specify a non-standard directory to
                    write the QC outputs to; will be used for all
                    projects that are processed (default is 'qc')
          report_html: (optional) specify the name for the output
                    HTML QC report (default is '<QC_DIR>_report.html')
          run_multiqc: if True then run MultiQC at the end of the
                    QC run (default)

        Returns:
          UNIX-style integer returncode: 0 = successful termination,
          non-zero indicates an error occurred.

        """
        # Check QC script version
        compatible_versions = ('1.3.0','1.3.1')
        print "Getting QC script information"
        status,qc_script_info = applications.Command('illumina_qc.sh',
                                                     '--version').subprocess_check_output()
        print "Using QC script %s" % qc_script_info.strip()
        version = qc_script_info.strip().split()[-1]
        if version not in compatible_versions:
            logging.error("QC script version is %s, needs %s" %
                          (version,'/'.join(compatible_versions)))
            return 1
        # Process project pattern matching
        if projects is None:
            project_pattern = '*'
            sample_pattern = '*'
        else:
            project_pattern = projects.split('/')[0]
            try:
                sample_pattern = projects.split('/')[1]
            except IndexError:
                sample_pattern = '*'
        # Get project dir data
        projects = self.get_analysis_projects(project_pattern)
        # Check we have projects
        if len(projects) == 0:
            logging.warning("No projects found for QC analysis")
            return 1
        # Set up runner
        if runner is not None:
            qc_runner = fetch_runner(runner)
        else:
            qc_runner = self.settings.runners.qc
        # Set up a simple scheduler
        sched = simple_scheduler.SimpleScheduler(runner=qc_runner,
                                                 max_concurrent=max_jobs)
        sched.start()
        # Look for samples with no/invalid QC outputs and populate
        # pipeline with the associated fastq.gz files
        for project in projects:
            print "*** Setting up QC for %s ***" % project.name
            # Set up qc directory
            project_qc_dir = project.setup_qc_dir(qc_dir,fastq_dir=fastq_dir)
            project.use_qc_dir(project_qc_dir)
            print "Using QC directory %s" % project.qc_dir
            # Sort out fastq source
            fastq_dir = project.qc_info(project_qc_dir).fastq_dir
            project.use_fastq_dir(fastq_dir)
            print "Set fastq_dir to %s" % project.fastq_dir
            # Set up the logs directory
            log_dir = os.path.join(project_qc_dir,'logs')
            if not os.path.exists(log_dir):
                print "Making QC logs directory: %s" % log_dir
                bcf_utils.mkdir(log_dir,mode=0775)
            # Loop over samples and queue up those where the QC
            # isn't validated
            samples = project.get_samples(sample_pattern)
            if len(samples) == 0:
                logging.warning("No samples found for QC analysis in project '%s'" %
                                project.name)
            groups = []
            for sample in samples:
                group = None
                print "Examining files in sample %s" % sample.name
                for fq in sample.fastq:
                    if utils.AnalysisFastq(fq).is_index_read:
                        # Reject index read Fastqs
                        logging.warning("Ignoring index read: %s" %
                                        os.path.basename(fq))
                        continue
                    if sample.verify_qc(project_qc_dir,fq):
                        logging.debug("\t%s: QC verified" % fq)
                    else:
                        print "\t%s: setting up QC run" % os.path.basename(fq)
                        # Create a group if none exists for this sample
                        if group is None:
                            group = sched.group("%s.%s" % (project.name,sample.name),
                                                log_dir=log_dir)
                        # Create and submit a QC job
                        fastq = os.path.join(project.dirn,'fastqs',fq)
                        label = "illumina_qc.%s.%s" % \
                                (project.name,str(utils.AnalysisFastq(fq)))
                        qc_cmd = applications.Command('illumina_qc.sh',fastq)
                        if ungzip_fastqs:
                            qc_cmd.add_args('--ungzip-fastqs')
                        if fastq_screen_subset is None:
                            fastq_screen_subset = 0
                        qc_cmd.add_args(
                            '--threads',nthreads,
                            '--subset',fastq_screen_subset,
                            '--qc_dir',project_qc_dir)
                        job = group.add(qc_cmd,name=label,wd=project.dirn)
                        print "Job: %s" %  job
                # Indicate no more jobs to add
                if group:
                    group.close()
                    groups.append(group.name)
            # Add MultiQC job (if requested)
            if run_multiqc:
                multiqc_out = "multi%s_report.html" % \
                              os.path.basename(project_qc_dir)
                if (not os.path.exists(multiqc_out)) or groups:
                    multiqc_cmd = applications.Command(
                        'multiqc',
                        '--title','%s/%s' % (self.run_name,project.name),
                        '--filename','./%s' % multiqc_out,
                        '--force',
                        project_qc_dir)
                    print "Running %s" % multiqc_cmd
                    label = "multiqc.%s" % project.name
                    job = sched.submit(multiqc_cmd,
                                       name=label,
                                       wd=project.dirn,
                                       log_dir=log_dir,
                                       wait_for=groups)
                else:
                    print "MultiQC report '%s': already exists" % multiqc_out
        # Wait for the scheduler to run all jobs
        sched.wait()
        sched.stop()
        # Verify the outputs and generate QC reports
        failed_projects = []
        for project in projects:
            if not project.verify_qc(qc_dir=qc_dir):
                failed_projects.append(project)
            else:
                qc_base = os.path.basename(qc_dir)
                if report_html is None:
                    out_file = '%s_report.html' % qc_base
                else:
                    out_file = report_html
                if not os.path.isabs(out_file):
                    out_file = os.path.join(project.dirn,out_file)
                title = "%s/%s" % (self.run_name,
                                   project.name)
                if fastq_dir is not None:
                    title = "%s (%s)" % (title,fastq_dir)
                title = "%s: QC report" % title
                print "QC okay, generating report for %s" % project.name
                project.qc_report(qc_dir=qc_dir,
                                  title=title,
                                  report_html=out_file)
            if run_multiqc:
                multiqc_report = os.path.join(project.dirn,multiqc_out)
                if not os.path.exists(multiqc_report):
                    logging.warning("Missing MultiQC report for %s" %
                                    project.name)
                    failed_projects.append(project)
        # Report failed projects
        if failed_projects:
            logging.error("QC failed for one or more samples in following projects:")
            for project in failed_projects:
                logging.error("- %s" % project.name)
            return 1
        # Finish
        return 0

    def copy_to_archive(self,archive_dir=None,platform=None,year=None,dry_run=False,
                        chmod=None,group=None,include_bcl2fastq=False,
                        read_only_fastqs=True,force=False,runner=None):
        """Copy the analysis directory and contents to an archive area

        Copies the contents of the analysis directory to an archive
        area, which can be on a local or remote system.

        The archive directory is constructed in the form

        <TOP_DIR>/<YEAR>/<PLATFORM>/<DIR>/...

        The YEAR and PLATFORM can be overriden using the appropriate
        arguments.

        By default the 'bcl2fastq' directory is omitted from the
        archive, unless the fastq files in any projects are links to
        the data. Inclusion of this directory can be forced by
        setting the appropriate argument.

        The fastqs will be switched to be read-only in the archive
        by default.

        'rsync' is used to perform the transfer.

        Arguments:
          archive_dir: top level archive directory, of the form
            '[[user@]host:]dir'; if not set then use the value from
            the settings.ini file.
          platform: set the value of the <PLATFORM> level in the
            archive.
          year: set the value of the <YEAR> level in the archive;
            if not set then defaults to the current year (4 digits)
          dry_run: report what would be done but don't perform any
            operations.
          chmod: change the mode of the destination files and
            directories according to the supplied argument (e.g.
            'g+w'); if not set then use the value from
            the settings.ini file.
          group: set the group of the destination files to the
            supplied argument; if not set then use the value from
            the settings.ini file.
          include_bcl2fastq: if True then force inclusion of the
            'bcl2fastq' subdirectory; otherwise only include it if
            fastq files in project subdirectories are symlinks.
          read_only_fastqs: if True then make the fastqs read-only
            in the destination directory; otherwise keep the original
            permissions.
          force: if True then do archiving even if key metadata items
            are not set; otherwise abort archiving operation.
          runner: (optional) specify a non-default job runner to use
            for primary data rsync

        Returns:
          UNIX-style integer returncode: 0 = successful termination,
          non-zero indicates an error occurred.

        """
        # Return value
        retval = 0
        # Check first: are there any projects?
        projects = self.get_analysis_projects()
        if not projects:
            raise Exception("No project directories found, nothing to archive")
        # Check metadata
        check_metadata = self.check_metadata(('source','run_number'))
        if not check_metadata:
            if not force:
                logging.error("Some metadata items not set, stopping")
                return
            logging.warning("Some metadata items not set, proceeding")
        # Fetch archive location
        if archive_dir is None:
            archive_dir = self.settings.archive.dirn
        if archive_dir is None:
            raise Exception, "No archive directory specified (use --archive_dir option?)"
        # Construct subdirectory structure i.e. platform and year
        if platform is None:
            platform = self.metadata.platform
        if platform is None:
            raise Exception, "No platform specified (use --platform option?)"
        if year is None:
            year = time.strftime("%Y")
        archive_dir = os.path.join(archive_dir,year,platform)
        print "Copying to archive directory: %s" % archive_dir
        print "Platform: %s" % platform
        print "Year    : %s" % year
        # Determine which directories to exclude
        excludes = ['--exclude=primary_data',
                    '--exclude=save.*',
                    '--exclude=*.bak',
                    '--exclude=tmp.*']
        if not include_bcl2fastq:
            # Determine whether bcl2fastq dir should be included implicitly
            # because there are links from the analysis directories
            for project in projects:
                if project.fastqs_are_symlinks:
                    print "Found at least one project with fastq symlinks (%s)" % project.name
                    include_bcl2fastq = True
                    break
        if not include_bcl2fastq:
            print "Excluding '%s' directory from archive" % self.params.unaligned_dir
            excludes.append('--exclude=%s' % self.params.unaligned_dir)
        # 10xgenomics products to exclude
        excludes.append('--exclude=*.mro')
        excludes.append('--exclude="%s*"' %
                        tenx_genomics_utils.flow_cell_id(self.run_name))
        # Log dir
        self.set_log_dir(self.get_log_subdir('archive'))
        # Set up runner
        if runner is not None:
            runner = fetch_runner(runner)
        else:
            runner = self.settings.runners.rsync
        runner.set_log_dir(self.log_dir)
        # Setup a scheduler for multiple rsync jobs
        sched = simple_scheduler.SimpleScheduler(
            runner=runner,
            max_concurrent=self.settings.general.max_concurrent_jobs)
        sched.start()
        # If making fastqs read-only then transfer them separately
        if read_only_fastqs:
            rsync_fastqs = applications.general.rsync(self.analysis_dir,
                                                      archive_dir,
                                                      prune_empty_dirs=True,
                                                      dry_run=dry_run,
                                                      chmod='ugo-w',
                                                      extra_options=(
                                                          '--include=*/',
                                                          '--include=fastqs/**',
                                                          '--exclude=*',))
            print "Running %s" % rsync_fastqs
            rsync_fastqs_job = sched.submit(rsync_fastqs,
                                            name="rsync.archive_fastqs")
            # Exclude fastqs from main rsync
            excludes.append('--exclude=fastqs')
            wait_for = [rsync_fastqs_job.job_name]
        else:
            rsync_fastqs_job = None
            wait_for = None
        # Main rsync command
        rsync = applications.general.rsync(self.analysis_dir,archive_dir,
                                           prune_empty_dirs=True,
                                           dry_run=dry_run,
                                           chmod=chmod,
                                           extra_options=excludes)
        print "Running %s" % rsync
        rsync_job = sched.submit(rsync,name="rsync.archive",
                                 wait_for=wait_for)
        # Wait for scheduler to complete
        sched.wait()
        sched.stop()
        # Check exit status on job(s)
        if rsync_fastqs_job:
            exit_code = rsync_fastqs_job.exit_code
            print "rsync of FASTQs completed: exit code %s" % exit_code
            if exit_code != 0:
                logging.error("Failed to copy FASTQs to archive location "
                              "(non-zero exit code returned)")
        else:
            exit_code = 0
        print "rsync of data completed: exit code %s" % rsync_job.exit_code
        if rsync_job.exit_code != 0:
            logging.error("Failed to copy data to archive location "
                          "(non-zero exit code returned)")
        # Set returncode
        retval = exit_code if exit_code != 0 else rsync_job.exit_code
        # Set the group
        if group is not None:
            print "Setting group of archived files to '%s'" % group
            fileops.set_group(
                group,
                os.path.join(archive_dir,
                             os.path.basename(self.analysis_dir)))
        # Finish
        return retval

    def log_analysis(self):
        # Add a record of the analysis to the logging file
        raise NotImplementedError

    def publish_qc(self,projects=None,location=None,ignore_missing_qc=False,
                   regenerate_reports=False,force=False):
        """
        Copy the QC reports to the webserver

        Looks for and copies various QC reports and outputs to
        a 'QC server' directory, and generates an HTML index.

        The reports include:

        - processing QC report
        - 'cellranger mkfastq' QC report
        - barcode analysis report

        Also if the analysis includes project directories then for
        each Fastq set in each project:

        - QC report for standard QC
        - MultiQC report

        Also if a project comprises ICell8 or 10xGenomics Chromium
        data:

        - ICell8 processing report, or
        - 'cellranger count' reports for each sample

        Raises an exception if:

        - 'source' and 'run_number' metadata items are not set
        - a subset of projects don't have associated QC outputs
          (unless 'ignore_missing_qc' is True)

        Arguments:
          projects (str): specify a glob-style pattern to match one
            or more projects to publish the reports for (default is
            to publish all reports)
          location (str): override the target location specified in
            the settings; can be of the form
            '[[user@]server:]directory'
          ignore_missing_qc (bool): if True then skip directories
            with missing QC data or reports (default is to raise
            an exception if projects have missing QC)
          regenerate_reports (bool): if True then try to create
            reports even when they already exist (default is to
            use existing reports)
          force (bool): if True then force QC report (re)generation
            even if QC is unverified (default is to raise an
            exception if projects cannot be verified)
        """
        # Turn off saving of parameters etc
        self._save_params = False
        self._save_metadata = False
        # Check metadata
        check_metadata = self.check_metadata(('source','run_number'))
        if not check_metadata:
            raise Exception("Some metadata items not set, stopping")
        # Process pattern matching
        if projects is None:
            project_pattern = '*'
        else:
            project_pattern = projects
        # Get location to publish qc reports to
        if location is None:
            location = fileops.Location(self.settings.qc_web_server.dirn)
        else:
            location = fileops.Location(location)
        # Check the settings
        if location.is_remote:
            print "Copying QC to remote directory"
            print "user:\t%s" % location.user
            print "host:\t%s" % location.server
            print "dirn:\t%s" % location.path
        else:
            print "Copying QC to local directory"
            print "dirn:\t%s" % location.path
        if location.path is None:
            raise Exception, "No target directory specified"
        dirn = os.path.join(str(location),
                            os.path.basename(self.analysis_dir))
        # Get general data
        analysis_dir = utils.AnalysisDir(self.analysis_dir)
        # Get project data
        projects = analysis_dir.get_projects(project_pattern)
        # Check QC situation for each project
        print "Checking QC status for each project:"
        project_qc = {}
        for project in projects:
            print "* Project '%s':" % project.name
            print "  Verifying QC directories:"
            # Get putative qc dirs and verify
            project_qc[project.name] = []
            for qc_dir in project.qc_dirs:
                # Set the source Fastqs dir
                fastq_dir = project.qc_info(qc_dir).fastq_dir
                project.use_fastq_dir(fastq_dir)
                if project.verify_qc(
                        qc_dir=os.path.join(project.dirn,qc_dir)):
                    status = "ok"
                    project_qc[project.name].append(qc_dir)
                else:
                    status = "failed"
                    if force:
                        project_qc[project.name].append(qc_dir)
                print "-- %s...%s" % (qc_dir,status)
            # Check status of reports for verified dirs
            print "  Checking QC reports:"
            for qc_dir in project_qc[project.name]:
                fastq_dir = project.qc_info(qc_dir).fastq_dir
                project.use_fastq_dir(fastq_dir)
                generate_report = regenerate_reports
                qc_zip = os.path.join(project.dirn,
                                      "%s_report.%s.%s.zip" %
                                      (qc_dir,project.name,
                                       os.path.basename(self.analysis_dir)))
                if os.path.isfile(qc_zip):
                    status = "existing report"
                else:
                    status = "no QC report"
                    generate_report = True
                print "-- %s/%s...%s" % (project.name,
                                         qc_dir,
                                         status)
                # (Re)create report
                if generate_report:
                    try:
                        project.qc_report(qc_dir=qc_dir,
                                          force=force)
                    except Exception as ex:
                        logging.error("publish_qc: failed to generate "
                                      "QC report for %s (%s): %s" %
                                      (project.name,qc_dir,ex))
                        project_qc[project.name].remove(qc_dir)
        # Projects with no QC
        no_qc_projects = filter(lambda p: not project_qc[p.name],
                                projects)
        # Final results
        if no_qc_projects:
            # Failed to generate results for some projects
            err_msg = "No QC reports for projects: %s" % \
                      ', '.join([x.name for x in no_qc_projects])
            if not ignore_missing_qc:
                # Fatal error
                raise Exception(err_msg)
            # Proceed with a warning
            logging.warning(err_msg)
        # Remove the 'bad' projects from the list before proceeding
        for project in no_qc_projects:
            print "Project %s will be skipped" % project.name
            projects.remove(project)
        if not projects:
            logging.warning("No projects with QC results to publish")
        # Collect barcode analysis artefacts
        print "Checking for barcode analysis"
        if self.params.barcode_analysis_dir is not None:
            barcode_analysis_dir = self.params.barcode_analysis_dir
        else:
            barcode_analysis_dir = 'barcode_analysis'
        if not os.path.isabs(barcode_analysis_dir):
            barcode_analysis_dir = os.path.join(self.analysis_dir,
                                                barcode_analysis_dir)
        barcodes_files = []
        if os.path.exists(barcode_analysis_dir):
            print "- Barcode analysis dir: %s" % barcode_analysis_dir
            for filen in ('barcodes.report',
                          'barcodes.xls',
                          'barcodes.html'):
                filen = os.path.join(barcode_analysis_dir,filen)
                if os.path.exists(filen):
                    print "...found %s" % os.path.basename(filen)
                    barcodes_files.append(filen)
        if not barcodes_files:
            print "...no barcode analysis found"
        # Make a directory for the QC reports
        fileops.mkdir(dirn)
        # Start building an index page
        title = "QC reports for %s" % os.path.basename(self.analysis_dir)
        index_page = htmlpagewriter.HTMLPageWriter(title)
        # Add CSS rules
        index_page.addCSSRule("h1 { background-color: #42AEC2;\n"
                              "     color: white;\n"
                              "     padding: 5px 10px; }")
        index_page.addCSSRule("table { margin: 10 10;\n"
                              "        border: solid 1px grey;\n"
                              "        background-color: white; }")
        index_page.addCSSRule("th    { background-color: grey;\n"
                              "        color: white;\n"
                              "        padding: 2px 5px; }")
        index_page.addCSSRule("td    { text-align: left;\n"
                              "        vertical-align: top;\n"
                              "        padding: 2px 5px;\n"
                              "        border-bottom: solid 1px lightgray; }")
        index_page.addCSSRule("td.param { background-color: grey;\n"
                              "           color: white;\n"
                              "           padding: 2px 5px;\n"
                              "           font-weight: bold; }")
        index_page.addCSSRule("p.footer { font-style: italic;\n"
                              "           font-size: 70%; }")
        # Build the page
        index_page.add("<h1>%s</h1>" % title)
        # General info
        index_page.add("<h2>General information</h2>")
        index_page.add("<table>")
        index_page.add("<tr><td class='param'>Run name</td><td>%s</td></tr>" %
                       os.path.basename(self.analysis_dir))
        index_page.add("<tr><td class='param'>Run number</td><td>%s</td></tr>" %
                       self.metadata.run_number)
        index_page.add("<tr><td class='param'>Platform</td><td>%s</td></tr>" %
                       self.metadata.platform)
        index_page.add("<tr><td class='param'>Endedness</td><td>%s</td></tr>" %
                       ('Paired end' if analysis_dir.paired_end else 'Single end'))
        try:
            bcl2fastq_software = ast.literal_eval(
                self.metadata.bcl2fastq_software)
        except ValueError:
            bcl2fastq_software = None
        index_page.add("<tr><td class='param'>Bcl2fastq</td><td>%s</td></tr>" %
                       ('unspecified' if not bcl2fastq_software else
                       "%s %s" % (bcl2fastq_software[1],
                                  bcl2fastq_software[2])))
        index_page.add("<tr><td class='param'>Reference</td><td>%s</td></tr>" %
                       self.run_reference_id)
        index_page.add("</table>")
        # Add link to processing statistics
        processing_qc_html = os.path.join(self.analysis_dir,
                                          "processing_qc.html")
        if os.path.exists(processing_qc_html):
            index_page.add("<h2>Processing Statistics</h2>")
            index_page.add("<a href='%s'>Processing QC report</a>" %
                           os.path.basename(processing_qc_html))
        else:
            processing_qc_html = None
        # Add to link to 10xGenomics cellranger QC summaries
        cellranger_qc_html = filter(lambda f: os.path.isfile(f) and
                                    os.path.basename(f).startswith(
                                        "cellranger_qc_summary")
                                    and f.endswith(".html"),
                                    [os.path.join(self.analysis_dir,f)
                                     for f in os.listdir(self.analysis_dir)])
        if cellranger_qc_html:
            index_page.add("<h2>QC summary: cellranger mkfastq</h2>")
            for qc_html in cellranger_qc_html:
                # Check for lane list at tail of QC summary
                # e.g. cellranger_qc_summary_45.html
                # This might not be present
                lanes = qc_html.split('.')[0].split('_')[-1]
                if all(c in string.digits for c in lanes):
                    lanes = ','.join(lanes)
                else:
                    lanes = None
                index_page.add("<a href='%s'>QC summary for %s</a>" %
                               (qc_html,
                                ("all lanes" if lanes is None
                                 else "lanes %s" % lanes)))
        # Barcode analysis
        if barcodes_files:
            # Create section
            index_page.add("<h2>Barcode analysis</h2>")
            index_page.add("<p>Plain text report: "
                           "<a href='barcodes/barcodes.report'>barcodes.report</a> "
                           " | XLS: "
                           "<a href='barcodes/barcodes.xls'>barcodes.xls</a>"
                           " | HTML: "
                           "<a href='barcodes/barcodes.html'>barcodes.html</a>"
                           "</p>")
            # Create subdir and copy files
            barcodes_dirn = os.path.join(dirn,'barcodes')
            fileops.mkdir(barcodes_dirn)
            for filen in barcodes_files:
                fileops.copy(filen,barcodes_dirn)
        if projects:
            # Table of projects
            index_page.add("<h2>QC Reports</h2>")
            index_page.add("<table>")
            index_page.add("<tr><th>Project</th><th>User</th><th>Library</th><th>Organism</th><th>PI</th><th>Samples</th><th>#Samples</th><th>Reports</th></tr>")
            # Set the string to represent "null" table entries
            null_str = '&nbsp;'
            # Deal with QC for each project
            for project in projects:
                # Reset source Fastqs dir
                project.use_fastq_dir()
                # Get local versions of project information
                info = project.info
                project_user = null_str if info.user is None else info.user
                library_type = null_str if info.library_type is None else info.library_type
                organism = null_str if info.organism is None else info.organism
                PI = null_str if info.PI is None else info.PI
                # Generate line in the table of projects
                index_page.add("<tr>")
                index_page.add("<td>%s</td>" % project.name)
                index_page.add("<td>%s</td>" % project_user)
                index_page.add("<td>%s</td>" % library_type)
                index_page.add("<td>%s</td>" % organism)
                index_page.add("<td>%s</td>" % PI)
                index_page.add("<td>%s</td>" % project.prettyPrintSamples())
                index_page.add("<td>%d</td>" % len(project.samples))
                # Locate and copy QC reports
                report_html = []
                for qc_dir in project_qc[project.name]:
                    qc_zip = os.path.join(project.dirn,
                                          "%s_report.%s.%s.zip" %
                                          (qc_dir,
                                           project.name,
                                           os.path.basename(self.analysis_dir)))
                    print qc_zip
                    assert(os.path.isfile(qc_zip))
                    qc_report_copied = True
                    try:
                        fileops.copy(qc_zip,dirn)
                        fileops.unzip(os.path.join(
                            dirn,
                            os.path.basename(qc_zip)),
                                      fileops.Location(dirn).path)
                    except Exception, ex:
                        print "Failed to copy QC report: %s" % ex
                        qc_report_copied = False
                    # Append info to the index page
                    if qc_report_copied:
                        qc_base = "%s_report" % qc_dir
                        fastq_dir = project.qc_info(qc_dir).fastq_dir
                        if fastq_dir != project.info.primary_fastq_dir:
                            fastq_set = fastq_dir
                        else:
                            fastq_set = None
                        # Index
                        report_html.append(
                            "<a href='%s.%s.%s/%s.html'>[Report%s]</a>"
                            % (qc_base,
                               project.name,
                               os.path.basename(self.analysis_dir),
                               qc_base,
                               (" (%s)" % fastq_set
                                if fastq_set is not None else "")))
                        # Zip file
                        report_html.append(
                            "<a href='%s'>[Zip%s]</a>"
                            % (os.path.basename(qc_zip),
                               (" (%s)" % fastq_dir
                                if fastq_set is not None else "")))
                    # MultiQC
                    multiqc_report = os.path.join(project.dirn,
                                                  "multi%s.html" % qc_base)
                    if os.path.exists(multiqc_report):
                        print "Found MultiQC report: %s" % multiqc_report
                        final_multiqc = "multi%s.%s.html" % (qc_base,project.name)
                        try:
                            fileops.copy(multiqc_report,
                                         os.path.join(dirn,final_multiqc))
                        except Exception, ex:
                            print "Failed to copy MultiQC report: %s" % ex
                            multiqc_report = None
                    else:
                        print "No MultiQC report found for %s" % \
                            os.path.basename(qc_dir)
                        multiqc_report = None
                    if multiqc_report:
                        report_html.append("<a href='%s'>[MultiQC%s]</a>"
                                           % (final_multiqc,
                                              (" (%s)" % fastq_dir
                                               if fastq_dir != 'fastqs' else "")))
                # Check there is something to add
                if not report_html:
                    report_html.append("QC reports not available")
                # Locate and copy ICell8 processing reports
                icell8_processing_zip = os.path.join(
                    project.dirn,
                    "icell8_processing.%s.%s.zip" %
                    (project.name,
                     os.path.basename(self.analysis_dir)))
                if os.path.isfile(icell8_processing_zip):
                    print icell8_processing_zip
                    icell8_report_copied = True
                    try:
                        fileops.copy(icell8_processing_zip,dirn)
                        fileops.unzip(os.path.join(
                            dirn,
                            os.path.basename(icell8_processing_zip)),
                                      fileops.Location(dirn).path)
                    except Exception as ex:
                        print "Failed to copy ICell8 report: %s" % ex
                        icell8_report_copied = False
                else:
                    # No ICell8 processing report to copy
                    icell8_report_copied = False
                # Append info to the index page
                if icell8_report_copied:
                    report_html.append(
                        "<a href='icell8_processing.%s.%s/" \
                        "icell8_processing.html'>" \
                        "[Icell8 processing]</a>" % \
                        (project.name,
                         os.path.basename(self.analysis_dir)))
                    report_html.append("<a href='%s'>[Zip]</a>" % \
                                       os.path.basename(icell8_processing_zip))
                # Locate and copy cellranger count reports
                cellranger_zip = os.path.join(project.dirn,
                                      "cellranger_count_report.%s.%s.zip" %
                                      (project.name,
                                       os.path.basename(self.analysis_dir)))
                if os.path.isfile(cellranger_zip):
                    print cellranger_zip
                    cellranger_report_copied = True
                    try:
                        fileops.copy(cellranger_zip,dirn)
                        fileops.unzip(os.path.join(
                            dirn,
                            os.path.basename(cellranger_zip)),
                                      fileops.Location(dirn).path)
                    except Exception as ex:
                        print "Failed to copy cellranger report: %s" % ex
                        cellranger_report_copied = False
                else:
                    # No cellranger count data to copy
                    cellranger_report_copied = False
                # Append info to the index page
                if cellranger_report_copied:
                    report_html.append(
                        "<a href='cellranger_count_report.%s.%s/" \
                        "cellranger_count_report.html'>" \
                        "[Cellranger count]</a>" % \
                        (project.name,
                         os.path.basename(self.analysis_dir)))
                    report_html.append("<a href='%s'>[Zip]</a>" % \
                                       os.path.basename(cellranger_zip))
                # Add to the index
                index_page.add("<td>%s</td>"
                               % null_str.join(report_html))
                index_page.add("</tr>")
            index_page.add("</table>")
        # Finish index page
        index_page.add("<p class='footer'>Generated by auto_process.py %s on %s</p>" % \
                       (get_version(),time.asctime()))
        # Copy to server
        index_html = os.path.join(self.tmp_dir,'index.html')
        index_page.write(index_html)
        fileops.copy(index_html,dirn)
        if processing_qc_html:
            fileops.copy(processing_qc_html,dirn)
        for qc_html in cellranger_qc_html:
            fileops.copy(qc_html,dirn)
        # Print the URL if given
        if self.settings.qc_web_server.url is not None:
            print "QC published to %s" % self.settings.qc_web_server.url
        
    def report(self,logging=False,summary=False,full=False,projects=False):
        # Report the contents of the run in various formats
        # Turn off saving of parameters etc
        self._save_params = False
        self._save_metadata = False
        report = None
        if logging:
            # Short form "logging"-style report
            report = self.report_logging_format()
        if summary:
            # Summary report listing projects one per line
            report = self.report_summary_format()
        if full:
            # More extensive report for one or 
            report = self.report_full_format()
        if projects:
            # Report projects one per line for injection into
            # spreadsheet
            report = self.report_projects()
        # Generate the report
        if report is not None:
            print report
        else:
            # Generate a verbose general report on the contents
            # of the data directory
            print "Run ref      : %s" % self.run_reference_id
            print "Directory    : %s" % self.analysis_dir
            print "Platform     : %s" % self.metadata.platform
            print "Unaligned dir: %s" % self.params.unaligned_dir
            if self.readme:
                print "README.txt found: %s" % self.readme
            if self.params.unaligned_dir is not None or \
               not os.path.exists(self.params.unaligned_dir):
                try:
                    illumina_data = self.load_illumina_data()
                    print "\nSummary of data in '%s' dir:\n" % self.params.unaligned_dir
                    for project in illumina_data.projects:
                        print "- %s" % IlluminaData.describe_project(project)
                except IlluminaData.IlluminaDataError,ex:
                    print "Failed to load data from %s:" % self.params.unaligned_dir
                    print "%s" % ex
            else:
                print "No information on source fastq data (no unaligned dir found)"
            print "\nAnalysis projects:"
            for project in self.get_analysis_projects():
                print "\n- %s" % project.name
                print "  %s" % ('-'*len(project.name),)
                print "  User    : %s" % project.info.user
                print "  PI      : %s" % project.info.PI
                print "  Library : %s" % project.info.library_type
                print "  Organism: %s" % project.info.organism
                print "  Dir     : %s" % os.path.basename(project.dirn)
                print "  #samples: %s" % len(project.samples)
                print "  Samples : %s" % project.prettyPrintSamples()
                print "  QC      : %s" % ('ok' if project.verify_qc() else 'not verified')
                print "  Comments: %s" % (project.info.comments)

    def report_logging_format(self):
        """Generate one-line report suitable for logging

        Generates a one-line report that can be used for logging, for
        example:

        Paired end: 'PJB': Peter Briggs, Mouse ChIP-seq (PI: P Briggs) (6 samples); ...

        The report is based on project directories that are located
        in the analysis directory, and not from other information
        (e.g. contents of the 'bcl2fastq' directory or other outputs
        from processing).

        Returns:
          String with the report text.

        """
        report = []
        analysis_dir = utils.AnalysisDir(self.analysis_dir)
        for p in analysis_dir.projects:
            report.append("'%s': %s, %s %s (PI: %s) (%d sample%s)" % \
                          (p.name,
                           p.info.user,
                           p.info.organism,
                           p.info.library_type,
                           p.info.PI,
                           len(p.samples),
                           's' if len(p.samples) > 1 else ''
                       ))
        report = '; '.join(report)
        # Paired end run?
        if analysis_dir.paired_end:
            endedness = "Paired end"
        else:
            endedness = "Single end"
        report = "%s: %s" % (endedness,report)
        return report

    def report_summary_format(self):
        """Generate summary report suitable for bioinformaticians

        Generates a multi-line report which gives general information
        about the run, plus one-line summaries for each project, plus
        any additional information that has been recorded.

        The general information includes:

        - Platform
        - Run name
        - Run reference id
        - Processing software
        - Assay (i.e. sequencing kit)

        For each project:

        - Project subdirectory
        - Researcher (aka user)
        - PI
        - Application (aka library type)
        - Organism
        - Number of samples

        Returns:
          String with the report text.

        """
        # Gather information
        analysis_dir = utils.AnalysisDir(self.analysis_dir)
        datestamp = None
        instrument = None
        run_number = None
        run_name = self.run_name
        try:
            datestamp,instrument,run_number = IlluminaData.split_run_name(run_name)
        except Exception, ex:
            logging.warning("Unable to extract information from run name '%s'" \
                            % run_name)
            logging.warning("Exception: %s" % ex)
        if self.metadata.platform is not None:
            platform = self.metadata.platform.upper()
        else:
            platform = 'unknown'
        if self.metadata.run_number is not None:
            run_number = self.metadata.run_number
        try:
            bcl2fastq_software = ast.literal_eval(
                self.metadata.bcl2fastq_software)
        except ValueError:
            bcl2fastq_software = None
        if self.metadata.assay is not None:
            assay = self.metadata.assay
        else:
            assay = ''
        # Generate report text
        report = []
        # Report header
        if datestamp and instrument and run_number:
            title = "%s run #%s datestamped %s" % (platform,
                                                   run_number,
                                                   datestamp)
        else:
            title = "%s" % os.path.basename(self.analysis_dir)
        report.append("%s\n%s" % (title,'='*len(title)))
        report.append("Run name : %s" % run_name)
        report.append("Reference: %s" % self.run_reference_id)
        report.append("Platform : %s" % platform)
        report.append("Directory: %s" % self.params.analysis_dir)
        report.append("Endedness: %s" % \
                      ('Paired end' if analysis_dir.paired_end else 'Single end'))
        report.append("Bcl2fastq: %s" %
                      ('Unknown' if not bcl2fastq_software else
                       "%s %s" % (bcl2fastq_software[1],
                                  bcl2fastq_software[2])))
        report.append("Assay    : %s" % assay)
        report.append("")
        # Projects
        report.append("%d project%s:" % (analysis_dir.n_projects,
                                         '' if analysis_dir.n_projects == 1 else 's'))
        rows = []
        comments = bcf_utils.OrderedDictionary()
        for project in analysis_dir.projects:
            project_data = dict(project=project.name)
            for item in ('user','PI','library_type','organism'):
                value = project.info[item]
                project_data[item] = value if value not in ('.','?') else \
                                    '<unspecified %s>' % item.lower()
            rows.append(("- '%s':" % project_data['project'],
                         project_data['user'],
                         project_data['organism'],
                         project_data['library_type'],
                         "%d sample%s" % (len(project.samples),
                                          's' if len(project.samples) > 1 else ''),
                         "(PI %s)" % project_data['PI']))
            if project.info.comments:
                comments[project.name] = project.info.comments
        report.append(utils.pretty_print_rows(rows))
        # Additional comments/notes
        if comments:
            width = max([len(x) for x in comments])
            report.append("")
            report.append("Additional notes/comments:")
            for project in comments:
                first_line = True
                for line in bcf_utils.split_into_lines(comments[project],70-width):
                    if first_line:
                        report.append("- %s%s: %s" % (project,
                                                      ' '*(width-len(project)),
                                                      line))
                        first_line = False
                    else:
                        report.append("  %s  %s" % (' '*width,line))
        return '\n'.join(report)

    def report_full_format(self):
        # Generate long form "full"-style report suitable for sending
        # to bioinformaticians
        # e.g. Peter Briggs data is now available at
        #
        # /path/to/data/140204_SQ12345_0001_AB12CDXYZ/PJB/
        #
        # The samples are:
        #
        # PJB1, PJBA1-4 (5 paired end samples, multiple fastqs per sample)
        #
        # Additionally:
        # Some extra information.
        #
        # Acquire data
        illumina_data = self.load_illumina_data()
        project_metadata = self.load_project_metadata(self.params.project_metadata)
        # Generate report text
        report = []
        for p in project_metadata:
            project = utils.AnalysisProject(p['Project'],
                                            os.path.join(self.params.analysis_dir,
                                                         p['Project']))
            # Title
            title = "%s %s %s data from %s run %s" % \
                          (project.info.user,
                           project.info.library_type,
                           project.info.organism,
                           self.params.platform.upper(),
                           os.path.basename(self.params.analysis_dir).split('_')[0])
            report.append("%s\n%s\n" % (title,'-'*len(title)))
            # Location
            report.append("The data for %(user)s's %(org)s %(lib)s is now "
                          "available at\n\n%(dirn)s\n" % \
                          dict(user=project.info.user,
                               dirn=project.dirn,
                               org=project.info.organism,
                               lib=project.info.library_type))
            # Samples
            report.append("The samples are:\n\n%s (%d%s sample%s%s)" % \
                          (project.prettyPrintSamples(),
                           len(project.samples),
                           " paired end" if project.info.paired_end else '',
                           's' if len(project.samples) > 1 else '',
                           ", multiple fastqs per sample" if project.multiple_fastqs else ''))
            # Additional information
            report.append("\nAdditional information:\n")
            report.append("Endedness:\t%s" % \
                          ('Paired end' if illumina_data.paired_end else 'Single end'))
            report.append("Comments :\t%s" % project.info.comments)
        report = '\n'.join(report)
        return report

    def report_projects(self):
        """Generate one line reports suitable for pasting into spreadsheet

        Generate one-line report for each each project with tab-separated
        data items, suitable for injection into a spreadsheet.

        Each line has the following information:

        - Run id e.g. HISEQ_140328
        - Run number
        - Source
        - Date
        - User
        - PI
        - Application
        - Organism
        - Platform
        - #Samples
        - PE (yes/no)
        - Samples
        
        Returns:
          String with the report text.

        """
        # Acquire data
        analysis_dir = utils.AnalysisDir(self.analysis_dir)
        # General information
        run_name = self.run_name
        try:
            datestamp,instrument,run_number = IlluminaData.split_run_name(run_name)
            run_number = run_number.lstrip('0')
        except Exception, ex:
            logging.warning("Unable to extract information from run name '%s'" \
                            % run_name)
            logging.warning("Exception: %s" % ex)
            date_stamp = ''
            run_number = ''
        if self.metadata.platform is not None:
            platform = self.metadata.platform.upper()
        else:
            platform = ''
        run_id = self.run_reference_id
        if self.metadata.run_number is not None:
            run_number = self.metadata.run_number
        if self.metadata.source is not None:
            data_source = self.metadata.source
        else:
            data_source = ''
        paired_end = 'yes' if analysis_dir.paired_end else 'no'
        # Generate report, one line per project
        report = []
        nprojects = len(analysis_dir.projects)
        if nprojects == 0:
            report.append("No projects found")
        else:
            report.append("%s project%s found" % (nprojects,
                                                  ('' if nprojects == 1
                                                   else 's')))
        for project in analysis_dir.projects:
            project_line = [run_id,str(run_number),data_source,'']
            info = project.info
            project_line.append('' if not info.user else info.user)
            project_line.append('' if not info.PI else info.PI)
            project_line.append('' if not info.library_type else info.library_type)
            project_line.append('' if not info.organism else info.organism)
            project_line.append(platform)
            project_line.append(str(len(project.samples)))
            project_line.append(paired_end)
            project_line.append(project.prettyPrintSamples())
            report.append('\t'.join(project_line))
        report = '\n'.join(report)
        return report

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
           utils.AnalysisProject(project_name,
                                 os.path.join(self.analysis_dir,
                                              project_name)).exists:
            raise Exception("Project called '%s' already exists" %
                            project_name)
        # Load target as a project
        project = utils.AnalysisProject(project_name,project_dir)
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
