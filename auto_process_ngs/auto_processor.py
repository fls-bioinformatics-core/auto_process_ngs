#!/usr/bin/env python
#
#     auto_processor.py: core module for automated processing of Illumina sequence data
#     Copyright (C) University of Manchester 2013-16 Peter Briggs
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
import bcftbx.IlluminaData as IlluminaData
import bcftbx.platforms as platforms
import bcftbx.TabFile as TabFile
import bcftbx.utils as bcf_utils
import bcftbx.htmlpagewriter as htmlpagewriter
from bcftbx.JobRunner import fetch_runner
import config
import applications
import utils
import simple_scheduler
import bcl2fastq_utils
import settings
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
            if self.has_parameter_file and self.params.data_dir is not None:
                # Get run name from data directory
                self.metadata['run_name'] = os.path.basename(
                    self.params.data_dir)
            elif self.analysis_dir.endswith('_analysis'):
                # Guess from analysis dir name
                self.metadata['run_name'] = os.path.basename(
                    self.analysis_dir[:len('_analysis')])
            else:
                # Unknown run name
                logging.warning("Unable to identify or guess run name")
        # Instrument-related metadata
        if self.metadata.instrument_name is None or \
           self.metadata.instrument_datestamp is None or \
           self.metadata.instrument_run_number is None:
            print "Attempting to set missing instrument metadata items"
            if self.metadata.run_name is not None:
                # Extract from run name
                try:
                    datestamp,instrument,run_number = \
                        IlluminaData.split_run_name(self.metadata.run_name)
                    if self.metadata.instrument_name is None:
                        self.metadata['instrument_name'] = instrument
                    if self.metadata.instrument_datestamp is None:
                        self.metadata['instrument_datestamp'] = datestamp
                    if self.metadata.instrument_run_number is None:
                        self.metadata['instrument_run_number'] = run_number
                except Exception, ex:
                    logging.warning("Unable to extract information from "
                                    "run name")
            else:
                # Unable to get missing data items
                logging.warning("Unable to set missing instrument metadata")

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
                fp.write("%s\n%s\n\n\n--" % (title,'='*len(title)))
        else:
            logging.warning("'%s' already exists" % self.readme_file)

    def edit_readme(self,editor='vi'):
        """
        Bring up README in an editor
        """
        if self.readme_file is None:
            logging.error("No README file to edit")
            return
        editor = editor.split(' ')
        edit_cmd = applications.Command(editor[0],*editor[1:])
        edit_cmd.add_args(self.readme_file)
        edit_cmd.run_subprocess()

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
        # Get the highest current step number
        # from the names of the existing subdirs
        i = 0
        for d in bcf_utils.list_dirs(self.log_dir):
            try:
                i = max(i,int(d.split('_')[0]))
            except ValueError:
                pass
        # Return the name
        return "%03d_%s" % (i+1,str(name))

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
        run_name = os.path.basename(self.analysis_dir)
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
            if run_name.endswith('_analysis'):
                # Strip trailing _analysis
                run_id = run_name[:-len('_analysis')]
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
                sample_sheet = os.path.join(data_dir,
                                            'Data','Intensities',
                                            'BaseCalls','SampleSheet.csv')
                tmp_sample_sheet = os.path.join(self.tmp_dir,'SampleSheet.csv')
            else:
                tmp_sample_sheet = os.path.join(self.tmp_dir,os.path.basename(sample_sheet))
            rsync = applications.general.rsync(sample_sheet,self.tmp_dir)
            print "%s" % rsync
            status = rsync.run_subprocess(log=self.log_path('rsync.sample_sheet.log'))
            if status != 0:
                logging.error("Failed to rsync sample sheet '%s' (status %s)"
                              % (sample_sheet,status))
                return status
            custom_sample_sheet = os.path.join(self.analysis_dir,
                                               'custom_SampleSheet.csv')
            sample_sheet = bcl2fastq_utils.make_custom_sample_sheet(
                tmp_sample_sheet,custom_sample_sheet)
            print "Keeping copy of original sample sheet"
            original_sample_sheet = os.path.join(self.analysis_dir,
                                                 'SampleSheet.orig.csv')
            os.rename(tmp_sample_sheet,original_sample_sheet)
            # Set the permissions for the original SampleSheet
            os.chmod(original_sample_sheet,0664)
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
        # Print the predicted ouputs
        projects = sample_sheet.predict_output()
        print "Predicted output from sample sheet:"
        print "Project\tSample\tFastq"
        for project in projects:
            project_name = project[8:]
            sample_names = []
            for sample in projects[project]:
                sample_name = sample[7:]
                for fastq_base in projects[project][sample]:
                    print "%s\t%s\t%s" % (project_name,sample_name,fastq_base)
                sample_names.append(sample_name)
        # Store the parameters
        self.params['data_dir'] = data_dir
        self.params['analysis_dir'] = self.analysis_dir
        self.params['sample_sheet'] = custom_sample_sheet
        self.params['bases_mask'] = bases_mask
        # Store the metadata
        self.metadata['run_name'] = os.path.basename(data_dir)
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
        clone_unaligned_dir = os.path.join(clone_dir,
                                           os.path.basename(self.params.unaligned_dir))
        if not copy_fastqs:
            # Link to unaligned dir
            print "Symlinking %s" % clone_unaligned_dir
            os.symlink(unaligned_dir,clone_unaligned_dir)
        else:
            # Copy unaligned dir
            print "Copying %s" % clone_unaligned_dir
            shutil.copytree(unaligned_dir,clone_unaligned_dir)
        # Copy additional files, if found
        for f in (self.params.sample_sheet,
                  self.params.stats_file,
                  self.params.project_metadata):
            srcpath = os.path.join(self.analysis_dir,f)
            if os.path.exists(srcpath):
                shutil.copy(srcpath,clone_dir)
        # Basic set of subdirectories
        d = AutoProcess(analysis_dir=clone_dir)
        d.log_dir
        d.script_code_dir

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
        # By default returns all projects within the analysis
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
        
    def make_fastqs(self,ignore_missing_bcl=False,ignore_missing_stats=False,
                    skip_rsync=False,remove_primary_data=False,
                    nprocessors=None,require_bcl2fastq_version=None,
                    unaligned_dir=None,sample_sheet=None,
                    bases_mask=None,no_lane_splitting=None,
                    minimum_trimmed_read_length=None,
                    mask_short_adapter_reads=None,
                    generate_stats=True,stats_file=None,
                    per_lane_stats_file=None,
                    report_barcodes=False,barcodes_file=None,
                    skip_bcl2fastq=False,only_fetch_primary_data=False,
                    runner=None):
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
          unaligned_dir       : if set then use this as the output directory for
                                bcl-to-fastq conversion. Default is 'bcl2fastq' (unless
                                an alternative is already specified in the config file)
          sample_sheet        : if set then use this as the input samplesheet
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
          skip_bcl2fastq      : if True then don't perform fastq generation
          only_fetch_primary_data: if True then fetch primary data, don't do anything else
          runner              : (optional) specify a non-default job runner to use for
                                fastq generation

        """
        #
        # Check for pre-existing bcl2fastq outputs
        if self.verify_bcl_to_fastq(unaligned_dir=unaligned_dir):
            print "Bcl to fastq outputs already present"
            # Check for project metadata file
            self.make_project_metadata_file()
            # (Re)generate stats?
            if generate_stats:
                self.generate_stats(stats_file=stats_file,
                                    per_lane_stats_file=per_lane_stats_file,
                                    unaligned_dir=unaligned_dir,
                                    nprocessors=nprocessors,
                                    runner=runner)
            return
        # Log dir
        self.set_log_dir(self.get_log_subdir('make_fastqs'))
        # Fetch primary data
        if not skip_rsync:
            if self.get_primary_data() != 0:
                logging.error("Failed to acquire primary data")
                raise Exception, "Failed to acquire primary data"
            if only_fetch_primary_data:
                return
        # Run bcl_to_fastq
        if not skip_bcl2fastq:
            try:
                self.bcl_to_fastq(require_bcl2fastq=require_bcl2fastq_version,
                                  unaligned_dir=unaligned_dir,
                                  sample_sheet=sample_sheet,
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
            if not self.verify_bcl_to_fastq():
                raise Exception("Bcl2fastq failed to produce expected outputs")
        # Generate statistics
        if generate_stats:
            self.generate_stats(stats_file=stats_file,
                                per_lane_stats_file=per_lane_stats_file,
                                unaligned_dir=unaligned_dir,
                                nprocessors=nprocessors,
                                runner=runner)
        # Count and report barcode sequences
        if report_barcodes:
            self.report_barcodes(barcodes_file)
        # Make a 'projects.info' metadata file
        self.make_project_metadata_file()
        # Remove primary data
        if remove_primary_data:
            self.remove_primary_data()

    def get_primary_data(self):
        """Acquire the primary sequencing data (i.e. BCL files)

        Copies the primary sequencing data (bcl files etc) to a local area
        using rsync.

        """
        data_dir = self.params.data_dir
        self.params["primary_data_dir"] = self.add_directory('primary_data')
        try:
            rsync = applications.general.rsync(data_dir,self.params.primary_data_dir,
                                               prune_empty_dirs=True,
                                               extra_options=('--copy-links',
                                                              '--include=*/',
                                                              '--include=Data/**',
                                                              '--include=RunInfo.xml',
                                                              '--include=SampleSheet.csv',
                                                              '--exclude=*'))
            print "Running %s" % rsync
            status = rsync.run_subprocess(log=self.log_path('rsync.primary_data.log'))
        except Exception, ex:
            logging.error("Exception getting primary data: %s" % ex)
            status = -1
        if status != 0:
            logging.error("Failed to acquire primary data (status %s)" % status)
        return status

    def bcl_to_fastq(self,require_bcl2fastq=None,unaligned_dir=None,
                     sample_sheet=None,bases_mask=None,
                     ignore_missing_bcl=False,ignore_missing_stats=False,
                     no_lane_splitting=None,minimum_trimmed_read_length=None,
                     mask_short_adapter_reads=None,nprocessors=None,
                     runner=None):
        """Generate FASTQ files from the raw BCL files

        Performs FASTQ generation from raw BCL files produced by an Illumina
        sequencer, by running the external 'bclToFastq.py' program (which
        wraps the 'configureBclToFastq' and 'make' steps).

        Arguments:
          require_bcl2fastq: if set then should be a string of the form
            '1.8.4' or '>2.0' explicitly specifying the version of
            bcl2fastq to use. (Default to use specifications from the
            settings)
          unaligned_dir: if set then use this as the output directory for
            bcl-to-fastq conversion. Default is 'bcl2fastq' (unless an
            alternative is already specified in the settings)
          sample_sheet: if set then use this as the input sample sheet file
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
        if unaligned_dir is not None:
            self.params['unaligned_dir'] = unaligned_dir
        elif self.params['unaligned_dir'] is None:
            self.params['unaligned_dir'] = 'bcl2fastq'
        # Examine primary data
        primary_data = os.path.join(self.params.primary_data_dir,
                                    os.path.basename(self.params.data_dir))
        if not os.path.isdir(primary_data):
            raise Exception("Missing primary data directory: %s" % primary_data)
        # Bases mask
        if bases_mask is None:
            bases_mask = self.params.bases_mask
        else:
            self.params['bases_mask'] = bases_mask
        # Check for basic information needed to do bcl2fastq conversion
        if self.params.data_dir is None:
            raise Exception, "No source data directory"
        if bases_mask is None:
            raise Exception, "No bases mask"
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
            print "Found available bcl2fastq packages:"
            for i,package in enumerate(bcl2fastq):
                bcl2fastq_info = bcl2fastq_utils.bcl_to_fastq_info(package)
                print "%s %s\t%s\t%s" % (('*' if i == 0 else ' '),
                                         bcl2fastq_info[0],
                                         bcl2fastq_info[1],
                                         bcl2fastq_info[2])
            bcl2fastq_exe = bcl2fastq[0]
            bcl2fastq_info = bcl2fastq_utils.bcl_to_fastq_info(bcl2fastq_exe)
        else:
            raise Exception("No appropriate bcl2fastq software located")
        # Store info on bcl2fastq package
        self.metadata['bcl2fastq_software'] = bcl2fastq_info
        # Sample sheet
        if sample_sheet is None:
            sample_sheet = self.params.sample_sheet
        if not os.path.isabs(sample_sheet):
            sample_sheet = os.path.join(self.analysis_dir,sample_sheet)
        if not os.path.isfile(sample_sheet):
            raise Exception("Missing sample sheet '%s'" % sample_sheet)
        self.params['sample_sheet'] = sample_sheet
        sample_sheet = self.params.sample_sheet
        # Generate temporary sample sheet with required format
        fmt = bcl2fastq_utils.get_required_samplesheet_format(bcl2fastq_info[2])
        timestamp = time.strftime("%Y%m%d%H%M%S")
        tmp_sample_sheet = os.path.join(self.tmp_dir,
                                        "SampleSheet.%s.csv" % timestamp)
        print "Generating '%s' format sample sheet: %s" % (fmt,tmp_sample_sheet)
        IlluminaData.SampleSheet(sample_sheet).write(tmp_sample_sheet,fmt=fmt)
        # Put a copy in the log directory
        shutil.copy(tmp_sample_sheet,self.log_dir)
        # Create bcl2fastq directory
        bcl2fastq_dir = self.add_directory(self.params.unaligned_dir)
        # Get info on the run
        illumina_run = IlluminaData.IlluminaRun(primary_data)
        # Determine initial number of mismatches
        nmismatches = bcl2fastq_utils.get_nmismatches(bases_mask)
        # Check for barcode collisions
        collisions = bcl2fastq_utils.check_barcode_collisions(sample_sheet,
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
                    sample_sheet,
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
        print "Primary data dir      : %s" % primary_data
        print "%s" % illumina_run.run_dir
        print "Platform              : %s" % illumina_run.platform
        print "Bcl format            : %s" % illumina_run.bcl_extension
        print "Source sample sheet   : %s" % sample_sheet
        print "Bases mask            : %s" % bases_mask
        print "Nmismatches           : %d" % nmismatches
        print "Nprocessors           : %s" % nprocessors
        print "Ignore missing bcl    : %s" % ignore_missing_bcl
        print "Ignore missing stats  : %s" % ignore_missing_stats
        print "No lane splitting     : %s" % no_lane_splitting
        print "Min trimmed read len  : %s" % minimum_trimmed_read_length
        print "Mask short adptr reads: %s" % mask_short_adapter_reads
        print "Output dir            : %s" % bcl2fastq_dir
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
        if minimum_trimmed_read_length:
            bcl2fastq.add_args('--minimum-trimmed-read-length',
                               minimum_trimmed_read_length)
        if mask_short_adapter_reads:
            bcl2fastq.add_args('--mask-short-adapter-reads',
                               mask_short_adapter_reads)
        bcl2fastq.add_args('--bcl2fastq_path',
                           bcl2fastq_exe,
                           primary_data,
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
        print "bcl2fastq completed"
        # Verify outputs
        try:
            illumina_data = self.load_illumina_data()
        except IlluminaData.IlluminaDataError,ex:
            logging.error("Unable to find outputs from bclToFastq (%s)" % ex)
            return
        if not IlluminaData.verify_run_against_sample_sheet(illumina_data,
                                                            sample_sheet):
            logging.error("Failed to verify bcl to fastq outputs against sample sheet")
            return

    def generate_stats(self,stats_file=None,per_lane_stats_file=None,
                       unaligned_dir=None,nprocessors=None,
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
                                                '--nprocessors',nprocessors,
                                                '--force')
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
        self.params['stats_file'] = stats_file
        self.params['per_lane_stats_file'] = per_lane_stats_file
        print "Statistics generation completed: %s" % self.params.stats_file

    def analyse_barcodes(self,unaligned_dir=None,lanes=None,truncate_barcodes=None,
                         nprocessors=None,runner=None):
        """Analyse the barcode sequences for FASTQs for each specified lane

        Run 'count_barcodes.py' for one or more lanes, to analyse the
        barcode index sequences in each lane.

        Arguments:
          unaligned_dir: if set then use this as the output directory for
            bcl-to-fastq conversion. Default is 'bcl2fastq' (unless an
            alternative is already specified in the settings)
          lanes: a list of lane numbers (integers) to perform the analysis
            for. Default is to analyse all lanes.
          truncate_barcodes: if set then truncate barcode sequences to the
            specified length for analysis
          runner: set a non-default job runner.
        
        """
        # Sort out parameters
        if unaligned_dir is not None:
            self.params['unaligned_dir'] = unaligned_dir
        elif self.params['unaligned_dir'] is None:
            self.params['unaligned_dir'] = 'bcl2fastq'
        # Load data and determine lanes
        illumina_data = self.load_illumina_data(unaligned_dir=unaligned_dir)
        lane_numbers = []
        for s in illumina_data.undetermined.samples:
            for f in s.fastq_subset(read_number=1,full_path=True):
                lane = str(IlluminaData.IlluminaFastq(f).lane_number)
            if lane not in lane_numbers:
                lane_numbers.append(lane)
        lane_numbers.sort()
        print "Found lanes: %s" % ','.join(lane_numbers)
        if lanes is None:
            lanes = lane_numbers
        else:
            lanes.sort()
        # Check there are some lanes to examine
        if len(lanes) == 0:
            print "No lanes specified/found: nothing to do"
            return
        # Create a subdirectory for barcode analysis
        output_dir = self.add_directory('barcode_analysis')
        # Base name for output counts
        counts_base = os.path.join(output_dir,'counts')
        # Log dir
        self.set_log_dir(self.get_log_subdir('analyse_barcodes'))
        # Set up runner
        # Use the stats runner for now
        if runner is not None:
            runner = fetch_runner(runner)
        else:
            runner = self.settings.runners.stats
        runner.set_log_dir(self.log_dir)
        # Schedule the jobs
        sched = simple_scheduler.SimpleScheduler(
            runner=runner,
            max_concurrent=self.settings.general.max_concurrent_jobs)
        sched.start()
        for lane in lanes:
            print "Starting analysis of barcodes for lane %s" % lane
            # Initial command
            barcode_cmd = applications.Command('count_barcodes.py',
                                               '-l',lane,
                                               '-r',os.path.join(output_dir,
                                                                 'report.lane%s' % lane),
                                               '-s',self.params.sample_sheet,
                                               '-t',0.01)
            # Truncate barcodes
            if truncate_barcodes is not None:
                barcode_cmd.add_args('-T',truncate_barcodes)
            # Multicore
            if nprocessors is not None:
                barcode_cmd.add_args('-N',nprocessors)
            # Look for an existing counts file
            counts_file = "%s.lane%s" % (counts_base,lane)
            if os.path.exists(counts_file):
                print "Using counts from existing file %s" % counts_file
                barcode_cmd.add_args('-c',counts_file)
            else:
                barcode_cmd.add_args('-o',counts_base)
            # Finish command and submit
            barcode_cmd.add_args(os.path.join(self.analysis_dir,
                                              unaligned_dir))
            sched.submit(barcode_cmd,
                         name='analyse_barcodes.lane%s' % lane)
        # Wait for the scheduler to run all jobs
        sched.wait()
        sched.stop()

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

    def verify_bcl_to_fastq(self,unaligned_dir=None):
        """Check that bcl to fastq outputs match sample sheet predictions

        Arguments:
          unaligned_dir

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
        # Try to create an IlluminaData object
        try:
            illumina_data = IlluminaData.IlluminaData(self.analysis_dir,
                                                      unaligned_dir=unaligned_dir)
        except IlluminaData.IlluminaDataError, ex:
            # Failed to initialise
            logging.debug("Failed to get information from %s: %s" % (bcl_to_fastq_dir,ex))
            return False
        # Do check
        return IlluminaData.verify_run_against_sample_sheet(illumina_data,
                                                            self.params.sample_sheet)

    def merge_fastq_dirs(self,primary_unaligned_dir,dry_run=True):
        # Combine multiple output directories from bcl2fastq into
        # a single directory
        # It is intended to combine the output from multiple runs
        # of bcl2fastq into a single 'unaligned'-equivalent directory
        # It operates in an automatic mode and should detect
        # additional 'unaligned' dirs on its own.
        # 
        # primary_unaligned_dir: the 'unaligned' dir that data from
        #                        from all others will be put into
        # dry_run: if True then just report operations
        #
        if primary_unaligned_dir is None:
            raise Exception,"Primary unaligned dir not defined"
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
                else:
                    print "* %s" % dirn
                    unaligned_dirs[dirn] = illumina_data
            except Exception, ex:
                logging.debug("Rejecting %s: %s" % (dirn,ex))
        # Check primary unaligned dir
        if primary_illumina_data is None:
            raise Exception, "Primary dir '%s' doesn't exist, or doesn't contain data?" % \
                primary_unaligned_dir
        # Is there anything to do?
        if not unaligned_dirs:
            print "No extra bcl2fastq output directories found, nothing to do"
            return
        # Make a directory to move replaced data to
        if not dry_run:
            unaligned_backup = self.add_directory("save.%s" % 
                                                  os.path.basename(primary_unaligned_dir))
        # Examine each additional directory and move data as required
        for unaligned_dir in unaligned_dirs:
            # Deal with projects
            print "Importing projects from %s:" % unaligned_dir
            illumina_data = unaligned_dirs[unaligned_dir]
            for project in illumina_data.projects:
                try:
                    primary_project = primary_illumina_data.get_project(project.name)
                    print "- '%s' will be replaced by version from %s" % (project.name,
                                                                          unaligned_dir)
                    if not dry_run:
                        backup_dirn = os.path.join(unaligned_backup,
                                                   "save.%s" % os.path.basename(project.dirn))
                        print "- moving %s to %s" % (primary_project.dirn,
                                                     backup_dirn)
                        shutil.move(primary_project.dirn,backup_dirn)
                except IlluminaData.IlluminaDataError:
                    print "- '%s' will be imported from %s"  % (project.name,
                                                                unaligned_dir)
                if not dry_run:
                    print "- moving %s to %s" % (project.dirn,
                                                 primary_unaligned_dir)
                    shutil.move(project.dirn,primary_unaligned_dir)
            # Deal with undetermined indices
            if illumina_data.undetermined is not None:
                print "Importing undetermined indices data from %s:" % unaligned_dir
                for lane in illumina_data.undetermined.samples:
                    primary_lane = None
                    for lane0 in primary_illumina_data.undetermined.samples:
                        if lane.name == lane0.name:
                            primary_lane = lane0
                            break
                    # Finished looking for existing data
                    if primary_lane is not None:
                        print "- '%s' will be replaced by version from %s" % (lane.name,
                                                                              unaligned_dir)
                        if not dry_run:
                                backup_dirn = os.path.join(
                                    unaligned_backup,
                                    "save.%s" % os.path.basename(
                                        primary_illumina_data.undetermined.dirn))
                                print "- moving %s to %s" % (primary_lane.dirn,
                                                             backup_dirn)
                                shutil.move(primary_lane.dirn,backup_dirn)
                    else:
                        print "- '%s' will be imported from %s"  % (lane.name,
                                                                    unaligned_dir)
                    if not dry_run:
                        print "- moving %s to %s" % \
                            (lane.dirn,primary_illumina_data.undetermined.dirn)
                        shutil.move(lane.dirn,primary_illumina_data.undetermined.dirn)
            else:
                print "No undetermined indices found"
        if dry_run:
            return
        # Reset the bcl2fastq dir
        self.params['unaligned_dir'] = primary_unaligned_dir
        # Make a 'projects.merged.info' metadata file
        merged_project_metadata_file='projects.merged.info'
        self.make_project_metadata_file(merged_project_metadata_file)

    def setup_analysis_dirs(self,ignore_missing_metadata=False,
                            short_fastq_names=False,
                            link_to_fastqs=False):
        # Construct and populate the analysis directories for each project
        # ignore_missing_metadata: if set True then make projects even if
        #                          metadata hasn't been set (defaults to False
        #                          i.e. stop if metadata isn't set)
        if self.params.unaligned_dir is None:
            logging.error("No unaligned directory, cannot build analysis directories")
            raise Exception,"Cannot build analysis directories"
        illumina_data = self.load_illumina_data()
        # Project metadata file
        project_metadata_file = self.params.project_metadata
        if project_metadata_file is None:
            project_metadata_file = 'projects.info'
        if not os.path.exists(os.path.join(self.params.analysis_dir,project_metadata_file)):
            logging.warning("No project metadata file '%s' found, attempting to create"
                            % project_metadata_file)
            self.make_project_metadata_file(project_metadata_file)
            logging.warning("Update '%s' and rerun" % project_metadata_file)
            return
        project_metadata = self.load_project_metadata(project_metadata_file=project_metadata_file,
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
                raise Exception, "Missing project metadata"
        # Create the projects
        n_projects = 0
        for line in project_metadata:
            # Acquire the run name
            if self.params.data_dir is not None:
                run_name = os.path.basename(self.params.data_dir)
            else:
                run_name = os.path.basename(self.params.analysis_dir)
            # Look up project data
            project_name = line['Project']
            user = line['User']
            PI = line['PI']
            organism = line['Organism']
            library_type = line['Library']
            comments = line['Comments']
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
                logging.warning("Project '%s' already exists, skipping" % project.name)
                continue
            print "Creating project: '%s'" % project_name
            project.create_directory(illumina_data.get_project(project_name),
                                     short_fastq_names=short_fastq_names,
                                     link_to_fastqs=link_to_fastqs)
            n_projects += 1
        # Tell us how many were made
        print "Created %d project%s" % (n_projects,'s' if n_projects != 1 else '')
        # Also set up analysis directory for undetermined reads
        undetermined = illumina_data.undetermined
        if illumina_data.undetermined is not None:
            undetermined = utils.AnalysisProject('undetermined',
                                                 os.path.join(self.analysis_dir,
                                                              'undetermined'),
                                                 run=run_name,
                                                 comments="Analysis of reads "
                                                 "with undetermined indices",
                                                 platform=self.metadata.platform)
            if not undetermined.exists:
                print "Creating directory 'undetermined' for analysing reads " \
                "with undetermined indices"
                undetermined.create_directory(illumina_data.undetermined,
                                              link_to_fastqs=link_to_fastqs)
            else:
                logging.warning("'undetermined' directory already exists, skipping")

    def run_qc(self,projects=None,max_jobs=4,ungzip_fastqs=False,
               fastq_screen_subset=1000000,runner=None):
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
                    (default is 1000000, set to zero or None to use
                    all reads)
          runner:   (optional) specify a non-default job runner to
                    use for the QC.
        Returns:
          UNIX-style integer returncode: 0 = successful termination,
          non-zero indicates an error occurred.

        """
        # Check QC script version
        compatible_versions = ('1.3.0',)
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
            # Make the qc directory if it doesn't exist
            qc_dir = os.path.join(project.dirn,'qc')
            if not os.path.exists(qc_dir):
                print "Making 'qc' subdirectory"
                bcf_utils.mkdir(qc_dir,mode=0775)
            # Set up the logs directory
            log_dir = os.path.join(qc_dir,'logs')
            if not os.path.exists(log_dir):
                print "Making 'qc/logs' subdirectory"
                bcf_utils.mkdir(log_dir,mode=0775)
            # Loop over samples and queue up those where the QC
            # isn't validated
            samples = project.get_samples(sample_pattern)
            if len(samples) == 0:
                logging.warning("No samples found for QC analysis in project '%s'" %
                                project.name)
            for sample in samples:
                group = None
                print "Examining files in sample %s" % sample.name
                for fq in sample.fastq:
                    if sample.verify_qc(qc_dir,fq):
                        logging.debug("\t%s: QC verified" % fq)
                    else:
                        print "\t%s: setting up QC run" % os.path.basename(fq)
                        # Create a group if none exists for this sample
                        if group is None:
                            group = sched.group("%s.%s" % (project.name,sample.name),
                                                log_dir=log_dir)
                        # Create and submit a QC job
                        fastq = os.path.join(project.dirn,'fastqs',fq)
                        label = "illumina_qc.%s" % str(utils.AnalysisFastq(fq))
                        qc_cmd = applications.Command('illumina_qc.sh',fastq)
                        if ungzip_fastqs:
                            qc_cmd.add_args('--ungzip-fastqs')
                        if fastq_screen_subset is None:
                            fastq_screen_subset = 0
                        qc_cmd.add_args('--subset',fastq_screen_subset)
                        job = group.add(qc_cmd,name=label,wd=project.dirn)
                        print "Job: %s" %  job
                # Indicate no more jobs to add
                if group: group.close()
        # Wait for the scheduler to run all jobs
        sched.wait()
        sched.stop()
        # Verify the outputs and generate QC reports
        failed_projects = []
        for project in projects:
            if not project.verify_qc():
                failed_projects.append(project)
            else:
                print "QC okay, generating report for %s" % project.name
                project.qc_report
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
                        read_only_fastqs=True,force=False):
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
                    '--exclude=save.*']
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
        # Log dir
        self.set_log_dir(self.get_log_subdir('archive'))
        # If making fastqs read-only then transfer them separately
        if read_only_fastqs:
            try:
                rsync_fastqs = applications.general.rsync(self.analysis_dir,archive_dir,
                                                          prune_empty_dirs=True,
                                                          dry_run=dry_run,
                                                          chmod='ugo-w',
                                                          extra_options=('--include=*/',
                                                                         '--include=fastqs/**',
                                                                         '--exclude=*',))
                print "Running %s" % rsync_fastqs
                status = rsync_fastqs.run_subprocess(
                    log=self.log_path('rsync.archive_fastqs.log'))
            except Exception, ex:
                logging.error("Exception rsyncing fastqs to archive: %s" % ex)
                status = -1
            if status != 0:
                logging.error("Failed to rsync fastqs to archive (returned status %d)" %
                              status)
            # Exclude from main rsync
            excludes.append('--exclude=fastqs')
            # Set returncode
            retval = status if status != 0 else retval
        # Run the rsync command
        try:
            rsync = applications.general.rsync(self.analysis_dir,archive_dir,
                                               prune_empty_dirs=True,
                                               dry_run=dry_run,
                                               chmod=chmod,
                                               extra_options=excludes)
            print "Running %s" % rsync
            status = rsync.run_subprocess(log=self.log_path('rsync.archive.log'))
        except Exception, ex:
            logging.error("Exception rsyncing to archive: %s" % ex)
            status = -1
        if status != 0:
            logging.error("Failed to rsync to archive (returned status %d)" % status)
        # Set returncode
        retval = status if status != 0 else retval
        # Set the group (local copies only)
        if group is not None:
            user,server,dirn = utils.split_user_host_dir(archive_dir)
            if user is None and server is None:
                # Local archive
                print "Setting group of archived files to '%s'" % group
                gid = bcf_utils.get_gid_from_group(group)
                if gid is None:
                    logging.error("Failed to get gid for group '%s'" % group)
                    retval = 1
                else:
                    for f in bcf_utils.walk(
                            os.path.join(dirn,os.path.basename(self.analysis_dir)),
                            include_dirs=True):
                        logging.debug("Updating group for %s" % f)
                        if not dry_run:
                            os.lchown(f,-1,gid)
        # Finish
        return retval

    def log_analysis(self):
        # Add a record of the analysis to the logging file
        raise NotImplementedError

    def publish_qc(self,projects=None,location=None,ignore_missing_qc=False,
                   regenerate_reports=False):
        # Copy the QC reports to the webserver
        #
        # projects: specify a pattern to match one or more projects to
        #           publish the reports for (default is to publish all reports)
        # location: override the target location specified in the settings
        #           can be of the form '[[user@]server:]directory'
        # ignore_missing_qc: if True then skip directories with missing QC data
        #           or reports (otherwise raises an exception)
        # regenerate_reports: if True then try to create reports even when they
        #           already exist
        #
        # Turn off saving of parameters etc
        self._save_params = False
        self._save_metadata = False
        # Check metadata
        check_metadata = self.check_metadata(('source','run_number'))
        if not check_metadata:
            logging.error("Some metadata items not set, stopping")
            return
        # Process pattern matching
        if projects is None:
            project_pattern = '*'
        else:
            project_pattern = projects
        # Get location to publish qc reports to
        if location is None:
            user,server,dirn = utils.split_user_host_dir(
                self.settings.qc_web_server.dirn)
        else:
            user,server,dirn = utils.split_user_host_dir(location)
        if server is not None:
            remote = True
        else:
            remote = False
        # Check the settings
        if remote:
            print "Copying QC to remote directory"
            print "user:\t%s" % user
            print "host:\t%s" % server
            print "dirn:\t%s" % dirn
        else:
            print "Copying QC to local directory"
            print "dirn:\t%s" % dirn
        if dirn is None:
            raise Exception, "No target directory specified"
        dirn = os.path.join(dirn,os.path.basename(self.analysis_dir))
        # Get general data
        analysis_dir = utils.AnalysisDir(self.analysis_dir)
        # Get project data
        projects = analysis_dir.get_projects(project_pattern)
        # Check QC situation for each project
        print "Checking QC status for each project:"
        no_qc_projects = []
        for project in projects:
            if project.qc is None:
                # No QC available, can't even report
                logging.warning("No QC available for %s" % project.name)
                no_qc_projects.append(project)
            else:
                # QC is available, check status of reports
                generate_report = regenerate_reports
                qc_zip = os.path.join(project.dirn,"%s.zip" % project.qc.report_name)
                if os.path.isfile(qc_zip):
                    print "Existing QC report found for %s" % project.name
                else:
                    # No QC report available
                    print "No QC report found for %s" % project.name
                    generate_report = True
                # (Re)create report
                if generate_report:
                    if project.qc.verify():
                        print "Generating report and zip file"
                        try:
                            project.qc.zip()
                        except Exception, ex:
                            logging.error("Failed to generate QC report for %s" %
                                          project.name)
                            qc_zip = None
                    else:
                        logging.error("Unable to verify QC for %s" % project.name)
                        qc_zip = None
                if qc_zip is None:
                    logging.error("Failed to make QC report for %s" % project.name)
                    no_qc_projects.append(project)
        # Final results
        if no_qc_projects:
            logging.error("QC reports missing for projects: %s" %
                          ', '.join([x.name for x in no_qc_projects]))
            if not ignore_missing_qc:
                raise Exception, "QC reports missing for projects: %s" % \
                    ', '.join([x.name for x in no_qc_projects])
        # Remove the 'bad' projects from the list before proceeding
        for project in no_qc_projects:
            print "Project %s will be skipped" % project.name
            projects.remove(project)
        if not projects:
            logging.error("No projects with QC results to publish")
            raise Exception, "No projects with QC results to publish"
        # Make a directory for the QC reports
        if not remote:
            # Local directory
            bcf_utils.mkdir(dirn)
        else:
            # Remote directory
            try:
                mkdir_cmd = applications.general.ssh_command(user,server,('mkdir',dirn))
                print "Running %s" % mkdir_cmd
                mkdir_cmd.run_subprocess()
            except Exception, ex:
                raise Exception, "Exception making remote directory for QC reports: %s" % ex
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
        # Table of projects
        index_page.add("<h2>QC Reports</h2>")
        index_page.add("<table>")
        index_page.add("<tr><th>Project</th><th>User</th><th>Library</th><th>Organism</th><th>PI</th><th>Samples</th><th>#Samples</th><th colspan='2'>Reports</th></tr>")
        # Set the string to represent "null" table entries
        null_str = '&nbsp;'
        # Deal with QC for each project
        for project in projects:
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
            # Locate and copy QC report
            qc_zip = os.path.join(project.dirn,"%s.zip" % project.qc.report_name)
            assert(os.path.isfile(qc_zip))
            report_copied = True
            if not remote:
                # Local directory
                shutil.copy(qc_zip,dirn)
                # Unpack
                unzip_cmd = applications.Command('unzip','-q','-o','-d',dirn,qc_zip)
                print "Running %s" % unzip_cmd
                unzip_cmd.run_subprocess()
            else:
                try:
                    # Remote directory
                    scp = applications.general.scp(user,server,qc_zip,dirn)
                    print "Running %s" % scp
                    scp.run_subprocess()
                    # Unpack at the other end
                    unzip_cmd = applications.general.ssh_command(
                        user,server,
                        ('unzip','-q','-o','-d',dirn,
                         os.path.join(dirn,os.path.basename(qc_zip))))
                    print "Running %s" % unzip_cmd
                    unzip_cmd.run_subprocess()
                except Exception, ex:
                    print "Failed to copy QC report: %s" % ex
                    report_copied = False
            # Append info to the index page
            if report_copied:
                index_page.add("<td><a href='%s/qc_report.html'>[Report]</a></td>"
                               % project.qc.report_name)
                index_page.add("<td><a href='%s'>[Zip]</a></td>"
                               % os.path.basename(qc_zip))
            else:
                # QC not available
                index_page.add("<td colspan='2'>QC reports not available</td>")
            # Finish table row for this project
            index_page.add("</tr>")
        index_page.add("</table>")
        # Add table of statistics
        stats_table = self.stats_as_html_table()
        if stats_table is not None:
            index_page.add("<h2>Statistics for Fastq files</h2>")
            index_page.add(stats_table)
        # Add per-lane statistics
        per_lane_stats_table = self.per_lane_stats_as_html_table()
        if per_lane_stats_table is not None:
            index_page.add("<h2>Per-lane statistics</h2>")
            index_page.add(per_lane_stats_table)
        # Finish index page
        index_page.add("<p class='footer'>Generated by auto_process.py %s on %s</p>" % \
                       (get_version(),time.asctime()))
        # Copy to server
        index_html = os.path.join(self.tmp_dir,'index.html')
        index_page.write(index_html)
        if not remote:
            # Local directory
            shutil.copy(index_html,dirn)
        else:
            # Remote directory
            scp = applications.general.scp(user,server,index_html,dirn)
            print "Running %s" % scp
            scp.run_subprocess()
        # Print the URL if given
        if self.settings.qc_web_server.url is not None:
            print "QC published to %s" % self.settings.qc_web_server.url

    def stats_as_html_table(self):
        # Return statistics as string containing HTML table, or None if no
        # stats file was found
        stats_file = self.params.stats_file
        if stats_file is None or not os.path.exists(os.path.join(self.analysis_dir,
                                                                 stats_file)):
            logging.warning("No statistics file found")
            return None
        else:
            stats_file = os.path.join(self.analysis_dir,stats_file)
        # Load the statistics and dump as HTML table
        header = ('Project','Sample','Fastq','Size','Nreads')
        table = ['<table class="stats">']
        stats = TabFile.TabFile(stats_file,first_line_is_header=True)
        html_line = ['<tr>']
        for field in header:
            html_line.append("<th>%s</th>" % field)
        html_line.append('</tr>')
        table.append(''.join(html_line))
        last_project = None
        last_sample = None
        for line in stats:
            html_line = ['<tr>']
            for field in header:
                if field == 'Project':
                    if  line[field] == last_project:
                        data = '&nbsp;'
                    else:
                        data = line[field]
                        last_project = data
                elif field == 'Sample':
                    if line[field] == last_sample:
                        data = '&nbsp;'
                    else:
                        data = line[field]
                        last_sample = data
                else:
                    data = line[field]
                html_line.append("<td>%s</td>" % data)
            html_line.append('</tr>')
            table.append(''.join(html_line))
        table.append('</table>')
        return '\n'.join(table)

    def per_lane_stats_as_html_table(self):
        # Return per-lane statistics as string containing HTML table, or None if
        # no per-lane stats file was found
        per_lane_stats_file = self.params.per_lane_stats_file
        if not per_lane_stats_file:
            per_lane_stats_file = 'per_lane_stats.info'
        if not os.path.exists(os.path.join(self.analysis_dir,
                                           per_lane_stats_file)):
            logging.warning("No per-lane statistics file found")
            return None
        per_lane_stats_file = os.path.join(self.analysis_dir,
                                           per_lane_stats_file)
        # Load the statistics and dump as HTML table
        table = ['<table class="per_lane_stats">']
        with open(per_lane_stats_file,'r') as fp:
            for line in fp:
                line = line.rstrip('\n')
                if not line:
                    continue
                html_line = ['<tr>']
                if line.startswith('Lane '):
                    html_line.append('<td colspan="4">%s</td>' % line)
                elif line.startswith('Total reads = '):
                    html_line.append('<td colspan="4">%s</td>' % line)
                elif line.startswith('- '):
                    items = line[2:].split('\t')
                    project,sample = items[0].split('/')
                    nreads = items[1]
                    percent_reads = items[2]
                    html_line.append('<td>%s</td>' % project)
                    html_line.append('<td>%s</td>' % sample)
                    html_line.append('<td>%s</td>' % nreads)
                    html_line.append('<td>%s</td>' % percent_reads)
                html_line.append('<tr>')
                table.append(''.join(html_line))
        table.append('</table>')
        return '\n'.join(table)
        
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
        if self.params.data_dir is not None:
            run_name = os.path.basename(self.params.data_dir)
            try:
                datestamp,instrument,run_number = IlluminaData.split_run_name(run_name)
            except Exception, ex:
                logging.warning("Unable to extract information from run name '%s'" \
                                % run_name)
                logging.warning("Exception: %s" % ex)
        else:
            run_name = os.path.basename(self.analysis_dir)
            if run_name.endswith('_analysis'):
                # Strip trailing _analysis
                run_name = run_name[:-len('_analysis')]
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
        run_name = os.path.basename(self.analysis_dir)
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
