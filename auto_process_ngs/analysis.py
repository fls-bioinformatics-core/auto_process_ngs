#!/usr/bin/env python
#
#     analysis: classes & funcs for handling analysis dirs and projects
#     Copyright (C) University of Manchester 2018-2025 Peter Briggs
#
########################################################################
#
# analysis.py
#
#########################################################################

"""
Classes and functions for handling analysis directories and
projects.

Classes:

- AnalysisFastq: extract information from Fastq name
- AnalysisDir: API for sequencing run analysis directory
- AnalysisProject: API for project in an analysis directory
- AnalysisSample: API for sample in an analysis project

Functions:

- run_id: fetch run ID for sequencing run
- run_reference_id: fetch run reference ID for sequencing run
- split_sample_name: split sample name into components
- split_sample_reference: split sample reference ID into components
- match_run_id: check if directory matches run identifier
- locate_run: search for an analysis directory by ID
- locate_project: search for an analysis project by ID
- locate_project_info_file: search for 'project.info' file
- copy_analysis_project: make copy of an AnalysisDir
"""

#######################################################################
# Imports
#######################################################################

import os
import re
import fnmatch
import logging
import bcftbx.IlluminaData as IlluminaData
import bcftbx.Pipeline as Pipeline
import bcftbx.utils as bcf_utils
from bcftbx.Md5sum import md5sum
from .metadata import AnalysisDirMetadata
from .metadata import AnalysisDirParameters
from .metadata import AnalysisProjectInfo
from .metadata import ProjectMetadataFile
from .metadata import AnalysisProjectQCDirInfo
from .fastq_utils import BaseFastqAttrs
from .fastq_utils import IlluminaFastqAttrs
from .utils import sort_sample_names
from functools import reduce

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Constants
#######################################################################

# Regular expression for matching SRA-style Fastq names
SRA_REGEX = re.compile(
    r"(?P<sample_name>(SRR|ERR)[0-9]+)(_(?P<read_number>[1-2]))?$")

#######################################################################
# Classes
#######################################################################

class AnalysisFastq(BaseFastqAttrs):
    """
    Class for extracting information about Fastq files

    Given the name of a Fastq file, extract data about the sample name,
    barcode sequence, lane number, read number and set number.

    Uses the IlluminaFastqAttrs class to handle Fastq filenames
    which consist of a valid Fastq name as defined by
    IlluminaFastqAttrs, but with additional elements appended.

    Instances of this class have the following attributes
    (defined in the base class):

    * fastq:            the original fastq file name
    * basename:         basename with NGS extensions stripped
    * extension:        full extension e.g. '.fastq.gz'
    * sample_name:      name of the sample
    * sample_number:    integer (or None if no sample number)
    * barcode_sequence: barcode sequence (string or None)
    * lane_number:      integer (or None if no lane number)
    * read_number:      integer (or None if no read number)
    * set_number:       integer (or None if no set number)
    * is_index_read:    boolean (True if index read, False if not)

    There are four additional attributes:

    * format: string identifying the format of the Fastq name
      ('Illumina', 'SRA', or None)
    * implicit_read_number: flag indicating whether the read
      number was implied (i.e. doesn't appear explicitly in the
      name)
    * canonical_name: the 'canonical' part of the name (string,
      or None if no canonical part could be extracted)
    * extras: the 'extra' part of the name (string, or None if
      there was no trailing extra part)

    Arguments:
      fastq (str): path or name of Fastq file
    """
    def __init__(self,fastq):
        """
        Create a new AnalysisFastq instance
        """
        # Initialise the base class
        BaseFastqAttrs.__init__(self,fastq)
        # Additional attributes
        self.format = None
        self.extras = None
        self.implicit_read_number = False
        # Check for SRA format
        fq = SRA_REGEX.match(self.basename)
        if fq is not None:
            fq_attrs = fq.groupdict()
            self.format = "SRA"
            self.sample_name = fq_attrs['sample_name']
            if fq_attrs['read_number'] is not None:
                self.read_number = int(fq_attrs['read_number'])
            else:
                # Assume that it's implicitly read 1
                self.read_number = 1
                # Set flag to indicate that the read number
                # is implicit (i.e. not in the name)
                self.implicit_read_number = True
            return
        # Try to identify a canonical component of the
        # name by removing trailing parts and testing
        fastq_base = self.basename
        extras = []
        fq_attrs = IlluminaFastqAttrs(fastq_base)
        while fq_attrs.sample_name == fq_attrs.basename:
            # Remove more trailing components and try again
            extras.append(fastq_base.split('_')[-1])
            fastq_base = '_'.join(fastq_base.split('_')[:-1])
            if not fastq_base:
                # No more components to remove
                fq_attrs = IlluminaFastqAttrs(self.basename)
                extras = None
                break
            else:
                fq_attrs = IlluminaFastqAttrs(fastq_base)
        # Copy attributes across
        self.sample_name = fq_attrs.sample_name
        self.sample_number = fq_attrs.sample_number
        self.barcode_sequence = fq_attrs.barcode_sequence
        self.lane_number = fq_attrs.lane_number
        self.read_number = fq_attrs.read_number
        self.set_number = fq_attrs.set_number
        self.is_index_read = fq_attrs.is_index_read
        self.format = "Illumina"
        # Store extras
        if extras:
            self.extras = "_%s" % '_'.join(extras[::-1])
        else:
            self.extras = None
    @property
    def canonical_name(self):
        """
        Return the 'canonical' part of the name
        """
        if self.format == "Illumina":
            illumina_fastq_attrs = IlluminaFastqAttrs(self.fastq)
            illumina_fastq_attrs.sample_name = self.sample_name
            illumina_fastq_attrs.sample_number = self.sample_number
            illumina_fastq_attrs.barcode_sequence = self.barcode_sequence
            illumina_fastq_attrs.lane_number = self.lane_number
            illumina_fastq_attrs.read_number = self.read_number
            illumina_fastq_attrs.set_number = self.set_number
            illumina_fastq_attrs.is_index_read = self.is_index_read
            if illumina_fastq_attrs.sample_name != \
               illumina_fastq_attrs.basename:
                return str(illumina_fastq_attrs)
        elif self.format == "SRA":
            return "%s%s" % (self.sample_name,
                             "_%s" % self.read_number
                             if self.read_number and
                             not self.implicit_read_number
                             else '')
        return None
    def __repr__(self):
        """
        Implement __repr__ built-in
        """
        if self.canonical_name:
            return "%s%s" % (self.canonical_name,
                             self.extras if self.extras else '')
        else:
            return self.basename

class AnalysisDir:
    """
    Class describing an analysis directory

    Conceptually an analysis directory maps onto a sequencing run.
    It consists of one or more sets of samples from that run,
    which are represented by subdirectories.

    It is also possible to have one or more subdirectories containing
    outputs from the CASAVA or bclToFastq processing software.

    Properties:

    * analysis_dir: full path to the directory
    * run_name: the name of the parent sequencing run
    * run_id: the ID assigned to the run
    * run_reference_id: the reference ID assigned to the run
    * metadata: metadata items associated with the run
    * projects: list of AnalysisProject objects
    * undetermined: AnalysisProject object for 'undetermined'
    * sequencing_data: list of IlluminaData objects
    * projects_metadata: metadata from the 'projects.info' file
    * datestamp: datestamp extracted from run name
    * instrument_name: instrument name extracted from run name
    * instrument_run_number: run number extracted from run name
    * n_projects: number of projects
    * n_sequencing_data: number of sequencing data directories
    * paired_end: whether data are paired ended

    Arguments:
      analysis_dir (str): name (and path) to analysis directory
    """
    def __init__(self,analysis_dir):
        """
        Create a new AnalysisDir instance
        """
        # Store location
        self._analysis_dir = os.path.abspath(analysis_dir)
        self._name = os.path.basename(analysis_dir)
        self._bcl2fastq_dirs = []
        self._project_dirs = []
        self._extra_dirs = []
        self.sequencing_data = []
        self.projects = []
        self.undetermined = None
        # Metadata
        self.metadata = AnalysisDirMetadata()
        metadata_file = os.path.join(self._analysis_dir,
                                     "metadata.info")
        if os.path.exists(metadata_file):
            try:
                self.metadata.load(metadata_file)
            except Exception as ex:
                logger.debug("Failed to load metadata from %s: %s" %
                             (metadata_file,ex))
        if not self.metadata:
            # No metadata loaded, try parameter file instead
            logger.debug("Attempting to load parameter file")
            parameter_file = os.path.join(self._analysis_dir,
                                          "auto_process.info")
            params = AnalysisDirParameters()
            if os.path.exists(parameter_file):
                try:
                    params.load(parameter_file,strict=False)
                except Exception as ex:
                    logger.warning("Failed to load parameters from %s: %s"
                                   % (parameter_file,ex))
            if params:
                # Attempt to acquire values from parameters
                for param in ('platform','run_number','source',):
                    if param not in params:
                        logger.debug("-- %s: missing" % param)
                        continue
                    logger.debug("-- %s: setting to '%s'" % (param,
                                                             params[param]))
                    self.metadata[param] = params[param]
            else:
                # No metadata or parameter data loaded
                raise Exception("Failed to load parameters: %s is not an "
                                "analysis directory?" % self._analysis_dir)
        # Projects metadata
        try:
            self.projects_metadata = ProjectMetadataFile(
                os.path.join(self._analysis_dir,"projects.info"))
        except Exception as ex:
            logger.warning("Failed to load projects metadata: %s" % ex)
            self.projects_metadata = None
        # Run name
        try:
            self.run_name = self.metadata.run_name
        except AttributeError:
            self.run_name = None
        if not self.run_name:
            self.run_name = self._analysis_dir[0:-len('_analysis')]
        self.run_name = os.path.basename(self.run_name)
        self.date_stamp,\
            self.instrument_name,\
            self.instrument_run_number = IlluminaData.split_run_name(
                self.run_name)
        # Run ID
        try:
            self.run_id = self.metadata.run_id
        except AttributeError:
            self.run_id = None
        if self.run_id is None:
            self.run_id = run_id(
                self.run_name,
                platform=self.metadata.platform,
                facility_run_number=self.metadata.run_number,
                analysis_number=self.metadata.analysis_number)
        # Run reference ID
        try:
            self.run_reference_id = self.metadata.run_reference_id
        except AttributeError:
            self.run_reference_id = None
        if self.run_reference_id is None:
            self.run_reference_id = run_reference_id(
                self.run_id,
                flow_cell_mode=self.metadata.flow_cell_mode)
        # Look for outputs from bclToFastq and analysis projects
        logger.debug("Examining subdirectories of %s" % self._analysis_dir)
        for dirn in bcf_utils.list_dirs(self._analysis_dir):
            # Look for sequencing data
            try:
                data = IlluminaData.IlluminaData(self._analysis_dir,
                                                 unaligned_dir=dirn)
                logger.debug("- %s: sequencing data" % dirn)
                self._bcl2fastq_dirs.append(dirn)
                self.sequencing_data.append(data)
                continue
            except IlluminaData.IlluminaDataError:
                pass
            except Exception as ex:
                logger.warning("Exception when attempting to load "
                               "subdir '%s' as CASAVA/bcl2fastq output "
                               "(ignored): %s" % (dirn,ex))
            # Look for analysis data
            data = AnalysisProject(dirn,os.path.join(self._analysis_dir,dirn))
            if data.is_analysis_dir:
                if dirn == 'undetermined':
                    logger.debug("- %s: undetermined indexes" % dirn)
                    self.undetermined = data
                else:
                    # Check against projects.info, if possible
                    try:
                        if dirn not in self.projects_metadata:
                            logger.debug("- %s: not in projects.info" % dirn)
                            self._extra_dirs.append(dirn)
                            continue
                    except AttributeError:
                        pass
                    logger.debug("- %s: project directory" % dirn)
                    self._project_dirs.append(dirn)
                    self.projects.append(data)
                continue
            else:
                # Unidentified contents
                self._extra_dirs.append(dirn)
                logger.debug("- %s: unknown" % dirn)

    @property
    def analysis_dir(self):
        """
        Return the path to the analysis directory
        """
        return self._analysis_dir

    @property
    def n_projects(self):
        """
        Return number of projects found
        """
        return len(self.projects)

    @property
    def n_sequencing_data(self):
        """
        Return number of sequencing data dirs found
        """
        return len(self.sequencing_data)

    @property
    def paired_end(self):
        """
        Return True if run is paired end, False if single end
        """
        return reduce(lambda x,y: x and y.info.paired_end,self.projects,True)

    def get_projects(self,pattern=None,include_undetermined=True):
        """
        Return the analysis projects in a list

        By default returns all projects within the analysis
        
        If the 'pattern' is not None then it should be a simple pattern
        used to match against available names to select a subset of
        projects (see bcf_utils.name_matches).

        If 'include_undetermined' is True then the undetermined
        project will also be included; otherwise it will be omitted.

        Arguments:
          pattern (str): (optional) glob-style pattern to match
            project names against
          include_undetermined (bool): if True (the default) then
            include the 'Undetermined' project

        Returns:
          List: list of AnalysisProject instances.
        """
        projects = [p for p in self.projects]
        if include_undetermined and self.undetermined:
            projects.append(self.undetermined)
        # Filter on pattern
        if pattern is not None:
            projects = list(filter(lambda p:
                                   fnmatch.fnmatch(p.name,pattern),
                                   projects))
        return projects
        
class AnalysisProject:
    """
    Class describing an analysis project

    Conceptually an analysis project consists of a set of samples
    from a single sequencing experiment, plus associated data e.g.
    QC results.

    Practically an analysis project is represented by a directory
    with a set of fastq files.

    Provides the following properties:

    * name: name of the project
    * dirn: associated directory (full path)
    * fastq_dirs: list of all subdirectories with fastq files (relative
      to dirn)
    * fastq_dir: directory with 'active' fastq file set (full path)
    * fastqs: list of fastq files in fastq_dir
    * read_numbers: list of read numbers from the Fastqs
    * samples: list of AnalysisSample objects generated from fastq_dir
    * multiple_fastqs: True if at least one sample has more than one fastq
      file per read associated with it
    * fastq_format: either 'fastqgz' or 'fastq'

    There is also an 'info' property with the following additional
    properties:

    * run: run name
    * user: user name
    * PI: PI name
    * library_type: library type, either None or e.g. 'RNA-seq' etc
    * single_cell_platform: single cell prep platform, either None or
      '10xGenomics Chromium 3'' etc
    * number of cells: number of cells in single cell projects
    * organism: organism, either None or e.g. 'Human' etc
    * platform: sequencing platform, either None or e.g. 'miseq' etc
    * comments: additional comments, either None or else string of text
    * paired_end: True if data is paired end, False if not
    * primary_fastq_dir: subdirectory holding the 'primary' fastq set
    * sequencer_model: model of sequencer used to generate the data

    It is possible for a project to have multiple sets of associated
    fastq files, held within separate subdirectories of the project
    directory. A list of subdirectory names with fastq sets can be
    accessed via the 'fastq_dirs' property.

    The 'active' fastq set defaults to the 'primary' set (taken from the
    'primary_fastq_dir' info property). An alternative active set can be
    specified using the 'fastq_dir' argument when instantiating the
    AnalysisProject; the active fastq set can also be switched for an
    existing AnalysisProject using the 'use_fastq_dir' method.

    The directory holding the primary fastq set is taken from the
    'fastq_dir' argument of the 'create_project' method when creating
    the project directory (by default this is the 'fastqs' subdirectory
    of the project directory). It can be changed using the
    'set_primary_fastq_dir' method.

    Arguments:
      name (str): name of the project (or path to project
        directory, if 'dirn' not supplied)
      dirn (str): optional, project directory (can be full or
        relative path)
      user (str): optional, specify name of the user
        PI (str): optional, specify name of the principal
        investigator(s)
      library_type (str): optional, specify library type e.g.
        'RNA-seq', 'miRNA' etc
      single_cell_platform (str): optional, specify single cell
        preparation platform e.g. '10xGenomics'
      organism (str): optional, specify organism e.g. 'Human',
        'Mouse' etc (separate multiple organisms with ';',
        use '?' if organism is not known)
      platform (str): optional, specify sequencing platform
        e.g 'miseq'
      run (str): optional, name of the run
      comments (str): optional, free text comments associated
        with the run (separate multiple commenst with ';')
      fastq_attrs (BaseFastqAttrs): optional, specify a class
        to use to get attributes from a Fastq file name (e.g.
      sample name, read number etc). Defaults to
        'AnalysisFastq'.
      fastq_dir (str): optional, explicitly specify the
        subdirectory holding the set of Fastq files to load;
        defaults to 'fastq' (if present) or to the top-level
        of the project directory (if absent).
    """
    def __init__(self,name,dirn=None,user=None,PI=None,library_type=None,
                 single_cell_platform=None,organism=None,run=None,
                 comments=None,platform=None,sequencer_model=None,
                 fastq_attrs=None,fastq_dir=None):
        """
        Create a new AnalysisProject instance
        """
        if dirn is not None:
            self.dirn = os.path.abspath(dirn)
        else:
            self.dirn = os.path.abspath(name)
        self.fastq_dir = None
        self.fastq_dirs = []
        self.fastq_format = None
        self.samples = []
        self._qc_dir = None
        self.info = AnalysisProjectInfo()
        self.info_file = os.path.join(self.dirn,"README.info")
        # Function to use for getting Fastq information
        if fastq_attrs is None:
            self.fastq_attrs = AnalysisFastq
        else:
            self.fastq_attrs = fastq_attrs
        # Populate from the directory contents
        self.populate(fastq_dir=fastq_dir)
        # Set project name
        if dirn is not None:
            # Name was explicitly supplied
            self.name = name
        elif self.info.name is not None:
            # Get name from metadata
            self.name = self.info.name
        else:
            # Default name from directory
            self.name = os.path.basename(self.dirn)
        # (Re)set metadata
        self.info['name'] = self.name
        if run is not None:
            self.info['run'] = run
        if user is not None:
            self.info['user'] = user
        if PI is not None:
            self.info['PI'] = PI
        if library_type is not None:
            self.info['library_type'] = library_type
        if single_cell_platform is not None:
            self.info['single_cell_platform'] = single_cell_platform
        if organism is not None:
            self.info['organism'] = organism
        if platform is not None:
            self.info['platform'] = platform
        if sequencer_model is not None:
            self.info['sequencer_model'] = sequencer_model
        if comments is not None:
            self.info['comments'] = comments

    def populate(self,fastq_dir=None):
        """
        Populate data structure from directory contents

        Arguments:
          fastq_dir (str): (optional) specify the
            subdirectory with Fastq files to use for
            populating the 'AnalysisProject'
        """
        if not os.path.exists(self.dirn):
            # Nothing to do, yet
            return
        # Get data from info file, if present
        if os.path.isfile(self.info_file):
            self.info.load(self.info_file)
        # Identify possible fastq subdirectories
        fastq_dirs = []
        for d in bcf_utils.list_dirs(self.dirn):
            fq_dir = os.path.join(self.dirn,d)
            fastqs = self.find_fastqs(fq_dir)
            if fastqs:
                fastq_dirs.append(d)
        # Also check top-level dir
        if self.find_fastqs(self.dirn):
            fastq_dirs.append('.')
        self.fastq_dirs = fastq_dirs
        logger.debug("Possible fastq dirs: %s" %
                     ','.join(self.fastq_dirs))
        # Set primary fastq file directory
        if not self.fastq_dirs:
            logger.debug("No fastq dirs located for %s" % self.dirn)
            return
        if self.info.primary_fastq_dir is None:
            if 'fastqs' in self.fastq_dirs:
                self.info['primary_fastq_dir'] = 'fastqs'
            else:
                self.info['primary_fastq_dir'] = self.fastq_dirs[0]
        if fastq_dir is None:
            fastq_dir = self.info.primary_fastq_dir
        else:
            if fastq_dir.startswith("%s%s" % (self.dirn,os.sep)):
                fastq_dir_ = os.path.relpath(fastq_dir,self.dirn)
            else:
                fastq_dir_ = fastq_dir
            if fastq_dir_ not in self.fastq_dirs:
                logger.warning("Requested fastqs dir '%s' not in list "
                               "of possible dirs %s" %
                               (fastq_dir,
                                ', '.join(self.fastq_dirs)))
        self.fastq_dir = os.path.normpath(
            os.path.join(self.dirn,fastq_dir))
        # Collect fastq files
        fastqs = self.find_fastqs(self.fastq_dir)
        if fastqs:
            self.fastq_format = self.determine_fastq_format(fastqs[0])
        logger.debug("Assigning fastqs to samples...")
        self.samples = []
        for fq in fastqs:
            name = self.fastq_attrs(fq).sample_name
            try:
                sample = self.get_sample(name)
            except KeyError:
                sample = AnalysisSample(name,
                                        fastq_attrs=self.fastq_attrs)
                self.samples.append(sample)
            sample.add_fastq(os.path.normpath(
                os.path.join(self.fastq_dir,fq)))
        # Sort samples by name
        self.samples = sorted(self.samples,
                              key=lambda s: split_sample_name(s.name))
        logger.debug("Listing samples and files:")
        for sample in self.samples:
            logger.debug("* %s: %s" % (sample.name,sample.fastq))
        # Set paired_end flag for project
        self.info['paired_end'] = self.determine_paired_end()
        # Set the QC output dir, if not already set
        if self.qc_dir is None:
            self.use_qc_dir('qc')

    def find_fastqs(self,dirn):
        """
        Return list of Fastq files found in directory

        Arguments:
          dirn (str): path to directory to search

        Returns:
          List: list of Fastq file names.
        """
        logger.debug("Searching '%s' for fastqs" % dirn)
        fastq_tuples = Pipeline.GetFastqGzFiles(dirn)
        if not fastq_tuples:
            logger.debug("No fastq.gz files found")
            fastq_tuples = Pipeline.GetFastqFiles(dirn)
            if not fastq_tuples:
                logger.debug("No fastq files found")
        # Unpack tuples
        # see https://stackoverflow.com/a/952952/579925
        fastqs = [item for sublist in fastq_tuples for item in sublist]
        return fastqs

    def determine_fastq_format(self,fastq):
        """
        Return type for Fastq file ('fastq' or 'fastqgz')

        Arguments:
          fastq (str): path or name of Fastq file

        Returns:
          String: either 'fastqgz' or 'fastq'.
        """
        if fastq.endswith('.gz'):
            return 'fastqgz'
        else:
            return 'fastq'

    def set_primary_fastq_dir(self,new_primary_fastq_dir):
        """
        Update the primary fastq directory for the project

        This sets the primary fastq directory (aka primary
        fastq set) to the specified name, which must be a
        subdirectory of the project directory.

        Updating the primary fastq directory also causes
        the 'samples' metadata item for the project to be
        updated.

        Relative paths are assumed to be subdirectories
        of the project directory.

        Note that it doesn't change the active fastq set;
        use the 'use_fastq_dir' method to do this.

        Arguments:
          new_primary_fastq_dir (str): path to the
            (sub)directory to be treated as the primary
            Fastq directory for the project

        Raises:
          Exception: if specified directory doesn't exist.
        """
        if not os.path.isabs(new_primary_fastq_dir):
            full_fastq_dir = os.path.join(self.dirn,
                                          new_primary_fastq_dir)
        else:
            full_fastq_dir = new_primary_fastq_dir
        if full_fastq_dir in [os.path.join(self.dirn,d)
                              for d in self.fastq_dirs]:
            self.info['primary_fastq_dir'] = new_primary_fastq_dir
            self.info['samples'] = self.sample_summary()
            self.info.save(self.info_file)
        else:
            raise Exception("Can't update primary fastq dir to '%s' "
                            "for project '%s': directory doesn't exist"
                            % (new_primary_fastq_dir,self.name))

    def use_fastq_dir(self,fastq_dir=None,strict=True):
        """
        Switch fastq directory and repopulate

        Switch to a specified source fastq dir, or to the
        primary fastq dir if none is supplied.

        Relative paths are assumed to be subdirectories
        of the project directory.

        Arguments:
          fastq_dir (str): path to the fastq dir to switch
            to; must be a subdirectory of the project,
            otherwise an exception is raised unless 'strict'
            is set to False
          strict (bool): if True (the default) then
            'fastq_dir' must resolve to a subdirectory of
            the project; otherwise an exception is raised.
            Setting 'strict' to False allows the fastq dir
            to be outside of the project

        Raises:
          Exception: if specified directory is not a
            subdirectory of the project.
        """
        if fastq_dir is None:
            fastq_dir = self.info.primary_fastq_dir
        if fastq_dir is None:
            if 'fastqs' in self.fastq_dirs:
                fastq_dir = 'fastqs'
            elif '.' in self.fastq_dirs:
                fastq_dir = '.'
            else:
                fastq_dir = self.fastq_dirs[0]
        if not os.path.isabs(fastq_dir):
            full_fastq_dir = os.path.join(self.dirn,fastq_dir)
        else:
            full_fastq_dir = fastq_dir
        full_fastq_dir = os.path.normpath(full_fastq_dir)
        if full_fastq_dir not in [os.path.normpath(os.path.join(self.dirn,d))
                                  for d in self.fastq_dirs]:
            msg = "Fastq dir '%s' not found in project '%s' (%s)" % \
                  (fastq_dir,self.name,self.dirn)
            if strict:
                raise Exception(msg)
            else:
                logger.warning(msg)
        self.populate(fastq_dir=fastq_dir)

    def setup_qc_dir(self,qc_dir=None,fastq_dir=None):
        """
        Set up a QC outputs directory

        Creates a QC outputs directory with a metadata
        file 'qc.info'.

        Arguments:
          qc_dir (str): path to QC outputs directory
            to set up. If a relative path is supplied then
            is assumed to be relative to the analysis
            project directory. If 'None' then defaults to
            the current 'qc_dir' for the project.
          fastq_dir (str): set the associated source Fastq
            directory (optional). If 'None' then defaults
            to the previously associated fastq_dir for the
            QC dir (or the current 'fastq_dir' for the
            project if that isn't set).

        Returns:
          String: full path to the QC directory.

        Raises:
          Exception: if previously stored Fastq source dir
            doesn't match the one supplied via 'fastq_dir'.
        """
        print("Setting up QC directory")
        if qc_dir is None:
            qc_dir = os.path.relpath(self.qc_dir,self.dirn)
            print("Assuming default QC dir: %s" % qc_dir)
        if not os.path.isabs(qc_dir):
            qc_dir = os.path.join(self.dirn,qc_dir)
        if not os.path.exists(qc_dir):
            print("Creating QC dir: %s" % qc_dir)
            bcf_utils.mkdir(qc_dir,mode=0o775)
        else:
            print("QC dir already exists: %s" % qc_dir)
        # Set up metadata
        qc_info = self.qc_info(qc_dir)
        print("qc_dir            : %s" % qc_dir)
        print("Supplied fastq_dir: %s" % fastq_dir)
        print("Stored fastq_dir  : %s" % qc_info.fastq_dir)
        if fastq_dir is None:
            if qc_info.fastq_dir is not None:
                fastq_dir = qc_info.fastq_dir
                print("Using stored Fastq dir for this QC dir")
            else:
                fastq_dir = os.path.relpath(self.fastq_dir,self.dirn)
                print("Assuming default Fastq dir: %s" % fastq_dir)
        if qc_info.fastq_dir is not None:
            if qc_info.fastq_dir != fastq_dir:
                raise Exception("Project '%s': supplied Fastq dir ('%s') "
                                "differs from stored dir ('%s') for QC "
                                "dir '%s'" % (self.name,
                                              fastq_dir,
                                              qc_info.fastq_dir,
                                              qc_dir))
        print("Setting associated Fastq dir: %s" % fastq_dir)
        qc_info['fastq_dir'] = fastq_dir
        qc_info.save()
        # Return the path to the QC directory
        return qc_dir

    def use_qc_dir(self,qc_dir):
        """
        Switch the default QC outputs directory

        Arguments:
          qc_dir (str): path to new default QC outputs
            directory. If a relative path is supplied then
            is assumed to be relative to the analysis
            project directory.
        """
        self._qc_dir = qc_dir
        if not os.path.isabs(self._qc_dir):
            self._qc_dir = os.path.join(self.dirn,
                                        self._qc_dir)

    def qc_info(self,qc_dir):
        """
        Fetch the metadata object for with QC dir

        Arguments:
          qc_dir (str): path to QC outputs directory

        Returns:
          AnalysisProjectQCDirInfo: metadata object
            with the metadata for the QC directory.
        """
        if not os.path.isabs(qc_dir):
            qc_dir = os.path.join(self.dirn,qc_dir)
        qc_info = os.path.join(qc_dir,"qc.info")
        return AnalysisProjectQCDirInfo(filen=qc_info)

    def create_directory(self,illumina_project=None,fastqs=None,
                         fastq_dir=None,
                         short_fastq_names=False,
                         link_to_fastqs=False):
        """Create and populate analysis directory for an IlluminaProject

        Creates a new directory corresponding to the AnalysisProject
        object, and optionally also populates with links to FASTQ files
        from a supplied IlluminaProject object.

        The directory structure it creates is:

        ::

            dir/
               fastqs/
               logs/
               ScriptCode/

        It also creates an info file with metadata about the project.

        Arguments:
          illumina_project (IlluminaProject): (optional) populated
            IlluminaProject object from which the analysis directory
            will be populated
          fastqs (list): (optional) list of Fastq files to import
          fastq_dir (str): (optional) name of subdirectory to put
            Fastq files into; defaults to 'fastqs'
          short_fastq_names (bool): (optional) if True then transform
            Fastq file names to be the shortest possible unique names;
            if False (default) then use the original Fastq names
          link_to_fastqs (bool): (optional) if True then make symbolic
            links to the Fastq files; if False (default) then make hard
            links
        """
        logger.debug("Creating analysis directory for project '%s'" % self.name)
        # Check for & create directory
        if os.path.exists(self.dirn):
            logger.warning("Directory %s already exists" % self.dirn)
        else:
            logger.debug("Making analysis directory %s" % self.dirn)
            bcf_utils.mkdir(self.dirn,mode=0o775)
        # Make a 'ScriptCode' directory
        scriptcode_dir = os.path.join(self.dirn,"ScriptCode")
        bcf_utils.mkdir(scriptcode_dir,mode=0o775)
        # Put a file in ScriptCode to make sure it's
        # not pruned on subsequent rsync operations
        fp = open(os.path.join(self.dirn,'ScriptCode','README.txt'),'w')
        fp.write("The ScriptCode directory is a place to put custom scripts and programs")
        fp.close()
        # Make a 'fastqs' directory
        if fastq_dir is None:
            fastq_dir = "fastqs"
        fastq_dir = os.path.join(self.dirn,fastq_dir)
        bcf_utils.mkdir(fastq_dir,mode=0o775)
        # Check for & create links to fastq files
        if fastqs is None:
            # Make a list of fastqs to import from the supplied
            # IlluminaProject object
            fastqs = []
            if illumina_project is not None:
                for sample in illumina_project.samples:
                    for fastq in sample.fastq:
                        fastqs.append(os.path.join(sample.dirn,fastq))
        if short_fastq_names:
            # Get mapping to (shortened) unique names
            fastq_names = IlluminaData.get_unique_fastq_names(fastqs)
        else:
            # Use full names
            fastq_names = {}
            for fq in fastqs:
                fastq_names[fq] = os.path.basename(fq)
        for fastq in fastqs:
            target_fq = os.path.join(fastq_dir,fastq_names[fastq])
            if os.path.exists(target_fq):
                logger.warning("Target '%s' already exists" % target_fq)
            else:
                if link_to_fastqs:
                    logger.debug("Making symlink to %s" % fastq)
                    bcf_utils.mklink(fastq,target_fq,relative=True)
                else:
                    logger.debug("Making hard link to %s" % fastq)
                    os.link(fastq,target_fq)
        # Populate
        self.populate(fastq_dir=os.path.basename(fastq_dir))
        # Update metadata: primary fastq dir
        self.info['primary_fastq_dir'] = os.path.relpath(fastq_dir,
                                                         self.dirn)
        # Update metadata: sample summary
        self.info['samples'] = self.sample_summary()
        # Save metadata
        self.info.save(self.info_file)

    def sample_summary(self):
        """
        Generate a summary of the sample names

        Generates a description string which summarises
        the number and names of samples in the project.

        The description is of the form:

        ::

            2 samples (PJB1, PJB2)

        Returns:
          String: summary of sample names.
        """
        # Get Fastqs
        fq_dir = os.path.join(self.dirn,
                              self.info.primary_fastq_dir)
        fastqs = self.find_fastqs(fq_dir)
        # Assign Fastqs to sample names
        samples = dict()
        for fq in fastqs:
            fq = self.fastq_attrs(fq)
            name = fq.sample_name
            try:
                samples[name].append(fq)
            except KeyError:
                samples[name] = [fq,]
        # Reduce Fastqs for each sample to the minimum set
        multiple_fastqs = False
        for sample_name in samples:
            fastqs = samples[sample_name]
            reduced_fastqs = [fastqs[0],]
            for fq in fastqs[1:]:
                for attr in ('sample_number',
                             'barcode_sequence',
                             'lane_number',
                             'set_number'):
                    if getattr(reduced_fastqs[0],attr) != \
                       getattr(fq,attr):
                        reduced_fastqs.append(fq)
                        break
            # Remove index reads
            reduced_fastqs = list(filter(lambda fq:
                                         not fq.is_index_read,
                                         reduced_fastqs))
            if len(reduced_fastqs) > 1:
                multiple_fastqs = True
            samples[sample_name] = reduced_fastqs
        # Generate description
        sample_names = sort_sample_names(samples.keys())
        n_samples = len(sample_names)
        if n_samples == 0:
            samples_description = "No samples"
        else:
            samples_description = "%s %s" % \
                                 (n_samples,
                                  'sample' if n_samples == 1 else 'samples')
            samples_description += " (%s" % ', '.join(sample_names)
            if multiple_fastqs:
                samples_description += ", multiple fastqs per sample"
            samples_description += ")"
        return samples_description

    @property
    def exists(self):
        """
        Check if analysis project directory already exists
        """
        return os.path.exists(self.dirn)

    @property
    def is_analysis_dir(self):
        """
        Determine if directory really is an analysis project

        This is a strict test:

        - the project must contain Fastqs
        - the project must contain a valid metadata file
        """
        if not self.fastqs:
            # No Fastqs
            return False
        if not os.path.exists(self.info_file):
            # No metadata file
            return False
        try:
            AnalysisProjectInfo().load(self.info_file,
                                       fail_on_error=True)
        except Exception:
            # Failed to load valid metadata file
            return False
        # All tests passed
        return True

    @property
    def qc_dir(self):
        """
        Return path to default QC outputs directory
        """
        return self._qc_dir

    @property
    def qc_dirs(self):
        """
        List QC output directories
        """
        qc_dirs = []
        for d in bcf_utils.list_dirs(self.dirn):
            if d.startswith("qc") and \
               os.path.exists(os.path.join(self.dirn,d,'qc.info')):
                qc_dirs.append(d)
        return qc_dirs

    @property
    def multiple_fastqs(self):
        """
        Determine if there are multiple Fastqs per sample
        """
        if not len(self.samples):
            return False
        else:
            return reduce(lambda x,y: x and y,
                          [len(s.fastq_subset(read_number=1)) > 1 for s in self.samples])

    @property
    def fastqs(self):
        """
        Return a list of Fastqs
        """
        fastqs = []
        for s in self.samples:
            fastqs.extend(s.fastq)
        return fastqs

    @property
    def read_numbers(self):
        """
        List the read numbers from the Fastqs

        Returns:
          List: a list of integer read numbers from the
            associated Fastqs; for example if there are
            R1 and R2 reads then [1, 2] will be returned.
        """
        read_numbers = set()
        for s in self.samples:
            for n in s.read_numbers:
                read_numbers.add(n)
        return sorted(list(read_numbers))

    @property
    def fastqs_are_symlinks(self):
        """
        Return True if Fastq files are symbolic links, False if not
        """
        for s in self.samples:
            if s.fastqs_are_symlinks:
                return True
        return False

    def get_sample(self,name):
        """
        Return sample that matches 'name'

        Arguments:
          name (str): name of a sample

        Returns:
          AnalysisSample: sample object with the matching name

        Raises
          KeyError: if no match is found.
        """
        for sample in self.samples:
            if sample.name == name: return sample
        raise KeyError("No matching sample for '%s'" % name)

    def get_samples(self,pattern):
        """
        Return list of sample matching pattern

        Arguments:
          pattern (str): simple 'glob' style pattern

        Returns:
          List: list of samples with names matching the supplied
            pattern (or an empty list if no names match).
        """
        samples = []
        for sample in self.samples:
            if bcf_utils.name_matches(sample.name,pattern):
                samples.append(sample)
        return samples

    def prettyPrintSamples(self):
        """
        Return a nicely formatted string describing the sample names

        Wraps a call to 'pretty_print_names' function.

        Returns:
          String: pretty description of sample names.
        """
        return bcf_utils.pretty_print_names(self.samples)

    def determine_paired_end(self):
        """
        Return whether or not project has paired end samples
        """
        paired_end = True
        for sample in self.samples:
            paired_end = (paired_end and sample.paired_end)
        return paired_end

    def __repr__(self):
        return "AnalysisProject(%s)" % self.name

class AnalysisSample:
    """
    Class describing an analysis sample

    An analysis sample consists of a set of Fastqs files
    corresponding to a single sample.

    AnalysisSample has the following properties:

    * name: name of the sample
    * fastq: list of Fastq files associated with the sample
    * paired_end: True if sample is paired end, False if not
    * read_numbers: list of read numbers from the Fastqs
    * fastqs_are_symlinks: True if associated Fastqs are symlinks

    Note that the 'fastq' list will include any index read fastqs
    (i.e. I1/I2) as well as R1/R2 fastqs.

    It also provides the following methods:

    * add_fastq: associates a Fastq file with the sample
    * fastq_subset: returns a subset of associated Fastqs

    Arguments:
      name (str): sample name
      fastq_attrs (BaseFastqAttrs): optional, specify a
        class to use to get attributes from a Fastq file name
        (e.g. sample name, read number etc). Defaults to
        'AnalysisFastq'.
    """
    def __init__(self,name,fastq_attrs=None):
        """
        Create a new AnalysisSample instance
        """
        self.name = name
        self.fastq = []
        self.paired_end = False
        # Function to use for getting Fastq information
        if fastq_attrs is None:
            self.fastq_attrs = AnalysisFastq
        else:
            self.fastq_attrs = fastq_attrs

    def add_fastq(self,fastq):
        """
        Add a reference to a Fastq file in the sample

        Arguments:
          fastq (str): full path for the Fastq file
        """
        assert(os.path.isabs(fastq))
        self.fastq.append(fastq)
        # Sort fastq's into order
        self.fastq.sort()
        # Check paired-end status
        if not self.paired_end:
            fq = self.fastq_attrs(fastq)
            if fq.read_number == 2:
                self.paired_end = True

    def fastq_subset(self,read_number=None):
        """
        Return a subset of Fastq files from the sample

        Note that only R1/R2 files will be returned; index read
        fastqs (i.e. I1/I2) are excluded regardless of read number.

        Arguments:
          read_number (int): select subset based on read_number
            (1 or 2)

        Returns:
          List: list of full paths to Fastq files matching the
            selection criteria.
        """
        # Build list of fastqs that match the selection criteria
        fastqs = []
        for fastq in self.fastq:
            fq = self.fastq_attrs(fastq)
            if fq.is_index_read:
                logger.debug("Rejecting index read %s" % fastq)
                continue
            if fq.read_number is None:
                logger.debug("Unable to determine read number for %s,"
                             "assume R1" % fastq)
                fq_read_number = 1
            else:
                fq_read_number = fq.read_number
            if fq_read_number == read_number:
                fastqs.append(fastq)
        # Sort into dictionary order and return
        fastqs.sort()
        return fastqs

    @property
    def read_numbers(self):
        """
        List the read numbers from the Fastqs

        Returns:
          List: a list of integer read numbers from the
            associated Fastqs; for example if there are
            R1 and R2 reads then [1, 2] will be returned.
        """
        read_numbers = set()
        for fastq in self.fastq:
            fq = self.fastq_attrs(fastq)
            if fq.is_index_read:
                continue
            if fq.read_number is None:
                read_numbers.add(1)
            else:
                read_numbers.add(fq.read_number)
        return sorted(list(read_numbers))

    @property
    def fastqs_are_symlinks(self):
        """
        Return True if Fastq files are symlinked, False if not
        """
        for fastq in self.fastq:
            if os.path.islink(fastq):
                return True
        return False

    def __repr__(self):
        """
        Implement __repr__ built-in

        Returns:
          String: representation for the sample (i.e. the
            sample name).
        """
        return str(self.name)

#######################################################################
# Functions
#######################################################################

def run_id(run_name,platform=None,facility_run_number=None,
           analysis_number=None):
    """
    Return a run ID e.g. 'HISEQ_140701/242#22'

    The run ID is a code that identifies the sequencing run, and has
    the general form:

    ``PLATFORM_DATESTAMP[/INSTRUMENT_RUN_NUMBER]#FACILITY_RUN_NUMBER[.ANALYSIS_NUMBER]``

    - PLATFORM is always uppercased e.g. HISEQ, MISEQ, GA2X
    - DATESTAMP is the YYMMDD code e.g. 140701
    - INSTRUMENT_RUN_NUMBER is the run number that forms part of the
      run directory e.g. for '140701_SN0123_0045_000000000-A1BCD'
      it is '45'
    - FACILITY_RUN_NUMBER is the run number that has been assigned
      by the facility
    - ANALYSIS_NUMBER is an optional number assigned to the analysis
      to distinguish it from other analysis attempts (for example, if
      a run is reprocessed at a later date with updated software)

    Note that the instrument run number is only used if it differs
    from the facility run number.

    If the platform isn't supplied then the instrument name is
    used instead, e.g.:

    ``SN0123_140701/242#22``

    If the run name can't be split into components then the
    general form will be:

    ``[PLATFORM_]RUN_NAME[#FACILITY_RUN_NUMBER]``

    depending on whether platform and/or facility run number have
    been supplied. For example for a run called 'rag_05_2017':

    ``MISEQ_rag_05_2017#90``

    Arguments:
      run_name (str): the run name (can be a path)
      platform (str): the platform name (optional)
      facility_run_number (int): the run number assigned by the
        local facility (can be different from the instrument
        run number) (optional)
      analysis_number (int): number assigned to this analysis
        to distinguish it from other analysis attempts
        (optional)

    Returns:
      String: run ID.
    """
    # Extract information from run name
    run_name = os.path.basename(os.path.normpath(run_name))
    try:
        datestamp,instrument,run_number = IlluminaData.split_run_name(run_name)
    except Exception as ex:
        logger.warning("Unable to extract information from run name '%s'" \
                       % run_name)
        logger.warning("Exception: %s" % ex)
        instrument = None
        date_stamp = None
        run_number = None
    # Platform
    if platform is not None:
        platform = platform.upper()
    else:
        platform = instrument
    # Run number
    if run_number is not None:
        run_number = int(run_number)
    # Facility run number
    if facility_run_number is not None:
        try:
            facility_run_number = int(facility_run_number)
        except ValueError:
            facility_run_number = None
    else:
        facility_run_number = None
    # Construct the run id
    if platform is not None:
        run_id = platform
        if datestamp is not None:
            run_id += "_%s" % datestamp
        else:
            run_id += "_%s" % run_name
    else:
        run_id = run_name
    if run_number is not None:
        try:
            if run_number != facility_run_number:
                run_id += "/%s" % run_number
        except ValueError:
            run_id += "/%s" % run_number
    if facility_run_number is not None:
        run_id += "#%s" % facility_run_number
    if analysis_number is not None:
        run_id += ".%s" % analysis_number
    return run_id

def run_reference_id(run_id,flow_cell_mode=None):
    """
    Return a run reference ID (e.g. 'NOVASEQ6000_230419/74#22_SP')

    The reference is constructed from the run ID
    pluse the following additional items, if defined:

    - flow cell mode
    """
    return "%s%s" % (run_id,
                     '' if not flow_cell_mode
                     else "_%s" % flow_cell_mode)

def split_sample_name(s):
    """
    Split sample name into numerical and non-numerical parts

    Utility function which splits the supplied sample name
    into numerical (i.e. integer) and non-numerical (i.e. all
    other types of character) parts, and returns the parts as
    a list.

    For example:

    >>> split_sample_name("PJB_01-123_T004")
    ['PJB_',1,'-',123,'_T',4]

    Arguments:
      s (str): the sample name to be split

    Returns:
      List: list with the numerical and non-numerical parts
        of the name.
    """
    parts = []
    current = ''
    for c in str(s):
        if (current.isdigit() and not c.isdigit()) or \
           (not current.isdigit() and c.isdigit()):
            parts.append(current)
            current = ''
        current += c
    if current:
        parts.append(current)
    for i,part in enumerate(parts):
        if part.isdigit():
            parts[i] = int(part)
    return parts

def split_sample_reference(s):
    """
    Split a '[RUN][:PROJECT[/SAMPLE]]' reference id

    Decomposes a reference id of the form:

    [RUN][:PROJECT[/SAMPLE]]

    where:

    - RUN is a run identifier (either a run name e.g.
      '201027_SN00284_0000161_AHXXJHJH', a reference id e.g.
      'HISEQ_201027/161#122', or a path to an analysis
      directory)
    - PROJECT is the name of a project within the run, and
    - SAMPLE is a sample name.

    A subset of elements can be present, in which case
    the missing components will be returned as None.

    Arguments:
      s (str): sample reference identifier

    Returns:
      Tuple: tuple of the form (RUN,PROJECT,SAMPLE)
        extracted from the supplied reference id;
        missing elements are set to None.
    """
    # Try to extract run name first
    if ':' in s:
        run,name = s.split(':')
    else:
        # Either a run or a sample specification
        if '#' in s:
            # Definitely a run reference id
            # Note: this test is quite weak!
            run = s
            name = None
        else:
            # Assume it's PROJECT[/SAMPLE]
            run = None
            name = s
    # Extract project and sample names
    if name:
        if '/' in name:
            # Split PROJECT/SAMPLE
            project,sample = name.split('/')
        elif run is None:
            # There was no leading RUN
            # so assume this is just SAMPLE
            project = None
            sample = name
        else:
            # There was a leading RUN
            # so assume this is just PROJECT
            project = name
            sample = None
    else:
        # RUN but no PROJECT and/or SAMPLE
        project = None
        sample = None
    # Normalise outputs: make sure empty values
    # are replaced with None
    run = run if run else None
    project = project if project else None
    sample = sample if sample else None
    return (run,project,sample)

def match_run_id(run,d):
    """
    Check if a directory matches a run identifier

    The identifier can be any one of:

    - a path to an analysis directory (e.g.
      '/path/to/201029_SN01234_0000123_AHXXXX_analysis')
    - the name of an analysis directory (e.g.
      '201029_SN01234_0000123_AHXXXX_analysis')
    - the name of a sequencing run associated with an
      analysis directory (e.g. '201029_SN01234_0000123_AHXXXX')
    - a run identifier (e.g. 'HISEQ_201029#123')

    If the identifier is a wildcard ('*') then any valid
    analysis directory will be a match.

    Arguments:
      run (str): run identifier
      d (str): path to a run to check against
        the supplied identifier

    Returns:
      Boolean: True if directory matches run ID, False
        if not.
    """
    # Check if run specification is actually a path
    logger.debug("%s: checking full path" % d)
    if d == run:
        return True
    # Check if name matches
    if os.path.basename(d) == run:
        return True
    # Attempt to load as an analysis dir
    try:
        analysis_dir = AnalysisDir(d)
        # Check run name
        logger.debug("%s: run name = %s" % (d,analysis_dir.run_name))
        if analysis_dir.run_name == run or run == '*':
            return d
        # Check run ID
        run_id_ = run_id(
            analysis_dir.run_name,
            platform=analysis_dir.metadata.platform,
            facility_run_number=analysis_dir.metadata.run_number,
            analysis_number=analysis_dir.metadata.analysis_number)
        logger.debug("%s: run ID = %s" % (d,run_id_))
        if run_id_ == run:
            return True
    except Exception as ex:
        # Not an analysis directory
        logger.debug("%s: exception: %s" % (d,ex))
        pass
    return False

def locate_run(run,start_dir=None,ascend=False):
    """
    Locate an analysis directory

    Searches the file system to locate an analysis
    directory which matches the supplied run
    identifier.

    The identifier can be any one of:

    - a path to an analysis directory (e.g.
      '/path/to/201029_SN01234_0000123_AHXXXX_analysis')
    - the name of an analysis directory (e.g.
      '201029_SN01234_0000123_AHXXXX_analysis')
    - the name of a sequencing run associated with an
      analysis directory (e.g. '201029_SN01234_0000123_AHXXXX')
    - a run identifier (e.g. 'HISEQ_201029#123')

    If the run identifier is a wildcard ('*') then
    the first valid analysis directory that is
    encountered on the search path will be matched;
    note that this may result in non-deterministic
    behaviour.

    Arguments:
      run (str): identifier for the run to locate
      start_dir (str): optional path to start
        searching from (defaults to the current
        directory)
      ascend (bool): if True then search by
        ascending into parent directories of
        'start_dir' (default is to search by
        descending into its subdirectories)

    Returns:
      String: path to the analysis directory, or
        None if the specified analysis directory can't
        be located.
    """
    # Set starting directory
    if start_dir is None:
        d = os.getcwd()
    else:
        d = os.path.abspath(start_dir)
    # Check if 'run' is actually a path
    logger.debug("%s: checking full path" % d)
    if match_run_id(run,d):
        return d
    # Current directory doesn't match
    if not ascend:
        # Search the subdirectories recursively
        for subdir in [os.path.join(d,sd) for sd in os.listdir(d)
                       if os.path.isdir(os.path.join(d,sd))]:
            logger.debug("Going into %s" % subdir)
            d = locate_run(run,start_dir=subdir)
            if d:
                # Found in subdirectory
                return d
    else:
        # Check the subdirectories of the current directory
        for subdir in [os.path.join(d,sd) for sd in os.listdir(d)
                       if os.path.isdir(os.path.join(d,sd))]:
            logger.debug("Checking %s" % subdir)
            if match_run_id(run,subdir):
                # Found in subdirectory
                return subdir
        # Move up a level and try again
        if d != os.sep:
            parent_dir = os.path.dirname(d)
            logger.debug("%s: ascending to %s" % (d,parent_dir))
            d = locate_run(run,start_dir=parent_dir,ascend=True)
            if d:
                # Found in parent
                return d
        else:
            # Reached filesystem root
            logger.debug("Reached filesystem root")
    # Run not found
    return None

def locate_project(project_id,start_dir=None,ascend=False):
    """
    Locate an analysis project

    Searches the file system to locate an analysis
    project which matches the supplied project
    identifier.

    The identifier can be either:

    - a path to an analysis project (e.g.
      '/path/to/201029_SN01234_0000123_AHXXXX_analysis/AB'),
      or
    - a valid project identifier (e.g.
      '201029_SN01234_0000123_AHXXXX_analysis:AB',
      'HISEQ_201029#123:AB' etc)
    - a project name

    Arguments:
      project_id (str): identifier for the project
        to locate
      start_dir (str): optional path to start
        searching from (defaults to the current
        directory)
      ascend (bool): if True then search by
        ascending into parent directories of
        'start_dir' (default is to search by
        descending into its subdirectories)

    Returns:
      AnalysisProject: path to the analysis project,
        or None if the specified project can't be
        located.
    """
    # Set starting directory
    if start_dir is None:
        d = os.getcwd()
    else:
        d = os.path.abspath(start_dir)
    # Check if 'project' is actually a path
    project = os.path.join(d,project_id)
    logger.debug("%s: checking full path" % project)
    if os.path.exists(project):
        return AnalysisProject(project)
    # Assume it's a run/project/sample specification
    run,project,sample = split_sample_reference(project_id)
    # No run identifier supplied
    if run is None:
        run = '*'
        # Also: if a project name was supplied
        # without a sample then it will have been
        # assigned to 'sample'
        if sample and project is None:
            project = sample
            sample = None
    # Locate the run
    run = locate_run(run,start_dir=start_dir,ascend=ascend)
    if not run:
        return None
    # Look for the project
    try:
        for p in AnalysisDir(run).projects:
            if p.name == project:
                return p
    except Exception:
        d = os.path.join(run,project)
        if os.path.exists(d):
            return AnalysisProject(d)
    return None

def locate_project_info_file(start_dir):
    """
    Locate project metadata file

    Searches the current directory and its parents for a
    project metadata file ('README.info'), ascending up
    directory levels until either a valid metadata file is
    found, or the root of the filesystem is reached.

    Arguments:
      start_dir (str): path of directory to start
        searching from

    Returns:
      String: path to the metadata file, or 'None' if
        no file can be located.
    """
    d = os.path.abspath(start_dir)
    while True:
        info_file = os.path.join(d,"README.info")
        if os.path.exists(info_file):
            try:
                # Try to load metadata
                AnalysisProjectInfo().load(info_file,
                                           fail_on_error=True)
                return info_file
            except Exception as ex:
                # Failed to load valid metadata file
                pass
        # Try next level up
        d = os.path.dirname(d)
        if d == os.path.sep:
            # Run out of directories
            return None

def copy_analysis_project(project,fastq_dir=None):
    """
    Make a copy of an AnalysisProject instance

    Arguments:
      project (AnalysisProject): project intance to copy
      fastq_dir (str): if set then specifies the Fastq
        subdirectory to use in the new instance

    Returns:
      AnalysisProject: new AnalysisProject instance which
        is a copy of the one supplied on input.
    """
    if fastq_dir is None:
        fastq_dir = project.fastq_dir
    return AnalysisProject(project.name,
                           project.dirn,
                           user=project.info.user,
                           PI=project.info.PI,
                           library_type=project.info.library_type,
                           single_cell_platform=
                           project.info.single_cell_platform,
                           organism=project.info.organism,
                           run=project.info.run,
                           comments=project.info.comments,
                           platform=project.info.platform,
                           sequencer_model=project.info.sequencer_model,
                           fastq_dir=fastq_dir,
                           fastq_attrs=
                           project.fastq_attrs)
