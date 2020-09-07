#!/usr/bin/env python
#
#     analysis: classes & funcs for handling analysis dirs and projects
#     Copyright (C) University of Manchester 2018-2020 Peter Briggs
#
########################################################################
#
# analysis.py
#
#########################################################################

"""
analysis.py

Classes and functions for handling analysis directories and
projects.

Classes:

- AnalysisFastq:
- AnalysisDir:
- AnalysisProject:
- AnalysisSample:

Functions:

- run_reference_id:
- split_sample_name:
- copy_analysis_project:
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

    fastq:            the original fastq file name
    basename:         basename with NGS extensions stripped
    extension:        full extension e.g. '.fastq.gz'
    sample_name:      name of the sample
    sample_number:    integer (or None if no sample number)
    barcode_sequence: barcode sequence (string or None)
    lane_number:      integer (or None if no lane number)
    read_number:      integer (or None if no read number)
    set_number:       integer (or None if no set number)
    is_index_read:    boolean (True if index read, False if not)

    There are three additional attributes:

    format:           string identifying the format of the
                      Fastq name ('Illumina', 'SRA', or None)
    canonical_name:   the 'canonical' part of the name (string,
                      or None if no canonical part could be
                      extracted)
    extras:           the 'extra' part of the name (string,
                      or None if there was no trailing extra
                      part)
    """
    def __init__(self,fastq):
        # Initialise the base class
        BaseFastqAttrs.__init__(self,fastq)
        # Additional attributes
        self.format = None
        self.extras = None
        # Check for SRA format
        fq = SRA_REGEX.match(self.basename)
        if fq is not None:
            fq_attrs = fq.groupdict()
            self.format = "SRA"
            self.sample_name = fq_attrs['sample_name']
            if fq_attrs['read_number'] is not None:
                self.read_number = int(fq_attrs['read_number'])
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
                             "_%s" % self.read_number if self.read_number
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

class AnalysisDir(object):
    """Class describing an analysis directory

    Conceptually an analysis directory maps onto a sequencing run.
    It consists of one or more sets of samples from that run,
    which are represented by subdirectories.

    It is also possible to have one or more subdirectories containing
    outputs from the CASAVA or bclToFastq processing software.

    """
    def __init__(self,analysis_dir):
        """Create a new AnalysisDir instance for a specified directory

        Arguments:
          analysis_dir: name (and path) to analysis directory

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
        try:
            metadata_file = os.path.join(self._analysis_dir,
                                         "metadata.info")
            self.metadata.load(metadata_file)
        except Exception as ex:
            logger.warning("Failed to load metadata file %s: %s" %
                           (metadata_file,ex))
            logger.warning("Attempting to load parameter file")
            try:
                params = AnalysisDirParameters()
                parameter_file = os.path.join(self._analysis_dir,
                                         "auto_process.info")
                params.load(parameter_file,strict=False)
                # Attempt to acquire values from parameters
                for param in ('platform','run_number','source','assay'):
                    if param not in params:
                        print("-- %s: missing" % param)
                        continue
                    print("-- %s: setting to '%s'" % (param,
                                                      params[param]))
                    self.metadata[param] = params[param]
            except Exception as ex:
                # No parameter file either
                logger.warning("Failed to load parameters: %s" % ex)
                logger.warning("Perhaps this is not an auto_process project?")
                raise ex
        # Projects metadata
        try:
            self.projects_metadata = ProjectMetadataFile(
                os.path.join(self._analysis_dir,"projects.info"))
        except Exception as ex:
            logger.warning("Failed to load projects metadata: %s" % ex)
            self.projects_metadata = None
        # Run name
        try:
            self.run_name = self.metadata.run
        except AttributeError:
            self.run_name = self._analysis_dir[0:-len('_analysis')]
        self.run_name = os.path.basename(self.run_name)
        self.date_stamp,\
            self.instrument_name,\
            self.instrument_run_number = IlluminaData.split_run_name(
                self.run_name)
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
                        if not self.projects_metadata.lookup('Project',dirn):
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
        """Return the path to the analysis directory

        """
        return self._analysis_dir

    @property
    def n_projects(self):
        """Return number of projects found

        """
        return len(self.projects)

    @property
    def n_sequencing_data(self):
        """Return number of sequencing data dirs found

        """
        return len(self.sequencing_data)

    @property
    def paired_end(self):
        """Return True if run is paired end, False if single end

        """
        return reduce(lambda x,y: x and y.info.paired_end,self.projects,True)

    def get_projects(self,pattern=None,include_undetermined=True):
        """Return the analysis projects in a list

        By default returns all projects within the analysis
        
        If the 'pattern' is not None then it should be a simple pattern
        used to match against available names to select a subset of
        projects (see bcf_utils.name_matches).

        If 'include_undetermined' is True then the undetermined
        project will also be included; otherwise it will be omitted.

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
        
class AnalysisProject(object):
    """Class describing an analysis project

    Conceptually an analysis project consists of a set of samples
    from a single sequencing experiment, plus associated data e.g.
    QC results.

    Practically an analysis project is represented by a directory
    with a set of fastq files.

    Provides the following properties:

    name        : name of the project
    dirn        : associated directory (full path)
    fastq_dirs  : list of all subdirectories with fastq files (relative
                  to dirn)
    fastq_dir   : directory with 'active' fastq file set (full path)
    fastqs      : list of fastq files in fastq_dir
    samples     : list of AnalysisSample objects generated from fastq_dir
    multiple_fastqs: True if at least one sample has more than one fastq
                     file per read associated with it
    fastq_format: either 'fastqgz' or 'fastq'

    There is also an 'info' property with the following additional
    properties:

    run         : run name
    user        : user name
    PI          : PI name
    library_type: library type, either None or e.g. 'RNA-seq' etc
    single_cell_platform: single cell prep platform, either None or 'ICell8' etc
    number of cells: number of cells in single cell projects
    ICELL8 well list: well list file for ICELL8 single cell projects
    organism    : organism, either None or e.g. 'Human' etc
    platform    : sequencing platform, either None or e.g. 'miseq' etc
    comments    : additional comments, either None or else string of text
    paired_end  : True if data is paired end, False if not
    primary_fastq_dir: subdirectory holding the 'primary' fastq set

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
    """
    def __init__(self,name,dirn=None,user=None,PI=None,library_type=None,
                 single_cell_platform=None,organism=None,run=None,
                 comments=None,platform=None,fastq_attrs=None,
                 fastq_dir=None):
        """Create a new AnalysisProject instance

        Arguments:
          name: name of the project (or path to project directory,
            if 'dirn' not supplied)
          dirn: optional, project directory (can be full or relative
            path)
          user: optional, specify name of the user
          PI: optional, specify name of the principal investigator
          library_type: optional, specify library type e.g. 'RNA-seq',
            'miRNA' etc
          single_cell_platform: optional, specify single cell
            preparation platform e.g. 'Icell8', '10xGenomics' etc
          organism: optional, specify organism e.g. 'Human', 'Mouse'
            etc
          platform: optional, specify sequencing platform e.g 'miseq'
          run: optional, name of the run
          comments: optional, free text comments associated with the
            run
          fastq_attrs: optional, specify a class to use to get
            attributes from a Fastq file name (e.g. sample name, read
            number etc). The supplied class should be a subclass of
            BaseFastqAttrs; defaults to 'AnalysisFastq'.
          fastq_dir: optional, explicitly specify the subdirectory
            holding the set of Fastq files to load; defaults to
            'fastq' (if present) or to the top-level of the project
            directory (if absent).
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
        if comments is not None:
            self.info['comments'] = comments

    def populate(self,fastq_dir=None):
        """Populate data structure from directory contents

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
        paired_end = True
        for sample in self.samples:
            paired_end = (paired_end and sample.paired_end)
        self.info['paired_end'] = paired_end
        # Set the QC output dir, if not already set
        if self.qc_dir is None:
            self.use_qc_dir('qc')

    def find_fastqs(self,dirn):
        """
        Return list of Fastq files found in directory
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

    def use_fastq_dir(self,fastq_dir=None):
        """
        Switch fastq directory and repopulate

        Switch to a specified source fastq dir, or to the
        primary fastq dir if none is supplied.

        Relative paths are assumed to be subdirectories
        of the project directory.
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
            raise Exception("Fastq dir '%s' not found in "
                            "project '%s' (%s)" %
                            (fastq_dir,self.name,self.dirn))
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

        dir/
           fastqs/
           logs/
           ScriptCode/

        It also creates an info file with metadata about the project.

        Arguments:
          illumina_project: (optional) populated IlluminaProject object
            from which the analysis directory will be populated
          fastqs: (optional) list of fastq files to import
          fastq_dir: (optional) name of subdirectory to put fastq files
            into; defaults to 'fastqs'
          short_fastq_names: (optional) if True then transform fastq file
            names to be the shortest possible unique names; if False
            (default) then use the original fastq names
          link_to_fastqs: (optional) if True then make symbolic links to
            to the fastq files; if False (default) then make hard links
    
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
        """Generate a summary of the sample names

        Generates a description string which summarises
        the number and names of samples in the project.

        The description is of the form:

        2 samples (PJB1, PJB2)
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
        sample_names = sorted(samples.keys())
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
        """Check if analysis project directory already exists

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
            if d.startswith("qc"):
                qc_dirs.append(d)
        return qc_dirs

    @property
    def qc(self):
        """
        Return QCReporter instance for QC outputs

        Not supported after version 0.15.0

        Raises exception if invoked
        """
        raise Exception("%s.qc method no longer supported"
                        % self.__class__)

    def qc_report(self,title=None,report_html=None,qc_dir=None,
                  force=False):
        """
        Report QC outputs for project

        Not supported after version 0.15.0

        Use the 'report_qc' function from qc.utils instead.
        """
        raise Exception("%s.qc_report method no longer supported"
                        % self.__class__)

    @property
    def multiple_fastqs(self):
        # Determine if there are multiple fastqs per sample
        if not len(self.samples):
            return False
        else:
            return reduce(lambda x,y: x and y,
                          [len(s.fastq_subset(read_number=1)) > 1 for s in self.samples])

    @property
    def fastqs(self):
        """Return a list of fastqs

        """
        fastqs = []
        for s in self.samples:
            fastqs.extend(s.fastq)
        return fastqs

    @property
    def fastqs_are_symlinks(self):
        """Return True if fastq files are symbolic links, False if not

        """
        for s in self.samples:
            if s.fastqs_are_symlinks:
                return True
        return False

    def verify_qc(self,qc_dir=None):
        """
        Check if QC run has completed successfully

        Not supported after version 0.15.0

        Use the 'verify_qc' function from the qc.utils module
        instead.
        """
        raise Exception("%s.qc_report method no longer supported"
                        % self.__class__)

    def get_sample(self,name):
        """Return sample that matches 'name'

        Arguments:
          name: name of a sample

        Returns:
          AnalysisSample object with the matching name; raises
          KeyError exception if no match is found.

        """
        for sample in self.samples:
            if sample.name == name: return sample
        raise KeyError("No matching sample for '%s'" % name)

    def get_samples(self,pattern):
        """Return list of sample matching pattern

        Arguments:
          pattern: simple 'glob' style pattern

        Returns:
          Python list of samples with names matching the supplied
          pattern (or an empty list if no names match).

        """
        samples = []
        for sample in self.samples:
            if bcf_utils.name_matches(sample.name,pattern):
                samples.append(sample)
        return samples

    def prettyPrintSamples(self):
        """Return a nicely formatted string describing the sample names

        Wraps a call to 'pretty_print_names' function.

        """
        return bcf_utils.pretty_print_names(self.samples)

class AnalysisSample(object):
    """Class describing an analysis sample

    An analysis sample consists of a set of fastqs file corresponding
    to single sample.

    AnalysisSample has the following properties:

    name      : name of the sample
    fastq     : list of fastq files associated with the sample
    paired_end: True if sample is paired end, False if not

    Note that the 'fastq' list will include any index read fastqs
    (i.e. I1/I2) as well as R1/R2 fastqs.

    """

    def __init__(self,name,fastq_attrs=None):
        """Create a new AnalysisSample instance

        Arguments:
          name: sample name
          fastq_attrs: optional, specify a class to use to get
            attributes from a Fastq file name (e.g. sample name, read
            number etc). The supplied class should be a subclass of
            BaseFastqAttrs; defaults to 'AnalysisFastq'.

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
        """Add a reference to a fastq file in the sample

        Arguments:
          fastq: full path for the fastq file

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
        """Return a subset of fastq files from the sample

        Note that only R1/R2 files will be returned; index read
        fastqs (i.e. I1/I2) are excluded regardless of read number.

        Arguments:
          read_number: select subset based on read_number (1 or 2)

        Returns:
          List of full paths to fastq files matching the selection criteria.

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
    def fastqs_are_symlinks(self):
        """Return True if fastq files are symlinked, False if not

        """
        for fastq in self.fastq:
            if os.path.islink(fastq):
                return True
        return False

    def qc_sample(self):
        """Fetch QCSample object for this sample

        Not supported after version 0.15.0

        Raises exception if invoked
        """
        raise Exception("%s.qc_sample method no longer supported"
                        % self.__class__)

    def verify_qc(self,qc_dir,fastq):
        """Check if QC completed for a fastq file

        Not supported after version 0.15.0

        Raises exception if invoked
        """
        raise Exception("%s.verify_qc method no longer supported"
                        % self.__class__)

    def __repr__(self):
        """Implement __repr__ built-in

        Return string representation for the sample -
        i.e. the sample name.

        """
        return str(self.name)

#######################################################################
# Functions
#######################################################################

def run_reference_id(run_name,platform=None,facility_run_number=None):
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

    If the platform isn't supplied then the instrument name is
    used instead, e.g.:

    'SN0123_140701/242#22'

    If the run name can't be split into components then the
    general form will be:

    [PLATFORM_]RUN_NAME[#FACILITY_RUN_NUMBER]

    depending on whether platform and/or facility run number have
    been supplied. For example for a run called 'rag_05_2017':

    'MISEQ_rag_05_2017#90'

    Arguments:
      run_name (str): the run name (can be a path)
      platform (str): the platform name (optional)
      facility_run_number (int): the run number assigned by the
        local facility (can be different from the instrument
        run number) (optional)
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
    # Construct the reference id
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
    return run_id

def split_sample_name(s):
    """Split sample name into numerical and non-numerical parts

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
                           fastq_dir=fastq_dir,
                           fastq_attrs=
                           project.fastq_attrs)
