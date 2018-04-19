#!/usr/bin/env python
#
#     utils: utility classes & funcs for auto_process_ngs module
#     Copyright (C) University of Manchester 2013-2018 Peter Briggs
#
########################################################################
#
# utils.py
#
#########################################################################

__version__ = "0.0.18"

"""utils

Utility classes and functions to support auto_process_ngs module.

Ultimately these should be relocated in the main 'genomics' code
tree at some point.

Classes:

- AnalysisFastq:
- AnalysisDir:
- AnalysisProject:
- AnalysisSample:
- ZipArchive:
- OutputFiles:
- ProgressChecker:

Functions:

- bases_mask_is_paired_end:
- split_user_host_dir:
- get_numbered_subdir:
- find_executables:
- parse_version:
- pretty_print_rows:
- write_script_file:
- edit_file:
- paginate:

"""

#######################################################################
# Imports
#######################################################################

import sys
import os
import fnmatch
import logging
import zipfile
import pydoc
import tempfile
import operator
import applications
import bcftbx.IlluminaData as IlluminaData
import bcftbx.JobRunner as JobRunner
import bcftbx.Pipeline as Pipeline
import bcftbx.utils as bcf_utils
from bcftbx.Md5sum import md5sum
from .metadata import AnalysisDirMetadata
from .metadata import AnalysisDirParameters
from .metadata import AnalysisProjectInfo
from .metadata import ProjectMetadataFile
from .metadata import AnalysisProjectQCDirInfo
from .fastq_utils import IlluminaFastqAttrs
from qc.reporting import QCReporter
from qc.reporting import QCSample
from qc.illumina_qc import expected_qc_outputs
from qc.illumina_qc import check_qc_outputs

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Classes
#######################################################################

class AnalysisFastq(IlluminaFastqAttrs):
    """
    Wrapper for IlluminaFastqAttrs
    """

class AnalysisDir:
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
                        print "-- %s: missing" % param
                        continue
                    print "-- %s: setting to '%s'" % (param,
                                                      params[param])
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
                logging.warning("Exception when attempting to load "
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
            projects = filter(lambda p: fnmatch.fnmatch(p.name,pattern),
                              projects)
        return projects
        
class AnalysisProject:
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
    def __init__(self,name,dirn,user=None,PI=None,library_type=None,
                 single_cell_platform=None,organism=None,run=None,
                 comments=None,platform=None,fastq_attrs=None,
                 fastq_dir=None):
        """Create a new AnalysisProject instance

        Arguments:
          name: name of the project
          dirn: project directory (can be full or relative path)
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
        self.name = name
        self.dirn = os.path.abspath(dirn)
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
        # (Re)set metadata
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
            if fastq_dir not in self.fastq_dirs:
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

        Note that it doesn't change the active fastq set;
        use the 'use_fastq_dir' method to do this.
        """
        if new_primary_fastq_dir in self.fastq_dirs:
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
        primary fastq dir if none is supplied
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
        elif fastq_dir not in self.fastq_dirs:
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
        print "Setting up QC directory"
        if qc_dir is None:
            qc_dir = os.path.relpath(self.qc_dir,self.dirn)
            print "Assuming default QC dir: %s" % qc_dir
        if not os.path.isabs(qc_dir):
            qc_dir = os.path.join(self.dirn,qc_dir)
        if not os.path.exists(qc_dir):
            print "Creating QC dir: %s" % qc_dir
            bcf_utils.mkdir(qc_dir,mode=0775)
        else:
            print "QC dir already exists: %s" % qc_dir
        # Set up metadata
        qc_info = self.qc_info(qc_dir)
        print "qc_dir            : %s" % qc_dir
        print "Supplied fastq_dir: %s" % fastq_dir
        print "Stored fastq_dir  : %s" % qc_info.fastq_dir
        if fastq_dir is None:
            if qc_info.fastq_dir is not None:
                fastq_dir = qc_info.fastq_dir
                print "Using stored Fastq dir for this QC dir"
            else:
                fastq_dir = os.path.relpath(self.fastq_dir,self.dirn)
                print "Assuming default Fastq dir: %s" % fastq_dir
        if qc_info.fastq_dir is not None:
            if qc_info.fastq_dir != fastq_dir:
                raise Exception("Project '%s': supplied Fastq dir ('%s') "
                                "differs from stored dir ('%s') for QC "
                                "dir '%s'" % (self.name,
                                              fastq_dir,
                                              qc_info.fastq_dir,
                                              qc_dir))
        print "Setting associated Fastq dir: %s" % fastq_dir
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
            bcf_utils.mkdir(self.dirn,mode=0775)
        # Make a 'ScriptCode' directory
        scriptcode_dir = os.path.join(self.dirn,"ScriptCode")
        bcf_utils.mkdir(scriptcode_dir,mode=0775)
        # Put a file in ScriptCode to make sure it's
        # not pruned on subsequent rsync operations
        fp = open(os.path.join(self.dirn,'ScriptCode','README.txt'),'w')
        fp.write("The ScriptCode directory is a place to put custom scripts and programs")
        fp.close()
        # Make a 'fastqs' directory
        if fastq_dir is None:
            fastq_dir = "fastqs"
        fastq_dir = os.path.join(self.dirn,fastq_dir)
        bcf_utils.mkdir(fastq_dir,mode=0775)
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
            reduced_fastqs = filter(lambda fq: not fq.is_index_read,
                                    reduced_fastqs)
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
        """Determine if directory really is an analysis project

        """
        return len(self.samples) > 0

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
        """
        return QCReporter(self)

    def qc_report(self,title=None,report_html=None,qc_dir=None,
                  force=False):
        """
        Report QC outputs for project

        Generates HTML and zipped QC reports.

        Arguments:
          title (str): title text for the report
          report_html (str): path for output HTML report file
          qc_dir (str): path for QC output dir (if None then
            use default QC directory)
          force (bool): if True then force reports to be
            regenerated (by default reports will not be
            regenerated if they already exist)

        Returns:
          String: name of zip file, or None if there was a
            problem.
        """
        if qc_dir is None:
            qc_dir = self._qc_dir
        elif not os.path.isabs(qc_dir):
            qc_dir = os.path.join(self.dirn,qc_dir)
        if not (force or self.verify_qc(qc_dir=qc_dir)):
            logger.debug("Failed to generate QC report for %s: QC "
                          "not verified and force not specified"
                          % self.name)
            return None
        # Create HTML report
        logger.debug("Creating HTML QC report for %s" % self.name)
        try:
            if not title:
                if self.info.run is not None:
                    title = "%s/%s: QC report" % (self.info.run,self.name)
                else:
                    title = "%s: QC report" % self.name
            if report_html is None:
                report_html = os.path.join(self.dirn,"qc_report.html")
            self.qc.report(title=title,
                           filename=report_html,
                           qc_dir=qc_dir,
                           relative_links=True)
        except Exception as ex:
            logger.error("Exception trying to generate QC report "
                         "for %s: %s" % (self.name,ex))
            return None
        # Get a name for the zip file derived from the HTML
        # report filename
        zip_name = os.path.splitext(os.path.basename(report_html))[0]
        # Create zip file
        logger.debug("Creating zip archive of QC report for %s" %
                      self.name)
        try:
            analysis_dir = os.path.basename(os.path.dirname(self.dirn))
            report_zip = os.path.join(self.dirn,
                                      "%s.%s.%s.zip" %
                                      (zip_name,
                                       self.name,
                                       analysis_dir))
            zip_file = ZipArchive(report_zip,relpath=self.dirn,
                                  prefix="%s.%s.%s" %
                                  (zip_name,
                                   self.name,
                                   analysis_dir))
            # Add the HTML report
            zip_file.add_file(report_html)
            # Add the FastQC and screen files
            for sample in self.qc.samples:
                for fastqs in sample.fastq_pairs:
                    for fq in fastqs:
                        logger.debug("Adding QC outputs for %s" % fq)
                        for f in expected_qc_outputs(fq,qc_dir):
                            if f.endswith('.zip'):
                                # Exclude .zip file
                                continue
                            if os.path.exists(f):
                                zip_file.add(f)
            # Finished
            return report_zip
        except Exception as ex:
            logger.error("Exception trying to generate zip archive "
                         "of QC report for %s: %s" % (self.name,ex))
            return None

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

        Arguments:
          qc_dir (str): path for QC output dir (if None then
            use default QC directory)

        Returns:
          Boolean: True if QC run is completed, False
            if QC couldn't be verified.
        """
        if qc_dir is None:
            qc_dir = self._qc_dir
        elif not os.path.isabs(qc_dir):
            qc_dir = os.path.join(self.dirn,qc_dir)
        try:
            return self.qc.verify(qc_dir=qc_dir)
        except AttributeError:
            return False

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
        raise KeyError, "No matching sample for '%s'" % name

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

class AnalysisSample:
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

        Returns:
          Populated QCSample object.

        """
        return QCSample(self)

    def verify_qc(self,qc_dir,fastq):
        """Check if QC completed for a fastq file

        Arguments:
          qc_dir: name of the QC directory
          fastq : fastq file to get the QC information for

        Returns:
          True if QC completed correctly, False otherwise.

        """
        present,missing = check_qc_outputs(fastq,qc_dir)
        if missing:
            return False
        return True

    def __repr__(self):
        """Implement __repr__ built-in

        Return string representation for the sample -
        i.e. the sample name.

        """
        return str(self.name)

class OutputFiles:
    """Class for managing multiple output files

    Usage:

    Create a new OutputFiles instance:
    >>> fp = OutputFiles()

    Set up files against keys:
    >>> fp.open('file1','first_file.txt')
    >>> fp.open('file2','second_file.txt')

    Write content to files:
    >>> fp.write('file1','some content for first file')
    >>> fp.write('file2','content for\nsecond file')

    Append content to an existing file:
    >>> fp.open('file3','third_file.txt',append=True)
    >>> fp.write('file2','appended content')

    Check if key exists and associated file handle is
    available for writing:
    >>> 'file1' in fp
    True
    >>> 'file3' in fp
    False

    Finish and close all open files
    >>> fp.close()

    Reopen and append to a previously opened and closed
    file:
    >>> fp.open('file4','fourth_file.txt')
    >>> fp.write('file4','some content')
    >>> fp.close('file4')
    >>> fp.open('file4',append=True)
    >>> fp.write('file4','more content')

    """
    def __init__(self,base_dir=None):
        """Create a new OutputFiles instance

        Arguments:
          base_dir (str): optional 'base' directory
            which files will be created relative to

        """
        self._fp = dict()
        self._file = dict()
        self._base_dir = base_dir

    def open(self,name,filen=None,append=False):
        """Open a new output file

        'name' is the handle used to reference the
        file when using the 'write' and 'close' methods.

        'filen' is the name of the file, and is unrelated
        to the handle. If not supplied then 'name' must
        be associated with a previously closed file (which
        will be reopened).

        If 'append' is True then append to an existing
        file rather than overwriting (i.e. use mode 'a'
        instead of 'w').

        If 'append' is True then append to an existing
        file rather than overwriting (i.e. use mode 'a'
        instead of 'w').

        """
        if append:
            mode = 'a'
        else:
            mode = 'w'
        if filen is None:
            filen = self.file_name(name)
        elif self._base_dir is not None:
            filen = os.path.join(self._base_dir,filen)
        else:
            filen = os.path.abspath(filen)
        self._file[name] = filen
        self._fp[name] = open(filen,mode)

    def write(self,name,s):
        """Write content to file (newline-terminated)

        Writes 's' as a newline-terminated string to the
        file that is referenced with the handle 'name'.

        """
        self._fp[name].write("%s\n" % s)

    def file_name(self,name):
        """Get the file name associated with a handle

        NB the file name will be available even if the
        file has been closed.

        Raises KeyError if the key doesn't exist.

        """
        return self._file[name]

    def close(self,name=None):
        """Close one or all open files

        If a 'name' is specified then only the file matching
        that handle will be closed; with no arguments all
        open files will be closed.

        """
        if name is not None:
            self._fp[name].close()
            del(self._fp[name])
        else:
            names = self._fp.keys()
            for name in names:
                self.close(name)

    def __contains__(self,name):
        return name in self._fp

    def __len__(self):
        return len(self._fp.keys())

class ZipArchive(object):
    """
    Utility class for creating .zip archive files

    Example usage:

    >>> z = ZipArchive('test.zip',relpath='/data')
    >>> z.add('/data/file1') # Add a single file
    >>> z.add('/data/dir2/') # Add a directory and all contents
    >>> z.close()  # to write the archive

    """
    def __init__(self,zip_file,contents=None,relpath=None,prefix=None):
        """
        Make an new zip archive instance

        Arguments:
          zip_file (str): path to the zip file to be created
          contents (list): list of file and/or directory paths
            which will be added to the zip file
          relpath (str): optional, if specified then this path
            will be stripped from the leading path for each item
            before being written (see also 'prefix')
          prefix (str): optional, if specified then this path
            will be prepended to the names of the items written
            to the archive. The prepending takes place after the
            relpath argument has been applied

        """
        self._zipfile = zipfile.ZipFile(zip_file,'w',
                                        allowZip64=True)
        self._relpath = relpath
        self._prefix = prefix
        if contents is not None:
            for item in contents:
                self.add(item)

    def add(self,item):
        """
        Add an item (file or directory) to the zip archive
        """
        item = os.path.abspath(item)
        if os.path.isfile(item):
            # Add file
            self.add_file(item)
        elif os.path.isdir(item):
            # Add directory and contents
            self.add_dir(item)
        else:
            raise Exception("ZipArchive: unknown item type for '%s'"
                            % item)

    def add_file(self,filen):
        """
        Add a file to the zip archive
        """
        if self._relpath:
            zip_pth = os.path.relpath(filen,self._relpath)
        else:
            zip_pth = filen
        if self._prefix:
            zip_pth = os.path.join(self._prefix,zip_pth)
        self._zipfile.write(filen,zip_pth)

    def add_dir(self,dirn):
        """
        Recursively add a directory and its contents
        """
        for item in os.listdir(dirn):
            f = os.path.join(dirn,item)
            if os.path.isdir(f):
                self.add_dir(f)
            else:
                self.add_file(f)

    def close(self):
        self._zipfile.close()

    def __del__(self):
        self.close()

class ProgressChecker(object):
    """
    Check if an index is a multiple of a value or percentage

    Utility class to help with reporting progress of iterations
    over large numbers of items.

    Typically progress would only be reported after a certain
    number or percentage of items have been consumed; the
    ProgressChecker can be used to check if this number or
    percentage has been reached.

    Example usage: to report after every 100th item:

    >>> progress = ProgressChecker(every=100)
    >>> for i in range(10000):
    >>>    if progress.check(i):
    >>>       print "Item %d" % i

    To report every 5% of items:

    >>> nitems = 10000
    >>> progress = ProgressChecker(percent=5,total=nitems)
    >>> for i in range(nitems):
    >>>    if progress.check(i):
    >>>       print "Item %d (%.2f%%)" % (i,progress.percent(i))
    """
    def __init__(self,every=None,percent=None,total=None):
        """
        Create a new ProgressChecker instance

        Arguments:
          every (int): specify interval number of items
            for reporting
          percent (float): specify a percentage interval
          total (int): total number of items (must be
            provided if using `percent`)
        """
        if every is None:
            every = max(int(float(total)*float(percent)/100.0),1)
        self._every = int(every)
        self._total = total

    def check(self,i):
        """
        Check index to see if it matches the interval

        Arguments:
          i (int): index to check

        Returns:
          Boolean: True if index matches the interval,
            False if not.
        """
        return (i%self._every == 0)

    def percent(self,i):
        """
        Convert index to a percentage

        Arguments:
          i (int): index to convert

        Returns:
          Float: index expressed as a percentage of the
            total number of items.
        """
        return float(i)/float(self._total)*100.0

#######################################################################
# Functions
#######################################################################

def bases_mask_is_paired_end(bases_mask):
    # Determine if run is paired end based on bases mask string
    non_index_reads = []
    for read in bases_mask.split(','):
        try:
            read.index('I')
        except ValueError:
            non_index_reads.append(read)
    if len(non_index_reads) == 2:
        # Paired end
        return True
    elif len(non_index_reads) < 2:
        # Single end
        return False
    else:
        # An error?
        raise Exception, "Bad bases mask '%s'?" % bases_mask

def split_user_host_dir(location):
    # Split a location of the form [[user@]host:]dir into its
    # user, hostname and directory components
    try:
        location = location.strip()
    except AttributeError:
        # Not a string?
        logger.error("Bad input to split_user_host_dir: '%s'" % location)
        return (None,None,None)
    if not location:
        return (None,None,None)
    try:
        location.index(':')
        location,dirn = location.split(':')
        try:
            location.index('@')
            user,host = location.split('@')
        except ValueError:
            user = None
            host = location
    except ValueError:
        user = None
        host = None
        dirn = location
    return (user,host,dirn)

def get_numbered_subdir(name,parent_dir=None,full_path=False):
    """
    Return a name for a new numbered log subdirectory

    Generates the name for a numbered subdirectory.

    Subdirectories are named as NNN_<name>  e.g.
    001_setup, 002_make_fastqs etc.

    'Gaps' are ignored, so the number associated with
    the new name will be one plus the highest index
    that already exists.

    **Note that a directory is not created** - this
    must be done by the calling subprogram. As a
    result there is the possibility of a race
    condition.

    Arguments:
      name (str): name for the subdirectory
        (typically the name of the processing
        stage that will produce logs to be
        written to the subdirs
      parent_dir (str): path to the parent
        directory where the indexed directory
        would be created; defaults to CWD if
        not set
      full_path (bool): if True then return the
        full path for the new subdirectory;
        default is to return the name relative
        to the parent directory

    Returns:
      String: name for the new log subdirectory
        (will be the full path if 'full_path'
        was specified).
    """
    # Sort out parent directory
    if parent_dir is None:
        parent_dir = os.getcwd()
    parent_dir = os.path.abspath(parent_dir)
    # Get the highest number from the names of
    # any other existing numbered subdirs
    i = 0
    for d in bcf_utils.list_dirs(parent_dir):
        try:
            i = max(i,int(d.split('_')[0]))
        except ValueError:
            pass
    # Generate and return name/path
    subdir = "%03d_%s" % (i+1,str(name))
    if full_path:
        subdir = os.path.join(parent_dir,subdir)
    return subdir

def find_executables(names,info_func,reqs=None,paths=None):
    """
    List available executables matching list of names

    By default searches the PATH for the executables listed
    in 'names', using the supplied 'info_func' to acquire
    package names and versions of each, returns a list of
    executables with the full path, package and version.

    'info_func' is a function that must be supplied by the
    calling subprogram. Its signature should look like:

    >>> def info_func(p):
    ...   # Determine full_path, package_name and
    ...   # version
    ...   # Then return these as a tuple
    ...   return (full_path,package_name,version)

    The 'reqs' argument allows a specific version or range
    of versions to be requested; in this case the returned
    list will only contain those packages which satisfy
    the requested versions.

    A range of version specifications can be requested by
    separating multiple specifiers with a comma - for
    example '>1.8.3,<2.16'.

    The full set of operators is:

    - ==, >, >=, <=, <

    If no versions are requested then the packages will
    be returned in PATH order; otherwise they will be
    returned in version order (highest to lowest).

    Arguments:
      names (list): list of executable names to look for.
        These can be full paths or executables with no
        leading paths
      info_func (function): function to use to get tuples
        of (full_path,package_name,version) for an
        executable
      reqs (str): optional version requirement expression
        (for example '>=1.8.4'). If supplied then only
        executables fulfilling the requirement will be
        returned. If no operator is supplied then '=='
        is implied.
      paths (list): optional set of directory paths to
        search when looking for executables. If not
        supplied then the set of paths specified in
        the PATH environment variable will be searched.

    Returns:
      List: full paths to executables matching specified
        criteria.

    """
    # Search paths
    if paths is None:
        paths = os.environ['PATH'].split(os.pathsep)
    # Search for executables
    available_exes = []
    for path in paths:
        if os.path.isfile(path):
            path = os.path.dirname(path)
        for name in names:
            prog_path = os.path.abspath(os.path.join(path,name))
            if bcf_utils.PathInfo(prog_path).is_executable:
                available_exes.append(prog_path)
    # Filter on requirement
    if reqs:
        # Loop over ranges
        for req in reqs.split(','):
            logging.debug("Filtering on expression: %s" % req)
            # Determine operator and version
            req_op = None
            req_version = None
            for op in ('==','>=','<=','>','<'):
                if req.startswith(op):
                    req_op = op
                    req_version = req[len(op):].strip()
                    break
            if req_version is None:
                req_op = '=='
                req_version = req.strip()
            logging.debug("Required version: %s %s" % (req_op,req_version))
            if req_op == '==':
                op = operator.eq
            elif req_op == '>=':
                op = operator.ge
            elif req_op == '>':
                op = operator.gt
            elif req_op == '<':
                op = operator.lt
            elif req_op == '<=':
                op = operator.le
            # Filter the available executables on version
            logging.debug("Pre filter: %s" % available_exes)
            logging.debug("Versions  : %s" % [info_func(x)[2]
                                              for x in available_exes])
            available_exes = filter(lambda x: op(
                parse_version(info_func(x)[2]),
                parse_version(req_version)),
                                    available_exes)
            logging.debug("Post filter: %s" % available_exes)
        # Sort into version order, highest to lowest
        available_exes.sort(
            key=lambda x: parse_version(info_func(x)[2]),
            reverse=True)
        logging.debug("Post sort: %s" % available_exes)
        logging.debug("Versions : %s" % [info_func(x)[2]
                                 for x in available_exes])
    if available_exes:
        print "%s: found available packages:" % ','.join(names)
        for i,package in enumerate(available_exes):
            package_info = info_func(package)
            print "%s %s\t%s\t%s" % (('*' if i == 0 else ' '),
                                     package_info[0],
                                     package_info[1],
                                     package_info[2])
    else:
        print "%s: No packages found" % ','.join(names)
    return available_exes

def parse_version(s):
    """
    Split a version string into a tuple for comparison

    Given a version string of the form e.g. "X.Y.Z",
    return a tuple of the components e.g. (X,Y,Z)

    Where possible components will be coverted to
    integers.

    If the version string is empty then the version
    number will be set to an arbitrary negative
    integer.

    Typically the result from this function would not
    be used directly, instead it is used to compare
    two versions, for example:

    >>> parse_version("2.17") < parse_version("1.8")
    False

    Arguments:
      s (str): version string

    Returns:
      Tuple: tuple of the version string
    """
    if s == "":
        # Essentially versionless; set to an
        # arbitrarily small integer
        s = "-99999"
    items = []
    for i in s.split('.'):
        try:
            i = int(i)
        except ValueError:
            pass
        items.append(i)
    return tuple(items)

def pretty_print_rows(data,prepend=False):
    """Format row-wise data into 'pretty' lines

    Given 'row-wise' data (in the form of a list of lists),
    for example:

    [['hello','A salutation'],[goodbye','The End']]

    formats into a string of text where lines are
    newline-separated and the 'fields' are padded with
    spaces so that they line up left-justified in
    columns, for example:

    hello   A salutation
    goodbye The End

    Arguments:
      data: row-wise data as a list of lists
      prepend: (optional), if True then columns
        are right-justified (i.e. padding is
        added before each value).

    """
    # Get maximum field widths for each column
    widths = []
    for row in data:
        for i in range(len(row)):
            width = len(str(row[i]))
            try:
                widths[i] = max(width,widths[i])
            except IndexError:
                widths.append(width)
    # Build output
    output = []
    for row in data:
        line = []
        for item,width in zip([str(x) for x in row],widths):
            padding = ' '*(width-len(item))
            if prepend:
                line.append(padding + item)
            else:
                line.append(item + padding)
        output.append(' '.join(line))
    return '\n'.join(output)

def write_script_file(script_file,contents,append=False,shell=None):
    """Write command to file

    Arguments:
      script_file (str): path of file to write command to
      contents (str): content to write to the file
      append (bool): optional, if True and script_file exists
        then append content (default is to overwrite existing
        contents) 
      shell: optional, if set then defines the shell to
        specify after '!#'

    """
    if append:
        mode = 'a'
    else:
        mode = 'w'
    with open(script_file,mode=mode) as fp:
        if (not append) and (shell is not None):
            fp.write("#!%s\n" % shell)
        fp.write("%s\n" % contents)
    os.chmod(script_file,0775)

def edit_file(filen,editor="vi",append=None):
    """
    Send a file to an editor

    Creates a temporary copy of a file and opens an
    editor to allow the user to make changes. Any
    edits are saved back to the original file.

    Arguments:
      filen (str): path to the file to be edited
      editor (str): optional, editor command to be used
        (will be overriden by user's EDITOR environment
        variable even if set). Defaults to 'vi'.
      append (str): optional, if set then append the
        supplied text to the end of the file before
        editing. NB the text will only be kept if the
        user saves a change to the file in the editor.

    """
    # Acquire an editor command
    try:
        editor = os.environ["EDITOR"]
    except KeyError:
        pass
    if editor is None:
        logger.critical("No editor specified!")
        return
    # Make a temporary copy for editing
    f,tmpfile = tempfile.mkstemp()
    os.fdopen(f).close()
    with open(tmpfile,'w') as fp:
        if os.path.exists(filen):
            fp.write(open(filen,'r').read())
        else:
            fp.write()
        if append:
            fp.write("%s\n" % str(append))
    checksum = md5sum(tmpfile)
    # Build command line to run the editor
    editor = str(editor).split(' ')
    edit_cmd = applications.Command(editor[0],*editor[1:])
    edit_cmd.add_args(tmpfile)
    edit_cmd.run_subprocess()
    # Finished
    if md5sum(tmpfile) != checksum:
        with open(filen,'w') as fp:
            fp.write(open(tmpfile,'r').read())
            os.remove(tmpfile)
    else:
        logger.warning("no changes to write")

def paginate(text):
    """
    Send text to stdout with pagination

    If the function detects that the stdout is an interactive
    terminal then the supplied text will be piped via a
    paginator command.

    The pager command will be the default for ``pydoc``, but
    can be over-ridden by the ``PAGER`` environment variable.

    If stdout is not a terminal (for example if it's being
    set to a file, or piped to another command) then the
    pagination is skipped.

    Arguments:
      text (str): text to be printed using pagination

    """
    # If stdout is a terminal
    if os.isatty(sys.stdout.fileno()):
        # Acquire a pager command
        try:
            pager = os.environ["PAGER"]
        except KeyError:
            pager = None
        # Output the prediction with paging
        if pager is not None:
            pydoc.pipepager(text,cmd=pager)
        else:
            pydoc.pager(text)
    else:
        # Stdout not a terminal
        print text
