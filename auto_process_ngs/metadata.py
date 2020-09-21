#!/usr/bin/env python
#
#     metadata: classes for storing metadata on analysis objects
#     Copyright (C) University of Manchester 2018-2020 Peter Briggs
#
########################################################################
#
# metadata.py
#
#########################################################################

"""
metadata

Classes for storing, accessing and updating metadata for analysis
directories, projects and so on.

Classes:

- MetadataDict:
- AnalysisDirParameters:
- AnalysisDirMetadata:
- AnalysisProjectInfo:
- ProjectMetadataFile:
- AnalysisProjectQCDirInfo:
"""

#######################################################################
# Imports
#######################################################################

import os
import uuid
import logging
import bcftbx.TabFile as TabFile
import bcftbx.utils as bcf_utils

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Classes
#######################################################################

class MetadataDict(bcf_utils.AttributeDictionary):
    """Class for storing metadata in an analysis project

    Provides storage for arbitrary data items in the form of
    key-value pairs, which can be saved to and loaded from
    an external file.

    The data items are defined on instantiation via a dictionary
    supplied to the 'attributes' argument. For example:

    Create a new metadata object:
    >>> metadata = MetadataDict(attributes={'salutation':'Salutation',
    ...                                     'valediction': 'Valediction'})
 
    The dictionary keys correspond to the keys in the MetadataDict
    object; the corresponding values are the keys that are used
    when saving and loading the data to and from a file.

    Set attributes:
    >>> metadata['salutation'] = 'hello'
    >>> metadata['valediction'] = 'goodbye'

    Retrieve values:
    >>> print("Salutation is %s" % metadata.salutation)

    Save to file:
    >>> metadata.save('metadata.tsv')

    Load data from a file:
    >>> metadata = MetadataDict('metadata.tsv')
    or
    >>> metadata = MetadataDict()
    >>> metadata.load('metadata.tsv')

    List items with 'null' values:
    >>> metadata.null_items()

    The external file storage is intended to be readable by
    humans so longer names are used to describe the keys; also
    Python None values are stored as '.', and True and False
    values are stored as 'Y' and 'N' respectively. These values
    are automatically converted back to the Python equivalents
    on reload.

    """

    def __init__(self,attributes=dict(),order=None,filen=None):
        """Create a new MetadataDict object

        By default an empty metadata object is created
        i.e. all attributes will have be None.

        If an input file is specified then the attributes
        will be assigned values according to the key-value
        pairs in that file.

        Arguments:
          attributes: dictionary defining metadata items
          filen: (optional) name of the tab-delimited file
            with key-value pairs to load in.

        """
        bcf_utils.AttributeDictionary.__init__(self)
        self.__filen = filen
        # Set up empty metadata attributes
        self.__attributes = attributes
        for key in self.__attributes:
            self[key] = None
        if self.__filen:
            # Load data from external file
            if os.path.exists(self.__filen):
                self.load(self.__filen)
        # Set up order of keys for output
        if order is None:
            self.__key_order = sorted(list(self.__attributes.keys()))
        else:
            # Use supplied key order
            self.__key_order = []
            for key in order:
                if key in self.__attributes:
                    self.__key_order.append(key)
                else:
                    raise KeyError("Key '%s' not defined in attributes")
            # Append keys not explicitly listed in the order
            extra_keys = []
            for key in self.__attributes:
                if key not in self.__key_order:
                    extra_keys.append(key)
            if extra_keys:
                extra_keys.sort()
                self.__key_order.extend(extra_keys)

    def __iter__(self):
        return iter(self.__key_order)

    def load(self,filen,strict=True,fail_on_error=False):
        """Load key-value pairs from a tab-delimited file
        
        Loads the key-value pairs from a previously created
        tab-delimited file written by the 'save' method.

        Note that this overwrites any existing values
        already assigned to keys within the metadata object.

        Arguments:
          filen (str): name of the tab-delimited file with
            key-value pairs
          strict (bool): if True (default) then discard
            items in the input file which are missing from
            the definition; if False then add them to the
            definition.
          fail_on_error (bool): if True then raise an
            exception if the file contains invalid content
            (if 'strict' is also specified then this
            includes any unrecognised keys); default is
            to warn and then ignore these errors.

        """
        self.__filen = filen
        metadata = TabFile.TabFile(filen)
        for line in metadata:
            try:
                # Get data from file and convert special values
                # to Python equivalents
                attr,value = line[0],line[1]
                if value == '.' or value == 'None':
                    value = None
                elif value == 'Y' or value == 'True':
                    value = True
                elif value == 'N' or value == 'False':
                    value = False
                # Locate dictionary key matching file key
                found_key = False
                for key in self.__attributes:
                    if self.__attributes[key] == attr:
                        self[key] = value
                        found_key = True
                        break
                if not found_key:
                    if strict:
                        logger.warning("Unrecognised key in %s: %s"
                                       % (filen,attr))
                        if fail_on_error:
                            raise Exception("%s: failed to load: bad key"
                                            % filen)
                    else:
                        logger.debug("Adding key from %s: %s"
                                     % (filen,attr))
                        self.__attributes[attr] = attr
                        self.__key_order.append(attr)
                        self[attr] = value
            except IndexError:
                logger.warning("Bad line in %s: %s" % (filen,line))
                if fail_on_error:
                    raise Exception("%s: failed to load: bad line"
                                    % filen)

    def save(self,filen=None):
        """Save metadata to tab-delimited file

        Writes key-value paires to a tab-delimited file.
        The data can be recovered using the 'load' method.
 
        Note that if the specified file already exists then
        it will be overwritten.

        Arguments:
          filen: name of the tab-delimited file with key-value
            pairs; if None then the file specified when the
            object was instantiated will be used instead.

        """
        metadata = TabFile.TabFile()
        for key in self.__key_order:
            # Retrieve value and convert to appropriate
            # format for persistent storage
            value = self[key]
            if value is None:
                value = '.'
            elif value is True:
                value = 'Y'
            elif value is False:
                value = 'N'
            # Get the equivalent file key
            attr = self.__attributes[key]
            # Store in the file
            metadata.append(data=(attr,value))
        # Write data to temporary file
        if filen is not None:
            self.__filen = filen
        if self.__filen is None:
            # Nowehere to write to
            logger.warning("Cannot save metadata: no destination file "
                           "defined")
            return
        tmp_filen = os.path.join(
            os.path.dirname(self.__filen),
            "%s.%s.tmp" % (os.path.basename(self.__filen),
                           uuid.uuid4()))
        metadata.write(tmp_filen)
        # Move to final destination
        os.rename(tmp_filen,self.__filen)

    def null_items(self):
        """
        Return a list of data items with 'null' values

        """
        null_items = []
        for key in self.__key_order:
            if self[key] is None:
                null_items.append(key)
        return null_items

class AnalysisDirParameters(MetadataDict):
    """Class for storing parameters in an analysis directory

    Provides a set of data items representing parameters for
    the current analysis, which are loaded from and saved to
    an external file.

    The parameter data items are:

    analysis_dir: path to the analysis directory
    data_dir: path to the directory holding the raw sequencing data
    platform: sequencing platform e.g. 'miseq'
    sample_sheet: path to the customised SampleSheet.csv file
    bases_mask: bases mask string
    project_metadata: name of the project metadata file
    primary_data_dir: directory used to hold copies of primary data
    acquired_primary_data: whether primary data has been copied
    unaligned_dir: output directory for bcl2fastq conversion
    barcode_analysis_dir: directory holding barcode analysis outputs
    stats_file: name of file with per-fastq statistics
    per_lane_stats_file: name of file with per-lane statistics

    """
    def __init__(self,filen=None):
        """Create a new AnalysisDirParameters object

        Arguments:
          filen (str): (optional) name of the tab-delimited
            file with key-value pairs to load in.

        """
        MetadataDict.__init__(self,
                              attributes = {
                                  'analysis_dir':'analysis_dir',
                                  'data_dir':'data_dir',
                                  'sample_sheet':'sample_sheet',
                                  'bases_mask':'bases_mask',
                                  'project_metadata':'project_metadata',
                                  'primary_data_dir':'primary_data_dir',
                                  'acquired_primary_data':'acquired_primary_data',
                                  'unaligned_dir':'unaligned_dir',
                                  'barcode_analysis_dir':'barcode_analysis_dir',
                                  'stats_file':'stats_file',
                                  'per_lane_stats_file':'per_lane_stats_file',
                              },
                              filen=filen)

class AnalysisDirMetadata(MetadataDict):
    """Class for storing metadata about an analysis directory

    Provides a set of data items representing metadata about
    the current analysis, which are loaded from and saved to
    an external file.

    The metadata items are:

    run_name: name of the run
    run_number: run number assigned by local facility
    source: source of the data (e.g. local facility)
    platform: sequencing platform e.g. 'miseq'
    assay: the 'assay' from the IEM SampleSheet e.g. 'Nextera XT'
    bcl2fastq_software: info on the Bcl conversion software used
    cellranger_software: info on the 10xGenomics cellranger software
      used
    instrument_name: name/i.d. for the sequencing instrument
    instrument_datestamp: datestamp from the sequencing instrument
    instrument_run_number: the run number from the sequencing
      instrument
    instrument_flow_cell_id: the flow cell ID from the sequencing
      instrument
    sequencer_model: the model of the sequencing instrument

    """
    def __init__(self,filen=None):
        """Create a new AnalysisDirMetadata object

        Arguments:
          filen (str): (optional) name of the tab-delimited
            file with key-value pairs to load in.

        """
        MetadataDict.__init__(self,
                              attributes = {
                                  'run_name':'run_name',
                                  'run_number': 'run_number',
                                  'source': 'source',
                                  'platform':'platform',
                                  'assay': 'assay',
                                  'bcl2fastq_software': 'bcl2fastq_software',
                                  'cellranger_software': 'cellranger_software',
                                  'instrument_name': 'instrument_name',
                                  'instrument_datestamp': 'instrument_datestamp',
                                  'instrument_flow_cell_id': 'instrument_flow_cell_id',
                                  'instrument_run_number': 'instrument_run_number',
                                  'sequencer_model': 'sequencer_model',
                              },
                              order=('run_name',
                                     'platform',
                                     'run_number',
                                     'source',),
                              filen=filen)

class AnalysisProjectInfo(MetadataDict):
    """Class for storing metadata in an analysis project

    Provides a set of metadata items which are loaded from
    and saved to an external file.

    The data items are:

    name: the project name
    run: the name of the sequencing run
    platform: the sequencing platform name e.g. 'miseq'
    sequencer_model: the sequencer model e.g. 'MiSeq'
    user: the user associated with the project
    PI: the principal investigator associated with the project
    organism: the organism associated with the project
    library_type: the library type e.g. 'RNA-seq'
    single_cell_platform: the single cell preparation platform
    number of cells: number of cells in single cell projects
    ICELL8 well list: well list file for ICELL8 single cell projects
    paired_end: True if the data is paired end, False if not
    primary_fastq_dir: the primary subdir with FASTQ files
    samples: textual description of the samples in the project
    comments: free-text comments

    """
    def __init__(self,filen=None):
        """Create a new AnalysisProjectInfo object

        Arguments:
          filen: (optional) name of the tab-delimited file
            with key-value pairs to load in.

        """
        MetadataDict.__init__(self,
                              attributes = {
                                  'name':'Project name',
                                  'run':'Run',
                                  'platform':'Platform',
                                  'sequencer_model':'Sequencer model',
                                  'user':'User',
                                  'PI':'PI',
                                  'organism':'Organism',
                                  'library_type':'Library type',
                                  'single_cell_platform':'Single cell platform',
                                  'number_of_cells':'Number of cells',
                                  'icell8_well_list':'ICELL8 well list',
                                  'paired_end':'Paired_end',
                                  'primary_fastq_dir':'Primary fastqs',
                                  'samples':'Samples',
                                  'comments':'Comments',
                              },
                              order = (
                                  'name',
                                  'run',
                                  'platform',
                                  'user',
                                  'PI',
                                  'organism',
                                  'library_type',
                                  'single_cell_platform',
                                  'number_of_cells',
                                  'icell8_well_list',
                                  'paired_end',
                                  'primary_fastq_dir',
                                  'samples',
                                  'sequencer_model',
                                  'comments',
                              ),
                              filen=filen)

class ProjectMetadataFile(TabFile.TabFile):
    """File containing metadata about multiple projects in analysis dir
 
    The file consists of a header line plus one line per project
    with the following tab-delimited fields:

    Project: name of the project
    Samples: list/description of sample names
    User: name(s) of the associated user(s)
    Library: the library type
    Single cell platform: single-cell preparation platform (e.g. 'ICELL8')
    Organism: name(s) of the organism(s)
    PI: name(s) of the associated principal investigator(s)
    Comments: free text containing additional information
              about the project

    Any fields set to None will be written to file with a '.'
    placeholder.
    """
    def __init__(self,filen=None):
        """Create a new ProjectsMetadataFile instance

        Arguments:
          filen: (optional) name of an existing file to read
            projects in from.

        """
        # List of expected fields
        # Add new fields to this list
        self._default_fields = ('Project',
                                'Samples',
                                'User',
                                'Library',
                                'SC_Platform',
                                'Organism',
                                'PI',
                                'Comments')
        # Map keywords to column names
        self._kwmap = { 'Project': 'project_name',
                        'Samples': 'sample_names',
                        'User': 'user',
                        'Library': 'library_type',
                        'SC_Platform': 'sc_platform',
                        'Organism': 'organism',
                        'PI' : 'PI',
                        'Comments': 'comments', }
        # List of default values
        self._default_values = { }
        # Optional file to read from
        self.__filen = filen
        if self.__filen is None:
            # No existing file so set the default
            # fields to write to the file
            self._fields = self._default_fields
        else:
            # Get columns from existing file
            with open(self.__filen,'r') as fp:
                header = fp.readline()
                self._fields = header.rstrip('\n').lstrip('#').split('\t')
        # Open the file
        TabFile.TabFile.__init__(self,filen=self.__filen,
                                 column_names=self._fields,
                                 skip_first_line=True,
                                 convert=False,
                                 keep_commented_lines=True)
        # Add any missing columns
        for field in self._default_fields:
            if field not in self._fields:
                self.appendColumn(field)

    def add_project(self,project_name,sample_names,**kws):
        """Add information about a project into the file

        Arguments:
          project_name (str): name of the new project
          sample_names (list): Python list of sample names
          user (str): (optional) user name(s)
          library_type (str): (optional) library type
          sc_platform (str): (optional) single-cell prep platform
          organism (str): (optional) organism(s)
          PI (str): (optional) principal investigator name(s)
          comments (str): (optional) additional information
            about the project
        """
        # Check project name doesn't already exist
        if project_name in self:
            raise Exception("Project '%s' already exists" %
                            project_name)
        kws['project_name'] = project_name
        kws['sample_names'] = ','.join(sample_names)
        # Create an empty data line for the project
        project = self.append()
        # Assign the data
        for field in self._fields:
            # Identify the keyword parameter for this field
            try:
                kw = self._kwmap[field]
            except KeyError as ex:
                raise ex
            # Look up the assigned value
            try:
                value = kws[kw]
            except KeyError:
                value = None
            # Set value in the data line
            if value is None:
                project[field] = '.'
            else:
                project[field] = value

    def update_project(self,project_name,**kws):
        """Update information about a project in the file

        Arguments:
          project_name (str): name of the project to update
          sample_names (list): (optional) Python list of
            new sample names
          user (str): (optional) new user name(s)
          library_type (str): (optional) new library type
          sc_platform (str): (optional) single-cell prep platform
          organism (str): (optional) new organism(s)
          PI (str): (optional) new principal investigator
            name(s)
          comments (str): (optional) new additional
            information about the project
        """
        # Fetch data line for existing project
        project = self.lookup(project_name)
        # Set project name
        kws['project_name'] = project_name
        # Set sample names, if supplied
        try:
            kws['sample_names'] = ','.join(kws['sample_names'])
        except KeyError:
            pass
        # Update the data
        for field in self._fields:
            # Identify the keyword parameter for this field
            try:
                kw = self._kwmap[field]
            except KeyError as ex:
                raise ex
            # Assign the new values
            if kw not in kws:
                continue
            value = kws[kw]
            # Set value in the data line
            if value is None:
                project[field] = '.'
            else:
                project[field] = value

    def lookup(self,project_name):
        """
        Return data for line with specified project name

        Leading comment characters (i.e. '#') are ignored
        when performing the lookup.
        """
        pname = project_name.lstrip('#')
        for name in (pname,"#%s" % pname):
            try:
                return TabFile.TabFile.lookup(self,
                                              "Project",
                                              name)[0]
            except IndexError:
                pass
        raise KeyError("'%s': project not found in metadata file"
                       % project_name)

    def save(self,filen=None):
        """Save the data back to file

        Arguments:
          filen: name of the file to save to (if not specified then
            defaults to the same file as data was read in from)

        """
        # Sort into project name order
        self.sort(lambda line: str(line['Project']).lstrip('#'))
        # Write to file
        if filen is not None:
            self.__filen = filen
        self.write(filen=self.__filen,include_header=True)

    def __contains__(self,name):
        """
        Internal: check if field 'name' exists

        'name' should be the name of a project; if either the
        supplied name and/or the project name in the file are
        commented (i.e. preceeded by '#' symbol), then the
        leading '#' is ignored.
        """
        return (name.lstrip('#') in
                [p[self._fields[0]].lstrip('#') for p in self])

class AnalysisProjectQCDirInfo(MetadataDict):
    """Class for storing metadata for a QC output directory

    Provides a set of metadata items which are loaded from
    and saved to an external file.

    The data items are:

    fastq_dir: the name of the associated Fastq subdirectory
    protocol: the QC protocol used
    organism: the organism(s) that the QC was run with
    cellranger_refdata: reference datasets used with cellranger
    """
    def __init__(self,filen=None):
        """Create a new AnalysisProjectQCDirInfo instance

        Arguments:
          filen: (optional) name of the tab-delimited file
            with key-value pairs to load in.
        """
        MetadataDict.__init__(self,
                              attributes = {
                                  'fastq_dir':'Fastq dir',
                                  'protocol':'QC protocol',
                                  'organism':'Organism',
                                  'cellranger_refdata':
                                  'Cellranger reference datasets',
                              },
                              order = (
                                  'fastq_dir',
                                  'protocol',
                                  'organism',
                              ),
                              filen=filen)
