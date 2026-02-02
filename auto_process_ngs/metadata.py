#!/usr/bin/env python
#
#     metadata: classes for storing metadata on analysis objects
#     Copyright (C) University of Manchester 2018-2026 Peter Briggs
#
########################################################################
#
# metadata.py
#
#########################################################################

"""
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

    In 'strict' mode when data items are loaded from a file,
    any items which are not in the list of attributes will be
    discarded: they will not be accessible as attributes or
    keys, and will not be preserved on save.

    If 'strict' mode is turned off then any extra data items
    loaded from the file will be preserved on save, however
    they will still be inaccessible as attributes or keys.

    If 'include_undefined' is turned on then any extra data items
    will become accessible as attributes and keys, as well as
    being preserved on save. (Note that this functionality only
    works if 'strict' mode is turned off, otherwise an exception
    is raised when file load is attempted.)

    The attribute/key names for these extra data items within the
    MetadataDict object are created by replacing spaces in the
    item names to underscores and converting to lowercase (for
    example 'Bioinformatics Analyst' -> 'bioinformatics_analyst').
    """

    def __init__(self, attributes=dict(), order=None, filen=None,
                 strict=True, fail_on_error=False, enable_fallback=False,
                 include_undefined=False):
        """Create a new MetadataDict object

        By default an empty metadata object is created
        i.e. all attributes will be None.

        If an input file is specified then the attributes
        will be assigned values according to the key-value
        pairs in that file.

        Arguments:
          attributes: dictionary defining metadata items
          filen: (optional) name of the tab-delimited file
            with key-value pairs to load in
          strict (bool): if True then by default discard
            items from the input file which aren't defined
            in the 'attributes' dictionary (default: True)
          fail_on_error (bool): if True then raise an
            exception if the input file contains invalid
            content (if 'strict' is also specified then this
            includes any unrecognised keys) (default:
            False, errors will be ignored)
          enable_fallback (bool): if True then try matching
            keys directly if lookup fails when reading file
            (default: False, don't enable fallback)
        """
        bcf_utils.AttributeDictionary.__init__(self)
        self.__filen = filen
        self.__strict = bool(strict)
        self.__fail_on_error = bool(fail_on_error)
        self.__enable_fallback = bool(enable_fallback)
        self.__include_undefined = bool(include_undefined)
        # Set up empty metadata attributes
        self.__attributes = attributes.copy()
        for key in self.__attributes:
            self[key] = None
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
                    raise KeyError("Key '%s' not defined in attributes" % key)
            # Append keys not explicitly listed in the order
            extra_keys = []
            for key in self.__attributes:
                if key not in self.__key_order:
                    extra_keys.append(key)
            if extra_keys:
                extra_keys.sort()
                self.__key_order.extend(extra_keys)
        # Store any undefined keys
        self.__undefined_items = []
        # Load data from external file
        self.__file_keys = list()
        if self.__filen:
            if os.path.exists(self.__filen):
                self.load(self.__filen)

    def __iter__(self):
        return iter(self.__key_order)

    def load(self, filen, strict=None, fail_on_error=None,
             enable_fallback=None, include_undefined=None):
        """Load key-value pairs from a tab-delimited file
        
        Loads the key-value pairs from a previously created
        tab-delimited file written by the 'save' method.

        Note that this overwrites any existing values
        already assigned to keys within the metadata object.

        Arguments:
          filen (str): name of the tab-delimited file with
            key-value pairs
          strict (bool): if True then discard items in the
            input file which are missing from the definition;
            if False then add them to the definition. Defaults
            to the value supplied on creation (or True if not
            supplied)
          fail_on_error (bool): if True then raise an
            exception if the file contains invalid content
            (if 'strict' is also specified then this
            includes any unrecognised keys). Defaults to the
            value supplied on creation (or False if not
            supplied)
          enable_fallback (bool): if True then try matching
            keys directly if lookup fails when reading file.
            Defaults to the value supplied on creation (or
            False if not supplied)
        """
        self.__filen = filen
        if strict is None:
            strict = self.__strict
        if fail_on_error is None:
            fail_on_error = self.__fail_on_error
        if enable_fallback is None:
            enable_fallback = self.__enable_fallback
        if include_undefined is None:
            include_undefined = self.__include_undefined
        if include_undefined and strict:
            # Can't have both together
            raise Exception("MetadataDict.load: 'include_undefined' is "
                            "incompatible with 'strict'")
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
                found_key = None
                for key in self.__attributes:
                    if self.__attributes[key] == attr:
                        self[key] = value
                        found_key = key
                        break
                # Key wasn't found
                if found_key is None and enable_fallback:
                    # Fallback to matching keys directly
                    if attr in self.__attributes:
                        self[attr] = value
                        found_key = attr
                # Key still not found (even after fallback
                # was possibly attempted)
                if found_key is None:
                    if strict:
                        if fail_on_error:
                            # Raise an exception
                            raise Exception("%s: failed to load: bad key "
                                            "'%s'"% (filen, attr))
                        else:
                            # Warn and ignore; the item will not be
                            # added and will be discarded on save
                            logger.warning("Unrecognised key in %s: '%s'"
                                           % (filen, attr))

                    else:
                        if include_undefined:
                            # Add the item; it will be preserved on save
                            logger.debug("Adding key from %s: %s"
                                        % (filen,attr))
                            # Construct attribute name for save
                            name = name_to_item(attr)
                            # Add the undefined item
                            self.__attributes[name] = attr
                            self[name] = value
                            self.__key_order.append(name)
                            found_key = name
                        else:
                            # Store undefined item (so it can be preserved
                            # on output) but don't make it available
                            self.__undefined_items.append((attr, value))
                # Store keys found in file
                if found_key:
                    self.__file_keys.append(found_key)
            except IndexError:
                # Unable to parse the line
                if fail_on_error:
                    # Fatal error
                    raise Exception("%s: failed to load: bad line"
                                    % filen)
                else:
                    # Warn and continue
                    logger.warning("Bad line in %s (ignored): %s" %
                                   (filen,line))

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
        # Add undefined data from input
        for item, value in self.__undefined_items:
            metadata.append(data=(item, value))
        # Write data to temporary file
        if filen is not None:
            self.__filen = filen
        if self.__filen is None:
            # Nowhere to write to
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

    def keys_in_file(self):
        """
        Return a list of the key names found explicitly in the file

        """
        return [k for k in self.__file_keys]

    def null_items(self):
        """
        Return a list of data items with 'null' values

        """
        null_items = []
        for key in self.__key_order:
            if self[key] is None:
                null_items.append(key)
        return null_items

    def __len__(self):
        return len([k for k in self.__key_order if self[k] is not None])

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
    analysis_number: arbitrary number assigned to analysis (to
      distinguish it from other analysis attempts)
    source: source of the data (e.g. local facility)
    platform: sequencing platform e.g. 'miseq'
    processing_software: dictionary of software packages used in
      in the processing
    bcl2fastq_software: info on the Bcl conversion software used
      (deprecated)
    cellranger_software: info on the 10xGenomics cellranger software
      used (deprecated)
    instrument_name: name/i.d. for the sequencing instrument
    instrument_datestamp: datestamp from the sequencing instrument
    instrument_run_number: the run number from the sequencing
      instrument
    instrument_flow_cell_id: the flow cell ID from the sequencing
      instrument
    sequencer_model: the model of the sequencing instrument
    flow_cell_mode: the flow cell configuration (if present in the
      run parameters)
    run_configuration: read names & lengths derived from RunInfo.xml
    default_bases_mask: default bases mask derived from RunInfo.xml

    """
    def __init__(self,filen=None):
        """Create a new AnalysisDirMetadata object

        Arguments:
          filen (str): (optional) name of the tab-delimited
            file with key-value pairs to load in.

        """
        MetadataDict.__init__(self,
                              attributes = {
                                  'run_name': 'run_name',
                                  'run_number': 'run_number',
                                  'run_id': 'run_id',
                                  'run_reference_id': 'run_reference_id',
                                  'source': 'source',
                                  'platform':'platform',
                                  'processing_software': 'processing_software',
                                  'bcl2fastq_software': 'bcl2fastq_software',
                                  'cellranger_software': 'cellranger_software',
                                  'instrument_name': 'instrument_name',
                                  'instrument_datestamp': 'instrument_datestamp',
                                  'instrument_flow_cell_id': 'instrument_flow_cell_id',
                                  'instrument_run_number': 'instrument_run_number',
                                  'sequencer_model': 'sequencer_model',
                                  'flow_cell_mode': 'flow_cell_mode',
                                  'run_configuration': 'run_configuration',
                                  'default_bases_mask': 'default_bases_mask',
                                  'analysis_number': 'analysis_number',
                              },
                              order=('run_name',
                                     'platform',
                                     'run_number',
                                     'analysis_number',
                                     'source',),
                              filen=filen)

class AnalysisProjectInfo(MetadataDict):
    """Class for storing metadata in an analysis project

    Provides a set of metadata items which are loaded from
    and saved to an external file.

    The core data items are:

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
    paired_end: True if the data is paired end, False if not
    primary_fastq_dir: the primary subdir with FASTQ files
    samples: textual description of the samples in the project
    biological_samples: comma-separated sample names with biological data
    multiplexed_samples: comma-separated names of multiplexed samples
    comments: free-text comments

    Additional user-defined metadata items can be specified
    via the 'custom_items' argument; if supplied then this
    should be a list of metadata item names (e.g.
    'order_numbers').

    Custom metadata item names:

     - can only contain alphanumeric characters (i.e. letters
       and numbers)
     - cannot start with a number

    These names are used to access and update the metadata
    through the object's attributes and keys.

    When the custom items are saved to file then the names are
    converted as follows:

     - underscores become spaces
     - for all-lowercase items, the first letter is capitalized
       (e.g. "order_numbers" becomes "Order numbers")
     - for items containing uppercase letters, the case of the
       elements are preserved (e.g. "EOL_date" becomes "EOL date")

    Arguments:
      filen (str): name of a tab-delimited file with key-value
        pairs to load in
      custom_items (list): list of additional custom metadata
        items to add
    """
    def __init__(self, filen=None, custom_items=None):
        # Core metadata
        data_items = {
            'name': 'Project name',
            'run': 'Run',
            'platform': 'Platform',
            'sequencer_model': 'Sequencer model',
            'user': 'User',
            'PI': 'PI',
            'organism': 'Organism',
            'library_type': 'Library type',
            'single_cell_platform': 'Single cell platform',
            'number_of_cells': 'Number of cells',
            'paired_end': 'Paired_end',
            'primary_fastq_dir': 'Primary fastqs',
            'samples': 'Samples',
            'biological_samples': 'Biological samples',
            'multiplexed_samples': 'Multiplexed samples',
            'comments': 'Comments'
        }
        order = ['name',
                 'run',
                 'platform',
                 'user',
                 'PI',
                 'organism',
                 'library_type',
                 'single_cell_platform',
                 'number_of_cells',
                 'paired_end',
                 'primary_fastq_dir',
                 'samples',
                 'biological_samples',
                 'multiplexed_samples',
                 'sequencer_model',
                 'comments']
        # Additional custom items
        if custom_items:
            for item in custom_items:
                # Check custom item name
                if item[0].isdigit():
                    raise Exception(f"'{item}': metadata items must not start with a number")
                if any([not (c.isalnum() or c == "_") for c in item]):
                    raise Exception(f"'{item}': metadata items must only contain letters and underscores")
                if item[0].isupper() and not any([c.isupper() if c.isalpha() else False for c in item[1:]]):
                    raise Exception(f"'{item}': metadata items cannot be capitalized (use '{item.lower()}' instead)")
                # Create a name for writing to file
                name = item_to_name(item)
                data_items[item] = name
                order.append(item)
        MetadataDict.__init__(self,
                              attributes=data_items,
                              order=order,
                              filen=filen,
                              strict=False,
                              include_undefined=True)

class ProjectMetadataFile(TabFile.TabFile):
    """File containing metadata about multiple projects in analysis dir
 
    The file consists of a header line plus one line per project
    with the following tab-delimited fields:

    Project: name of the project
    Samples: list/description of sample names
    User: name(s) of the associated user(s)
    Library: the library type
    Single cell platform: single-cell preparation platform
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
                raise KeyError("Unrecognised field: '%s'" % field)
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
                raise KeyError("Unrecognised field: '%s'" % field)
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
    fastqs: list of the Fastq files (without leading paths)
    protocol: name of the QC protocol used
    organism: the organism(s) that the QC was run with
    seq_data_samples: samples with sequence (i.e. biological) data
    cellranger_version: version of cellranger/10x software used
    cellranger_refdata: reference datasets used with cellranger
    cellranger_probeset: probe set used with cellranger
    fastq_screens: names of panels used with fastq_screen
    star_index: index used by STAR
    annotation_bed: BED file with gene annotation
    annotation_gtf: GTF file with gene annotation
    protocol_summary: free-text summary of the QC protocol
    protocol_specification: full QC protocol specification
    """
    def __init__(self,filen=None):
        """Create a new AnalysisProjectQCDirInfo instance

        Arguments:
          filen: (optional) name of the tab-delimited file
            with key-value pairs to load in.
        """
        MetadataDict.__init__(
            self,
            attributes = {
                'fastq_dir' :'Fastq dir',
                'fastqs': 'Fastqs',
                'protocol': 'QC protocol',
                'organism': 'Organism',
                'seq_data_samples': 'Sequence data samples',
                'cellranger_version': 'Cellranger version',
                'cellranger_refdata': 'Cellranger reference datasets',
                'cellranger_probeset': 'Cellranger reference probe set',
                'fastq_screens': 'FastqScreen panels',
                'star_index': 'STAR index',
                'annotation_bed': 'BED gene annotation file',
                'annotation_gtf': 'GTF gene annotation file',
                'protocol_summary': 'Protocol summary',
                'protocol_specification': 'Protocol specification',
                'fastqs_split_by_lane': 'Fastqs split by lane',
            },
            order = (
                'fastq_dir',
                'protocol',
                'organism',
                'fastqs',
            ),
            filen=filen)


def item_to_name(item):
    """
    Convert a metadata item to its name when writing to file

    Replaces underscores with spaces, then:

     - if item is all lowercase then capitalize
     - otherwise preserve case

    For example:

    - 'order_name' converts to 'Order name'
    - 'EOL_date_stamp' converts to 'EOL date stamp'

    The operation can be reversed by calling 'name_to_item'.

    Arguments:
      item: metadata item name

    Returns:
      String: converted name for the metadata item.
    """
    name = str(item)
    if name.islower():
        # Replace underscores with spaces and capitalize
        return name.replace("_", " ").capitalize()
    else:
        # Only replace spaces
        return name.replace("_", " ")


def name_to_item(name):
    """
    Convert a metadata item name for accessing in objects

    Replaces spaces with underscores, then:

     - if starts with capital letter but is otherwise
       lowercase then convert to lowercase
     - otherwise preserve case

    For example:

    - 'Order name' converts to 'order_name'
    - 'EOL date stamp' converts to 'EOL_date_stamp'

    The operation can be reversed by calling 'item_to_name'.

    Arguments:
      name: metadata item name used in file

    Returns:
      String: converted item name for internal use.
    """
    item = str(name)
    if name[0].isupper() and all([c.islower() if c.isalpha() else True for c in name[1:]]):
        # Convert to lowercase
        name = name.lower()
    # Replace spaces with underscores
    return name.replace(" ", "_")
