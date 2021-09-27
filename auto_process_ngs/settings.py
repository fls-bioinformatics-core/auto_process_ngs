#!/bin/env python
#
#     settings.py: handle configuration settings for autoprocessing
#     Copyright (C) University of Manchester 2014-2021 Peter Briggs
#
#########################################################################
#
# settings.py
#
#########################################################################

"""
Classes and functions for handling the collection of configuration settings
for automated processing.

The settings are stored in a '.ini'-formatted file (by default called
'auto_process.ini'; this file can be created by making a copy of the
'auto_process.ini.sample' file).

The simplest usage example is:

>>> from settings import Settings
>>> s = Settings()

The values of the configuration parameters can then be accessed using
e.g.

>>> s.general.max_concurrent_jobs
4

To print the values of all parameters use

>>> s.report_settings()

To import values from a non-standard named file use e.g.

>>> s = Settings('my_auto_process.ini')

The 'locate_settings_file' function is used implicitly to locate the
settings file if none is given; it can also automatically create a settings
file if none is found but there is a sample version on the search path.

To update values once the settings have been read in do e.g.

>>> s.set('general.max_concurrent_jobs',4)

To update the configuration file use the save method e.g.

>>> s.save()

"""

#######################################################################
# Imports
#######################################################################

import os
import sys
import logging
import bcftbx.JobRunner as JobRunner
from bcftbx.utils import AttributeDictionary
from .config import Config
from .config import NoSectionError

#######################################################################
# Classes
#######################################################################

class Settings(object):
    """
    Load parameter values from an external config file

    The input file should be in '.ini' format and contain
    sections and values consistent with the sample
    'auto_process.ini' file.

    """
    def __init__(self,settings_file=None):
        """
        Create new Settings instance

        If 'settings_file' is specified then this should be the
        full path to an appropriately formatted '.ini' file.

        Otherwise the class will attempt to locate an appropriate
        file to use: by default this will be a file called
        'auto_process.ini' which will exist somewhere in the
        search path defined by the 'locate_settings_file'
        function; if no file with this name can be found then
        the class will fallback to looking for a file with the
        older 'settings.ini' file name.
        
        """
        # Initialise list of sections
        self._sections = []
        # Locate settings file
        if settings_file is None:
            # Look for default
            self.settings_file = locate_settings_file(
                name="auto_process.ini",create_from_sample=False)
            if self.settings_file is None:
                # Fallback to old name
                self.settings_file = locate_settings_file(
                    name="settings.ini",create_from_sample=False)
        else:
            self.settings_file = os.path.abspath(settings_file)
        # Import site-specific settings from local version
        config = Config()
        if self.settings_file:
            config.read(self.settings_file)
        else:
            # Look for sample settings file
            config.read(os.path.join(get_config_dir(),
                                     'auto_process.ini.sample'))
        # General parameters
        self.add_section('general')
        default_runner = config.get('general','default_runner',
                                    'SimpleJobRunner')
        self.general['default_runner'] = config.getrunner('general',
                                                          'default_runner',
                                                          'SimpleJobRunner')
        self.general['max_concurrent_jobs'] = config.getint('general',
                                                            'max_concurrent_jobs',12)
        self.general['max_cores'] = config.getint('general','max_cores')
        self.general['max_batches'] = config.getint('general','max_batches')
        self.general['poll_interval'] = config.getfloat('general',
                                                        'poll_interval',5)
        # modulefiles
        self.add_section('modulefiles')
        self.modulefiles['make_fastqs'] = config.get('modulefiles',
                                                     'make_fastqs')
        self.modulefiles['bcl2fastq'] = config.get('modulefiles','bcl2fastq')
        self.modulefiles['cellranger_mkfastq'] = config.get('modulefiles',
                                                    'cellranger_mkfastq')
        self.modulefiles['cellranger_atac_mkfastq'] = config.get('modulefiles',
                                                    'cellranger_atac_mkfastq')
        self.modulefiles['cellranger_arc_mkfastq'] = config.get('modulefiles',
                                                    'cellranger_arc_mkfastq')
        self.modulefiles['spaceranger_mkfastq'] = config.get('modulefiles',
                                                    'spaceranger_mkfastq')
        self.modulefiles['run_qc'] = config.get('modulefiles','run_qc')
        self.modulefiles['publish_qc'] = config.get('modulefiles','publish_qc')
        self.modulefiles['process_icell8'] = config.get('modulefiles','process_icell8')
        self.modulefiles['illumina_qc'] = config.get('modulefiles','illumina_qc')
        self.modulefiles['fastq_strand'] = config.get('modulefiles','fastq_strand')
        self.modulefiles['cellranger'] = config.get('modulefiles','cellranger')
        self.modulefiles['report_qc'] = config.get('modulefiles','report_qc')
        self.modulefiles['cutadapt'] = config.get('modulefiles','cutadapt')
        self.modulefiles['fastq_screen'] = config.get('modulefiles','fastq_screen')
        # conda
        self.add_section('conda')
        self.conda['enable_conda'] = config.getboolean('conda',
                                                       'enable_conda',
                                                       False)
        # bcl2fastq
        self.add_section('bcl2fastq')
        self.bcl2fastq = self.get_bcl2fastq_config('bcl2fastq',config)
        # qc
        self.add_section('qc')
        self.qc['nprocessors'] = config.getint('qc','nprocessors',None)
        self.qc['fastq_screen_subset'] = config.getint('qc',
                                                       'fastq_screen_subset',
                                                       100000)
        # Organisms
        self.add_section('organisms')
        for section in filter(lambda x: x.startswith('organism:'),
                              config.sections()):
            organism = section.split(':')[1]
            self.organisms[organism] = self.get_organism_config(
                section,config)
        # Handle legacy STAR index specifications (fastq_strand_indexes)
        try:
            for organism,index_file in config.items('fastq_strand_indexes'):
                if organism not in self.organisms:
                    self.organisms[organism] = self.get_organism_config()
                self['organisms'][organism]['star_index'] = index_file
            logging.warning("Added STAR index information from "
                            "deprecated 'fastq_strand_indexes' section (use "
                            "'organism:ORGANISM' sections instead)")
        except NoSectionError:
            pass
        # Legacy 10xgenomics transcriptome references
        try:
            for organism,reference in config.items('10xgenomics_transcriptomes'):
                if organism not in self.organisms:
                    self.organisms[organism] = self.get_organism_config()
                self['organisms'][organism]['cellranger_reference'] = reference
            logging.warning("Added cellranger references from deprecated "
                            "'10xgenomics_transcriptomes' section (use "
                            "'organism:ORGANISM' sections instead)")
        except NoSectionError:
            pass
        # Legacy 10xgenomics snRNA-seq pre-mRNA references
        try:
            for organism,reference in config.items('10xgenomics_premrna_references'):
                if organism not in self.organisms:
                    self.organisms[organism] = self.get_organism_config()
                self['organisms'][organism]['cellranger_premrna_reference'] = reference
            logging.warning("Added cellranger pre-mRNA references from "
                            "deprecated '10xgenomics_premrna_references' "
                            "section (use 'organism:ORGANISM' sections "
                            "instead)")
        except NoSectionError:
            pass
        # Legacy 10xgenomics scATAC-seq genome references
        try:
            for organism,reference in config.items('10xgenomics_atac_genome_references'):
                if organism not in self.organisms:
                    self.organisms[organism] = self.get_organism_config()
                self['organisms'][organism]['cellranger_atac_reference'] = reference
            logging.warning("Added cellranger-atac references from deprecated "
                            "'10xgenomics_atac_genome_references' section "
                            "(use 'organism:ORGANISM' sections instead)")
        except NoSectionError:
            pass
        # Legacy 10xGenomics cellranger ARC single cell multiome references
        try:
            for organism,reference in config.items('10xgenomics_multiome_references'):
                if organism not in self.organisms:
                    self.organisms[organism] = self.get_organism_config()
                self['organisms'][organism]['cellranger_arc_reference'] = reference
            logging.warning("Added cellranger-arc references from deprecated "
                            "'10xgenomics_multiome_references' section "
                            "(use 'organism:ORGANISM' sections instead)")
        except NoSectionError:
            pass
        # Sequencers
        self.add_section('sequencers')
        for section in filter(lambda x: x.startswith('sequencer:'),
                              config.sections()):
            instrument = section.split(':')[1]
            self.sequencers[instrument] = self.get_sequencer_config(
                section,config)
        # Add any settings legacy 'sequencers' section
        try:
            for instrument,platform in config.items('sequencers'):
                if instrument not in self.sequencers:
                    self['sequencers'][instrument] = \
                        AttributeDictionary(platform=None,
                                            model=None)
                self['sequencers'][instrument]['platform'] = platform
            logging.warning("Added sequencer information from "
                            "deprecated 'sequencers' section (use "
                            "'sequencer:INSTRUMENT' sections "
                            "instead)")
        except NoSectionError:
            pass
        # Sequencing platform-specific defaults
        self.add_section('platform')
        for section in filter(lambda x: x.startswith('platform:'),
                              config.sections()):
            platform = section.split(':')[1]
            self.platform[platform] = self.get_bcl2fastq_config(section,config)
        # Handle deprecated bcl2fastq settings
        for platform in ('hiseq','miseq','nextseq'):
            if config.has_option('bcl2fastq',platform):
                logging.warning("Deprecated setting in [bcl2fastq]: '%s'"
                                % platform)
            try:
                bcl2fastq = self.platform[platform]['bcl2fastq']
            except KeyError:
                bcl2fastq = config.get('bcl2fastq',platform)
                if bcl2fastq is None:
                    continue
                logging.warning("Setting 'bcl2fastq' in '[platform:%s]' to '%s'"
                                % (platform,bcl2fastq))
                if platform not in self.platform:
                    self.platform[platform] = AttributeDictionary()
                self.platform[platform]['bcl2fastq'] = bcl2fastq
        # Metadata defaults
        self.add_section('metadata')
        self.metadata['default_data_source'] = config.get('metadata',
                                                          'default_data_source')
        # icell8
        self.add_section('icell8')
        self.icell8['aligner'] = config.get('icell8','aligner')
        self.icell8['batch_size'] = config.getint('icell8','batch_size',5000000)
        self.icell8['mammalian_conf_file'] = config.get('icell8',
                                                        'mammalian_conf_file')
        self.icell8['contaminants_conf_file'] = config.get('icell8',
                                                           'contaminants_conf_file')
        self.icell8['nprocessors_contaminant_filter'] = config.getint('icell8','nprocessors_contaminant_filter',None)
        self.icell8['nprocessors_statistics'] = config.getint('icell8','nprocessors_statistics',None)
        # 10xgenomics
        self.add_section('10xgenomics')
        self['10xgenomics']['cellranger_jobmode'] = config.get('10xgenomics',
                                                               'cellranger_jobmode',
                                                               'local')
        self['10xgenomics']['cellranger_maxjobs'] = config.getint('10xgenomics','cellranger_maxjobs',24)
        self['10xgenomics']['cellranger_mempercore'] = config.getint('10xgenomics','cellranger_mempercore',5)
        self['10xgenomics']['cellranger_jobinterval'] = config.getint('10xgenomics','cellranger_jobinterval',100)
        self['10xgenomics']['cellranger_localmem'] = config.getint('10xgenomics','cellranger_localmem',5)
        self['10xgenomics']['cellranger_localcores'] = config.getint('10xgenomics','cellranger_localcores',None)
        # fastq_stats
        self.add_section('fastq_stats')
        self.fastq_stats['nprocessors'] = config.getint('fastq_stats','nprocessors',None)
        # Define runners for specific jobs
        self.add_section('runners')
        for name in ('bcl2fastq',
                     'qc',
                     'stats',
                     'rsync',
                     'icell8',
                     'icell8_contaminant_filter',
                     'icell8_statistics',
                     'icell8_report',
                     'cellranger',):
            self.runners[name] = config.getrunner('runners',name,
                                                  default_runner)
        # Information for archiving analyses
        # dirn should be a directory in the form [[user@]host:]path]
        self.add_section('archive')
        self.archive['dirn'] = config.get('archive','dirn',None)
        self.archive['log'] = config.get('archive','log',None)
        self.archive['group'] = config.get('archive','group',None)
        self.archive['chmod'] = config.get('archive','chmod',None)
        # Information for uploading QC reports
        # dirn should be a directory in the form [[user@]host:]path]
        self.add_section('qc_web_server')
        self.qc_web_server['dirn'] = config.get('qc_web_server','dirn',None)
        self.qc_web_server['url'] = config.get('qc_web_server','url',None)
        self.qc_web_server['use_hierarchy'] = config.getboolean(
            'qc_web_server','use_hierarchy')
        self.qc_web_server['exclude_zip_files'] = config.getboolean(
            'qc_web_server','exclude_zip_files')
        # Templates for reporting project data
        self.add_section('reporting_templates')
        try:
            for template,fields in config.items('reporting_templates'):
                self['reporting_templates'][template] = fields
        except NoSectionError:
            logging.debug("No reporting templates defined")
        # Destinations for data transfer
        self.add_section('destination')
        for section in filter(lambda x: x.startswith('destination:'),
                              config.sections()):
            dest = section.split(':')[1]
            self.destination[dest] = self.get_destination_config(
                section,config)

    def get_bcl2fastq_config(self,section,config):
        """
        Retrieve bcl2fastq configuration options from .ini file

        Given the name of a section (e.g. 'bcl2fastq',
        'platform:miseq'), fetch the bcl2fastq settings and return
        in an AttributeDictionary object.

        The options that can be extracted are:

        - default_version
        - bcl2fastq
        - nprocessors
        - no_lane_splitting
        - create_empty_fastqs

        Arguments:
          section (str): name of the section to retrieve the
            settings from
          config (Config): Config object with settings loaded

        Returns:
          AttributeDictionary: dictionary of option:value pairs.

        """
        values = AttributeDictionary()
        if section == 'bcl2fastq':
            values['default_version'] = config.get(section,'default_version',
                                                   None)
            values['nprocessors'] = config.getint(section,'nprocessors',None)
            values['no_lane_splitting'] = config.getboolean(section,'no_lane_splitting',
                                                            False)
            values['create_empty_fastqs'] = config.getboolean(
                section,
                'create_empty_fastqs',
                True)
        else:
            values['bcl2fastq'] = config.get(section,'bcl2fastq',None)
            values['nprocessors'] = config.getint(section,'nprocessors',None)
            values['no_lane_splitting'] = config.getboolean(section,'no_lane_splitting',
                                                            None)
            values['create_empty_fastqs'] = config.getboolean(
                section,
                'create_empty_fastqs',
                None)
        return values

    def get_destination_config(self,section,config):
        """
        Retrieve 'destination' configuration options from .ini file

        Given the name of a section (e.g. 'destination:webserver'),
        fetch the associated data transfer settings and return
        in an AttributeDictionary object.

        The options that can be extracted are:

        - directory (compulsory, str)
        - subdir (optional, str, default 'None')
        - readme_template (optional, str, default 'None')
        - url (optional, str, default 'None')
        - include_downloader (optional, boolean, default 'False')
        - include_qc_report (optional, boolean, default 'False')
        - hard_links (optional, boolean, default 'False')

        Arguments:
          section (str): name of the section to retrieve the
            settings from
          config (Config): Config object with settings loaded

        Returns:
          AttributeDictionary: dictionary of option:value pairs.

        """
        values = AttributeDictionary()
        values['directory'] = config.get(section,'directory',None)
        values['subdir'] = config.get(section,'subdir',None)
        values['readme_template'] = config.get(section,'readme_template',None)
        values['url'] = config.get(section,'url',None)
        values['include_downloader'] = config.getboolean(
            section,'include_downloader',False)
        values['include_qc_report'] = config.getboolean(
            section,'include_qc_report',False)
        values['hard_links'] = config.getboolean(section,'hard_links',False)
        return values

    def get_organism_config(self,section=None,config=None):
        """
        Retrieve 'organism' configuration options from .ini file

        Given the name of a section (e.g. 'organism:Human'),
        fetch the data association with the organism and return in
        an AttributeDictionary object.

        The items that can be extracted are:

        - star_index (str, path to STAR index)
        - bowtie_index (str, path to Bowtie index)
        - cellranger_reference (str)
        - cellranger_premrna_reference (str)
        - cellranger_atac_reference (str)
        - cellranger_arc_reference (str)

        Arguments:
          section (str): name of the section to retrieve the
            settings from
          config (Config): Config object with settings loaded

        Returns:
          AttributeDictionary: dictionary of option:value pairs.
        """
        values = AttributeDictionary()
        for param in (
                'star_index',
                'bowtie_index',
                'cellranger_reference',
                'cellranger_premrna_reference',
                'cellranger_atac_reference',
                'cellranger_arc_reference'):
            if section and config:
                values[param] = config.get(section,param,None)
            else:
                values[param] = None
        return values

    def get_sequencer_config(self,section,config):
        """
        Retrieve 'sequencer' configuration options from .ini file

        Given the name of a section (e.g. 'sequencer:SN7001250'),
        fetch the data associated with the sequencer instrument
        and return in an AttributeDictionary object.

        The items that can be extracted are:

        - platform (compulsory, str)
        - model (str, default 'None')

        Arguments:
          section (str): name of the section to retrieve the
            settings from
          config (Config): Config object with settings loaded

        Returns:
          AttributeDictionary: dictionary of option:value pairs.
        """
        values = AttributeDictionary()
        values['platform'] = config.get(section,'platform',None)
        values['model'] = config.get(section,'model',None)
        if values['platform'] is None:
            raise Exception("%s: missing required 'platform'" % section)
        if values['model']:
            # Strip quotes
            model = values['model']
            while model[0] in ('"','\'',) and model[-1] in ('"','\'',):
                model = model[1:-1]
            values['model'] = model
        return values

    def set(self,param,value):
        """
        Update a configuration parameter value

        Arguments:
          param (str): an identifier of the form
            SECTION[:SUBSECTION].ATTR which specifies the
            parameter to update
          value (str): the new value of the parameter

        """
        section,attr = param.split('.')
        try:
            section,subsection = section.split(':')
            getattr(self,section)[subsection][attr] = value
        except ValueError:
            getattr(self,section)[attr] = value

    def add_section(self,section):
        """
        Add a new section

        Arguments:
          section (str): an identifier of the form
            SECTION[:SUBSECTION] which specifies the
            section to add

        """
        try:
            section,subsection = section.split(':')
            if section not in self._sections:
                self.add_section(section)
            getattr(self,section)[subsection] = AttributeDictionary()
        except ValueError:
            self._sections.append(section)
            setattr(self,section,AttributeDictionary())

    def __getitem__(self,section):
        """
        Implement __getitem__ to enable s[SECTION]
        """
        return getattr(self,section)

    def has_subsections(self,section):
        """
        Check if section contains subsections

        Arguments:
          section (str): name of the section to check

        """
        for item in getattr(self,section):
            if not isinstance(getattr(self,section)[item],
                              AttributeDictionary):
                return False
        return True

    def save(self):
        """
        Save the current configuration to the config file

        If no config file was specified on initialisation then
        this method doesn't do anything.

        """
        config = Config()
        if self.settings_file:
            for section in self._sections:
                if not self.has_subsections(section):
                    name = section
                    values = getattr(self,section)
                    config.add_section(name)
                    for attr in values:
                        config.set(name,attr,values[attr])
                else:
                    for subsection in getattr(self,section):
                        name = "%s:%s" % (section,subsection)
                        values = getattr(self,section)[subsection]
                        config.add_section(name)
                        for attr in values:
                            config.set(name,attr,values[attr])
            config.write(open(self.settings_file,'w'))
        else:
            logging.warning("No settings file found, nothing saved")
    
    def report_settings(self):
        """
        Report the settings read from the config file
        """
        text = []
        if self.settings_file:
            text.append("Settings from %s" % self.settings_file)
        else:
            logging.warning("No settings file found, reporting built-in "
                            "defaults")
        for section in self._sections:
            if section == 'sequencers':
                display_name = 'sequencer'
            else:
                display_name = section
            if self.has_subsections(section):
                for subsection in getattr(self,section):
                    text.append(
                        show_dictionary('%s:%s' % (display_name,subsection),
                                        getattr(self,section)[subsection]))
            else:
                text.append(show_dictionary(display_name,
                                            getattr(self,section)))
        return '\n'.join(text)

#######################################################################
# Functions
#######################################################################

def get_install_dir():
    """
    Return location of top-level directory of installation

    This is a directory one or more level above the location of this
    module which contains a 'config' subdir with an
    'auto_process.ini.sample' file, for example: if this file is
    located in

    /opt/auto_process/lib/python2.7/site-packages/auto_process_ngs

    then each level will be searched until a matching 'config' dir is
    located.

    If it can't be located then the directory of this module is
    returned.

    """
    path = os.path.dirname(__file__)
    while path != os.sep:
        if os.path.isdir(os.path.join(path,'config')) and \
           os.path.isfile(os.path.join(path,'config',
                                       'auto_process.ini.sample')):
            logging.debug("Found install dir: %s" % path)
            return os.path.abspath(os.path.normpath(path))
        path = os.path.dirname(path)
    return os.path.dirname(__file__)

def get_config_dir():
    """
    Return location of config directory

    Returns the path to the 'config' directory, or None if it doesn't
    exist.

    """
    path = os.path.join(get_install_dir(),'config')
    logging.debug("Putative config dir: %s" % path)
    if os.path.isdir(path):
        return path
    else:
        return None

def locate_settings_file(name='auto_process.ini',create_from_sample=True):
    """
    Locate configuration settings file

    Look for a configuration settings file (default name
    'auto_process.ini'). The search path is:

    1. file specified by the AUTO_PROCESS_CONF environment
       variable (if it exists)
    2. current directory
    3. 'config' subdir of installation location
    4. top-level installation location

    The first file with a matching name is returned.

    If no matching file is located but one of the locations
    contains a file with the correct name ending in
    '.sample', and if the 'create_from_sample' argument is
    set, then use this to make a settings file in the same
    location.

    Returns the path to a settings file, or None if one isn't
    found.

    """
    # Check for environment variable
    try:
        settings_file = os.environ['AUTO_PROCESS_CONF']
        if os.path.exists(settings_file):
            return settings_file
    except KeyError:
        pass
    # Check locations
    install_dir = get_install_dir()
    config_dir = get_config_dir()
    config_file_dirs = (os.getcwd(),
                        config_dir,
                        install_dir,)
    settings_file = None
    sample_settings_file = None
    for path in config_file_dirs:
        settings_file = os.path.join(path,name)
        if os.path.exists(settings_file):
            # Located settings file
            break
        # No settings file here, look for a sample version
        if sample_settings_file is None:
            sample_settings_file = settings_file + '.sample'
            if not os.path.exists(sample_settings_file):
                sample_settings_file = None
        # Reset settings file to keep looking
        settings_file = None
    # No settings file found anywhere on search path
    if settings_file is None:
        logging.debug("No local settings file found in %s" %
                      ', '.join(config_file_dirs))
        if sample_settings_file is not None and create_from_sample:
            logging.warning("Attempting to make a copy from sample "
                            "settings file")
            settings_file = os.path.splitext(sample_settings_file)[0]
            try:
                with open(settings_file,'w') as fp:
                    with open(sample_settings_file,'r') as fpp:
                        fp.write(fpp.read())
                logging.warning("Created new file %s" % settings_file)
            except Exception as ex:
                raise Exception("Failed to create %s: %s" %
                                (settings_file,ex))
    # Finish
    return settings_file

def show_dictionary(name,d):
    """
    Print the contents of a dictionary
    """
    text = ["[%s]" % name]
    for key in d:
        text.append("\t%s = %s" % (key,(d[key] if d[key] is not None
                                        else '<Not set>')))
    return '\n'.join(text)
