#!/bin/env python
#
#     settings.py: handle configuration settings for autoprocessing
#     Copyright (C) University of Manchester 2014-2019 Peter Briggs
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
'settings.ini'; this file can be created by making a copy of the
'settings.ini.sample' file).

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

>>> s = Settings('my_settings.ini')

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
from config import Config
from config import NoSectionError

#######################################################################
# Classes
#######################################################################

class Settings(object):
    """
    Load parameter values from an external config file

    The input file should be in '.ini' format and contain
    sections and values consistent with the sample
    settings.ini file.

    """
    def __init__(self,settings_file=None):
        """
        Create new Settings instance

        If 'settings_file' is specified then this should be the
        full path to an appropriately formatted '.ini' file.

        Otherwise the class will attempt to locate an appropriate
        file to use.
        
        """
        # Initialise list of sections
        self._sections = []
        # Locate settings file
        if settings_file is None:
            self.settings_file = locate_settings_file(create_from_sample=False)
        else:
            self.settings_file = os.path.abspath(settings_file)
        # Import site-specific settings from local version
        config = Config()
        if self.settings_file:
            config.read(self.settings_file)
        else:
            # Look for sample settings file
            config.read(os.path.join(get_config_dir(),'settings.ini.sample'))
        # General parameters
        self.add_section('general')
        default_runner = config.get('general','default_runner',
                                    'SimpleJobRunner')
        self.general['default_runner'] = config.getrunner('general',
                                                          'default_runner',
                                                          'SimpleJobRunner')
        self.general['max_concurrent_jobs'] = config.getint('general',
                                                            'max_concurrent_jobs',12)
        self.general['poll_interval'] = config.getfloat('general',
                                                        'poll_interval',5)
        # modulefiles
        self.add_section('modulefiles')
        self.modulefiles['make_fastqs'] = config.get('modulefiles','make_fastqs')
        self.modulefiles['run_qc'] = config.get('modulefiles','run_qc')
        self.modulefiles['publish_qc'] = config.get('modulefiles','publish_qc')
        self.modulefiles['process_icell8'] = config.get('modulefiles','process_icell8')
        self.modulefiles['process_10xgenomics'] = config.get('modulefiles',
                                                             'process_10xgenomics')
        self.modulefiles['illumina_qc'] = config.get('modulefiles','illumina_qc')
        self.modulefiles['fastq_strand'] = config.get('modulefiles','fastq_strand')
        self.modulefiles['cellranger'] = config.get('modulefiles','cellranger')
        self.modulefiles['report_qc'] = config.get('modulefiles','report_qc')
        # bcl2fastq
        self.add_section('bcl2fastq')
        self.bcl2fastq = self.get_bcl2fastq_config('bcl2fastq',config)
        # qc
        self.add_section('qc')
        self.qc['nprocessors'] = config.getint('qc','nprocessors',1)
        self.qc['fastq_screen_subset'] = config.getint('qc',
                                                       'fastq_screen_subset',
                                                       100000)
        # fastq_strand indexes
        self.add_section('fastq_strand_indexes')
        try:
            for genome,conf_file in config.items('fastq_strand_indexes'):
                self.fastq_strand_indexes[genome] = conf_file
        except NoSectionError:
            logging.debug("No strand stats conf files defined")
        # Sequencers
        self.add_section('sequencers')
        try:
            for instrument,platform in config.items('sequencers'):
                self['sequencers'][instrument] = platform
        except NoSectionError:
            logging.debug("No sequencers defined")
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
        self.icell8['nprocessors_contaminant_filter'] = config.getint('icell8','nprocessors_contaminant_filter',1)
        self.icell8['nprocessors_statistics'] = config.getint('icell8','nprocessors_statistics',1)
        # 10xgenomics
        self.add_section('10xgenomics')
        self['10xgenomics']['cellranger_jobmode'] = config.get('10xgenomics',
                                                               'cellranger_jobmode',
                                                               'local')
        self['10xgenomics']['cellranger_mempercore'] = config.getint('10xgenomics','cellranger_mempercore',5)
        self['10xgenomics']['cellranger_jobinterval'] = config.getint('10xgenomics','cellranger_jobinterval',100)
        self['10xgenomics']['cellranger_localmem'] = config.getint('10xgenomics','cellranger_localmem',5)
        self['10xgenomics']['cellranger_localcores'] = config.getint('10xgenomics','cellranger_localcores',1)
        # 10xgenomics transcriptomes
        self.add_section('10xgenomics_transcriptomes')
        try:
            for genome,transcriptome in config.items('10xgenomics_transcriptomes'):
                self['10xgenomics_transcriptomes'][genome] = transcriptome
        except NoSectionError:
            logging.debug("No 10xgenomics transcriptomes defined")
        # 10xgenomics snRNA-seq pre-mRNA references
        self.add_section('10xgenomics_premrna_references')
        try:
            for genome,reference in config.items('10xgenomics_premrna_references'):
                self['10xgenomics_premrna_references'][genome] = reference
        except NoSectionError:
            logging.debug("No 10xgenomics snRNA-seq pre-mRNA references defined")
        # 10xgenomics scATAC-seq genome references
        self.add_section('10xgenomics_atac_genome_references')
        try:
            for genome,reference in config.items('10xgenomics_atac_genome_references'):
                self['10xgenomics_atac_genome_references'][genome] = reference
        except NoSectionError:
            logging.debug("No 10xgenomics scATAC-seq genome references defined")
        # fastq_stats
        self.add_section('fastq_stats')
        self.fastq_stats['nprocessors'] = config.getint('fastq_stats','nprocessors',1)
        # Define runners for specific jobs
        self.add_section('runners')
        for name in ('bcl2fastq',
                     'qc',
                     'stats',
                     'rsync',
                     'icell8',
                     'icell8_contaminant_filter',
                     'icell8_statistics',
                     'icell8_report',):
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
            values['nprocessors'] = config.getint(section,'nprocessors',1)
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
            if self.has_subsections(section):
                for subsection in getattr(self,section):
                    text.append(
                        show_dictionary('%s:%s' % (section,subsection),
                                        getattr(self,section)[subsection]))
            else:
                text.append(show_dictionary(section,getattr(self,section)))
        return '\n'.join(text)

#######################################################################
# Functions
#######################################################################

def get_install_dir():
    """
    Return location of top-level directory of installation

    This is a directory one or more level above the location of this
    module which contains a 'config' subdir with a 'settings.ini.sample'
    file, for example:

    If this file is located in:

    /opt/auto_process/lib/python2.7/site-packages/auto_process_ngs

    then each level will be searched until a matching 'config' dir is
    located.

    If it can't be located then the directory of this module is
    returned.

    """
    path = os.path.dirname(__file__)
    while path != os.sep:
        if os.path.isdir(os.path.join(path,'config')) and \
           os.path.isfile(os.path.join(path,'config','settings.ini.sample')):
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

def locate_settings_file(name='settings.ini',create_from_sample=True):
    """
    Locate configuration settings file

    Look for a configuration settings file (default name
    'settings.ini'). The search path is:

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
        logging.debug("No local settings file found in %s" % ', '.join(config_file_dirs))
        if sample_settings_file is not None and create_from_sample:
            logging.warning("Attempting to make a copy from sample settings file")
            settings_file = os.path.splitext(sample_settings_file)[0]
            try:
                open(settings_file,'w').write(open(sample_settings_file,'r').read())
                logging.warning("Created new file %s" % settings_file)
            except Exception as ex:
                raise Exception("Failed to create %s: %s" % (settings_file,ex))
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
