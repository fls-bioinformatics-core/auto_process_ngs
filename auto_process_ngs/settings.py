#!/bin/env python
#
#     settings.py: handle configuration settings for autoprocessing
#     Copyright (C) University of Manchester 2014-15 Peter Briggs
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

#######################################################################
# Classes
#######################################################################

class Settings:
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
        # modulefiles
        self.add_section('modulefiles')
        self.modulefiles['make_fastqs'] = config.get('modulefiles','make_fastqs')
        self.modulefiles['run_qc'] = config.get('modulefiles','run_qc')
        # bcl2fastq
        self.add_section('bcl2fastq')
        self.bcl2fastq = self.get_bcl2fastq_config('bcl2fastq',config)
        # qc
        self.add_section('qc')
        self.qc['nprocessors'] = config.getint('qc','nprocessors',1)
        self.qc['fastq_screen_subset'] = config.getint('qc',
                                                       'fastq_screen_subset',
                                                       100000)
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
        # fastq_stats
        self.add_section('fastq_stats')
        self.fastq_stats['nprocessors'] = config.getint('fastq_stats','nprocessors',1)
        # Define runners for specific jobs
        self.add_section('runners')
        for name in ('bcl2fastq','qc','stats','rsync'):
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

    def get_bcl2fastq_config(self,section,config):
        """
        Retrieve bcl2fastq configuration options from .ini file

        Given the name of a section (e.g. 'blc2fastq',
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
        if self.settings_file:
            print "Settings from %s:" % self.settings_file
        else:
            logging.warning("No settings file found, reporting built-in defaults")
        for section in self._sections:
            if self.has_subsections(section):
                for subsection in getattr(self,section):
                    show_dictionary('%s:%s' % (section,subsection),
                                    getattr(self,section)[subsection])
            else:
                show_dictionary(section,getattr(self,section))

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

    1. current directory
    2. 'config' subdir of installation location
    3. top-level installation location

    The first file with a matching name is returned.

    If no matching file is located but one of the locations
    contains a file with the correct name ending in
    '.sample', and if the 'create_from_sample' argument is
    set, then use this to make a settings file in the same
    location.

    Returns the path to a settings file, or None if one isn't
    found.

    """
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
            except Exception,ex:
                raise Exception("Failed to create %s: %s" % (settings_file,ex))
    # Finish
    return settings_file

def show_dictionary(name,d):
    """
    Print the contents of a dictionary

    """
    print "[%s]" % name
    for key in d:
        print "\t%s = %s" % (key,(d[key] if d[key] is not None
                                  else '<Not set>'))
