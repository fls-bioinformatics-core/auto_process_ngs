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
        # Locate settings file
        if settings_file is None:
            self.settings_file = locate_settings_file(create_from_sample=False)
        else:
            self.settings_file = os.path.abspath(settings_file)
        # Import site-specific settings from local version
        config = Config()
        if self.settings_file:
            config.read(self.settings_file)
        # General parameters
        self.general = AttributeDictionary()
        self.general['default_runner'] = config.get('general','default_runner',
                                                    'SimpleJobRunner')
        self.general['max_concurrent_jobs'] = config.getint('general',
                                                            'max_concurrent_jobs',12)
        # modulefiles
        self.modulefiles = AttributeDictionary()
        self.modulefiles['make_fastqs'] = config.get('modulefiles','make_fastqs')
        self.modulefiles['run_qc'] = config.get('modulefiles','run_qc')
        # bcl2fastq
        self.bcl2fastq = AttributeDictionary()
        self.bcl2fastq['nprocessors'] = config.getint('bcl2fastq','nprocessors',1)
        # fastq_stats
        self.fastq_stats = AttributeDictionary()
        self.fastq_stats['nprocessors'] = config.getint('fastq_stats','nprocessors',1)
        # Define runners for specific jobs
        self.runners = AttributeDictionary()
        for name in ('bcl2fastq','qc','stats',):
            self.runners[name] = config.getrunner('runners',name,
                                                  self.general.default_runner)
        # Information for archiving analyses
        # dirn should be a directory in the form [[user@]host:]path]
        self.archive = AttributeDictionary()
        self.archive['dirn'] = config.get('archive','dirn',None)
        self.archive['log'] = config.get('archive','log',None)
        self.archive['group'] = config.get('archive','group',None)
        self.archive['chmod'] = config.get('archive','chmod',None)
        # Information for uploading QC reports
        # dirn should be a directory in the form [[user@]host:]path]
        self.qc_web_server = AttributeDictionary()
        self.qc_web_server['dirn'] = config.get('qc_web_server','dirn',None)
        self.qc_web_server['url'] = config.get('qc_web_server','url',None)

    def set(self,param,value):
        """
        Update a configuration parameter value

        Arguments:
          param (str): an identifier of the form SECTION.ATTR
            which specifies the parameter to update
          value (str): the new value of the parameter

        """
        section,attr = param.split('.')
        getattr(self,section)[attr] = value

    def save(self):
        """
        Save the current configuration to the config file

        If no config file was setting on initialisation then
        this method doesn't do anything.

        """
        config = Config()
        if self.settings_file:
            config.read(self.settings_file)
            for section in config.sections():
                values = getattr(self,section)
                for attr in values:
                    config.set(section,attr,values[attr])
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
        show_dictionary('general',self.general)
        show_dictionary('modulefiles',self.modulefiles)
        show_dictionary('bcl2fastq',self.bcl2fastq)
        show_dictionary('runners',self.runners)
        show_dictionary('archive',self.archive)
        show_dictionary('qc_web_server',self.qc_web_server)

#######################################################################
# Functions
#######################################################################

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
    install_dir = os.path.abspath(os.path.normpath(
        os.path.join(os.path.dirname(sys.argv[0]),'..')))
    config_file_dirs = (os.getcwd(),
                        os.path.join(install_dir,'config'),
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
