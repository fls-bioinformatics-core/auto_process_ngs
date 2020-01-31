#!/usr/bin/env python
#
#     config: classes & funcs for configuring auto_process_ngs
#     Copyright (C) University of Manchester 2014-2019 Peter Briggs
#
########################################################################
#
# config.py
#
#########################################################################

"""config

Classes and functions to support configuration of the auto_process_ngs
module.

"""

from configparser import ConfigParser
from configparser import NoOptionError
from configparser import NoSectionError
from bcftbx.JobRunner import fetch_runner

#######################################################################
# Classes
#######################################################################

class Config(ConfigParser):
    """Wraps ConfigParser to set defaults for missing options

    Implements a wrapper for ConfigParser:

    - 'get' and 'getint' methods take a 'default' argument, which
      is returned if the specified option is missing from the file
    - implements a 'getrunner' method that returns a JobRunner
      instance based on a specification string.

    Example usage (read in a file and get a parameter value from a
    section):

    >>> c = Config()
    >>> c.read(conf_file)
    >>> c.get('section1','parameter2')

    See also the configparser documentation at
    https://docs.python.org/3/library/configparser.html

    """
    def __init__(self):
        ConfigParser.__init__(self)
        self.optionxform = str
    def get(self,section,option,default=None,**kwargs):
        try:
            value = super(Config,self).get(section,option,**kwargs)
            if value == 'None' or value == '':
                return default
            else:
                return value
        except NoOptionError:
            return default
        except NoSectionError:
            return default
    def getint(self,section,option,default=None):
        try:
            return super(Config,self).getint(section,option,
                                             fallback=default)
        except TypeError:
            return default
    def getfloat(self,section,option,default=None):
        try:
            return super(Config,self).getfloat(section,option,
                                               fallback=default)
        except (TypeError,AttributeError):
            return default
    def getboolean(self,section,option,default=None):
        try:
            return super(Config,self).getboolean(section,option,
                                                 fallback=default)
        except (TypeError,AttributeError):
            return default
    def getrunner(self,section,option,default='SimpleJobRunner'):
        try:
            return fetch_runner(self.get(section,option,default))
        except Exception:
            return None
