#!/usr/bin/env python
#
#     config: classes & funcs for configuring auto_process_ngs
#     Copyright (C) University of Manchester 2014-2018 Peter Briggs
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

from ConfigParser import SafeConfigParser
from ConfigParser import NoOptionError
from ConfigParser import NoSectionError
from bcftbx.JobRunner import fetch_runner

#######################################################################
# Classes
#######################################################################

class Config(SafeConfigParser):
    """Wraps SafeConfigParser to set defaults for missing options

    Implements a wrapper for SafeConfigParser:

    - 'get' and 'getint' methods take a 'default' argument, which
      is returned if the specified option is missing from the file
    - implements a 'getrunner' method that returns a JobRunner
      instance based on a specification string.

    Example usage (read in a file and get a parameter value from a
    section):

    >>> c = Config()
    >>> c.read(conf_file)
    >>> c.get('section1','parameter2')

    See also the ConfigParser documentation at
    https://docs.python.org/2/library/configparser.html

    """
    def __init__(self):
        SafeConfigParser.__init__(self)
        self.optionxform = str
    def get(self,section,option,default=None):
        try:
            value = SafeConfigParser.get(self,section,option)
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
            return SafeConfigParser.getint(self,section,option)
        except TypeError:
            return default
    def getboolean(self,section,option,default=None):
        try:
            return SafeConfigParser.getboolean(self,section,option)
        except (TypeError,AttributeError):
            return default
    def getrunner(self,section,option,default='SimpleJobRunner'):
        try:
            return fetch_runner(self.get(section,option,default))
        except Exception:
            return None
