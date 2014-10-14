#!/bin/env python
#
#     config: classes & funcs for configuring auto_process_ngs
#     Copyright (C) University of Manchester 2014 Peter Briggs
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

from ConfigParser import ConfigParser,NoOptionError
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

    """
    def __init__(self):
        ConfigParser.__init__(self)
    def get(self,section,option,default=None):
        try:
            return ConfigParser.get(self,section,option)
        except NoOptionError:
            return default
    def getint(self,section,option,default):
        try:
            return ConfigParser.getint(self,section,option)
        except NoOptionError:
            return default
    def getrunner(self,section,option,default):
        return fetch_runner(self.get(section,option,default))
