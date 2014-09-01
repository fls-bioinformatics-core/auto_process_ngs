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
from bcftbx import JobRunner

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

#######################################################################
# Functions
#######################################################################

def fetch_runner(definition):
    """Return job runner instance based on a definition string

    Given a definition string, returns an appropriate runner
    instance.

    Definitions are of the form:

      RunnerName[(args)]

    RunnerName can be 'SimpleJobRunner' or 'GEJobRunner'.
    If '(args)' are also supplied then these are passed to
    the job runner on instantiation (only works for
    GE runners).

    """
    if definition.startswith('SimpleJobRunner'):
        return JobRunner.SimpleJobRunner(join_logs=True)
    elif definition.startswith('GEJobRunner'):
        if definition.startswith('GEJobRunner(') and definition.endswith(')'):
            ge_extra_args = definition[len('GEJobRunner('):len(definition)-1].split(' ')
            return JobRunner.GEJobRunner(ge_extra_args=ge_extra_args)
        else:
            return JobRunner.GEJobRunner()
    raise Exception,"Unrecognised runner definition: %s" % definition
