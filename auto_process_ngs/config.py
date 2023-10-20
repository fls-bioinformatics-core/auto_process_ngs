#!/usr/bin/env python
#
#     config: utilities for managing auto_process_ngs config files
#     Copyright (C) University of Manchester 2014-2023 Peter Briggs
#
########################################################################
#
# config.py
#
#########################################################################

"""
Classes to support handling configuration files within the
'auto_process_ngs' library.

The following classes are provided:

- Config: modified application-specific version of ConfigParser
- NullValue: represents a null (unset) configuration parameter
"""

from configparser import ConfigParser
from configparser import NoOptionError
from configparser import NoSectionError
from bcftbx.JobRunner import BaseJobRunner
from bcftbx.JobRunner import fetch_runner

#######################################################################
# Classes
#######################################################################

class Config(ConfigParser):
    """
    Modified application-specific version of ConfigParser

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
        self.nullvalue = NullValue()

    def get(self,section,option,**kwargs):
        """
        Fetch value of option from config file

        Wrapper for 'get' method of superclass

        Arguments:
          section (str): section name in config file
          option (str): name of option
          default (object): value to return if option is
            to set to empty string or to string 'None'
            in config file (default: None)
          kwargs (mapping): additional keyword arguments
            to pass directly to 'get' method of superclass

        Returns:
          String: value of specified option (or default),
            or 'NullValue' instance (if section or option
            are not found and not default specified).
        """
        # Sort out whether default is in use
        use_default,default_value,kwargs = self.__use_default(kwargs)
        # Try to extract the value
        try:
            value = super(Config,self).get(section,option,**kwargs)
            if value == 'None' or value == '':
                # Explicit 'None' or empty string
                if use_default:
                    value = default_value
                else:
                    value = None
            return value
        except (NoOptionError,NoSectionError):
            # Section or option not found
            if use_default:
                return default_value
            return self.nullvalue

    def getint(self,section,option,**kwargs):
        """
        Fetch integer value option from config file

        Wrapper for 'getint' method of superclass

        Arguments:
          section (str): section name in config file
          option (str): name of option
          default (object): value to return if option is
            to set to empty string or to string 'None'
            in config file (default: None)

        Returns:
          Integer: value of specified option (or default),
            or 'NullValue' instance (if section or option
            are not found and not default specified).
        """
        # Sort out whether default is in use
        use_default,default_value,kwargs = self.__use_default(kwargs)
        if kwargs:
            raise TypeError("unexpected extra keyword arguments")
        # Extract the value
        try:
            return super(Config,self).getint(section,option)
        except TypeError:
            # NB could also indicate that value is None?
            if use_default:
                return default_value
            return self.nullvalue

    def getfloat(self,section,option,**kwargs):
        """
        Fetch float value option from config file

        Wrapper for 'getfloat' method of superclass

        Arguments:
          section (str): section name in config file
          option (str): name of option
          default (object): value to return if option is
            to set to empty string or to string 'None'
            in config file (default: None)

        Returns:
          Float: value of specified option (or default),
            or 'NullValue' instance (if section or option
            are not found and not default specified).
        """
        # Sort out whether default is in use
        use_default,default_value,kwargs = self.__use_default(kwargs)
        if kwargs:
            raise TypeError("unexpected extra keyword arguments")
        # Extract the value
        try:
            return super(Config,self).getfloat(section,option)
        except (TypeError,AttributeError):
            # NB could also indicate that value is None?
            if use_default:
                return default_value
            return self.nullvalue

    def getboolean(self,section,option,**kwargs):
        """
        Fetch boolean value option from config file

        Wrapper for 'getboolean' method of superclass

        Arguments:
          section (str): section name in config file
          option (str): name of option
          default (object): value to return if option is
            to set to empty string or to string 'None'
            in config file (default: None)

        Returns:
          Booelan: value of specified option (or default),
            or 'NullValue' instance (if section or option
            are not found and not default specified).
        """
        # Sort out whether default is in use
        use_default,default_value,kwargs = self.__use_default(kwargs)
        if kwargs:
            raise TypeError("unexpected extra keyword arguments")
        # Extract the value
        try:
            return super(Config,self).getboolean(section,option)
        except (TypeError,AttributeError):
            # E.g. if it is not recognised boolean format?
            if use_default:
                return default_value
            return self.nullvalue

    def getrunner(self,section,option,**kwargs):
        """
        Return job runner instance from option in config file

        Fetches string representation of runner from
        config file and returns a runner instance.

        Arguments:
          section (str): section name in config file
          option (str): name of option
          default (object): value to return if option is
            to set to empty string or to string 'None'
            in config file (default: None)

        Returns:
          String: value of specified option (or default),
            or 'NullValue' instance (if section or option
            are not found and not default specified).
        """
        # Sort out whether default is in use
        use_default,default_value,kwargs = self.__use_default(kwargs)
        if kwargs:
            raise TypeError("unexpected extra keyword arguments")
        # Extract the value
        runner = self.get(section,option)
        if runner is self.nullvalue or runner is None:
            if not use_default:
                return runner
            else:
                runner = default_value
        # Return existing runner instance?
        if isinstance(runner,BaseJobRunner):
            return runner
        # Get a runner instance
        return fetch_runner(runner)

    def __use_default(self,kwargs):
        # Internal: extract 'default' option from keywords
        # Return tuple:
        # - use_default: True/False (default keyword was/n't found)
        # - default_value: value of default (if found, None otherwise)
        # - kwargs: updated keywords with 'default' removed
        if 'default' in kwargs:
            use_default = True
            default_value = kwargs['default']
            kwargs = { x: kwargs[x] for x in kwargs if x != 'default' }
        else:
            use_default = False
            default_value = self.nullvalue
        return (use_default,default_value,kwargs)

class NullValue:
    """
    Represents a null (i.e. unset) value for a config setting
    """
    def __init__(self):
        pass

    def __bool__(self):
        # NullValue is always False
        return False

    def __eq__(self,x):
        # NullValue instances are all equivalent
        # to each other
        return x.__class__ == self.__class__

    def __repr__(self):
        return "<NullValue>"
