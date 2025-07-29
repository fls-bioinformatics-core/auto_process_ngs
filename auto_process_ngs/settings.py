#!/bin/env python
#
#     settings.py: handle configuration settings for autoprocessing
#     Copyright (C) University of Manchester 2014-2025 Peter Briggs
#

"""
Classes and functions for handling the collection of configuration
settings for automated processing from 'ini'-style files.

The core class is called ``GenericSettings`` and provides the base
functionality for defining and handling arbitrary configuration files.
This should be subclassed to handle specific configurations, by
defining parameters in different sections along with their types and
default values (as required).

For example, a simple ``.ini`` file could look like:

::

    [general]
    max_cores = 12
    #max_jobs = 24
    #poll_interval = 0.5

A simple settings subclass to handle this might look like:

::

    class SimpleSettings(GenericSettings):
        def __init__(self, settings_file):
            GenericSettings.__init__(
                self,
                settings = {
                    "general": { "max_cores": int,
                                 "max_jobs": int,
                                 "poll_interval": float },
                },
                defaults = {
                    "general.max_cores": 8,
                    "general.max_jobs": 24,
                    "general.poll_interval": 0.5 },
                settings_file = settings_file)

The values of the parameters can then be accessed either via a chain
of attributes, for example:

>>> s.general.max_cores
12

or via keys, for example:

>>> s["general"]["max_cores"]
12

To update values programmatically after they have been set, use the
``set`` method with the fully-qualified parameter name, for example:

>>> s.set("general.max_cores", 8)

To print the values of all parameters use the ``report_settings``
method:

>>> s.report_settings()

and to write the settings back to file use the ``save`` method:

>>> s.save()

More advanced configuration settings can also be defined, including
named and user-defined subsections, and user-defined parameters.

Named subsections, for example:

::

    [general:cores]
    max_cores = 12

    [general:jobs]
    #max_jobs = 24
    #poll_interval = 0.5

could be handled using by:

::

    class SimpleSettings(GenericSettings):
        def __init__(self, settings_file):
            GenericSettings.__init__(
                self,
                settings = {
                    "general:cores": { "max_cores": int },
                    "general:jobs": { "max_jobs": int,
                                      "poll_interval": float },
                },
                defaults = {
                    "general:cores.max_cores": 8,
                    "general:jobs.max_jobs": 24,
                    "general:jobs.poll_interval": 0.5 },
                settings_file = settings_file)

User-defined subsections can be managed using a wild card for the
section name, for example:

::

    [organism:human]
    fasta = /data/hg38.fasta

    [organism:mouse]
    fasta = /data/mm10.fasta

could be handled by:

::

    class SimpleSettings(GenericSettings):
        def __init__(self, settings_file):
            GenericSettings.__init__(
                self,
                settings = {
                    "organism:*": { "fasta": str },
                },
                settings_file = settings_file)

Alternatively a section can include user-defined parameters in a
similar way:

::

    [templates]
    default = "a,b,c"
    mine = "d,e,f"

::

    class SimpleSettings(GenericSettings):
        def __init__(self, settings_file):
            GenericSettings.__init__(
                self,
                settings = {
                    "templates": { "*": str },
                },
                settings_file = settings_file)

Optionally it is also possible to specify "fallbacks" (parameters
from which other unassigned parameters should take their values)
and legacy parameters (deprecated parameters which can be used to
supply values to other paramters, as a mechanism for backwards
compatibility).

The ``Settings`` class is built using the ``GenericSettings`` class
specifically to handle the ``auto_process.ini`` configuration file
(by default called 'auto_process.ini').

The simplest usage example is:

>>> from settings import Settings
>>> s = Settings()

The 'locate_settings_file' function is used implicitly to locate the
settings file if none is given; it can also automatically create a
settings file if none is found.
"""

#######################################################################
# Imports
#######################################################################

import os
import sys
from bcftbx.JobRunner import BaseJobRunner
from bcftbx.JobRunner import fetch_runner
from bcftbx.utils import AttributeDictionary
from fnmatch import fnmatch
from .config import Config
from .config import NullValue
from .config import NoSectionError

# Module specific logger
import logging
logger = logging.getLogger(__name__)

#######################################################################
# Classes
#######################################################################

class GenericSettings:
    """
    Base class for handling .ini configuration files

    Arguments:
      settings (dict): mapping of section names to
        dictionaries defining parameters names and types
      defaults (dict): dictionary of fully qualified
        parameter names mapped to default settings
      fallbacks (dict): dictionary of fully qualified
        parameter names mapped to fallback parameters
      legacy (dict): dictionary of fully qualified parameter
        names mapped to "legacy" fallback parameter names
      expand_vars (list): list of fully qualified parameter
        names where the values can be expanded by
        substituting environment variables
      aliases (dict): dictionary defining "aliases" for
        section names
      settings_file (str): path to an .ini format file
        to load values from
      resolve_undefined (bool): if True (default) then
        assign values to "null" parameters by checking
        fallback parameters and default values
    """
    def __init__(self, settings, defaults={}, fallbacks={},
                 legacy={}, expand_vars=[], aliases={},
                 settings_file=None, resolve_undefined=True):
        self._settings = settings
        self._defaults = defaults
        self._fallbacks = fallbacks
        self._legacy_fallbacks = legacy
        self._expand_vars = expand_vars
        self._aliases = aliases
        self._sections = []
        self._settings_file = None
        # Null value indicates no setting assigned
        self.nullvalue = NullValue()
        # Build the initial structure with all parameters
        # assigned to 'null'
        self._create_sections(self._settings)
        # Load data from file
        if settings_file:
            self._settings_file = os.path.abspath(settings_file)
            self._load_from_file(self._settings_file, self._settings)
        # Sort out legacy parameters
        self._legacy_settings = {}
        for p in self._legacy_fallbacks:
            legacy_param = self._legacy_fallbacks[p]
            section, param = legacy_param.split(".")
            if section not in self._legacy_settings:
                self._legacy_settings[section] = {}
            # Get the type from the current parameter definition
            param_type = self._fetch_parameter_type(p)
            self._legacy_settings[section][param] = param_type
        if self._legacy_settings and self._settings_file:
            # Load legacy settings into a separate object
            self._legacy = GenericSettings(
                self._legacy_settings,
                settings_file=self._settings_file,
                resolve_undefined=False)
        else:
            # No legacy settings
            self._legacy = None
        # Set defaults
        if resolve_undefined:
            self.resolve_undefined_params()

    def __getitem__(self,section):
        """
        Implement __getitem__ to enable s[SECTION]
        """
        return getattr(self, self.section_name(section))

    def __contains__(self, item):
        """
        Implement __contains__ to enable 'SECTION[:SUBSECTION][.NAME] in s'
        """
        # Break up item into section, subsection and name
        section, subsection, name = self._split_parameter(item)
        # Get the internal name for the section
        section = self.section_name(section)
        # Check for existence of item components
        if section not in self._sections:
            # Section not found
            return False
        s = getattr(self, section)
        if subsection is not None:
            if subsection not in getattr(self, section):
                # Subsection not found
                return False
            s = getattr(self, section)[subsection]
        if name is not None and name != self.nullvalue:
            # Return existence or otherwise of attr
            return name in s
        else:
            # All checks passed
            return True

    def _split_parameter(self, p):
        """
        Split paramter into section, subsection and name

        Parameter should be of the form
        ``SECTION[:SUBSECTION][.NAME]``.

        Missing elements are returned as ``None``.

        Arguments:
          p (str): parameter to split

        Returns:
          Tuple: tuple of (SECTION, SUBSECTION, NAME).
        """
        try:
            section, name = p.split('.')
        except ValueError:
            section = p
            name = None
        try:
            section, subsection = section.split(':')
        except ValueError:
            subsection = None
        return (section, subsection, name)

    def set(self, param, value):
        """
        Update a configuration parameter value

        NB parameters are referenced by the names that
        appear in the config file (rather than the
        internal representation within the Settings
        instance).

        Arguments:
          param (str): an identifier of the form
            SECTION[:SUBSECTION].NAME which specifies the
            parameter to update
          value (str): the new value of the parameter
        """
        section, subsection, name = self._split_parameter(param)
        if subsection:
            section = self.section_name(section)
            getattr(self, section)[subsection][name] = value
            logger.debug("set: %s:%s.%s -> %r" % (section,
                                                  subsection,
                                                  name,
                                                  value))
        else:
            getattr(self, section)[name] = value
            logger.debug("set: %s.%s -> %r" % (section,
                                               name,
                                               value))

    def _create_sections(self, settings, sections=None):
        """
        Create the sections defined in the settings

        Builds the initial (empty) structure with the values
        for all parameters assigned to 'null'

        Arguments:
          settings (dict):mapping of section names to
            dictionaries defining parameters names and
            types
        """
        for param in settings:
            section, subsection, name = self._split_parameter(param)
            self._add_section(section)
            if param == section:
                for var in settings[param]:
                    self[section][var] = self.nullvalue

    def _add_section(self, section):
        """
        Add a new section

        Arguments:
          section (str): an identifier of the form
            SECTION[:SUBSECTION] which specifies the
            section to add
        """
        section, subsection, name = self._split_parameter(section)
        if subsection:
            section = self.section_name(section)
            if section not in self._sections:
                self._add_section(section)
            if subsection not in self[section]:
                getattr(self, section)[subsection] = AttributeDictionary()
        else:
            if section not in self._sections:
                self._sections.append(section)
                setattr(self, section, AttributeDictionary())
            for alias in self._aliases:
                # FIXME only needs to be done once when
                # FIXME section is first added?
                if section == self._aliases[alias]:
                    setattr(self, alias, self[section])

    def _load_from_file(self, settings_file, settings):
        """
        Load settings data from .ini file

        Arguments:
          settings_file (str): path to the .ini file
          settings (dict):mapping of section names to
            dictionaries defining parameters names and
            types
        """
        logger.debug(f"Loading values from '{settings_file}'")
        config = Config()
        self.nullvalue = config.nullvalue
        config.read(settings_file)
        for name in settings:
            logger.debug(f"- loading data for section '{name}'")
            if name.endswith(":*"):
                # Handle pseudo-sections e.g. [organism:human]
                section = name.split(":")[0]
                for conf in filter(lambda x: x.startswith(name[:-1]),
                                   config.sections()):
                    # Look up matching settings in the data file
                    subsection = conf.split(":")[1]
                    if subsection not in self[section]:
                        # Add new subsection
                        self[section][subsection] = AttributeDictionary()
                    for var in settings[name]:
                        # For all parameters defined in the settings
                        # instance, look for matching values in the
                        # configuration file
                        type = settings[name][var]
                        self[section][subsection][var] = self.update_value(
                            config.get(conf, var), type)
            else:
                # Handle standard sections with no subsections
                for var in settings[name]:
                    # For all parameters defined in the settings
                    # instance, look for matching values in the
                    # configuration file
                    logger.debug(f"-- locating data for parameter '{var}'")
                    type = settings[name][var]
                    if var != "*":
                        logger.debug(f"-- setting '{name}.{var}'")
                        self[name][var] = self.update_value(
                            config.get(name, var), type)
                    else:
                        # Special case: new settings are added
                        # for all parameters in the equivalent
                        # section in the configuration file
                        try:
                            for param, value in config.items(name):
                                getattr(self, name)[param] = \
                                    self.update_value(value)
                        except NoSectionError:
                            logger.debug(f"{name}: no section in config file")

    def has_subsections(self, section):
        """
        Check if section contains subsections

        Arguments:
          section (str): name of the section to check
        """
        for item in getattr(self, section):
            if not isinstance(getattr(self, section)[item],
                              AttributeDictionary):
                return False
        return True

    def fetch_value(self, param):
        """
        Return the value stored against a parameter

        NB parameters are referenced by the names that
        appear in the config file (rather than the
        internal representation within the Settings
        instance).

        Arguments:
          param (str): parameter name of the form
            SECTION[:SUBSECTION][.NAME]
        """
        # Break up param into section, subsection and attr
        section, subsection, name = self._split_parameter(param)
        # Get the internal name for the section
        section = self.section_name(section)
        if subsection is None:
            s = getattr(self,section)
        else:
            try:
                s = getattr(self,section)[subsection]
            except KeyError as ex:
                print(f"Exception fetching '{param}': {ex}")
                raise ex
        return s[name]

    def list_params(self,pattern=None,exclude_undefined=False):
        """
        Return (yield) all the stored parameters

        NB parameters are referenced by the names that
        appear in the config file (rather than the
        internal representation within the Settings
        instance).

        Arguments:
          pattern (str): optional glob-style pattern;
            if supplied then only parameters matching
            the pattern will be returned
          exclude_undefined (bool): if True then only
            undefined parameters (i.e. those with null
            values) will be returned

        Yields:
          String: parameter names of the form
            SECTION[:SUBSECTION][.NAME]
        """
        if pattern:
            if '.' not in pattern:
                pattern += '.*'
        for section in self._sections:
            if not self.has_subsections(section):
                name = self.section_name(section)
                values = getattr(self, section)
                for v in values:
                    param = "%s.%s" % (name, v)
                    if pattern and not fnmatch(param, pattern):
                        continue
                    if exclude_undefined and \
                       self.fetch_value(param) == self.nullvalue:
                        continue
                    else:
                        yield param
            else:
                subsections = getattr(self, section)
                for subsection in subsections:
                    name = f"{self.section_name(section)}:{subsection}"
                    values = getattr(self,section)[subsection]
                    for v in values:
                        param = "%s.%s" % (name, v)
                        if pattern and not fnmatch(param,pattern):
                            continue
                        if exclude_undefined and \
                           self.fetch_value(param) == self.nullvalue:
                            continue
                        else:
                            yield param

    def resolve_undefined_params(self):
        """
        Set non-null values for all parameters which are null

        Resolution comprises the following stages:

        - Legacy parameters: if a parameter is unset (i.e. null)
          and has a "legacy" fallback defined, then it will be
          set to the value of the legacy parameter;
        - Fallback parameters: if a parameter is unset and
          has a standard fallback defined, then it will be set
          to the value of the fallback;
        - Default paramaters: if a parameter is unset after the
          fallbacks are exhausted but has a defined default, then
          it will be set to the default value.

        A final round of checking fallback parameters is then
        performed, in cases fallbacks that were previously unset
        have been assigned a default value.

        Any parameters that are still undefined at this point
        are then assigned the value of 'None'.
        """
        logger.debug(f"Resolving undefined values")
        # Legacy fallback resolution
        if self._legacy:
            self._resolve_fallbacks(self._legacy_fallbacks,
                                    self._legacy)
        # Initial round of fallback resolution
        self._resolve_fallbacks(self._fallbacks)
        # Populate incomplete subsections
        self._populate_subsections()
        # Set defaults
        self._resolve_defaults(self._defaults)
        # Second round of fallback resolution
        # (in case fallbacks which were previously unset
        # are now set to their defaults)
        self._resolve_fallbacks(self._fallbacks)
        # Set remaining undefined parameters to 'None'
        for param in self.list_params():
            if self.fetch_value(param) == self.nullvalue:
                logger.debug(f"Updating undefined parameter '{param}' "
                             "to None")
                self.set(param, None)
        # Expand variables
        for param in self._expand_vars:
            value = self.fetch_value(param)
            if value:
                self.set(param,os.path.expandvars(value))

    def _resolve_fallbacks(self, fallbacks, source_settings=None):
        """
        Set values for unset parameters to fallback values

        Arguments:
          fallbacks (dict): dictionary mapping parameters
            to their fallbacks
          source_settings (GenericSettings): source to get
            fallback values from (default: self)
        """
        logger.debug("Resolving fallbacks")
        if source_settings is None:
            source_settings = self
        logger.debug(f"Fallbacks: {fallbacks}")
        logger.debug(f"Source settings: {list(source_settings.list_params())}")
        for param in fallbacks:
            fallback_param = fallbacks[param]
            if not ("*" in fallback_param and
                    "*" in param):
                # Isn't a transformation, skip
                continue
            logger.debug(f"- fallback for '{param}' is '{fallback_param}'")
            # Expand the wildcard to get a list of actual
            # fallback parameters
            for p in source_settings.list_params(fallback_param):
                if "*" in p:
                    # Ignore unexpanded wildcards
                    # FIXME can this happen?
                    continue
                # Transform the fallback parameter into the
                # actual parameter
                transformed_param = self.transform_parameter(p,
                                                             fallback_param,
                                                             param)
                logger.debug(f"- transformed '{p}' to '{transformed_param}'")
                if transformed_param not in self or \
                   self.fetch_value(transformed_param) == self.nullvalue:
                    # Ensure the subsection is present
                    section, subsection, name = self._split_parameter(
                        transformed_param)
                    self._add_section(f"{section}:{subsection}")
                    # Add/update the parameter
                    value = source_settings.fetch_value(p)
                    self.set(transformed_param, value)
                    logger.debug(
                        f"- '{transformed_param}' now set to "
                        f"'{self[self.section_name(section)][subsection][name]}'")
        # Deal with standard fallbacks
        fallback_params = {}
        for param in fallbacks:
            # Skip parameters which require transformation
            if "*" in fallbacks[param] and "*" in param:
                # Isn't a standard fallback
                continue
            # Make expanded set of fallback params by resolving
            # any wildcards
            for p in self.list_params(param):
                fallback_params[p] = fallbacks[param]
        for param in fallback_params:
            if param not in self or self.fetch_value(param) == self.nullvalue:
                # Handle specific parameter
                logger.debug(f"- checking fallbacks for {param}")
                tried_params = set()
                fallback_param = fallback_params[param]
                while fallback_param and fallback_param not in tried_params:
                    try:
                        value = source_settings.fetch_value(fallback_param)
                    except AttributeError:
                        logger.warning(f"Fallback '{fallback_param}' not "
                                       f"found")
                        value = self.nullvalue
                    if value != self.nullvalue:
                        self.set(param, value)
                        logger.debug(f"- '{param}' now set to "
                                     f"'{self.fetch_value(param)}'")
                        fallback_param = None
                    elif fallback_param in fallback_params:
                        tried_params.add(fallback_param)
                        fallback_param = fallback_params[fallback_param]
                    else:
                        fallback_param = None

    def _resolve_defaults(self, defaults):
        """
        Set values of unset parameters to their defaults

        Arguments:
          defaults (dict): dictionary mapping parameters
            to their default values
        """
        for param in defaults:
            for p in self.list_params(pattern=param):
                if self.fetch_value(p) == self.nullvalue:
                    default_value = defaults[param]
                    logger.debug(f"updating '%s' with default value %r" %
                                 (p, default_value))
                    param_type = self._fetch_parameter_type(p)
                    self.set(p, self.update_value(default_value,
                                                  param_type))

    def _populate_subsections(self):
        """
        Adds missing parameters to incomplete wildcard sections

        Updates sections defined with wildcard subsections
        (e.g. 'organism:*:...') by adding any missing
        parameters.
        """
        # Fully populate wildcard sections
        logger.debug("Populating incomplete wildcard sections")
        for s in self._settings:
            section, subsection, name = self._split_parameter(s)
            if subsection == "*":
                section = s.split(":")[0]
                for name in self[section]:
                    if name == "*":
                        # Skip wildcard name
                        continue
                    for var in self._settings[f"{section}:*"]:
                        if f"{section}:{name}.{var}" not in self:
                            logger.debug(f"- adding missing "
                                         f"{section}:{name}.{var}")
                            self.set(f"{section}:{name}.{var}",
                                     self.nullvalue)

    def _fetch_parameter_type(self, param):
        """
        Fetch the type assigned to a parameter
        """
        section, subsection, name = self._split_parameter(param)
        if subsection:
            try:
                return self._settings[f"{section}:{subsection}"][name]
            except KeyError:
                return self._settings[f"{section}:*"][name]
        else:
             return self._settings[section][name]

    def save(self, out_file=None, exclude_undefined=True):
        """
        Save the current configuration to the config file

        If no config file was specified on initialisation then
        this method doesn't do anything.

        Arguments:
          out_file (str): specify output file (default:
            overwrite initial config file)
          exclude_undefined (bool): if True then parameters
            with null values will not be written to the
            output config file (default)
        """
        if not out_file:
            out_file = self._settings_file
        if out_file:
            out_file = os.path.abspath(out_file)
            config = Config()
            # Output sections with defined parameters
            for param in self.list_params(exclude_undefined=exclude_undefined):
                section, subsection, name = self._split_parameter(param)
                if name == "*":
                    # Skip wild-card definitions
                    continue
                if subsection:
                    section = f"{self.section_name(section)}:{subsection}"
                else:
                    section = self.section_name(section)
                if not config.has_section(section):
                    config.add_section(section)
                value = self.fetch_value(param)
                if exclude_undefined and value == self.nullvalue:
                    continue
                config.set(section, name, str(value))
            # Ensure user-defined sections (i.e. '[section:name]')
            # are output, even if they don't contain any defined
            # parameters
            for section in self._sections:
                if self.has_subsections(section):
                    for subsection in getattr(self, section):
                        name = f"{self.section_name(section)}:{subsection}"
                        if not config.has_section(name):
                            config.add_section(name)
            with open(out_file,'wt') as fp:
                config.write(fp)
        else:
            logger.warning("No output file, nothing saved")

    def report_settings(self, exclude_undefined=False):
        """
        Report the settings read from the config file

        Arguments:
          exclude_undefined (bool): if True then parameters
            with null values will not be shown (default: show
            all parameters)

        Returns:
          String: report of the settings.
        """
        if exclude_undefined:
            exclude_value = self.nullvalue
        else:
            exclude_value = None
        text = []
        if self._settings_file:
            text.append("Settings from %s" % self._settings_file)
        else:
            logger.warning("No settings file found, reporting built-in "
                           "defaults")
        reported = set()
        for param in self.list_params(exclude_undefined=exclude_undefined):
            section, subsection, name = self._split_parameter(param)
            if subsection:
                # Handle subsection
                sname = f"{section}:{subsection}"
                if sname in reported:
                    # Already reported, skip
                    continue
                # Generate content
                content = show_dictionary(
                        getattr(self, section)[subsection],
                        exclude_value=exclude_value)
                if content:
                    text.append(f"[{sname}]")
                    text.append(content)
                # Add to set of reported sections
                reported.add(sname)
            else:
                # Handle standard section
                sname = f"{section}"
                if sname in reported:
                    # Already reported, skip
                    continue
                # Generate content
                content = show_dictionary(
                    getattr(self, section),
                    exclude_value=exclude_value)
                if content:
                    text.append(f"[{sname}]")
                    text.append(content)
                # Add to set of reported sections
                reported.add(sname)
        return '\n'.join(text)

    def update_value(self, value, param_type=None):
        """
        Update raw value by stripping quotes and converting to type

        'None' values are converted to "null"; other non-null
        values will have any surrounding quotes removed and
        then are converted to type.

        Boolean values will also have 'yes' and 'True' converted
        to True and 'no' and 'False' converted to False.

        Arguments:
          value (object): raw value
          param_type (function): type conversion function
        """
        if value is None:
            # Set None to null
            value = self.nullvalue
        elif value != self.nullvalue:
            # Strip quotes
            if str(value)[0] in ("\"'"):
                if str(value)[0] == str(value)[-1]:
                    value = str(value)[1:-1]
            # Handle non-null values
            if param_type is bool:
                # Special handling of boolean types
                if str(value).lower() in ("true", "yes"):
                    value = True
                elif str(value).lower() in ("false", "no"):
                    value = False
                else:
                    raise TypeError(f"{value}: invalid value for boolean")
            elif param_type:
                # All other types
                value = param_type(value)
        return value

    def section_name(self, name):
        """
        Return the internal name for a section
        """
        try:
            return self._aliases[name]
        except KeyError:
            return name

    def transform_parameter(self, param, source_pattern, target_pattern):
        """
        Transforms parameter by mapping from a source to a target pattern

        Given a parameter, map this to another parameter using a
        source pattern (e.g. ``section.*``) and a target pattern
        (e.g. ``new_section:*.value``).

        The value in the ``*`` wildcard position in the target is
        replaced by that in the wildcard position from the source
        (e.g. transforming ``section.name`` to
        ``new_section:name.value``).

        Arguments:
          param (str): fully-qualified source paramter name to
            transform
          source_pattern (str): pattern for source parameter
          target_pattern (str): pattern to use for transformation

        Returns:
           String: transformed parameter name.
        """
        # Break the parameter up into components
        section, param = param.split(".")
        try:
            section, subsection = section.split(":")
        except ValueError:
            subsection = None
        # Identify wildcard in source pattern
        if ":*." in source_pattern:
            # Matches subsection
            sub = subsection
        elif ".*" in source_pattern:
            # Matches trailing parameter
            sub = param
        # Substitute wildcard into the target_pattern
        section, param = target_pattern.split(".")
        if ":*" in section:
            section = section.split(":")[0]
            return f"{section}:{sub}.{param}"
        elif param == "*":
            return f"{section}.{sub}"


class Settings(GenericSettings):
    """
    Handle local settings for ``auto_process`` parameters

    Defines a set of configuration parameters and provides
    an interface for loading and accessing local settings
    defined in an ``.ini`` file.

    If a config file isn't explicitly specified then the
    instance will will attempt to locate one by searching
    various locations (as defined within the
    ``locate_settings_file`` function) first using the
    name ``auto_process.ini``, then with the legacy name
    ``settings.ini``.

    If a config file still cannot be found then the
    parameters will be set to any default values defined
    within the ``Settings`` class.

    Arguments:
      settings_file (str): optional, path to .ini file to
        load parameters from
      resolve_undefined (bool): if True (default) then
        assign values to "null" parameters by checking
        fallback parameters and default values
    """

    def __init__(self, settings_file=None, resolve_undefined=True):
        GenericSettings.__init__(
            self,
            # Define the sections, parameters and types
            settings = {
                "general": { "default_runner": jobrunner,
                             "max_concurrent_jobs": int,
                             "max_cores": int,
                             "max_batches": int,
                             "poll_interval": float },
                "modulefiles": {
                    'make_fastqs': str,
                    'bcl2fastq': str,
                    'bcl_convert': str,
                    'cellranger_mkfastq': str,
                    'cellranger_atac_mkfastq': str,
                    'cellranger_arc_mkfastq': str,
                    'spaceranger_mkfastq': str,
                    'run_qc': str,
                    'publish_qc': str,
                    'process_icell8': str,
                    'fastqc': str,
                    'fastq_screen': str,
                    'fastq_strand': str,
                    'cellranger': str,
                    'report_qc': str,
                    'cutadapt': str },
                "conda": { "enable_conda": bool,
                           "env_dir": str },
                "bcl_conversion": { "bcl_converter": str,
                                    "nprocessors": int,
                                    "no_lane_splitting": bool,
                                    "create_empty_fastqs": bool,
                                    "use_conda_for_bcl2fastq": bool },
                "qc": { "nprocessors": int,
                        "fastq_screens": str,
                        "fastq_subset_size": int,
                        "split_undetermined_fastqs": bool,
                        "use_legacy_screen_names": bool },
                "screen:*": { "conf_file": str },
                "organism:*": { "star_index": str,
                                "bowtie_index": str,
                                "cellranger_reference": str,
                                "annotation_bed": str,
                                "annotation_gtf": str,
                                "cellranger_premrna_reference": str,
                                "cellranger_atac_reference": str,
                                "cellranger_arc_reference": str,
                                "cellranger_probe_set": str },
                "sequencer:*": { "platform": str,
                                 "model": str },
                "platform:*": {  "bcl_converter": str,
                                 "nprocessors": int,
                                 "no_lane_splitting": bool,
                                 "create_empty_fastqs": bool },
                "metadata": { "default_data_source": str },
                "icell8" : { "aligner": str,
                             "batch_size": int,
                             "mammalian_conf_file": str,
                             "contaminants_conf_file": str,
                             "nprocessors_contaminant_filter": int,
                             "nprocessors_statistics": int },
                "10xgenomics": { "cellranger_jobmode": str,
                                 "cellranger_maxjobs": int,
                                 "cellranger_mempercore": int,
                                 "cellranger_jobinterval": int,
                                 "cellranger_localmem": int,
                                 "cellranger_localcores": int },
                "fastq_stats": { "nprocessors": int },
                "runners": { 'barcode_analysis': jobrunner,
                             'bcl2fastq': jobrunner,
                             'bcl_convert': jobrunner,
                             'cellranger': jobrunner,
                             'cellranger_count': jobrunner,
                             'cellranger_mkfastq': jobrunner,
                             'cellranger_multi': jobrunner,
                             'fastqc': jobrunner,
                             'fastq_screen': jobrunner,
                             'icell8': jobrunner,
                             'icell8_contaminant_filter': jobrunner,
                             'icell8_statistics': jobrunner,
                             'icell8_report': jobrunner,
                             'merge_fastqs': jobrunner,
                             'picard': jobrunner,
                             'publish_qc': jobrunner,
                             'qc': jobrunner,
                             'qualimap': jobrunner,
                             'rseqc': jobrunner,
                             'rsync': jobrunner,
                             'star': jobrunner,
                             'stats': jobrunner },
                "archive": { "dirn": str,
                             "log": str,
                             "group": str,
                             "chmod": str },
                "qc_web_server": { "dirn": str,
                                   "url": str,
                                   "use_hierarchy": bool,
                                   "exclude_zip_files": bool },
                "reporting_templates": { "*": str },
                "destination:*": { "directory": str,
                                   "subdir": str,
                                   "zip_fastqs": bool,
                                   "max_zip_size": str,
                                   "readme_template": str,
                                   "url": str,
                                   "include_downloader": bool,
                                   "include_qc_report": bool,
                                   "hard_links": bool }
            },
            # Aliases for sections
            aliases = { "organisms": "organism",
                        "screens": "screen",
                        "sequencers": "sequencer" },
            # Default values
            defaults =  {
                "general.default_runner": "SimpleJobRunner",
                "general.max_concurrent_jobs": 12,
                "general.poll_interval": 5,
                "conda.enable_conda": False,
                "bcl_conversion.bcl_converter": "bcl2fastq>=2.20",
                "bcl_conversion.no_lane_splitting": False,
                "bcl_conversion.create_empty_fastqs": False,
                "bcl_conversion.use_conda_for_bcl2fastq": False,
                "destination:*.include_qc_report": False,
                "destination:*.hard_links": False,
                "destination:*.zip_fastqs": False,
                "destination:*.include_downloader": False,
                "icell8.batch_size": 5000000,
                "qc.fastq_subset_size": 100000,
                "qc.split_undetermined_fastqs": True,
                "qc.use_legacy_screen_names": False,
                "qc_web_server.use_hierarchy": False,
                "qc_web_server.exclude_zip_files": False,
                "10xgenomics.cellranger_jobmode": "local",
                "10xgenomics.cellranger_maxjobs": 24,
                "10xgenomics.cellranger_mempercore": 5,
                "10xgenomics.cellranger_jobinterval": 100,
                "10xgenomics.cellranger_localmem": 5,
                "10xgenomics.cellranger_localcores": 1 },
            # Fallbacks
            fallbacks = {
                "runners.*": "general.default_runner" },
            # Legacy fallbacks
            legacy = {
                # Legacy bcl conversion parameters
                "bcl_conversion.bcl_converter": "bcl2fastq.default_version",
                "bcl_conversion.nprocessors": "bcl2fastq.nprocessors",
                "bcl_conversion.no_lane_splitting":
                "bcl2fastq.no_lane_splitting",
                "bcl_conversion.create_empty_fastqs":
                "bcl2fastq.create_empty_fastqs",
                # Legacy module files
                "modulefiles.fastqc": "modulefiles.illumina_qc",
                "modulefiles.fastq_screen": "modulefiles.illumina_qc",
                # Legacy platform settings
                "platform:*.bcl_converter": "platform:*.bcl2fastq",
                # Legacy QC subsetting
                "qc.fastq_subset_size": "qc.fastq_screen_subset",
                # Legacy organism-specific reference data
                "organism:*.star_index": "fastq_strand_indexes.*",
                "organism:*.cellranger_reference":
                "10xgenomics_transcriptomes.*",
                "organism:*.cellranger_premrna_reference":
                "10xgenomics_premrna_references.*",
                "organism:*.cellranger_atac_reference":
                "10xgenomics_atac_genome_references.*",
                "organism:*.cellranger_arc_reference":
                "10xgenomics_multiome_references.*",
                # Legacy sequencers
                "sequencer:*.platform": "sequencers.*",
                # Legacy runners
                "runners.cellranger_mkfastq": "runners.cellranger",
                "runners.cellranger_count": "runners.cellranger",
                "runners.cellranger_multi": "runners.cellranger_count",
                "runners.fastqc": "runners.qc",
                "runners.fastq_screen": "runners.qc",
                "runners.picard": "runners.qc",
                "runners.qualimap": "runners.qc",
                "runners.rseqc": "runners.qc",
                "runners.star": "runners.qc" },
            # Parameters where variables can be expanded
            # e.g. $HOME/auto_process -> /home/user/auto_process
            expand_vars = ["conda.env_dir"],
            settings_file=self._find_config(settings_file),
            resolve_undefined=resolve_undefined)
        # Post-process
        self._check_sequencer_platforms()
        self._update_default_bcl_converter()
        self._update_platform_bcl_converter()

    def _find_config(self, settings_file=None):
        """
        Locate settings file to load values from

        Arguments:
          settings_file (str): if set then use this
            as the settings file; otherwise attempt
            locate a file by searching the default
            locations
        """
        if settings_file is None:
            # Look for default
            settings_file = locate_settings_file(
                name="auto_process.ini",
                create_from_sample=False)
            if settings_file is None:
                # Fallback to old name
                settings_file = locate_settings_file(
                    name="settings.ini",
                    create_from_sample=False)
        return settings_file

    def _check_sequencer_platforms(self):
        """
        Check 'platform' set for all defined sequencers

        Raises ``Exception`` if any sequencer doesn't have
        a defined ``platform`` parameter.
        """
        for sequencer in self.sequencers:
            if not self.sequencers[sequencer].platform:
                raise Exception(f"{sequencer}: platform not set")

    def _update_default_bcl_converter(self):
        """
        Update default 'bcl_conversion.bcl_converter'

        If the default BCL converter was set from a legacy
        fallback then it will probably be of the form e.g.
        ">=2.0". If so then the value will be updated to
        prepend the value with ``bcl2fastq``.
        """
        default_bcl_converter = self.bcl_conversion.bcl_converter
        if default_bcl_converter and default_bcl_converter[0] in ("><="):
            self.set(f"bcl_conversion.bcl_converter",
                     f"bcl2fastq{default_bcl_converter}")

    def _update_platform_bcl_converter(self):
        """
        Update 'platform:*.bcl_converter' for legacy platforms

        If the BCL converter for any platform was set from a
        legacy fallback then it will probably be of the form
        e.g. ">=2.0". If so then the value will be updated to
        prepend the value with ``bcl2fastq``.
        """
        for platform in self.platform:
            bcl_converter = self.platform[platform].bcl_converter
            if bcl_converter and bcl_converter[0] in ("><="):
                self.set(f"platform:{platform}.bcl_converter",
                         f"bcl2fastq{bcl_converter}")


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

    /opt/auto_process/lib/python3.6/site-packages/auto_process_ngs

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
            logger.debug("Found install dir: %s" % path)
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
    logger.debug("Putative config dir: %s" % path)
    if os.path.isdir(path):
        return path
    else:
        return None

def locate_settings_file(name='auto_process.ini',create_from_sample=False):
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
        logger.debug("No local settings file found in %s" %
                     ', '.join(config_file_dirs))
        if sample_settings_file is not None and create_from_sample:
            logger.warning("Attempting to make a copy from sample "
                           "settings file")
            settings_file = os.path.splitext(sample_settings_file)[0]
            try:
                with open(settings_file,'w') as fp:
                    with open(sample_settings_file,'r') as fpp:
                        fp.write(fpp.read())
                logger.warning("Created new file %s" % settings_file)
            except Exception as ex:
                raise Exception("Failed to create %s: %s" %
                                (settings_file,ex))
    # Finish
    return settings_file

def fetch_reference_data(s,name):
    """
    Fetch specific reference data for all organisms

    Given the name of a reference data item (e.g.
    'star_index'), extracts the relevant data
    for all organisms from the supplied settings
    object and returns it as a dictionary.

    Arguments:
      s (Settings): populated Settings instance
      name (str): name of the reference data
        items required (e.g. 'star_index')

    Returns:
      Dictionary: keys are organism names and
        values are the corresponding reference
        data.
    """
    refdata = dict()
    for organism in s.organisms:
        data_item = s.organisms[organism][name]
        print(f"fetch_reference_data: {organism}/{name} = {data_item}")
        print("%r" % data_item)
        if data_item:
            refdata[organism] = data_item
    return refdata

def show_dictionary(d,indent='   ',exclude_value=None,
                    exclude_wildcards=True):
    """
    Print the contents of a dictionary

    Arguments:
      d (str): dictionary instance to show
      exclude_value (object): optional, if not 'None'
        then don't include entries which match this
        value
      exclude_wildcards (bool): optional, if True
        (default) then don't include keys which are
        wildcards (i.e. "*")
    """
    text = []
    for key in d:
        if exclude_wildcards and key == "*":
            continue
        if exclude_value is not None and d[key] is exclude_value:
            continue
        text.append("%s%s = %s" % (indent,
                                   key,
                                   (d[key] if d[key] is not None
                                    else '<Not set>')))
    return '\n'.join(text)

def jobrunner(runner):
    """
    Implement a 'type' function for JobRunner instances

    Arguments:
      runner (object): JobRunner instance or string
        definition

    Returns:
      JobRunner: instance created from the supplied
        definition.
    """
    if not isinstance(runner, BaseJobRunner):
        return fetch_runner(runner)
    else:
        return runner
