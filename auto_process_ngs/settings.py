#!/bin/env python
#
#     settings.py: handle configuration settings for autoprocessing
#     Copyright (C) University of Manchester 2014-2023 Peter Briggs
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
from bcftbx.JobRunner import fetch_runner
from bcftbx.utils import AttributeDictionary
from fnmatch import fnmatch
from .config import Config
from .config import NoSectionError

# Module specific logger
import logging
logger = logging.getLogger(__name__)

#######################################################################
# Classes
#######################################################################

class Settings:
    """
    Load parameter values from an external config file

    The input file should be in '.ini' format and contain
    sections and values consistent with the sample
    'auto_process.ini' file.

    """
    # Look up table mapping names of sections used internally
    # within the Settings class to the names that appear in
    # the config file, when these differ
    __SECTIONS_INTERNAL_TO_CONFIG_NAMES = {
        "sequencers": "sequencer",
        "organisms": "organism",
        "screens": "screen"
    }

    # Runners
    __RUNNERS = (
        'barcode_analysis',
        'bcl2fastq',
        'bcl_convert',
        'cellranger',
        'cellranger_count',
        'cellranger_mkfastq',
        'cellranger_multi',
        'fastqc',
        'fastq_screen',
        'icell8',
        'icell8_contaminant_filter',
        'icell8_statistics',
        'icell8_report',
        'merge_fastqs',
        'picard',
        'publish_qc',
        'qc',
        'qualimap',
        'rseqc',
        'rsync',
        'star',
        'stats',
    )

    # Environment modules
    __ENVIRONMENT_MODULES = (
        'make_fastqs',
        'bcl2fastq',
        'bcl_convert',
        'cellranger_mkfastq',
        'cellranger_atac_mkfastq',
        'cellranger_arc_mkfastq',
        'spaceranger_mkfastq',
        'run_qc',
        'publish_qc',
        'process_icell8',
        'fastqc',
        'fastq_screen',
        'fastq_strand',
        'cellranger',
        'report_qc',
        'cutadapt',
        'illumina_qc',
    )

    # Default values
    # These will be used if the parameter was not
    # explicitly defined in the config file and
    # there is no fallback (or the fallback was
    # also not defined)
    __DEFAULTS = {
        "general.default_runner": fetch_runner("SimpleJobRunner"),
        "general.max_concurrent_jobs": 12,
        "general.poll_interval": 5,
        "conda.enable_conda": False,
        "qc.fastq_subset_size": 100000,
        "qc.use_legacy_screen_names": False,
        "10xgenomics.cellranger_jobmode": "local",
        "10xgenomics.cellranger_maxjobs": 24,
        "10xgenomics.cellranger_mempercore": 5,
        "10xgenomics.cellranger_jobinterval": 100,
        "10xgenomics.cellranger_localmem": 5,
        "10xgenomics.cellranger_localcores": 1,
        "destination:*.include_zip_fastqs": False,
        "destination:*.include_downloader": False,
        "destination:*.include_qc_report": False,
        "destination:*.hard_links": False,
        "icell8.batch_size": 5000000,
    }

    # Parameters where variables can be expanded
    # e.g. $HOME/auto_process -> /home/user/auto_process
    __EXPAND_VARS = (
        "conda.env_dir",
    )

    # Fallback parameters
    # These will be used if the parameter was not
    # explicitly defined in the config file but the
    # fallback parameter was
    __FALLBACKS = {
        "qc.fastq_subset_size": "qc.fastq_screen_subset",
        "modulefiles.fastqc": "modulefiles.illumina_qc",
        "modulefiles.fastq_screen": "modulefiles.illumina_qc",
        "runners.cellranger_mkfastq": "runners.cellranger",
        "runners.cellranger_count": "runners.cellranger",
        "runners.cellranger_multi": "runners.cellranger_count",
        "runners.fastqc": "runners.qc",
        "runners.fastq_screen": "runners.qc",
        "runners.picard": "runners.qc",
        "runners.qualimap": "runners.qc",
        "runners.rseqc": "runners.qc",
        "runners.star": "runners.qc",
    }

    # Deprecated parameters
    __DEPRECATED = (
        "bcl2fastq",
        "modulefiles.illumina_qc",
        "fastq_strand_indexes",
        "10xgenomics_transcriptomes",
        "10xgenomics_premrna_references",
        "10xgenomics_atac_genome_references",
        "10xgenomics_multiome_references",
    )

    def __init__(self,settings_file=None,resolve_undefined=True):
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

        Arguments:
          settings_file (str): path to file to read config from
          resolve_undefined (bool): if True then resolve values
            of undefined parameters from fallbacks and defaults
            (default)
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
        self.nullvalue = config.nullvalue
        if self.settings_file:
            config.read(self.settings_file)
        else:
            # Look for sample settings file
            config.read(os.path.join(get_config_dir(),
                                     'auto_process.ini.sample'))
        # General parameters
        self.add_section('general')
        default_runner = config.get('general','default_runner')
        self.general['default_runner'] = config.getrunner('general',
                                                          'default_runner')
        self.general['max_concurrent_jobs'] = config.getint('general',
                                                            'max_concurrent_jobs')
        self.general['max_cores'] = config.getint('general','max_cores')
        self.general['max_batches'] = config.getint('general','max_batches')
        self.general['poll_interval'] = config.getfloat('general',
                                                        'poll_interval')
        # Environment modulefiles
        self.add_section('modulefiles')
        for module_file in self.__ENVIRONMENT_MODULES:
            self.modulefiles[module_file] = config.get('modulefiles',
                                                       module_file)
        # conda
        self.add_section('conda')
        self.conda['enable_conda'] = config.getboolean('conda',
                                                       'enable_conda')
        self.conda['env_dir'] = config.get('conda','env_dir')
        # bcl_conversion
        self.add_section('bcl_conversion')
        # Add settings from legacy bcl2fastq section first
        self.bcl_conversion = self.get_bcl_converter_config('bcl2fastq',
                                                            config)
        # Update with settings from bcl_conversion section
        self.get_bcl_converter_config('bcl_conversion',
                                      config,
                                      self.bcl_conversion)
        # qc
        self.add_section('qc')
        self.qc['nprocessors'] = config.getint('qc','nprocessors')
        self.qc['fastq_screens'] = config.get('qc','fastq_screens')
        self.qc['fastq_screen_subset'] = config.getint('qc',
                                                       'fastq_screen_subset')
        self.qc['fastq_subset_size'] = config.getint(
            'qc',
            'fastq_subset_size')
        self.qc['use_legacy_screen_names'] = config.getboolean(
            'qc',
            'use_legacy_screen_names')
        # Fastq screens
        self.add_section('screens')
        for section in filter(lambda x: x.startswith('screen:'),
                              config.sections()):
            screen = section.split(':')[1]
            self.screens[screen] = AttributeDictionary(conf_file=None)
            self.screens[screen]['conf_file'] = config.get(section,
                                                           'conf_file')
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
            logger.warning("Added STAR index information from "
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
            logger.warning("Added cellranger references from deprecated "
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
            logger.warning("Added cellranger pre-mRNA references from "
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
            logger.warning("Added cellranger-atac references from deprecated "
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
            logger.warning("Added cellranger-arc references from deprecated "
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
            logger.warning("Added sequencer information from "
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
            self.platform[platform] = self.get_bcl_converter_config(section,
                                                                    config)
        # Handle deprecated bcl2fastq settings
        for platform in ('hiseq','miseq','nextseq'):
            if config.has_option('bcl2fastq',platform):
                logger.warning("Deprecated setting in [bcl2fastq]: '%s'"
                               % platform)
            try:
                bcl2fastq = self.platform[platform]['bcl2fastq']
            except KeyError:
                bcl2fastq = config.get('bcl2fastq',platform)
                if bcl2fastq is None or bcl2fastq is self.nullvalue:
                    continue
                logger.warning("Setting 'bcl2fastq' in '[platform:%s]' to '%s'"
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
        self.icell8['batch_size'] = config.getint('icell8','batch_size')
        self.icell8['mammalian_conf_file'] = config.get('icell8',
                                                        'mammalian_conf_file')
        self.icell8['contaminants_conf_file'] = config.get('icell8',
                                                           'contaminants_conf_file')
        self.icell8['nprocessors_contaminant_filter'] = config.getint('icell8','nprocessors_contaminant_filter')
        self.icell8['nprocessors_statistics'] = config.getint('icell8','nprocessors_statistics')
        # 10xgenomics
        self.add_section('10xgenomics')
        tenx = self['10xgenomics']
        tenx['cellranger_jobmode'] = config.get('10xgenomics',
                                                'cellranger_jobmode')
        tenx['cellranger_maxjobs'] = config.getint('10xgenomics',
                                                   'cellranger_maxjobs')
        tenx['cellranger_mempercore'] = config.getint('10xgenomics',
                                                      'cellranger_mempercore')
        tenx['cellranger_jobinterval'] = config.getint('10xgenomics',
                                                       'cellranger_jobinterval')
        tenx['cellranger_localmem'] = config.getint('10xgenomics',
                                                    'cellranger_localmem')
        tenx['cellranger_localcores'] = config.getint('10xgenomics',
                                                      'cellranger_localcores')
        # fastq_stats
        self.add_section('fastq_stats')
        self.fastq_stats['nprocessors'] = config.getint('fastq_stats',
                                                        'nprocessors')
        # Define runners for specific jobs
        self.add_section('runners')
        for name in self.__RUNNERS:
            self.runners[name] = config.getrunner('runners',name)
        # Information for archiving analyses
        # dirn should be a directory in the form [[user@]host:]path]
        self.add_section('archive')
        self.archive['dirn'] = config.get('archive','dirn')
        self.archive['log'] = config.get('archive','log')
        self.archive['group'] = config.get('archive','group')
        self.archive['chmod'] = config.get('archive','chmod')
        # Information for uploading QC reports
        # dirn should be a directory in the form [[user@]host:]path]
        self.add_section('qc_web_server')
        self.qc_web_server['dirn'] = config.get('qc_web_server','dirn')
        self.qc_web_server['url'] = config.get('qc_web_server','url')
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
            logger.debug("No reporting templates defined")
        # Destinations for data transfer
        self.add_section('destination')
        for section in filter(lambda x: x.startswith('destination:'),
                              config.sections()):
            dest = section.split(':')[1]
            self.destination[dest] = self.get_destination_config(
                section,config)
        # Set defaults
        if resolve_undefined:
            self.resolve_undefined_params()

    def get_bcl_converter_config(self,section,config,attr_dict=None):
        """
        Retrieve BCL conversion configuration options from .ini file

        Given the name of a section (e.g. 'bcl_conversion',
        'platform:miseq'), fetch the BCL converter settings and return
        in an AttributeDictionary object.

        The options that can be extracted are:

        - bcl_converter
        - nprocessors
        - no_lane_splitting
        - create_empty_fastqs

        There are also some legacy options:

        - default_version
        - bcl2fastq

        Arguments:
          section (str): name of the section to retrieve the
            settings from
          config (Config): Config object with settings loaded
          attr_dict (AttributeDictionary): optional, existing
            AttributeDictionary which will be added to

        Returns:
          AttributeDictionary: dictionary of option:value pairs.

        """
        if attr_dict:
            values = attr_dict
        else:
            values = AttributeDictionary()
        if section == 'bcl2fastq':
            # Deprecated [bcl2fastq] section
            value = config.get(section,'default_version')
            if value:
                values['bcl_converter'] = "bcl2fastq%s" % value
        else:
            # [bcl_conversion] and [platform:...] sections
            bcl2fastq = config.get(section,'bcl2fastq')
            value = config.get(section,'bcl_converter')
            if value:
                values['bcl_converter'] = value
            elif bcl2fastq:
                values['bcl_converter'] = "bcl2fastq%s" % bcl2fastq
            elif 'bcl_converter' not in values:
                values['bcl_converter'] = self.nullvalue
        # Common settings
        value = config.getint(section,'nprocessors')
        if value or 'nprocessors' not in values:
            values['nprocessors'] = value
        value = config.getboolean(section,'no_lane_splitting')
        if value is not self.nullvalue or 'no_lane_splitting' not in values:
            values['no_lane_splitting'] = value
        value = config.getboolean(section,'create_empty_fastqs')
        if value is not self.nullvalue or 'create_empty_fastqs' not in values:
            values['create_empty_fastqs'] = value
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
        - zip_fastqs (optional, boolean, default 'False')
        - max_zip_size (option, str, default 'None')
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
        values['directory'] = config.get(section,'directory')
        values['subdir'] = config.get(section,'subdir')
        values['zip_fastqs'] = config.getboolean(section,'zip_fastqs')
        values['max_zip_size'] = config.get(section,'max_zip_size')
        values['readme_template'] = config.get(section,'readme_template')
        values['url'] = config.get(section,'url')
        values['include_downloader'] = config.getboolean(section,
                                                         'include_downloader')
        values['include_qc_report'] = config.getboolean(section,
                                                        'include_qc_report')
        values['hard_links'] = config.getboolean(section,'hard_links')
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
        - annotation_bed (str, path to BED file with annotation)
        - annotation_gtf (str, path to GTF file with annotation)
        - cellranger_reference (str)
        - cellranger_premrna_reference (str)
        - cellranger_atac_reference (str)
        - cellranger_arc_reference (str)
        - cellranger_probe_set (str)

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
                'annotation_bed',
                'annotation_gtf',
                'cellranger_premrna_reference',
                'cellranger_atac_reference',
                'cellranger_arc_reference',
                'cellranger_probe_set'):
            if section and config:
                values[param] = config.get(section,param)
            else:
                values[param] = self.nullvalue
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
        values['platform'] = config.get(section,'platform')
        values['model'] = config.get(section,'model')
        if values['platform'] is None or values['platform'] is self.nullvalue:
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

        NB parameters are referenced by the names that
        appear in the config file (rather than the
        internal representation within the Settings
        instance).

        Arguments:
          param (str): an identifier of the form
            SECTION[:SUBSECTION].ATTR which specifies the
            parameter to update
          value (str): the new value of the parameter
        """
        section,attr = param.split('.')
        try:
            section,subsection = section.split(':')
            section = self._section_internal_name(section)
            getattr(self,section)[subsection][attr] = value
            ##print("set: %s:%s.%s -> %r" % (section,
            ##                               subsection,
            ##                               attr,
            ##                               value))
        except ValueError:
            getattr(self,section)[attr] = value
            ##print("set: %s.%s -> %r" % (section,
            ##                            attr,
            ##                            value))

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
            section = self._section_internal_name(section)
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

    def __contains__(self,item):
        """
        Implement __contains__ to enable 'SECTION[:SUBSECTION][.ATTR] in s'
        """
        # Break up item into section, subsection and attr
        try:
            section,attr = item.split('.')
        except ValueError:
            section = item
            attr = None
        try:
            section,subsection = section.split(':')
        except ValueError:
            subsection = None
        # Get the internal name for the section
        section = self._section_internal_name(section)
        # Check for existence of item components
        if section not in self._sections:
            # Section not found
            return False
        s = getattr(self,section)
        if subsection is not None:
            if subsection not in getattr(self,section):
                # Subsection not found
                return False
            s = getattr(self,section)[subsection]
        if attr is not None:
            # Return existence or otherwise of attr
            return attr in s
        else:
            # All checks passed
            return True

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

    def fetch_value(self,param):
        """
        Return the value stored against a parameter

        NB parameters are referenced by the names that
        appear in the config file (rather than the
        internal representation within the Settings
        instance).

        Arguments:
          param (str): parameter name of the form
            SECTION[:SUBSECTION][.ATTR]
        """
        # Break up param into section, subsection and attr
        try:
            section,attr = param.split('.')
        except ValueError:
            section = param
            attr = None
        try:
            section,subsection = section.split(':')
        except ValueError:
            subsection = None
        # Get the internal name for the section
        section = self._section_internal_name(section)
        if subsection is None:
            s = getattr(self,section)
        else:
            s = getattr(self,section)[subsection]
        return s[attr]

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
            SECTION[:SUBSECTION][.ATTR]
        """
        if pattern:
            if '.' not in pattern:
                pattern += '.*'
        for section in self._sections:
            if not self.has_subsections(section):
                name = self._section_config_name(section)
                values = getattr(self,section)
                for attr in values:
                    param = "%s.%s" % (name,attr)
                    if pattern and not fnmatch(param,pattern):
                        continue
                    if exclude_undefined and \
                       self.fetch_value(param) is self.nullvalue:
                        continue
                    else:
                        yield param
            else:
                for subsection in getattr(self,section):
                    name = "%s:%s" % (self._section_config_name(section),
                                      subsection)
                    values = getattr(self,section)[subsection]
                    for attr in values:
                        param = "%s.%s" % (name,attr)
                        if pattern and not fnmatch(param,pattern):
                            continue
                        if exclude_undefined and \
                           self.fetch_value(param) is self.nullvalue:
                            continue
                        else:
                            yield param

    def resolve_undefined_params(self):
        """
        Set non-null values for all parameters which are null

        Resolution has three stages:

        - Fallback parameters: if a parameter is undefined and
          has a defined fallback that value is assigned
        - Default paramaters: if a parameter is undefined after
          fallbacks are exhausted and has a defined default then
          that value is assigned
        - Default runner: any undefined runner parameters are
          assigned the value of the default runner, if set

        Any parameters that are still undefined at this point
        are then assigned the value of 'None'.
        """
        # Fallbacks
        for param in self.__FALLBACKS:
            if param not in self or self.fetch_value(param) is self.nullvalue:
                tried_params = set()
                fallback_param = self.__FALLBACKS[param]
                while fallback_param and fallback_param not in tried_params:
                    value = self.fetch_value(fallback_param)
                    if value is not self.nullvalue:
                        if fallback_param in self.__DEPRECATED:
                            logger.warning("Setting '%s' parameter "
                                           "using value from deprecated "
                                           "'%s' parameter" %
                                           (param,fallback_param))
                        else:
                            print("fallback: setting '%s' to value of '%s' "
                                  "(%r)" % (param,fallback_param,value))
                        self.set(param,value)
                        fallback_param = None
                    elif fallback_param in self.__FALLBACKS:
                        tried_params.add(fallback_param)
                        fallback_param = self.__FALLBACKS[fallback_param]
                    else:
                        fallback_param = None
        # Defaults
        for param in self.__DEFAULTS:
            for p in self.list_params(pattern=param):
                if self.fetch_value(p) is self.nullvalue:
                    print("updating '%s' with default value %r" %
                          (p,self.__DEFAULTS[param]))
                    self.set(p,self.__DEFAULTS[param])
        # Set default runners
        default_runner = self.fetch_value("general.default_runner")
        if default_runner:
            for param in self.list_params(pattern="runners"):
                if self.fetch_value(param) is self.nullvalue:
                    print("Updating '%s' to default runner" % param)
                    self.set(param,default_runner)
        # Set remaining undefined parameters to 'None'
        for param in self.list_params():
            if self.fetch_value(param) is self.nullvalue:
                print("Updating undefined parameter '%s' to None" % param)
                self.set(param,None)
        # Expand variables
        for param in self.__EXPAND_VARS:
            value = self.fetch_value(param)
            if value:
                self.set(param,os.path.expandvars(value))

    def save(self,out_file=None,exclude_undefined=True):
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
            out_file = self.settings_file
        if out_file:
            out_file = os.path.abspath(out_file)
            config = Config()
            for param in self.list_params(exclude_undefined=exclude_undefined):
                name,attr = param.split('.')
                if ':' in name:
                    name = "%s:%s" % (
                        self._section_config_name(name.split(':')[0]),
                        name.split(':')[1])
                else:
                    name = self._section_config_name(name)
                if not config.has_section(name):
                    config.add_section(name)
                value = self.fetch_value(param)
                if exclude_undefined and value is self.nullvalue:
                    continue
                config.set(name,attr,str(value))
            with open(out_file,'wt') as fp:
                config.write(fp)
        else:
            logger.warning("No settings file found, nothing saved")
    
    def report_settings(self,exclude_undefined=False):
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
        if self.settings_file:
            text.append("Settings from %s" % self.settings_file)
        else:
            logger.warning("No settings file found, reporting built-in "
                           "defaults")
        for section in self._sections:
            display_name = self._section_config_name(section)
            if self.has_subsections(section):
                for subsection in getattr(self,section):
                    text.append(
                        show_dictionary('%s:%s' % (display_name,subsection),
                                        getattr(self,section)[subsection],
                                        exclude_value=exclude_value))
            else:
                text.append(show_dictionary(display_name,
                                            getattr(self,section),
                                            exclude_value=exclude_value))
        return '\n'.join(text)

    def _section_config_name(self,internal_name):
        """
        Internal: look up section name used in config file

        Given the internal name for a section (e.g. 'sequencers')
        returns the name used for that section in the config file
        (e.g. 'sequencer').
        """
        try:
            return self.__SECTIONS_INTERNAL_TO_CONFIG_NAMES[internal_name]
        except KeyError:
            return internal_name

    def _section_internal_name(self,config_name):
        """
        Internal: look up section name for 'Settings' class

        Given the name for a section as it appears in the config
        file (e.g. 'sequencer'), returns the name used internally
        for that section within the 'Settings' class (e.g.
        'sequencers').
        """
        for name in self.__SECTIONS_INTERNAL_TO_CONFIG_NAMES:
            if self.__SECTIONS_INTERNAL_TO_CONFIG_NAMES[name] == \
               config_name:
                return name
        return config_name

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
        if data_item:
            refdata[organism] = data_item
    return refdata

def show_dictionary(name,d,exclude_value=None):
    """
    Print the contents of a dictionary

    Arguments:
      name (str): name of dictionary
      d (str): dictionary instance to show
      exclude_value (object): optional, if not 'None'
        then don't include entries which match this
        value
    """
    text = ["[%s]" % name]
    for key in d:
        if exclude_value is not None and d[key] is exclude_value:
            continue
        text.append("\t%s = %s" % (key,(d[key] if d[key] is not None
                                        else '<Not set>')))
    return '\n'.join(text)
