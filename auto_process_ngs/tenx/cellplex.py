#!/usr/bin/env python
#
#     tenx/cellplex.py: utilities for handling 10xGenomics Cellplex data
#     Copyright (C) University of Manchester 2023-2025 Peter Briggs
#

"""
Utilities for working with 10x Genomics single cell multiplexing
(CellPlex) pipelines:

- CellrangerMultiConfigCsv
"""

#######################################################################
# Imports
#######################################################################

import os
import re
from bcftbx.utils import pretty_print_names

# Initialise logging
import logging
logger = logging.getLogger(__name__)

#######################################################################
# Data
#######################################################################

# Internal representation of known feature types that can
# appear in library definitions
KNOWN_FEATURE_TYPES = (
    "gene_expression",
    "multiplexing_capture",
    "antibody_capture",
    "vdj_b",
    "vdj_t",
)

#######################################################################
# Classes
#######################################################################

class CellrangerMultiConfigCsv:
    """
    Class to handle cellranger multi 'config.csv' files

    See https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/multi#cellranger-multi

    Provides the following properties:

    - sample_names: list of multiplexed sample names
    - sections: list of the sections in the config
    - reference_data_path: path to the reference dataset
    - probe_set_path: path to the probe set
    - feature_reference_path: path to the feature reference
    - vdj_reference_path: path to the V(D)J-compatible reference
    - gex_libraries: list of Fastq IDs associated
      with GEX data
    - physical_sample: physical sample name extracted from
      the config file name if present (otherwise None)
    - is_valid: indicates whether the file appears to be valid

    Provides the following methods:

    - sample: returns information on a specific multiplexed
      sample
    - gex_library: returns information on a specific GEX
      library
    - fastq_dirs: returns mapping of library names to the
      associated Fastq directory paths
    - pretty_print_samples: returns a string with a 'nice'
      description of the multiplexed sample names
    - get_errors: returns a list of error messages (if any)
      indicating problems with the config.csv file

    By default data from the config.csv file is read in
    'strict' mode; any errors detected in formatting will
    cause an exception to be raised. If the file is read
    with 'strict' turned off then the 'is_valid' property
    can be used to check if the file is corrected formatted,
    and any errors can be accessed via the 'get_errors'
    method.
    """
    def __init__(self, filen, strict=True):
        """
        Create new CellrangerMultiConfigCsv instance

        Arguments:
          filen (str): path to cellranger multi config.csv
            file
          strict (bool): if True (the default) then raise
            an exception if problems are detected with
            the config file; otherwise skip the offending
            lines
        """
        self._filen = os.path.abspath(filen)
        self._sections = []
        self._samples = {}
        self._physical_sample = None
        self._reference_data_path = None
        self._probe_set_path = None
        self._feature_reference_path = None
        self._vdj_reference_path = None
        self._libraries = {}
        self._fastq_dirs = {}
        self._errors = []
        self._read_config_csv(strict=strict)

    def _read_config_csv(self, strict=True):
        """
        Internal: read in data from a multiplex 'config.csv' file
        """
        logger.debug("Reading data from '%s'" % self._filen)
        psample = ".".join(os.path.basename(self._filen).split(".")[1:-1])
        if psample:
            self._physical_sample = psample
        sections = set()
        cmo_list = set()
        feature_type_list = set()
        with open(self._filen,'rt') as config_csv:
            current_section = None
            for lineno, line in enumerate(config_csv, start=1):
                if current_section:
                    sections.add(current_section)
                line = line.rstrip('\n')
                if line == "[samples]":
                    current_section = "samples"
                    continue
                elif line == "[gene-expression]":
                    current_section = "gene-expression"
                    continue
                elif line == "[feature]":
                    current_section = "feature"
                    continue
                elif line == "[vdj]":
                    current_section = "vdj"
                    continue
                elif line == "[libraries]":
                    current_section = "libraries"
                    continue
                elif not line:
                    # Blank line ends section
                    current_section = None
                    continue
                if current_section == "samples":
                    if line.startswith('sample_id,cmo_ids') or \
                       line.startswith('sample_id,probe_barcode_ids'):
                        # Header line, skip
                        continue
                    else:
                        # Extract sample name
                        values = [x.strip() for x in line.split(',')]
                        if len(values) < 2:
                            self._error("samples",
                                        f"L{lineno}: bad line: {line}")
                            continue
                        sample = values[0]
                        cmos = values[1]
                        for cmo in cmos.split("|"):
                            if not re.fullmatch("[A-Za-z0-9_]+", cmo):
                                self._error("samples",
                                            f"L{lineno}: bad CMO or probe ID "
                                            f"'{cmo}': {line}")
                                continue
                            elif cmo in cmo_list:
                                self._error("samples",
                                            f"L{lineno}: CMO or probe ID  "
                                            f"'{cmo}'(provided for "
                                            f"'{sample}') already in use")
                                continue
                            cmo_list.add(cmo)
                        if len(values) > 2:
                            desc = values[2]
                        else:
                            desc = ""
                        logger.debug("Found sample '%s'" % sample)
                        if sample in self._samples:
                            self._error("samples",
                                        f"L{lineno}: duplicate entry for "
                                        f"sample name '{sample}': {line}")
                            continue
                        self._samples[sample] = { 'cmo': cmos,
                                                  'description': desc }
                elif current_section == "gene-expression":
                    if line.startswith('reference,'):
                        # Extract reference dataset
                        self._reference_data_path = ','.join(
                            line.split(',')[1:]).strip()
                    elif line.startswith('probe-set,'):
                        # Extract probe set
                        self._probe_set_path = ','.join(
                            line.split(',')[1:]).strip()
                elif current_section == "feature":
                    if line.startswith('reference,'):
                        # Extract feature reference file
                        self._feature_reference_path = ','.join(
                            line.split(',')[1:]).strip()
                elif current_section == "vdj":
                    if line.startswith('reference,'):
                        # Extract feature reference file
                        self._vdj_reference_path = ','.join(
                            line.split(',')[1:]).strip()
                elif current_section == "libraries":
                    if line.startswith('fastq_id,fastqs,'):
                        # Header line, skip
                        continue
                    else:
                        # Extract data
                        values = [x.strip() for x in line.split(',')]
                        if len(values) == 6:
                            name,fastqs,lanes,library_id,feature_type,\
                                subsample_rate = values
                        elif len(values) == 3:
                            name,fastqs,feature_type = values
                            lanes = "any"
                            library_id = None
                            subsample_rate = ""
                        else:
                            self._error("libraries",
                                        f"L{lineno}: bad line: {line}")
                            continue
                        # Check feature type
                        feature_name = self._feature_name(feature_type)
                        if feature_name not in KNOWN_FEATURE_TYPES:
                            self._error("libraries",
                                        f"L{lineno}: '{feature_type}': "
                                        "unrecognised feature type for "
                                        f"sample '{name}'")
                        elif feature_type.lower() in feature_type_list:
                            self._error("libraries",
                                        f"L{lineno}: '{feature_type}': "
                                        "duplicated feature type for "
                                        f"sample '{name}'")
                        # Store Fastq dir
                        self._fastq_dirs[name] = fastqs
                        # Store library
                        self._libraries[name] = {
                            'fastqs': fastqs,
                            'lanes': lanes,
                            'library_id': library_id,
                            'feature_type': feature_type,
                            'subsample_rate': subsample_rate
                        }
                        feature_type_list.add(feature_type.lower())
        self._sections = sorted(list(sections))
        if not self.is_valid:
            if strict:
                # Report errors and raise exception
                for err in self.get_errors():
                    logger.critical(f"[{err[0]}]: {err[1]}")
                raise Exception(f"Errors encountered reading in {self._filen}")
            else:
                logger.warning(f"Errors encountered reading in {self._filen}")

    @property
    def is_valid(self):
        """
        Indicate whether config.csv file is valid

        Returns:
          Boolean: True if no errors were encountered reading in
            the file, False if not.
        """
        return not bool(self._errors)

    def get_errors(self):
        """
        Return errors detected on reading in the config.csv file

        Returns:
          List: list of error messages that were encountered;
            will be empty if there were no errors.
        """
        return [err for err in self._errors]

    @property
    def sample_names(self):
        """
        Return the multiplexed sample names from config.csv

        Samples are listed in the '[samples]' section.
        """
        return sorted(list(self._samples.keys()))

    @property
    def physical_sample(self):
        """
        Return the physical sample from config.csv name

        Physical sample name is extracted from config file
        names of the form:

        10x_multi_config[.SAMPLE].csv

        If no physical sample name is present in the name
        then returns 'None'.
        """
        return self._physical_sample

    @property
    def sections(self):
        """
        Return the list of sections in the config.csv file
        """
        return self._sections

    @property
    def reference_data_path(self):
        """
        Return the path to the reference dataset from config.csv
        """
        return self._reference_data_path

    @property
    def probe_set_path(self):
        """
        Return the path to the probe set file from config.csv
        """
        return self._probe_set_path

    @property
    def feature_reference_path(self):
        """
        Return the path to the feature reference file from config.csv
        """
        return self._feature_reference_path

    @property
    def vdj_reference_path(self):
        """
        Return the path to the V(D)J reference file from config.csv
        """
        return self._vdj_reference_path

    @property
    def gex_libraries(self):
        """
        Return the library names associated with GEX data from config.csv

        Libraries are listed in the '[libraries]' section
        """
        return self.libraries("gene_expression")

    @property
    def fastq_dirs(self):
        """
        Return mapping of library names to Fastq directories
        """
        return { k: self._fastq_dirs[k] for k in self._fastq_dirs }

    def sample(self,sample_name):
        """
        Return dictionary of values associated with multiplexed sample

        Keys include 'cmo' (list of CMO ids) and 'description'
        (description text) associated with the sample in the
        '[samples]' section of the config.csv file.

        Arguments:
          sample_name (str): name of the sample of interest
        """
        return self._samples[sample_name]

    @property
    def feature_types(self):
        """
        Return list of feature types defined in config file

        Feature type names are returned converted to lower case.
        """
        return sorted(list(set([self._libraries[n]['feature_type'].lower()
                                for n in self._libraries])))

    def libraries(self,feature_type):
        """
        Return library names associated with specified feature type
        """
        return sorted(list(self._libraries_for_feature_type(feature_type)))

    def _feature_name(self,feature_type):
        """
        Convert feature type to internal feature name

        Converts the feature type as it appears in the config
        file (e.g. 'Gene Expression', 'VDJ-T') to an internal
        version, by converting to lower case and replacing
        spaces and hyphens with underscores (e.g.
        'gene_expression', 'vdj_t')
        """
        return str(feature_type).lower().replace(' ','_').replace('-','_')

    def _libraries_for_feature_type(self,feature_type):
        """
        Return subset of libraries matching specified feature type

        Subset is returned as a dictionary.

        Arguments:
          feature_type (str): feature to return list of (e.g.
            'Gene expression', 'VDJ-T')

        Returns:
          Dictionary: keys are library names with the matching
            feature type.

        Raises:
          KeyError: if the feature type is not found.
        """
        feature_name = self._feature_name(feature_type)
        try:
            return {
                n: self._libraries[n]
                for n in self._libraries
                if self._feature_name(self._libraries[n]['feature_type'])
                == feature_name
            }
        except KeyError:
            raise KeyError("'%s': feature type not found" % feature_type)

    def _error(self, section, err):
        """
        Register an error in the config file

        Arguments:
          section (str): name of the section that the error
            occurred in
          err (str): the details of the error
        """
        self._errors.append((section, err))

    def library(self,feature_type,name):
        """
        Return dictionary of values associated with library

        Keys include:

        - 'fastqs' (path to Fastqs)
        - 'lanes' (associated lanes)
        - 'library_id' (physical library ID)
        - 'feature_type' (e.g. 'Gene Expression')
        - 'subsample_rate' (the associated subsampling rate)

        Arguments:
          feature_type (str): feature type of the library of
            interest (e.g. 'Gene Expression')
          name (str): name of the library of interest
        """
        return self._libraries_for_feature_type(feature_type)[name]

    def gex_library(self,name):
        """
        Return dictionary of values associated with GEX library

        Arguments:
          name (str): name of the sample of interest
        """
        return self.library("gene expression",name)

    def pretty_print_samples(self):
        """
        Return string describing the multiplexed sample names

        Wraps a call to 'pretty_print_names' function.

        Returns:
          String: pretty description of multiplexed sample names.
        """
        return pretty_print_names(self.sample_names)
