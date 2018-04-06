#!/usr/bin/env python
#
#     samplesheet_utils.py: utilities for handling samplesheet files
#     Copyright (C) University of Manchester 2016-2018 Peter Briggs
#
########################################################################
#
# samplesheet_utils.py
#
#########################################################################

"""
samplesheet_utils.py

Utilities for handling SampleSheet files:

- SampleSheetLinter: core class which provides methods for checking
  sample sheet contents for potential problems
- predict_outputs: generate expected outputs in human-readable form
- check_and_warn: check sample sheet for problems and issue warnings

Helper functions:

- has_invalid_characters: check if a file contains invalid characters
- get_close_names: return closely matching names from a list

"""

#######################################################################
# Imports
#######################################################################

import logging
import difflib
import re
from bcftbx.IlluminaData import SampleSheet
from bcftbx.IlluminaData import SampleSheetPredictor

#######################################################################
# Classes
#######################################################################

class SampleSheetLinter(SampleSheetPredictor):
    """
    Class for checking sample sheets for problems

    Provides the following methods for checking different aspects
    of a sample sheet:

    - close_project_names: check if sample sheet projects look similar
    - samples_with_multiple_barcodes: check for samples with multiple
      barcodes
    - samples_in_multiple_projects: check for samples assigned to
      multiple projects
    - has_invalid_lines: check for invalid sample sheet lines
    - has_invalid_characters: check if sample sheet contains invalid
      characters

    Example usage:

    Initialise linter:
    >>> linter = SampleSheetLinter(sample_sheet_file="SampleSheet.txt")

    Get closely-matching names:
    >>> linter.close_project_names()
    ...

    """
    def __init__(self,sample_sheet=None,sample_sheet_file=None,fp=None):
        """
        Create a new SampleSheetLinter instance

        Arguments:
          sample_sheet (SampleSheet): a SampleSheet instance to use
            for prediction (if None then must provide a file via
            the `sample_sheet_file` argument; if both are provided
            then `sample_sheet` takes precedence)
          sample_sheet_file (str): path to a sample sheet file, if
            `sample_sheet` argument is None
          fp (File): File-like object opened for reading; if this
            is not None then the SampleSheet object will be populated
            from this in preference to `sample_sheet`

        """
        # Initialise
        self._fp = fp
        self._sample_sheet_file = sample_sheet_file
        self._sample_sheet = sample_sheet
        if self._fp is None:
            if self._sample_sheet is None:
                self._sample_sheet = SampleSheet(sample_sheet_file)
        else:
            self._sample_sheet = SampleSheet(fp=self._fp)
        SampleSheetPredictor.__init__(self,
                                      sample_sheet=self._sample_sheet)


    def walk(self):
        """
        Traverse the list of projects and samples

        Generator that yields tuples consisting of
        (SampleSheetProject,SampleSheetSample) pairs
        
        Yields:
          Tuple: SampleSheetProject, SampleSheetSample pair

        """
        for project in [self.get_project(name)
                        for name in self.project_names]:
            for sample in [project.get_sample(idx)
                           for idx in project.sample_ids]:
                yield (project,sample)
        
    def close_project_names(self):
        """
        Return list of closely-matching project names

        Returns:
          Dictionary: keys are project names which have at least one
            close match; the values for each key are lists with the
            project names which are close matches.

        """
        return get_close_names(self.project_names)

    def samples_with_multiple_barcodes(self):
        """
        Return list of samples which have multiple associated barcodes

        Returns:
          Dictionary: keys are sample IDs which have more than one
          associated barcode; the values for each key are lists of
          the associated barcodes.

        """
        # Look for samples with multiple barcodes
        multiple_barcodes = {}
        for project,sample in self.walk():
            if len(sample.barcode_seqs) > 1:
                multiple_barcodes[sample.sample_id] = \
                    [s for s in sample.barcode_seqs]
        return multiple_barcodes

    def samples_in_multiple_projects(self):
        """
        Return list of samples which are in multiple projects

        Returns:
          Dictionary: dictionary with sample IDs which appear in
            multiple projects as keys; the associated values are
            lists with the project names.

        """
        # Look for samples with multiple projects
        samples = {}
        for project,sample in self.walk():
            if sample.sample_id not in samples:
                samples[sample.sample_id] = []
            samples[sample.sample_id].append(project.name)
        multiple_projects = {}
        for sample in samples:
            if len(samples[sample]) > 1:
                multiple_projects[sample] = samples[sample]
        return multiple_projects

    def has_invalid_lines(self):
        """
        Return list of samplesheet lines which are invalid

        Returns:
          List: list of lines which are invalid (i.e. missing
            required data) in the sample sheet.

        """
        # Convience variables
        sample_id = self._sample_sheet.sample_id_column
        sample_name = self._sample_sheet.sample_name_column
        sample_project = self._sample_sheet.sample_project_column
        # Look at first line to see which items have been provided
        line = self._sample_sheet.data[0]
        has_sample_id = line[sample_id] != ''
        has_sample_name = (sample_name is not None) and \
                          (line[sample_name] != '')
        has_project = line[sample_project] != ''
        # Look for invalid data lines
        invalid_lines = []
        for line in self._sample_sheet.data:
            if self._sample_sheet.has_lanes and line['Lane'] == '':
                invalid_lines.append(line)
            elif has_sample_id and line[sample_id] == '':
                invalid_lines.append(line)
            elif has_sample_name and line[sample_name] == '':
                invalid_lines.append(line)
            elif has_project and line[sample_project] == '':
                invalid_lines.append(line)
        return invalid_lines

    def has_invalid_barcodes(self):
        """
        Return list of lines with invalid barcodes

        Returns:
          List: list of lines which contain invalid barcode
            sequences in the sample sheet.
        """
        invalid_lines = list()
        indices = list()
        for indx in ('index','index2'):
            if indx in self._sample_sheet.data.header():
                indices.append(indx)
        if indices:
            for line in self._sample_sheet.data:
                for indx in indices:
                    if not barcode_is_valid(line[indx]):
                        invalid_lines.append(line)
                        continue
        return invalid_lines

    def has_invalid_characters(self):
        """
        Check if text file contains any 'invalid' characters

        In this context a character is 'invalid' if:
        - it is non-ASCII (decimal code > 127), or
        - it is a non-printing ASCII character (code < 32)

        Returns:
          Boolean: True if file contains at least one invalid
            character, False if all characters are valid.

        """
        return has_invalid_characters(text=self._sample_sheet.show())

#######################################################################
# Functions
#######################################################################

def barcode_is_valid(s):
    """
    Check if a sample sheet barcode sequence is valid

    Valid barcodes must consist of only the letters A,T,G or C
    in any order, and always uppercase.

    10xGenomics sample set IDs of the form e.g. 'SI-P03-C9' or
    'SI-GA-B3' are also considered to be valid.

    Arguments:
      s (str): barcode sequence to validate

    Returns:
      Boolean: True if barcode is valid, False if not.

    """
    return bool(re.match(r'^[ATGC]*$',s) or
                re.match(r'^SI\-[A-Z0-9]+\-[A-Z0-9]+$',s))

def predict_outputs(sample_sheet=None,sample_sheet_file=None):
    """
    Generate expected sample sheet output in human-readable form

    Arguments:
      sample_sheet (SampleSheet): if supplied then must be a
        populated ``SampleSheet`` instance (if ``None`` then
        data will be loaded from file specified by
        ``sample_sheet_file``)
      sample_sheet_file (str): if ``sample_sheet`` is ``None``
        then read data from the file specified by this argument

    Returns:
      String: text describing the expected projects, sample names,
        indices, barcodes and lanes.

    """
    # Set up linter
    linter = SampleSheetLinter(sample_sheet=sample_sheet,
                               sample_sheet_file=sample_sheet_file)
    # Do checks
    close_names = linter.close_project_names()
    # Generate prediction report
    prediction = []
    title = "Predicted projects:"
    prediction.append("%s\n%s" % (title,('='*len(title))))
    for project_name in linter.project_names:
        if project_name in close_names:
            warning_flag = " *"
        else:
            warning_flag = ""
        prediction.append("- %s%s" % (project_name,
                                      warning_flag))
    for project_name in linter.project_names:
        project = linter.get_project(project_name)
        title = "%s (%d sample%s)" % (project_name,
                                      len(project.sample_ids),
                                      ('s' if len(project.sample_ids) != 1
                                       else ''))
        prediction.append("\n%s\n%s" % (title,('-'*len(title))))
        for sample_id in project.sample_ids:
            sample = project.get_sample(sample_id)
            for barcode in sample.barcode_seqs:
                lanes = sample.lanes(barcode)
                if lanes:
                    lanes = "L%s" % (','.join([str(l)
                                               for l in lanes]))
                else:
                    lanes = "L*"
                line = [sample_id,
                        "S%d" % sample.s_index,
                        barcode,
                        lanes]
                prediction.append("%s" % '\t'.join([str(i) for i in line]))
    return '\n'.join(prediction)

def check_and_warn(sample_sheet=None,sample_sheet_file=None):
    """
    Check for sample sheet problems and issue warnings

    The following checks are performed:

    - closely matching project names
    - samples with more than one barcode assigned
    - samples associated with more than one project
    - invalid lines
    - invalid characters

    Arguments:
      sample_sheet (SampleSheet): if supplied then must be a
        populated ``SampleSheet`` instance (if ``None`` then
        data will be loaded from file specified by
        ``sample_sheet_file``)
      sample_sheet_file (str): if ``sample_sheet`` is ``None``
        then read data from the file specified by this argument

    Returns:
      Boolean: True if problems were identified, False otherwise.

    """
    # Acquire sample sheet linter instance
    linter = SampleSheetLinter(sample_sheet=sample_sheet,
                               sample_sheet_file=sample_sheet_file)
    # Do checks
    warnings = False
    if linter.close_project_names():
        logging.warning("Some projects have similar names: check for typos")
        warnings = True
    if linter.samples_with_multiple_barcodes():
        logging.warning("Some samples have more than one barcode assigned")
        warnings = True
    if linter.samples_in_multiple_projects():
        logging.warning("Some samples appear in more than one project")
        warnings = True
    if linter.has_invalid_characters():
        logging.warning("Sample sheet file contains invalid characters "
                        "(non-printing ASCII or non-ASCII)")
        warnings = True
    if linter.has_invalid_lines():
        logging.warning("Sample sheet has one or more invalid lines")
        warnings = True
    return warnings

def has_invalid_characters(filen=None,text=None):
    """
    Check if text file contains any 'invalid' characters

    In this context a character is 'invalid' if:
    - it is non-ASCII (decimal code > 127), or
    - it is a non-printing ASCII character (code < 32)

    Returns:
      Boolean: True if file contains at least one invalid
        character, False if all characters are valid.

    """
    if filen is not None:
        with open(filen,'r') as fp:
            for line in fp:
                for c in set(line.replace('\n','').replace('\t','')):
                    if ord(c) > 127 or ord(c) < 32:
                        return True
    else:
        for c in set(text.replace('\n','').replace('\t','')):
            if ord(c) > 127 or ord(c) < 32:
                return True
    return False

def get_close_names(names):
    """
    Given a list of names, find pairs which are similar

    Returns:
      Dictionary: keys are names which have at least one close
        match; the values for each key are lists with the
        close matches.

    """
    close_names = {}
    for i,name in enumerate(names[:-1]):
        matches = difflib.get_close_matches(name,names[i+1:])
        if len(matches):
            #print "*** %s: matches %s ***" % (name,matches)
            close_names[name] = matches
            for name1 in matches:
                if name1 not in close_names:
                    close_names[name1] = []
                close_names[name1].append(name)
    return close_names
