#!/usr/bin/env python
#
#     samplesheet_utils.py: utilities for handling samplesheet files
#     Copyright (C) University of Manchester 2016 Peter Briggs
#
########################################################################
#
# samplesheet_utils.py
#
#########################################################################

"""
samplesheet_utils.py

Utilities for handling SampleSheet files:

- predict_outputs: generate expected outputs in human-readable form
- check_and_warn: check sample sheet for problems and issue warnings
- close_project_names: check if sample sheet projects look similar
- samples_with_multiple_barcodes: check for samples with multiple barcodes
- samples_in_multiple_projects: check for samples assigned to multiple projects
- has_invalid_lines: check for invalid sample sheet lines
- has_invalid_characters: check if a file contains invalid characters
- get_close_names: return closely matching names from a list

"""

#######################################################################
# Imports
#######################################################################

import logging
import difflib
from bcftbx.IlluminaData import SampleSheet
from bcftbx.IlluminaData import SampleSheetPredictor

#######################################################################
# Functions
#######################################################################

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
    # Set up predictor instance
    if sample_sheet is None:
        sample_sheet = SampleSheet(sample_sheet_file)
    predictor = SampleSheetPredictor(sample_sheet=sample_sheet)
    # Do checks
    close_names = get_close_names(predictor.project_names)
    # Generate prediction report
    prediction = []
    title = "Predicted projects:"
    prediction.append("%s\n%s" % (title,('='*len(title))))
    for project_name in predictor.project_names:
        if project_name in close_names:
            warning_flag = " *"
        else:
            warning_flag = ""
        prediction.append("- %s%s" % (project_name,
                                      warning_flag))
    for project_name in predictor.project_names:
        project = predictor.get_project(project_name)
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
    # Acquire sample sheet instance
    if sample_sheet is None:
        sample_sheet = SampleSheet(sample_sheet_file)
    # Do checks
    warnings = False
    if close_project_names(sample_sheet=sample_sheet):
        logging.warning("Some projects have similar names: check for typos")
        warnings = True
    if samples_with_multiple_barcodes(sample_sheet=sample_sheet):
        logging.warning("Some samples have more than one barcode assigned")
        warnings = True
    if samples_in_multiple_projects(sample_sheet=sample_sheet):
        logging.warning("Some samples appear in more than one project")
        warnings = True
    if has_invalid_characters(text=sample_sheet.show()):
        logging.warning("Sample sheet file contains invalid characters "
                        "(non-printing ASCII or non-ASCII)")
        warnings = True
    if has_invalid_lines(sample_sheet=sample_sheet):
        logging.warning("Sample sheet has one or more invalid lines")
        warnings = True
    return warnings

def close_project_names(sample_sheet=None,sample_sheet_file=None):
    """
    Return list of closely-matching project names in samplesheet

    Arguments:
      sample_sheet (SampleSheet): if supplied then must be a
        populated ``SampleSheet`` instance (if ``None`` then
        data will be loaded from file specified by
        ``sample_sheet_file``)
      sample_sheet_file (str): if ``sample_sheet`` is ``None``
        then read data from the file specified by this argument

    Returns:
      Dictionary: keys are project names which have at least one
        close match; the values for each key are lists with the
        project names which are close matches.

    """
    # Set up predictor instance
    if sample_sheet is None:
        sample_sheet = SampleSheet(sample_sheet_file)
    predictor = SampleSheetPredictor(sample_sheet=sample_sheet)
    # Return close names
    return get_close_names(predictor.project_names)

def samples_with_multiple_barcodes(sample_sheet=None,sample_sheet_file=None):
    """
    Return list of samples which have multiple associated barcodes

    Arguments:
      sample_sheet (SampleSheet): if supplied then must be a
        populated ``SampleSheet`` instance (if ``None`` then
        data will be loaded from file specified by
        ``sample_sheet_file``)
      sample_sheet_file (str): if ``sample_sheet`` is ``None``
        then read data from the file specified by this argument

    Returns:
      Dictionary: keys are sample IDs which have more than one
        associated barcode; the values for each key are lists of
        the associated barcodes.

    """
    # Set up predictor instance
    if sample_sheet is None:
        sample_sheet = SampleSheet(sample_sheet_file)
    predictor = SampleSheetPredictor(sample_sheet=sample_sheet)
    logging.debug(predict_outputs(sample_sheet=sample_sheet))
    # Look for samples with multiple barcodes
    multiple_barcodes = {}
    for project in [predictor.get_project(name)
                    for name in predictor.project_names]:
        for sample in [project.get_sample(idx)
                       for idx in project.sample_ids]:
            if len(sample.barcode_seqs) > 1:
                multiple_barcodes[sample.sample_id] = \
                            [s for s in sample.barcode_seqs]
    return multiple_barcodes

def samples_in_multiple_projects(sample_sheet=None,
                                 sample_sheet_file=None):
    """
    Return list of samples which are in multiple projects

    Arguments:
      sample_sheet (SampleSheet): if supplied then must be a
        populated ``SampleSheet`` instance (if ``None`` then
        data will be loaded from file specified by
        ``sample_sheet_file``)
      sample_sheet_file (str): if ``sample_sheet`` is ``None``
        then read data from the file specified by this argument

    Returns:
      Dictionary: dictionary with sample IDs which appear in
        multiple projects as keys; the associated values are
        lists with the project names.

    """
    # Set up predictor instance
    if sample_sheet is None:
        sample_sheet = SampleSheet(sample_sheet_file)
    predictor = SampleSheetPredictor(sample_sheet=sample_sheet)
    # Look for samples with multiple projects
    samples = {}
    for project in [predictor.get_project(name)
                    for name in predictor.project_names]:
        for sample in [project.get_sample(name)
                       for name in project.sample_ids]:
            if sample.sample_id not in samples:
                samples[sample.sample_id] = []
            samples[sample.sample_id].append(project.name)
    multiple_projects = {}
    for sample in samples:
        if len(samples[sample]) > 1:
            multiple_projects[sample] = samples[sample]
    return multiple_projects

def has_invalid_lines(sample_sheet=None,sample_sheet_file=None):
    """
    Return list of samplesheet lines which are invalid

    Arguments:
      sample_sheet (SampleSheet): if supplied then must be a
        populated ``SampleSheet`` instance (if ``None`` then
        data will be loaded from file specified by
        ``sample_sheet_file``)
      sample_sheet_file (str): if ``sample_sheet`` is ``None``
        then read data from the file specified by this argument

    Returns:
      List: list of lines which are invalid (i.e. missing
        required data) in the sample sheet.

    """
    # Look for invalid data lines
    if sample_sheet is None:
        sample_sheet = SampleSheet(sample_sheet_file)
    invalid_lines = []
    for line in sample_sheet.data:
        if sample_sheet.has_lanes and line['Lane'] == '':
            invalid_lines.append(line)
    return invalid_lines

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
