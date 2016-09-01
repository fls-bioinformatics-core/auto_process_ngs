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
- get_close_names: return closely matching names from a list
- has_invalid_characters: check if a file contains invalid characters

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

def warn_close_names(sample_sheet=None,sample_sheet_file=None):
    """
    Issue a warning for each set of closely matching project names

    Arguments:
      sample_sheet (SampleSheet): if supplied then must be a
        populated ``SampleSheet`` instance (if ``None`` then
        data will be loaded from file specified by
        ``sample_sheet_file``)
      sample_sheet_file (str): if ``sample_sheet`` is ``None``
        then read data from the file specified by this argument

    Returns:
      Boolean: True if there is at least one pair of closely
        matching project names, False if not.

    """
    # Set up predictor instance
    if sample_sheet is None:
        sample_sheet = SampleSheet(sample_sheet_file)
    predictor = SampleSheetPredictor(sample_sheet=sample_sheet)
    # Look for close names
    close_names = get_close_names(predictor.project_names)
    for name in close_names:
        logging.warning("Project '%s': similar to %d others (%s)" %
                        (name,
                         len(close_names[name]),
                         ','.join(close_names[name])))
    # Return status indicates whether or not there were
    # closely matching names
    if close_names:
        return True
    else:
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

def has_invalid_characters(filen):
    """
    Check if text file contains any 'invalid' characters

    In this context a character is 'invalid' if:
    - it is non-ASCII (decimal code > 127), or
    - it is a non-printing ASCII character (code < 32)

    Returns:
      Boolean: True if file contains at least one invalid
        character, False if all characters are valid.

    """
    with open(filen,'r') as fp:
        for line in fp:
            for c in line.replace('\n','').replace('\t',''):
                if ord(c) > 127 or ord(c) < 32:
                    return True
    return False
