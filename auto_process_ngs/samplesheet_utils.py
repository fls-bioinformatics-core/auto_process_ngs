#!/usr/bin/env python
#
#     samplesheet_utils.py: utility functions for handling samplesheet files
#     Copyright (C) University of Manchester 2016 Peter Briggs
#
########################################################################
#
# samplesheet_utils.py
#
#########################################################################

"""
samplesheet_utils.py

Utility functions for handling SampleSheet files:

- predict_outputs: generate expected outputs in human-readable form

"""

#######################################################################
# Imports
#######################################################################

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
    prediction = []
    if sample_sheet is None:
        sample_sheet = SampleSheet(sample_sheet_file)
    predictor = SampleSheetPredictor(sample_sheet=sample_sheet)
    title = "Predicted projects:"
    prediction.append("%s\n%s" % (title,('='*len(title))))
    for project_name in predictor.project_names:
        prediction.append("- %s" % project_name)
    for project_name in predictor.project_names:
        project = predictor.get_project(project_name)
        title = "%s (%d samples)" % (project_name,
                                     len(project.sample_ids))
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
