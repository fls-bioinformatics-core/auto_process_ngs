#!/usr/bin/env python
#
#     samplesheet_cmd.py: implement 'samplesheet' command
#     Copyright (C) University of Manchester 2022 Peter Briggs
#
#########################################################################

#######################################################################
# Imports
#######################################################################

import os
import shutil
import json
import logging
from bcftbx.IlluminaData import SampleSheet
from ..samplesheet_utils import check_and_warn
from ..samplesheet_utils import predict_outputs
from ..samplesheet_utils import set_samplesheet_column
from ..utils import edit_file
from ..utils import paginate

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Constants
#######################################################################

class SampleSheetOperation:
    SET_PROJECT = 0
    SET_SAMPLE_NAME = 1
    SET_SAMPLE_ID = 2
    VIEW = 4
    PREDICT = 5
    EDIT = 6

#######################################################################
# Command functions
#######################################################################

def samplesheet(ap,cmd,*args,**kws):
    """
    Various sample sheet manipulations

    Arguments:
      ap (AutoProcessor): autoprocessor pointing to the
        analysis directory to operate on
      cmd (int): sample sheet operation to perform
      args (list): positional arguments specific to the
        command
      kws (mapping): keyword arguments specific to the
        command
    """
    if cmd == SampleSheetOperation.SET_PROJECT:
        # Set the sample project
        set_project(ap,*args,**kws)
    elif cmd == SampleSheetOperation.SET_SAMPLE_NAME:
        # Set the sample names
        set_sample_name(ap,*args,**kws)
    elif cmd == SampleSheetOperation.SET_SAMPLE_ID:
        # Set the sample IDs
        set_sample_id(ap,*args,**kws)
    elif cmd == SampleSheetOperation.VIEW:
        # Show raw sample sheet
        view_samplesheet(ap)
    elif cmd == SampleSheetOperation.PREDICT:
        # Predict the outputs
        predict_samplesheet_outputs(ap)
    elif cmd == SampleSheetOperation.EDIT:
        # Show raw sample sheet
        edit_samplesheet(ap)

def set_project(ap,new_project,lanes=None,where=None):
    """
    Update the project names in the SampleSheet

    Arguments:
      ap (AutoProcessor): autoprocessor pointing to the
        analysis directory to operate on
      new_project (str): new project name
      lanes (list): optional list of lane numbers to
        apply project name update to
      where (tuple): optional tuple '(COLUMN,PATTERN)'
        to only apply update to lines where value in
        COLUMN matches glob-style PATTERN
    """
    s = set_samplesheet_column(
        sample_sheet=ap.params.sample_sheet,
        column='SAMPLE_PROJECT',
        new_value=new_project,
        lanes=lanes,
        where=where)
    s.write(ap.params.sample_sheet)

def set_sample_name(ap,new_name,lanes=None,where=None):
    """
    Update the sample names in the SampleSheet

    Arguments:
      ap (AutoProcessor): autoprocessor pointing to the
        analysis directory to operate on
      new_name (str): new sample name
      lanes (list): optional list of lane numbers to
        apply project name update to
      where (tuple): optional tuple '(COLUMN,PATTERN)'
        to only apply update to lines where value in
        COLUMN matches glob-style PATTERN
    """
    s = set_samplesheet_column(
        sample_sheet=ap.params.sample_sheet,
        column='SAMPLE_NAME',
        new_value=new_name,
        lanes=lanes,
        where=where)
    s.write(ap.params.sample_sheet)

def set_sample_id(ap,new_id,lanes=None,where=None):
    """
    Update the sample IDs in the SampleSheet

    Arguments:
      ap (AutoProcessor): autoprocessor pointing to the
        analysis directory to operate on
      new_id (str): new sample ID
      lanes (list): optional list of lane numbers to
        apply project name update to
      where (tuple): optional tuple '(COLUMN,PATTERN)'
        to only apply update to lines where value in
        COLUMN matches glob-style PATTERN
    """
    s = set_samplesheet_column(
        sample_sheet=ap.params.sample_sheet,
        column='SAMPLE_ID',
        new_value=new_id,
        lanes=lanes,
        where=where)
    s.write(ap.params.sample_sheet)

def view_samplesheet(ap):
    """
    Show the raw SampleSheet content

    Arguments:
      ap (AutoProcessor): autoprocessor pointing to the
        analysis directory to operate on
    """
    with open(ap.params.sample_sheet,'rt') as fp:
        paginate(fp.read())

def predict_samplesheet_outputs(ap):
    """
    Predict the outputs from the SampleSheet

    Arguments:
      ap (AutoProcessor): autoprocessor pointing to the
        analysis directory to operate on
    """
    paginate(predict_outputs(
        sample_sheet_file=ap.params.sample_sheet))

def edit_samplesheet(ap):
    """
    Bring up SampleSheet in an editor

    Arguments:
      ap (AutoProcessor): autoprocessor pointing to the
        analysis directory to operate on
    """
    # Fetch the sample sheet
    sample_sheet_file = ap.params.sample_sheet
    if sample_sheet_file is None:
        logging.error("No sample sheet file to edit")
        return
    edit_file(sample_sheet_file)
    # Check updated sample sheet and issue warnings
    if check_and_warn(sample_sheet_file=sample_sheet_file):
        logger.error("Sample sheet may have problems, see warnings above")
        return 1
    return 0

