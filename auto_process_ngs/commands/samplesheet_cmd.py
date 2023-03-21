#!/usr/bin/env python
#
#     samplesheet_cmd.py: implement 'samplesheet' command
#     Copyright (C) University of Manchester 2023 Peter Briggs
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
from ..bcl2fastq.utils import make_custom_sample_sheet
from ..samplesheet_utils import check_and_warn
from ..samplesheet_utils import predict_outputs
from ..samplesheet_utils import set_samplesheet_column
from ..utils import edit_file
from ..utils import fetch_file
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
    IMPORT = 7

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
    elif cmd == SampleSheetOperation.IMPORT:
        # Replace with new sample sheet
        import_samplesheet(ap,*args,**kws)

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
        content = "Current samplesheet '%s'\n\n%s" % (
            ap.params.sample_sheet,
            fp.read())
        paginate(content)

def predict_samplesheet_outputs(ap):
    """
    Predict the outputs from the SampleSheet

    Arguments:
      ap (AutoProcessor): autoprocessor pointing to the
        analysis directory to operate on
    """
    content = "Current samplesheet '%s'\n\n%s" % (
        ap.params.sample_sheet,
        predict_outputs(sample_sheet_file=ap.params.sample_sheet))
    paginate(content)

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
        logger.error("Sample sheet has problems, see warnings above")

def import_samplesheet(ap,new_sample_sheet):
    """
    Update the SampleSheet contents from a file

    Arguments:
      ap (AutoProcessor): autoprocessor pointing to the
        analysis directory to operate on
      sample_sheet (str): path or URL pointing to the
        SampleSheet file to import the contents from
    """
    # Fetch the new sample sheet
    tmp_sample_sheet = os.path.join(ap.tmp_dir,
                                    os.path.basename(new_sample_sheet))
    try:
        fetch_file(new_sample_sheet,tmp_sample_sheet)
    except Exception as ex:
        raise Exception("Error importing sample sheet data "
                        "from '%s': %s" % (new_sample_sheet,ex))
    # Keep a copy of the imported sample sheet
    imported_sample_sheet = os.path.join(ap.analysis_dir,
                                         'SampleSheet.imported.csv')
    print("Copying imported sample sheet to %s" %
          imported_sample_sheet)
    shutil.copyfile(tmp_sample_sheet,imported_sample_sheet)
    os.chmod(imported_sample_sheet,0o664)
    # Process acquired sample sheet
    custom_sample_sheet = os.path.join(ap.analysis_dir,
                                       'custom_%s' %
                                       os.path.basename(
                                           imported_sample_sheet))
    print("Generating custom version and writing to '%s'" %
          custom_sample_sheet)
    make_custom_sample_sheet(tmp_sample_sheet,custom_sample_sheet)
    # Update metadata to make this the default sample sheet
    print("Updating the default sample sheet")
    ap.params['sample_sheet'] = custom_sample_sheet
    # Generate and print predicted outputs
    print(predict_outputs(sample_sheet=SampleSheet(custom_sample_sheet)))
    # Check the sample sheet for problems
    if check_and_warn(sample_sheet_file=custom_sample_sheet):
        logger.warning("Imported sample sheet has problems, see "
                       "warnings above")
