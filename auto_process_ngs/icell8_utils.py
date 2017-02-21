#!/usr/bin/env python
#
#     icell8_utils.py: utility functions for handling Wafergen iCell8 data
#     Copyright (C) University of Manchester 2017 Peter Briggs
#

"""
icell8_utils.py

Utility classes and functions for processing the outputs from Wafergen's
iCell8 platform:

- ICell8WellList: class representing iCell8 well list file

"""

#######################################################################
# Imports
#######################################################################

from bcftbx.TabFile import TabFile

######################################################################
# Classes
######################################################################

class ICell8WellList(object):
    """
    Class representing an iCell8 well list file

    The file is tab-delimited and consists of an uncommented header
    line which lists the fields ('Row','Col','Candidate',...),
    followed by lines of data.

    The key columns are 'Sample' (gives the cell type) and 'Barcode'
    (the inline barcode sequence).
    """
    def __init__(self,well_list_file):
        self._data = TabFile(filen=well_list_file,
                             first_line_is_header=True)
    def barcodes(self):
        """
        Return a list of barcodes
        """
        return [x['Barcode'] for x in self._data]
    def sample(self,barcode):
        """
        Return sample (=cell type) corresponding to barcode
        """
        samples = self._data.lookup('Barcode',barcode)
        try:
            return samples[0]['Sample']
        except IndexError:
            raise KeyError("Failed to locate sample for '%s'" % barcode)
