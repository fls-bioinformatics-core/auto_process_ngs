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
- ICell8ReadPair: class representing an iCell8 R1/R2 read-pair
- ICell8FastqIterator: class for iterating over iCell8 R1/R2 FASTQ-pair
"""

#######################################################################
# Imports
#######################################################################

from itertools import izip
from collections import Iterator
from bcftbx.FASTQFile import FastqIterator
from bcftbx.TabFile import TabFile

######################################################################
# Magic numbers
######################################################################

INLINE_BARCODE_LENGTH = 11
UMI_LENGTH = 10

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
        """
        Create a new ICell8WellList instance.

        Arguments:
          well_list_file (str): path to the well list
            file.
        """
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

class ICell8ReadPair(object):
    """
    Class representing an iCell8 R1/R2 read-pair
    """
    def __init__(self,r1,r2):
        """
        Create a new ICell8ReadPair instance.

        Arguments:
          r1 (FastqRead): the R1 read in the pair
          r2 (FasrqRead): the R2 read
        """
        if not r1.seqid.is_pair_of(r2.seqid):
            raise Exception("Reads are not paired")
        self._r1 = r1
        self._r2 = r2

    @property
    def r1(self):
        """
        R1 read from the pair
        """
        return self._r1

    @property
    def r2(self):
        """
        R2 read from the pair
        """
        return self._r2

    @property
    def barcode(self):
        """
        Inline barcode sequence extracted from the R1 read
        """
        return self._r1.sequence[0:INLINE_BARCODE_LENGTH]

    @property
    def umi(self):
        """
        UMI sequence extracted from the R1 read
        """
        return self._r1.sequence[INLINE_BARCODE_LENGTH:
                                 INLINE_BARCODE_LENGTH+UMI_LENGTH]

    @property
    def min_barcode_quality(self):
        """
        Minimum inline barcode quality score

        The score is encoded as a character e.g. '/' or 'A'.
        """
        return min(self._r1.quality[0:INLINE_BARCODE_LENGTH])

    @property
    def min_umi_quality(self):
        """
        Minimum UMI sequence quality score

        The score is encoded as a character e.g. '/' or 'A'.
        """
        return min(self._r1.quality[INLINE_BARCODE_LENGTH:
                                    INLINE_BARCODE_LENGTH+UMI_LENGTH])

class ICell8FastqIterator(Iterator):
    """
    Class for iterating over an iCell8 R1/R2 FASTQ-pair

    The iterator returns a set of ICell8ReadPair
    instances, for example:

    >>> for pair in ICell8FastqIterator(fq1,fq2):
    >>>   print "-- R1: %s" % pair.r1
    >>>   print "   R2: %s" % pair.r2
    """
    def __init__(self,fqr1,fqr2):
        """
        Create a new ICell8FastqIterator instance

        Arguments:
          fqr1 (str): path to the R1 FASTQ file
          fqr2 (str): path to the R2 FASTQ
        """
        self._fqr1 = FastqIterator(fqr1)
        self._fqr2 = FastqIterator(fqr2)
    def next(self):
        return ICell8ReadPair(self._fqr1.next(),
                              self._fqr2.next())
