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
- ICell8Read1: class representing an iCell8 R1 read
- ICell8ReadPair: class representing an iCell8 R1/R2 read-pair
- ICell8FastqIterator: class for iterating over iCell8 R1/R2 FASTQ-pair
- ICell8Stats: class for gathering stats from iCell8 FASTQ pairs
- collect_fastq_stats: get barcode and distince UMI counts for Fastq
"""

#######################################################################
# Imports
#######################################################################

import time
import logging
from itertools import izip
from collections import Iterator
from multiprocessing import Pool
from bcftbx.FASTQFile import FastqIterator
from bcftbx.TabFile import TabFile
from .fastq_utils import pair_fastqs

# Initialise logging
import logging
logger = logging.getLogger(__name__)

######################################################################
# Magic numbers
######################################################################

INLINE_BARCODE_LENGTH = 11
UMI_LENGTH = 10

######################################################################
# Functions
######################################################################

def collect_fastq_stats(fastq):
    """
    Get barcode and distince UMI counts for Fastq file

    Used by Icell8Stats to collect counts for each file
    supplied.

    Arguments:
      fastq (str): path to Fastq file

    Returns:
      Tuple: tuple consisting of (fastq,counts,umis)
        where 'fastq' is the path to the input Fastq
        file, 'counts' is a dictionary with barcodes
        as keys and read counts as values, and 'umis'
        is a dictionary with barcodes as keys and
        sets of UMIs as values.
    """
    counts = {}
    umis = {}
    for r in FastqIterator(fastq):
        r = ICell8Read1(r)
        barcode = r.barcode
        try:
            counts[barcode] += 1
        except KeyError:
            counts[barcode] = 1
        umi = r.umi
        try:
            umis[barcode].add(umi)
        except KeyError:
            umis[barcode] = set((umi,))
    return (fastq,counts,umis)

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

    def samples(self):
        """
        Return a list of samples
        """
        samples = set([x['Sample'] for x in self._data])
        return sorted(list(samples))

    def sample(self,barcode):
        """
        Return sample (=cell type) corresponding to barcode
        """
        samples = self._data.lookup('Barcode',barcode)
        try:
            return samples[0]['Sample']
        except IndexError:
            raise KeyError("Failed to locate sample for '%s'" % barcode)

class ICell8Read1(object):
    """
    Class representing an iCell8 R1 read
    """
    def __init__(self,fastq_read):
        """
        Create a new ICell8Read1 instance.

        Arguments:
          fastq_read (FastqRead): the R1 read in the pair
        """
        self._barcode_length = 11
        self._umi_length = 10
        self._read = fastq_read
    @property
    def read(self):
        """
        R1 read
        """
        return self._read
    @property
    def barcode(self):
        """
        Inline barcode sequence extracted from the R1 read
        """
        return self._read.sequence[0:INLINE_BARCODE_LENGTH]
    @property
    def umi(self):
        """
        UMI sequence extracted from the R1 read
        """
        return self._read.sequence[INLINE_BARCODE_LENGTH:
                                   INLINE_BARCODE_LENGTH+UMI_LENGTH]
    @property
    def barcode_quality(self):
        """
        Inline barcode sequence quality extracted from the R1 read
        """
        return self._read.quality[0:INLINE_BARCODE_LENGTH]
    @property
    def umi_quality(self):
        """
        UMI sequence quality extracted from the R1 read
        """
        return self._read.quality[INLINE_BARCODE_LENGTH:
                                  INLINE_BARCODE_LENGTH+UMI_LENGTH]
    @property
    def min_barcode_quality(self):
        """
        Minimum inline barcode quality score

        The score is encoded as a character e.g. '/' or 'A'.
        """
        return min(self.barcode_quality)
    @property
    def min_umi_quality(self):
        """
        Minimum UMI sequence quality score

        The score is encoded as a character e.g. '/' or 'A'.
        """
        return min(self.umi_quality)

class ICell8ReadPair(ICell8Read1):
    """
    Class representing an iCell8 R1/R2 read-pair
    """
    def __init__(self,r1,r2):
        """
        Create a new ICell8ReadPair instance.

        Arguments:
          r1 (FastqRead): the R1 read in the pair
          r2 (FasrqRead): the matching R2 read
        """
        if not r1.seqid.is_pair_of(r2.seqid):
            raise Exception("Reads are not paired")
        ICell8Read1.__init__(self,r1)
        self._read2 = r2

    @property
    def r1(self):
        """
        R1 read from the pair
        """
        return self.read

    @property
    def r2(self):
        """
        R2 read from the pair
        """
        return self._read2

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
        self._read_count = 0
        self._fqr1 = FastqIterator(fqr1)
        self._fqr2 = FastqIterator(fqr2)
    def next(self):
        self._read_count += 1
        r1 = self._fqr1.next()
        r2 = self._fqr2.next()
        try:
            return ICell8ReadPair(r1,r2)
        except Exception as ex:
            print "Failed to create read pair:"
            print "-- Read pair number: %d" % self._read_count
            print "-- Read 1:\n%s" % r1
            print "-- Read 2:\n%s" % r2
            logging.critical("Failed to create read pair: %s" % ex)
            raise ex

class ICell8Stats(object):
    """
    Class for gathering statistics on iCell8 FASTQ R1 files

    Given a set of paths to FASTQ R1 files (from Icell8
    Fastq file pairs), collects statistics on the number of
    reads, barcodes and distinct UMIs.

    NB the list of distinct UMIs are where each UMI
    appears only once. Each UMI may appear multiple times
    across the FASTQ files.

    """
    def __init__(self,*fastqs,**kws):
        """
        Create a new ICell8Stats instance

        Arguments:
          fastqs: set of paths to ICell8 R1/R2 FASTQ
            file pairs to be processed
          nprocs (int): number of cores to use for
            statistics generation (default: 1)
        """
        # Handle keywords
        nprocs = 1
        for kw in kws:
            if kw not in ('nprocs',):
                raise TypeError("%s got an unexpected keyword "
                                "argument '%s'" %
                                (self.__class__.__name__,kw))
            if kw == 'nprocs':
                nprocs = int(kws['nprocs'])
        print "#procs = %s" % nprocs
        # Collect statistics for each file
        print "Collecting stats..."
        if nprocs > 1:
            # Multiple cores
            pool = Pool(nprocs)
            results = pool.map(collect_fastq_stats,fastqs)
            pool.close()
            pool.join()
        else:
            # Single core
            results = map(collect_fastq_stats,fastqs)
        # Combine results
        print "Merging stats from each Fastq:"
        self._counts = {}
        self._umis = {}
        for fq,fq_counts,fq_umis in results:
            print "%s" % fq
            for barcode in fq_counts:
                try:
                    self._counts[barcode] += fq_counts[barcode]
                except KeyError:
                    self._counts[barcode] = fq_counts[barcode]
                try:
                    self._umis[barcode].update(fq_umis[barcode])
                except KeyError:
                    self._umis[barcode] = set(fq_umis[barcode])
        # Create unique sorted UMI lists
        for barcode in self._umis:
            self._umis[barcode] = sorted(list(self._umis[barcode]))

    def barcodes(self):
        """
        Return list of barcodes from the FASTQs
        """
        return sorted([b for b in self._counts.keys()])

    def nreads(self,barcode=None):
        """
        Return total number of reads, or per barcode

        Invoked without arguments, returns the
        total number of reads analysed. If a barcode
        is specified then returns the number of reads
        with that barcode.

        Arguments:
          barcode (str): optional, specify barcode
            for which the read count will be returned.

        Returns:
          Integer: number of reads.
        """
        if barcode is not None:
            return self._counts[barcode]
        else:
            nreads = 0
            for b in self.barcodes():
                nreads += self.nreads(b)
            return nreads

    def distinct_umis(self,barcode=None):
        """
        Return all distinct UMIs, or by barcode

        Invoked without arguments, returns a list
        of distinct UMIs found across the files. If a
        barcode is specified then returns a list of
        UMIs associated with that barcode.

        Arguments:
          barcode (str): optional, specify barcode
            for which the list of distinct UMIs will be
            returned.

        Returns:
          List: list of distinct UMI sequences.
        """
        if barcode is not None:
            return sorted(list(self._umis[barcode]))
        else:
            umis = set()
            for b in self.barcodes():
                umis.update(self.distinct_umis(b))
            return sorted(list(umis))
