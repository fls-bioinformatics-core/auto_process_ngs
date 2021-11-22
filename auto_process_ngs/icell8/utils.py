#!/usr/bin/env python
#
#     icell8.utils.py: utility functions for handling ICELL8 data
#     Copyright (C) University of Manchester 2017-2021 Peter Briggs
#

"""
icell8.utils.py

Utility classes and functions for processing the outputs from the ICELL8
single-cell platform.

Classes:

- ICell8WellList: class representing ICELL8 well list file
- ICell8Read1: class representing an ICELL8 R1 read
- ICell8ReadPair: class representing an ICELL8 R1/R2 read-pair
- ICell8FastqIterator: class for iterating over ICELL8 R1/R2 FASTQ-pair
- ICell8Stats: class for gathering stats from ICELL8 FASTQ pairs

Functions:

- get_batch_size: get optimal size for batches of reads
- batch_fastqs: split reads into batches
- normalize_sample_name: replace special characters in well list sample names
- get_bases_mask_icell8: generate bases mask for ICELL8 run
- get_bases_mask_icell8_atac: generate bases mask for ICELL8 ATAC-seq run
"""

#######################################################################
# Imports
#######################################################################

import os
import time
import logging
from collections import Iterator
from multiprocessing import Pool
from builtins import range
from bcftbx.FASTQFile import FastqIterator
from bcftbx.IlluminaData import SampleSheet
from bcftbx.IlluminaData import samplesheet_index_sequence
from bcftbx.IlluminaData import fix_bases_mask
from bcftbx.TabFile import TabFile
from ..bcl2fastq.utils import get_bases_mask
from ..command import Command
from ..fastq_utils import FastqReadCounter
from ..fastq_utils import pair_fastqs
from ..fastq_utils import get_read_count
from ..fastq_utils import get_read_number
from ..utils import ProgressChecker

# Initialise logging
import logging
logger = logging.getLogger(__name__)

#######################################################################
# Constants
#######################################################################

from . import constants
INLINE_BARCODE_LENGTH = constants.INLINE_BARCODE_LENGTH
UMI_LENGTH = constants.UMI_LENGTH
MAXIMUM_BATCH_SIZE = constants.MAXIMUM_BATCH_SIZE
SAMPLENAME_ILLEGAL_CHARS = constants.SAMPLENAME_ILLEGAL_CHARS

######################################################################
# Functions
######################################################################

def get_batch_size(fastqs,min_batches=1,
                   max_batch_size=MAXIMUM_BATCH_SIZE,
                   incr_function=None):
    """
    Determine number of reads per batch

    Given a maximum batch size (i.e. number of reads per
    batch), determine the number of batches and actual
    batch size.

    Arguments:
      fastqs (list): list of paths to one or more Fastq
        files to take reads from
      min_batches (int): initial minimum number of batches
        to try
      max_batch_size (int): the maxiumum batch size
      incr_function (Function): optional function to use
        to generate new number of batches to try

    Returns:
      Tuple: tuple of (batch_size,nbatches).
    """
    # Count the total number of reads
    print("Fetching read counts")
    nreads = get_read_count(fastqs)
    print("Total reads: %d" % nreads)

    # Default incrementer function: add the initial
    # number of batches on to get new number
    if incr_function is None:
        incr_function = lambda n: n + min_batches

    # Determine batch size
    batch_size = nreads//min_batches
    nbatches = min_batches
    print("Initial batch size: %d" % batch_size)
    print("Maximum batch size: %d" % max_batch_size)
    if max_batch_size > 0:
        while batch_size > max_batch_size or batch_size*nbatches < nreads:
            # Reset the number of batches
            nbatches = incr_function(nbatches)
            # Set the new batch size
            batch_size = nreads//nbatches
            if nreads%nbatches:
                batch_size += 1
            print("Trying %d batches: %d reads" % (nbatches,batch_size))
        logger.warning("Maximum batch size exceeded (%d), "
                       "increasing number of batches to %d"
                       % (max_batch_size,nbatches))
        logger.warning("New batch size: %d" % batch_size)
    print("Final batch size: %d" % batch_size)

    # Assert that all reads are covered with none left over
    # with no remainder
    assert(batch_size*nbatches >= nreads),"Batches won't cover all reads"
    return (batch_size,nbatches)

def batch_fastqs(fastqs,batch_size,basename="batched",
                 out_dir=None):
    """
    Splits reads from one or more Fastqs into batches

    Concatenates input Fastq files and then splits
    reads into smaller Fastqs using the external 'batch'
    utility.

    Arguments:
      fastqs (list): list of paths to one or more Fastq
        files to take reads from
      batch_size (int): number of reads to allocate to
        each batch
      basename (str): optional basename to use for the
        output Fastq files (default: 'batched')
      out_dir (str): optional path to a directory where
        the batched Fastqs will be written
    """
    # Determine number of batches
    nreads = get_read_count(fastqs)
    nbatches = nreads//batch_size
    if nbatches*batch_size < nreads:
        nbatches += 1
    print("Creating %d batches of %d reads" % (nbatches,
                                               batch_size))
    assert(batch_size*nbatches >= nreads)

    # Check if fastqs are compressed
    gzipped = fastqs[0].endswith('.gz')
    if gzipped:
        batch_cmd = Command('zcat')
    else:
        batch_cmd = Command('cat')

    # Get the read number
    read_number = get_read_number(fastqs[0])
    suffix = ".r%s.fastq" % read_number

    # Build and run the batching command
    batch_cmd.add_args(*fastqs)
    batch_cmd.add_args('|',
                       'split',
                       '-l',batch_size*4,
                       '-d',
                       '-a',3,
                       '--additional-suffix=%s' % suffix,
                       '-',
                       os.path.join(out_dir,"%s.B" % basename))
    batch_script = os.path.join(out_dir,"batch.sh")
    batch_cmd.make_wrapper_script("/bin/bash",
                                  batch_script)

    # Check for successful exit code
    retcode = Command("/bin/bash",
                      batch_script).run_subprocess(
                          working_dir=out_dir)
    if retcode != 0:
        raise Exception("Batching failed: exit code %s" % retcode)
    print("Batching completed")

    # Collect and return the batched Fastq names
    batched_fastqs = [os.path.join(out_dir,
                                   "%s.B%03d%s"
                                   % (basename,i,suffix))
                      for i in range(0,nbatches)]
    return batched_fastqs

def normalize_sample_name(s):
    """
    Clean up sample name from well list file

    Replaces 'illegal' characters in the supplied name
    with underscore characters and returns the
    normalized name.

    Arguments:
      s (str): sample name from well list file

    Returns:
      String: normalized sample name
    """
    name = []
    # Replace special characters with underscores
    for c in str(s):
        if c in SAMPLENAME_ILLEGAL_CHARS:
            name.append('_')
        else:
            name.append(c)
    return ''.join(name)

def get_bases_mask_icell8(bases_mask,sample_sheet=None):
    """
    Reset the supplied bases mask string so that only the
    bases containing the inline barcode and UMIs are kept,
    and any remaining bases are ignored.

    If a sample sheet is also supplied then an additional
    update will be made to ensure that the bases mask
    respects the barcode lengths given there.

    Arguments:
      bases_mask (str): initial bases mask string to update
      sample_sheet (str): path to optional sample sheet

    Returns:
      String: updated bases mask string
    """
    # Extract R1 mask
    bases_mask = bases_mask.split(',')
    r1_mask = bases_mask[0]
    # Update to restrict to 21 bases
    num_cycles = int(r1_mask[1:])
    icell8_inline_length = (INLINE_BARCODE_LENGTH + UMI_LENGTH)
    assert(num_cycles >= icell8_inline_length)
    discard_length = (num_cycles - icell8_inline_length)
    r1_mask = "y%d" % icell8_inline_length
    r1_mask += ("n%d" % discard_length if discard_length > 0 else "")
    bases_mask[0] = r1_mask
    # Rebuild full bases mask
    bases_mask = ','.join(bases_mask)
    # Handle sample sheet
    if sample_sheet is not None:
        index_seq = samplesheet_index_sequence(
            SampleSheet(sample_sheet).data[0])
        if index_seq is None:
            index_seq = ""
        bases_mask = fix_bases_mask(bases_mask,index_seq)
    return bases_mask

def get_bases_mask_icell8_atac(runinfo_xml):
    """
    Acquire a bases mask for ICELL8 scATAC-seq

    Generates an initial bases mask based on the run
    contents, and then updates this so that only the
    first 8 bases of each of index reads are used.

    Arguments:
      runinfo_xml (str): path to the RunInfo.xml for
        the sequencing run

    Returns:
      String: ICELL8 scATAC-seq bases mask string
    """
    # Get initial bases mask from run
    bases_mask = get_bases_mask(runinfo_xml).lower().split(',')
    # First read
    new_bases_mask = [bases_mask[0]]
    # Assume barcode is dual-index (8bp plus 8bp)
    # Update first index and second indices to restrict to 8 bases
    for index_mask in bases_mask[1:3]:
        num_cycles = int(index_mask[1:])
        if num_cycles < 8:
            raise Exception("Index read < 8 bases")
        new_bases_mask.append("I8%s" % ('n'*(num_cycles-8),))
    # Keep last read as is
    new_bases_mask.append(bases_mask[3])
    # Reassemble and return
    return ','.join(new_bases_mask)

def pass_quality_filter(s,cutoff):
    """
    Check if sequence passes quality filter cutoff

    Arguments:
      s (str): sequence quality scores (PHRED+33)
      cutoff (int): minimum quality value; all
        quality scores must be equal to or greater
        than this value for the filter to pass

    Returns:
      Boolean: True if the quality scores pass the
        filter cutoff, False if not.
    """
    cutoff = chr(cutoff + 33)
    for c in s:
        if c < cutoff:
            return False
    return True

######################################################################
# Classes
######################################################################

class ICell8WellList(object):
    """
    Class representing an ICELL8 well list file

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
    Class representing an ICELL8 R1 read
    """
    def __init__(self,fastq_read):
        """
        Create a new ICell8Read1 instance.

        Arguments:
          fastq_read (FastqRead): the R1 read in the pair
        """
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
    Class representing an ICELL8 R1/R2 read-pair
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
    Class for iterating over an ICELL8 R1/R2 FASTQ-pair

    The iterator returns a set of ICell8ReadPair
    instances, for example:

    >>> for pair in ICell8FastqIterator(fq1,fq2):
    >>>   print("-- R1: %s" % pair.r1)
    >>>   print("   R2: %s" % pair.r2)
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
    def __next__(self):
        self._read_count += 1
        r1 = self._fqr1.next()
        r2 = self._fqr2.next()
        try:
            return ICell8ReadPair(r1,r2)
        except Exception as ex:
            print("Failed to create read pair:")
            print("-- Read pair number: %d" % self._read_count)
            print("-- Read 1:\n%s" % r1)
            print("-- Read 2:\n%s" % r2)
            logging.critical("Failed to create read pair: %s" % ex)
            raise ex

class ICell8StatsCollector(object):
    """
    Class to collect ICELL8 barcode and UMI counts

    This class essentially wraps a single function
    which gets ICELL8 barcodes and distinct UMI
    counts from a Fastq file. It is used by the
    `Icell8Stats` class to collect counts for each
    file supplied.

    Example usage:

    >>> collector = ICell8StatsCollector()
    >>> fq,counts,barcodes = collector(fastq)

    By default the collection process is (relatively)
    quiet; more verbose output can be requested by
    setting the `verbose` argument to True on
    instantiation.

    The collector has been implemented as a callable
    class so that it can be used with both the built-in
    `map` function and `Pool.map` from the Python
    `multiprocessing` module. (Specifically, this
    implementation works around issues with
    `multiprocessing` being unable to pickle an
    instance method - otherwise we could use e.g.

    >>> pool.map(collector.collect_fastq_stats,fastqs)

    See the question at
    https://stackoverflow.com/q/1816958/579925
    and specifically the answer at
    https://stackoverflow.com/a/6975654/579925
    for more elaboration.)
    """
    def __init__(self,verbose=False):
        """
        Create a new ICell8StatsCollector instance

        Arguments:
          verbose (bool): if True then periodically
            reports progress to stdout (default:
            False)
        """
        self._verbose = bool(verbose)

    def __call__(self,fastq):
        return self.collect_fastq_stats(fastq)

    def collect_fastq_stats(self,fastq):
        """
        Get barcode and distinct UMI counts for Fastq file

        This method can be called directly, but is
        also invoked implicitly if its parent instance
        is called.

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
        print("collect_fastq_stats: started: %s" % fastq)
        try:
            n = FastqReadCounter.zcat_wc(fastq)
            print("%s: processing %d read%s" % (
                os.path.basename(fastq),
                n,('s' if n != 1 else '')))
            counts = {}
            umis = {}
            progress = ProgressChecker(percent=5,total=n)
            for i,r in enumerate(FastqIterator(fastq),start=1):
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
                if self._verbose:
                    if progress.check(i):
                        print("%s: %s: processed %d reads (%.1f%%)" %
                              (time.strftime("%Y%m%d.%H%M%S"),
                               os.path.basename(fastq),
                               i,progress.percent(i)))
        except Exception as ex:
            print("collect_fastq_stats: caught exception: '%s'" % ex)
            raise Exception("collect_fastq_stats: %s: caught exception "
                            "'%s'",(fastq,ex))
        print("collect_fastq_stats: returning: %s" % fastq)
        return (fastq,counts,umis)

class ICell8Stats(object):
    """
    Class for gathering statistics on ICELL8 FASTQ R1 files

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
          fastqs: set of paths to ICELL8 R1/R2 FASTQ
            file pairs to be processed
          nprocs (int): number of cores to use for
            statistics generation (default: 1)
          verbose (bool): if True then print additional
            output reporting progress of statistics
            gathering (default: don't report progress)
        """
        # Handle keywords
        nprocs = 1
        verbose = False
        for kw in kws:
            if kw not in ('nprocs','verbose'):
                raise TypeError("%s got an unexpected keyword "
                                "argument '%s'" %
                                (self.__class__.__name__,kw))
            if kw == 'nprocs':
                nprocs = int(kws['nprocs'])
            elif kw == 'verbose':
                verbose = bool(kws['verbose'])
        # Set up collector instance
        collector = ICell8StatsCollector(verbose=True)
        # Collect statistics for each file
        print("Collecting stats...")
        if nprocs > 1:
            # Multiple cores
            print("Multicore mode (%d cores)" % nprocs)
            pool = Pool(nprocs)
            results = pool.map(collector,fastqs)
            print("Processes completed, disposing of pool..")
            pool.close()
            pool.join()
            print("Pool disposal complete")
        else:
            # Single core
            print("Single core mode")
            results = map(collector,fastqs)
        # Combine results
        print("Merging stats from each Fastq:")
        self._counts = {}
        self._umis = {}
        for fq,fq_counts,fq_umis in results:
            print("- %s" % fq)
            nbarcodes = len(fq_counts)
            progress = ProgressChecker(percent=5,total=nbarcodes)
            print("  %d barcode%s" % (nbarcodes,
                                      ('s' if nbarcodes != 1
                                       else '')))
            for i,barcode in enumerate(fq_counts):
                try:
                    self._counts[barcode] += fq_counts[barcode]
                except KeyError:
                    self._counts[barcode] = fq_counts[barcode]
                try:
                    self._umis[barcode].update(fq_umis[barcode])
                except KeyError:
                    self._umis[barcode] = set(fq_umis[barcode])
                if verbose:
                    if progress.check(i):
                        print("  %d barcodes merged (%.1f%%)" %
                              (i,progress.percent(i)))
        nbarcodes = len(self._counts)
        print("Total %s barcode%s" % (nbarcodes,
                                      ('s' if nbarcodes != 1
                                       else '')))
        # Create unique sorted UMI lists
        print("Sorting UMI lists for each barcode")
        progress = ProgressChecker(percent=5,total=nbarcodes)
        for i,barcode in enumerate(self._umis):
            self._umis[barcode] = sorted(list(self._umis[barcode]))
            if verbose:
                if progress.check(i):
                    print("- UMIs sorted for %d barcodes (%.1f%%)" %
                          (i,progress.percent(i)))
        print("Finished stats collection")

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
