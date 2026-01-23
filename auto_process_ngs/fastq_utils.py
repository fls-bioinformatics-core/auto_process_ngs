#!/usr/bin/env python
#
#     fastq_utils.py: utility functions for operating on fastq files
#     Copyright (C) University of Manchester 2016-2026 Peter Briggs
#
########################################################################
#
# fastq_utils.py
#
#########################################################################

"""
Utility classes and functions for operating on Fastq files:

- BaseFastqAttrs: base class for extracting info from Fastq file
- IlluminaFastqAttrs: class for extracting info from Illumina Fastqs
- FastqReadCounter: implements various methods for counting reads
  in FASTQ files
- assign_barcodes_single_end: extract and assign inline barcodes
- get_read_number: get the read number (1 or 2) from a Fastq file
- get_read_count: count total reads across one or more Fastqs
- pair_fastqs: automagically pair up FASTQ files
- pair_fastqs_by_name: pair up FASTQ files based on their names
- group_fastqs_by_name: group FASTQ files based on their names
  (more general version of 'pair_fastqs_by_name' which can handle
  arbitrary collections of read IDs)
- remove_index_fastqs: remove index (I1/I2) Fastqs from a list
- build_custom_fastq_attrs_regex: make regular expression patterns
  and format templates for custom 'FastqAttrs' classes
- get_custom_fastq_attrs_class: create a custom 'FastqAttrs' class
  for handling non-canonical Fastq filenames.
"""

#######################################################################
# Imports
#######################################################################

import os
import re
import gzip
import subprocess
import logging
from bcftbx.FASTQFile import FastqIterator
from bcftbx.FASTQFile import get_fastq_file_handle
from bcftbx.FASTQFile import nreads

#######################################################################
# Classes
#######################################################################

class BaseFastqAttrs:
    """
    Base class for extracting information about a Fastq file

    Instances of this class provide the follow attributes:

    fastq:            the original fastq file name
    basename:         basename with path and NGS extensions stripped
    extension:        full extension e.g. '.fastq.gz'
    type:             file type ('fastq' or 'bam')
    compression:      compression type ('gz' or 'bz2')
    sample_name:      name of the sample
    sample_number:    integer (or None if no sample number)
    barcode_sequence: barcode sequence (string or None)
    lane_number:      integer (or None if no lane number)
    read_number:      integer (or None if no read number)
    set_number:       integer (or None if no set number)
    is_index_read:    boolean (True if index read, False if not)

    The 'fastq', 'basename', 'extension', 'type' and 'compression',
    attributes are set by this class; the remaining attributes
    (at minimum 'sample_name') should be set by the subclass.

    Subclasses should also implement the 'fastq_basename' and
    'bam_basename' methods.
    """
    def __init__(self, fastq):
        # Store name
        self.fastq = fastq
        # Basename and extensions
        self.basename = None
        self.extension = None
        self.type = None
        self.compression = None
        self._set_basename_and_extensions()
        # Values that should be derived from the name
        # (should be set by subclass)
        self.sample_name = None
        self.sample_number = None
        self.barcode_sequence = None
        self.lane_number = None
        self.read_number = None
        self.set_number = None
        self.is_index_read = False

    def _set_basename_and_extensions(self):
        """
        Set basename and extract extensions
        """
        basename = os.path.basename(self.fastq)
        for ext in (".gz", ".bz2"):
            # Strip compression extension
            if basename.endswith(ext):
                basename = basename[:-len(ext)]
                self.compression = ext[1:]
                break
        for ext in (".fastq", ".fq", ".bam"):
            if basename.endswith(ext):
                basename = basename[:-len(ext)]
                self.type = ext[1:]
                break
        if self.type in ("fastq", "fq"):
            self.type = "fastq"
        elif self.type in ("bam",):
            self.type = "bam"
        self.basename = basename
        self.extension = os.path.basename(self.fastq)[len(self.basename):]

    def fq_attrs(self):
        """
        Return attributes as a dictionary
        """
        return {
            "sample_name": self.sample_name,
            "sample_number": self.sample_number,
            "barcode_sequence": self.barcode_sequence,
            "lane_number": self.lane_number,
            "read_number": self.read_number,
            "set_number": self.set_number,
            "is_index_read": self.is_index_read,
        }

    def fastq_basename(self, fq_attrs=None):
        """
        Return the basename of the Fastq file
        """
        raise NotImplementedError("'fastq_basename' method not implemented")

    def bam_basename(self, fq_attrs=None):
        """
        Return basename for an associated BAM file

        Should be overridden by subclasses.
        """
        raise NotImplementedError("'bam_basename' method not implemented")

    def __repr__(self):
        try:
            return self.fastq_basename()
        except NotImplementedError:
            # Default
            return self.basename


class IlluminaFastqAttrs(BaseFastqAttrs):
    """Class for extracting information about Fastq files

    Given the name of a Fastq file, extract data about the sample name,
    barcode sequence, lane number, read number and set number.

    The name format can be a 'full' Fastq name as generated by CASAVA or
    bcl2fastq 1.8, which follows the general form:

    <sample_name>_<barcode_sequence>_L<lane_number>_R<read_number>_<set_number>.astq.gz

    e.g. for

    NA10831_ATCACG_L002_R1_001.fastq.gz

    sample_name = 'NA10831_ATCACG_L002_R1_001'
    barcode_sequence = 'ATCACG'
    lane_number = 2
    read_number = 1
    set_number = 1

    Alternatively it can be a full Fastq name as generated by bcl2fastq2,
    of the general form:

    <sample_name>_S<sample_number>_L<lane_number>_R<read_number>_001.fastq.gz

    e.g. for

    ES_exp1_S4_L003_R2_001.fastq.gz

    sample_name = 'ES_exp1'
    sample_number = 4
    lane_number = 3
    read_number = 2
    set_number = 1

    bcl2fastq can also produce 'index read' Fastq files where the
    R1/R2 is replaced by I1, e.g.:

    ES_exp1_S4_L003_I1_001.fastq.gz

    Alternatively it can be a 'reduced' version where one or more
    of the components has been omitted (typically because they are
    redundant in uniquely distinguishing a Fastq file within a
    set of Fastqs).

    The reduced formats are:

    <sample_name>
    <sample_name>_L<lane_number>
    <sample_name>_<barcode_sequence>
    <sample_name>_<barcode_sequence>_L<lane_number>

    with an optional suffix '_R<read_number>' for paired end sets.

    e.g.

    NA10831
    NA10831_L002
    NA10831_ATCACG
    NA10831_ATCACG_L002

    Finally, the name can be a non-standard name of the form:

    <sample_name>.r<read_number>

    or

    <sample_name>.<barcode_sequence>.r<read_number>

    In this case the sample_names are permitted to include dots.

    Provides the follow attributes:

    fastq:            the original fastq file name
    sample_name:      name of the sample (leading part of the name)
    sample_number:    integer (or None if no sample number)
    barcode_sequence: barcode sequence (string or None)
    lane_number:      integer (or None if no lane number)
    read_number:      integer (or None if no read number)
    set_number:       integer (or None if no set number)
    is_index_read:    boolean (True if index read, False if not)
    """

    def __init__(self,fastq):
        """Create and populate a new IlluminaFastqAttrs object

        Arguments:
          fastq: name of the fastq.gz (optionally can include leading path)
        """
        BaseFastqAttrs.__init__(self,fastq)
        # Base name for sample (no leading path or extension)
        fastq_base = self.basename
        # Determine if it's a non-standard (dot-separated) name
        #
        # These are names of the form e.g.
        # NH1.2.r2
        # or
        # NH1_ChIP-seq.ACAGTG.r2
        if fastq_base.count('.') > 0:
            self.delimiter = '.'
            fields = fastq_base.split('.')
            field = fields[-1]
            if len(field) == 2 and field.startswith('r'):
                # Read number
                self.read_number = int(field[1])
                fields = fields[:-1]
                field = fields[-1]
            if len(fields) > 1:
                # Barcode sequence
                # Determine separator for dual index barcodes ('-' or '+')
                for sep in ('-','+'):
                    if sep in field:
                        break
                is_tag = True
                for f in field.split(sep):
                    for c in f:
                        is_tag = is_tag and c in 'ACGTN'
                if is_tag:
                    self.barcode_sequence = field
                    fields = fields[:-1]
                    field = fields[-1]
            # Remaining fields are the sample name
            self.sample_name = '.'.join(fields)
            assert(self.sample_name != '')
            return
        # Some form of Illumina-derived name
        #
        # Full Illumina-style names are e.g.
        # NH1_ChIP-seq_Gli1_ACAGTG_L001_R1_001
        # or
        # NH1_ChIP-seq_Gli1_S4_L003_R2_001
        #
        # We have shorter name formats where redundant parts are
        # omitted, the patterns are:
        # NAME          e.g. NH1_ChIP-seq_Gli1
        # NAME+LANE     e.g. NH1_ChIP-seq_Gli1_L001
        # NAME+TAG      e.g. NH1_ChIP-seq_Gli1_ACAGTG
        # NAME+TAG+LANE e.g. NH1_ChIP-seq_Gli1_ACAGTG_L001
        #
        # Also read number (i.e. R1 or R2) is appended but only for
        # paired end samples
        #
        # The set number is never included, except for full names
        fields = fastq_base.split('_')
        self.delimiter = '_'
        # Deal with set number first e.g. 001
        field = fields[-1]
        ##logger.debug("Test for set number %s" % field)
        if len(field) == 3 and field.isdigit():
            self.set_number = int(field)
            fields = fields[:-1]
        # Deal with trailing read number e.g. R1
        field = fields[-1]
        ##logger.debug("Test for read number %s" % field)
        if len(field) == 2:
            if field.startswith('R'):
                self.read_number = int(field[1])
                fields = fields[:-1]
            elif field.startswith('I'):
                self.read_number = int(field[1])
                self.is_index_read = True
                fields = fields[:-1]
        # Deal with trailing lane number e.g. L001
        field = fields[-1]
        ##logger.debug("Test for lane number %s" % field)
        if len(field) == 4 and field.startswith('L') and field[1:].isdigit():
            self.lane_number = int(field[1:])
            fields = fields[:-1]
        # Deal with trailing index tag e.g. ATTGCT or ATTGCT-CCTAAG
        field = fields[-1]
        ##logger.debug("Test for barcode sequence %s" % field)
        if len(fields) > 1:
            # This mustn't be the last field: if it is then it's
            # not the tag - it's the name
            # Determine separator for dual index barcodes ('-' or '+')
            for sep in ('-','+'):
                if sep in field:
                    break
            is_tag = True
            for f in field.split(sep):
                for c in f:
                    is_tag = is_tag and c in 'ACGTN'
            if is_tag:
                self.barcode_sequence = field
                fields = fields[:-1]
                ##logger.debug("Identified barcode sequence as %s" % self.barcode_sequence)
            else:
                # Alternatively might be the sample number
                if field.startswith('S'):
                    try:
                        if field[1:].isdigit():
                            self.sample_number = int(field[1:])
                            fields = fields[:-1]
                    except IndexError:
                        pass
        # What's left is the name
        ##logger.debug("Remaining fields: %s" % fields)
        self.sample_name = '_'.join(fields)
        assert(self.sample_name != '')

    def fastq_basename(self, fq_attrs=None):
        """
        Reconstruct basename for Fastq files
        """
        if fq_attrs is None:
            fq_attrs = self.fq_attrs()
        delimiter = self.delimiter
        fq = ["%s" % fq_attrs["sample_name"]]
        if fq_attrs["sample_number"] is not None:
            fq.append("S%d" % fq_attrs["sample_number"])
        if fq_attrs["barcode_sequence"] is not None:
            fq.append("%s" % fq_attrs["barcode_sequence"])
        if fq_attrs["lane_number"] is not None:
            fq.append("L%03d" % fq_attrs["lane_number"])
        if fq_attrs["read_number"] is not None:
            if delimiter == '.':
                fq.append("r%d" % fq_attrs["read_number"])
            elif fq_attrs["is_index_read"]:
                fq.append("I%d" % fq_attrs["read_number"])
            else:
                fq.append("R%d" % fq_attrs["read_number"])
        if fq_attrs["set_number"] is not None:
            fq.append("%03d" % fq_attrs["set_number"])
        return delimiter.join(fq)

    def bam_basename(self, fq_attrs=None):
        """
        Construct basename for BAM files
        """
        if fq_attrs is None:
            fq_attrs = self.fq_attrs()
        else:
            fq_attrs = fq_attrs.copy()
        fq_attrs["read_number"] = None
        return self.fastq_basename(fq_attrs)


class FastqReadCounter:
    """
    Implements various methods for counting reads in FASTQ file

    The methods are:

    - simple: a wrapper for the FASTQFile.nreads() function
    - fastqiterator: counts reads using FASTQFile.FastqIterator
    - zcat_wc: runs 'zcat | wc -l' in the shell
    - reads_per_lane: counts reads by lane using FastqIterator

    """
    @staticmethod
    def simple(fastq=None,fp=None):
        """
        Return number of reads in a FASTQ file

        Uses the FASTQFile.nreads function to do the counting.

        Arguments:
          fastq: fastq(.gz) file
          fp: open file descriptor for fastq file

        Returns:
          Number of reads

        """
        return nreads(fastq=fastq,fp=fp)
    @staticmethod
    def fastqiterator(fastq=None,fp=None):
        """
        Return number of reads in a FASTQ file

        Uses the FASTQFile.FastqIterator class to do the
        counting.

        Arguments:
          fastq: fastq(.gz) file
          fp: open file descriptor for fastq file

        Returns:
          Number of reads

        """
        nreads = 0
        for r in FastqIterator(fastq_file=fastq,fp=fp):
            nreads += 1
        return nreads
    @staticmethod
    def zcat_wc(fastq=None,fp=None):
        """
        Return number of reads in a FASTQ file

        Uses a system call to run 'zcat FASTQ | wc -l' to do
        the counting (or just 'wc -l' if not a gzipped FASTQ).

        Note that this can only operate on fastq files (not
        on streams provided via the 'fp' argument; this will
        raise an exception).

        Arguments:
          fastq: fastq(.gz) file
          fp: open file descriptor for fastq file

        Returns:
          Number of reads

        """
        if fastq is None:
            raise Exception("zcat_wc: can only operate on a file")
        if fastq.endswith(".gz"):
            cmd = "zcat %s | wc -l" % fastq
        else:
            cmd = "wc -l %s | cut -d' ' -f1" % fastq
        output = subprocess.check_output(cmd,shell=True)
        try:
            return int(output)//4
        except Exception as ex:
            raise Exception("zcat_wc returned: %s" % output)
    @staticmethod
    def reads_per_lane(fastq=None,fp=None):
        """
        Return counts of reads in each lane of FASTQ file

        Uses the FASTQFile.FastqIterator class to do the
        counting, with counts split by lane.

        Arguments:
          fastq: fastq(.gz) file
          fp: open file descriptor for fastq file

        Returns:
          Dictionary where keys are lane numbers (as integers)
            and values are number of reads in that lane.

        """
        nreads = {}
        for r in FastqIterator(fastq_file=fastq,fp=fp):
            lane = int(r.seqid.flowcell_lane)
            try:
                nreads[lane] += 1
            except KeyError:
                nreads[lane] = 1
        return nreads

#######################################################################
# Functions
#######################################################################

def assign_barcodes_single_end(fastq_in,fastq_out,n=5):
    """
    Extract inline barcodes and assign to Fastq read headers

    Strips the first n bases from each read of the input
    FASTQ file and assigns it to the index sequence for that
    read in the output file.

    If the supplied output file name ends with '.gz' then it
    will be gzipped.

    Arguments:
      fastq_in (str): input FASTQ file (can be gzipped)
      fastq_out (str): output FASTQ file (will be gzipped if
        ending with '.gz')
      n (integer): number of bases to extract and assign as
        index sequence (default: 5)

    Returns:
      Integer: number of reads processed.

    """
    if fastq_out.endswith('.gz'):
        fp = gzip.GzipFile(filename=fastq_out,mode='wb')
    else:
        fp = open(fastq_out,'w')
    print("Processing reads from %s" % fastq_in)
    nread = 0
    for read in FastqIterator(fastq_in):
        # Extract new barcode sequence
        barcode = read.sequence[:n]
        # Truncate sequence and quality accordingly
        sequence = read.sequence[n:]
        quality = read.quality[n:]
        # Assign new values and write to output
        read.seqid.index_sequence = barcode
        read.sequence = sequence
        read.quality = quality
        fp.write("%s\n" % read)
        nread += 1
    fp.close()
    print("Finished (%d reads processed)" % nread)
    return nread

def get_read_number(fastq):
    """
    Get the read number (1 or 2) from a Fastq file

    Arguments:
      fastq (str): path to a Fastq file

    Returns:
      Integer: read number (1 or 2) extracted from the first read.
    """
    with get_fastq_file_handle(fastq) as fp:
        for r in FastqIterator(fp=fp):
            seq_id = r.seqid
            break
    return int(seq_id.pair_id)

def get_read_count(fastqs):
    """
    Get the total count of reads across multiple Fastqs

    Arguments:
      fastqs (list): lpaths to one or more Fastq files

    Returns:
      Integer: total number of reads across all files.
    """
    nreads = 0
    for fq in fastqs:
        n = FastqReadCounter.zcat_wc(fq)
        print("%s:\t%d" % (os.path.basename(fq),n))
        nreads += n
    return nreads

def pair_fastqs(fastqs):
    """
    Automagically pair up FASTQ files

    Given a list of FASTQ files, generate a list of R1/R2
    pairs by examining the header for the first read in
    each file.

    Arguments:
      fastqs (list): list of paths to FASTQ files which
        will be paired.

    Returns:
      Tuple: pair of lists of the form (paired,unpaired),
        where `paired` is a list of tuples consisting of
        FASTQ R1/R2 pairs and `unpaired` is a list of
        FASTQs which couldn't be paired.
    """
    fq_pairs = []
    seq_ids = {}
    bad_files = []
    for fq in [os.path.abspath(fq) for fq in fastqs]:
        # Get header from first read
        seq_id = None
        with get_fastq_file_handle(fq) as fp:
            for r in FastqIterator(fp=fp):
                seq_id = r.seqid
                break
        if seq_id is None:
            logging.debug("'Bad' file: %s" % fq)
            bad_files.append(fq)
            continue
        fq_pair = None
        for fq1 in seq_ids:
            if seq_id.is_pair_of(seq_ids[fq1]):
                # Found a pair
                if seq_id.pair_id == '1':
                    fq_pair = (fq,fq1)
                else:
                    fq_pair = (fq1,fq)
                fq_pairs.append(fq_pair)
                logging.debug("*** Paired: %s\n"
                              "          : %s" % fq_pair)
                # Remove paired fastq
                del(seq_ids[fq1])
                break
        if fq_pair is None:
            # Unable to pair, store for now
            logging.debug("Unpaired: %s" % fq)
            seq_ids[fq] = seq_id
    # Sort pairs into order
    fq_pairs = sorted(fq_pairs,key=lambda x: x[0])
    unpaired = sorted(list(seq_ids.keys()) + bad_files)
    # Return paired and upaired fastqs
    return (fq_pairs,unpaired)

def pair_fastqs_by_name(fastqs,fastq_attrs=IlluminaFastqAttrs):
    """
    Pair Fastq files based on their name

    Pairing is based on the read number for the supplied
    Fastq files being present in the file names; the file
    contents are not examined.

    Unpaired Fastqs (i.e. those for which a mate cannot be
    found) are returned as a "pair" where the equivalent R1
    or R2 mate is missing.

    Arguments:
      fastqs (list): list of Fastqs to pair
      fastq_attrs (BaseFastqAttrs): optional, class to use
        for extracting data from the filename (default:
        IlluminaFastqAttrs)

    Returns:
      List: list of tuples (R1,R2) with the R1/R2 pairs,
       or (R1,) or (R2,) for unpaired files.
    """
    pairs = []
    fastqs_r1 = sorted(list(filter(lambda f:
                                   fastq_attrs(f).read_number != 2,fastqs)))
    fastqs_r2 = sorted(list(filter(lambda f:
                                   fastq_attrs(f).read_number == 2,fastqs)))
    for fqr1 in fastqs_r1:
        # Split up R1 name
        logging.debug("fqr1 %s" % os.path.basename(fqr1))
        dir_path = os.path.dirname(fqr1)
        # Generate equivalent R2 file
        fqr2 = fastq_attrs(fqr1)
        fqr2.read_number = 2
        fqr2 = os.path.join(dir_path,"%s%s" % (fqr2,fqr2.extension))
        logging.debug("fqr2 %s" % os.path.basename(fqr2))
        if fqr2 in fastqs_r2:
            pairs.append((fqr1,fqr2))
        else:
            pairs.append((fqr1,))
    # Looking for unpaired R2 files
    for fqr2 in fastqs_r2:
        for pair in pairs:
            try:
                if fqr2 == pair[1]:
                    fqr2 = None
                    break
            except IndexError:
                pass
        if fqr2 is not None:
            pairs.append((fqr2,))
    pairs = sorted(pairs,key=lambda x: x[0])
    return pairs

def group_fastqs_by_name(fastqs,fastq_attrs=IlluminaFastqAttrs):
    """
    Group Fastq files based on their name

    Grouping is based on the read number and type for the
    supplied Fastq files being present in the file names; the
    file contents are not examined.

    Unpaired Fastqs (i.e. those for which a mate cannot be
    found) are returned as a "pair" where the equivalent R1
    or R2 mate is missing.

    Arguments:
      fastqs (list): list of Fastqs to pair
      fastq_attrs (BaseFastqAttrs): optional, class to use
        for extracting data from the filename (default:
        IlluminaFastqAttrs)

    Returns:
      List: list of tuples (R1,R2) with the R1/R2 pairs,
        or (R1,) or (R2,) for unpaired files.
    """
    # Get reads and put into groups by read ID
    reads = set()
    index_reads = set()
    fastq_sets = dict()
    for fastq in fastqs:
        fq = fastq_attrs(fastq)
        read_number = fq.read_number if fq.read_number else 1
        if not fq.is_index_read:
            read = "r%d" % read_number
            reads.add(read)
        else:
            read = "i%d" % read_number
            index_reads.add(read)
        try:
            fastq_sets[read].append(fastq)
        except KeyError:
            fastq_sets[read] = [fastq]
    reads = sorted(list(reads)) + sorted(list(index_reads))
    # Rearrange into groups
    groups = []
    for ii,read in enumerate(reads):
        for fastq in fastq_sets[read]:
            # Create reference Fastq name
            fq_ref = fastq_attrs(fastq)
            fq_ref.is_index_read = False
            fq_ref.read_number = 1
            # Initialise a new group
            group = [fastq]
            # Look for matches
            for r in reads[ii+1:]:
                unmatched_fastqs = list()
                for fastq1 in fastq_sets[r]:
                    fq1_ref = fastq_attrs(fastq1)
                    fq1_ref.is_index_read = False
                    fq1_ref.read_number = 1
                    if str(fq1_ref) == str(fq_ref):
                        group.append(fastq1)
                    else:
                        unmatched_fastqs.append(fastq1)
                fastq_sets[r] = unmatched_fastqs
            groups.append(group)
    groups = sorted(groups,key=lambda x: x[0])
    return groups

def remove_index_fastqs(fastqs,fastq_attrs=IlluminaFastqAttrs):
    """
    Remove index (I1/I2) Fastqs from list

    Arguments:
      fastqs (list): list of paths to Fastq files
      fastq_attrs (BaseFastqAttrs): class to use for
        extracting attributes from Fastq names
        (defaults to IlluminaFastqAttrs)

    Returns:
      List: input Fastq list with any index read
        Fastqs removed.
    """
    return list(filter(lambda fq:
                       not fastq_attrs(fq).is_index_read,
                       fastqs))


def build_custom_fastq_attrs_regex(pattern):
    """
    Build regex pattern and string templates for Fastq filenames

    Given a glob-like pattern describing a Fastq file name
    format, returns a tuple of (REGEX_PATTERN, FORMAT_STRING),
    which can be used to extract attributes from a Fastq file
    name and regenerate the name using those attributes.

    Patterns are strings which should include the elements
    '{SAMPLE}' and '{READ}' along with constant characters and
    wildcard elements '*'.

    For example:

    ::

        {SAMPLE}_*_{READ}

    would generate a regular expression and template for
    matching and regenerating file names of the form:

    ::

        PJB1_1.fastq

    where the sample name would be 'PJB1' and the read
    number would be '1'.

    Arguments:
        pattern (str): glob-like pattern describing a Fastq file
        name format

    Returns:
        Tuple: a tuple consisting of regular expression pattern
        and string template derived from the input pattern.
    """
    # Break the pattern into tokens and fixed elements
    parts = []
    for item in pattern.split("}"):
        if not item:
            continue
        elif item.startswith("{"):
            # Token e.g. '{READ}'
            parts.append(item + "}")
        elif "{" in item:
            # Fixed element + trailing token
            parts.extend([item.split("{")[0], "{" + item.split("{")[1] + "}"])
        else:
            # Fixed element without trailing token
            parts.append(item)
    # Build the regex pattern and format string
    re_pattern = []
    format_str = []
    idx = 0
    for part in parts:
        if part.startswith("{") and part.endswith("}"):
            # Token
            token = part[1:-1]
            if token == "SAMPLE":
                re_pattern.append("(?P<sample_name>.+)")
                format_str.append("{sample_name}")
            elif token == "READ":
                re_pattern.append("(?P<read_number>[1-3])")
                format_str.append("{read_number}")
            else:
                raise Exception(f"Unrecognised token '{item}'")
        else:
            # Interstitial string
            if "*" in part:
                # Contains wildcards - break into subparts
                # and replace with appropriate regex matches
                subparts = part.split("*")
                for subpart in subparts[:-1]:
                    idx += 1
                    re_pattern.append(f"{subpart}(?P<p{idx}>.*)")
                    format_str.append("%s{p%s}" % (subpart, idx))
                if subparts[-1]:
                    re_pattern.append(f"{subparts[-1]}")
                    format_str.append(subparts[-1])
            else:
                re_pattern.append(part)
                format_str.append(part)
    re_pattern = f"^{''.join(re_pattern)}$"
    format_str = "".join(format_str)
    return (re_pattern, format_str)


def get_custom_fastqattrs_class(pattern):
    """
    Create and return a custom FastqAttrs class

    This is a factory function that can be used to create
    custom FastqAttr-type classes (subclassed from the
    ``BaseFastqAttr`` class) to extract basic sample name
    and read number information from non-canonical Fastq
    file names.

    The function takes a single argument which is a basic
    glob-like pattern describing the expected filename.

    Patterns should include the strings ``{SAMPLE}`` in
    the position where the sample name is expected, and
    ``{READ}`` where the read number is expected. The
    rest of the pattern should consist of some combination
    of the "fixed" (i.e. invariant) parts of the name
    and the "*" for the variable parts which are not
    either sample name or read number.

    Some examples:

    - if the files are "S1-1_1.fastq", "S1-1_2.fastq", etc
      then the pattern "{SAMPLE}_{READ}" will enable
      name and read number extraction
    - if the files are "S1-1_R1.fastq", "S1-1_R2.fastq",
      etc then the pattern could be "{SAMPLE}_R{READ}"
    - for simple Illumina-style names (e.g.
      "SMP1-1_S1_L002_R2_001.fastq"
      the pattern "{SAMPLE}_S*_L*_R{READ}_001" could be
      a suitable pattern

    If the pattern includes a recognised Fastq extension
    then this will be ignored.

    Arguments:
      pattern (str): basic glob-like pattern to
        define sample and read number

    Returns:
      Class: custom subclass of ``BaseFastqAttrs`` for
      parsing files with the specified name structure.
    """
    # Fetch the Fastq regular expression and template
    fq_re_pattern, fq_format_str = build_custom_fastq_attrs_regex(pattern)

    # Transform Fastq regex pattern for BAM names
    bam_re_pattern = fq_re_pattern
    for read_pattern in ("_R(?P<read_number>[1-3])",
                         "_(?P<read_number>[1-3])",
                         "(?P<read_number>[1-3])"):
        # Remove the read pattern
        bam_re_pattern = bam_re_pattern.replace(read_pattern, "")
    # Strip any trailing 'separator' characters
    bam_re_pattern = bam_re_pattern.rstrip("_-.")

    # Transform Fastq format template for BAM names
    bam_format_str = fq_format_str
    for read_pattern in ("_R{read_number}",
                         "_{read_number}",
                         "{read_number}"):
        # Remove the read pattern
        bam_format_str = bam_format_str.replace(read_pattern, "")

    # Compile the regular expressions
    fq_re_pattern = re.compile(fq_re_pattern)
    bam_re_pattern = re.compile(bam_re_pattern)

    # Create the custom class
    class CustomFastqAttrs(BaseFastqAttrs):
        # Raw pattern
        fq_pattern = pattern
        # Regular expressions
        _fq_re_pattern = fq_re_pattern
        _bam_re_pattern = bam_re_pattern
        # Format strings
        _fq_format_string = fq_format_str
        _bam_format_string = bam_format_str

        def __init__(self, fastq):
            BaseFastqAttrs.__init__(self, fastq)
            if self.type == "bam":
                # Try matching as a BAM file name
                self._fq = self._bam_re_pattern.match(self.basename)
            else:
                # Default is to match as FASTQ name
                self._fq = self._fq_re_pattern.match(self.basename)
            if self._fq is not None:
                fq_attrs = self._fq.groupdict()
                self.sample_name = fq_attrs['sample_name']
                if "read_number" in fq_attrs and fq_attrs['read_number']:
                    self.read_number = int(fq_attrs['read_number'])

        def fastq_basename(self, fq_attrs=None):
            if self._fq is not None:
                # Get attributes
                if not fq_attrs:
                    fq_attrs = self.fq_attrs()
                # If no read number then fallback to BAM
                if fq_attrs["read_number"] is None:
                    return self.bam_basename(fq_attrs=fq_attrs)
                # Add in values for the interstitial components
                for name in self._fq.groupdict():
                    if name.startswith("p"):
                        fq_attrs[name] = self._fq.groupdict()[name]
                # Make basename
                return self._fq_format_string.format(**fq_attrs)
            else:
                return self.basename

        def bam_basename(self, fq_attrs=None):
            if self._fq is not None:
                if not fq_attrs:
                    fq_attrs = self.fq_attrs()
                # Add in values for the interstitial components
                for name in self._fq.groupdict():
                    if name.startswith("p"):
                        fq_attrs[name] = self._fq.groupdict()[name]
                # Make basename
                return self._bam_format_string.format(**fq_attrs)
            else:
                return self.basename

        def __repr__(self):
            if self._fq is not None:
                if self.type == "bam":
                    return self.bam_basename()
                else:
                    return self.fastq_basename()
            else:
                return self.basename

    return CustomFastqAttrs