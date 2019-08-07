#!/usr/bin/env python
#
#     barcodes/splitter.py: classes & functions to sort reads by barcode
#     Copyright (C) University of Manchester 2019 Peter Briggs
#
#########################################################################
#
# barcodes/splitter.py
#
#########################################################################

"""
Classes and functions for sorting reads from FASTQ files by barcode (i.e.
index sequences) from the read headers:

- HammingMetrics: class implementing Hamming distance calculators
- HammingLookup: class for fetching Hamming distance for two strings
- get_fastqs_from_dir: list Fastqs matching specific lane
- split_single_end: split reads from single ended data
- split_paired_end: split reads from pair ended data

"""

import itertools
import logging
import bcftbx.IlluminaData as IlluminaData
import bcftbx.FASTQFile as FASTQFile
from ..utils import OutputFiles

#########################################################################
# Classes
#########################################################################

class HammingMetrics(object):
    """
    Calculate Hamming distances between two strings

    Implements a set of functions (as class methods) fo
    calculating Hamming distances between two strings
    (sequences):

    - hamming_distance: for equal-length sequences
    - hamming_distance_truncate: for non-equal-length
      sequences
    - hamming_distance_with_N: for equal-length sequences
      where N's automatically mismatch

    For background on Hamming distance see e.g.:
    
    http://en.wikipedia.org/wiki/Hamming_distance
    
    """
    @classmethod
    def hamming_distance(self,s1,s2):
        """
        Get Hamming distance for equal-length sequences

        """
        if len(s1) != len(s2):
            raise ValueError("Undefined for sequences of unequal length")
        return sum(ch1 != ch2 for ch1, ch2 in itertools.izip(s1, s2))

    @classmethod
    def hamming_distance_truncate(self,s1,s2):
        """
        Get Hamming distance for non-equal-length sequences

        """
        if len(s1) != len(s2):
            l = min(len(s1),len(s2))
            s1 = s1[:l]
            s2 = s2[:l]
        return sum(ch1 != ch2 for ch1, ch2 in itertools.izip(s1, s2))

    @classmethod
    def hamming_distance_with_N(self,s1,s2):
        """
        Get Hamming distance for equal-length sequences (no match for Ns)

        """
        # modified version where N's don't match anything (not even other N's)
        if len(s1) != len(s2):
            raise ValueError("Undefined for sequences of unequal length")
        return sum((ch1 != ch2 or ch1 == 'N' or ch2 == 'N') \
                   for ch1, ch2 in itertools.izip(s1, s2))

class HammingLookup(object):
    """
    Class for handling Hamming distances for multiple sequences

    
    
    """
    def __init__(self,hamming_func=HammingMetrics.hamming_distance):
        """
        Create a new HammingLookup instance

        Arguments:
          hamming_func (Function): the function to use
            for calculating Hamming distances
        """
        self._hamming = hamming_func
        self._distances = dict()

    def dist(self,s1,s2):
        """
        Return the Hamming distance for two sequences

        If the sequences are new then calculates,
        caches and returns the Hamming distance for this
        pair (calculated using the function supplied on
        instantiation); otherwise returns the cached distance
        calculated previously.

        Arguments:
          s1 (str): first sequence
          s2 (str): second sequence

        Returns:
          Integer: Hamming distance for the two sequences.
        """
        try:
            return self._distances[s1][s2]
        except KeyError:
            try:
                return self._distances[s2][s1]
            except KeyError:
                pass
        if s1 not in self._distances:
            self._distances[s1] = dict()
        self._distances[s1][s2] = self._hamming(s1,s2)
        return self._distances[s1][s2]

class BarcodeMatcher(object):
    """
    Class to match barcode sequences against reference set
    """
    def __init__(self,index_seqs,max_dist=0):
        """
        Create a new BarcodeMatcher instance

        Argument:
          index_seqs (list): list of reference sequences
            to match against
          max_dist (int): optional, sequence pair must be
            at or below this distance to be counted as a
            match (defaults to 0 i.e. exact matches only)
        """
        self._hamming = HammingLookup(
            hamming_func=HammingMetrics.hamming_distance_truncate)
        for i,seq1 in enumerate(index_seqs):
            for seq2 in index_seqs[i+1:]:
                if self._hamming.dist(seq1,seq2) <= max_dist:
                    raise Exception("Index sequence ambiguity: '%s' "
                                    "and '%s' are too similar (differ "
                                    "by %d bases, must be > %d)" %
                                    (seq1,
                                     seq2,
                                     self._hamming.dist(seq1,seq2),
                                     max_dist))
        self._index_seqs = sorted(list(index_seqs))
        self._max_dist = max_dist

    @property
    def sequences(self):
        """
        Return list of stored index sequences
        """
        return self._index_seqs

    def match(self,seq):
        """
        Find matching index sequence for supplied sequence

        Arguments:
          seq (str): sequence to match against indexes

        Returns:
          String: matching index sequence, or None if no match
            was found.
        """
        for index_seq in self._index_seqs:
            logging.debug("%s, %s" % (index_seq,seq))
            if self._hamming.dist(index_seq,seq) <= self._max_dist:
                return index_seq
        return None

#########################################################################
# Functions
#########################################################################

def get_fastqs_from_dir(dirn,lane,unaligned_dir=None):
    """
    Collect Fastq files for specified lane

    Arguments:
      dirn (str): path to directory to collect Fastq
        files from
      lane (int): lane Fastqs must have come from
      unaligned_dir (str): subdirectory of 'dirn' with
        outputs from bcl2fastq

    Returns:
      List: list of Fastqs (for single ended data) or of
        Fastq pairs (for pair ended data).
    """
    try:
        illumina_data = IlluminaData.IlluminaData(dirn,
                                                  unaligned_dir=unaligned_dir)
    except Exception as ex:
        raise Exception("Unable to read fastqs from %s: %s\n" % (dirn,ex))
    paired_end = illumina_data.paired_end
    fastqs_r1 = []
    fastqs_r2 = []
    for project in illumina_data.projects:
        for sample in project.samples:
            for fastq in sample.fastq_subset(read_number=1,full_path=True):
                if IlluminaData.IlluminaFastq(fastq).lane_number == lane:
                    fastqs_r1.append(fastq)
            for fastq in sample.fastq_subset(read_number=2,full_path=True):
                if IlluminaData.IlluminaFastq(fastq).lane_number == lane:
                    fastqs_r2.append(fastq)
    if illumina_data.undetermined:
        for sample in illumina_data.undetermined.samples:
            for fastq in sample.fastq_subset(read_number=1,full_path=True):
                if IlluminaData.IlluminaFastq(fastq).lane_number == lane:
                    fastqs_r1.append(fastq)
            for fastq in sample.fastq_subset(read_number=2,full_path=True):
                if IlluminaData.IlluminaFastq(fastq).lane_number == lane:
                    fastqs_r2.append(fastq)
    if not paired_end:
        return fastqs_r1
    fastqs = []
    fastqs_r1.sort()
    fastqs_r2.sort()
    for fq1,fq2 in zip(fastqs_r1,fastqs_r2):
        fastqs.append("%s,%s" % (fq1,fq2))
    return fastqs

def split_single_end(matcher,fastqs,base_name=None,output_dir=None):
    """
    Split reads from single ended data

    For each fastq file in 'fastqs', check reads against the index
    sequences in the BarcodeMatcher 'matcher' and write to an
    appropriate file.

    Arguments:
      matcher (BarcodeMatcher): barcoder matcher instance
      fastqs (list): list of Fastqs to split
      base_name (str): optional, base name to use for output
        Fastq files
      output_dir (str): optional, path to directory to write
        output Fastqs to

    """
    if base_name is None:
        base_name = ''
    else:
        base_name = "%s." % base_name
    fp = OutputFiles(base_dir=output_dir)
    for barcode in matcher.sequences:
        fp.open(barcode,"%s%s.fastq" % (base_name,barcode))
    fp.open('undetermined',"%sundetermined.fastq" % base_name)
    # Filter reads
    nread = 0
    for fastq in fastqs:
        print("Processing reads from %s" % fastq)
        for read in FASTQFile.FastqIterator(fastq):
            nread += 1
            seq = read.seqid.index_sequence
            if not seq:
                raise Exception("%s: no index sequence for read %d" %
                                (fastq,nread))
            assigned_index = matcher.match(seq)
            # Read not assigned
            if assigned_index is None:
                assigned_index = 'undetermined'
            logging.debug("Assigned read #%d to %s" % (nread,assigned_index))
            fp.write(assigned_index,read)
    print("Finished (%d reads processed)" % nread)

def split_paired_end(matcher,fastq_pairs,base_name=None,output_dir=None):
    """
    Split reads from paired end data

    For each fastq file pair in 'fastqs', check reads against the
    index sequences in the BarcodeMatcher 'matcher' and write to an
    appropriate file.

    Arguments:
      matcher (BarcodeMatcher): barcoder matcher instance
      fastqs (list): list of Fastq pairs to split
      base_name (str): optional, base name to use for output
        Fastq files
      output_dir (str): optional, path to directory to write
        output Fastqs to

    """
    if base_name is None:
        base_name = ''
    else:
        base_name = "%s." % base_name
    fp = OutputFiles(base_dir=output_dir)
    for barcode in matcher.sequences:
        fp.open((barcode,'R1'),"%s%s_R1.fastq" % (base_name,barcode))
        fp.open((barcode,'R2'),"%s%s_R2.fastq" % (base_name,barcode))
    fp.open(('undetermined','R1'),"%sundetermined_R1.fastq" % base_name)
    fp.open(('undetermined','R2'),"%sundetermined_R2.fastq" % base_name)
    # Filter reads
    nread = 0
    for fq_r1,fq_r2 in fastq_pairs:
        print("Processing reads from fastq pair %s %s" % (fq_r1,fq_r2))
        for read1,read2 in itertools.izip(FASTQFile.FastqIterator(fq_r1),
                                          FASTQFile.FastqIterator(fq_r2)):
            nread += 1
            seq = read1.seqid.index_sequence
            if not seq:
                raise Exception("%s: no index sequence for read %d" %
                                (fq_r1,nread))
            if seq != read2.seqid.index_sequence:
                raise Exception("Index sequence mismatch between R1 and "
                                "R2 reads")
            assigned_index = matcher.match(seq)
            # Read not assigned
            if assigned_index is None:
                assigned_index = 'undetermined'
            logging.debug("Assigned read #%d to %s" % (nread,assigned_index))
            fp.write((assigned_index,'R1'),read1)
            fp.write((assigned_index,'R2'),read2)
    print("Finished (%d read pairs processed)" % nread)

