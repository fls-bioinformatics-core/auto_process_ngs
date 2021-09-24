#!/usr/bin/env python
#
# sequence lengths library
import os
import json
import logging
from bcftbx.FASTQFile import FastqIterator

# Module specific logger
logger = logging.getLogger(__name__)

class SeqLens(object):
    """
    Wrapper class for handling sequence length data

    The ``SeqLens`` class wraps data handling for sequence
    length data output from the 'get_sequence_lengths'
    function.
    """
    def __init__(self,json_file):
        """
        Initialise new SeqLens class

        Arguments:
          json_file (str): path to JSON file with
            sequence length data output from the
            ``get_sequence_lengths`` function
        """
        self._json_file = os.path.abspath(json_file)
        self._data = dict()
        try:
            with open(self._json_file,'rt') as fp:
                self._data = json.load(fp)
        except Exception as ex:
            logger.warn("Failed to load data from '%s': %s"
                        % (self._json_file,ex))

    @property
    def data(self):
        """
        Return the raw data dictionary
        """
        return self._data

    @property
    def fastq(self):
        """
        Return the path to the associated Fastq file
        """
        return self.data['fastq']

    @property
    def nreads(self):
        """
        Return the number of reads in the Fastq
        """
        return self.data['nreads']

    @property
    def min_length(self):
        """
        Return the minimum sequence length
        """
        return self.data['min_length']

    @property
    def max_length(self):
        """
        Return the maximum sequence length
        """
        return self.data['max_length']

    @property
    def mean(self):
        """
        Return the mean sequence length
        """
        return self.data['mean_length']

    @property
    def range(self):
        """
        Return the range of sequence lengths as a string
        """
        if self.min_length != self.max_length:
            return "%s-%s" % (self.min_length,self.max_length)
        else:
            return "%s" % self.max_length

    @property
    def nmasked(self):
        """
        Return the number of masked reads in the Fastq
        """
        return self.data['nreads_masked']

    @property
    def frac_masked(self):
        """
        Return the fraction of masked reads
        """
        return self.data['frac_reads_masked']

    @property
    def npadded(self):
        """
        Return the number of padded reads in the Fastq
        """
        return self.data['nreads_padded']

    @property
    def frac_padded(self):
        """
        Return the fraction of padded reads
        """
        return self.data['frac_reads_padded']

    @property
    def dist(self):
        """
        Return the sequence length distribution

        This is a dictionary where the keys are
        sequence lengths and the values are the
        corresponding number of sequences
        """
        dist = self.data['seq_lengths_dist']
        return { int(i):dist[i] for i in dist }

    @property
    def masked_dist(self):
        """
        Return the masked sequence length distribution

        This is a dictionary where the keys are
        sequence lengths and the values are the
        corresponding number of masked sequences
        """
        dist = self.data['seq_lengths_masked_dist']
        return { int(i):dist[i] for i in dist }

    def __nonzero__(self):
        return bool(self._data)

def get_sequence_lengths(fastq,outfile=None,show_progress=False,
                         limit=None):
    """
    Get sequence lengths and masking statistics for Fastq

    Returns a dictionary with the following keys:

    - fastq: the Fastq file that metrics were calculated
      from
    - nreads: total number of reads processed
    - nreads_masked: number of reads that are completely
      masked (i.e. consist only of 'N's)
    - nreads_padded: number of partially masked reads
      (i.e. contain trailing 'N's)
    - frac_reads_masked: fraction of the processed reads
      which are masked
    - frac_reads_padded: fraction of the processed reads
      which are padded
    - min_length: minimum read length
    - max_length: maximum read length
    - mean_length: mean read length
    - median_length: median read length
    - seq_lengths_dist: distribution of lengths for all
      reads
    - seq_lengths_masked_dist: distribution of lengths
      for masked reads
    - seq_lengths_padded_dist: distribution of lengths
      for padded reads

    The distributions are each themselves dictionaries
    where the keys are read lengths and the values are
    the number of reads with the matching length; note
    that only lengths with a non-zero number of reads
    are included (zeroes are implied for all other
    lengths).

    Arguments:
      fastq (str): path to Fastq file
      outfile (str): optional, path to output JSON file
      show_progress (bool): if True then print message
        to stdout every 100000 reads indicating progress
        (default: operate silently)
      limit (int): if set then only process this number
        of reads from the head of the Fastq and return
        stats based on these (default: process all reads
        in the file)

    Returns:
      Dictionary: containing the metrics for the Fastq.
    """
    nreads = 0
    nreads_masked = 0
    nreads_padded = 0
    sequence_length = dict()
    reads_all_n = dict()
    reads_padded = dict()
    if show_progress:
        print("\n%s" % fastq)
    for r in FastqIterator(fastq):
        nreads += 1
        if show_progress and nreads%100000 == 0:
            print("...%d reads" % nreads)
        if limit and nreads == limit:
            logging.info("Stopping at limit: %d" % limit)
            break
        seqlen = r.seqlen
        try:
            sequence_length[seqlen] += 1
        except KeyError:
            sequence_length[seqlen] = 1
        nns = 0
        seqlen_no_ns = len(r.sequence.strip('N'))
        if seqlen_no_ns == 0:
            nreads_masked += 1
            try:
                reads_all_n[seqlen] += 1
            except KeyError:
                reads_all_n[seqlen] = 1
        elif seqlen_no_ns < seqlen:
            nreads_padded += 1
            try:
                reads_padded[seqlen] += 1
            except KeyError:
                reads_padded[seqlen] = 1
    # Get statistics
    seqlens = sorted(sequence_length.keys())
    min_len = min(seqlens)
    max_len = max(seqlens)
    mean_len = float(sum([l*sequence_length[l] for l in seqlens]))\
               /float(nreads)
    median_read = nreads//2
    median_len = None
    read_count = 0
    for l in seqlens:
        read_count += sequence_length[l]
        if read_count >= median_read:
            median_len = l
            break
    frac_reads_masked = float(nreads_masked)/nreads*100.0
    frac_reads_padded = float(nreads_padded)/nreads*100.0
    # Build dictionary for output
    stats = dict(fastq=fastq,
                 nreads=nreads,
                 nreads_masked=nreads_masked,
                 nreads_padded=nreads_padded,
                 frac_reads_masked=frac_reads_masked,
                 frac_reads_padded=frac_reads_padded,
                 min_length=min_len,
                 max_length=max_len,
                 mean_length=mean_len,
                 median_length=median_len,
                 seq_lengths_dist=sequence_length,
                 seq_lengths_masked_dist=reads_all_n,
                 seq_lengths_padded_dist=reads_padded)
    # Output to file
    if outfile:
        with open(outfile,'wt') as fp:
            json.dump(stats,fp,sort_keys=True,indent=4)
    # Return stats
    return stats
