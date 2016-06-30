#!/usr/bin/env python
#
# Fastq statistics utilities
from bcftbx.FASTQFile import FastqIterator
from .fastqc import FastqcData

class FastqQualityStats:
    """
    Class for storing per-base quality stats from a FASTQ

    This class acts as a container for per-base quality
    statistics from a FASTQ file; the statistics can be
    accessed via various properties:

    mean: mean of all quality scores at each position
    median: mean of all quality scores
    q25: lower quartile
    q75: upper quartile
    p10: 10th percentile
    p90: 90th percentile

    Each of these properties is a Python list, so to
    access e.g. the mean for the 5th base::

    >>> mean = fqstats.mean[4]

    (as the lists are indexed starting from zero).

    The ``nbases`` property returns the number of bases
    for which data is stored.

    The statistics can be populated directly from a FASTQ
    file, for example::

    >>> stats = FastqQualityStats()
    >>> stats.from_fastq('example.fq')

    Alternatively they can be loaded from a
    ``fastqc_data.txt`` file output from the FastQC program:

    >>> stats.from_fastqc_data('example_fastqc/fastqc_data.txt')

    """
    def __init__(self):
        self.mean = []
        self.median = []
        self.q25 = []
        self.q75 = []
        self.p10 = []
        self.p90 = []

    @property
    def nbases(self):
        """
        Return the number of bases for which data is stored

        """
        return len(self.mean)

    def from_fastq(self,fastq):
        """
        Get statistics from a FASTQ file

        Generates and stores statistics from a FASTQ file.

        Arguments:
          fastq (str): path to a FASTQ file (can be gzipped)

        """
        # Initialise using first read for sequence length
        quality_per_base = []
        for read in FastqIterator(fastq):
            for i in xrange(read.seqlen):
                quality_per_base.append({})
                for j in xrange(ord('!'),ord('I')+1):
                    quality_per_base[i][chr(j)] = 0
            break

        # Iterate through fastq file and count quality scores
        nreads = 0
        for read in FastqIterator(fastq):
            nreads += 1
            for pos,q in enumerate(read.quality):
                quality_per_base[pos][q] += 1
        #print quality_per_base

        # Median etc positions
        # FIXME these are not correct if the list has an odd number of values!
        median_pos = nreads/2
        q25_pos = median_pos/2
        q75_pos = median_pos + q25_pos
        # For percentiles see http://stackoverflow.com/a/2753343/579925
        # FIXME 10th/90th percentiles is a fudge here
        p10_pos = nreads/10
        p90_pos = nreads - p10_pos

        # For each base position determine stats
        for pos,counts in enumerate(quality_per_base):
            #print "Position: %d" % pos
            # Expand to a list
            scores = ''
            for q in counts:
                scores += q*counts[q]
            # Sort into order
            scores = ''.join(sorted(scores))
            #print scores
            # Get the mean (scores are Phred+33 encoded)
            self.mean.append(float(sum([(ord(q)-33)
                                        for q in scores]))/nreads)
            #print "Mean: %.2f" % self.mean[pos]
            # Get the median etc
            self.median.append(ord(scores[median_pos])-33)
            self.q25.append(ord(scores[q25_pos])-33)
            self.q75.append(ord(scores[q75_pos])-33)
            self.p10.append(ord(scores[p10_pos])-33)
            self.p90.append(ord(scores[p90_pos])-33)
            #print "Median: %d" % self.median[pos]
            #print "Q25   : %d" % self.q25[pos]
            #print "Q75   : %d" % self.q75[pos]
            #print "P10   : %d" % self.p10[pos]
            #print "P90   : %d" % self.p90[pos]

    def from_fastqc_data(self,fastqc_data):
        """
        Get statistics from a FastQC data file

        Reads in the per-base quality statistics from the
        ``Per base sequence quality`` module of the FastQC
        program, which are stored in the ``fastqc_data.txt``
        output file.

        Arguments:
          fastqc_data (str): path to a FastQC
            ``fastqc_data.txt`` file

        """
        for line in FastqcData(fastqc_data).data('Per base sequence quality'):
            if line.startswith('#'):
                continue
            i,mean,median,q25,q75,p10,p90 = line.strip().split('\t')
            self.mean.append(int(float(mean)))
            self.median.append(int(float(median)))
            self.q25.append(int(float(q25)))
            self.q75.append(int(float(q75)))
            self.p10.append(int(float(p10)))
            self.p90.append(int(float(p90)))
