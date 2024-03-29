#!/usr/bin/env python
#
# fastq strand library
import os
import logging
from bcftbx.TabFile import TabFile
from bcftbx.utils import AttributeDictionary

# Module specific logger
logger = logging.getLogger(__name__)

"""
Example fastq_strand file for 0.0.1:

#fastq_strand version: 0.0.1	#Aligner: STAR	#Reads in subset: 3
#Genome	1st forward	2nd reverse
Genome1	13.13	93.21
Genome2	13.13	93.21

Example with count data:

"#fastq_strand version: 0.0.1	#Aligner: STAR	#Reads in subset: 3
#Genome	1st forward	2nd reverse	Unstranded	1st read strand aligned	2nd read strand aligned
Genome1	13.13	93.21	391087	51339	364535
"""

class Fastqstrand:
    """
    Class representing data from a fastq_strand.py run

    """
    def __init__(self,fastq_strand_out):
        """
        Create a new Fastqstrand instance

        """
        self._fastq_strand_out = os.path.abspath(fastq_strand_out)
        self._version = None
        self._genomes = AttributeDictionary()
        # Read in data
        tabfile = None
        with open(self._fastq_strand_out,'r') as fp:
            for line in fp:
                line = line.strip()
                if line.startswith('#fastq_strand version:'):
                    self._version = line.split()[2]
                    continue
                elif line.startswith('#Genome'):
                    tabfile = TabFile(column_names=line[1:].split('\t'))
                    continue
                tabfile.append(tabdata=line)
        # Check there is some data
        if tabfile is None:
            raise Exception("Unable to extract fastq_strand data from %s" %
                            self._fastq_strand_out)
        # Copy data to main object
        for line in tabfile:
            # Store the data
            data = AttributeDictionary()
            self._genomes[line['Genome']] = data
            data['forward'] = line['1st forward']
            data['reverse'] = line['2nd reverse']
            # Additional processing
            if data.reverse > 0.0:
                ratio = float(data.forward)/float(data.reverse)
            elif data.forward > 0.0:
                ratio = float("+inf")
            else:
                ratio = None
            if ratio is not None:
                if ratio < 0.2:
                    strandedness = "reverse"
                elif ratio > 5 or ratio == float("+inf"):
                    strandedness = "forward"
                else:
                    strandedness = "unstranded?"
            else:
                strandedness = "undetermined"
            data['ratio'] = ratio
            data['strandedness'] = strandedness

    @property
    def txt(self):
        """
        Path of the fastq_strand output txt file
        """
        return self._fastq_strand_out

    @property
    def version(self):
        """
        Version of fastq_strand which produced the stats
        """
        return self._version

    @property
    def genomes(self):
        """
        List of genome names with strand stats
        """
        return sorted([g for g in self._genomes])

    @property
    def stats(self):
        """
        Data associated with the genomes
        """
        return self._genomes

def build_fastq_strand_conf(organisms,indexes):
    """
    Construct a conf file for input into fastq_strand

    Arguments:
      organisms (list): list of organism names to
        include
      indexes (dict): mapping of organism names to
        paths to STAR indexes

    Returns:
      String: contents of a fastq_strand conf file,
        or None if there were no matching indexes.
    """
    if not organisms:
        return None
    conf = list()
    for organism in organisms:
        try:
            conf.append('\t'.join([organism,indexes[organism]]))
        except KeyError:
            pass
    if conf:
        return '\n'.join(conf)
    else:
        return None
