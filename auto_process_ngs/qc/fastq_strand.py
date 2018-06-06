#!/usr/bin/env python
#
# fastq screens library
import os
from bcftbx.TabFile import TabFile
from bcftbx.utils import AttributeDictionary

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

class Fastqstrand(object):
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
            ratio = float(data.forward)/float(data.reverse)
            if ratio < 0.1:
                strandedness = "reverse"
            elif ratio > 10:
                strandedness = "forward"
            else:
                strandedness = "unstranded"
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
