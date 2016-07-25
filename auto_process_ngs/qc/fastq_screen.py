#!/usr/bin/env python
#
# fastq screens library
import os
from bcftbx.TabFile import TabFile

"""
Example screen file for v0.4.1:

#Fastq_screen version: 0.4.1
Library	%Unmapped	%One_hit_one_library	%Multiple_hits_one_library	%One_hit_multiple_libraries	%Multiple_hits_multiple_libraries
hg19	98.10	0.02	0.27	0.55	1.06
mm9	35.92	47.46	10.18	3.56	2.88
rn4	93.01	0.18	0.17	3.87	2.77
dm3	99.97	0.00	0.00	0.01	0.02
ws200	99.98	0.00	0.00	0.00	0.02
ecoli	96.52	0.43	3.05	0.00	0.00
saccer	99.97	0.00	0.00	0.00	0.03
PhiX	99.33	0.67	0.00	0.00	0.00
Vectors	99.99	0.00	0.00	0.01	0.00
SpR6	100.00	0.00	0.00	0.00	0.00

%Hit_no_libraries: 30.80

Example screen file for v0.4.2:

#Fastq_screen version: 0.4.2	#Reads in subset: 1000000
Library	#Reads_processed	#Unmapped	%Unmapped	#One_hit_one_library	%One_hit_one_library	#Multiple_hits_one_library	%Multiple_hits_one_library	#One_hit_multiple_libraries	%One_hit_multiple_libraries	Multiple_hits_multiple_libraries	%Multiple_hits_multiple_libraries
hg19	89393	89213	99.80	1	0.00	0	0.00	11	0.01	168	0.19
mm9	89393	89157	99.74	11	0.01	5	0.01	2	0.00	218	0.24
rn4	89393	89170	99.75	2	0.00	1	0.00	8	0.01	212	0.24
dm3	89393	89391	100.00	0	0.00	0	0.00	0	0.00	2	0.00
ws200	89393	89391	100.00	0	0.00	0	0.00	1	0.00	1	0.00
ecoli	89393	89393	100.00	0	0.00	0	0.00	0	0.00	0	0.00
saccer	89393	89392	100.00	0	0.00	0	0.00	1	0.00	0	0.00
PhiX	89393	89393	100.00	0	0.00	0	0.00	0	0.00	0	0.00
Vectors	89393	89393	100.00	0	0.00	0	0.00	0	0.00	0	0.00
SpR6	89393	89393	100.00	0	0.00	0	0.00	0	0.00	0	0.00

%Hit_no_libraries: 99.73

"""

class Fastqscreen(TabFile):
    """
    Class representing data from a FastqScreen run

    """
    def __init__(self,screen_file):
        """
        Create a new FastqscreenData instance

        """
        TabFile.__init__(self,
                         column_names=('Library',
                                       '%Unmapped',
                                       '%One_hit_one_library',
                                       '%Multiple_hits_one_library',
                                       '%One_hit_multiple_libraries',
                                       '%Multiple_hits_multiple_libraries',))
        self._screen_file = os.path.abspath(screen_file)
        self._version = None
        self._no_hits = None
        # Read in data
        with open(self._screen_file,'r') as fp:
            for line in fp:
                line = line.strip()
                if line.startswith('#Fastq_screen version:'):
                    self._version = line.split()[2]
                    continue
                elif line.startswith('Library') or line.startswith('Genome'):
                    tabfile = TabFile(column_names=line.split())
                    continue
                elif line.startswith('%Hit_no_libraries:') or \
                     line.startswith('%Hit_no_genomes:'):
                    self._no_hits = float(line.split()[-1])
                    continue
                elif not line or \
                   line.startswith('#') or \
                   line.startswith('%'):
                    continue
                tabfile.append(tabdata=line)
        # Handle different terminology for different versions
        if tabfile.header()[0] == 'Library':
            library = 'Library'
            unmapped = '%Unmapped'
            one_hit_one_library = '%One_hit_one_library'
            multiple_hits_one_library = '%Multiple_hits_one_library'
            one_hit_multiple_libraries = '%One_hit_multiple_libraries'
            multiple_hits_multiple_libraries = '%Multiple_hits_multiple_libraries'
        elif tabfile.header()[0] == 'Genome':
            library = 'Genome'
            unmapped = '%Unmapped'
            one_hit_one_library = '%One_hit_one_genome'
            multiple_hits_one_library = '%Multiple_hits_one_genome'
            one_hit_multiple_libraries = '%One_hit_multiple_genomes'
            multiple_hits_multiple_libraries = '%Multiple_hits_multiple_genomes'
        # Copy data to main object
        for line in tabfile:
            data = [line[library],
                    line[unmapped],
                    line[one_hit_one_library],
                    line[multiple_hits_one_library],
                    line[one_hit_multiple_libraries],
                    line[multiple_hits_multiple_libraries]]
            self.append(data=data)

    @property
    def txt(self):
        """
        Path of the fastq_screen.txt file
        """
        return self._screen_file

    @property
    def png(self):
        """
        Path of the fastq_screen.png file
        """
        return os.path.splitext(self._screen_file)[0]+'.png'

    @property
    def version(self):
        """
        Version of fastq_screen which produced the screens
        """
        return self._version

    @property
    def libraries(self):
        """
        List of library names used in the screen
        """
        return [lib['Library'] for lib in self]

    @property
    def no_hits(self):
        """
        Percentage of reads with no hits on any library
        """
        return self._no_hits
