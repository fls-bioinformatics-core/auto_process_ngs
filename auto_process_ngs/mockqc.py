#     mockqc.py: module providing mock Illumina QC data for testing
#     Copyright (C) University of Manchester 2016-2019 Peter Briggs
#
########################################################################

"""
mockqc.py

Provides classes for mocking up examples of QC outputs for the
process pipeline, to be used in testing.

These include:

- MockQCOutputs

"""

#######################################################################
# Import modules that this module depends on
#######################################################################

import os
import base64
import bcftbx.utils
from bcftbx.mock import MockIlluminaData
from . import mockqcdata

#######################################################################
# Class definitions
#######################################################################

class MockQCOutputs(object):
    """
    Utility class for creating mock auto-process QC outputs
    """
    @classmethod
    def fastq_basename(self,fastq):
        """
        Return the basename for a FASTQ file
        """
        basename = os.path.basename(fastq)
        if basename.endswith('.gz'):
            basename = '.'.join(basename.split('.')[:-1])
        if basename.endswith('.fastq'):
            basename = '.'.join(basename.split('.')[:-1])
        return basename

    @classmethod
    def fastqc_v0_11_2(self,fastq,qc_dir):
        """
        Create mock outputs from FastQC v0.11.2
        """
        # Basename for FastQC products
        basename = self.fastq_basename(fastq)
        # FastQC output directory:
        # - <basename>_fastqc
        fastqc_dir = os.path.join(qc_dir,basename+'_fastqc')
        os.mkdir(fastqc_dir)
        # Make fake HTML and zip files:
        # - <basename>_fastqc.html
        # - <basename>_fastqc.zip
        for ext in ('.html','.zip'):
            with open(fastqc_dir+ext,'w') as fp:
                fp.write('')
        # Make summary.txt file
        with open(os.path.join(fastqc_dir,'summary.txt'),'w') as fp:
            for name,status in (('Basic Statistics','PASS'),
                                ('Per base sequence quality','FAIL'),
                                ('Per tile sequence quality','PASS'),
                                ('Per sequence quality scores','PASS'),
                                ('Per sequence quality scores','PASS'),
                                ('Per base sequence content','FAIL'),
                                ('Per sequence GC content','FAIL'),
                                ('Per base N content','WARN'),
                                ('Sequence Length Distribution','PASS'),
                                ('Sequence Duplication Levels','PASS'),
                                ('Overrepresented sequences','FAIL'),
                                ('Adapter Content','FAIL'),
                                ('Kmer Content','PASS')):
                fp.write("%s\t%s\t%s.fastq\n" % (status,name,basename))
        # Make fastqc_data.txt file
        with open(os.path.join(fastqc_dir,'fastqc_data.txt'),'w') as fp:
            fp.write(mockqcdata.FASTQC_0_11_3['fastqc_data.txt'] %
                     { 'fastq': ('%s.fastq' % basename) })
        # Make fake files:
        for f in ('fastqc.fo',
                  'fastqc_report.html'):
            with open(os.path.join(fastqc_dir,f),'w') as fp:
                fp.write('')
        # Make and populate fake subdirectories
        icons_dir = os.path.join(fastqc_dir,'Icons')
        os.mkdir(icons_dir)
        for f in ('error.png',
                  'fastqc_icon.png',
                  'tick.png',
                  'warning.png'):
            with open(os.path.join(icons_dir,f),'w') as fp:
                fp.write('')
        # Write fake PNG images
        images_dir = os.path.join(fastqc_dir,'Images')
        os.mkdir(images_dir)
        for f in ('adapter_content.png',
                  'duplication_levels.png',
                  'kmer_profiles.png',
                  'per_base_n_content.png',
                  'per_base_quality.png',
                  'per_base_sequence_content.png',
                  'per_sequence_gc_content.png',
                  'per_sequence_quality.png',
                  'per_tile_quality.png',
                  'sequence_length_distribution.png'):
            with open(os.path.join(images_dir,f),'wb') as fp:
                fp.write(base64.b64decode(mockqcdata.BASE64_PNG_DATA))

    @classmethod
    def fastq_screen_v0_9_2(self,fastq,qc_dir,screen_name=None):
        """
        Create mock outputs from Fastq_screen v0.9.2
        """
        # Basename for FastQC products
        basename = self.fastq_basename(fastq)
        if screen_name is not None:
            screen = "_%s" % screen_name
        else:
            screen = ''
        screen_basename = os.path.join(qc_dir,
                                       "%s%s_screen." % (basename,screen))
        # Raw data
        with open(screen_basename+'txt','w') as fp:
            fp.write(mockqcdata.FASTQ_SCREEN_V0_9_2['screen.txt'])
        # Fake PNG
        with open(screen_basename+'png','wb') as fp:
            fp.write(base64.b64decode(mockqcdata.BASE64_PNG_DATA))

    @classmethod
    def fastq_strand_v0_0_4(self,fastq,qc_dir):
        """
        Create mock outputs from fastq_strand.py v0.0.4
        """
        # Basename for Fastq_strand products
        basename = self.fastq_basename(fastq)
        fastq_strand_txt = os.path.join(qc_dir,
                                        "%s_fastq_strand.txt" %
                                        basename)
        with open(fastq_strand_txt,'w') as fp:
            fp.write(mockqcdata.FASTQ_STRAND_V0_0_4['fastq_strand.txt'])
