#     mockqc.py: module providing mock Illumina QC data for testing
#     Copyright (C) University of Manchester 2016-2022 Peter Briggs
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
from .tenx_genomics_utils import CellrangerMultiConfigCsv
from . import mockqcdata
from . import mock10xdata

#######################################################################
# Class definitions
#######################################################################

class MockQCOutputs:
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
    def fastq_screen_v0_9_2(self,fastq,qc_dir,screen_name=None,
                            legacy=False):
        """
        Create mock outputs from Fastq_screen v0.9.2
        """
        # Basename for FastQC products
        basename = self.fastq_basename(fastq)
        if screen_name is not None:
            screen = "_%s" % screen_name
        else:
            screen = ''
        if not legacy:
            screen_basename = "%s_screen%s." % (basename,screen)
        else:
            screen_basename = "%s%s_screen." % (basename,screen)
        screen_basename = os.path.join(qc_dir,screen_basename)
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

    @classmethod
    def seqlens(self,fastq,qc_dir):
        """
        Create mock outputs from sequence length pipeline task
        """
        # Basename for sequence length outputs
        basename = self.fastq_basename(fastq)
        seq_lens_json = os.path.join(qc_dir,
                                     "%s_seqlens.json" %
                                     basename)
        with open(seq_lens_json,'wt') as fp:
            fp.write(mockqcdata.SEQ_LENS_JSON % { 'fastq': fastq })

    @classmethod
    def multiqc_v1_8(self,qc_dir,multiqc_name=None):
        """
        Create mock output from MultiQC v1.8
        """
        if multiqc_name is None:
            multiqc_name = "multi%s_report" % os.path.basename(qc_dir)
        multiqc_html = os.path.join(os.path.dirname(qc_dir),
                                    "%s.html" % multiqc_name)
        with open(multiqc_html,'wt') as fp:
            fp.write("   <a href=\"http://multiqc.info\" "
                     "target=\"_blank\">MultiQC v1.8</a>\n")

    @classmethod
    def cellranger_count(self,sample,qc_dir,cellranger='cellranger',
                         version=None,reference_data_path=
                         "/data/refdata-cellranger-1.2.0",
                         prefix="cellranger_count"):
        """
        Create mock outputs for 'cellranger[-atac|-arc] count'

        Arguments:
          sample (str): sample name to create mock outputs
            for
          qc_dir (str): path to top level QC directory
          cellranger (str): pipeline name; one of 'cellranger',
            'cellranger-atac' or 'cellranger-arc' (default:
            'cellranger')
          version (str): explicit version of 10x pipeline to
            associate with mock outputs (default: determined
            from pipeline name)
          reference_data_path (str): explicit path to
            reference dataset (doesn't have to exist)
          prefix (str): relative path to QC directory to put
            mock outputs into (default: 'cellranger_count')
        """
        # Set internals based on 10x pipeline
        if cellranger == 'cellranger':
            if not version:
                version = '6.1.2'
            cmdline = "cellranger --transcriptome %s" \
                      % reference_data_path
            metrics_data = mock10xdata.METRICS_SUMMARY
            cellranger_output_files = ("web_summary.html",
                                       "metrics_summary.csv")
        elif cellranger == 'cellranger-atac':
            if not version:
                version = '2.0.0'
            cmdline = "cellranger-atac --reference %s" \
                      % reference_data_path
            metrics_data = mock10xdata.ATAC_SUMMARY_2_0_0
            cellranger_output_files = ("web_summary.html",
                                       "summary.csv")
        elif cellranger == 'cellranger-arc':
            if not version:
                version = '2.0.0'
            cmdline = "cellranger-arc --reference %s" \
                      % reference_data_path
            metrics_data = mock10xdata.MULTIOME_SUMMARY_2_0_0
            cellranger_output_files = ("web_summary.html",
                                       "summary.csv")
        # Create and populate the directory for the sample
        sample_dir = os.path.join(qc_dir,
                                  prefix,
                                  sample)
        os.makedirs(sample_dir)
        os.makedirs(os.path.join(sample_dir,"outs"))
        for f in cellranger_output_files:
            with open(os.path.join(sample_dir,"outs",f),'wt') as fp:
                fp.write(metrics_data)
        for f in ("_cmdline",):
            with open(os.path.join(sample_dir,f),'wt') as fp:
                fp.write(cmdline)

    @classmethod
    def cellranger_multi(self,samples,qc_dir,config_csv=None,
                         prefix='cellranger_multi'):
        """
        Create mock outputs for 'cellranger multi'

        Arguments:
          samples (list): sample names to create mock outputs
            for
          qc_dir (str): path to top level QC directory
          config_csv (str): path to an associated CSV config
            file
          prefix (str): relative path to QC directory to put
            mock outputs into (default: 'cellranger_multi')
        """
        # Read in multiplexing config
        if config_csv:
            config = CellrangerMultiConfigCsv(config_csv)
            reference_data_path = config.reference_data_path
            cmdline = "cellranger --csv %s" % config_csv
        else:
            cmdline = "cellranger"
        # Per sample outputs
        per_sample_output_files = ("web_summary.html",
                                   "metrics_summary.csv")
        for sample in samples:
            sample_dir = os.path.join(qc_dir,
                                      prefix,
                                      "outs",
                                      "per_sample_outs",
                                      sample)
            os.makedirs(sample_dir)
            for f in per_sample_output_files:
                with open(os.path.join(sample_dir,f),'wt') as fp:
                    fp.write(mock10xdata.CELLPLEX_METRICS_SUMMARY)
        # Multiplexing analysis outputs
        multiplexing_output_files = ("assignment_confidence_table.csv",
                                     "cells_per_tag.json",
                                     "tag_calls_per_cell.csv",
                                     "tag_calls_summary.csv")
        multiplexing_dir = os.path.join(qc_dir,
                                        prefix,
                                        "outs",
                                        "multi",
                                        "multiplexing_analysis")
        os.makedirs(multiplexing_dir)
        for f in multiplexing_output_files:
            with open(os.path.join(multiplexing_dir,f),'wt') as fp:
                fp.write("")
        # Top-level outputs
        for f in ("_cmdline",):
            with open(os.path.join(qc_dir,prefix,f),'wt') as fp:
                fp.write(cmdline)
