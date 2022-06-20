#     mockqc.py: module providing mock Illumina QC data for testing
#     Copyright (C) University of Manchester 2016-2022 Peter Briggs
#
########################################################################

"""
mockqc.py

Provides classes and functions for mocking up examples of QC artefacts
from the process pipeline, to be used in testing.

The core class is:

- MockQCOutputs

which provides a set of static methods for creating mock outputs from
different internal and third-party QC components (e.g. fastqc,
fastq_screen etc).

In addition there is a factory function for creating mock QC
directories in different configurations:

- make_mock_qc_dir: create and populate a mock QC directory

"""

#######################################################################
# Import modules that this module depends on
#######################################################################

import os
import base64
import bcftbx.utils
from bcftbx.mock import MockIlluminaData
from .analysis import AnalysisFastq
from .fastq_utils import group_fastqs_by_name
from .metadata import AnalysisProjectQCDirInfo
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
    def picard_collect_insert_size_metrics(self,fq,organism,qc_dir):
        """
        Create mock outputs from Picard CollectInsertSizeMetrics
        """
        # Basename for insert size metrics outputs
        basename = AnalysisFastq(fq)
        basename.read_number = None
        basename = str(basename)
        out_dir = os.path.join(qc_dir,
                               "picard",
                               organism)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        for ext in ('.insert_size_metrics.txt',
                    '.insert_size_histogram.pdf'):
            f = os.path.join(out_dir,"%s%s" % (basename,ext))
            with open(f,'wt') as fp:
                if f.endswith('.txt'):
                    fp.write(mockqcdata.PICARD_COLLECT_INSERT_SIZE_METRICS)
                else:
                    fp.write("Placeholder\n")

    @classmethod
    def rseqc_infer_experiment(self,fq,organism,qc_dir):
        """
        Create mock outputs from RSeQC infer_experiment.py
        """
        # Basename for mock infer_experiment.py log file
        basename = AnalysisFastq(fq)
        basename.read_number = None
        basename = str(basename)
        out_dir = os.path.join(qc_dir,
                               "rseqc_infer_experiment",
                               organism)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        f = os.path.join(out_dir,
                         "%s.infer_experiment.log" % basename)
        with open(f,'wt') as fp:
            fp.write("""This is PairEnd Data
Fraction of reads failed to determine: 0.0172
Fraction of reads explained by "1++,1--,2+-,2-+": 0.4903
Fraction of reads explained by "1+-,1-+,2++,2--": 0.4925
""")

    @classmethod
    def rseqc_genebody_coverage(self,name,organism,qc_dir):
        """
        Create mock outputs from RSeQC geneBody_coverage.py
        """
        out_dir = os.path.join(qc_dir,
                               "rseqc_genebody_coverage",
                               organism)
        os.makedirs(out_dir)
        for ext in ('.geneBodyCoverage.curves.png',
                    '.geneBodyCoverage.r',
                    '.geneBodyCoverage.txt'):
            f = os.path.join(out_dir,"%s%s" % (name,ext))
            with open(f,'wt') as fp:
                if f.endswith('.txt'):
                    fp.write(mockqcdata.RSEQC_GENEBODY_COVERAGE_TXT)
                else:
                    fp.write("Placeholder\n")

    @classmethod
    def qualimap_rnaseq(self,fq,organism,qc_dir):
        """
        Create mock outputs from Qualimap 'rnaseq' function
        """
        # Basename for Qualimap rnaseq outputs
        basename = AnalysisFastq(fq)
        basename.read_number = None
        basename = str(basename)
        out_dir = os.path.join(qc_dir,
                               "qualimap-rnaseq",
                               organism,
                               basename)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        for f in ('qualimapReport.html',
                  'rnaseq_qc_results.txt'):
            ff = os.path.join(out_dir,f)
            with open(ff,'wt') as fp:
                if f == 'rnaseq_qc_results.txt':
                    fp.write(mockqcdata.QUALIMAP_RNASEQ_RESULTS)
                else:
                    fp.write("Placeholder\n")

    @classmethod
    def multiqc(self,dirn,multiqc_html=None,version="1.8"):
        """
        Create mock output from MultiQC
        """
        if multiqc_html is None:
            multiqc_html = "multiqc_report.html"
        multiqc_html = os.path.join(dirn,multiqc_html)
        with open(multiqc_html,'wt') as fp:
            fp.write("   <a href=\"http://multiqc.info\" "
                     "target=\"_blank\">MultiQC v%s</a>\n" %
                     version)

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

#######################################################################
# Functions
#######################################################################

def make_mock_qc_dir(qc_dir,fastq_names,fastq_dir=None,
                     protocol=None,
                     project_name=None,
                     screens=('model_organisms',
                              'other_organisms',
                              'rRNA',),
                     organisms=('human',),
                     cellranger_pipelines=('cellranger',),
                     cellranger_samples=None,
                     cellranger_multi_samples=None,
                     include_fastqc=True,
                     include_fastq_screen=True,
                     include_strandedness=True,
                     include_seqlens=True,
                     include_rseqc_infer_experiment=False,
                     include_rseqc_genebody_coverage=False,
                     include_picard_insert_size_metrics=False,
                     include_qualimap_rnaseq=False,
                     include_multiqc=False,
                     include_cellranger_count=False,
                     include_cellranger_multi=False,
                     legacy_screens=False,
                     legacy_cellranger_outs=False):
    """
    Create a mock QC directory with QC artefacts

    Arguments:
      qc_dir (str): path to the mock QC directory to create
      fastq_names (list): Fastq names to make outputs for
      fastq_dir (str): optional, set a non-standard
        directory for the Fastq files
      protocol (str): QC protocol to emulate
      project_name (str): optional, specify the project name
      screens (list): optional, list of non-standard
        FastqScreen panel names
      organisms (list): optional, list of organism names for
        extended QC metrics
      cellranger_pipelines (list): list of 10xGenomics pipelines
        to make mock outputs for (e.g. 'cellranger',
        'cellranger-atac' etc)
      cellranger_samples (list): list of sample names to
        produce 'cellranger count' outputs for
      include_fastqc (bool): include outputs from Fastqc
      include_fastq_screen (bool): include outputs from
        FastqScreen
      include_strandedness (bool): include outputs from
        strandedness
      include_seqlens (bool): include sequence length metrics
      include_rseqc_infer_experiment (bool): include RSeQC
        infer_experiment.py outputs
      include_rseqc_genebody_coverage (bool): include RSeQC
        geneBody_coverage.py outputs
      include_picard_insert_size_metrics (bool): include Picard
        CollectInsertSizeMetrics outputs
      include_qualimap_rnaseq (bool): include Qualimap 'rnaseq'
        outputs
      include_multiqc (bool): include MultiQC outputs
      include_celllranger_count (bool): include 'cellranger
        count' outputs
      include_cellranger_multi (bool): include 'cellranger
        multi' outputs
      cellranger_multi_samples (list): list of sample names to
        produce 'cellranger multi' outputs for
      legacy_screens (bool): if True then use legacy naming
        convention for FastqScreen outputs
      legacy_cellranger_outs (bool): if True then use legacy
        naming convention for 10xGenomics pipeline outputs

    Returns:
      String: path to the mock QC directory that was created.
    """
    # Make an empty QC dir
    if not os.path.exists(qc_dir):
        os.mkdir(qc_dir)
    # Project name
    if project_name is None:
        project_name = os.path.basename(os.path.dirname(qc_dir))
    # QC metadata
    qc_info = AnalysisProjectQCDirInfo()
    qc_info['protocol'] = protocol
    if organisms:
        qc_info['organism'] = ','.join(organisms)
    # Populate with fake QC products
    for fq in fastq_names:
        # FastQC
        if include_fastqc:
            MockQCOutputs.fastqc_v0_11_2(fq,qc_dir)
        # Fastq_screen
        if include_fastq_screen:
            for screen in screens:
                MockQCOutputs.fastq_screen_v0_9_2(
                    fq,qc_dir,screen,legacy=legacy_screens)
            qc_info['fastq_screens'] = ','.join(screens)
        # Sequence lengths
        if include_seqlens:
            MockQCOutputs.seqlens(fq,qc_dir)
    for fq_group in group_fastqs_by_name(fastq_names):
        if protocol in ('10x_scRNAseq',
                        '10x_snRNAseq',
                        '10x_Multiome_GEX',
                        '10x_CellPlex',
                        '10x_Visium',):
            fq = fq_group[1]
        else:
            fq = fq_group[0]
        # Strandedness
        if include_strandedness:
            MockQCOutputs.fastq_strand_v0_0_4(fq,qc_dir)
        # RSeQC infer_experiment.py
        if include_rseqc_infer_experiment:
            for organism in organisms:
                MockQCOutputs.rseqc_infer_experiment(
                    fq,organism,qc_dir)
        # Picard insert size metrics
        if include_picard_insert_size_metrics:
            for organism in organisms:
                MockQCOutputs.picard_collect_insert_size_metrics(
                    fq,organism,qc_dir)
        # Qualimap rnaseq
        if include_qualimap_rnaseq:
            for organism in organisms:
                MockQCOutputs.qualimap_rnaseq(fq,organism,qc_dir)
    # Version file for RSeQC infer_experiment.py
    if include_rseqc_infer_experiment:
        for organism in organisms:
            with open(os.path.join(qc_dir,
                                   "rseqc_infer_experiment",
                                   organism,
                                   "__versions"),'wt') as fp:
                fp.write("rseqc:infer_experiment\t4.0.0\n")
    # Extra files for insert sizes
    if include_picard_insert_size_metrics:
        # Collated insert sizes
        for organism in organisms:
            with open(os.path.join(
                    qc_dir,
                    "insert_sizes.%s.tsv" % organisms),'wt') as fp:
                fp.write("Placeholder\n")
        # Picard version
        for organism in organisms:
            with open(os.path.join(qc_dir,
                                   "picard",
                                   organism,
                                   "__versions"),'wt') as fp:
                fp.write("picard\t2.27.1\n")
    # Version file for Qualimap rnaseq
    if include_qualimap_rnaseq:
        for organism in organisms:
            with open(os.path.join(qc_dir,
                                   "qualimap-rnaseq",
                                   organism,
                                   "__versions"),'wt') as fp:
                fp.write("qualimap\tv.2.2.2\n")
    # Strandedness conf file
    if include_strandedness:
        with open(os.path.join(qc_dir,
                               "fastq_strand.conf"),'wt') as fp:
            fp.write("Placeholder\n")
    # RSeQC gene body coverage
    if include_rseqc_genebody_coverage:
        for organism in organisms:
            MockQCOutputs.rseqc_genebody_coverage(project_name,
                                                  organism,
                                                  qc_dir)
            with open(os.path.join(qc_dir,
                                   "rseqc_genebody_coverage",
                                   organism,
                                   "__versions"),'wt') as fp:
                fp.write("rseqc:genebody_coverage\t4.0.0\n")
    # MultiQC
    if include_multiqc:
        out_file = "multi%s_report.html" % os.path.basename(qc_dir)
        MockQCOutputs.multiqc(os.path.dirname(qc_dir),
                              multiqc_html=out_file,
                              version="1.8")
    # Cellranger count
    if include_cellranger_count:
        for cellranger in cellranger_pipelines:
            # Set defaults
            if cellranger == "cellranger":
                version = "6.1.2"
                refdata = "/data/refdata-cellranger-2020-A"
            elif cellranger == "cellranger-atac":
                version = "2.0.0"
                refdata = "/data/refdata-cellranger-atac-2020-A"
            elif cellranger == "cellranger-arc":
                version = "2.0.0"
                refdata = "/data/refdata-cellranger-arc-2020-A"
            # Set top-level output dir
            if not legacy_cellranger_outs:
                count_dir = os.path.join("cellranger_count",
                                         version,
                                         os.path.basename(refdata))
            else:
                count_dir = "cellranger_count"
            # Make pipeline outputs
            for sample in cellranger_samples:
                MockQCOutputs.cellranger_count(
                    sample,
                    qc_dir,
                    cellranger=cellranger,
                    version=version,
                    reference_data_path=refdata,
                    prefix=count_dir)
                if cellranger == "cellranger-arc":
                    multiome_config = os.path.join(qc_dir,
                                                   "libraries.%s.csv" %
                                                   sample)
                    with open(multiome_config,'wt') as fp:
                        fp.write("Placeholder\n")
            if not legacy_cellranger_outs:
                qc_info['cellranger_version'] = version
            qc_info['cellranger_refdata'] = refdata
    # Cellranger multi
    if include_cellranger_multi:
        # Make cellranger multi config.csv file
        multi_config = os.path.join(qc_dir,"10x_multi_config.csv")
        with open(multi_config,'wt') as fp:
            if not fastq_dir:
                fastq_dir = os.path.join(os.path.dirname(qc_dir),
                                         "fastqs")
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-2020-A

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB1_GEX,%s,any,PJB1,gene expression,
PJB2_MC,%s,any,PJB2,Multiplexing Capture,

[samples]
sample_id,cmo_ids,description
PJB_CML1,CMO301,CML1
PJB_CML2,CMO302,CML2
""" % (fastq_dir,fastq_dir))
        # Cellranger version
        version = "6.1.2"
        # Set top-level output dir
        multi_dir = os.path.join("cellranger_multi",
                                 version,
                                 "refdata-cellranger-2020-A")
        # Make outputs
        MockQCOutputs.cellranger_multi(cellranger_multi_samples,
                                       qc_dir,
                                       config_csv=multi_config,
                                       prefix=multi_dir)
        qc_info['cellranger_version'] = version
    # Additional metadata items
    star_index = "/data/star/hg38"
    annotation_bed = "/data/annotation/hg38.bed"
    annotation_gtf = "/data/annotation/hg38.gtf"
    if include_picard_insert_size_metrics:
        qc_info['star_index'] = star_index
    if include_rseqc_genebody_coverage or \
       include_rseqc_infer_experiment:
        qc_info['star_index'] = star_index
        qc_info['annotation_bed'] = annotation_bed
    if include_qualimap_rnaseq:
        qc_info['star_index'] = star_index
        qc_info['annotation_bed'] = annotation_bed
        qc_info['annotation_gtf'] = annotation_gtf
    # Write out metadata file
    qc_info.save(os.path.join(qc_dir,"qc.info"))
    return qc_dir
