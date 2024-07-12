#!/usr/bin/env python
#
#     qc/outputs: utilities to predict and check QC pipeline outputs
#     Copyright (C) University of Manchester 2019-2024 Peter Briggs
#
"""
Provides utility classes and functions for QC outputs.

Provides the following classes:

- QCOutputs: detect and characterise QC outputs
- ExtraOutputs: helper class for reading 'extra_outputs.tsv' file

Provides the following functions:

- fastq_strand_output: get name for fastq_strand.py output
- picard_collect_insert_size_metrics_output: get names for Picard
  CollectInsertSizeMetrics output
- rseqc_genebody_coverage_output: get names for RSeQC geneBody_coverage.py
  output
- qualimap_rnaseq_output: get names for Qualimap 'rnaseq' output
- check_fastq_strand_outputs: fetch Fastqs without fastq_strand.py outputs
"""

#######################################################################
# Imports
#######################################################################

import os
import logging
from bcftbx.qc.report import strip_ngs_extensions
from bcftbx.utils import AttributeDictionary
from ..analysis import AnalysisFastq
from ..analysis import split_sample_name
from ..fastq_utils import group_fastqs_by_name
from ..fastq_utils import remove_index_fastqs
from ..metadata import AnalysisProjectQCDirInfo
from ..tenx.cellplex import CellrangerMultiConfigCsv
from .fastq_strand import Fastqstrand
from .cellranger import CellrangerMulti
from .modules import QCDir
from .modules.cellranger_count import CellrangerCount
from .modules.fastqc import Fastqc
from .modules.fastq_screen import FastqScreen
from .modules.sequence_lengths import SequenceLengths

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Classes
#######################################################################

class QCOutputs:
    """
    Class to detect and characterise QC outputs

    On instantiation this class scans the supplied
    directory to identify and classify various
    artefacts which are typically produced by the
    QC pipeline.

    The following attributes are available:

    - fastqs: sorted list of Fastq names
    - reads: list of reads (e.g. 'r1', 'r2', 'i1' etc)
    - samples: sorted list of sample names extracted
      from Fastqs
    - seq_data_samples: sorted list of samples with
      biological data (rather than e.g. feature
      barcodes)
    - bams: sorted list of BAM file names
    - organisms: sorted list of organism names
    - fastq_screens: sorted list of screen names
    - cellranger_references: sorted list of reference
      datasets used with 10x pipelines
    - cellranger_probe_sets: sorted list of probe set
      files used with 10x pipelines
    - multiplexed_samples: sorted list of sample names
      for multiplexed samples (e.g. 10x CellPlex)
    - outputs: list of QC output categories detected (see
      below for valid values)
    - output_files: list of absolute paths to QC output
      files
    - software: dictionary with information on the
      QC software packages
    - stats: AttrtibuteDictionary with useful stats from
      across the project
    - config_files: list of QC configuration files found
      in the QC directory (see below for valid values)

    The following are valid values for the 'outputs'
    property:

    - 'fastqc_[r1...]'
    - 'screens_[r1...]'
    - 'strandedness'
    - 'sequence_lengths'
    - 'picard_insert_size_metrics'
    - 'rseqc_genebody_coverage'
    - 'rseqc_infer_experiment'
    - 'qualimap_rnaseq'
    - 'icell8_stats'
    - 'icell8_report'
    - 'cellranger_count'
    - 'cellranger_multi'
    - 'cellranger-atac_count'
    - 'cellranger-arc_count'
    - 'multiqc'
    - 'extra_outputs'

    The following are valid values for the 'config_files'
    property:

    - fastq_strand.conf
    - 10x_multi_config.csv
    - 10x_multi_config.<SAMPLE>.csv
    - libraries.<SAMPLE>.csv

    Arguments:
      qc_dir (str): path to directory to examine
      fastq_attrs (BaseFastqAttrs): (optional) class for
        extracting data from Fastq names
    """
    def __init__(self,qc_dir,fastq_attrs=None):
        # Directory to examine
        self.qc_dir = os.path.abspath(qc_dir)
        # How to handle Fastq names
        if fastq_attrs is not None:
            self.fastq_attrs = fastq_attrs
        else:
            self.fastq_attrs = AnalysisFastq
        # QC metadata
        self.qc_info = None
        qc_info_file = os.path.join(self.qc_dir,"qc.info")
        if os.path.exists(qc_info_file):
            self.qc_info = AnalysisProjectQCDirInfo(qc_info_file)
        # Properties
        self.fastqs = set()
        self.samples = []
        self.seq_data_samples = []
        self.bams = set()
        self.reads = []
        self.organisms = set()
        self.stats = AttributeDictionary(
            max_seqs=None,
            min_sequence_length_read={},
            max_sequence_length_read={},
            min_sequence_length=None,
            max_sequence_length=None,
        )
        self.software = {}
        self.outputs = set()
        self.output_files = []
        self.config_files = []
        # Raw QC output data
        self._qc_data = {}
        # Scan the directory for QC outputs
        self._collect_qc_outputs()
        # Report
        print("Collected %d output files" % len(self.output_files))
        print("Outputs:")
        for output in self.outputs:
            print("- %s" % output)
        print("Fastqs:")
        for fq in self.fastqs:
            print("- %s" % fq)
        print("BAMs:")
        for b in self.bams:
            print("- %s" % b)
        print("Samples:")
        for s in self.samples:
            print("- %s" % s)
        print("Reads:")
        for r in self.reads:
            print("- %s" % r)
        print("Organisms:")
        for o in self.organisms:
            print("- %s" % o)
        print("Biological samples:")
        for s in self.seq_data_samples:
            print("- %s" % s)
        print("Software:")
        for sw in self.software:
            print("- %s: %s" % (sw,','.join(self.software[sw])))
        print("Screens:")
        for scrn in self.fastq_screens:
            print("- %s" % scrn)
        print("QC config files:")
        for cfg in self.config_files:
            print("- %s" % cfg)

    def data(self,name):
        """
        Return the 'raw' data associated with a QC output

        Arguments:
          name (str): name identifier for a QC output
            (e.g. 'fastqc')

        Returns:
          AttributeDictionary containing the raw data
            associated with the named QC output.

        Raises:
          KeyError: if the name doesn't match a stored
            QC output.
        """
        return self._qc_data[name]

    def _collect_qc_outputs(self):
        """
        Internal: determine which QC outputs are present
        """
        # QC directory object
        qcdir = QCDir(self.qc_dir,fastq_attrs=self.fastq_attrs)
        # Get list of files
        print("Scanning contents of %s" % self.qc_dir)
        files = qcdir.file_list
        print("\t- %d objects found" % len(files))
        logger.debug("files: %s" % files)
        # Collect QC outputs from modules
        for qc_data in (
                self._collect_fastq_screens(qcdir),
                self._collect_fastqc(qcdir),
                self._collect_fastq_strand(files),
                self._collect_seq_lengths(qcdir),
                self._collect_picard_insert_size_metrics(self.qc_dir),
                self._collect_rseqc_genebody_coverage(self.qc_dir),
                self._collect_rseqc_infer_experiment(self.qc_dir),
                self._collect_qualimap_rnaseq(self.qc_dir),
                self._collect_icell8(self.qc_dir),
                self._collect_cellranger_count(qcdir),
                self._collect_cellranger_multi(self.qc_dir),
                self._collect_multiqc(self.qc_dir),
                self._collect_extra_outputs(self.qc_dir)):
            self._add_qc_outputs(qc_data)
        # Collect QC config files
        self.config_files = self._collect_qc_config_files(files)
        # Fastq screens
        fastq_screen = self.data('fastq_screen')
        self.fastq_screens = sorted(fastq_screen.screen_names)
        # Sequence lengths
        seq_lengths = self.data('sequence_lengths')
        if "sequence_lengths" in seq_lengths.tags:
            # Maximum number of sequences
            self.stats['max_seqs'] = seq_lengths.max_seqs
            # Maximum and minimum sequence lengths for each read
            for read in seq_lengths.reads:
                self.stats.min_sequence_length_read[read] = \
                    seq_lengths.min_seq_length[read]
                self.stats.max_sequence_length_read[read] = \
                    seq_lengths.max_seq_length[read]
            # Maximum and minimum sequence lengths for all reads
            min_lens = [seq_lengths.min_seq_length[r]
                        for r in seq_lengths.min_seq_length]
            if min_lens:
                self.stats['min_sequence_length'] = min(min_lens)
            max_lens = [seq_lengths.max_seq_length[r]
                        for r in seq_lengths.max_seq_length]
            if max_lens:
                self.stats['max_sequence_length'] = max(max_lens)
        # Cellranger count
        cellranger_count = self.data('cellranger_count')
        # Single library analyses reference data
        self.cellranger_references = [ref for ref in
                                      cellranger_count.references]
        # Cellranger multi
        cellranger_multi = self.data('cellranger_multi')
        # Multiplexed samples
        self.multiplexed_samples = sorted(cellranger_multi.multiplexed_samples,
                                          key=lambda s: split_sample_name(s))
        for reference_data in cellranger_multi.references:
            if reference_data not in self.cellranger_references:
                self.cellranger_references.append(reference_data)
        self.cellranger_references = sorted(self.cellranger_references)
        self.cellranger_probe_sets = cellranger_multi.probe_sets
        # Fastqs
        self.fastqs = sorted(list(self.fastqs))
        # BAMs
        self.bams = sorted(list(self.bams))
        # Samples
        samples = set([self.fastq_attrs(fq).sample_name
                       for fq in self.fastqs])
        for bam in self.bams:
            samples.add(self.fastq_attrs(bam).sample_name)
        for s in cellranger_count.samples:
            samples.add(s)
        self.samples = sorted(list(samples),
                              key=lambda s: split_sample_name(s))
        # Samples with biological sequence data
        if self.qc_info is not None and self.qc_info.seq_data_samples:
            # Explicitly defined in QC metadata
            self.seq_data_samples.extend(
                [s for s in self.qc_info.seq_data_samples.split(',')
                 if s in samples])
        for cf in self.config_files:
            # Implicitly defined in 10x multi config file(s) as GEX data
            if cf.startswith('10x_multi_config.'):
                cf = CellrangerMultiConfigCsv(os.path.join(self.qc_dir,cf))
                self.seq_data_samples.extend([s for s in cf.gex_libraries
                                              if s in samples])
        if not self.seq_data_samples:
            # Nothing defined - assume all samples are biological data
            self.seq_data_samples = [s for s in self.samples]
        self.seq_data_samples = sorted(list(set(self.seq_data_samples)))
        # Determine reads
        reads = set()
        for fastq in self.fastqs:
            fq = self.fastq_attrs(fastq)
            if fq.is_index_read:
                reads.add("i%s" % (fq.read_number
                                   if fq.read_number is not None else '1'))
            else:
                reads.add("r%s" % (fq.read_number
                                   if fq.read_number is not None else '1'))
        self.reads = sorted(list(reads))
        # Organisms
        self.organisms = sorted(list(self.organisms))
        # Sort outputs
        self.outputs = sorted(self.outputs)
        # Sort output files
        self.output_files = sorted(self.output_files)

    def _add_qc_outputs(self,data):
        """
        Internal: store 'raw' data for a specific QC output

        The data will be an AttributeDictionary instance
        returned from one of the ``_collect_...`` methods,
        which should contain at minimum the following
        attributes:

        - name: identifier for the QC output (e.g. 'fastqc')
        - fastqs: list of Fastq names associated with the
          output
        - software: a dictionary compromising the names of
          associated software packages as keys, and a list
          of versions as the values
        - output_files: a list of associated output files
        - tags: a list of output 'classes' associated with
          the QC output (e.g. 'fastqc_r1')

        The data can be accessed using the name identifier,
        via the 'data' method; the information about Fastqs,
        software, output files and tags are added to the
        global properties.
        """
        print("- adding data for '%s'" % data.name)
        if data.name in self._qc_data:
            raise KeyError("'%s': data already stored against "
                           "this name")
        ##print(data)
        # Store raw data
        self._qc_data[data.name] = data
        # Output files
        for f in data.output_files:
            self.output_files.append(f)
        # Fastqs
        try:
            for fq in data.fastqs:
                self.fastqs.add(fq)
        except AttributeError:
            pass
        # BAMs
        try:
            for bam in data.bam_files:
                self.bams.add(bam)
        except AttributeError:
            pass
        # Organisms
        try:
            for organism in data.organisms:
                self.organisms.add(organism)
        except AttributeError:
            pass
        # Software versions
        for sw in data.software.keys():
            self.software[sw] = [str(v) if v else '?'
                                 for v in data.software[sw]]
        # Tags
        for tag in data.tags:
            self.outputs.add(tag)

    def _read_versions_file(self,f,pkgs=None):
        """
        Internal: extract software info from 'versions' file

        'versions' files (typically named ``_versions``)
        should consist of one or more lines of text, with
        each line comprising a software package name and
        a version number, separated by a tab character.

        Returns a dictionary where package names are keys,
        and the corresponding values are lists of versions.

        If an existing dictionary is supplied via the 'pkgs'
        argument then any package information is added to this
        dictionary; otherwise an empty dictionary is created
        and populated.
        """
        if pkgs is not None:
            software = pkgs
        else:
            software = dict()
        if not os.path.exists(f):
            return software
        try:
            with open(f,'rt') as fp:
                for line in fp:
                    try:
                        pkg,version = line.strip().split('\t')
                    except IndexError:
                        continue
                    try:
                        software[pkg].append(version)
                    except KeyError:
                        software[pkg] = [version]
        except Exception as ex:
            logger.warning("%s: unable to extract versions: %s" %
                           (f,ex))
        return software

    def _collect_fastq_screens(self,qcdir):
        """
        Collect information on FastqScreen outputs

        Returns an AttributeDictionary with the following
        attributes:

        - name: set to 'fastq_screen'
        - software: dictionary of software and versions
        - screen_names: list of associated panel names
        - fastqs: list of associated Fastq names
        - fastqs_for_screen: dictionary of panel names
          and lists of Fastq names associated with each
          panel
        - output_files: list of associated output files
        - tags: list of associated output classes

        Arguments:
          qcdir (QCDir): QC directory object to examine
        """
        return FastqScreen.collect_qc_outputs(qcdir)

    def _collect_fastqc(self,qcdir):
        """
        Collect information on FastQC outputs

        Returns an AttributeDictionary with the following
        attributes:

        - name: set to 'fastqc'
        - software: dictionary of software and versions
        - fastqs: list of associated Fastq names
        - output_files: list of associated output files
        - tags: list of associated output classes

        Arguments:
          qcdir (QCDir): QC directory object to examine
        """
        return Fastqc.collect_qc_outputs(qcdir)

    def _collect_fastq_strand(self,files):
        """
        Collect information on FastqStrand outputs

        Returns an AttributeDictionary with the following
        attributes:

        - name: set to 'fastq_strand'
        - software: dictionary of software and versions
        - fastqs: list of associated Fastq names
        - output_files: list of associated output files
        - tags: list of associated output classes

        Arguments:
          files (list): list of file names to examine.
        """
        versions = set()
        output_files = list()
        fastqs = set()
        tags = set()
        # Look for fastq_strand outputs
        fastq_strand = list(filter(lambda f:
                                   f.endswith("_fastq_strand.txt"),
                                   files))
        logger.debug("fastq_strand: %s" % fastq_strand)
        print("\t- %d fastq_strand files" % len(fastq_strand))
        if fastq_strand:
            tags.add("strandedness")
            for f in fastq_strand:
                fq = self.fastq_attrs(os.path.splitext(f)[0])
                fastqs.add(
                    os.path.basename(
                        os.path.splitext(f)[0])[:-len("_fastq_strand")])
                versions.add(Fastqstrand(f).version)
            # Store the fastq_strand files
            output_files.extend(fastq_strand)
        # Return collected information
        if versions:
            software = { 'fastq_strand': sorted(list(versions)) }
        else:
            software = {}
        return AttributeDictionary(
            name='fastq_strand',
            software=software,
            fastqs=sorted(list(fastqs)),
            output_files=output_files,
            tags=sorted(list(tags))
        )

    def _collect_seq_lengths(self,qcdir):
        """
        Collect information on sequence length outputs

        Returns an AttributeDictionary with the following
        attributes:

        - name: set to 'sequence_lengths'
        - software: dictionary of software and versions
        - max_seqs: maximum number of sequences found in
          a single Fastq
        - min_seq_length: dictionary with minimum sequence
          lengths for each read
        - max_seq_length: dictionary with maximum sequence
          lengths for each read
        - reads: list of read IDs (e.g. 'r1', 'i2')
        - fastqs: list of associated Fastq names
        - output_files: list of associated output files
        - tags: list of associated output classes

        Arguments:
          qcdir (QCDir): QC directory object to examine
        """
        return SequenceLengths.collect_qc_outputs(qcdir)

    def _collect_picard_insert_size_metrics(self,qc_dir):
        """
        Collect information on Picard CollectInsertSizeMetrics outputs

        Returns an AttributeDictionary with the following
        attributes:

        - name: set to 'picard_collect_insert_size_metrics'
        - software: dictionary of software and versions
        - organisms: list of organisms with associated
          outputs
        - bam_files: list of associated BAM file names
        - output_files: list of associated output files
        - tags: list of associated output classes

        Arguments:
          qc_dir (str): top-level directory to look under.
        """
        software = {}
        bam_files = set()
        output_files = list()
        tags = set()
        # Look for Picard CollectInsertSizeMetrics outputs
        organisms = set()
        picard_dir = os.path.join(qc_dir,"picard")
        if os.path.isdir(picard_dir):
            # Look for subdirs with organism names
            for d in filter(
                    lambda dd:
                    os.path.isdir(os.path.join(picard_dir,dd)),
                    os.listdir(picard_dir)):
                # Check for outputs
                for f in filter(
                        lambda ff:
                        ff.endswith(".insert_size_metrics.txt"),
                        os.listdir(os.path.join(picard_dir,d))):
                    name = f[:-len(".insert_size_metrics.txt")]
                    outputs = picard_collect_insert_size_metrics_output(
                        name,
                        prefix=os.path.join(picard_dir,d))
                    if all([os.path.exists(f) for f in outputs]):
                        # All outputs present
                        organisms.add(d)
                        bam_files.add(name)
                        output_files.extend(outputs)
                # Check for software information
                software = self._read_versions_file(
                    os.path.join(picard_dir,d,"_versions"),
                    software)
        if organisms:
            tags.add("picard_insert_size_metrics")
            if not software:
                software['picard'] = [None]
        # Look for collated insert sizes files
        for f in filter(
                lambda ff:
                ff.startswith("insert_sizes.") and ff.endswith(".tsv"),
                os.listdir(qc_dir)):
            tags.add("collated_insert_sizes")
            output_files.append(os.path.join(qc_dir,f))
            organisms.add(f[len("insert_sizes."):-len(".tsv")])
        # Return collected information
        return AttributeDictionary(
            name='picard_collect_insert_size_metrics',
            software=software,
            bam_files=sorted(list(bam_files)),
            organisms=sorted(list(organisms)),
            output_files=output_files,
            tags=sorted(list(tags))
        )

    def _collect_rseqc_genebody_coverage(self,qc_dir):
        """
        Collect information on RSeQC geneBody_coverage.py outputs

        Returns an AttributeDictionary with the following
        attributes:

        - name: set to 'rseqc_genebody_coverage'
        - software: dictionary of software and versions
        - organisms: list of organisms with associated
          outputs
        - output_files: list of associated output files
        - tags: list of associated output classes

        Arguments:
          qc_dir (str): top-level directory to look under.
        """
        software = {}
        output_files = list()
        tags = set()
        # Look for RSeQC geneBody_coverage.py outputs
        organisms = set()
        genebody_cov_dir = os.path.join(qc_dir,
                                        "rseqc_genebody_coverage")
        if os.path.isdir(genebody_cov_dir):
            # Look for subdirs with organism names
            for d in filter(
                    lambda dd:
                    os.path.isdir(os.path.join(genebody_cov_dir,dd)),
                    os.listdir(genebody_cov_dir)):
                # Check for outputs
                for f in filter(
                        lambda ff:
                        ff.endswith(".geneBodyCoverage.txt"),
                        os.listdir(os.path.join(genebody_cov_dir,d))):
                    name = f[:-len(".geneBodyCoverage.txt")]
                    outputs = rseqc_genebody_coverage_output(
                        name,
                        prefix=os.path.join(genebody_cov_dir,d))
                    if all([os.path.exists(f) for f in outputs]):
                        # All outputs present
                        organisms.add(d)
                        output_files.extend(outputs)
                # Check for software information
                software = self._read_versions_file(
                    os.path.join(genebody_cov_dir,
                                 d,"_versions"),
                    software)
        if organisms:
            tags.add("rseqc_genebody_coverage")
            if not software:
                software['rseqc:genebody_coverage'] = [None]
        # Return collected information
        return AttributeDictionary(
            name='rseqc_genebody_coverage',
            software=software,
            organisms=sorted(list(organisms)),
            output_files=output_files,
            tags=sorted(list(tags))
        )

    def _collect_rseqc_infer_experiment(self,qc_dir):
        """
        Collect information on RSeQC infer_experiment.py outputs

        Returns an AttributeDictionary with the following
        attributes:

        - name: set to 'rseqc_infer_experiment'
        - software: dictionary of software and versions
        - organisms: list of organisms with associated
          outputs
        - bam_files: list of associated BAM file names
        - output_files: list of associated output files
        - tags: list of associated output classes

        Arguments:
          qc_dir (str): top-level directory to look under.
        """
        software = {}
        bam_files = set()
        output_files = list()
        tags = set()
        # Look for RSeQC infer_experiment.py outputs
        organisms = set()
        infer_experiment_dir = os.path.join(qc_dir,
                                            "rseqc_infer_experiment")
        if os.path.isdir(infer_experiment_dir):
            # Look for subdirs with organism names
            for d in filter(
                    lambda dd:
                    os.path.isdir(os.path.join(infer_experiment_dir,dd)),
                    os.listdir(infer_experiment_dir)):
                # Check for outputs
                for f in filter(
                        lambda ff:
                        ff.endswith(".infer_experiment.log"),
                        os.listdir(os.path.join(infer_experiment_dir,d))):
                    name = f[:-len(".infer_experiment.log")]
                    organisms.add(d)
                    bam_files.add(name)
                    output_files.append(os.path.join(infer_experiment_dir,
                                                     d,f))
                # Check for software information
                software = self._read_versions_file(
                    os.path.join(infer_experiment_dir,
                                 d,"_versions"),
                    software)
        if organisms:
            tags.add("rseqc_infer_experiment")
            if not software:
                software['rseqc:infer_experiment'] = [None]
        # Return collected information
        return AttributeDictionary(
            name='rseqc_infer_experiment',
            software=software,
            bam_files=sorted(list(bam_files)),
            organisms=sorted(list(organisms)),
            output_files=output_files,
            tags=sorted(list(tags))
        )

    def _collect_qualimap_rnaseq(self,qc_dir):
        """
        Collect information on Qualimap 'rnaseq' outputs

        Returns an AttributeDictionary with the following
        attributes:

        - name: set to 'qualimap_rnaseq'
        - software: dictionary of software and versions
        - organisms: list of organisms with associated
          outputs
        - bam_files: list of associated BAM file names
        - output_files: list of associated output files
        - tags: list of associated output classes

        Arguments:
          qc_dir (str): top-level directory to look under.
        """
        software = {}
        bam_files = set()
        output_files = list()
        tags = set()
        # Look for Qualimap 'rnaseq' outputs
        organisms = set()
        qualimap_dir = os.path.join(qc_dir,"qualimap-rnaseq")
        if os.path.isdir(qualimap_dir):
            # Look for subdirs with organism names
            for d in filter(
                    lambda dd:
                    os.path.isdir(os.path.join(qualimap_dir,dd)),
                    os.listdir(qualimap_dir)):
                # Look for subdirs with BAM file names
                organism_dir = os.path.join(qualimap_dir,d)
                for bam in filter(
                        lambda dd:
                        os.path.isdir(os.path.join(organism_dir,dd)),
                        os.listdir(organism_dir)):
                    # Check for Qualimap rnaseq outputs
                    for f in filter(
                            lambda ff:
                            ff == "qualimapReport.html",
                            os.listdir(os.path.join(organism_dir,bam))):
                        outputs = [os.path.join(organism_dir,bam,ff)
                                   for ff in ('qualimapReport.html',
                                              'rnaseq_qc_results.txt')]
                        if all([os.path.exists(ff) for ff in outputs]):
                            # All outputs present
                            organisms.add(d)
                            bam_files.add(bam)
                            output_files.extend(outputs)
                            # Add additional outputs (CSS, images etc)
                            for subdir in ('css',
                                           'images_qualimapReport',
                                           'raw_data_qualimapReport',):
                                dd = os.path.join(organism_dir,bam,subdir)
                                if os.path.exists(dd):
                                    extra_files = [os.path.join(dd,ff)
                                                   for ff in os.listdir(dd)]
                                    output_files.extend(extra_files)
                # Check for software information
                software = self._read_versions_file(
                    os.path.join(organism_dir,"_versions"),
                    software)
        if organisms:
            tags.add("qualimap_rnaseq")
            if not software:
                software['qualimap'] = [None]
        # Return collected information
        return AttributeDictionary(
            name='qualimap_rnaseq',
            software=software,
            bam_files=sorted(list(bam_files)),
            organisms=sorted(list(organisms)),
            output_files=output_files,
            tags=sorted(list(tags))
        )

    def _collect_icell8(self,qc_dir):
        """
        Collect information on ICell8 reports

        Returns an AttributeDictionary with the following
        attributes:

        - name: set to 'icell8'
        - software: dictionary of software and versions
        - fastqs: list of associated Fastq names
        - output_files: list of associated output files
        - tags: list of associated output classes

        Arguments:
          qc_dir (str): top-level directory to look under.
        """
        output_files = list()
        tags = set()
        # Look for ICELL8 outputs
        icell8_top_dir = os.path.dirname(qc_dir)
        print("Checking for ICELL8 reports in %s/stats" %
              icell8_top_dir)
        icell8_stats_xlsx = os.path.join(icell8_top_dir,
                                         "stats",
                                         "icell8_stats.xlsx")
        if os.path.exists(icell8_stats_xlsx):
            tags.add("icell8_stats")
            output_files.append(icell8_stats_xlsx)
        icell8_report_html = os.path.join(icell8_top_dir,
                                          "icell8_processing.html")
        if os.path.exists(icell8_report_html):
            tags.add("icell8_report")
            output_files.append(icell8_report_html)
        # Return collected information
        return AttributeDictionary(
            name='icell8',
            software={},
            fastqs=[],
            output_files=output_files,
            tags=sorted(list(tags))
        )

    def _collect_cellranger_count(self,qcdir):
        """
        Collect information on Cellranger count outputs

        Returns an AttributeDictionary with the following
        attributes:

        - name: set to 'cellranger_count'
        - software: dictionary of software and versions
        - references: list of associated reference datasets
        - fastqs: list of associated Fastq names
        - samples: list of associated sample names
        - pipelines: list of tuples defining 10x pipelines
          in the form (name,version,reference)
        - samples_by_pipeline: dictionary with lists of
          sample names associated with each 10x pipeline
          tuple
        - output_files: list of associated output files
        - tags: list of associated output classes

        Arguments:
          qcdir (QCDir): QC directory object to examine
        """
        return CellrangerCount.collect_qc_outputs(qcdir)

    def _collect_cellranger_multi(self,qc_dir):
        """
        Collect information on Cellranger multi outputs

        Returns an AttributeDictionary with the following
        attributes:

        - name: set to 'cellranger_multi'
        - software: dictionary of software and versions
        - references: list of associated reference datasets
        - probe_sets: list of associated probe sets
        - fastqs: list of associated Fastq names
        - multiplexed_samples: list of associated multiplexed
          sample names
        - pipelines: list of tuples defining 10x pipelines
          in the form (name,version,reference)
        - samples_by_pipeline: dictionary with lists of
          multiplexed sample names associated with each 10x
          pipeline tuple
        - output_files: list of associated output files
        - tags: list of associated output classes

        Arguments:
          qc_dir (str): top-level directory to look under.
        """
        software = {}
        output_files = list()
        multiplexed_samples = set()
        cellranger_references = set()
        cellranger_probe_sets = set()
        samples_by_pipeline = dict()
        tags = set()
        # Look for cellranger multi outputs
        cellranger_multi_dir = os.path.join(qc_dir,
                                            "cellranger_multi")
        print("Checking for cellranger multi outputs under %s" %
              cellranger_multi_dir)
        if os.path.isdir(cellranger_multi_dir):
            cellranger_name = None
            versions = set()
            cellranger_multi_samples = {}
            for ver in filter(
                    lambda f:
                    os.path.isdir(os.path.join(cellranger_multi_dir,f)),
                    os.listdir(cellranger_multi_dir)):
                cellranger_multi_samples[ver] = {}
                for ref in filter(
                        lambda f:
                        os.path.isdir(os.path.join(cellranger_multi_dir,ver,f)),
                        os.listdir(os.path.join(cellranger_multi_dir,ver))):
                    # Check putative reference dataset names
                    cellranger_multi_samples[ver][ref] = []
                    cellranger_multi = CellrangerMulti(
                        os.path.join(
                            cellranger_multi_dir,
                            ver,
                            ref))
                    for smpl in cellranger_multi.sample_names:
                        cellranger_multi_samples[ver][ref].append(smpl)
                        try:
                            output_files.append(cellranger_multi.web_summary(smpl))
                            output_files.append(cellranger_multi.metrics_csv(smpl))
                            cellranger_name = cellranger_multi.pipeline_name
                            if cellranger_name is None:
                                cellranger_name = 'cellranger'
                            if cellranger_multi.reference_data:
                                cellranger_references.add(
                                    cellranger_multi.reference_data)
                            if cellranger_multi.probe_set:
                                cellranger_probe_sets.add(
                                    cellranger_multi.probe_set)
                        except OSError:
                            pass
                    # Add outputs, samples and version
                    if cellranger_multi_samples[ver][ref]:
                        tags.add("cellranger_multi")
                        versions.add(ver)
                    for smpl in cellranger_multi_samples[ver][ref]:
                        multiplexed_samples.add(smpl)
            # Store sample lists associated with pipeline,
            # version and reference dataset
            for version in cellranger_multi_samples:
                for reference in cellranger_multi_samples[version]:
                    pipeline_key = (cellranger_name,version,reference)
                    samples_by_pipeline[pipeline_key] = \
                        [s for s in
                         cellranger_multi_samples[version][reference]]
            # Store cellranger versions
            if cellranger_name and versions:
                if cellranger_name not in software:
                    software[cellranger_name] = list(versions)
                else:
                    for v in list(versions):
                        if v not in software[cellranger_name]:
                            software[cellranger_name].append(v)
                software[cellranger_name] = sorted(software[cellranger_name])
        # Return collected information
        return AttributeDictionary(
            name='cellranger_multi',
            software=software,
            references=sorted(list(cellranger_references)),
            probe_sets=sorted(list(cellranger_probe_sets)),
            fastqs=[],
            multiplexed_samples=sorted(list(multiplexed_samples)),
            pipelines=sorted([p for p in samples_by_pipeline]),
            samples_by_pipeline=samples_by_pipeline,
            output_files=output_files,
            tags=sorted(list(tags))
        )

    def _collect_multiqc(self,qc_dir):
        """
        Collect information on MultiQC reports

        Returns an AttributeDictionary with the following
        attributes:

        - name: set to 'multiqc'
        - software: dictionary of software and versions
        - fastqs: list of associated Fastq names
        - output_files: list of associated output files
        - tags: list of associated output classes

        Arguments:
          qc_dir (str): top-level directory to look under.
        """
        version = None
        output_files = list()
        tags = set()
        # Look for MultiQC report
        multiqc_dir = os.path.dirname(qc_dir)
        print("Checking for MultiQC report in %s" % multiqc_dir)
        multiqc_report = os.path.join(multiqc_dir,
                                      "multi%s_report.html"
                                      % os.path.basename(qc_dir))
        if os.path.isfile(multiqc_report):
            tags.add("multiqc")
            output_files.append(multiqc_report)
            # Try to locate version from HTML file
            # Look for line like e.g.
            # <a href="http://multiqc.info" target="_blank">MultiQC v1.8</a>
            with open(multiqc_report,'rt') as fp:
                for line in fp:
                    if line.strip().startswith("<a href=\"http://multiqc.info\" target=\"_blank\">MultiQC v"):
                        try:
                            version = line.strip().split()[3][1:-4]
                        except Exception as ex:
                            logger.warning("Failed to extract MultiQC version "
                                           "from '%s': %s" % (line,ex))
                        break
        # Return collected information
        if version:
            software = { 'multiqc': [ version ] }
        else:
            software = {}
        return AttributeDictionary(
            name='multiqc',
            software=software,
            fastqs=[],
            output_files=output_files,
            tags=sorted(list(tags))
        )

    def _collect_extra_outputs(self,qc_dir):
        """
        Collect information on additional outputs

        Returns an AttributeDictionary with the following
        attributes:

        - name: set to 'extra_outputs'
        - software: dictionary of software and versions
        - output_files: list of associated output files
        - tags: list of associated output classes

        Arguments:
          qc_dir (str): top-level directory to look under.
        """
        version = None
        software = dict()
        output_files = set()
        tags = set()
        # Look for extra_outputs.tsv
        extra_outputs_tsv = os.path.join(qc_dir,"extra_outputs.tsv")
        for output in ExtraOutputs(extra_outputs_tsv).outputs:
            output_files.add(os.path.join(qc_dir,output.file_path))
            if output.additional_files:
                # Add in any additional files or dirs
                for f in output.additional_files:
                    output_files.add(os.path.join(qc_dir,f))
        if output_files:
            tags.add("extra_outputs")
        # Return collected information
        return AttributeDictionary(
            name='extra_outputs',
            software=software,
            fastqs=[],
            output_files=sorted(list(output_files)),
            tags=sorted(list(tags))
        )

    def _collect_qc_config_files(self,files):
        """
        Collect information on QC config files

        Returns a list of QC configuration (if any) found
        in the supplied list (with the leading paths
        removed).

        The following configuration files are checked for:

        - fastq_strand.conf
        - 10x_multi_config.csv
        - 10x_multi_config.<SAMPLE>.csv
        - libraries.<SAMPLE>.csv

        Arguments:
          files (list): list of file names to examine.
        """
        config_files = set()
        # Strand and cellranger multi configs
        for f in ("fastq_strand.conf",
                  "10x_multi_config.csv",):
            if os.path.join(self.qc_dir,f) in files:
                config_files.add(f)
        # Cellranger-arc configs ('libraries.SAMPLE.csv') and per-sample
        # cellranger multi configs ('10x_multi_config.SAMPLE.csv')
        for f in list(
                filter(lambda f:
                       f.endswith(".csv") and
                       (os.path.basename(f).startswith("libraries.") or
                        os.path.basename(f).startswith("10x_multi_config.")),
                       files)):
            config_files.add(os.path.basename(f))
        return sorted(list(config_files))

class ExtraOutputs:
    """
    Class for handling files specifying external QC outputs

    Reads data from the supplied tab-delimited (TSV) file
    specifying one or more external QC output files.

    Each line in the file should have up to three items
    separated by tabs:

    - file or directory (relative to the qc dir)
    - text description (used in HTML)
    - optionally, comma-separated list of additional files
      or directories to include in the final ZIP archive
      (relative to the qc dir)

    Blank lines and lines starting with the '#' comment
    character are ignored.

    The data from each line of the file is then available
    via the 'outputs' attribute, which provides a list
    of 'AttributeDictionary' objects with the following
    properties:

    - 'file_path': relative path to the output file
    - 'description': associated description
    - 'additional_files': list of the associated files

    Arguments:
      tsv_file (str): path to the input TSV file
    """
    def __init__(self,tsv_file):
        self.tsv_file = os.path.abspath(tsv_file)
        self.outputs = list()
        self.__load()

    def __load(self):
        """
        Internal: loads data from the supplied TSV file
        """
        if not os.path.exists(self.tsv_file):
            return
        with open(self.tsv_file,'rt') as fp:
            for line in fp:
                line = line.strip()
                if line.startswith('#') or not line:
                    continue
                items = line.split('\t')
                if len(items) > 3:
                    raise IndexError("Bad line (too many items): '%s'" %
                                     line)
                try:
                    file_path = items[0]
                    description = items[1]
                except IndexError:
                    raise IndexError("Bad line (not enough items): '%s'" %
                                     line)
                try:
                    additional_files = items[2].split(',')
                except IndexError:
                    additional_files = None
                self.outputs.append(
                    AttributeDictionary(file_path=file_path,
                                        description=description,
                                        additional_files=additional_files))

#######################################################################
# Functions
#######################################################################

def fastq_strand_output(fastq):
    """
    Generate name for fastq_strand.py output

    Given a Fastq file name, the output from fastq_strand.py
    will look like:

    - {FASTQ}_fastq_strand.txt

    Arguments:
       fastq (str): name of Fastq file

    Returns:
       tuple: fastq_strand.py output (without leading paths)

    """
    return "%s_fastq_strand.txt" % strip_ngs_extensions(
        os.path.basename(fastq))

def picard_collect_insert_size_metrics_output(filen,prefix=None):
    """
    Generate names of Picard CollectInsertSizeMetrics output

    Given a Fastq or BAM file name, the output from Picard's
    CollectInsertSizeMetrics function will look like:

    - {PREFIX}/{FASTQ}.insert_size_metrics.txt
    - {PREFIX}/{FASTQ}.insert_size_histogram.pdf

    Arguments:
      filen (str): name of Fastq or BAM file
      prefix (str): optional directory to prepend to
        outputs

    Returns:
      tuple: CollectInsertSizeMetrics output (without leading
        paths)

    """
    outputs = []
    basename = os.path.basename(filen)
    while basename.split('.')[-1] in ('bam',
                                      'fastq',
                                      'gz'):
        basename = '.'.join(basename.split('.')[:-1])
    for ext in ('.insert_size_metrics.txt',
                '.insert_size_histogram.pdf'):
        outputs.append("%s%s" % (basename,ext))
    if prefix is not None:
        outputs = [os.path.join(prefix,f) for f in outputs]
    return tuple(outputs)

def rseqc_genebody_coverage_output(name,prefix=None):
    """
    Generate names of RSeQC geneBody_coverage.py output

    Given a basename, the output from geneBody_coverage.py
    will look like:

    - {PREFIX}/{NAME}.geneBodyCoverage.curves.png
    - {PREFIX}/{NAME}.geneBodyCoverage.r
    - {PREFIX}/{NAME}.geneBodyCoverage.txt

    Arguments:
      name (str): basename for output files
      prefix (str): optional directory to prepend to
        outputs

    Returns:
      tuple: geneBody_coverage.py output (without leading paths)

    """
    outputs = []
    for ext in ('.geneBodyCoverage.curves.png',
                '.geneBodyCoverage.r',
                '.geneBodyCoverage.txt'):
        outputs.append("%s%s" % (name,ext))
    if prefix is not None:
        outputs = [os.path.join(prefix,f) for f in outputs]
    return tuple(outputs)

def qualimap_rnaseq_output(prefix=None):
    """
    Generate names of Qualimap 'rnaseq' output

    The output from Qualimap 'rnaseq' are always::

    - {PREFIX}/qualimapReport.html
    - {PREFIX}/rnaseq_qc_results.txt
    - {PREFIX}/...

    Arguments:
      prefix (str): optional directory to prepend to
        outputs

    Returns:
      tuple: Qualimap 'rnaseq' output (without leading paths)

    """
    outputs = ['qualimapReport.html',
               'rnaseq_qc_results.txt']
    if prefix is not None:
        outputs = [os.path.join(prefix,f) for f in outputs]
    return tuple(outputs)

def check_fastq_strand_outputs(project,qc_dir,fastq_strand_conf,
                               fastqs=None,read_numbers=None):
    """
    Return Fastqs missing QC outputs from fastq_strand.py

    Returns a list of the Fastqs from a project for which
    one or more associated outputs from `fastq_strand.py`
    don't exist in the specified QC directory.

    Arguments:
      project (AnalysisProject): project to check the
        QC outputs for
      qc_dir (str): path to the QC directory (relative
        path is assumed to be a subdirectory of the
        project)
      fastq_strand_conf (str): path to a fastq_strand
        config file; strandedness QC outputs will be
        included unless the path is `None` or the
        config file doesn't exist. Relative path is
        assumed to be a subdirectory of the project
      fastqs (list): optional list of Fastqs to check
        against (defaults to Fastqs from the project)
      read_numbers (list): read numbers to predict
        outputs for

    Returns:
      List: list of Fastq file "pairs" with missing
        outputs; pairs are (R1,R2) tuples, with 'R2'
        missing if only one Fastq is used for the
        strandedness determination.
    """
    # Sort out QC directory
    if not os.path.isabs(qc_dir):
        qc_dir = os.path.join(project.dirn,qc_dir)
    # Sort out fastq_strand config file
    if fastq_strand_conf is not None:
        if not os.path.isabs(fastq_strand_conf):
            fastq_strand_conf = os.path.join(project.dirn,
                                             fastq_strand_conf)
    if not os.path.exists(fastq_strand_conf):
        # No conf file, nothing to check
        return list()
    if not fastqs:
        fastqs_in = project.fastqs
    else:
        fastqs_in = fastqs
    fastq_pairs = set()
    for fq_group in group_fastqs_by_name(
            remove_index_fastqs(fastqs_in,
                                project.fastq_attrs),
            fastq_attrs=project.fastq_attrs):
        # Assemble Fastq pairs based on read numbers
        if read_numbers:
            fq_pair = tuple([fq_group[r-1] for r in read_numbers])
        else:
            fq_pair = tuple([fq_group[0]])
        # Strand stats output
        output = os.path.join(qc_dir,
                              fastq_strand_output(fq_pair[0]))
        if not os.path.exists(output):
            fastq_pairs.add(fq_pair)
    return sorted(list(fastq_pairs))
