#!/usr/bin/env python
#
#     qc/outputs: utilities to predict and check QC pipeline outputs
#     Copyright (C) University of Manchester 2019-2022 Peter Briggs
#
"""
Provides utility classes and functions for QC outputs.

Provides the following class:

- QCOutputs: detect and characterise QC outputs

Provides the following functions:

- fastq_screen_output: get names for fastq_screen outputs
- fastqc_output: get names for FastQC outputs
- fastq_strand_output: get name for fastq_strand.py output
- picard_collect_insert_size_metrics_output: get names for Picard
  CollectInsertSizeMetrics output
- rseqc_genebody_coverage_output: get names for RSeQC geneBody_coverage.py
  output
- qualimap_rnaseq_output: get names for Qualimap 'rnaseq' output
- cellranger_count_output: get names for cellranger count output
- cellranger_atac_count_output: get names for cellranger-atac count output
- cellranger_arc_count_output: get names for cellranger-arc count output
- cellranger_multi_output: get names for cellranger multi output
- check_fastq_strand_outputs: fetch Fastqs without fastq_strand.py outputs
- check_cellranger_count_outputs: fetch sample names without cellranger
  count outputs
- check_cellranger_atac_count_outputs: fetch sample names without
  cellranger-atac count outputs
- check_cellranger_arc_count_outputs: fetch sample names without
  cellranger-arc count outputs
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
from ..tenx_genomics_utils import CellrangerMultiConfigCsv
from .constants import FASTQ_SCREENS
from .fastqc import Fastqc
from .fastq_screen import Fastqscreen
from .fastq_strand import Fastqstrand
from .cellranger import CellrangerCount
from .cellranger import CellrangerMulti
from .protocols import get_read_numbers
from .seqlens import SeqLens

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
    - bams: sorted list of BAM file names
    - organisms: sorted list of organism names
    - fastq_screens: sorted list of screen names
    - cellranger_references: sorted list of reference
      datasets used with 10x pipelines
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

    The following are valid values for the 'config_files'
    property:

    - fastq_strand.conf
    - 10x_multi_config.csv
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
        # Properties
        self.fastqs = set()
        self.samples = []
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
        # Get list of files
        print("Scanning contents of %s" % self.qc_dir)
        files = [os.path.join(self.qc_dir,f)
                 for f in os.listdir(self.qc_dir)]
        print("\t- %d objects found" % len(files))
        logger.debug("files: %s" % files)
        # Collect QC outputs
        for qc_data in (
                self._collect_fastq_screens(files),
                self._collect_fastqc(files),
                self._collect_fastq_strand(files),
                self._collect_seq_lengths(files),
                self._collect_picard_insert_size_metrics(self.qc_dir),
                self._collect_rseqc_genebody_coverage(self.qc_dir),
                self._collect_rseqc_infer_experiment(self.qc_dir),
                self._collect_qualimap_rnaseq(self.qc_dir),
                self._collect_icell8(self.qc_dir),
                self._collect_cellranger_count(self.qc_dir),
                self._collect_cellranger_multi(self.qc_dir),
                self._collect_multiqc(self.qc_dir),):
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
        self.multiplexed_samples = cellranger_multi.multiplexed_samples
        for reference_data in cellranger_multi.references:
            if reference_data not in self.cellranger_references:
                self.cellranger_references.append(reference_data)
        self.cellranger_references = sorted(self.cellranger_references)
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

        'versions' files (typically named ``__versions``)
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

    def _collect_fastq_screens(self,files):
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
          files (list): list of file names to examine.
        """
        versions = set()
        output_files = list()
        fastqs = set()
        fastqs_for_screen = dict()
        screen_names = set()
        tags = set()
        # Look for legacy screen files
        legacy_screens = list(filter(lambda f:
                                     f.endswith("_screen.txt") or
                                     f.endswith("_screen.png"),
                                     files))
        logger.debug("Screens (legacy): %s" % legacy_screens)
        # Look for new-style screen files
        screens = list(filter(lambda f:
                              "_screen_" in os.path.basename(f) and
                              (f.endswith(".txt") or
                               f.endswith(".png")),
                               files))
        logger.debug("Screens: %s" % screens)
        print("\t- %d fastq_screen files" % (len(legacy_screens) +
                                             len(screens)))
        if legacy_screens:
            for screen_name in FASTQ_SCREENS:
                # Explicitly check for legacy screens
                for screen in list(filter(lambda s:
                                          s.endswith("_%s_screen.txt" %
                                                     screen_name),
                                          legacy_screens)):
                    fq = self.fastq_attrs(screen[:-len("_screen.txt")])
                    fastq_name = str(fq)[:-len("_%s" % screen_name)]
                    tags.add("screens_%s%s" %
                             (('i' if fq.is_index_read else 'r'),
                              (fq.read_number
                               if fq.read_number is not None else '1')))
                    # Store general information
                    fastqs.add(fastq_name)
                    screen_names.add(screen_name)
                    versions.add(Fastqscreen(screen).version)
                    # Store Fastq against screen name
                    if screen_name not in fastqs_for_screen:
                        fastqs_for_screen[screen_name] = set()
                    fastqs_for_screen[screen_name].add(fastq_name)
            # Store the legacy screen files
            output_files.extend(legacy_screens)
        if screens:
            # Pull out the Fastq names from the .txt files
            for screen in list(filter(lambda s: s.endswith(".txt"),
                                      screens)):
                # Assume that names are 'FASTQ_screen_SCREENNAME.txt'
                fastq_name,screen_name = os.path.basename(screen)\
                                         [:-len(".txt")].\
                                         split("_screen_")
                fq = self.fastq_attrs(fastq_name)
                tags.add("screens_%s%s" %
                         (('i' if fq.is_index_read else 'r'),
                          (fq.read_number
                           if fq.read_number is not None else '1')))
                # Store general information
                fastqs.add(fastq_name)
                screen_names.add(screen_name)
                versions.add(Fastqscreen(screen).version)
                # Store Fastq against screen name
                if screen_name not in fastqs_for_screen:
                    fastqs_for_screen[screen_name] = set()
                fastqs_for_screen[screen_name].add(fastq_name)
            # Store the screen files
            output_files.extend(screens)
        # Return collected information
        if versions:
            software = { 'fastq_screen': sorted(list(versions)) }
        else:
            software = {}
        return AttributeDictionary(
            name='fastq_screen',
            software=software,
            screen_names=sorted(list(screen_names)),
            fastqs=sorted(list(fastqs)),
            fastqs_for_screen={ s: sorted(list(fastqs_for_screen[s]))
                                for s in fastqs_for_screen },
            output_files=output_files,
            tags=sorted(list(tags))
        )

    def _collect_fastqc(self,files):
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
          files (list): list of file names to examine.
        """
        versions = set()
        output_files = list()
        fastqs = set()
        tags = set()
        # Look for fastqc outputs
        fastqcs = list(
            filter(lambda f:
                   f.endswith("_fastqc") and
                   os.path.exists("%s.html" % f) and
                   os.path.exists(os.path.join(f,"summary.txt")) and
                   os.path.exists(os.path.join(f,"fastqc_data.txt")),
                   files))
        logger.debug("Fastqc: %s" % fastqcs)
        print("\t- %d fastqc files" % len(fastqcs))
        if fastqcs:
            versions = set()
            # Pull out the Fastq names from the Fastqc files
            for fastqc in fastqcs:
                f = os.path.basename(fastqc)[:-len("_fastqc")]
                fastqs.add(f)
                fq = self.fastq_attrs(f)
                tags.add("fastqc_%s%s" %
                         (('i' if fq.is_index_read else 'r'),
                          (fq.read_number
                           if fq.read_number is not None else '1')))
                versions.add(Fastqc(fastqc).version)
            # Store the fastqc files needed for reporting
            output_files.extend(["%s.html" % f for f in fastqcs])
            output_files.extend([os.path.join(f,"summary.txt")
                                 for f in fastqcs])
            output_files.extend([os.path.join(f,"fastqc_data.txt")
                                 for f in fastqcs])
            # Fastqc plot images
            for png in ("per_base_quality",):
                output_files.extend([os.path.join(f,
                                                  "Images",
                                                  "%s.png" % png)
                                     for f in fastqcs])
        # Return collected information
        if versions:
            software = { 'fastqc': sorted(list(versions)) }
        else:
            software = {}
        return AttributeDictionary(
            name='fastqc',
            software=software,
            fastqs=sorted(list(fastqs)),
            output_files=output_files,
            tags=sorted(list(tags))
        )

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

    def _collect_seq_lengths(self,files):
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
          files (list): list of file names to examine.
        """
        output_files = list()
        fastqs = set()
        reads = set()
        max_seqs = None
        min_seq_length = {}
        max_seq_length = {}
        tags = set()
        # Look for sequence length outputs
        seq_lens = list(filter(lambda f:
                               f.endswith("_seqlens.json"),
                               files))
        logger.debug("seq_lens: %s" % seq_lens)
        print("\t- %d sequence length files" % len(seq_lens))
        if seq_lens:
            tags.add("sequence_lengths")
            for f in seq_lens:
                fq = self.fastq_attrs(os.path.splitext(f)[0])
                read = '%s%s' % ('i' if fq.is_index_read else 'r',
                                 fq.read_number)
                fastqs.add(
                    os.path.basename(
                        os.path.splitext(f)[0])[:-len("_seqlens")])
                reads.add(read)
                seqlens_data = SeqLens(f)
                try:
                    # Try to extract the sequence lengths
                    max_seqs = max(max_seqs,seqlens_data.nreads)
                except Exception:
                    if not max_seqs:
                        max_seqs = seqlens_data.nreads
                if not read in min_seq_length:
                    min_seq_length[read] = seqlens_data.min_length
                else:
                    min_seq_length[read] = min(min_seq_length[read],
                                               seqlens_data.min_length)
                if not read in max_seq_length:
                    max_seq_length[read] = seqlens_data.max_length
                else:
                    max_seq_length[read] = max(max_seq_length[read],
                                               seqlens_data.max_length)
            # Store the sequence length files
            output_files.extend(seq_lens)
        # Return collected information
        return AttributeDictionary(
            name='sequence_lengths',
            software={},
            max_seqs=max_seqs,
            min_seq_length=min_seq_length,
            max_seq_length=max_seq_length,
            reads=sorted(list(reads)),
            fastqs=sorted(list(fastqs)),
            output_files=output_files,
            tags=sorted(list(tags))
        )

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
                    os.path.join(picard_dir,d,"__versions"),
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
                                 d,"__versions"),
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
                    print("- %s" % f)
                    name = f[:-len(".infer_experiment.log")]
                    organisms.add(d)
                    bam_files.add(name)
                    output_files.append(os.path.join(infer_experiment_dir,
                                                     d,f))
                # Check for software information
                software = self._read_versions_file(
                    os.path.join(infer_experiment_dir,
                                 d,"__versions"),
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
                                           'images_qualimapReport'):
                                dd = os.path.join(organism_dir,bam,subdir)
                                if os.path.exists(dd):
                                    extra_files = [os.path.join(dd,ff)
                                                   for ff in os.listdir(dd)]
                                    output_files.extend(extra_files)
                # Check for software information
                software = self._read_versions_file(
                    os.path.join(organism_dir,"__versions"),
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

    def _collect_cellranger_count(self,qc_dir):
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
          qc_dir (str): top-level directory to look under.
        """
        software = {}
        output_files = list()
        cellranger_samples = []
        cellranger_references = set()
        samples_by_pipeline = dict()
        tags = set()
        # Look for cellranger_count outputs
        cellranger_count_dir = os.path.join(qc_dir,
                                            "cellranger_count")
        print("Checking for cellranger* count outputs under %s" %
              cellranger_count_dir)
        cellranger_versioned_samples = {}
        if os.path.isdir(cellranger_count_dir):
            cellranger_name = None
            versions = set()
            # Old-style (unversioned)
            for d in filter(
                    lambda f:
                    os.path.isdir(os.path.join(cellranger_count_dir,f)),
                    os.listdir(cellranger_count_dir)):
                sample_dir = os.path.join(cellranger_count_dir,d)
                try:
                    cellranger = CellrangerCount(sample_dir)
                    output_files.append(cellranger.web_summary)
                    output_files.append(cellranger.metrics_csv)
                    output_files.append(cellranger.cmdline_file)
                    cellranger_samples.append(d)
                    cellranger_name = cellranger.pipeline_name
                    cellranger_references.add(cellranger.reference_data)
                    # Store as version '?'
                    ref = os.path.basename(cellranger.reference_data)
                    if cellranger_name not in cellranger_versioned_samples:
                        cellranger_versioned_samples[cellranger_name] = {}
                    if '?' not in cellranger_versioned_samples[cellranger_name]:
                        cellranger_versioned_samples[cellranger_name]['?'] = {}
                    if ref not in \
                       cellranger_versioned_samples[cellranger_name]['?']:
                        cellranger_versioned_samples[cellranger_name]['?'][ref] = []
                    cellranger_versioned_samples[cellranger_name]['?'][ref].append(d)
                    versions.add('?')
                except OSError:
                    pass
            if cellranger_samples:
                tags.add("%s_count" % cellranger_name)
            # New-style (versioned)
            cellranger_name = None
            for ver in filter(
                    lambda f:
                    os.path.isdir(os.path.join(cellranger_count_dir,f)),
                    os.listdir(cellranger_count_dir)):
                # Check putative version numbers
                for ref in filter(
                        lambda f:
                        os.path.isdir(os.path.join(cellranger_count_dir,ver,f)),
                        os.listdir(os.path.join(cellranger_count_dir,ver))):
                    # Check putative reference dataset names
                    samples = []
                    for smpl in filter(
                            lambda f:
                            os.path.isdir(os.path.join(cellranger_count_dir,
                                                       ver,ref,f)),
                            os.listdir(os.path.join(cellranger_count_dir,
                                                    ver,ref))):
                        sample_dir = os.path.join(cellranger_count_dir,
                                                  ver,ref,smpl)
                        cellranger_name = None
                        try:
                            cellranger = CellrangerCount(sample_dir)
                            output_files.append(cellranger.web_summary)
                            output_files.append(cellranger.metrics_csv)
                            output_files.append(cellranger.cmdline_file)
                            samples.append(smpl)
                            cellranger_name = cellranger.pipeline_name
                            cellranger_references.add(
                                cellranger.reference_data)
                        except OSError:
                            pass
                    # Add outputs, samples and version
                    if samples:
                        tags.add("%s_count" % cellranger_name)
                        if cellranger_name not in cellranger_versioned_samples:
                            cellranger_versioned_samples[cellranger_name] = {}
                        if ver not in cellranger_versioned_samples[cellranger_name]:
                            cellranger_versioned_samples[cellranger_name][ver] = {}
                        cellranger_versioned_samples[cellranger_name][ver][ref] = samples
                        versions.add(ver)
                        for smpl in cellranger_versioned_samples[cellranger_name][ver][ref]:
                            if smpl not in cellranger_samples:
                                cellranger_samples.append(smpl)
            # Store cellranger versions
            for cellranger_name in cellranger_versioned_samples:
                software[cellranger_name] = sorted(list(cellranger_versioned_samples[cellranger_name].keys()))
        # Store sample lists associated with pipeline,
        # version and reference dataset
        for name in cellranger_versioned_samples:
            for version in cellranger_versioned_samples[name]:
                for reference in cellranger_versioned_samples[name][version]:
                    pipeline_key = (name,version,reference)
                    samples_by_pipeline[pipeline_key] = \
                        [s for s in
                         cellranger_versioned_samples[name][version][reference]]
        # Return collected information
        return AttributeDictionary(
            name='cellranger_count',
            software=software,
            references=sorted(list(cellranger_references)),
            fastqs=[],
            samples=cellranger_samples,
            pipelines=sorted([p for p in samples_by_pipeline]),
            samples_by_pipeline=samples_by_pipeline,
            output_files=output_files,
            tags=sorted(list(tags))
        )

    def _collect_cellranger_multi(self,qc_dir):
        """
        Collect information on Cellranger multi outputs

        Returns an AttributeDictionary with the following
        attributes:

        - name: set to 'cellranger_multi'
        - software: dictionary of software and versions
        - references: list of associated reference datasets
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
                            cellranger_references.add(
                                cellranger_multi.reference_data)
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
            fastqs=[],
            multiplexed_samples=sorted(list(multiplexed_samples)),
            pipelines=sorted([p for p in samples_by_pipeline]),
            samples_by_pipeline=samples_by_pipeline,
            output_files=output_files,
            tags=sorted(list(tags))
        )

    def _collect_multiqc(self,qc_dir):
        """
        Collect information on MultQC reports

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

    def _collect_qc_config_files(self,files):
        """
        Collect information on QC config files

        Returns a list of QC configuration (if any) found
        in the supplied list (with the leading paths
        removed).

        The following configuration files are checked for:

        - fastq_strand.conf
        - 10x_multi_config.csv
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
        # Cellranger-arc configs ('libraries.SAMPLE.csv')
        for f in list(filter(lambda f:
                             f.endswith(".csv") and
                             os.path.basename(f).startswith("libraries."),
                             files)):
            config_files.add(os.path.basename(f))
        return sorted(list(config_files))

#######################################################################
# Functions
#######################################################################

def fastq_screen_output(fastq,screen_name,legacy=False):
    """
    Generate name of fastq_screen output files

    Given a Fastq file name and a screen name, the outputs from
    fastq_screen will look like:

    - {FASTQ}_screen_{SCREEN_NAME}.png
    - {FASTQ}_screen_{SCREEN_NAME}.txt

    "Legacy" screen outputs look like:

    - {FASTQ}_{SCREEN_NAME}_screen.png
    - {FASTQ}_{SCREEN_NAME}_screen.txt

    Arguments:
       fastq (str): name of Fastq file
       screen_name (str): name of screen
       legacy (bool): if True then use 'legacy' (old-style)
         naming convention (default: False)

    Returns:
       tuple: fastq_screen output names (without leading path)

    """
    if not legacy:
        base_name = "%s_screen_%s"
    else:
        base_name = "%s_%s_screen"
    base_name = base_name  % (strip_ngs_extensions(os.path.basename(fastq)),
                              str(screen_name))
    return (base_name+'.png',base_name+'.txt')

def fastqc_output(fastq):
    """
    Generate name of FastQC outputs

    Given a Fastq file name, the outputs from FastQC will look
    like:

    - {FASTQ}_fastqc/
    - {FASTQ}_fastqc.html
    - {FASTQ}_fastqc.zip

    Arguments:
       fastq (str): name of Fastq file

    Returns:
       tuple: FastQC outputs (without leading paths)

    """
    base_name = "%s_fastqc" % strip_ngs_extensions(os.path.basename(fastq))
    return (base_name,base_name+'.html',base_name+'.zip')

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

def cellranger_count_output(project,sample_name=None,
                            prefix="cellranger_count"):
    """
    Generate list of 'cellranger count' outputs

    Given an AnalysisProject, the outputs from 'cellranger
    count' will look like:

    - {PREFIX}/{SAMPLE_n}/outs/metrics_summary.csv
    - {PREFIX}/{SAMPLE_n}/outs/web_summary.html

    for each SAMPLE_n in the project.

    If a sample name is supplied then outputs are limited
    to those for that sample

    Arguments:
      project (AnalysisProject): project to generate
        output names for
      sample_name (str): sample to limit outputs to
      prefix (str): directory for outputs (defaults
        to "cellranger_count")

    Returns:
       tuple: cellranger count outputs (without leading paths)
    """
    outputs = []
    # Metrics and web summary files
    for sample in project.samples:
        if sample_name and sample_name != sample.name:
            continue
        sample_count_dir = os.path.join(prefix,
                                        sample.name)
        for f in ("metrics_summary.csv",
                  "web_summary.html"):
            outputs.append(os.path.join(sample_count_dir,
                                        "outs",f))
    return tuple(outputs)

def cellranger_atac_count_output(project,sample_name=None,
                                 prefix="cellranger_count"):
    """
    Generate list of 'cellranger-atac count' outputs

    Given an AnalysisProject, the outputs from 'cellranger-atac
    count' will look like:

    - {PREFIX}/{SAMPLE_n}/outs/summary.csv
    - {PREFIX}/{SAMPLE_n}/outs/web_summary.html

    for each SAMPLE_n in the project.

    If a sample name is supplied then outputs are limited
    to those for that sample

    Arguments:
      project (AnalysisProject): project to generate
        output names for
      sample_name (str): sample to limit outputs to
      prefix (str): directory for outputs (defaults
        to "cellranger_count")

    Returns:
       tuple: cellranger count outputs (without leading paths)
    """
    outputs = []
    # Metrics and web summary files
    for sample in project.samples:
        if sample_name and sample_name != sample.name:
            continue
        sample_count_dir = os.path.join(prefix,
                                        sample.name)
        for f in ("summary.csv",
                  "web_summary.html"):
            outputs.append(os.path.join(sample_count_dir,
                                        "outs",f))
    return tuple(outputs)

def cellranger_arc_count_output(project,sample_name=None,
                                prefix="cellranger_count"):
    """
    Generate list of 'cellranger-arc count' outputs

    Given an AnalysisProject, the outputs from 'cellranger-arc
    count' will look like:

    - {PREFIX}/{SAMPLE_n}/outs/summary.csv
    - {PREFIX}/{SAMPLE_n}/outs/web_summary.html

    for each SAMPLE_n in the project.

    If a sample name is supplied then outputs are limited
    to those for that sample

    Arguments:
      project (AnalysisProject): project to generate
        output names for
      sample_name (str): sample to limit outputs to
      prefix (str): directory for outputs (defaults
        to "cellranger_count")

    Returns:
       tuple: cellranger count outputs (without leading paths)
    """
    outputs = []
    # Metrics and web summary files
    for sample in project.samples:
        if sample_name and sample_name != sample.name:
            continue
        sample_count_dir = os.path.join(prefix,
                                        sample.name)
        for f in ("summary.csv",
                  "web_summary.html"):
            outputs.append(os.path.join(sample_count_dir,
                                        "outs",f))
    return tuple(outputs)

def cellranger_multi_output(project,config_csv,sample_name=None,
                            prefix="cellranger_multi"):
    """
    Generate list of 'cellranger multi' outputs

    Given an AnalysisProject, the outputs from 'cellranger
    multi' will look like:

    - {PREFIX}/outs/multi/multiplexing_analysis/tag_calls_summary.csv

    and

    - {PREFIX}/outs/per_sample_outs/{SAMPLE_n}/metrics_summary.csv
    - {PREFIX}/outs/per_sample_outs/{SAMPLE_n}/web_summary.html

    for each multiplexed SAMPLE_n defined in the config.csv file
    (nb these are not equivalent to the 'samples' defined by the
    Fastq files in the project).

    If a sample name is supplied then outputs are limited
    to those for that sample; if the supplied config.csv file isn't
    found then no outputs will be returned.

    Arguments:
      project (AnalysisProject): project to generate
        output names for
      config_csv (str): path to the cellranger multi
        config.csv file
      sample_name (str): multiplexed sample to limit outputs
        to (optional)
      prefix (str): directory for outputs (optional, defaults
        to "cellranger_multi")

    Returns:
       tuple: cellranger multi outputs (without leading paths)
    """
    outputs = []
    # Check that config.csv file exists
    if not os.path.isfile(config_csv):
        return outputs
    # Per-sample metrics and web summary files
    for sample in CellrangerMultiConfigCsv(config_csv).sample_names:
        if sample_name and sample_name != sample:
            continue
        sample_dir = os.path.join(prefix,
                                  "outs",
                                  "per_sample_outs",
                                  sample)
        for f in ("metrics_summary.csv",
                  "web_summary.html"):
            outputs.append(os.path.join(sample_dir,f))
    # Multiplexing outputs
    multi_analysis_dir = os.path.join(prefix,
                                      "outs",
                                      "multi",
                                      "multiplexing_analysis")
    for f in ("tag_calls_summary.csv",):
        outputs.append(os.path.join(multi_analysis_dir,f))
    return tuple(outputs)

def check_fastq_screen_outputs(project,qc_dir,screen,qc_protocol=None,
                               legacy=False):
    """
    Return Fastqs missing QC outputs from FastqScreen

    Returns a list of the Fastqs from a project for which
    one or more associated outputs from FastqScreen
    don't exist in the specified QC directory.

    Arguments:
      project (AnalysisProject): project to check the
        QC outputs for
      qc_dir (str): path to the QC directory (relative
        path is assumed to be a subdirectory of the
        project)
      screen (str): screen name to check
      qc_protocol (str): QC protocol to predict outputs
        for; if not set then defaults to standard QC
        based on ended-ness
      legacy (bool): if True then check for 'legacy'-style
         names (defult: False)

    Returns:
      List: list of Fastq files with missing outputs.
    """
    if not os.path.isabs(qc_dir):
        qc_dir = os.path.join(project.dirn,qc_dir)
    fastqs = set()
    for fastq in remove_index_fastqs(project.fastqs,
                                     project.fastq_attrs):
        read_numbers = get_read_numbers(qc_protocol).seq_data
        if project.fastq_attrs(fastq).read_number not in read_numbers:
            # Ignore non-data reads
            continue
        for output in [os.path.join(qc_dir,f)
                       for f in fastq_screen_output(fastq,screen,
                                                    legacy=legacy)]:
            if not os.path.exists(output):
                fastqs.add(fastq)
    return sorted(list(fastqs))

def check_fastqc_outputs(project,qc_dir,qc_protocol=None):
    """
    Return Fastqs missing QC outputs from FastQC

    Returns a list of the Fastqs from a project for which
    one or more associated outputs from FastQC don't exist
    in the specified QC directory.

    Arguments:
      project (AnalysisProject): project to check the
        QC outputs for
      qc_dir (str): path to the QC directory (relative
        path is assumed to be a subdirectory of the
        project)
      qc_protocol (str): QC protocol to predict outputs
        for; if not set then defaults to standard QC
        based on ended-ness

    Returns:
      List: list of Fastq files with missing outputs.
    """
    if not os.path.isabs(qc_dir):
        qc_dir = os.path.join(project.dirn,qc_dir)
    fastqs = set()
    for fastq in remove_index_fastqs(project.fastqs,
                                     project.fastq_attrs):
        read_numbers = get_read_numbers(qc_protocol).qc
        if project.fastq_attrs(fastq).read_number not in read_numbers:
            # Ignore non-QC reads
            continue
        # FastQC outputs
        for output in [os.path.join(qc_dir,f)
                       for f in fastqc_output(fastq)]:
            if not os.path.exists(output):
                fastqs.add(fastq)
    return sorted(list(fastqs))

def check_fastq_strand_outputs(project,qc_dir,fastq_strand_conf,
                               qc_protocol=None):
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
      qc_protocol (str): QC protocol to predict outputs
        for; if not set then defaults to standard QC
        based on ended-ness

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
    read_numbers = get_read_numbers(qc_protocol).seq_data
    fastq_pairs = set()
    for fq_group in group_fastqs_by_name(
            remove_index_fastqs(project.fastqs,
                                project.fastq_attrs),
            fastq_attrs=project.fastq_attrs):
        # Assemble Fastq pairs based on read numbers
        fq_pair = tuple([fq_group[r-1] for r in read_numbers])
        # Strand stats output
        output = os.path.join(qc_dir,
                              fastq_strand_output(fq_pair[0]))
        if not os.path.exists(output):
            fastq_pairs.add(fq_pair)
    return sorted(list(fastq_pairs))

def check_cellranger_count_outputs(project,qc_dir=None,
                                   prefix="cellranger_count"):
    """
    Return samples missing QC outputs from 'cellranger count'

    Returns a list of the samples from a project for which
    one or more associated outputs from `cellranger count`
    don't exist in the specified QC directory.

    Arguments:
      project (AnalysisProject): project to check the
        QC outputs for
      qc_dir (str): path to QC directory (if not the default
        QC directory for the project)
      prefix (str): directory for outputs (defaults
        to "cellranger_count")

    Returns:
      List: list of sample names with missing outputs
    """
    if qc_dir is None:
        qc_dir = project.qc_dir
    qc_dir = os.path.abspath(qc_dir)
    samples = set()
    for sample in project.samples:
        for output in cellranger_count_output(project,
                                              sample.name,
                                              prefix):
            if not os.path.exists(os.path.join(qc_dir,output)):
                samples.add(sample.name)
    return sorted(list(samples))

def check_cellranger_atac_count_outputs(project,qc_dir=None,
                                        prefix="cellranger_count"):
    """
    Return samples missing QC outputs from 'cellranger-atac count'

    Returns a list of the samples from a project for which
    one or more associated outputs from `cellranger-atac count`
    don't exist in the specified QC directory.

    Arguments:
      project (AnalysisProject): project to check the
        QC outputs for
      qc_dir (str): path to QC directory (if not the default
        QC directory for the project)
      prefix (str): directory for outputs (defaults
        to "cellranger_count")

    Returns:
      List: list of sample names with missing outputs
    """
    if qc_dir is None:
        qc_dir = project.qc_dir
    qc_dir = os.path.abspath(qc_dir)
    samples = set()
    for sample in project.samples:
        for output in cellranger_atac_count_output(project,
                                                   sample.name,
                                                   prefix):
            if not os.path.exists(os.path.join(qc_dir,output)):
                samples.add(sample.name)
    return sorted(list(samples))

def check_cellranger_arc_count_outputs(project,qc_dir=None,
                                       prefix="cellranger_count"):
    """
    Return samples missing QC outputs from 'cellranger-arc count'

    Returns a list of the samples from a project for which
    one or more associated outputs from `cellranger-arc count`
    don't exist in the specified QC directory.

    Arguments:
      project (AnalysisProject): project to check the
        QC outputs for
      qc_dir (str): path to QC directory (if not the default
        QC directory for the project)
      prefix (str): directory for outputs (defaults
        to "cellranger_count")

    Returns:
      List: list of sample names with missing outputs
    """
    if qc_dir is None:
        qc_dir = project.qc_dir
    qc_dir = os.path.abspath(qc_dir)
    samples = set()
    for sample in project.samples:
        if not os.path.exists(os.path.join(qc_dir,
                                           "libraries.%s.csv"
                                           % sample.name)):
            # Skip if there is no libraries.csv for the sample
            continue
        for output in cellranger_arc_count_output(project,
                                                  sample.name,
                                                  prefix):
            if not os.path.exists(os.path.join(qc_dir,output)):
                samples.add(sample.name)
    return sorted(list(samples))
