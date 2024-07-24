#!/usr/bin/env python3
#
#     fastqc: implements 'sequence_lengths' QC module
#     Copyright (C) University of Manchester 2024 Peter Briggs

"""
Implements the 'sequence_lengths' QC module:

* SequenceLengths: core QCModule class
* GetSeqLengthStats: pipeline task to generate sequence length data
"""

#######################################################################
# Imports
#######################################################################

import os
import logging
from bcftbx.utils import AttributeDictionary
from . import QCModule
from ..seqlens import SeqLens
from ..seqlens import get_sequence_lengths
from ..utils import filter_fastqs
from ...fastq_utils import remove_index_fastqs
from ...pipeliner import PipelineFunctionTask

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Core class
#######################################################################

class SequenceLengths(QCModule):
    """
    Class for handling the 'sequence_lengths' QC module
    """
    name = "sequence_lengths"
    mapped_metrics = False
    require_bam_files = False
    runners = ("fastqc_runner",)
    
    def __init__(self):
        QCModule.__init__(self)

    @classmethod
    def collect_qc_outputs(self,qc_dir):
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
          qc_dir (QCDir): QC directory to examine
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
                               qc_dir.file_list))
        logger.debug("seq_lens: %s" % seq_lens)
        print("\t- %d sequence length files" % len(seq_lens))
        if seq_lens:
            tags.add("sequence_lengths")
            for f in seq_lens:
                fq = qc_dir.fastq_attrs(f[:-len("_seqlens.json")])
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

    @classmethod
    def verify(self,params,qc_outputs):
        """
        Verify 'sequence_lengths' QC module against outputs

        Returns one of 3 values:

        - True: outputs verified ok
        - False: outputs failed to verify
        - None: verification not possible

        Arguments:
          params (AttributeDictionary): values of parameters
            used as inputs
          qc_outputs (AttributeDictionary): QC outputs returned
            from the 'collect_qc_outputs' method
        """
        if not params.fastqs:
            # Nothing to check
            return None
        try:
            # Filter Fastq names
            fastqs = filter_fastqs(params.qc_reads,params.fastqs)
            # Check that outputs exist for every Fastq
            for fq in fastqs:
                if fq not in qc_outputs.fastqs:
                    return False
            return True
        except KeyError:
            # No sequence length outputs present
            return False

    @classmethod
    def add_to_pipeline(self,p,project_name,project,qc_dir,
                        read_numbers,fastqs,require_tasks=[],
                        compute_runner=None):
        """
        Adds tasks for 'sequence_lengths' module to pipeline

        Arguments:
          p (Pipeline): pipeline to extend
          project_name (str): name of project
          project (AnalysisProject): project to run module on
          qc_dir (str): path to QC directory
          read_numbers (list): read numbers to include
          fastqs (list): Fastqs to run the module on
          require_tasks (list): list of tasks that the module
            needs to wait for
          compute_runner (JobRunner): runner to use for
            computation
        """
        get_seq_lengths = GetSeqLengthStats(
            "%s: get sequence length statistics" %
            project_name,
            project,
            qc_dir,
            read_numbers=read_numbers,
            fastqs=fastqs,
            fastq_attrs=project.fastq_attrs)
        p.add_task(get_seq_lengths,
                   requires=require_tasks,
                   runner=compute_runner)
        return get_seq_lengths

#######################################################################
# Pipeline tasks
#######################################################################

class GetSeqLengthStats(PipelineFunctionTask):
    """
    Get data on sequence lengths, masking and padding
    for Fastqs in a project, and write the data to
    JSON files.
    """
    def init(self,project,qc_dir,read_numbers=None,fastqs=None,
             fastq_attrs=None):
        """
        Initialise the GetSeqLengthStats task

        Arguments:
          project (AnalysisProject): project with Fastqs
            to get the sequence length data from
          qc_dir (str): directory for QC outputs (defaults
            to subdirectory 'qc' of project directory)
          read_numbers (sequence): list of read numbers to
            include (or None to include all reads)
          fastqs (list): optional, list of Fastq files
            (overrides Fastqs in project)
          fastq_attrs (BaseFastqAttrs): class to use for
            extracting data from Fastq names
        """
        self._fastqs = list()
    def setup(self):
        # Input Fastqs
        if self.args.fastqs:
            fastqs_in = self.args.fastqs
        else:
            fastqs_in = self.args.project.fastqs
        # Remove index Fastqs
        self._fastqs = remove_index_fastqs(
            fastqs_in,
            fastq_attrs=self.args.fastq_attrs)
        # Get sequence length data for Fastqs
        for fastq in self._fastqs:
            if self.args.read_numbers and \
               self.args.fastq_attrs(fastq).read_number \
               not in self.args.read_numbers:
                continue
            outfile = os.path.join(self.args.qc_dir,
                                   "%s_seqlens.json" %
                                   self.args.fastq_attrs(fastq))
            if os.path.exists(outfile):
                continue
            self.add_call(
                "Get read lengths for %s" % os.path.basename(fastq),
                get_sequence_lengths,
                fastq,
                outfile=outfile)
    def finish(self):
        for result in self.result():
            # Fastq name
            fastq = result['fastq']
            # Check output file exists
            outfile = os.path.join(self.args.qc_dir,
                                   "%s_seqlens.json" %
                                   self.args.fastq_attrs(fastq))
            if not os.path.exists(outfile):
                raise Exception("Missing sequence length file: %s"  %
                                outfile)
