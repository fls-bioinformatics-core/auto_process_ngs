#!/usr/bin/env python
#
#     illumina_qc: running and validating QC
#     Copyright (C) University of Manchester 2018 Peter Briggs
#
"""
illumina_qc.py
==============

Provides utility classes and functions for generating QC commands
and checking outputs from them.

Provides the following core class:

- IlluminaQC: QC command generator

Provides the following supporting functions:

- determine_qc_protocol: get QC protocol for a project
- fastq_screen_output: get names for fastq_screen outputs
- fastqc_output: get names for FastQC outputs
- fastq_strand_output: get name for fastq_strand.py output

Provides the following constants:

- FASTQ_SCREENS: tuple of screen names
- PROTOCOLS: tuple of QC protocol names
"""

#######################################################################
# Imports
#######################################################################

import os
import logging
from bcftbx.qc.report import strip_ngs_extensions
from ..applications import Command
from ..fastq_utils import IlluminaFastqAttrs
from ..fastq_utils import pair_fastqs_by_name

FASTQ_SCREENS = ('model_organisms',
                 'other_organisms',
                 'rRNA',)

PROTOCOLS = ('standardPE',
             'standardSE',
             'singlecell')

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Classes
#######################################################################

class IlluminaQC(object):
    """
    Utility class for running 'illumina_qc.sh'
    """
    def __init__(self,protocol="standardPE",
                 fastq_screen_subset=None,nthreads=1,
                 fastq_strand_conf=None,
                 ungzip_fastqs=False):
        """
        Create a new IlluminaQC instance

        Arguments:
          protocol (str): QC protocol, must one of
            those listed in PROTOCOLS. If 'None' then
            defaults to "standardPE"
          fastq_screen_subset (int): subset of reads
            to use when running Fastq_screen ('None'
            uses the script default)
          nthreads (int): number of cores (threads)
            to run the QC using (default: 1)
          fastq_strand_conf (str): path to a config
            file with STAR indexes to use for strand
            determination
          ungzip_fastqs (bool): if True then also
            ungzip the source Fastqs (if gzipped)
            (default is not to uncompress the Fastqs)
        """
        if protocol not in PROTOCOLS:
            raise Exception("IlluminaQC: unrecognised protocol '%s'" %
                            protocol)
        self.protocol = protocol
        self.fastq_screen_subset = fastq_screen_subset
        self.nthreads = nthreads
        self.fastq_strand_conf = fastq_strand_conf
        self.ungzip_fastqs = ungzip_fastqs

    def version(self):
        """
        Return version of QC script
        """
        status,qc_script_info = Command(
            'illumina_qc.sh',
            '--version').subprocess_check_output()
        if status == 0:
            return qc_script_info.strip().split()[-1]

    def commands(self,fastqs,qc_dir=None):
        """
        Generate commands for running QC scripts

        Note that index reads (e.g. I1 Fastqs) will
        not have commands generated for them.

        Arguments:
          fastqs (list): list of paths to Fastq files
            to run the QC script on
          qc_dir (str): path to the directory which
            will hold the outputs from the QC script

        Returns:
          List: list of `Command` instances (one per
            Fastq) for running the `illumina_qc.sh`
            script on.
        """
        if self.protocol == "standardSE":
            commands = self._commands_for_single_end
        elif self.protocol == "standardPE":
            commands = self._commands_for_paired_end
        elif self.protocol == "singlecell":
            commands = self._commands_for_single_cell
        return commands(fastqs,qc_dir=qc_dir)

    def _commands_for_paired_end(self,fastqs,qc_dir=None):
        """
        Generate commands for running paired-end QC scripts

        Note that index reads (e.g. I1 Fastqs) will
        not have commands generated for them.

        Arguments:
          fastqs (list): list of paths to Fastq files
            to run the QC script on
          qc_dir (str): path to the directory which
            will hold the outputs from the QC script

        Returns:
          List: list of `Command` instances (one per
            Fastq) for running the `illumina_qc.sh`
            script on.
        """
        cmds = list()
        # Convert to list and filter out index reads
        fastqs = self._remove_index_reads(self._to_list(fastqs))
        # Generate QC commands for individual Fastqs
        for fastq in fastqs:
            # Build command
            cmd = Command('illumina_qc.sh',fastq)
            if self.ungzip_fastqs:
                cmd.add_args('--ungzip-fastqs')
            cmd.add_args('--threads',self.nthreads)
            if self.fastq_screen_subset is not None:
                cmd.add_args('--subset',self.fastq_screen_subset)
            if qc_dir is not None:
                cmd.add_args('--qc_dir',os.path.abspath(qc_dir))
            cmds.append(cmd)
        # Generate pair-wise QC commands
        for fq_pair in self._fastq_pairs(fastqs):
            # Strandedness for this pair
            if self.fastq_strand_conf is not None:
                cmd = Command('fastq_strand.py',
                              '-n',self.nthreads,
                              '--conf',self.fastq_strand_conf,
                              '--outdir',os.path.abspath(qc_dir),
                              *fq_pair)
                cmds.append(cmd)
        return cmds

    def _commands_for_single_end(self,fastqs,qc_dir=None):
        """
        Generate commands for running single-end QC scripts

        Note that index reads (e.g. I1 Fastqs) will
        not have commands generated for them.

        Arguments:
          fastqs (list): list of paths to Fastq files
            to run the QC script on
          qc_dir (str): path to the directory which
            will hold the outputs from the QC script

        Returns:
          List: list of `Command` instances (one per
            Fastq) for running the `illumina_qc.sh`
            script on.
        """
        cmds = list()
        # Convert to list and filter out index reads
        fastqs = self._remove_index_reads(self._to_list(fastqs))
        # Generate QC commands for individual Fastqs
        for fastq in fastqs:
            # Build command
            cmd = Command('illumina_qc.sh',fastq)
            if self.ungzip_fastqs:
                cmd.add_args('--ungzip-fastqs')
            cmd.add_args('--threads',self.nthreads)
            if self.fastq_screen_subset is not None:
                cmd.add_args('--subset',self.fastq_screen_subset)
            if qc_dir is not None:
                cmd.add_args('--qc_dir',os.path.abspath(qc_dir))
            cmds.append(cmd)
            # Strandedness
            if self.fastq_strand_conf is not None:
                cmd = Command('fastq_strand.py',
                              '-n',self.nthreads,
                              '--conf',self.fastq_strand_conf,
                              '--outdir',os.path.abspath(qc_dir),
                              fastq)
                cmds.append(cmd)
        return cmds

    def _commands_for_single_cell(self,fastqs,qc_dir=None):
        """
        Generate commands for running single cell QC scripts

        Note that index reads (e.g. I1 Fastqs) will
        not have commands generated for them.

        Arguments:
          fastqs (list): list of paths to Fastq files
            to run the QC script on
          qc_dir (str): path to the directory which
            will hold the outputs from the QC script

        Returns:
          List: list of `Command` instances (one per
            Fastq) for running the `illumina_qc.sh`
            script on.
        """
        cmds = list()
        # Convert to list and filter out index reads
        fastqs = self._remove_index_reads(self._to_list(fastqs))
        # Generate QC commands for individual Fastqs
        for fastq in fastqs:
            # Build command
            cmd = Command('illumina_qc.sh',fastq)
            if self.ungzip_fastqs:
                cmd.add_args('--ungzip-fastqs')
            cmd.add_args('--threads',self.nthreads)
            if self.fastq_screen_subset is not None:
                cmd.add_args('--subset',self.fastq_screen_subset)
            if qc_dir is not None:
                cmd.add_args('--qc_dir',os.path.abspath(qc_dir))
            cmds.append(cmd)
            # Strandedness (R2 only)
            if self.fastq_strand_conf is not None:
                if IlluminaFastqAttrs(fastq).read_number == 2:
                    cmd = Command('fastq_strand.py',
                                  '-n',self.nthreads,
                                  '--conf',self.fastq_strand_conf,
                                  '--outdir',os.path.abspath(qc_dir),
                                  fastq)
                    cmds.append(cmd)
        return cmds

    def expected_outputs(self,fastqs,qc_dir):
        """
        Generate expected outputs for input Fastq

        Arguments:
          fastqs (str/list): either path to a single
            Fastq file, or a list of paths to multiple
            Fastqs
          qc_dir (str): path to the directory which
            will hold the outputs from the QC script

        Returns:
          List: list of expected output files from
            the QC for the supplied Fastq.
        """
        if self.protocol == "standardSE":
            outputs = self._outputs_for_single_end
        elif self.protocol == "standardPE":
            outputs = self._outputs_for_paired_end
        elif self.protocol == "singlecell":
            outputs = self._outputs_for_single_cell
        return outputs(fastqs,qc_dir)

    def _outputs_for_paired_end(self,fastqs,qc_dir):
        """
        Generate expected outputs for standard PE protocol

        Arguments:
          fastqs (str/list): either path to a single
            Fastq file, or a list of paths to multiple
            Fastqs
          qc_dir (str): path to the directory which
            will hold the outputs from the QC script

        Returns:
          List: list of expected output files from
            the QC for the supplied Fastq.
        """
        qc_dir = os.path.abspath(qc_dir)
        fastqs = self._remove_index_reads(self._to_list(fastqs))
        expected = []
        # Expected outputs for single Fastqs
        for fastq in fastqs:
            # FastQC outputs
            expected.extend([os.path.join(qc_dir,f)
                             for f in fastqc_output(fastq)])
            # Fastq_screen outputs
            for name in FASTQ_SCREENS:
                expected.extend([os.path.join(qc_dir,f)
                                 for f in fastq_screen_output(fastq,name)])
        # Pair-wise outputs
        for fq_pair in self._fastq_pairs(fastqs):
            # Strand stats output
            if self.fastq_strand_conf:
                expected.append(os.path.join(
                    qc_dir,fastq_strand_output(fq_pair[0])))
        return expected

    def _outputs_for_single_end(self,fastqs,qc_dir):
        """
        Generate expected outputs for standard SE protocol

        Arguments:
          fastqs (str/list): either path to a single
            Fastq file, or a list of paths to multiple
            Fastqs
          qc_dir (str): path to the directory which
            will hold the outputs from the QC script

        Returns:
          List: list of expected output files from
            the QC for the supplied Fastq.
        """
        qc_dir = os.path.abspath(qc_dir)
        fastqs = self._remove_index_reads(self._to_list(fastqs))
        expected = []
        # Expected outputs for single Fastqs
        for fastq in fastqs:
            # FastQC outputs
            expected.extend([os.path.join(qc_dir,f)
                             for f in fastqc_output(fastq)])
            # Fastq_screen outputs
            for name in FASTQ_SCREENS:
                expected.extend([os.path.join(qc_dir,f)
                                 for f in fastq_screen_output(fastq,name)])
            # Strand stats output
            if self.fastq_strand_conf:
                expected.append(os.path.join(
                    qc_dir,fastq_strand_output(fastq)))
        return expected

    def _outputs_for_single_cell(self,fastqs,qc_dir):
        """
        Generate expected outputs for single-cell protocol

        Arguments:
          fastqs (str/list): either path to a single
            Fastq file, or a list of paths to multiple
            Fastqs
          qc_dir (str): path to the directory which
            will hold the outputs from the QC script

        Returns:
          List: list of expected output files from
            the QC for the supplied Fastq.
        """
        qc_dir = os.path.abspath(qc_dir)
        fastqs = self._remove_index_reads(self._to_list(fastqs))
        expected = []
        # Expected outputs for single Fastqs
        for fastq in fastqs:
            # FastQC outputs
            expected.extend([os.path.join(qc_dir,f)
                             for f in fastqc_output(fastq)])
            # Fastq_screen outputs
            for name in FASTQ_SCREENS:
                expected.extend([os.path.join(qc_dir,f)
                                 for f in fastq_screen_output(fastq,name)])
            # Strand stats output (R2 only)
            if self.fastq_strand_conf:
                if IlluminaFastqAttrs(fastq).read_number == 2:
                    expected.append(os.path.join(
                        qc_dir,fastq_strand_output(fastq)))
        return expected

    def check_outputs(self,fastqs,qc_dir):
        """
        Check QC outputs for input Fastq

        Arguments:
          fastqs (str/list): either path to a single
            Fastq file, or a list of paths to multiple
            Fastqs
          qc_dir (str): path to the directory which
            will hold the outputs from the QC script

        Returns:
          Tuple: tuple (present,missing), where
            'present' is a list of outputs which were
            found, and 'missing' is a list of those
            which were not.
        """
        qc_dir = os.path.abspath(qc_dir)
        present = []
        missing = []
        # Check that outputs exist
        for output in self.expected_outputs(fastqs,qc_dir):
            if os.path.exists(output):
                present.append(output)
            else:
                missing.append(output)
        return (present,missing)

    def _remove_index_reads(self,fastqs):
        """
        Internal: remove index (I1/I2) Fastqs from list
        """
        return filter(lambda fq:
                      not IlluminaFastqAttrs(fq).is_index_read,
                      fastqs)

    def _fastq_pairs(self,fastqs):
        """
        Internal: remove singletons and return paired Fastqs
        """
        return filter(lambda p: len(p) == 2,
                      pair_fastqs_by_name(fastqs))

    def _to_list(self,args):
        """
        Internal: convert arg to appropriate list
        """
        try:
            if len(args[0]) == 1:
                return (args,)
        except IndexError:
            pass
        return args

#######################################################################
# Functions
#######################################################################

def determine_qc_protocol(project):
    """
    Determine the QC protocol for a project

    Arguments:
      project (AnalysisProject): project instance

    Return:
      String: QC protocol for the project
    """
    if project.info.paired_end:
        protocol = "standardPE"
    else:
        protocol = "standardSE"
    if project.info.single_cell_platform is not None:
        protocol = "singlecell"
    return protocol

def fastq_screen_output(fastq,screen_name):
    """
    Generate name of fastq_screen output files

    Given a Fastq file name and a screen name, the outputs from
    fastq_screen will look like:

    - {FASTQ}_{SCREEN_NAME}_screen.png
    - {FASTQ}_{SCREEN_NAME}_screen.txt

    Arguments:
       fastq (str): name of Fastq file
       screen_name (str): name of screen

    Returns:
       tuple: fastq_screen output names (without leading path)

    """
    base_name = "%s_%s_screen" % (strip_ngs_extensions(os.path.basename(fastq)),
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
