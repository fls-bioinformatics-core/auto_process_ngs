#!/usr/bin/env python
#
#     illumina_qc: running and validating QC
#     Copyright (C) University of Manchester 2018 Peter Briggs
#

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

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Classes
#######################################################################

class IlluminaQC(object):
    """
    Utility class for running 'illumina_qc.sh'
    """
    def __init__(self,fastq_screen_subset=None,nthreads=1,
                 fastq_strand_conf=None,ungzip_fastqs=False):
        """
        Create a new IlluminaQC instance

        Arguments:
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
        Generate commands for running QC script

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
        # Filter out index reads
        fastqs = filter(lambda fq:
                        not IlluminaFastqAttrs(fq).is_index_read,
                        fastqs)
        # Generate QC commands for individual Fastqs
        for fastq in fastqs:
            # Skip index reads (i.e. I1)
            if IlluminaFastqAttrs(fastq).is_index_read:
                continue
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
        for fq_pair in pair_fastqs_by_name(fastqs):
            if len(fq_pair) != 2:
                continue
            # Strandedness for this pair
            if self.fastq_strand_conf is not None:
                cmd = Command('fastq_strand.py',
                              '-n',self.nthreads,
                              '--conf',self.fastq_strand_conf,
                              '--outdir',os.path.abspath(qc_dir),
                              *fq_pair)
                cmds.append(cmd)
        return cmds

    def expected_outputs(self,fastq,qc_dir):
        """
        Generate expected outputs for input Fastq

        Arguments
          fastq (str): path to a Fastq file
          qc_dir (str): path to the directory which
            will hold the outputs from the QC script

        Returns:
          List: list of expected output files from
            the QC for the supplied Fastq.
        """
        qc_dir = os.path.abspath(qc_dir)
        expected = []
        # Skip index reads (i.e. I1)
        if IlluminaFastqAttrs(fastq).is_index_read:
            return expected
        # FastQC outputs
        expected.extend([os.path.join(qc_dir,f)
                         for f in fastqc_output(fastq)])
        # Fastq_screen outputs
        for name in FASTQ_SCREENS:
            expected.extend([os.path.join(qc_dir,f)
                             for f in fastq_screen_output(fastq,name)])
        return expected

    def check_outputs(self,fastq,qc_dir):
        """
        Check QC outputs for input Fastq

        Arguments:
          fastq (str): path to a Fastq file
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
        for output in self.expected_outputs(fastq,qc_dir):
            if os.path.exists(output):
                present.append(output)
            else:
                missing.append(output)
        return (present,missing)

#######################################################################
# Functions
#######################################################################

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
