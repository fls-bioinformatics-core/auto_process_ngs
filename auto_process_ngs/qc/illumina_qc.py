#!/usr/bin/env python
#
#     illumina_qc: running and validating QC
#     Copyright (C) University of Manchester 2018-2019 Peter Briggs
#
"""
Overview
--------

This module provides utility classes and functions for generating
QC commands and checking the associated outputs.

It provides the following core class:

- IlluminaQC: QC command generator

It also provides the following supporting functions:

- determine_qc_protocol: get QC protocol for a project
- fastq_screen_output: get names for fastq_screen outputs
- fastqc_output: get names for FastQC outputs
- fastq_strand_output: get name for fastq_strand.py output

and the following module-level constants:

- QC_MODULES: tuple of QC module names
- FASTQ_SCREENS: tuple of screen names
- PROTOCOLS: tuple of QC protocol names

Usage
-----

The procedure for using IlluminaQC instances is as follows:

- Create an IlluminaQC instance for the appropriate QC
  protocol and settings;

- Generate a dictionary of Fastq groups with missing QC
  outputs with each QC 'module' under this protocol, by
  invoking the `fastqs_missing_qc` method and supplying
  a list of Fastq files (typically all Fastqs in a
  project);

- Pass this dictionary of Fastq groups into the `commands`
  method, to generate a list of commands that should be
  run to generate the missing QC outputs.

- The QC for a project or other set of Fastqs can also be
  checked by invoking the `fastqs_missing_qc` method and
  examining whether any Fastq groups are returned for any
  of the QC modules.
"""

#######################################################################
# Imports
#######################################################################

import os
import logging
from bcftbx.qc.report import strip_ngs_extensions
from ..applications import Command
from ..analysis import AnalysisFastq
from ..fastq_utils import pair_fastqs_by_name

QC_MODULES = ('illumina_qc',
              'fastq_strand',)

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
                 ungzip_fastqs=False,
                 fastq_attrs=None):
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
          fastq_attrs (BaseFastqAttrs): optional,
            class to use for extracting data from the
            Fastq names (default: AnalysisFastq)
        """
        if protocol not in PROTOCOLS:
            raise Exception("IlluminaQC: unrecognised protocol '%s'" %
                            protocol)
        self.protocol = protocol
        self.fastq_screen_subset = fastq_screen_subset
        self.nthreads = nthreads
        self.fastq_strand_conf = fastq_strand_conf
        self.ungzip_fastqs = ungzip_fastqs
        if fastq_attrs is None:
            self.fastq_attrs = AnalysisFastq
        else:
            self.fastq_attrs = fastq_attrs

    def version(self):
        """
        Return version of QC script
        """
        status,qc_script_info = Command(
            'illumina_qc.sh',
            '--version').subprocess_check_output()
        if status == 0:
            return qc_script_info.strip().split()[-1]

    def commands(self,fastq_groups,qc_dir=None):
        """
        Generate commands for running QC scripts

        Takes a dictionary output from the
        `fastqs_missing_qc` method (where the keys
        are QC module names and the associated values
        are lists of Fastq groupings - either a single
        Fastq or a list or iterable with multiple
        Fastqs - which are missing expected QC
        outputs from that module).

        Returns a list of commands that will generate
        the missing QC outputs when executed.

        Arguments:
          fastq_groups (dictionary): Fastq groupings
            for each QC module which have missing
            associated QC outputs
          qc_dir (str): path to the directory which
            will hold the outputs from the QC script

        Returns:
          List: list of `Command` instances (one per
            Fastq group) for running the QC commands.
        """
        if self.protocol == "standardSE":
            commands = self._commands_for_single_end
        elif self.protocol == "standardPE":
            commands = self._commands_for_paired_end
        elif self.protocol == "singlecell":
            commands = self._commands_for_single_cell
        return commands(fastq_groups,qc_dir=qc_dir)

    def _commands_for_paired_end(self,fastq_groups,
                                 qc_dir=None):
        """
        Generate commands for running paired-end QC scripts

        Takes a dictionary output from the
        `fastqs_missing_qc` method (where the keys
        are QC module names and the associated values
        are lists of Fastq groupings - either a single
        Fastq or a list or iterable with multiple
        Fastqs - which are missing expected QC
        outputs from that module).

        Returns a list of commands that will generate
        the missing QC outputs when executed.

        Arguments:
          fastq_groups (dictionary): Fastq groupings
            for each QC module which have missing
            associated QC outputs
          qc_dir (str): path to the directory which
            will hold the outputs from the QC script

        Returns:
          List: list of `Command` instances (one per
            Fastq group) for running the QC commands.
        """
        cmds = list()
        # Iterate over QC modules
        for qc_module in filter(lambda s: s in fastq_groups,
                                QC_MODULES):
            # Iterate over Fastq groupings in each module
            for fqs in fastq_groups[qc_module]:
                fastqs = self._remove_index_reads(self._to_list(fqs))
                if qc_module == 'illumina_qc':
                    # illumina_qc.sh
                    for fastq in fastqs:
                        cmd = Command('illumina_qc.sh',fastq)
                        if self.ungzip_fastqs:
                            cmd.add_args('--ungzip-fastqs')
                        cmd.add_args('--threads',self.nthreads)
                        if self.fastq_screen_subset is not None:
                            cmd.add_args('--subset',self.fastq_screen_subset)
                        if qc_dir is not None:
                            cmd.add_args('--qc_dir',os.path.abspath(qc_dir))
                        cmds.append(cmd)
                elif qc_module == 'fastq_strand':
                    # fastq_strand.py
                    if self.fastq_strand_conf is not None:
                        for fq_pair in self._fastq_pairs(fastqs):
                            print fq_pair
                            cmd = Command('fastq_strand.py',
                                          '-n',self.nthreads,
                                          '--conf',self.fastq_strand_conf,
                                          '--outdir',os.path.abspath(qc_dir),
                                          *fq_pair)
                            cmds.append(cmd)
        # Return list of commands
        return cmds

    def _commands_for_single_end(self,fastq_groups,
                                 qc_dir=None):
        """
        Generate commands for running single-end QC scripts

        Takes a dictionary output from the
        `fastqs_missing_qc` method (where the keys
        are QC module names and the associated values
        are lists of Fastq groupings - either a single
        Fastq or a list or iterable with multiple
        Fastqs - which are missing expected QC
        outputs from that module).

        Returns a list of commands that will generate
        the missing QC outputs when executed.

        Arguments:
          fastq_groups (dictionary): Fastq groupings
            for each QC module which have missing
            associated QC outputs
          qc_dir (str): path to the directory which
            will hold the outputs from the QC script

        Returns:
          List: list of `Command` instances (one per
            Fastq group) for running the QC commands.
        """
        cmds = list()
        # Iterate over QC modules
        for qc_module in filter(lambda s: s in fastq_groups,
                                QC_MODULES):
            # Iterate over Fastqs in each module
            for fqs in fastq_groups[qc_module]:
                fastqs = self._remove_index_reads(self._to_list(fqs))
                if qc_module == 'illumina_qc':
                    # illumina_qc.sh
                    for fastq in fastqs:
                        cmd = Command('illumina_qc.sh',fastq)
                        if self.ungzip_fastqs:
                            cmd.add_args('--ungzip-fastqs')
                        cmd.add_args('--threads',self.nthreads)
                        if self.fastq_screen_subset is not None:
                            cmd.add_args('--subset',self.fastq_screen_subset)
                        if qc_dir is not None:
                            cmd.add_args('--qc_dir',os.path.abspath(qc_dir))
                        cmds.append(cmd)
                elif qc_module == 'fastq_strand':
                    # fastq_strand.py
                    if self.fastq_strand_conf is not None:
                        for fastq in fastqs:
                            cmd = Command('fastq_strand.py',
                                          '-n',self.nthreads,
                                          '--conf',self.fastq_strand_conf,
                                          '--outdir',os.path.abspath(qc_dir),
                                          fastq)
                            cmds.append(cmd)
        # Return list of commands
        return cmds

    def _commands_for_single_cell(self,fastq_groups,
                                  qc_dir=None):
        """
        Generate commands for running single cell QC scripts

        Takes a dictionary output from the
        `fastqs_missing_qc` method (where the keys
        are QC module names and the associated values
        are lists of Fastq groupings - either a single
        Fastq or a list or iterable with multiple
        Fastqs - which are missing expected QC
        outputs from that module).

        Returns a list of commands that will generate
        the missing QC outputs when executed.

        Arguments:
          fastq_groups (dictionary): Fastq groupings
            for each QC module which have missing
            associated QC outputs
          qc_dir (str): path to the directory which
            will hold the outputs from the QC script

        Returns:
          List: list of `Command` instances (one per
            Fastq group) for running the QC commands.
        """
        cmds = list()
        # Iterate over QC modules
        for qc_module in filter(lambda s: s in fastq_groups,
                                QC_MODULES):
            # Iterate over Fastq groupings in each module
            for fqs in fastq_groups[qc_module]:
                fastqs = self._remove_index_reads(self._to_list(fqs))
                if qc_module == 'illumina_qc':
                    # illumina_qc.sh
                    for fastq in fastqs:
                        cmd = Command('illumina_qc.sh',fastq)
                        if self.ungzip_fastqs:
                            cmd.add_args('--ungzip-fastqs')
                        cmd.add_args('--threads',self.nthreads)
                        if self.fastq_screen_subset is not None:
                            cmd.add_args('--subset',self.fastq_screen_subset)
                        if qc_dir is not None:
                            cmd.add_args('--qc_dir',os.path.abspath(qc_dir))
                        if self.fastq_attrs(fastq).read_number == 1:
                            # Screens for R2 only
                            cmd.add_args('--no-screens')
                        cmds.append(cmd)
                elif qc_module == 'fastq_strand':
                    # fastq_strand.py
                    if self.fastq_strand_conf is not None:
                        for fastq in fastqs:
                            if self.fastq_attrs(fastq).read_number == 2:
                                cmd = Command(
                                    'fastq_strand.py',
                                    '-n',self.nthreads,
                                    '--conf',self.fastq_strand_conf,
                                    '--outdir',os.path.abspath(qc_dir),
                                    fastq)
                                cmds.append(cmd)
        return cmds

    def expected_outputs(self,fastqs,qc_dir):
        """
        Generate expected outputs for input Fastq

        Returns a dictionary where keys are the names
        of QC modules (e.g. 'illumina_qc'), and the
        associated values are lists of tuples, which
        in turn consist of Fastq-output pairs i.e.

        (fastqs,outputs)

        `fastqs` can be a single Fastq file or many
        Fastqs in a list or iterable; `outputs` is a
        list or iterable of the associated QC outputs.

        Arguments:
          fastqs (str/list): either path to a single
            Fastq file, or a list of paths to multiple
            Fastqs
          qc_dir (str): path to the directory which
            will hold the outputs from the QC script

        Returns:
          Dictionary: dictionary of QC modules with
            lists of tuples with Fastqs and associated
            QC output file names.
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

        Returns a dictionary where keys are the names
        of QC modules (e.g. 'illumina_qc'), and the
        associated values are lists of tuples, which
        in turn consist of Fastq-output pairs i.e.

        `(fastqs,outputs)`

        `fastqs` can be a single Fastq file or many
        Fastqs in a list or iterable; `outputs` is a
        list or iterable of the associated QC outputs.

        Arguments:
          fastqs (str/list): either path to a single
            Fastq file, or a list of paths to multiple
            Fastqs
          qc_dir (str): path to the directory which
            will hold the outputs from the QC script

        Returns:
          Dictionary: dictionary of QC modules with
            lists of tuples with Fastqs and associated
            QC output file names.
        """
        qc_dir = os.path.abspath(qc_dir)
        fastqs = self._remove_index_reads(self._to_list(fastqs))
        outputs = dict()
        for qc_module in QC_MODULES:
            output = []
            if qc_module == 'illumina_qc':
                # Expected outputs for single Fastqs
                for fastq in fastqs:
                    expected = []
                    # FastQC outputs
                    expected.extend([os.path.join(qc_dir,f)
                                     for f in fastqc_output(fastq)])
                    # Fastq_screen outputs
                    for name in FASTQ_SCREENS:
                        expected.extend([os.path.join(qc_dir,f)
                                         for f in
                                         fastq_screen_output(fastq,name)])
                    output.append((fastq,expected))
            elif qc_module == 'fastq_strand':
                # Pair-wise outputs
                if self.fastq_strand_conf:
                    for fq_pair in self._fastq_pairs(fastqs):
                        # Strand stats output
                        expected = [os.path.join(
                            qc_dir,fastq_strand_output(fq_pair[0]))]
                        output.append((fq_pair,expected))
            else:
                logger.warning("Ignoring unimplemented QC module '%s'"
                               % qc_module)
            if output:
                outputs[qc_module] = output
        return outputs

    def _outputs_for_single_end(self,fastqs,qc_dir):
        """
        Generate expected outputs for standard SE protocol

        Returns a dictionary where keys are the names
        of QC modules (e.g. 'illumina_qc'), and the
        associated values are lists of tuples, which
        in turn consist of Fastq-output pairs i.e.

        (fastqs,outputs)

        `fastqs` can be a single Fastq file or many
        Fastqs in a list or iterable; `outputs` is a
        list or iterable of the associated QC outputs.

        Arguments:
          fastqs (str/list): either path to a single
            Fastq file, or a list of paths to multiple
            Fastqs
          qc_dir (str): path to the directory which
            will hold the outputs from the QC script

        Returns:
          Dictionary: dictionary of QC modules with
            lists of tuples with Fastqs and associated
            QC output file names.
        """
        qc_dir = os.path.abspath(qc_dir)
        fastqs = self._remove_index_reads(self._to_list(fastqs))
        outputs = dict()
        for qc_module in QC_MODULES:
            output = []
            if qc_module == 'illumina_qc':
                # Expected outputs for single Fastqs
                for fastq in fastqs:
                    expected = []
                    # FastQC outputs
                    expected.extend([os.path.join(qc_dir,f)
                                     for f in fastqc_output(fastq)])
                    # Fastq_screen outputs
                    for name in FASTQ_SCREENS:
                        expected.extend([os.path.join(qc_dir,f)
                                         for f in
                                         fastq_screen_output(fastq,name)])
                    output.append((fastq,expected))
            elif qc_module == 'fastq_strand':
                # Strand stats output
                if self.fastq_strand_conf:
                    expected = [os.path.join(
                        qc_dir,fastq_strand_output(fastq))]
                    output.append((fastq,expected))
            else:
                logger.warning("Ignoring unimplemented QC module '%s'"
                               % qc_module)
            if output:
                outputs[qc_module] = output
        return outputs

    def _outputs_for_single_cell(self,fastqs,qc_dir):
        """
        Generate expected outputs for single-cell protocol

        Returns a dictionary where keys are the names
        of QC modules (e.g. 'illumina_qc'), and the
        associated values are lists of tuples, which
        in turn consist of Fastq-output pairs i.e.

        (fastqs,outputs)

        `fastqs` can be a single Fastq file or many
        Fastqs in a list or iterable; `outputs` is a
        list or iterable of the associated QC outputs.

        Arguments:
          fastqs (str/list): either path to a single
            Fastq file, or a list of paths to multiple
            Fastqs
          qc_dir (str): path to the directory which
            will hold the outputs from the QC script

        Returns:
          Dictionary: dictionary of QC modules with
            lists of tuples with Fastqs and associated
            QC output file names.
        """
        qc_dir = os.path.abspath(qc_dir)
        fastqs = self._remove_index_reads(self._to_list(fastqs))
        outputs = dict()
        for qc_module in QC_MODULES:
            output = []
            if qc_module == 'illumina_qc':
                # Expected outputs for single Fastqs
                for fastq in fastqs:
                    expected = []
                    # FastQC outputs
                    expected.extend([os.path.join(qc_dir,f)
                                     for f in fastqc_output(fastq)])
                    # Fastq_screen outputs (R2 only)
                    if self.fastq_attrs(fastq).read_number == 2:
                        for name in FASTQ_SCREENS:
                            expected.extend([os.path.join(qc_dir,f)
                                             for f in
                                             fastq_screen_output(fastq,name)])
                    output.append((fastq,expected))
            elif qc_module == 'fastq_strand':
                # Strand stats output (R2 only)
                if self.fastq_strand_conf:
                    for fastq in fastqs:
                        if self.fastq_attrs(fastq).read_number == 2:
                            expected = [os.path.join(
                                qc_dir,fastq_strand_output(fastq))]
                            output.append((fastq,expected))
            else:
                logger.warning("Ignoring unimplemented QC module '%s'"
                               % qc_module)
            if output:
                outputs[qc_module] = output
        return outputs

    def fastqs_missing_qc(self,fastqs,qc_dir):
        """
        Filter Fastqs with missing QC outputs

        Returns a dictionary where the keys are
        QC modules and the associated values are
        lists of Fastq groupings (either a single
        Fastq or a list or iterable with multiple
        Fastqs) which are missing expected QC
        outputs from that module.

        The output from this method is the correct
        format for input into the `commands`
        method (to generate commands that will
        produce the missing outputs when executed).

        Arguments:
          fastqs (str/list): either path to a single
            Fastq file, or a list of paths to multiple
            Fastqs
          qc_dir (str): path to the directory which
            will hold the outputs from the QC script

        Returns:
          Dictionary: dictionary where keys are
            QC modules and values are lists of
            Fastq groupings which are missing
            expected QC outputs.
        """
        result = dict()
        qc_dir = os.path.abspath(qc_dir)
        outputs = self.check_outputs(fastqs,qc_dir)
        for qc_module in outputs:
            unverified = list()
            for fqs,present,missing in outputs[qc_module]:
                if missing:
                    unverified.append(fqs)
            if unverified:
                result[qc_module] = unverified
        return result

    def check_outputs(self,fastqs,qc_dir):
        """
        Check QC outputs for input Fastq

        Returns a dictionary with QC modules as
        keys and the corresponding values being
        lists of (fastqs,present,missing) tuples
        where:

        - `fastqs` are one or more Fastqs
        - `present` is a list of the output files
           that are found on the system
        - `missing` is a list of outputs files
           that aren't found.

        Arguments:
          fastqs (str/list): either path to a single
            Fastq file, or a list of paths to multiple
            Fastqs
          qc_dir (str): path to the directory which
            will hold the outputs from the QC script

        Returns:
          Dictionary: dictionary with QC modules as
            keys and corresponding values being
            lists of (fastqs,present,missing) tuples.
        """
        qc_dir = os.path.abspath(qc_dir)
        result = dict()
        # Check whether outputs exist
        expected = self.expected_outputs(fastqs,qc_dir)
        for qc_module in expected:
            for fqs,outputs in expected[qc_module]:
                present = []
                missing = []
                for output in outputs:
                    if os.path.exists(output):
                        present.append(output)
                    else:
                        missing.append(output)
                try:
                    result[qc_module].append((fqs,present,missing))
                except KeyError:
                    result[qc_module] = [(fqs,present,missing),]
        return result

    def _remove_index_reads(self,fastqs):
        """
        Internal: remove index (I1/I2) Fastqs from list
        """
        return filter(lambda fq:
                      not self.fastq_attrs(fq).is_index_read,
                      fastqs)

    def _fastq_pairs(self,fastqs):
        """
        Internal: remove singletons and return paired Fastqs
        """
        return filter(lambda p: len(p) == 2,
                      pair_fastqs_by_name(fastqs,
                                          fastq_attrs=self.fastq_attrs))

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
