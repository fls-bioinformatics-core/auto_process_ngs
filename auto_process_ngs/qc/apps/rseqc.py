#!/usr/bin/env python
#
#     qc/rseqc: utilities for handling RSeQC outputs
#     Copyright (C) University of Manchester 2024-2025 Peter Briggs
#
"""
Provides utility classes and functions for handling RSeQC outputs.

Provides the following classes:

- InferExperiment: wrapper for handling outputs from 'infer_experiment.py'

Provides the following functions:

- fastqc_output_files: generates names of FastQC outputs files
"""

#######################################################################
# Imports
#######################################################################

import os
import logging

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Classes
#######################################################################

"""
Example output from infer_experiment.py:

> This is PairEnd Data
> Fraction of reads failed to determine: 0.0172
> Fraction of reads explained by "1++,1--,2+-,2-+": 0.4903
> Fraction of reads explained by "1+-,1-+,2++,2--": 0.4925

Example output for poor data:

> Unknown Data type
"""

class InferExperiment(object):
    """
    Wrapper class for handling outputs from infer_experiment.py

    The ``InferExperiment`` object gives access to various
    aspects of the outputs of the RSeQC ``infer_experiment.py``
    utility.

    The following properties are available:

    - log_file (str): path to the source log file
    - paired_end (bool): True if data are paired, False
      if data are single end, None for unknown data type
    - forward (float): fraction of 'forward' aligned reads
    - reverse (float): fraction of 'reverse' aligned reads
    - unstranded (float): fraction of aligned reads neither
      'forward' nor 'reverse'
    - known_data_type (bool): True if data type could be
      identified, False if not
    """
    def __init__(self,infer_experiment_log):
        """
        Create a new InferExperiment instance

        Arguments:
          infer_experiment_log (str): path to the log file
            from a run of ``infer_experiment.py``
        """
        # Initialise
        self._infer_experiment_log = os.path.abspath(infer_experiment_log)
        self._paired_end = None
        self._unstranded = None
        self._forward = None
        self._reverse = None
        self._unknown = False
        # Process the log contents
        with open(self._infer_experiment_log,'rt') as fp:
            for line in fp:
                line = line.strip()
                if line == "This is PairEnd Data":
                    self._paired_end = True
                elif line == "This is SingleEnd Data":
                    self._paired_end = False
                elif line.startswith("Fraction of reads failed to determine:"):
                    self._unstranded = float(line.split()[-1])
                elif line.startswith("Fraction of reads explained by \"1++,1--,2+-,2-+\":") or \
                     line.startswith("Fraction of reads explained by \"++,--\":"):
                    self._forward = float(line.split()[-1])
                elif line.startswith("Fraction of reads explained by \"1+-,1-+,2++,2--\":") or \
                     line.startswith("Fraction of reads explained by \"+-,-+\":"):
                    self._reverse = float(line.split()[-1])
                elif line == "Unknown Data type":
                    self._unknown = True
                    break
                else:
                    pass

    @property
    def log_file(self):
        """
        Path to source log file
        """
        return self._infer_experiment_log

    @property
    def paired_end(self):
        """
        True if data are paired end, False if single end
        """
        return self._paired_end

    @property
    def forward(self):
        """
        Fraction of 'forward' aligned reads
        """
        return self._forward

    @property
    def reverse(self):
        """
        Fraction of 'reverse' aligned reads
        """
        return self._reverse

    @property
    def unstranded(self):
        """
        Fraction of aligned reads neither 'forward' nor 'reverse'
        """
        return self._unstranded

    @property
    def unknown(self):
        """
        True if data type could be identified, False if not
        """
        return self._unknown

#######################################################################
# Functions
#######################################################################

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
