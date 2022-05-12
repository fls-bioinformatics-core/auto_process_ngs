#!/usr/bin/env python
#
# rseqc library
import os

"""
Example output from infer_experiment.py:

This is PairEnd Data
Fraction of reads failed to determine: 0.0172
Fraction of reads explained by "1++,1--,2+-,2-+": 0.4903
Fraction of reads explained by "1+-,1-+,2++,2--": 0.4925
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
      if data are single end
    - forward (float): fraction of 'forward' aligned reads
    - reverse (float): fraction of 'reverse' aligned reads
    - unstranded (float): fraction of aligned reads neither
      'forward' nor 'reverse'

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
