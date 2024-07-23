#!/usr/bin/env python3
#
#     multiqc: implements 'multiqc' QC module
#     Copyright (C) University of Manchester 2024 Peter Briggs

"""
Implements the 'multiqc' QC module:

* Multiqc: core QCModule class
"""

#######################################################################
# Imports
#######################################################################

import os
import logging
from bcftbx.utils import AttributeDictionary
from . import QCModule

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Core class
#######################################################################

class Multiqc(QCModule):
    """
    Class for handling the 'multiqc' QC module
    """
    name = "multiqc"
    mapped_metrics = False
    runners = ()
    envmodules = ()

    def __init__(self):
        QCModule.__init__(self)

    @classmethod
    def collect_qc_outputs(self,qc_dir):
        """
        Collect information multiqc outputs

        Returns an AttributeDictionary with the following
        attributes:

        - name: set to 'multiqc'
        - software: dictionary of software and versions
        - fastqs: list of associated Fastq names
        - output_files: list of associated output files
        - tags: list of associated output classes

        Arguments:
          qc_dir (QCDir): QC directory to examine
        """
        version = None
        output_files = list()
        tags = set()
        # Look for MultiQC report
        multiqc_dir = os.path.dirname(qc_dir.path)
        print("Checking for MultiQC report in %s" % multiqc_dir)
        multiqc_report = os.path.join(multiqc_dir,
                                      "multi%s_report.html"
                                      % os.path.basename(qc_dir.path))
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
            name=self.name,
            software=software,
            fastqs=[],
            output_files=output_files,
            tags=sorted(list(tags))
        )

    @classmethod
    def verify(self,params,qc_outputs):
        """
        Verify 'multiqc' QC module against outputs

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
        return None
