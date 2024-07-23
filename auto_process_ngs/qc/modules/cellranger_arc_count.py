#!/usr/bin/env python3
#
#     cellranger_arc_count: implements 'cellranger-arc_count' QC module
#     Copyright (C) University of Manchester 2024 Peter Briggs

"""
Implements the 'cellranger-arc_count' QC module:

* CellrangerArcCount: core QCModule class

Currently only the 'verify' method of the 'CellrangerArcCount' class
is implemented here; the remainder of the methods are implicitly
delegated to the base 'CellrangerCount' class.
"""

#######################################################################
# Imports
#######################################################################

import logging
from .cellranger_count import CellrangerCount
from .cellranger_count import verify_10x_pipeline

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Core class
#######################################################################

class CellrangerArcCount(CellrangerCount):
    """
    Class for handling the 'cellranger-arc_count' QC module
    """
    name = "cellranger-arc_count"
    mapped_metrics = False
    runners = ("verify_runner",
               "cellranger_count_runner")
    envmodules = ("cellranger",)
    
    def __init__(self):
        CellrangerCount.__init__(self)

    @classmethod
    def verify(self,params,qc_outputs):
        """
        Verify 'cellranger-arc_count' QC module against outputs

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
        # Look for cellranger-arc config files and
        # make a list of expected samples
        expected_samples = []
        for cf in qc_outputs.config_files:
            if cf.startswith("libraries.") and \
               cf.endswith(".csv"):
                sample = '.'.join(cf.split('.')[1:-1])
                expected_samples.append(sample)
        if not expected_samples:
            # No libraries to check
            return True
        ##if "cellranger-arc_count" not in self.outputs:
        ##    # No cellranger-arc outputs present
        ##    return False
        # Check expected samples against actual samples
        # associated with specified version and dataset
        return verify_10x_pipeline(('cellranger-arc',
                                    params.cellranger_version,
                                    params.cellranger_refdata),
                                   params.samples,
                                   qc_outputs)
