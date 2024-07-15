#!/usr/bin/env python3
#
#     fastqc: implements 'cellranger_count' QC module
#     Copyright (C) University of Manchester 2024 Peter Briggs

"""
Implements the 'cellranger-atac_count' QC module:

* CellrangerAtacCount: core QCModule class

Currently only the 'verify' method of the 'CellrangerAtacCount' class
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

class CellrangerAtacCount(CellrangerCount):
    """
    Class for handling the 'cellranger-atac_count' QC module
    """
    name = "cellranger-atac_count"
    mapped_metrics = False
    runners = ("verify_runner",
               "cellranger_count_runner")
    envmodules = ("cellranger",)
    
    def __init__(self):
        CellrangerCount.__init__(self)

    @classmethod
    def verify(self,params,qc_outputs):
        """
        Verify 'cellranger_count' QC module against outputs

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
        if not params.samples:
            # No samples so cellranger-atac outputs not
            # expected
            return True
        if params.cellranger_refdata is None:
            # No reference data so cellranger-atac outputs not
            # expected
            return True
        # Check expected samples against actual samples
        # associated with specified version and dataset
        return verify_10x_pipeline(('cellranger-atac',
                                    params.cellranger_version,
                                    params.cellranger_refdata),
                                   params.samples,
                                   qc_outputs)
