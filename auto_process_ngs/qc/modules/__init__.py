#!/usr/bin/env python3
#
#     modules: base and utility classes for QC modules
#     Copyright (C) University of Manchester 2024 Peter Briggs


#######################################################################
# Imports
#######################################################################

import os
import logging
from ...metadata import AnalysisProjectQCDirInfo

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Classes
#######################################################################

class QCModule:
    """
    Base class for defining and handling QC modules
    """
    name = "base"
    mapped_metrics = False
    require_bam_files = False
    
    def __init__(self):
        pass

    @classmethod
    def collect_qc_outputs(self,qc_dir):
        """
        Collect information on QC module outputs

        Should be implemented by the subclass to
        handle the specific outputs from the QC
        module, and return an AttributeDictionary
        with the collected data.

        Arguments:
          qc_dir (QCDir): QC directory to examine
        """
        raise NotImplementedError("Subclass should implement "
                                  "'collect_qc_outputs'")

    @classmethod
    def verify(self,params,qc_outputs):
        """
        Verify the outputs of the QC module

        Should be implemented by the subclass to
        perform the verification based on the
        supplied parameters and the collected
        outputs.

        Should return one of 3 values:

        - True: outputs verified ok
        - False: outputs failed to verify
        - None: verification not possible

        Arguments:
          params (AttributeDictionary): values
            of parameters used as inputs to the
            QC module
          qc_outputs (AttributeDictionary):
            QC outputs returned from the
            'collect_qc_outputs' method
        """
        raise NotImplementedError("Subclass should implement "
                                  "'collect_qc_outputs'")


    def add_to_pipeline(self,p,*args,**kws):
        """
        Add the QC module tasks to a pipeline

        Should be implemented by the subclass
        to add the tasks required to execute the
        QC module to an existing pipeline,
        and return the final task instance.

        Arguments:
          p (Pipeline): pipeline instance
          args (list): positional arguments required
            by the QC module
          kws (mapping): additional optional arguments
            for configuring the module 
        """
        raise NotImplementedError("Subclass should implement "
                                  "'add_to_pipeline'")
        
class QCDir:
    """
    Helper class for interacting with a QC directory

    Provides the following attributes/properties:

    - path: full path to the QC directory
    - fastq_attrs: class for extracting information
      from Fastq filenames
    - info: AnalysisProjectQCDirInfo instance loaded
      with metadata from 'qc.info' (or None if there
      is no info file)
    - file_list: list of paths for files and
      directories under the the QC directory

    Arguments:
      qc_dir (str): path to the QC dir
      fastq_attrs (object): optional class for
        extracting info from Fastq filenames
        (default: AnalysisFastq)
    """
    def __init__(self,qc_dir,fastq_attrs=None):
        # QC directory
        self.path = os.path.abspath(qc_dir)
        # How to handle Fastq names
        if fastq_attrs is not None:
            self.fastq_attrs = fastq_attrs
        else:
            self.fastq_attrs = AnalysisFastq
        # QC metadata
        self.info = None
        info_file = os.path.join(self.path,"qc.info")
        if os.path.exists(info_file):
            self.info = AnalysisProjectQCDirInfo(info_file)
        # File list
        self._file_list = None

    @property
    def file_list(self):
        """
        Returns list of file and directory paths under QC dir

        The list contains the full path for each file
        and directory.
        """
        if self._file_list is None:
            self._file_list = [os.path.join(self.path,f)
                               for f in os.listdir(self.path)]
        return self._file_list
