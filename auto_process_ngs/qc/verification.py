#!/usr/bin/env python
#
#     verification: utilities for verification of QC outputs
#     Copyright (C) University of Manchester 2022 Peter Briggs
#

"""
Utilities for verifying QC pipeline outputs.

Provides the following classes:

- QCVerifier: enables verification of QC outputs against protocols

Provides the following functions:

- parse_qc_module_spec: process QC module specification string
- filter_fastqs: filter list of Fastqs based on read IDs
- filter_10x_pipelines: filter list of 10xGenomics pipeline tuples
- verify_project: check the QC outputs for a project
"""

#######################################################################
# Imports
#######################################################################

import os
import logging
from bcftbx.utils import AttributeDictionary
from ..analysis import AnalysisFastq
from ..metadata import AnalysisProjectQCDirInfo
from .protocols import fetch_protocol_definition
from .outputs import QCOutputs
from ..tenx_genomics_utils import CellrangerMultiConfigCsv

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Classes
#######################################################################

class QCVerifier(QCOutputs):
    """
    Class to perform verification of QC outputs

    The QCVerifier enables the QC outputs from a directory
    to be checked against arbitrary QC protocols via its
    ``verify`` method.

    For example:

    >>> project = AnalysisProject("/data/projects/PJB")
    >>> verifier = QCVerifier(project.qc_dir)
    >>> verifier.verify(project.fastqs,"standardPE")
    True

    Arguments:
      qc_dir (str): path to directory to examine
      fastq_attrs (BaseFastqAttrs): (optional) class for
        extracting data from Fastq names
    """
    def __init__(self,qc_dir,fastq_attrs=None):
        QCOutputs.__init__(self,qc_dir,fastq_attrs=fastq_attrs)

    def verify(self,fastqs,qc_protocol,fastq_screens=None,
               cellranger_version=None,cellranger_refdata=None):
        """
        Verify QC outputs for Fastqs against specified protocol

        Arguments:
          fastqs (list): list of Fastqs to verify outputs for
          qc_protocol (str): QC protocol to verify against
          fastq_screens (list): list of panel names to verify
            FastqScreen outputs against
          cellranger_version (str): specific version of 10x
            package to check for
          cellranger_refdata (str): specific 10x reference
            dataset to check for

        Returns:
          Boolean: True if all expected outputs are present,
            False otherwise.
        """
        # Look up protocol definition
        print("QC protocol: %s" % qc_protocol)
        reads,qc_modules = fetch_protocol_definition(qc_protocol)

        # Sample names
        samples = set()
        for fq in fastqs:
            samples.add(self.fastq_attrs(fq).sample_name)
        samples = sorted(list(samples))
        print("Samples: %s" % samples)

        # Default parameters for verification
        default_params = dict(
            fastqs=fastqs,
            samples=samples,
            data_reads=reads.data,
            qc_reads=reads.qc,
            fastq_screens=fastq_screens,
            cellranger_version=cellranger_version,
            cellranger_refdata=cellranger_refdata
        )

        # Perform verification
        verified = dict()

        for qc_module in qc_modules:

            # Handle QC module specification
            qc_module,module_params = parse_qc_module_spec(qc_module)

            # Initialise up parameters for this module
            params = AttributeDictionary(**default_params)

            # Override parameters from module definition
            # parameter list
            for p in module_params:
                params[p] = module_params[p]

            # Verify outputs for this QC module
            verified[qc_module] = self.verify_qc_module(qc_module,
                                                        **params)

        # Report status of checks
        print("-"*27)
        for name in verified:
            print("%21s: %s" % (name,('PASS' if verified[name]
                                      else 'FAIL')))
        status = all([verified[m] for m in verified])
        print("-"*27)
        print("%21s: %s" % ("QC STATUS",('PASS' if status else 'FAIL')))
        print("-"*27)

        # Return verification status
        return status

    def verify_qc_module(self,name,fastqs=None,samples=None,
                         data_reads=None,qc_reads=None,
                         fastq_screens=None,
                         cellranger_version=None,
                         cellranger_refdata=None):
        """
        Verify QC outputs for specific QC module

        Arguments:
          name (str): QC module name
          fastqs (list): list of Fastqs
          samples (list): list of sample names
          data_reads (list): list of reads containing
            sequence data
          qc_reads (list): list of reads to perform general
            QC on
          fastq_screens (list): list of panel names to verify
            FastqScreen outputs against
          cellranger_version (str): specific version of 10x
            package to check for
          cellranger_refdata (str): specific 10x reference
            dataset to check for

        Returns:
          Boolean: True if all outputs are present, False
            if one or more are missing.

        Raises:
          Exception: if the specified QC module name is not
            recognised.
        """
        print("%s: checking QC outputs" % name)

        # Perform checks based on QC module
        if name == "fastqc":
            if not fastqs:
                # Nothing to check
                return True
            try:
                # Filter Fastq names
                fastqs = self.filter_fastqs(qc_reads,fastqs)
                # Check that outputs exist for every Fastq
                for fq in fastqs:
                    if fq not in self.data('fastqc').fastqs:
                        return False
                return True
            except KeyError:
                # No Fastqc outputs present
                return False

        elif name == "fastq_screen":
            if not fastqs or not fastq_screens:
                # Nothing to check
                return True
            try:
                # Filter Fastq names
                fastqs = self.filter_fastqs(data_reads,fastqs)
                # Check outputs exist for each screen
                for screen in fastq_screens:
                    if screen not in self.data('fastq_screen').\
                       screen_names:
                        # No outputs associated with screen
                        return False
                    # Check outputs exist for each Fastq
                    for fq in fastqs:
                        if fq not in self.data('fastq_screen').\
                           fastqs_for_screen[screen]:
                            return False
                return True
            except KeyError as ex:
                # No FastqScreen outputs present
                return False

        elif name == "sequence_lengths":
            if not fastqs:
                # Nothing to check
                return True
            try:
                # Filter Fastq names
                fastqs = self.filter_fastqs(qc_reads,fastqs)
                # Check that outputs exist for every Fastq
                for fq in fastqs:
                    if fq not in self.data('sequence_lengths').fastqs:
                        return False
                return True
            except KeyError:
                # No sequence length outputs present
                return False

        elif name == "strandedness":
            if not fastqs or \
               "fastq_strand.conf" not in self.config_files:
                # No Fastqs or no conf file so strandedness
                # outputs not expected
                return True
            if "strandedness" not in self.outputs:
                # No strandedness outputs present
                return False
            # Filter Fastq names
            fastqs = self.filter_fastqs(data_reads[:1],fastqs)
            # Check that outputs exist for every Fastq
            for fq in fastqs:
                if fq not in self.data('fastq_strand').fastqs:
                    return False
            return True

        elif name == "multiqc":
            return ("multiqc" in self.outputs)

        elif name == "cellranger_count":
            print(self.data("cellranger_count"))
            if not samples:
                # No sample so cellranger outputs not expected
                return True
            if "cellranger_count" not in self.outputs:
                # No cellranger outputs present
                return False
            # Check expected samples against actual samples
            # associated with specified version and dataset
            pipelines = filter_10x_pipelines(
                ('cellranger',cellranger_version,cellranger_refdata),
                self.data('cellranger_count').pipelines)
            print("Matching cellranger pipelines: %s" % pipelines)
            for pipeline in pipelines:
                verified_pipeline = True
                for sample in samples:
                    if sample not in self.data("cellranger_count").\
                       samples_by_pipeline[pipeline]:
                        # At least one sample missing outputs from
                        # this pipeline, so move on to the next
                        verified_pipeline = False
                        break
                if verified_pipeline:
                    return True
            return False

        elif name == "cellranger-atac_count":
            if not samples:
                # No sample so cellranger-atac outputs not
                # expected
                return True
            if "cellranger-atac_count" not in self.outputs:
                # No cellranger-atac outputs present
                return False
            # Check expected samples against actual samples
            # associated with specified version and dataset
            pipelines = filter_10x_pipelines(
                ('cellranger-atac',cellranger_version,cellranger_refdata),
                self.data('cellranger_count').pipelines)
            for pipeline in pipelines:
                verified_pipeline = True
                for sample in samples:
                    if sample not in self.data("cellranger_count").\
                       samples_by_pipeline[pipeline]:
                        # At least one sample missing outputs from
                        # this pipeline, so move on to the next
                        verified_pipeline = False
                        break
                if verified_pipeline:
                    return True
            return False

        elif name == "cellranger-arc_count":
            # Look for cellranger-arc config files and
            # make a list of expected samples
            expected_samples = []
            for cf in self.config_files:
                if cf.startswith("libraries.") and \
                   cf.endswith(".csv"):
                    sample = '.'.join(cf.split('.')[1:-1])
                    expected_samples.append(sample)
            if not expected_samples:
                # No libraries to check
                return True
            if "cellranger-arc_count" not in self.outputs:
                # No cellranger-arc outputs present
                return False
            # Check expected samples against actual samples
            # associated with specified version and dataset
            pipelines = filter_10x_pipelines(
                ('cellranger-arc',cellranger_version,cellranger_refdata),
                self.data('cellranger_count').pipelines)
            for pipeline in pipelines:
                verified_pipeline = True
                for sample in samples:
                    if sample not in self.data("cellranger_count").\
                       samples_by_pipeline[pipeline]:
                        # At least one sample missing outputs from
                        # this pipeline, so move on to the next
                        verified_pipeline = False
                        break
                if verified_pipeline:
                    return True
            return False

        elif name == "cellranger_multi":
            if "10x_multi_config.csv" not in self.config_files:
                # No multi config file so no outputs expected
                return True
            if "cellranger_multi" not in self.outputs:
                return False
            # Get expected multiplexed sample names
            # from config file
            cf = os.path.join(self.qc_dir,"10x_multi_config.csv")
            expected_samples = CellrangerMultiConfigCsv(cf).\
                               sample_names
            if not expected_samples:
                # No samples to check outputs for
                return True
            # Check against actual samples
            # associated with specified version and dataset
            pipelines = filter_10x_pipelines(
                ('cellranger',cellranger_version,cellranger_refdata),
                self.data('cellranger_multi').pipelines)
            for pipeline in pipelines:
                verified_pipeline = True
                for sample in expected_samples:
                    if sample not in self.data("cellranger_multi").\
                       samples_by_pipeline[pipeline]:
                        # At least one sample missing outputs from
                        # this pipeline, so move on to the next
                        verified_pipeline = False
                        break
                if verified_pipeline:
                    return True
            return False

        else:
            raise Exception("unknown QC module: '%s'" % name)

    def filter_fastqs(self,reads,fastqs):
        """
        Filter list of Fastqs and return names matching reads

        Wrap external 'filter_fastqs' function
        """
        return filter_fastqs(reads,
                             fastqs,
                             fastq_attrs=self.fastq_attrs)

#######################################################################
# Functions
#######################################################################

def parse_qc_module_spec(module_spec):
    """
    Parse QC module spec into name and parameters

    Parse a QC module specification of the form
    ``NAME`` or ``NAME(KEY=VALUE;...)`` and return
    the module name and any additional parameters
    in the form of a dictionary.

    For example:

    >>> parse_qc_module_spec('NAME')
    ('NAME', {})
    >>> parse_qc_module_spec('NAME(K1=V1;K2=V2)')
    ('NAME', { 'K1':'V1', 'K2':'V2' })

    Arguments:
      module_spec (str): QC module specification

    Returns:
      Tuple: tuple of the form (name,params) where
        'name' is the QC module name and 'params'
        is a dictionary with the extracted key-value
        pairs.
    """
    # Handle module specification string of the form
    # 'NAME[(KEY=VALUE;...)]'
    items = module_spec.split('(')
    # Extract the module name and associated parameter list
    name = items[0]
    params = {}
    try:
        for item in items[1].rstrip(')').split(';'):
            key,value = item.split('=')
            params[key] = value
    except IndexError:
        pass
    return (name,params)

def filter_fastqs(reads,fastqs,fastq_attrs=AnalysisFastq):
    """
    Filter list of Fastqs and return names matching reads

    Arguments:
      reads (list): list of reads to filter ('r1',
        'i2' etc: '*' matches all reads, 'r*' matches
        all data reads, 'i*' matches all index reads)
      fastqs (list): list of Fastq files or names
        to filter
      fastq_attrs (BaseFastqAttrs): class for extracting
        attribute data from Fastq names

    Returns:
      List: matching Fastq names (i.e. no leading
        path or trailing extensions)
    """
    fqs = set()
    for read in reads:
        index_read = (read.startswith('i') or read == '*')
        if read == '*':
            # All reads
            for fastq in fastqs:
                fqs.add(fastq_attrs(fastq).basename)
            continue
        if read[1:] == '*':
            # All read numbers
            read_number = None
        else:
            # Specific read
            read_number = int(read[1:])
        for fastq in fastqs:
            fq = fastq_attrs(fastq)
            if (not index_read and fq.is_index_read) or \
               (index_read and not fq.is_index_read):
                # Skip index reads
                continue
            if fq.read_number == read_number or \
               not read_number:
                fqs.add(fq.basename)
    return sorted(list(fqs))

def filter_10x_pipelines(p,pipelines):
    """
    Filter list of 10x pipelines

    Pipelines are described using tuples of the form:

    (NAME,VERSION,REFERENCE)

    for example:

    ('cellranger','6.1.2','refdata-gex-2020')

    Only pipelines matching the specified name, version
    and reference data will be included in the returned
    list.

    Where the supplied version or reference dataset name
    are either None or '*', these will match any version
    and/or reference dataset.

    Arguments:
      p (tuple): tuple specifying pipeline(s) to match
        against
      pipelines (list): list of pipeline tuples to filter

    Returns:
      List: list of matching 10x pipeline tuples.
    """
    # Extract elements from pipeline pattern
    name = p[0]
    version = p[1]
    refdata = p[2]
    # Check for wildcard versions and reference data
    if version == "*":
        version = None
    if refdata == "*":
        refdata = None
    # Normalise reference dataset name
    refdata = (os.path.basename(refdata) if refdata else None)
    # Find all matching pipelines
    matching_pipelines = list()
    for pipeline in pipelines:
        if pipeline[0] != name:
            # Wrong 10x package name
            continue
        if version and pipeline[1] != version:
            # Wrong version
            continue
        if refdata and pipeline[2] != refdata:
            # Wrong reference dataset
            continue
        # Passed all filters
        matching_pipelines.append(pipeline)
    return matching_pipelines

def verify_project(project,qc_dir=None,qc_protocol=None):
    """
    Check the QC outputs are correct for a project

    Arguments:
      project (AnalysisProject): project to verify QC for
      qc_dir (str): path to the QC output dir; relative
        path will be treated as a subdirectory of the
        project being checked.
      qc_protocol (str): QC protocol to verify against
        (optional)

     Returns:
       Boolean: Returns True if all expected QC products
         are present, False if not.
    """
    logger.debug("verify: qc_dir (initial): %s" % qc_dir)
    if qc_dir is None:
        qc_dir = project.qc_dir
    else:
        if not os.path.isabs(qc_dir):
            qc_dir = os.path.join(project.dirn,
                                  qc_dir)
    logger.debug("verify: qc_dir (final)  : %s" % qc_dir)
    cellranger_version = None
    cellranger_refdata = None
    fastq_screens = None
    qc_info_file = os.path.join(qc_dir,"qc.info")
    if os.path.exists(qc_info_file):
        qc_info = AnalysisProjectQCDirInfo(filen=qc_info_file)
        if not qc_protocol:
            qc_protocol = qc_info['protocol']
        try:
            cellranger_refdata = qc_info['cellranger_refdata']
        except KeyError:
            pass
        try:
            cellranger_version = qc_info['cellranger_version']
        except KeyError:
            pass
        try:
            fastq_screens = qc_info['fastq_screens']
            if fastq_screens:
                fastq_screens = fastq_screens.split(',')
            elif 'fastq_screens' not in qc_info.keys_in_file():
                fastq_screens = FASTQ_SCREENS
        except KeyError:
            pass
    logger.debug("verify: cellranger reference data : %s" %
                 cellranger_refdata)
    logger.debug("verify: fastq screens : %s" % fastq_screens)
    verifier = QCVerifier(qc_dir,
                          fastq_attrs=project.fastq_attrs)
    return verifier.verify(project.fastqs,
                           qc_protocol,
                           fastq_screens=fastq_screens,
                           cellranger_version=cellranger_version,
                           cellranger_refdata=cellranger_refdata)
