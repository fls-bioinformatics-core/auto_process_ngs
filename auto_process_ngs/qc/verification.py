#!/usr/bin/env python
#
#     verification: utilities for verification of QC outputs
#     Copyright (C) University of Manchester 2022-2024 Peter Briggs
#

"""
Utilities for verifying QC pipeline outputs.

Provides the following classes:

- QCVerifier: enables verification of QC outputs against protocols

Provides the following functions:

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
from .fastq_screen import LEGACY_SCREENS
from .modules.cellranger_arc_count import CellrangerArcCount
from .modules.cellranger_atac_count import CellrangerAtacCount
from .modules.cellranger_count import CellrangerCount
from .modules.cellranger_multi import CellrangerMulti
from .modules.fastqc import Fastqc
from .modules.fastq_screen import FastqScreen
from .modules.multiqc import Multiqc
from .modules.picard_insert_size_metrics import PicardInsertSizeMetrics
from .modules.qualimap_rnaseq import QualimapRnaseq
from .modules.rseqc_genebody_coverage import RseqcGenebodyCoverage
from .modules.rseqc_infer_experiment import RseqcInferExperiment
from .modules.sequence_lengths import SequenceLengths
from .modules.strandedness import Strandedness
from .protocols import fetch_protocol_definition
from .protocols import parse_qc_module_spec
from .outputs import QCOutputs
from .utils import filter_fastqs
from .utils import get_bam_basename
from ..tenx.cellplex import CellrangerMultiConfigCsv
from ..utils import normalise_organism_name

# Module specific logger
logger = logging.getLogger(__name__)

######################################################################
# Supported QC module classes
######################################################################

QC_MODULES = (CellrangerCount,
              CellrangerAtacCount,
              CellrangerArcCount,
              CellrangerMulti,
              Fastqc,
              FastqScreen,
              Multiqc,
              PicardInsertSizeMetrics,
              QualimapRnaseq,
              RseqcGenebodyCoverage,
              RseqcInferExperiment,
              SequenceLengths,
              Strandedness)

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

    def verify(self,protocol,fastqs,organism=None,
               fastq_screens=None,star_index=None,
               annotation_bed=None,annotation_gtf=None,
               cellranger_version=None,cellranger_refdata=None,
               cellranger_use_multi_config=None,
               seq_data_samples=None):
        """
        Verify QC outputs for Fastqs against specified protocol

        Arguments:
          protocol (QCProtocol): QC protocol to verify against
          fastqs (list): list of Fastqs to verify outputs for
          organism (str): organism associated with outputs
          fastq_screens (list): list of panel names to verify
            FastqScreen outputs against
          star_index (str): path to STAR index
          annotation_bed (str): path to BED annotation file
          annotation_gtf (str): path to GTF annotation file
          cellranger_version (str): specific version of 10x
            package to check for
          cellranger_refdata (str): specific 10x reference
            dataset to check for
          cellranger_use_multi_config (bool): if True then
            cellranger count verification will attempt to
            use data (GEX samples and reference dataset) from
            the '10x_multi_config.csv' file
          seq_data_samples (list): list of sample names with
            sequence (i.e. biological) data

        Returns:
          Boolean: True if all expected outputs are present,
            False otherwise.
        """
        # Sample names
        samples = set()
        for fq in fastqs:
            samples.add(self.fastq_attrs(fq).sample_name)
        samples = sorted(list(samples))

        # Subsets of samples and Fastqs with sequence data
        # (i.e. biological data rather than feature barcodes etc)
        if not seq_data_samples:
            # No explicit list supplied so try to guess
            seq_data_samples = self.identify_seq_data(samples)
        seq_data_fastqs = [fq for fq in fastqs
                           if self.fastq_attrs(fq).sample_name
                           in seq_data_samples]

        # Seq data Fastqs defaults to all Fastqs
        if seq_data_fastqs is None:
            seq_data_fastqs = fastqs

        # Default parameters for verification
        default_params = dict(
            qc_dir=self.qc_dir,
            fastqs=fastqs,
            samples=samples,
            seq_data_fastqs=seq_data_fastqs,
            seq_data_samples=seq_data_samples,
            seq_data_reads=protocol.reads.seq_data,
            qc_reads=protocol.reads.qc,
            organism=organism,
            fastq_screens=fastq_screens,
            star_index=star_index,
            annotation_bed=annotation_bed,
            annotation_gtf=annotation_gtf,
            cellranger_version=cellranger_version,
            cellranger_refdata=cellranger_refdata,
            cellranger_use_multi_config=cellranger_use_multi_config
        )

        # Perform verification
        verified = dict()
        params_for_module = dict()

        for qc_module in protocol.qc_modules:

            # Handle QC module specification
            qc_module,module_params = parse_qc_module_spec(qc_module)

            # Store parameters for reporting
            params_for_module[qc_module] = dict(**module_params)

            # Initialise up parameters for this module
            params = AttributeDictionary(**default_params)

            # Override parameters from module definition
            # parameter list
            for p in module_params:
                params[p] = module_params[p]

            # Verify outputs for this QC module
            verified[qc_module] = self.verify_qc_module(qc_module,
                                                        params)
        # Make templates for parameter and status
        field_width = 21
        for qc_module in protocol.qc_modules:
            field_width = max(len(qc_module),field_width)
        parameter_template_str = "{parameter:%ss}: {value}" % field_width
        qc_module_template_str = "{name:%ss}: {status:4s}{params}" % field_width
        # Report parameters and status of checks
        print("-"*(10+len(self.qc_dir)))
        print("QC dir  : %s" % self.qc_dir)
        print("Protocol: %s" % protocol.name)
        print("-"*(10+len(self.qc_dir)))
        # Report parameters
        print("Parameters:")
        for p in default_params:
            if p == 'fastqs':
                fqs = ['.../%s' % os.path.basename(fq)
                       for fq in default_params[p]]
                if not fqs:
                    print(parameter_template_str.format(parameter=p,
                                                        value=''))
                else:
                    print(parameter_template_str.format(parameter=p,
                                                        value=fqs[0]))
                    for fq in fqs[1:]:
                        print(parameter_template_str.format(parameter='',
                                                            value=fq))
            elif p == 'samples':
                smpls = default_params[p]
                if not smpls:
                    print(parameter_template_str.format(parameter=p,
                                                        value=''))
                else:
                    print(parameter_template_str.format(parameter=p,
                                                        value=smpls[0]))
                    for smpl in smpls[1:]:
                        print(parameter_template_str.format(parameter='',
                                                            value=smpl))
            elif p == 'cellranger_refdata':
                refdata = default_params[p]
                print(parameter_template_str.format(
                    parameter=p,
                    value=('.../%s' % os.path.basename(refdata)
                           if refdata else refdata)))
            else:
                print(parameter_template_str.format(parameter=p,
                                                    value=default_params[p]))
        print("-"*(field_width+6))
        # Report the status for each module
        for name in verified:
            # Get string version of status
            if verified[name] is None:
                qc_module_status = '****'
            elif verified[name]:
                qc_module_status = 'PASS'
            else:
                qc_module_status = 'FAIL'
            # Report status of module
            print(qc_module_template_str.format(
                name=name,
                status=qc_module_status,
                params=(" %s" % params_for_module[name]
                        if params_for_module[name] else '')))
        # Status for QC as a whole
        status = all([verified[m] for m in verified
                      if verified[m] is not None])
        print("-"*(field_width+6))
        print(qc_module_template_str.format(
            name="QC STATUS",
            status=('PASS' if status else 'FAIL'),
            params=''))
        print("-"*(field_width+6))

        # Return verification status
        return status

    def verify_qc_module(self,name,params):
        """
        Verify QC outputs for specific QC module

        Arguments:
          name (str): QC module name
          params (AttributeDictionary): parameters to
            verify QC module using

        Returns:
          Boolean: True if all outputs are present, False
            if one or more are missing.

        Raises:
          Exception: if the specified QC module name is not
            recognised.
        """
        for m in QC_MODULES:
            if m.name == name:
                return m.verify(params,self.data(name))
        # No match
        raise Exception("unknown QC module: '%s'" % name)

    def identify_seq_data(self,samples):
        """
        Identify samples with sequence (biological) data

        Arguments:
          samples (list): list of all sample names

        Returns:
          List: subset of sample names with sequence data.
        """
        # Check for 10x_multi_config.csv
        if "10x_multi_config.csv" in self.config_files:
            # Get GEX sample names from multi config file
            cf = CellrangerMultiConfigCsv(
                os.path.join(self.qc_dir,"10x_multi_config.csv"))
            seq_data = [s for s in cf.gex_libraries if s in samples]
        else:
            seq_data = [s for s in samples]
        return seq_data

#######################################################################
# Functions
#######################################################################

def verify_project(project,qc_dir=None,qc_protocol=None,
                   fastqs=None):
    """
    Check the QC outputs are correct for a project

    Arguments:
      project (AnalysisProject): project to verify QC for
      qc_dir (str): path to the QC output dir; relative
        path will be treated as a subdirectory of the
        project being checked.
      qc_protocol (str): QC protocol name or specification
        to verify against (optional)
      fastqs (list): list of Fastqs to include (optional,
        defaults to Fastqs in the project)

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
    star_index = None
    annotation_bed=None
    annotation_gtf=None
    organism = None
    seq_data_samples = None
    fastq_screens = None
    qc_info_file = os.path.join(qc_dir,"qc.info")
    if os.path.exists(qc_info_file):
        # Get QC parameters from metadata file
        qc_info = AnalysisProjectQCDirInfo(filen=qc_info_file)
        # QC protocol
        if not qc_protocol:
            if qc_info['protocol_specification']:
                qc_protocol = qc_info['protocol_specification']
            else:
                qc_protocol = qc_info['protocol']
        # Organism
        try:
            organism = qc_info['organism']
        except KeyError:
            pass
        # Fastqs
        try:
            if not fastqs:
                if qc_info['fastqs']:
                    fastqs = qc_info['fastqs'].split(',')
        except KeyError:
            pass
        # Samples with biological data
        try:
            seq_data_samples = qc_info['seq_data_samples']
            if seq_data_samples:
                seq_data_samples = seq_data_samples.split(',')
        except KeyError:
            pass
        # Cellranger reference data and version
        try:
            cellranger_refdata = qc_info['cellranger_refdata']
        except KeyError:
            pass
        try:
            cellranger_version = qc_info['cellranger_version']
        except KeyError:
            pass
        # STAR index
        try:
            star_index = qc_info['star_index']
        except KeyError:
            pass
        # Reference annotation
        try:
            annotation_bed = qc_info['annotation_bed']
        except KeyError:
            pass
        try:
            annotation_gtf = qc_info['annotation_gtf']
        except KeyError:
            pass
        # Fastq screens
        try:
            fastq_screens = qc_info['fastq_screens']
            if fastq_screens:
                fastq_screens = fastq_screens.split(',')
            elif 'fastq_screens' not in qc_info.keys_in_file():
                fastq_screens = LEGACY_SCREENS
        except KeyError:
            pass
    logger.debug("verify: cellranger reference data : %s" %
                 cellranger_refdata)
    logger.debug("verify: fastq screens : %s" % (fastq_screens,))
    protocol = fetch_protocol_definition(qc_protocol)
    if fastqs:
        fastqs_in = fastqs
    else:
        fastqs_in = project.fastqs
    verifier = QCVerifier(qc_dir,
                          fastq_attrs=project.fastq_attrs)
    return verifier.verify(protocol,
                           fastqs_in,
                           organism=organism,
                           seq_data_samples=seq_data_samples,
                           fastq_screens=fastq_screens,
                           star_index=star_index,
                           annotation_bed=annotation_bed,
                           annotation_gtf=annotation_gtf,
                           cellranger_version=cellranger_version,
                           cellranger_refdata=cellranger_refdata)
