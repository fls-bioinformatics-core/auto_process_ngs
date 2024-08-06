#!/usr/bin/env python
#
#     qc_modules.py: list modules for QC pipeline
#     Copyright (C) University of Manchester 2024 Peter Briggs
#

"""
Provides the following constants for the QC pipeline:

- QC_MODULES: list of QC module classes
- QC_MODULE_NAMES: list of QC module names
"""

#######################################################################
# Imports
#######################################################################

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

#######################################################################
# Data
#######################################################################

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

QC_MODULE_NAMES = sorted([m.name for m in QC_MODULES])
