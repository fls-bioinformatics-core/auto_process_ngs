#!/usr/bin/env python
#
#     commands/__init__.py: core functions for autoprocess commands
#     Copyright (C) University of Manchester 2017-2019 Peter Briggs
#
#######################################################################
# Imports
#######################################################################

from .setup_cmd import setup
from .setup_analysis_dirs_cmd import setup_analysis_dirs
from .make_fastqs_cmd import make_fastqs
from .analyse_barcodes_cmd import analyse_barcodes
from .run_qc_cmd import run_qc
from .publish_qc_cmd import publish_qc
from .archive_cmd import archive
from .report_cmd import report
from .merge_fastq_dirs_cmd import merge_fastq_dirs
from .update_fastq_stats_cmd import update_fastq_stats
from .import_project_cmd import import_project
from .clone_cmd import clone
