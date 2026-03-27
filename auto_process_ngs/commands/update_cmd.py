#!/usr/bin/env python
#
#     update_cmd.py: implement 'update' command
#     Copyright (C) University of Manchester 2023-2026 Peter Briggs
#
#########################################################################

#######################################################################
# Imports
#######################################################################

import os
import logging
from ..qc.utils import report_qc

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Command functions
#######################################################################

def update(ap, update_paths=True, update_project_metadata=True,
           update_sync_projects=True, update_qc_reports=True):
    """
    Update metadata and artefacts in analysis directory

    Arguments:
      ap (AutoProcessor): autoprocessor pointing to the
        analysis directory to publish QC for
      update_paths (bool): whether to update analysis
        directory paths in metadata and parameter files
        (default: True)
      update_project_metadata (bool): whether to update
        metadata stored in 'projects.info' and in the
        project directories (default: True)
      update_sync_projects (bool): whether to update
        projects listed in 'projects.info' against
        project directories on the filesystem (default:
        True)
      update_qc_reports (bool): whether to update QC
        reports in projects where existing report is
        older than the project metadata file (default:
        True)
    """
    if not (update_paths or update_project_metadata or
            update_sync_projects or update_qc_reports):
        logger.warning("No updates requested")

    if update_paths:
        # Update paths if analysis dir has been moved or copied
        ap.update_paths()

    if update_sync_projects:
        # Synchronise projects listed in 'projects.info'
        # with contents of analysis directory
        ap.sync_project_metadata_file()

    if update_project_metadata:
        # Synchronise metadata in each project with the
        # contents of 'projects.info'
        ap.sync_project_metadata()

    if update_qc_reports:
        # Update QC reports that are older than project metadata
        for project in ap.get_analysis_projects():
            metadata_mtime = os.path.getmtime(project.info_file)
            for qc_dir in project.qc_dirs:
                qc_report = os.path.join(project.dirn,
                                         "%s_report.html" % qc_dir)
                if not os.path.exists(qc_report):
                    continue
                qc_report_mtime = os.path.getmtime(qc_report)
                if metadata_mtime > qc_report_mtime:
                    print("...regenerating QC report for '%s/%s'" %
                          (project.name,qc_dir))
                    report_status = report_qc(
                        project,
                        qc_dir=qc_dir,
                        multiqc=True,
                        force=True,
                        runner=ap.settings.runners.publish_qc,
                        log_dir=ap.tmp_dir)
                    if report_status == 0:
                        print("...ok")
                    else:
                        print("...failed")
