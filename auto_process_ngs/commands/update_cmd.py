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
from ..analysis import AnalysisProject
from ..metadata import AnalysisProjectInfo
from ..metadata import AnalysisProjectQCDirInfo
from ..qc.utils import report_qc
from ..utils import sort_sample_names

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
        # Update paths in the top-level parameter file
        # (if analysis dir has been moved or copied)
        if ap.params.analysis_dir != ap.analysis_dir:
            print("Updating analysis directory paths in parameter file")
            old_dir = ap.params.analysis_dir
            for p in ('analysis_dir',
                      'primary_data_dir',
                      'sample_sheet'):
                if not ap.params[p]:
                    continue
                print("...updating '%s'" % p)
                ap.params[p] = os.path.normpath(
                    os.path.join(ap.analysis_dir,
                                 os.path.relpath(ap.params[p],old_dir)))
            # Update paths in QC metadata in projects
            for project in ap.get_analysis_projects_from_dirs():
                # Iterate through all project directories
                for qc_dir in project.qc_dirs:
                    qc_info = AnalysisProjectQCDirInfo(
                        os.path.join(project.dirn,qc_dir,"qc.info"))
                    print("...updating QC info for %s/%s" % (project.name,
                                                             qc_dir))
                    qc_info['fastq_dir'] = os.path.normpath(
                        os.path.join(ap.analysis_dir,
                                     os.path.relpath(qc_info.fastq_dir,
                                                     old_dir)))
                    qc_info.save()
            # Save the updated parameter data
            ap.save_parameters(force=True)

    if update_sync_projects:
        # Load information from 'projects.info'
        project_metadata = ap.load_project_metadata(
            ap.params.project_metadata)
        # Comment out projects which don't exist on filesystem
        save_required = False
        for line in project_metadata:
            # Iterate through the named projects
            name = line['Project']
            if name.startswith('#'):
                # Commented out, ignore
                continue
            # Look for a matching project directory
            project_dir = os.path.join(ap.analysis_dir, name)
            if not os.path.exists(project_dir):
                print(f"Commenting out missing project '{name}'")
                line['Project'] = f"#{name}"
                save_required = True
        if save_required:
            project_metadata.save()
        # Add any project directories without entries
        save_required = False
        projects = [line['Project'] for line in project_metadata]
        for project in ap.get_analysis_projects_from_dirs():
            if project.name.endswith(".bak") \
                    or project.name.endswith(".orig") \
                    or project.name.endswith(".tmp"):
                # Skip directories with extensions indicating they
                # should be ignored
                print(f"Not adding entry for unlisted project '{project.name}'")
                continue
            elif project.name == "undetermined":
                # Skip undetermined
                continue
            elif project.name not in projects and f"#{project.name}" not in projects:
                # Add new entry
                print(f"Adding entry for unlisted project '{project.name}'")
                project_metadata.add_project(project.name,
                                             [s.name for s in project.samples],
                                             user=project.info.user,
                                             PI=project.info.PI,
                                             organism=project.info.organism,
                                             library_type=project.info.library_type,
                                             sc_platform=project.info.single_cell_platform,
                                             comments=project.info.comments)
                save_required = True
        # Save the updated project metadata if required
        if save_required:
            project_metadata.save()

    if update_project_metadata:
        # Update project metadata
        project_metadata = ap.load_project_metadata(
            ap.params.project_metadata)
        save_required = False
        for line in project_metadata:
            # Iterate through the named projects
            name = line['Project']
            if name.startswith('#'):
                # Commented out, ignore
                continue
            # Look for a matching project directory
            project_dir = os.path.join(ap.analysis_dir,name)
            if os.path.exists(project_dir):
                project = AnalysisProject(project_dir)
                print("Checking metadata for project '%s'" % name)
                # Synchronise metadata in projects with projects.info
                metadata_items = dict(
                    name=name,
                    user=line['User'],
                    PI=line['PI'],
                    organism=line['Organism'],
                    library_type=line['Library'],
                    single_cell_platform=line['SC_Platform'],
                    comments=line['Comments'],
                    samples=project.sample_summary()
                )
                # Only update items where values differ
                project_metadata_updated = False
                for item in metadata_items:
                    new_value = (metadata_items[item]
                                 if metadata_items[item] != '.' else None)
                    if project.info[item] != new_value:
                        print("...updating '%s' => %r" % (item,new_value))
                        project.info[item] = new_value
                        project_metadata_updated = True
                # Check paired-end info
                project_info = AnalysisProjectInfo(project.info_file)
                if project_info.paired_end != project.info.paired_end:
                    print("...updating paired-end info")
                    project_metadata_updated = True
                # Save the updated project metadata if required
                if project_metadata_updated:
                    print("...saving project metadata")
                    project.info.save()
                # Update list of sample names in projects.info
                sample_list = ','.join(sort_sample_names(
                    [s.name for s in project.samples]))
                if line['Samples'] != sample_list:
                    print("...updating sample list in projects.info")
                    line['Samples'] = sample_list
                    save_required = True
        # Save master project metadata
        if save_required:
            print("Saving projects.info")
            project_metadata.save()

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
