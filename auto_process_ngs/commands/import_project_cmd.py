#!/usr/bin/env python
#
#     import_project_cmd.py: implement auto process import_project command
#     Copyright (C) University of Manchester 2019 Peter Briggs
#
#########################################################################

#######################################################################
# Imports
#######################################################################

import os
import logging
from ..applications import general as applications_general
from ..analysis import AnalysisProject
from ..qc.utils import verify_qc
from ..qc.utils import report_qc

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Command functions
#######################################################################

def import_project(ap,project_dir):
    """
    Import a project directory into this analysis directory

    Arguments:
      ap (AutoProcessor): autoprocessor pointing to the parent
        analysis directory
      project_dir (str): path to project directory to be
        imported
    """
    # Check that target directory exists
    project_dir = os.path.abspath(project_dir)
    # Check that project doesn't already exist
    project_name = os.path.basename(project_dir)
    project_metadata = ap.load_project_metadata()
    if project_name in [p['Project'] for p in project_metadata] or \
       AnalysisProject(project_name,
                       os.path.join(ap.analysis_dir,
                                    project_name)).exists:
        raise Exception("Project called '%s' already exists" %
                        project_name)
    # Load target as a project
    project = AnalysisProject(project_name,project_dir)
    # Set up a log directory
    ap.set_log_dir(ap.get_log_subdir('import_project_%s' % project_name))
    # Rsync the project directory
    print "Importing project directory contents for '%s'" % project_name
    try:
        excludes = ['--exclude=tmp.*',
                    '--exclude=qc_report.*']
        rsync = applications_general.rsync(project_dir,
                                           ap.analysis_dir,
                                           extra_options=excludes)
        print "Running %s" % rsync
        status = rsync.run_subprocess(log=ap.log_path('import_project.rsync.log'))
    except Exception as ex:
        logger.error("Exception importing project: %s" % ex)
        raise ex
    if status != 0:
        raise Exception("Failed to import project from %s (status %s)" %
                        (project_dir,status))
    # Update the projects.info metadata file
    print "Updating projects.info file with imported project"
    project_metadata = ap.load_project_metadata()
    sample_names = [s.name for s in project.samples]
    project_metadata.add_project(project_name,
                                 sample_names,
                                 user=project.info.user,
                                 library_type=project.info.library_type,
                                 single_cell_platform=project.info.single_cell_platform,
                                 organism=project.info.organism,
                                 PI=project.info.PI,
                                 comments=project.info.comments)
    project_metadata.save()
    # Report
    print "Projects now in metadata file:"
    for p in project_metadata:
        print "- %s" % p['Project']
    # Update the QC report
    try:
        project = ap.get_analysis_projects(pattern=project_name)[0]
    except Exception as ex:
        logger.error("Exception when trying to acquire project %s: %s"
                      % (project_name,ex))
        return
    if project.qc_dir is None:
        print "No QC directory assigned for %s" % project_name
    else:
        print "QC directory for %s: %s" % (project_name,
                                           project.qc_dir)
        qc_info = project.qc_info(project.qc_dir)
        if qc_info.fastq_dir:
            print "Updating stored Fastq directory for QC"
            fastq_dir = os.path.join(project.dirn,
                                     os.path.relpath(qc_info.fastq_dir,
                                                     project_dir))
            print "Updated Fastq directory: %s" % fastq_dir
            qc_info['fastq_dir'] = fastq_dir
            qc_info.save()
        if verify_qc(project,log_dir=ap.log_dir):
            try:
                report_qc(project,
                          zip_outputs=True,
                          multiqc=True,
                          log_dir=ap.log_dir)
                print "Updated QC report for %s" % project_name
            except Exception as ex:
                raise Exception("Project '%s' imported but failed to "
                                "generate QC report: %s" % (project_name,
                                                            ex))
        else:
            print "Failed to verify QC: report not updated"
