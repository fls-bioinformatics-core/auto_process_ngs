#!/usr/bin/env python
#
#     import_project_cmd.py: implement auto process import_project command
#     Copyright (C) University of Manchester 2019-2021 Peter Briggs
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
from ..simple_scheduler import SchedulerJob

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Command functions
#######################################################################

def import_project(ap,project_dir,runner=None):
    """
    Import a project directory into an analysis directory

    Importing a project directory consists of the following
    actions:

    - Copying the project directory and contents to the
      analysis directory
    - Updating 'projects.info' to add in the data about the
      imported project
    - Updating the project metadata and the QC report

    Arguments:
      ap (AutoProcessor): autoprocessor pointing to the parent
        analysis directory
      project_dir (str): path to project directory to be
        imported
      runner (JobRunner): explicitly specify the job runner to
        send jobs to (overrides default runner set in the
        configuration)
    """
    # Load target directory as a project
    project_dir = os.path.abspath(project_dir)
    project = AnalysisProject(project_dir)
    # Check that project doesn't already exist
    project_metadata = ap.load_project_metadata()
    if project.name in [p['Project'] for p in project_metadata] or \
       AnalysisProject(os.path.join(ap.analysis_dir,
                                    project.name)).exists:
        raise Exception("Project called '%s' already exists" %
                        project.name)
    # Set up job runner
    if runner is None:
        runner = ap.settings.general.default_runner
    print("Using job runner '%s'" % runner)
    # Set up a log directory
    ap.set_log_dir(ap.get_log_subdir('import_project_%s' % project.name))
    # Rsync the project directory
    print("Importing project directory contents for '%s'" % project.name)
    try:
        excludes = ['--exclude=tmp.*',
                    '--exclude=qc_report.*',
                    '--exclude=multi_qc*']
        rsync = applications_general.rsync(project_dir,
                                           ap.analysis_dir,
                                           extra_options=excludes)
        print("Running %s" % rsync)
        rsync_job = SchedulerJob(runner,
                                 rsync.command_line,
                                 name="import_project.rsync.%s" % project.name,
                                 log_dir=ap.log_dir)
        rsync_job.start()
        try:
            rsync_job.wait()
        except KeyboardInterrupt as ex:
            logger.warning("Keyboard interrupt, terminating project import")
            rsync_job.terminate()
            raise ex
        status = rsync_job.exit_status
    except Exception as ex:
        logger.error("Exception importing project: %s" % ex)
        raise ex
    if status != 0:
        raise Exception("Failed to import project from %s (status %s)" %
                        (project_dir,status))
    # Update the projects.info metadata file
    print("Updating projects.info file with imported project")
    project_metadata = ap.load_project_metadata()
    sample_names = [s.name for s in project.samples]
    project_metadata.add_project(
        project.name,
        sample_names,
        user=project.info.user,
        library_type=project.info.library_type,
        single_cell_platform=project.info.single_cell_platform,
        organism=project.info.organism,
        PI=project.info.PI,
        comments=project.info.comments)
    project_metadata.save()
    # Report
    print("Projects now in metadata file:")
    for p in project_metadata:
        print("- %s" % p['Project'])
    # Update the project metadata and QC report
    try:
        project = ap.get_analysis_projects(pattern=project.name)[0]
    except Exception as ex:
        logger.error("Exception when trying to acquire project %s: %s"
                      % (project.name,ex))
        return
    project.info['run'] = ap.run_name
    project.info.save()
    print("Set run to %s for project '%s'" % (project.info.run,
                                              project.name))
    if project.qc_dir is None:
        print("No QC directory assigned for %s" % project.name)
    else:
        print("QC directory for %s: %s" % (project.name,
                                           project.qc_dir))
        qc_info = project.qc_info(project.qc_dir)
        if qc_info.fastq_dir:
            print("Updating stored Fastq directory for QC")
            fastq_dir = os.path.join(project.dirn,
                                     os.path.relpath(qc_info.fastq_dir,
                                                     project_dir))
            print("Updated Fastq directory: %s" % fastq_dir)
            qc_info['fastq_dir'] = fastq_dir
            qc_info.save()
        if verify_qc(project,log_dir=ap.log_dir,runner=runner):
            try:
                status = report_qc(project,
                                   zip_outputs=True,
                                   multiqc=True,
                                   runner=runner,
                                   log_dir=ap.log_dir)
                if status != 0:
                    raise Exception("Report generation returned non-zero "
                                    "exit code (%s)" % status)
                print("Updated QC report for %s" % project.name)
            except Exception as ex:
                raise Exception("Project '%s' imported but failed to "
                                "generate QC report: %s" % (project.name,
                                                            ex))
        else:
            print("Failed to verify QC: report not updated")
