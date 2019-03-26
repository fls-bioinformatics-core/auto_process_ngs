#!/usr/bin/env python
#
#     publish_qc_cmd.py: implement auto process publish_qc command
#     Copyright (C) University of Manchester 2017-2019 Peter Briggs
#
#########################################################################

#######################################################################
# Imports
#######################################################################

import os
import time
import string
import ast
import auto_process_ngs.fileops as fileops
from auto_process_ngs.simple_scheduler import SimpleScheduler
from auto_process_ngs.qc.utils import verify_qc
from auto_process_ngs.qc.utils import report_qc
from ..docwriter import Document
from ..docwriter import Table
from ..docwriter import Link
from ..docwriter import Para
import bcftbx.utils as bcf_utils
from auto_process_ngs import get_version

# Module specific logger
import logging
logger = logging.getLogger(__name__)

#######################################################################
# Module data
#######################################################################

CSS_RULES = \
"""
h1       { background-color: #42AEC2;
           color: white;
           padding: 5px 10px; }
table    { margin: 10 10;
           border: solid 1px grey;
           background-color: white; }
th       { background-color: grey;
           color: white;
           padding: 2px 5px; }
td       { text-align: left;
           vertical-align: top;
           padding: 2px 5px;
           border-bottom: solid 1px lightgray; }
td.param { background-color: grey;
           color: white;
           padding: 2px 5px;
           font-weight: bold; }
div.footer { font-style: italic;
             font-size: 70%; }
"""

#######################################################################
# Command functions
#######################################################################

def publish_qc(ap,projects=None,location=None,ignore_missing_qc=False,
               regenerate_reports=False,force=False,use_hierarchy=False,
               exclude_zip_files=False,legacy=False):
    """
    Copy the QC reports to the webserver

    Looks for and copies various QC reports and outputs to
    a 'QC server' directory, and generates an HTML index.

    The reports include:

    - processing QC report
    - 'cellranger mkfastq' QC report
    - barcode analysis report

    Also if the analysis includes project directories then for
    each Fastq set in each project:

    - QC report for standard QC

    Also if a project comprises ICell8 or 10xGenomics Chromium
    data:

    - ICell8 processing report, or
    - 'cellranger count' reports for each sample

    In 'legacy' mode, the top-level report will also contain
    explicit links for each project for the following (where
    appropriate):

    - MultiQC report

    (These reports should now be accessible from the per-project
    QC reports, regardless of whther 'legacy' mode is specified.)

    Raises an exception if:

    - 'source' and 'run_number' metadata items are not set
    - a subset of projects don't have associated QC outputs
      (unless 'ignore_missing_qc' is True)

    Arguments:
      ap (AutoProcessor): autoprocessor pointing to the
        analysis directory to publish QC for
      projects (str): specify a glob-style pattern to match one
        or more projects to publish the reports for (default is
        to publish all reports)
      location (str): override the target location specified in
        the settings; can be of the form '[[user@]server:]directory'
      ignore_missing_qc (bool): if True then skip directories
        with missing QC data or reports (default is to raise
        an exception if projects have missing QC)
      regenerate_reports (bool): if True then try to create
        reports even when they already exist (default is to
        use existing reports)
      force (bool): if True then force QC report (re)generation
        even if QC is unverified (default is to raise an
        exception if projects cannot be verified)
      use_hierarchy (bool): if True then publish to a
        YEAR/PLATFORM subdirectory under the target location
        (default is not to use the hierarchy)
      exclude_zip_files (bool): if True then exclude any ZIP
        archives from publication (default is to include ZIP
        files)
      legacy (bool): if True then operate in 'legacy' mode (i.e.
        explicitly include MultiQC reports for each project)
    """
    # Turn off saving of parameters etc
    ap._save_params = False
    ap._save_metadata = False
    # Check metadata
    check_metadata = ap.check_metadata(('source','run_number'))
    if not check_metadata:
        raise Exception("Some metadata items not set, stopping")
    # Process pattern matching
    if projects is None:
        project_pattern = '*'
    else:
        project_pattern = projects
    # Get location to publish qc reports to
    if location is None:
        location = fileops.Location(ap.settings.qc_web_server.dirn)
    else:
        location = fileops.Location(location)
    if use_hierarchy:
        datestamp = str(ap.metadata.instrument_datestamp)
        if len(datestamp) == 6:
            # Assume YYMMDD datestamp format
            year = "20%s" % datestamp[0:2]
        elif len(datestamp) == 8:
            # Assume YYYYMMDD datestamp format
            year = datestamp[0:4]
        else:
            raise Exception("Invalid datestamp '%s' (use "
                            "--year option)" % datestamp)
        platform = ap.metadata.platform
    # Check the settings
    if location.is_remote:
        print "Copying QC to remote directory"
        print "user:\t%s" % location.user
        print "host:\t%s" % location.server
        print "dirn:\t%s" % location.path
    else:
        print "Copying QC to local directory"
    if use_hierarchy:
        print "dirn\t%s" % os.path.join(location.path,
                                        year,platform)
    else:
        print "dirn:\t%s" % location.path
    if exclude_zip_files:
        print "ZIP archives will be excluded from publication"
    if location.path is None:
        raise Exception("No target directory specified")
    if not fileops.exists(str(location)):
        raise Exception("Target directory '%s' doesn't exist" %
                        location)
    # Collect processing statistics
    print "Checking for processing QC report"
    processing_qc_html = os.path.join(ap.analysis_dir,
                                      "processing_qc.html")
    if os.path.exists(processing_qc_html):
        print "...found %s" % os.path.basename(processing_qc_html)
    else:
        print "...no processing QC report found"
        processing_qc_html = None
    # Collect barcode analysis artefacts
    print "Checking for barcode analysis"
    if ap.params.barcode_analysis_dir is not None:
        barcode_analysis_dir = ap.params.barcode_analysis_dir
    else:
        barcode_analysis_dir = 'barcode_analysis'
    if not os.path.isabs(barcode_analysis_dir):
        barcode_analysis_dir = os.path.join(ap.analysis_dir,
                                            barcode_analysis_dir)
    barcodes_files = []
    if os.path.exists(barcode_analysis_dir):
        print "...found barcode analysis dir: %s" % barcode_analysis_dir
        for filen in ('barcodes.report',
                      'barcodes.xls',
                      'barcodes.html'):
            filen = os.path.join(barcode_analysis_dir,filen)
            if os.path.exists(filen):
                print "...found %s" % os.path.basename(filen)
                barcodes_files.append(filen)
    if not barcodes_files:
        print "...no barcode analysis found"
    # Collect 10xGenomics cellranger QC summaries
    print "Checking for 10xGenomics cellranger QC summaries"
    cellranger_qc_html = []
    for filen in os.listdir(ap.analysis_dir):
        if filen.startswith("cellranger_qc_summary") and \
           filen.endswith(".html"):
            print "...found %s" % filen
            cellranger_qc_html.append(
                os.path.join(ap.analysis_dir,filen))
    if not cellranger_qc_html:
        print "...no cellranger QC summaries found"
    # Collect QC for project directories
    print "Checking project directories"
    projects = ap.get_analysis_projects_from_dirs(pattern=project_pattern)
    if projects:
        ap.set_log_dir(ap.get_log_subdir('publish_qc'))
    project_qc = {}
    for project in projects:
        # Check qc subdirectories
        print "Checking project '%s':" % project.name
        project_qc[project.name] = bcf_utils.AttributeDictionary()
        project_qc[project.name]['qc_dirs'] = {}
        if not project.qc_dirs:
            print "...no QC directories found"
        for qc_dir in project.qc_dirs:
            print "...found QC dir '%s'" % qc_dir
            # Gather the QC artefacts
            qc_artefacts = bcf_utils.AttributeDictionary()
            project_qc[project.name].qc_dirs[qc_dir] = qc_artefacts
            # Base name for QC reports
            qc_base = "%s_report" % qc_dir
            # Set the source Fastqs dir
            fastq_dir = project.qc_info(qc_dir).fastq_dir
            project.use_fastq_dir(fastq_dir)
            print "...associated Fastq set '%s'" % \
                os.path.basename(fastq_dir)
            qc_protocol = project.qc_info(qc_dir).protocol
            print "...associated QC protocol '%s'" % qc_protocol
            if qc_protocol is None:
                qc_protocol = "standardPE"
                print "...assuming QC protocol '%s'" % qc_protocol
            # Verify the QC and check for report
            verified = verify_qc(
                project,
                fastq_dir=fastq_dir,
                qc_dir=os.path.join(project.dirn,qc_dir),
                qc_protocol=qc_protocol,
                log_dir=ap.log_dir)
            if verified:
                print "...%s: verified QC" % qc_dir
                # Check for an existing report
                qc_zip = os.path.join(
                    project.dirn,
                    "%s_report.%s.%s.zip" %
                    (qc_dir,project.name,
                     os.path.basename(ap.analysis_dir)))
                # Check if we need to (re)generate report
                if (regenerate_reports or
                    not os.path.exists(qc_zip)):
                    report_status = report_qc(project,
                                              qc_dir=qc_dir,
                                              multiqc=True,
                                              log_dir=ap.log_dir)
                    if report_status == 0:
                        print "...%s: (re)generated report" % qc_dir
                    else:
                        print "...%s: failed to (re)generate " \
                            "QC report" % qc_dir
                # Add to the list of verified QC dirs
                if os.path.exists(qc_zip):
                    qc_artefacts['qc_zip'] = qc_zip
                else:
                    print "...%s: missing QC report" % qc_dir
            else:
                # Not verified
                print "...%s: failed to verify QC" % qc_dir
            if legacy:
                # MultiQC report
                multiqc_report = os.path.join(project.dirn,
                                              "multi%s.html"
                                              % qc_base)
                if os.path.exists(multiqc_report):
                    print "...%s: found MultiQC report" % qc_dir
                    qc_artefacts['multiqc_report'] = multiqc_report
                else:
                    print "...%s: no MultiQC report" % qc_dir
        # ICell8 pipeline report
        icell8_zip = os.path.join(project.dirn,
                                  "icell8_processing.%s.%s.zip" %
                                  (project.name,
                                   os.path.basename(ap.analysis_dir)))
        if os.path.exists(icell8_zip):
            print "...%s: found ICell8 pipeline report" % project.name
            project_qc[project.name]['icell8_zip'] = icell8_zip
        # Cellranger count report
        cellranger_zip = os.path.join(project.dirn,
                                      "cellranger_count_report.%s.%s.zip" %
                                      (project.name,
                                       os.path.basename(ap.analysis_dir)))
        if os.path.exists(cellranger_zip):
            print "...%s: found cellranger count report" % project.name
            project_qc[project.name]['cellranger_zip'] = cellranger_zip
    # Projects with no QC
    no_qc_projects = filter(lambda p: not project_qc[p.name].qc_dirs,
                            projects)
    # Determine what projects are left and if we can proceed
    if no_qc_projects:
        # Failed to generate results for some projects
        err_msg = "No QC reports for projects: %s" % \
                  ', '.join([x.name for x in no_qc_projects])
        if not ignore_missing_qc:
            # Fatal error
            raise Exception(err_msg)
        # Proceed with a warning
        logging.warning(err_msg)
    # Remove the 'bad' projects from the list before proceeding
    for project in no_qc_projects:
        print "Project %s will be skipped" % project.name
        projects.remove(project)
    if not projects:
        logging.warning("No projects with QC results to publish")
    # Set up the final destination
    if use_hierarchy:
        dirn = os.path.join(str(location),year,platform)
    else:
        dirn = str(location)
    dirn = os.path.join(dirn,os.path.basename(ap.analysis_dir))
    # Create the publication directory
    fileops.mkdir(dirn,recursive=True)
    if not fileops.exists(dirn):
        raise Exception("Failed to create directory: %s" % dirn)
    # Do file transfer and unpacking
    if projects:
        # Make log directory and set up scheduler
        # to farm out the intensive operations to
        ap.set_log_dir(ap.get_log_subdir('publish_qc'))
        runner = ap.settings.general.default_runner
        runner.set_log_dir(ap.log_dir)
        sched = SimpleScheduler(
            runner=runner,
            max_concurrent=ap.settings.general.max_concurrent_jobs,
            poll_interval=ap.settings.general.poll_interval
        )
        sched.start()
    else:
        sched = None
    # Start building an index page
    title = "QC reports for %s" % os.path.basename(ap.analysis_dir)
    index_page = Document(title=title)
    index_page.add_css_rule(CSS_RULES)
    # General info
    general_info = index_page.add_section("General information")
    params_tbl = Table(("param","value",))
    params_tbl.no_header()
    params_tbl.add_css_classes("param",column="param")
    params_tbl.add_row(param="Run number",
                       value=ap.metadata.run_number)
    params_tbl.add_row(param="Platform",
                       value=ap.metadata.platform)
    params_tbl.add_row(param="Endedness",
                       value=('Paired end' if ap.paired_end
                              else 'Single end'))
    try:
        bcl2fastq_software = ast.literal_eval(
            ap.metadata.bcl2fastq_software)
    except ValueError:
        bcl2fastq_software = None
    params_tbl.add_row(param="Bcl2fastq",
                       value=('unspecified' if not bcl2fastq_software else
                              "%s %s" % (bcl2fastq_software[1],
                                         bcl2fastq_software[2])))
    try:
        cellranger_software = ast.literal_eval(
            ap.metadata.cellranger_software)
    except ValueError:
        cellranger_software = None
    params_tbl.add_row(param="Cellranger",
                       value=('unspecified' if not cellranger_software else
                              "%s %s" % (cellranger_software[1],
                                         cellranger_software[2])))
    params_tbl.add_row(param="Reference",
                       value=ap.run_reference_id)
    general_info.add(params_tbl)
    # Processing QC report
    if processing_qc_html:
        fileops.copy(processing_qc_html,dirn)
        processing_stats = index_page.add_section("Processing Statistics")
        processing_stats.add(Link("Processing QC report",
                                  os.path.basename(processing_qc_html)))
    # Barcode analysis
    if barcodes_files:
        # Create section
        barcodes = index_page.add_section("Barcode Analysis")
        barcodes.add(
            Para("Plain text report:",
                 Link("barcodes.report","barcodes/barcodes.report"),
                 "| XLS:",
                 Link("barcodes.xls","barcodes/barcodes.xls"),
                 "| HTML:",
                 Link("barcodes.html","barcodes/barcodes.html")))
        # Create subdir and copy files
        barcodes_dirn = os.path.join(dirn,'barcodes')
        fileops.mkdir(barcodes_dirn)
        for filen in barcodes_files:
            fileops.copy(filen,barcodes_dirn)
    # 10xGenomics cellranger QC summaries
    if cellranger_qc_html:
        cellranger_qc = index_page.add_section(
            "QC summary: cellranger mkfastq")
        for qc_html in cellranger_qc_html:
            # Check for optional lane list at tail of QC summary
            # e.g. cellranger_qc_summary_45.html
            # This might not be present
            lanes = qc_html.split('.')[0].split('_')[-1]
            if all(c in string.digits for c in lanes):
                lanes = ','.join(lanes)
            else:
                lanes = None
            fileops.copy(qc_html,dirn)
            cellranger_qc.add(Link("QC summary for %s" %
                                   ("all lanes" if lanes is None
                                    else "lanes %s" % lanes),
                                   os.path.basename(qc_html)))
    if projects:
        # Table of projects
        projects_summary = index_page.add_section("QC Reports")
        projects_tbl = Table(("Project",
                              "User",
                              "Library",
                              "Organism",
                              "PI",
                              "Samples",
                              "Nsamples",
                              "Reports"))
        projects_summary.add(projects_tbl)
        # Set the string to represent "null" table entries
        null_str = '&nbsp;'
        # Deal with QC for each project
        for project in projects:
            # Reset source Fastqs dir
            project.use_fastq_dir()
            # Get local versions of project information
            info = project.info
            project_user = null_str if info.user is None else info.user
            if info.library_type is not None:
                library_type = info.library_type
                if info.single_cell_platform is not None:
                    library_type += " (%s)" % info.single_cell_platform
            else:
                library_type = null_str
            organism = null_str if info.organism is None else info.organism
            PI = null_str if info.PI is None else info.PI
            # Generate line in the table of projects
            idx = projects_tbl.add_row(Project=project.name,
                                       User=project_user,
                                       Library=library_type,
                                       Organism=organism,
                                       PI=PI,
                                       Samples=project.prettyPrintSamples(),
                                       Nsamples=len(project.samples))
            # Copy QC reports and other artefacts
            report_html = Para()
            for qc_dir in project_qc[project.name].qc_dirs:
                qc_artefacts = project_qc[project.name].qc_dirs[qc_dir]
                qc_base = "%s_report" % qc_dir
                fastq_dir = project.qc_info(qc_dir).fastq_dir
                if fastq_dir.startswith("%s%s" % (project.dirn,os.sep)):
                    fastq_dir = os.path.relpath(fastq_dir,project.dirn)
                if fastq_dir != project.info.primary_fastq_dir:
                    fastq_set = fastq_dir
                else:
                    fastq_set = None
                # QC report
                try:
                    qc_zip = qc_artefacts.qc_zip
                    try:
                        # Copy and unzip QC report
                        copy_job = sched.submit(
                            fileops.copy_command(qc_zip,dirn),
                            name="copy.qc_report.%s.%s" % (project.name,
                                                           qc_dir))
                        unzip_job = sched.submit(
                            fileops.unzip_command(os.path.join(
                                dirn,
                                os.path.basename(qc_zip)),
                                                  fileops.Location(dirn).path),
                            name="unzip.qc_report.%s.%s" % (project.name,
                                                            qc_dir),
                            wait_for=(copy_job.name,))
                        sched.wait()
                        # Check the jobs completed ok
                        if copy_job.exit_status or unzip_job.exit_status:
                            raise Exception("copy and/or unzip job failed")
                        if exclude_zip_files:
                            print "Removing %s from server" % qc_zip
                            fileops.remove_file(os.path.join(
                                dirn,
                                os.path.basename(qc_zip)))
                        # Append info to the index page
                        fastq_set_name = (" (%s)" % fastq_set
                                          if fastq_set is not None
                                          else "")
                        report_html.add(
                            Link("[Report%s]" % fastq_set_name,
                                 "%s.%s.%s/%s.html"
                                 % (qc_base,
                                    project.name,
                                    os.path.basename(ap.analysis_dir),
                                    qc_base)))
                        if not exclude_zip_files:
                            report_html.add(
                                Link("[ZIP%s]" % fastq_set_name,
                                     os.path.basename(qc_zip)))
                    except Exception as ex:
                        print "Failed to copy QC report: %s" % ex
                except AttributeError:
                    # No QC report
                    pass
                # MultiQC
                try:
                    multiqc_report = qc_artefacts.multiqc_report
                    assert(os.path.isfile(multiqc_report))
                    final_multiqc = "multi%s.%s.html" % (qc_base,
                                                         project.name)
                    try:
                        fileops.copy(multiqc_report,
                                     os.path.join(dirn,final_multiqc))
                        report_html.add(
                            Link("[MultiQC%s]" % fastq_set_name,
                                 final_multiqc))
                    except Exception as ex:
                        print "Failed to copy MultiQC report: %s" % ex
                except AttributeError:
                    # No MultiQC report
                    pass
            # Check there is something to add
            if not report_html:
                report_html.add("QC reports not available")
            # ICell8 pipeline report
            try:
                icell8_zip = project_qc[project.name].icell8_zip
                try:
                    # Copy and unzip ICell8 report
                    copy_job = sched.submit(
                        fileops.copy_command(icell8_zip,dirn),
                        name="copy.icell8_report.%s" % project.name)
                    unzip_job = sched.submit(
                        fileops.unzip_command(os.path.join(
                            dirn,
                            os.path.basename(icell8_zip)),
                                              fileops.Location(dirn).path),
                        name="unzip.icell8_report.%s" % project.name,
                        wait_for=(copy_job.name,))
                    sched.wait()
                    # Check the jobs completed ok
                    if copy_job.exit_status or unzip_job.exit_status:
                        raise Exception("copy and/or unzip job failed")
                    if exclude_zip_files:
                        print "Removing %s from server" % icell8_zip
                        fileops.remove_file(os.path.join(
                            dirn,
                            os.path.basename(icell8_zip)))
                    # Append info to the index page
                    report_html.add(
                        Link("[Icell8 processing]",
                             "icell8_processing.%s.%s/"
                             "icell8_processing.html" %
                             (project.name,
                              os.path.basename(ap.analysis_dir))))
                    if not exclude_zip_files:
                        report_html.add(
                            Link("[ZIP]",
                                 os.path.basename(icell8_processing_zip)))
                except Exception as ex:
                    print "Failed to copy ICell8 report: %s" % ex
            except AttributeError:
                # No ICell8 report
                pass
            # Cellranger count reports
            try:
                cellranger_zip = project_qc[project.name].cellranger_zip
                try:
                    # Copy and unzip cellranger report
                    copy_job = sched.submit(
                        fileops.copy_command(cellranger_zip,dirn),
                        name="copy.cellranger_report.%s" % project.name)
                    unzip_job = sched.submit(
                        fileops.unzip_command(os.path.join(
                            dirn,
                            os.path.basename(cellranger_zip)),
                                              fileops.Location(dirn).path),
                        name="unzip.cellranger_report.%s" % project.name,
                        wait_for=(copy_job.name,))
                    sched.wait()
                    # Check the jobs completed ok
                    if copy_job.exit_status or unzip_job.exit_status:
                        raise Exception("copy and/or unzip job failed")
                    if exclude_zip_files:
                            print "Removing %s from server" % cellranger_zip
                            fileops.remove_file(os.path.join(
                                dirn,
                                os.path.basename(cellranger_zip)))
                    # Append info to the index page
                    report_html.add(
                        Link("[Cellranger count]",
                             "cellranger_count_report.%s.%s/"
                             "cellranger_count_report.html" %
                             (project.name,
                              os.path.basename(ap.analysis_dir))))
                    if not exclude_zip_files:
                        report_html.add(
                            Link("[ZIP]",
                                 os.path.basename(cellranger_zip)))
                except Exception as ex:
                    print "Failed to copy cellranger report: %s" % ex
            except AttributeError:
                # No cellranger count data to copy
                pass
            # Add to the index
            projects_tbl.set_value(idx,"Reports",report_html)
    # Finish index page
    footer = index_page.add_section(css_classes=("footer",))
    footer.add("Generated by auto_process.py %s on %s" %
               (get_version(),time.asctime()))
    # Copy to server
    index_html = os.path.join(ap.tmp_dir,'index.html')
    print "Writing index page: %s" % index_html
    index_page.write(index_html)
    fileops.copy(index_html,dirn)
    # Stop scheduler
    if sched is not None:
        sched.stop()
    # Print the URL if given
    if ap.settings.qc_web_server.url is not None:
        url = ap.settings.qc_web_server.url
        if use_hierarchy:
            url = os.path.join(url,
                               year,
                               platform,
                               os.path.basename(ap.analysis_dir),
                               "index.html")
        print "QC published to %s" % url
