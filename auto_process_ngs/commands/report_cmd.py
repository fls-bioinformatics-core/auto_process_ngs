#!/usr/bin/env python
#
#     report_cmd.py: implement auto process 'report' command
#     Copyright (C) University of Manchester 2018 Peter Briggs
#
#########################################################################

#######################################################################
# Imports
#######################################################################

import os
import ast
import logging
import bcftbx.IlluminaData as IlluminaData
import bcftbx.utils as bcf_utils
import auto_process_ngs.utils as utils

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Constants
#######################################################################

class ReportingMode(object):
    CONCISE = 0
    SUMMARY = 1
    PROJECTS = 2
    INFO = 4

#######################################################################
# Command functions
#######################################################################

def report(ap,mode=None):
    """
    Print a report on an analysis project

    Arguments:
      ap (AutoProcessor): autoprocessor pointing to the
        analysis directory to be reported on
      mode (int): reporting mode (concise, summary,
        projects or info)
    """
    # Turn off saving of parameters
    ap._save_params = False
    ap._save_metadata = False
    # Do the reporting
    if mode is None or mode == ReportingMode.INFO:
        f = report_info
    elif mode == ReportingMode.CONCISE:
        f = report_concise
    elif mode == ReportingMode.SUMMARY:
        f = report_summary
    elif mode == ReportingMode.PROJECTS:
        f = report_projects
    else:
        raise Exception("Unknown reporting mode")
    print f(ap)

def report_info(ap):
    """Generate a general report

    Generates an unstructured report on the contents
    of the analysis directory.

    Arguments:
      ap (AutoProcessor): autoprocessor pointing to the
        analysis directory to be reported on
        
    Returns:
      String with the report text.
    """
    report = []
    report.append("Run reference: %s" % ap.run_reference_id)
    report.append("Directory    : %s" % ap.analysis_dir)
    report.append("Platform     : %s" % ap.metadata.platform)
    report.append("Unaligned dir: %s" % ap.params.unaligned_dir)
    if ap.readme_file:
        report.append("README.txt found: %s" % ap.readme_file)
    if ap.params.unaligned_dir is not None or \
       not os.path.exists(ap.params.unaligned_dir):
        try:
            illumina_data = ap.load_illumina_data()
            report.append("\nSummary of data in '%s' dir:\n" %
                          ap.params.unaligned_dir)
            for project in illumina_data.projects:
                report.append("- %s" % IlluminaData.describe_project(project))
        except IlluminaData.IlluminaDataError as ex:
            report.append("Failed to load data from %s:" % ap.params.unaligned_dir)
            report.append("%s" % ex)
    else:
        report.append("No information on source fastq data (no unaligned dir "
                      "found)")
    projects = ap.get_analysis_projects()
    report.append("\n%d analysis project%s:" % (len(projects),
                                                "s" if len(projects) != 0
                                                else ""))
    for project in projects:
        info = project.info
        report.append("\n- %s" % project.name)
        report.append("  %s" % ('-'*len(project.name),))
        report.append("  User    : %s" % info.user)
        report.append("  PI      : %s" % info.PI)
        report.append("  Library : %s" % info.library_type)
        report.append("  SC Plat.: %s" % info.single_cell_platform)
        report.append("  Organism: %s" % info.organism)
        report.append("  Dir     : %s" % os.path.basename(project.dirn))
        report.append("  #samples: %s" % len(project.samples))
        report.append("  #cells  : %s" % default_value(info.number_of_cells))
        report.append("  Samples : %s" % project.prettyPrintSamples())
        report.append("  QC      : %s" % ('ok'
                                          if project.verify_qc()
                                          else 'not verified'))
        report.append("  Comments: %s" % (project.info.comments))
    return '\n'.join(report)

def report_concise(ap):
    """Generate one-line report suitable for logging

    Generates a one-line report that can be used for logging, for
    example:

    Paired end: 'PJB': Peter Briggs, Mouse ChIP-seq (PI: P Briggs) (6 samples); ...

    The report is based on project directories that are located
    in the analysis directory, and not from other information
    (e.g. contents of the 'bcl2fastq' directory or other outputs
    from processing).

    Arguments:
      ap (AutoProcessor): autoprocessor pointing to the
        analysis directory to be reported on

    Returns:
      String with the report text.
    """
    report = []
    analysis_dir = utils.AnalysisDir(ap.analysis_dir)
    for p in analysis_dir.projects:
        samples = "%d sample%s" % (len(p.samples),
                                   's' if len(p.samples) != 1
                                   else '')
        if p.info.number_of_cells is not None:
            samples += "/%d cell%s" % (p.info.number_of_cells,
                                       's' if p.info.number_of_cells != 1
                                       else '')
        report.append("'%s': %s, %s %s%s (PI: %s) (%s)" % \
                      (p.name,
                       p.info.user,
                       p.info.organism,
                       ('%s ' % p.info.single_cell_platform
                        if p.info.single_cell_platform else ''),
                       p.info.library_type,
                       p.info.PI,
                       samples
                   ))
    report = '; '.join(report)
    # Paired end run?
    if analysis_dir.paired_end:
        endedness = "Paired end"
    else:
        endedness = "Single end"
    report = "%s: %s" % (endedness,report)
    return report

def report_summary(ap):
    """Generate summary report suitable for bioinformaticians

    Generates a multi-line report which gives general information
    about the run, plus one-line summaries for each project, plus
    any additional information that has been recorded.

    The general information includes:

    - Platform
    - Run name
    - Run reference id
    - Processing software
    - Assay (i.e. sequencing kit)

    For each project:

    - Project subdirectory
    - Researcher (aka user)
    - PI
    - Application (aka library type)
    - Single cell prep platform (e.g. ICell8)
    - Organism
    - Number of samples

    Arguments:
      ap (AutoProcessor): autoprocessor pointing to the
        analysis directory to be reported on

    Returns:
      String with the report text.
    """
    # Gather information
    analysis_dir = utils.AnalysisDir(ap.analysis_dir)
    datestamp = None
    instrument = None
    run_number = None
    run_name = ap.run_name
    try:
        datestamp,instrument,run_number = IlluminaData.split_run_name(run_name)
    except Exception, ex:
        logger.warning("Unable to extract information from run name '%s'" \
                       % run_name)
        logger.warning("Exception: %s" % ex)
    if ap.metadata.platform is not None:
        platform = ap.metadata.platform.upper()
    else:
        platform = 'unknown'
    if ap.metadata.run_number is not None:
        run_number = ap.metadata.run_number
    try:
        bcl2fastq_software = ast.literal_eval(
            ap.metadata.bcl2fastq_software)
    except ValueError:
        bcl2fastq_software = None
    if ap.metadata.assay is not None:
        assay = ap.metadata.assay
    else:
        assay = ''
    # Generate report text
    report = []
    # Report header
    if datestamp and instrument and run_number:
        title = "%s run #%s datestamped %s" % (platform,
                                               run_number,
                                               datestamp)
    else:
        title = "%s" % os.path.basename(ap.analysis_dir)
    report.append("%s\n%s" % (title,'='*len(title)))
    report.append("Run name : %s" % run_name)
    report.append("Reference: %s" % ap.run_reference_id)
    report.append("Platform : %s" % platform)
    report.append("Directory: %s" % ap.params.analysis_dir)
    report.append("Endedness: %s" % \
                  ('Paired end' if analysis_dir.paired_end else 'Single end'))
    report.append("Bcl2fastq: %s" %
                  ('Unknown' if not bcl2fastq_software else
                   "%s %s" % (bcl2fastq_software[1],
                              bcl2fastq_software[2])))
    report.append("Assay    : %s" % assay)
    report.append("")
    # Projects
    report.append("%d project%s:" % (analysis_dir.n_projects,
                                     '' if analysis_dir.n_projects == 1 else 's'))
    data_items = ('user',
                  'PI',
                  'library_type',
                  'single_cell_platform',
                  'number_of_cells',
                  'organism')
    rows = []
    comments = bcf_utils.OrderedDictionary()
    for project in analysis_dir.projects:
        project_data = dict(project=project.name)
        for item in data_items:
            value = project.info[item]
            project_data[item] = value if value not in ('.','?') else \
                                 '<unspecified %s>' % item.lower()
        library = project_data['library_type']
        if project_data['single_cell_platform'] is not None:
            library += " (%s)" % project_data['single_cell_platform']
        samples = "%d sample%s" % (len(project.samples),
                                   's' if len(project.samples) != 1 else '')
        if project_data['number_of_cells'] is not None:
            samples += "/%d cell%s" % (int(project_data['number_of_cells']),
                                       's' if int(project_data['number_of_cells']) != 1 else '')
        rows.append(("- '%s':" % project_data['project'],
                     project_data['user'],
                     project_data['organism'],
                     library,
                     samples,
                     "(PI %s)" % project_data['PI']))
        if project.info.comments:
            comments[project.name] = project.info.comments
    report.append(utils.pretty_print_rows(rows))
    # Additional comments/notes
    if comments:
        width = max([len(x) for x in comments])
        report.append("")
        report.append("Additional notes/comments:")
        for project in comments:
            first_line = True
            for line in bcf_utils.split_into_lines(comments[project],70-width):
                if first_line:
                    report.append("- %s%s: %s" % (project,
                                                  ' '*(width-len(project)),
                                                  line))
                    first_line = False
                else:
                    report.append("  %s  %s" % (' '*width,line))
    return '\n'.join(report)

def report_projects(ap):
    """Generate one line reports suitable for pasting into spreadsheet

    Generate one-line report for each each project with tab-separated
    data items, suitable for injection into a spreadsheet.

    Each line has the following information:

    - Run id e.g. HISEQ_140328
    - Run number
    - Source
    - Date
    - User
    - PI
    - Application
    - Single Cell Platform
    - Organism
    - Platform
    - #Samples
    - #Cells
    - PE (yes/no)
    - Samples

    Arguments:
      ap (AutoProcessor): autoprocessor pointing to the
        analysis directory to be reported on
        
    Returns:
      String with the report text.
    """
    # Acquire data
    analysis_dir = utils.AnalysisDir(ap.analysis_dir)
    # General information
    run_name = ap.run_name
    try:
        datestamp,instrument,run_number = IlluminaData.split_run_name(run_name)
        run_number = run_number.lstrip('0')
    except Exception, ex:
        logger.warning("Unable to extract information from run name '%s'" \
                       % run_name)
        logger.warning("Exception: %s" % ex)
        date_stamp = ''
        run_number = ''
    if ap.metadata.platform is not None:
        platform = ap.metadata.platform.upper()
    else:
        platform = ''
    run_id = ap.run_reference_id
    if ap.metadata.run_number is not None:
        run_number = ap.metadata.run_number
    if ap.metadata.source is not None:
        data_source = ap.metadata.source
    else:
        data_source = ''
    paired_end = 'yes' if analysis_dir.paired_end else 'no'
    # Generate report, one line per project
    report = []
    nprojects = len(analysis_dir.projects)
    if nprojects == 0:
        report.append("No projects found")
    else:
        report.append("%s project%s found" % (nprojects,
                                              ('' if nprojects == 1
                                               else 's')))
    for project in analysis_dir.projects:
        project_line = [run_id,str(run_number),data_source,'']
        info = project.info
        project_line.append('' if not info.user else info.user)
        project_line.append('' if not info.PI else info.PI)
        project_line.append('' if not info.library_type
                            else info.library_type)
        project_line.append('' if not info.single_cell_platform
                            else info.single_cell_platform)
        project_line.append('' if not info.organism else info.organism)
        project_line.append(platform)
        project_line.append(str(len(project.samples)))
        project_line.append('' if not info.number_of_cells
                            else str(info.number_of_cells))
        project_line.append(paired_end)
        project_line.append(project.prettyPrintSamples())
        report.append('\t'.join(project_line))
    report = '\n'.join(report)
    return report

def default_value(s,default=""):
    """
    Returns supplied value, or default if supplied value is None
    """
    if s is None:
        return default
    return s
