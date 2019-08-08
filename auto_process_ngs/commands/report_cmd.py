#!/usr/bin/env python
#
#     report_cmd.py: implement auto process 'report' command
#     Copyright (C) University of Manchester 2018-2019 Peter Briggs
#
#########################################################################

#######################################################################
# Imports
#######################################################################

import sys
import os
import ast
import logging
import tempfile
import bcftbx.IlluminaData as IlluminaData
import bcftbx.utils as bcf_utils
from .. import analysis
from .. import utils
from .. import fileops
from ..qc.utils import verify_qc

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

def report(ap,mode=None,fields=None,out_file=None):
    """
    Print a report on an analysis project

    Arguments:
      ap (AutoProcessor): autoprocessor pointing to the
        analysis directory to be reported on
      mode (int): reporting mode (concise, summary,
        projects or info)
      fields (list): optional set of fields to report
        (only for 'projects' reporting mode)
      out_file (str): optional, path to a file to write
        the report to (default is to write to stdout)
    """
    # Turn off saving of parameters
    ap._save_params = False
    ap._save_metadata = False
    # Do the reporting
    if mode is None or mode == ReportingMode.INFO:
        f = report_info
        kws = {}
    elif mode == ReportingMode.CONCISE:
        f = report_concise
        kws = {}
    elif mode == ReportingMode.SUMMARY:
        f = report_summary
        kws = {}
    elif mode == ReportingMode.PROJECTS:
        f = report_projects
        kws = { 'fields': fields }
    else:
        raise Exception("Unknown reporting mode")
    # Generate and write the report
    report = f(ap,**kws)
    if out_file:
        fp,temp_file = tempfile.mkstemp()
        with os.fdopen(fp,'w') as fpp:
            fpp.write("%s\n" % report)
        fileops.copy(temp_file,out_file)
        os.remove(temp_file)
        print("Report written to %s" % out_file)
    else:
        print(report)

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
    try:
        projects = ap.get_analysis_projects()
        report.append("\n%d analysis project%s:" % (len(projects),
                                                    "s" if len(projects) != 0
                                                    else ""))
    except Exception as ex:
        projects = []
        report.append("\nNo analysis projects found")
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
                                          if verify_qc(project)
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
    analysis_dir = analysis.AnalysisDir(ap.analysis_dir)
    if analysis_dir.projects:
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
    else:
        # No projects - try loading data from unaligned dir
        try:
            illumina_data = ap.load_illumina_data()
            for p in illumina_data.projects:
                report.append("'%s' (%s sample%s)" %
                              (p.name,
                               len(p.samples),
                               's' if len(p.samples) != 1 else ''))
            report = ', '.join(report)
            report = "no projects found; contents of '%s' are: %s" % \
                     (ap.params.unaligned_dir,
                      report)
        except IlluminaData.IlluminaDataError as ex:
            report = "no projects found"
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
    # Default items to report
    report_items = ['Run name',
                    'Reference',
                    'Platform',
                    'Directory',
                    'Endedness',
                    'Bcl2fastq',]
    # Gather information
    analysis_dir = analysis.AnalysisDir(ap.analysis_dir)
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
    try:
        cellranger_software = ast.literal_eval(
            ap.metadata.cellranger_software)
        report_items.append('Cellranger')
    except ValueError:
        cellranger_software = None
    if ap.metadata.assay is not None:
        assay = ap.metadata.assay
        report_items.append('Assay')
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
    # General information
    field_width = max([len(i) for i in report_items])
    for item in report_items:
        # Get the value for each item
        if item == 'Run name':
            value = run_name
        elif item == 'Reference':
            value = ap.run_reference_id
        elif item == 'Platform':
            value = platform
        elif item == 'Directory':
            value = ap.params.analysis_dir
        elif item == 'Endedness':
            value = ('Paired end' if analysis_dir.paired_end else
                     'Single end')
        elif item == 'Bcl2fastq':
            value = ('Unknown' if not bcl2fastq_software else
                     "%s %s" % (bcl2fastq_software[1],
                                bcl2fastq_software[2]))
        elif item == 'Cellranger':
            value = ('Unknown' if not cellranger_software else
                     "%s %s" % (cellranger_software[1],
                                cellranger_software[2]))
        elif item == 'Assay':
            value = assay
        else:
            raise Exception("Unknown reporting item '%s'" % item)
        # Append a line reporting the value
        report.append("%s%s: %s" % (item,
                                    ' '*(field_width-len(item)),
                                    value))
    report.append("")
    # Projects
    rows = []
    comments = bcf_utils.OrderedDictionary()
    if analysis_dir.n_projects != 0:
        report.append("%d project%s:" % (analysis_dir.n_projects,
                                         '' if analysis_dir.n_projects == 1
                                         else 's'))
        data_items = ('user',
                      'PI',
                      'library_type',
                      'single_cell_platform',
                      'number_of_cells',
                      'organism')
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
                samples += "/%d cell%s" % (
                    int(project_data['number_of_cells']),
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
    else:
        # No projects - try loading data from unaligned dir
        try:
            illumina_data = ap.load_illumina_data()
            report.append("No projects found; '%s' directory contains "
                          "the following data:\n" %
                          ap.params.unaligned_dir)
            for project in illumina_data.projects:
                rows.append(("- '%s':" % project.name,
                             "%s sample%s" % (len(project.samples),
                                              's' if len(project.samples) != 1
                                              else '')))
            report.append(utils.pretty_print_rows(rows))
        except IlluminaData.IlluminaDataError as ex:
            report.append("No projects found")
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

def report_projects(ap,fields=None):
    """Generate one line reports suitable for pasting into spreadsheet

    Generate one-line report for each each project with tab-separated
    data items, suitable for injection into a spreadsheet.

    By default each line has the following information:

    - Run id e.g. HISEQ_140328
    - Run number
    - Source
    - Empty field (for user supplied date)
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

    Alternatively a custom list of fields can be specified via
    the 'fields' argument.

    Composite fields can be specified by joining two or more
    fields with '+' (e.g. 'project+run_id'); the resulting
    value will be the values of the individual fields joined by
    underscores.

    Arguments:
      ap (AutoProcessor): autoprocessor pointing to the
        analysis directory to be reported on
      fields (list): list or tuple of field names to
        output for each project
        
    Returns:
      String with the report text.
    """
    # Set default fields
    if fields is None:
        fields = ('run_id',
                  'run_number',
                  'source',
                  'null',
                  'user',
                  'PI',
                  'application',
                  'single_cell_platform',
                  'organism',
                  'sequencer_platform',
                  'no_of_samples',
                  'no_of_cells',
                  'paired_end',
                  'sample_names',)
    # Acquire data
    analysis_dir = analysis.AnalysisDir(ap.analysis_dir)
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
        project_line = []
        for field in fields:
            # Deal with composite fields (multiple field names
            # joined with '+')
            value = []
            for subfield in field.split('+'):
                value.append(fetch_value(ap,project,subfield))
            # Append to the output line
            project_line.append('_'.join(value))
        report.append('\t'.join(project_line))
    report = '\n'.join(report)
    return report

def fetch_value(ap,project,field):
    """
    Return the value of the supplied field

    Given a field name, return the value determined from
    the data in the supplied AutoProcessor and
    AnalysisProject instances.

    Arguments:
      ap (AutoProcessor): autoprocessor pointing to the
        analysis directory to be reported on
      project (AnalysisProject): project to report on
      field (str): field name to return value of

    Returns:
      String: value of supplied field.
    """
    # Convenience variable for project info
    try:
        info = project.info
    except AttributeError:
        info = None
    # Generate value for supplied field name
    if field == 'datestamp':
        return IlluminaData.split_run_name(ap.run_name)[0]
    elif field == 'run_id':
        return ap.run_reference_id
    elif field == 'run_number':
        return ('' if not ap.metadata.run_number
                else str(ap.metadata.run_number))
    elif field == 'source' or field == 'data_source':
        return ('' if not ap.metadata.source
                else ap.metadata.source)
    elif field == 'analysis_dir' or field == 'path':
        return ap.params.analysis_dir
    elif field == 'project' or field == 'project_name':
        return project.name
    elif field == 'user':
        return ('' if not info.user else info.user)
    elif field == 'PI' or field == 'pi':
        return ('' if not info.PI else info.PI)
    elif field == 'application' or field == 'library_type':
        return ('' if not info.library_type
                else info.library_type)
    elif field == 'single_cell_platform':
        return ('' if not info.single_cell_platform
                else info.single_cell_platform)
    elif field == 'organism':
        return ('' if not info.organism else info.organism)
    elif field == 'sequencer_platform' or field == 'platform':
        return ('' if not ap.metadata.platform
                else str(ap.metadata.platform).upper())
    elif field == 'no_of_samples' or field == '#samples':
        return str(len(project.samples))
    elif field == 'no_of_cells' or field == '#cells':
        return ('' if not info.number_of_cells
                else str(info.number_of_cells))
    elif field == 'paired_end':
        return ('yes' if ap.paired_end else 'no')
    elif field == 'sample_names' or field == 'samples':
        return project.prettyPrintSamples()
    elif field == 'null' or field == '':
        return ''
    else:
        raise KeyError("'%s': unrecognised field for reporting"
                       % field)

def default_value(s,default=""):
    """
    Returns supplied value, or default if supplied value is None
    """
    if s is None:
        return default
    return s
