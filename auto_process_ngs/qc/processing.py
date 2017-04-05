#!/usr/bin/env python
#
# Library for QC reporting of processing stages

import os
from bcftbx.TabFile import TabFile
from ..docwriter import Document
from ..docwriter import Table
from ..docwriter import List
from ..docwriter import Link
from ..docwriter import Img
from .. import css_rules
from .illumina_qc import pretty_print_reads
from .plots import ustackedbar

#######################################################################
# Functions
#######################################################################

def report_processing_qc(analysis_dir,html_file):
    """
    Generate HTML report for processing statistics

    Arguments:
      analysis_dir (AnalysisDir): 
      html_file (str): destination path and file name for
        HTML report
    """
    # Initialise the HTML report
    processing_qc = Document("Processing report for %s" %
                             os.path.basename(analysis_dir.analysis_dir))
    processing_qc.add_css_rule(css_rules.QC_REPORT_CSS_RULES)
    processing_qc.add_css_rule("table { font-size: 80%;\n"
                               "        font-family: sans-serif; }")
    processing_qc.add_css_rule("td { text-align: right; }")
    # Add table of contents
    toc = processing_qc.add_section("Contents",name="toc")
    toc_list = List()
    toc.add(toc_list)
    # Per-lane statistics
    per_lane_stats_file = analysis_dir.params.per_lane_stats_file
    if per_lane_stats_file is None:
        per_lane_stats_file = "per_lane_statistics.info"
    if os.path.exists(per_lane_stats_file):
        per_lane_stats = processing_qc.add_section(
            "Per-lane statistics",
            name="per_lane_stats")
        stats = TabFile(per_lane_stats_file,
                        first_line_is_header=True)
        tbl = Table(columns=stats.header())
        tbl.append_columns("Assigned/unassigned")
        for line in stats:
            n = tbl.add_row()
            for c in stats.header():
                if c in ("Total reads","Assigned reads","Unassigned reads"):
                    value = pretty_print_reads(line[c])
                else:
                    value = line[c]
                tbl.set_value(n,c,value)
            tbl.set_value(n,"Assigned/unassigned",
                          Img(ustackedbar((line["Assigned reads"],
                                           line["Unassigned reads"]),
                                          length=100,height=15,
                                          colors=('red','white'),
                                              inline=True)))
        per_lane_stats.add(tbl)
        toc_list.add_item(Link("Per-lane statistics",per_lane_stats))
    # Per lane by sample statistics
    per_lane_sample_stats_file = "per_lane_sample_stats.info"
    if os.path.exists(per_lane_sample_stats_file):
        per_lane_sample_stats = processing_qc.add_section(
            "Per-lane statistics by sample",
            name="per_lane_sample_stats")
        lane_toc_list = List()
        per_lane_sample_stats.add(lane_toc_list)
        with open("per_lane_sample_stats.info") as stats:
            for line in stats:
                if line.startswith("Lane "):
                    lane = int(line.split(' ')[-1])
                    s = per_lane_sample_stats.add_subsection(
                        "Lane %d" % lane,
                        name="per_lane_sample_stats_lane%d" % lane
                    )
                    preamble = s.add_subsection()
                    lane_toc_list.add_item(Link("Lane %d" % lane,s))
                    current_project = None
                    tbl = Table(columns=('pname','sname',
                                         'nreads','percreads',
                                         'barplot'),
                                pname='Project',
                                sname='Sample',
                                nreads='Nreads',
                                percreads='%reads',
                                barplot='',
                    )
                    s.add(tbl)
                elif line.startswith("Total reads = "):
                    preamble.add(line)
                    total_reads = int(line.split('=')[-1].strip())
                elif line.startswith("- "):
                    pname = line.split()[1].split('/')[0]
                    if pname == current_project:
                        pname = "&nbsp;"
                    else:
                        current_project = pname
                    sname = line.split()[1].split('/')[1]
                    nreads = int(line.split()[2])
                    percreads = line.split()[3]
                    barplot = ustackedbar((nreads,total_reads-nreads),
                                          length=100,height=5,
                                          colors=('black','lightgrey'),
                                          bbox=False,
                                          inline=True)
                    tbl.add_row(pname=pname,
                                sname=sname,
                                nreads=pretty_print_reads(nreads),
                                percreads=percreads,
                                barplot=Img(barplot))
        toc_list.add_item(Link("Per-lane statistics by sample",
                               per_lane_sample_stats),
                          lane_toc_list)
    # Per fastq statistics
    stats_file = "statistics_full.info"
    if not os.path.exists(stats_file):
        if analysis_dir.params.stats_file is not None:
            stats_file = analysis_dir.params.stats_file
        else:
            stats_file = "statistics.info"
    if os.path.exists(stats_file):
        per_file_stats = processing_qc.add_section(
            "Per-file statistics by project",
            name="per_file_stats")
        project_toc_list = List()
        per_file_stats.add(project_toc_list)
        stats = TabFile(stats_file,first_line_is_header=True)
        projects = sorted(list(set([d['Project'] for d in stats])))
        lanes = filter(lambda c: c.startswith('L'),stats.header())
        sample = None
        for project in projects:
            subset = filter(lambda d: d['Project'] == project,stats)
            subset_lanes = filter(lambda l:
                                  reduce(lambda x,y: x and bool(y),
                                         [d[l] for d in subset]),
                                  lanes)
            s = per_file_stats.add_subsection(
                "%s" % project,
                name="per_file_stats_%s" % project
            )
            project_toc_list.add_item(Link("%s" % project,s))
            tbl = Table(columns=('Sample','Fastq','Size'))
            if subset_lanes:
                tbl.append_columns(*subset_lanes)
            tbl.append_columns('Nreads')
            s.add(tbl)
            for line in subset:
                if sample == line['Sample']:
                    sname = "&nbsp;"
                else:
                    sample = line['Sample']
                    sname = sample
                data = { 'Sample': sname,
                         'Fastq': line['Fastq'],
                         'Size': line['Size'],
                         'Nreads': pretty_print_reads(line['Nreads'])}
                for l in subset_lanes:
                    data[l] = pretty_print_reads(line[l])
                tbl.add_row(**data)
        toc_list.add_item(Link("Per-file statistics by project",
                               per_file_stats),
                          project_toc_list)
    # Write the processing QC summary file
    processing_qc.write(html_file)
