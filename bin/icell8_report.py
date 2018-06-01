#!/usr/bin/env python
#
#     icell8_report.py: perform processing of Wafergen iCell8 data
#     Copyright (C) University of Manchester 2017 Peter Briggs
#
"""
icell8_report.py

Utility to report summary of ICell8 processing pipeline.
"""

######################################################################
# Imports
######################################################################

import os
import sys
import argparse
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from bcftbx.utils import AttributeDictionary
from auto_process_ngs.analysis import AnalysisProject
from auto_process_ngs.docwriter import Document
from auto_process_ngs.docwriter import Table
from auto_process_ngs.docwriter import List
from auto_process_ngs.docwriter import Link
from auto_process_ngs.docwriter import Img
from auto_process_ngs.utils import ZipArchive
import auto_process_ngs.css_rules as css_rules

# Initialise logging
import logging
logger = logging.getLogger(__name__)

######################################################################
# Main
######################################################################

if __name__ == "__main__":
    # Set plotting style
    matplotlib.style.use('ggplot')

    # Command line
    p = argparse.ArgumentParser()
    p.add_argument("project_dir",metavar="DIR",
                   nargs="?",default=None,
                   help="directory with ICell8 processing outputs")
    p.add_argument("-s","--stats_file",default=None,
                   help="ICell8 stats file (default: "
                   "DIR/stats/icell8_stats.tsv)")
    p.add_argument("-o","--out_file",action='store',
                   default="icell8_processing.html",
                   help="Output HTML file (default: "
                   "'icell8_processing.html')")
    p.add_argument("-n","--name",action="store",default=None,
                   help="specify a string to append to the zip "
                   "archive name and prefix")
    args = p.parse_args()

    # Project directory
    project_dir = args.project_dir
    if project_dir is None:
        project_dir = os.getcwd()
    else:
        project_dir = os.path.abspath(project_dir)
    print "Project dir: %s" % project_dir

    # Stats file
    if args.stats_file is None:
        stats_file = os.path.join(project_dir,
                                  "stats",
                                  "icell8_stats.tsv")
    else:
        stats_file = os.path.abspath(args.stats_file)
    print "Stats file: %s" % stats_file

    # Check for XLSX version of stats
    xlsx_file = os.path.splitext(stats_file)[0]+".xlsx"
    if not os.path.exists(xlsx_file):
        xlsx_file = None
    print "XLSX file: %s" % xlsx_file

    # Output file and directory for output data (images etc)
    out_file = args.out_file
    if not os.path.isabs(out_file):
        out_file = os.path.join(project_dir,out_file)
    out_dir = os.path.splitext(out_file)[0]+"_data"
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    print "Output file: %s" % out_file
    print "Data dir   : %s" % out_dir

    # Sort out title/name
    report_name = args.name
    report_title = ""
    if report_name is None:
        # Assume this is an analysis project
        project = AnalysisProject(os.path.basename(project_dir),
                                  project_dir)
        if project.is_analysis_dir:
            report_name = ".%s.%s" % (project.name,
                                      os.path.basename(
                                          os.path.dirname(
                                              project_dir)))

    # Load data from input file
    df = pd.read_csv(stats_file,sep='\t')

    # Check columns
    cols = list(df.columns.values)
    for c in ('Nreads',
              'Nreads_filtered',
              'Nreads_poly_g',
              'Nreads_trimmed',
              'Nreads_final'):
        if c not in cols:
            logging.critical("Missing column '%s' for stats file"
                             % c)
            sys.exit(1)
    # Special case: no filtered reads
    if 'Nreads_contaminant_filtered' in cols:
        contaminant_filtered = True
    else:
        logging.warning("No stats on contaminant filtering")
        contaminant_filtered = False

    # Rename the '#Barcodes' and '%reads_poly_g' columns
    df.rename(columns={'#Barcode':'Barcode',
                       '%reads_poly_g':'percent_poly_g'},
              inplace=True)
    print df.head()

    # Gather the data
    data = AttributeDictionary()

    # Total reads
    data['total_reads'] = df['Nreads'].sum()
    
    # Total assigned reads
    df = df.drop(df[df['Barcode'] == 'Unassigned'].index)
    data['total_assigned_reads'] = df['Nreads'].sum()
    # Mean and median reads per barcode
    data['median_read_count'] = df['Nreads'].median()
    data['mean_read_count'] = df['Nreads'].mean()
    data['std_read_count'] = df['Nreads'].std()
    # Number of barcodes (total and assigned)
    data['total_barcodes'] = len(df)
    data['assigned_barcodes'] = len(df[df['Nreads'] > 0])

    # Report data
    print "Total #reads         : %d" % data.total_reads
    print "Total assigned #reads: %d" % data.total_assigned_reads
    print "Median read count    : %d" % data.median_read_count
    print "Mean   read count    : %d" % data.mean_read_count
    print "Std    read count    : %d" % data.std_read_count
    print "Total #barcodes      : %d" % data.total_barcodes
    print "Assigned #barcodes   : %d" % data.assigned_barcodes

    # Make a HTML report
    report = Document(title="ICell8 processing summary")
    report.add_css_rule(css_rules.QC_REPORT_CSS_RULES)
    report.add_css_rule("body { border: solid 1px grey;\n"
                        "       font-size: 80%;\n"
                        "       font-family: sans-serif; }")

    # Table of contents
    toc = report.add_section("Contents",name="toc")
    toc_list = List()
    toc.add(toc_list)

    # Files
    files_info = report.add_section("Stats files")
    stats_files = "Statistics: %s" % Link("[TSV]",
                                          os.path.relpath(
                                              stats_file,
                                              os.path.dirname(out_file)))
    if xlsx_file is not None:
        stats_files += " | %s" % Link("[XLSX]",
                                      os.path.relpath(
                                          xlsx_file,
                                      os.path.dirname(out_file)))
    files_info.add(stats_files)
    toc_list.add_item(Link(files_info.title,files_info))

    # General info
    general_info = report.add_section("General info",
                                      name="general_info")
    tbl = Table(columns=('name','value'))
    tbl.no_header()
    data_items = (('Total #reads','total_reads'),
                  ('Total assigned #reads','total_assigned_reads'),
                  ('Total #barcodes','total_barcodes'),
                  ('Assigned #barcodes','assigned_barcodes'))
    for item in data_items:
        name,key = item
        tbl.add_row(name=name,value=data[key])
    general_info.add(tbl)

    # Reads at each stage
    cols = ['Nreads',
            'Nreads_filtered',
            'Nreads_trimmed']
    if contaminant_filtered:
        cols.append('Nreads_contaminant_filtered')
    plot_filen = os.path.join(out_dir,"reads_per_stage.png")
    reads_per_stage = \
        pd.Series([data.total_reads,],
                  ['Nreads_initial',]).append(
                      df[cols].sum())
    print reads_per_stage
    fig=plt.figure()
    plot = reads_per_stage.plot.bar(figsize=(6,4))
    labels = ['Initial',
              'Assigned',
              'Quality filtered',
              'Trimmed']
    if contaminant_filtered:
        labels.append('Uncontaminated')
    plot.set_xticklabels(labels,rotation=45)
    plot.set_ylabel("#read pairs")
    plot.get_figure().savefig(plot_filen,bbox_inches='tight')
    general_info.add(Img(os.path.relpath(plot_filen,
                                         os.path.dirname(out_file))))
    general_info.add(reads_per_stage.to_frame().to_html(header=False))
    toc_list.add_item(Link(general_info.title,general_info))
    
    # Low read counts
    low_read_threshold = data.total_assigned_reads/data.total_barcodes/10
    read_counts = report.add_section("Read counts",
                                     name="read_counts")
    low_read_count = df[['Barcode','Nreads']].query(
        "Nreads < %d" % low_read_threshold).query(
            "Nreads > 0")
    low_reads = len(low_read_count)
    no_reads = len(df.query("Nreads == 0"))
    tbl = Table(columns=('name','value'))
    tbl.no_header()
    data_items = (('Mean read count per barcode','mean_read_count'),
                  ('Median read count per barcode','median_read_count'))
    for item in data_items:
        name,key = item
        tbl.add_row(name=name,value=data[key])
    tbl.add_row(name="# barcodes with no reads",
                value="%d (%0.2f%%)" %
                (no_reads,
                 float(no_reads)/float(data.total_barcodes)*100.0))
    tbl.add_row(name="# barcodes with less than %d reads" %
                low_read_threshold,
                value="%d (%0.2f%%)" %
                (low_reads,
                 float(low_reads)/float(data.total_barcodes)*100.0))
    read_counts.add(tbl)

    # Histogram of initial and final read counts
    plot_filen = os.path.join(out_dir,"read_dist.png")
    fig = plt.figure(figsize=(12,4))
    ax = fig.add_subplot(1,2,1)
    df.query('Nreads > 0')[['Nreads']].plot.hist(ax=ax,
                                                 by='Nreads_final',
                                                 bins=100,
                                                 legend=False,
                                                 edgecolor='black')
    ax.set_title("Initial distribution")
    ax.set_xlabel("No of reads/barcode")
    ax.set_ylabel("No of barcodes")
    ax = fig.add_subplot(1,2,2)
    df.query('Nreads > 0')[['Nreads_final']].plot.hist(ax=ax,
                                                       by='Nreads_final',
                                                       bins=100,
                                                       legend=False,
                                                       edgecolor='black')
    ax.set_title("Final distribution")
    ax.set_xlabel("No of reads/barcode")
    ax.yaxis.label.set_visible(False)
    ax.get_figure().savefig(plot_filen,bbox_inches='tight')
    read_counts.add(Img(os.path.relpath(plot_filen,
                                        os.path.dirname(out_file))))

    # List barcodes with low read counts
    if low_reads:
        read_counts.add(low_read_count.to_html(index=False))
    print "#barcodes < %d reads: %d" % (low_read_threshold,
                                        len(low_read_count))
    toc_list.add_item(Link(read_counts.title,read_counts))

    # Sample info
    sample_info = report.add_section("Samples",
                                     name="sample_info")
    cols = ['Sample',
            'Nreads',
            'Nreads_filtered',
            'Nreads_trimmed',]
    if contaminant_filtered:
        cols.append('Nreads_contaminant_filtered')
    samples = df[cols].groupby('Sample').sum()
    sample_info.add(samples.to_html(na_rep='-'))
    toc_list.add_item(Link(sample_info.title,sample_info))

    # Group by sample
    plot_filen = os.path.join(out_dir,"samples.png")
    fig=plt.figure()
    plot = samples.transpose().plot.bar(stacked=True,
                                        figsize=(6,4))
    labels = ['Assigned',
              'Quality filtered',
              'Trimmed']
    if contaminant_filtered:
        labels.append('Uncontaminated')
    plot.set_xticklabels(labels,rotation=45)
    plot.set_ylabel("#read pairs")
    plot.legend(loc="upper left", bbox_to_anchor=(1,1))
    plot.get_figure().savefig(plot_filen,bbox_inches='tight')
    sample_info.add(Img(os.path.relpath(plot_filen,
                                        os.path.dirname(out_file))))
    
    # Poly-G regions
    poly_g_threshold = 20.0
    high_poly_g = df[['Barcode',
                      'Nreads_poly_g',
                      'percent_poly_g']].query("percent_poly_g > %f" %
                                               poly_g_threshold)
    median_poly_g = df['percent_poly_g'].median()
    mean_poly_g = df['percent_poly_g'].mean()
    total_poly_g = df['Nreads_poly_g'].sum()
    poly_g_info = report.add_section("Poly-G regions",
                                     name="poly_g_info")
    tbl = Table(columns=('name','value'))
    tbl.no_header()
    tbl.add_row(name="# assigned reads with poly-G regions",
                value="%d (%.2f%%)" %
                (total_poly_g,
                 float(total_poly_g)/float(data.total_assigned_reads)*100.0))
    tbl.add_row(name="Mean %reads with poly-G regions per barcode",
                value="%.2f%%" % mean_poly_g)
    tbl.add_row(name="Median %reads with poly-G regions per barcode",
                value="%.2f%%" % median_poly_g)
    tbl.add_row(name="# barcodes with > %.2f%% reads with poly-G regions" %
                poly_g_threshold,
                value="%d" % len(high_poly_g['Nreads_poly_g']))
    poly_g_info.add(tbl)
    print "#reads with poly-G regions: %d" % high_poly_g['Nreads_poly_g'].sum()
    print "Mean   poly-G percentage: %f%%" % mean_poly_g
    print "Median poly-G percentage: %f%%" % median_poly_g

    # Histogram of poly-g reads
    plot_filen = os.path.join(out_dir,"poly_g_dist.png")
    fig=plt.figure()
    plot = df[['percent_poly_g']].query("percent_poly_g > 0").plot.hist(
        by='percent_poly_g',
        bins=100,
        legend=False,
        edgecolor='black',
        figsize=(6,4))
    plot.set_title("Poly-G distribution")
    plot.set_xlabel("%reads with poly-G regions")
    plot.set_ylabel("No of barcodes")
    plot.get_figure().savefig(plot_filen,bbox_inches='tight')
    poly_g_info.add(Img(os.path.relpath(plot_filen,
                                        os.path.dirname(out_file))))
    # List of barcodes with high percentage of reads poly-G content
    if len(high_poly_g['Nreads_poly_g']):
        poly_g_info.add(high_poly_g.to_html(index=False))
    toc_list.add_item(Link(poly_g_info.title,poly_g_info))

    # Generate the HTML report
    report.write(out_file)

    # Collect everything into a zip archive
    parent_dir = os.path.dirname(out_file)
    zip_name = os.path.splitext(os.path.basename(out_file))[0]
    if report_name is not None:
        zip_name = "%s%s" % (zip_name,report_name)
    report_zip = os.path.join(parent_dir,"%s.zip" % zip_name)
    zip_file = ZipArchive(report_zip,
                          relpath=parent_dir,
                          prefix="%s" % zip_name)
    # Add the HTML report
    zip_file.add_file(out_file)
    # Add the data directory
    zip_file.add_dir(out_dir)
    # Add the stats file
    zip_file.add_file(stats_file)
    if xlsx_file is not None:
        zip_file.add_file(xlsx_file)
    zip_file.close()
    print "Wrote zip archive: %s" % report_zip
