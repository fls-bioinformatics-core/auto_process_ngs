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
import argparse
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from bcftbx.utils import AttributeDictionary
from auto_process_ngs.docwriter import Document
from auto_process_ngs.docwriter import Table
from auto_process_ngs.docwriter import List
from auto_process_ngs.docwriter import Link
from auto_process_ngs.docwriter import Img
from auto_process_ngs.utils import ZipArchive
import auto_process_ngs.css_rules as css_rules

######################################################################
# Main
######################################################################

if __name__ == "__main__":
    # Set plotting style
    matplotlib.style.use('ggplot')

    # Command line
    p = argparse.ArgumentParser()
    p.add_argument("stats_file",help="ICell8 stats file")
    p.add_argument("out_file",
                   nargs="?",default="icell8_processing.html",
                   help="Output HTML file (default: "
                   "'icell8_processing.html')")
    args = p.parse_args()

    # Output file name
    out_file = os.path.abspath(args.out_file)

    # Directory for output data (images etc)
    out_dir = os.path.splitext(out_file)[0]+"_data"
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    # Load data from input file
    df = pd.read_csv(args.stats_file,sep='\t')

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
    plot_filen = os.path.join(out_dir,"reads_per_stage.png")
    reads_per_stage = \
        pd.Series([data.total_reads,],
                  ['Nreads_initial',]).append(
                      df[['Nreads',
                          'Nreads_filtered',
                          'Nreads_trimmed',
                          'Nreads_contaminant_filtered']].sum())
    print reads_per_stage
    fig=plt.figure()
    plot = reads_per_stage.plot.bar()
    plot.set_xticklabels(['Initial',
                          'Assigned',
                          'Quality filtered',
                          'Trimmed',
                          'Uncontaminated'],
                         rotation=45)
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
    samples = df[['Sample',
                  'Nreads',
                  'Nreads_filtered',
                  'Nreads_trimmed',
                  'Nreads_contaminant_filtered']].groupby('Sample').sum()
    sample_info.add(samples.to_html(na_rep='-'))
    toc_list.add_item(Link(sample_info.title,sample_info))

    # Group by sample
    plot_filen = os.path.join(out_dir,"samples.png")
    fig=plt.figure()
    plot = samples.transpose().plot.bar(stacked=True)
    plot.set_xticklabels(['Assigned',
                          'Quality filtered',
                          'Trimmed',
                          'Uncontaminated'],
                         rotation=45)
    plot.set_ylabel("#read pairs")
    plot.legend(loc='best')
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
        edgecolor='black')
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
    report_zip = os.path.join(parent_dir,"%s.zip" % zip_name)
    zip_file = ZipArchive(report_zip,
                          relpath=parent_dir,
                          prefix="%s" % zip_name)
    # Add the HTML report
    zip_file.add_file(out_file)
    # Add the data directory
    zip_file.add_dir(out_dir)
    zip_file.close()
    print "Wrote zip archive: %s" % report_zip
