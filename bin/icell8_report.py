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
#from auto_process_ngs.docwriter import Target
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
    args = p.parse_args()

    # Initialise some parameters
    low_read_threshold = 10000
    poly_g_threshold = 5.0

    # Load data from input file
    df = pd.read_csv(args.stats_file,sep='\t')

    # Rename the '%reads_poly_g' column
    df.rename(columns={'%reads_poly_g':'percent_poly_g'},
              inplace=True)

    # Gather the data
    data = AttributeDictionary()
    
    # Total assigned reads
    data['total_assigned_reads'] = df['Nreads'].sum()
    # Mean and median reads per barcode
    data['median_read_count'] = df['Nreads'].median()
    data['mean_read_count'] = df['Nreads'].mean()
    # Number of barcodes (total and assigned)
    data['total_barcodes'] = len(df)
    data['assigned_barcodes'] = len(df[df['Nreads'] > 0])

    # Report data
    print "Total assigned #reads: %d" % data.total_assigned_reads
    print "Median read count    : %d" % data.median_read_count
    print "Mean   read count    : %d" % data.mean_read_count
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
    data_items = (('Total assigned #reads','total_assigned_reads'),
                  ('Total #barcodes','total_barcodes'),
                  ('Assigned #barcodes','assigned_barcodes'),
                  ('Mean read count per barcode','mean_read_count'),
                  ('Median read count per barcode','median_read_count'))
    for item in data_items:
        name,key = item
        tbl.add_row(name=name,value=data[key])
    general_info.add(tbl)

    # Reads at each stage
    reads_per_stage = df[['Nreads',
                          'Nreads_filtered',
                          'Nreads_trimmed',
                          'Nreads_contaminant_filtered']].sum()
    fig=plt.figure()
    plot = reads_per_stage.plot.bar()
    plot.set_xticklabels(['Assigned','Quality filtered','Trimmed','Uncontaminated'],
                         rotation=45)
    plot.set_ylabel("#read pairs")
    plot.get_figure().savefig("reads_per_stage.png",
                              bbox_inches='tight')
    general_info.add(Img("reads_per_stage.png"))
    general_info.add(reads_per_stage.to_frame().to_html())
    toc_list.add_item(Link(general_info.title,general_info))
    
    # Low read counts
    read_counts = report.add_section("Read counts",
                                     name="read_counts")
    #low_read_count = df.query("Nreads < %d" % low_read_threshold)
    low_read_count = df[['#Barcode','Nreads']].query("Nreads < %d" % low_read_threshold)
    n_low_reads = len(low_read_count)
    read_counts.add("%d barcodes with less than %d reads" %
                    (n_low_reads,
                     low_read_threshold))
    if n_low_reads:
        read_counts.add(low_read_count.to_html())
    print "#barcodes < %d reads: %d" % (low_read_threshold,
                                        len(low_read_count))
    toc_list.add_item(Link(read_counts.title,read_counts))

    # Histogram of initial and final read counts
    fig = plt.figure()
    ax = fig.add_subplot(1,2,1)
    df.query('Nreads > 0')[['Nreads']].plot.hist(ax=ax,
                                                 by='Nreads_final',
                                                 bins=100,
                                                 legend=False)
    ax.set_title("Initial")
    ax.set_xlabel("No of reads/barcode")
    ax.set_ylabel("No of barcodes")
    ax = fig.add_subplot(1,2,2)
    df.query('Nreads > 0')[['Nreads_final']].plot.hist(ax=ax,
                                                       by='Nreads_final',
                                                       bins=100,
                                                       legend=False,)
    ax.set_title("Final")
    ax.set_xlabel("No of reads/barcode")
    ax.set_ylabel(None)
    ax.get_figure().savefig("read_dist.png")
    read_counts.add(Img("read_dist.png"))

    # Sample info
    sample_info = report.add_section("Samples",
                                     name="sample_info")
    samples = df[['Sample',
                  'Nreads',
                  'Nreads_filtered',
                  'Nreads_trimmed',
                  'Nreads_contaminant_filtered']].groupby('Sample').sum()
    sample_info.add(samples.to_html())
    toc_list.add_item(Link(sample_info.title,sample_info))

    # Group by sample
    fig=plt.figure()
    plot = samples.transpose().plot.bar(stacked=True)
    plot.set_xticklabels(['Assigned',
                          'Quality filtered',
                          'Trimmed',
                          'Uncontaminated'],
                         rotation=45)
    plot.set_ylabel("#read pairs")
    plot.get_figure().savefig("samples.png",bbox_inches='tight')
    sample_info.add(Img("samples.png"))
    
    # Poly-G regions
    high_poly_g = df[['#Barcode',
                      'Nreads_poly_g',
                      'percent_poly_g']].query("percent_poly_g > %f" %
                                               poly_g_threshold)
    median_poly_g = df['percent_poly_g'].median()
    mean_poly_g = df['percent_poly_g'].mean()
    poly_g_info = report.add_section("Poly-G regions",
                                     name="poly_g_info")
    poly_g_info.add("Mean %%reads with poly-G regions per barcode: %.2f%%"
                    % mean_poly_g)
    poly_g_info.add("Median %%reads with poly-G regions per barcode: %.2f%%"
                    % median_poly_g)
    poly_g_info.add("%d barcodes with > %.2f%% reads with poly-G regions"
                    % (len(high_poly_g['Nreads_poly_g']),
                       poly_g_threshold))
    if len(high_poly_g['Nreads_poly_g']):
        poly_g_info.add(high_poly_g.to_html())
    print "#reads with > %f%% poly-G: %d" % (poly_g_threshold,
                                             high_poly_g['Nreads_poly_g'].sum())
    print "Mean   poly-G percentage: %f%%" % mean_poly_g
    print "Median poly-G percentage: %f%%" % median_poly_g

    # Histogram of poly-g reads
    fig=plt.figure()
    plot = df[['percent_poly_g']].plot.hist(by='percent_poly_g',
                                            bins=100,
                                            legend=False)
    plot.get_figure().savefig("poly_g_dist.png")
    poly_g_info.add(Img("poly_g_dist.png"))
    toc_list.add_item(Link(poly_g_info.title,poly_g_info))

    # Histogram of low read counts
    ##fig=plt.figure()
    ##plot = low_read_count[['Nreads']].plot.hist(by='Nreads',
    ##                                            bins=20,
    ##                                          legend=False)
    ##plot.get_figure().savefig("hist_low_counts.png")

    # Generate the HTML report
    report.write("icell8_processing.html")
