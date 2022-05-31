#!/usr/bin/env python
#
#     qc.constants.py: constants for QC pipeline
#     Copyright (C) University of Manchester 2019-2022 Peter Briggs
#

"""
Provides the following constants for the QC pipeline:

- FASTQ_SCREENS: tuple of screen names
- PROTOCOLS: tuple of QC protocol names
- SEQUENCE_DEDUP_CUTOFFS: dictionary of sequence deduplication thresholds
- QC_REPORT_CSS_STYLES: text of CSS styles for HTML QC reports
"""
######################################################################
# QC constants
######################################################################

# QC protocols
from .protocols import QC_PROTOCOLS
PROTOCOLS = tuple(QC_PROTOCOLS.keys())

# Screen names
FASTQ_SCREENS = ('model_organisms',
                 'other_organisms',
                 'rRNA',)

# Sequence deduplication thresholds
SEQUENCE_DEDUP_CUTOFFS = {
    "warn": 0.3,
    "fail": 0.2,
}

# CSS for QC report
QC_REPORT_CSS_STYLES = """/* General styles */
html { font-family: sans-serif; }
p { font-size: 85%;
    color: #808080; }
/* Headers */
h1, h2, h3, h4 { font-family: DejaVu Serif, serif; }
h1 { background-color: #42AEC2;
     color: white;\n
     padding: 5px 10px; }
h2 { background-color: #8CC63F;
     color: white;
     display: inline-block;
     padding: 5px 15px;
     margin: 0;
     border-top-left-radius: 20px;
     border-bottom-right-radius: 20px; }
h3, h4 { background-color: grey;
         color: white;
         display: block;
         padding: 5px 15px;
         margin: 0;
         border-top-left-radius: 20px;
         border-bottom-right-radius: 20px; }
.single_library_summary h3, .single_library_summary h4 {
     background-color: white;
     color: #808080; }
/* Summary section */
div.summary { margin: 10 10;
              border: solid 2px #8CC63F;
              padding: 0;
              border-top-left-radius: 25px; }
.info { padding: 5px 15px;
        float: left;
        font-size: 100%; }
.info h3 { margin: 5px; }
.info p { color: black; }
/* Samples and Fastqs */
.sample { margin: 10 10;
          border: solid 2px #8CC63F;
          padding: 0;
          background-color: #ffe;
          border-top-left-radius: 25px;
          border-bottom-right-radius: 25px; }
.fastqs { border: 1px solid grey;
          padding: 5px;
          margin: 5px 20px; }
.fastq { border: 2px solid lightgray;
         padding: 5px;
         margin: 5px;
         float: left; }
.strandedness { border: 2px solid lightgray;
                padding: 5px;
                margin: 5px;float: left; }
/* Metadata table */
table.metadata {
          margin: 10 10;
          border: solid 1px grey;
          background-color: white;
          font-size: 90%; }
table.metadata tr td:first-child {
          background-color: grey;
          color: white;
          padding: 2px 5px;
          font-weight: bold;
          vertical-align: top; }
/* Summary table */
table.summary { border: solid 1px grey;
                background-color: white;
                margin: 10 10;
                font-size: 80%; }
table.summary th { background-color: grey;
                   color: white;
                   padding: 2px 5px; }
table.summary td { text-align: center;
                   padding: 2px 5px;
                   border-bottom: solid 1px lightgray; }
table.summary td.single_library_analyses { text-align: left; }
table.summary tr td:first-child { text-align: right; }
table.summary tr td:first-child {
          background-color: grey;
          color: white;
          font-weight: bold; }
table.summary tr td:first-child a {
          color: white;
          font-weight: bold; }
/* Warnings section */
.warnings { padding: 2px;
            border: solid 3px red;
            color: red;
            background-color: #F5BCA9;
            font-weight: bold;
            margin: 10px; }
.warnings p   { color: red;
                font-size: 120%; }
.warnings img { vertical-align: middle; }
/* Toggle sections */
.toggle_button {
  background-color: #4CAF50; /* Green */
  border: none;
  color: white;
  padding: 15px 32px;
  text-align: center;
  text-decoration: none;
  display: block;
  font-size: 16px;
  cursor: pointer;
  margin: 5px 0px 0px 10px;
}
.toggle_section {
  border: solid 2px #4CAF50;
  margin: 0px 5px 10px 10px;
  padding: 5px;
}
/* Display control elements */
.clear { clear: both; }
.hide  { display: none; }
/* FastQC summary table */
table.fastqc_summary span.PASS { font-weight: bold;
                                 color: green; }
table.fastqc_summary span.WARN { font-weight: bold;
                                 color: orange; }
table.fastqc_summary span.FAIL { font-weight: bold;
                                 color: red; }
/* Program versions */
table.programs th { text-align: left;
                    background-color: grey;
                    color: white;
                    padding: 2px 5px; }
table.programs td { padding: 2px 5px;
                    border-bottom: solid 1px lightgray; }
/* Rules for printing */
@media print
{
a { color: black; text-decoration: none; }
.sample { page-break-before: always; }
table th { border-bottom: solid 1px lightgray; }
.no_print { display: none; }
}
"""

# Javascript for QC report

# Function to toggle a report section
# Based on code at
# https://www.w3schools.com/howto/howto_js_toggle_hide_show.asp
# https://www.w3schools.com/howto/howto_js_toggle_text.asp
JS_TOGGLE_FUNCTION = """
// Toggle the display style of an element
function toggleBlock(id,button_id,show_text,hide_text) {
  var x = document.getElementById(id);
  var b = document.getElementById(button_id);
  if (x.style.display === "none") {
    x.style.display = "block";
    b.innerHTML = hide_text;
  } else {
    x.style.display = "none";
    b.innerHTML = show_text;
  }
}
"""
