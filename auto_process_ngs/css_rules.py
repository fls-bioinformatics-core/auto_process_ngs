#!/usr/bin/env python
#
# CSS rules for styling HTML documents

#######################################################################
# Constants
#######################################################################

QC_REPORT_CSS_RULES = \
"""
h1    { background-color: #42AEC2;
        color: white;
        padding: 5px 10px; }
h2    { background-color: #8CC63F;
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
table { margin: 10 10;
        border: solid 1px grey;
        background-color: white; }
th    { background-color: grey;
        color: white;
        padding: 2px 5px; }
td    { text-align: left;
        vertical-align: top;
        padding: 2px 5px;
        border-bottom: solid 1px lightgray; }
td.param { background-color: grey;
           color: white;
           padding: 2px 5px;
           font-weight: bold; }
p.footer { font-style: italic;
           font-size: 70%; }
"""
