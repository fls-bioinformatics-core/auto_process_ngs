#!/usr/bin/env python
#
#     qcreporter2.py: generate report file for Illumina NGS qc runs
#     Copyright (C) University of Manchester 2015 Peter Briggs
#

#######################################################################
# Imports
#######################################################################

import sys
import os
import optparse
from auto_process_ngs.utils import AnalysisProject
from auto_process_ngs.qcreport.illumina import QCReporter
from auto_process_ngs import get_version

"""
qc_reporter2

"""

#######################################################################
# Main program
#######################################################################

def main():
    # Deal with command line
    p = optparse.OptionParser(usage="%prog DIR [DIR]",
                              version="%prog "+get_version(),
                              description="Generate QC report for each directory "
                              "DIR")
    p.add_option('--verify',action='store_true',dest='verify',
                 help="verify the QC products only (don't write the report)")
    opts,args = p.parse_args()
    if len(args) < 1:
        p.error("Need to supply at least one directory")

    # Examine projects i.e. supplied directories
    for d in args:
        project_name = os.path.basename(d)
        dir_path = os.path.abspath(d)
        p = AnalysisProject(project_name,dir_path)
        print "Project: %s" % p.name
        print "-"*(len('Project: ')+len(p.name))
        print "%d samples | %d fastqs" % (len(p.samples),len(p.fastqs))
        if opts.verify:
            if not QCReporter(p).verify():
                print "Verification: FAILED"
            else:
                print "Verification: OK"
        else:
            qc = QCReporter(p).report()

if __name__ == '__main__':
    main()
    
            
        
