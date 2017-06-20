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
from auto_process_ngs.qc.illumina_qc import QCReporter
from auto_process_ngs import get_version

"""
qc_reporter2

"""

#######################################################################
# Main program
#######################################################################

def main():
    # Deal with command line
    p = optparse.OptionParser(usage="%prog DIR [DIR...]",
                              version="%prog "+get_version(),
                              description="Generate QC report for each directory "
                              "DIR")
    p.add_option('--qc_dir',action='store',dest='qc_dir',default='qc',
                 help="explicitly specify QC output directory (nb if "
                 "supplied then the same QC_DIR will be used for each "
                 "DIR. Non-absolute paths are assumed to be relative to "
                 "DIR). Default: 'qc'")
    p.add_option('-f','--filename',action='store',dest='filename',
                 default=None,
                 help="file name for output QC report (default: "
                 "<DIR>/<QC_DIR>_report.html)")
    p.add_option('--verify',action='store_true',dest='verify',
                 help="verify the QC products only (don't write the "
                 "report)")
    opts,args = p.parse_args()
    if len(args) < 1:
        p.error("Need to supply at least one directory")

    # Examine projects i.e. supplied directories
    for d in args:
        project_name = os.path.basename(d)
        dir_path = os.path.abspath(d)
        p = AnalysisProject(project_name,dir_path)
        print "Project: %s" % p.name
        if opts.qc_dir is None:
            qc_dir = p.qc_dir
        else:
            qc_dir = opts.qc_dir
        if not os.path.isabs(qc_dir):
            qc_dir = os.path.join(p.dirn,qc_dir)
        print "QC output dir: %s" % qc_dir
        print "-"*(len('Project: ')+len(p.name))
        print "%d samples | %d fastqs" % (len(p.samples),len(p.fastqs))
        if opts.verify:
            if not QCReporter(p).verify(qc_dir=qc_dir):
                print "Verification: FAILED"
            else:
                print "Verification: OK"
        else:
            qc_base = os.path.basename(qc_dir)
            if opts.filename is None:
                out_file = '%s_report.html' % qc_base
            else:
                out_file = opts.filename
            if not os.path.isabs(out_file):
                out_file = os.path.join(p.dirn,out_file)
            print "Writing QC report to %s" % out_file
            qc = QCReporter(p).report(qc_dir=qc_dir,
                                      filename=out_file)

if __name__ == '__main__':
    main()
    
            
        
