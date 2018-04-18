#!/usr/bin/env python
#
#     reportqc.py: generate report file for Illumina NGS qc runs
#     Copyright (C) University of Manchester 2015-2018 Peter Briggs
#

#######################################################################
# Imports
#######################################################################

import sys
import os
import optparse
import logging
from bcftbx.utils import find_program
from auto_process_ngs.utils import AnalysisProject
from auto_process_ngs.utils import ZipArchive
from auto_process_ngs.applications import Command
from auto_process_ngs.qc.illumina_qc import QCReporter
from auto_process_ngs.qc.illumina_qc import expected_qc_outputs
from auto_process_ngs import get_version

# Module specific logger
logger = logging.getLogger(__name__)

__version__ = get_version()

"""
reportqc

Utility to verify and report on QC outputs from
auto_process pipeline.
"""

#######################################################################
# Functions
#######################################################################

def verify_qc(project,qc_dir=None):
    """
    Get list of fastqs in project failing verification

    Arguments:
      project (AnalysisProject): project object
      qc_dir (str): optional name of subdirectory
        containing QC outputs (defaults to default
        QC subdir from the project)

    Returns:
      List: list of Fastqs (including path) which
        don't pass the verification check.
    """
    if qc_dir is None:
        qc_dir = self.qc_dir
    fastqs = []
    for sample in project.samples:
        for fq in sample.fastq:
            if not sample.verify_qc(qc_dir,fq):
                fastqs.append(fq)
    return fastqs

def zip_report(project,report_html,qc_dir=None):
    """
    Create ZIP archive for a QC report

    Arguments:
      project (AnalysisProject): project object
      report_html (str): HTML QC report
      qc_dir (str): optional name of subdirectory
        containing QC outputs (defaults to default
        QC subdir from the project)

    Returns:
      String: path to the output ZIP file.
    """
    # Name for ZIP file
    basename = os.path.splitext(os.path.basename(report_html))[0]
    analysis_dir = os.path.basename(os.path.dirname(project.dirn))
    # Create ZIP archive
    report_zip = os.path.join(project.dirn,
                              "%s.%s.%s.zip" %
                              (basename,
                               project.name,
                               analysis_dir))
    zip_file = ZipArchive(report_zip,relpath=project.dirn,
                          prefix="%s.%s.%s" %
                          (basename,
                           project.name,
                           analysis_dir))
    # Get QC dir if not set
    if qc_dir is None:
        qc_dir = self.qc_dir
    # Add the HTML report
    zip_file.add_file(report_html)
    # Add the FastQC and screen files
    for sample in project.qc.samples:
        for fastqs in sample.fastq_pairs:
            for fq in fastqs:
                logger.debug("Adding QC outputs for %s" % fq)
                for f in expected_qc_outputs(fq,qc_dir):
                    if f.endswith('.zip'):
                        # Exclude .zip file
                        continue
                    if os.path.exists(f):
                        zip_file.add(f)
    # Finished
    return report_zip

#######################################################################
# Main program
#######################################################################

def main():
    # Deal with command line
    p = optparse.OptionParser(usage="%prog DIR [DIR...]",
                              version="%prog "+__version__,
                              description="Generate QC report for each directory "
                              "DIR")
    p.add_option('--fastq_dir',
                 action='store',dest='fastq_dir',default=None,
                 help="explicitly specify subdirectory of DIRs with "
                 "Fastq files to run the QC on.")
    p.add_option('--qc_dir',action='store',dest='qc_dir',default='qc',
                 help="explicitly specify QC output directory (nb if "
                 "supplied then the same QC_DIR will be used for each "
                 "DIR. Non-absolute paths are assumed to be relative to "
                 "DIR). Default: 'qc'")
    reporting = optparse.OptionGroup(p,'Reporting options')
    reporting.add_option('-t','--title',action='store',dest='title',
                         default=None,
                         help="title for output QC reports")
    reporting.add_option('-f','--filename',
                         action='store',dest='filename',default=None,
                 help="file name for output HTML QC report (default: "
                         "<DIR>/<QC_DIR>_report.html)")
    reporting.add_option('--zip',action='store_true',
                         dest='zip',default=False,
                         help="make ZIP archive for the QC report")
    reporting.add_option('--multiqc',action='store_true',
                         dest='multiqc',default=False,
                         help="generate MultiQC report")
    reporting.add_option('--force',action='store_true',
                         dest='force',default=False,
                         help="force generation of reports even if "
                         "verification fails")
    p.add_option_group(reporting)
    verification = optparse.OptionGroup(p,'Verification options')
    verification.add_option('--verify',action='store_true',dest='verify',
                            help="verify the QC products only (don't "
                            "write the report)")
    verification.add_option('-l','--list-unverified',action='store_true',
                            dest='list_unverified',default=False,
                            help="list the Fastqs that failed "
                            "verification")
    p.add_option_group(verification)
    opts,args = p.parse_args()
    if len(args) < 1:
        p.error("Need to supply at least one directory")

    # Report name and version
    print "%s version %s" % (os.path.basename(sys.argv[0]),__version__)

    # Check for MultiQC if required
    if opts.multiqc:
        if find_program("multiqc") is None:
            logger.critical("MultiQC report requested but 'multiqc' "
                            "not available")
            sys.exit(1)

    # Examine projects i.e. supplied directories
    retval = 0
    for d in args:
        dir_path = os.path.abspath(d)
        project_name = os.path.basename(dir_path)
        p = AnalysisProject(project_name,dir_path,
                            fastq_dir=opts.fastq_dir)
        print "Project: %s" % p.name
        print "Fastqs : %s" % p.fastq_dir
        if opts.qc_dir is None:
            qc_dir = p.qc_dir
        else:
            qc_dir = opts.qc_dir
        if not os.path.isabs(qc_dir):
            qc_dir = os.path.join(p.dirn,qc_dir)
        # Warning if there is a mismatch
        qc_info = p.qc_info(qc_dir)
        if qc_info.fastq_dir is not None and \
           os.path.join(p.dirn,qc_info.fastq_dir) != p.fastq_dir:
            logger.warning("Stored fastq dir mismatch (%s != %s)" %
                           (p.fastq_dir,qc_info.fastq_dir))
        print "QC output dir: %s" % qc_dir
        print "-"*(len('Project: ')+len(p.name))
        print "%d samples | %d fastqs" % (len(p.samples),len(p.fastqs))
        # Verification step
        unverified = verify_qc(p,qc_dir)
        if unverified:
            if opts.list_unverified:
                for fq in unverified:
                    print fq
            print "Verification: FAILED"
            retval = 1
            if not opts.force:
                continue
        else:
            print "Verification: OK"
            if opts.verify:
                continue
        # Report generation
        qc_base = os.path.basename(qc_dir)
        if opts.filename is None:
            out_file = '%s_report.html' % qc_base
        else:
            out_file = opts.filename
        if not os.path.isabs(out_file):
            out_file = os.path.join(p.dirn,out_file)
        print "Writing QC report to %s" % out_file
        report_html= QCReporter(p).report(qc_dir=qc_dir,
                                          title=opts.title,
                                          filename=out_file)
        # Generate ZIP archive
        if opts.zip:
            report_zip = zip_report(p,report_html,qc_dir)
            print "ZIP archive: %s" % report_zip
        # MultiQC report
        if opts.multiqc:
            multiqc_report = os.path.join(p.dirn,
                                          "multi%s_report.html" %
                                          qc_base)
            multiqc_cmd = Command(
                'multiqc',
                '--title','%s' % opts.title,
                '--filename','%s' % multiqc_report,
                '--force',
                qc_dir)
            print "Running %s" % multiqc_cmd
            multiqc_retval = multiqc_cmd.run_subprocess()
            if multiqc_retval == 0 and os.path.exists(multiqc_report):
                print "MultiQC: %s" % multiqc_report
            else:
                print "MultiQC: FAILED"
                retval += 1
    # Finish with appropriate exit code
    sys.exit(retval)

if __name__ == '__main__':
    main()
