#!/usr/bin/env python
#
#     reportqc.py: generate report file for Illumina NGS qc runs
#     Copyright (C) University of Manchester 2015-2020 Peter Briggs
#

#######################################################################
# Imports
#######################################################################

import sys
import os
import argparse
import logging
from bcftbx.utils import find_program
from auto_process_ngs.analysis import AnalysisProject
from auto_process_ngs.applications import Command
from auto_process_ngs.qc.constants import PROTOCOLS
from auto_process_ngs.qc.reporting import QCReporter
from auto_process_ngs.qc.outputs import expected_outputs
from auto_process_ngs.qc.utils import determine_qc_protocol
from auto_process_ngs import get_version

__version__ = get_version()

"""
reportqc

Utility to verify and report on QC outputs from
auto_process pipeline.
"""

#######################################################################
# Main program
#######################################################################

def main():
    # Deal with command line
    p = argparse.ArgumentParser(description="Generate QC report for each "
                                "directory DIR")
    p.add_argument('-v','--version',action='version',
                   version="%(prog)s "+__version__)
    p.add_argument('--protocol',
                   action='store',dest='qc_protocol',default=None,
                   help="explicitly specify QC protocol (must be one of "
                   "%s). Default is to determine the protocol "
                   "automatically (recommended)" %
                   str(','.join(["'%s'" % pr for pr in PROTOCOLS])))
    p.add_argument('--qc_dir',action='store',dest='qc_dir',default='qc',
                   help="explicitly specify QC output directory (nb if "
                   "supplied then the same QC_DIR will be used for each "
                   "DIR. Non-absolute paths are assumed to be relative to "
                   "DIR). Default: 'qc'")
    p.add_argument('--fastq_dir',
                   action='store',dest='fastq_dir',default=None,
                   help="explicitly specify subdirectory of DIRs with "
                   "Fastq files to run the QC on")
    reporting = p.add_argument_group('Reporting options')
    reporting.add_argument('-t','--title',action='store',dest='title',
                           default=None,
                           help="title for output QC reports")
    reporting.add_argument('-f','--filename',
                           action='store',dest='filename',default=None,
                           help="file name for output HTML QC report "
                           "(default: <DIR>/<QC_DIR>_report.html)")
    reporting.add_argument('--zip',action='store_true',
                           dest='zip',default=False,
                           help="make ZIP archive for the QC report")
    reporting.add_argument('--multiqc',action='store_true',
                           dest='multiqc',default=False,
                           help="generate MultiQC report")
    reporting.add_argument('--force',action='store_true',
                           dest='force',default=False,
                           help="force generation of reports even if "
                           "verification fails")
    verification = p.add_argument_group('Verification options')
    verification.add_argument('--verify',action='store_true',dest='verify',
                              help="verify the QC products only (don't "
                              "write the report); returns exit code 0 "
                              "if QC is verified, 1 if not")
    deprecated = p.add_argument_group('Deprecated options')
    deprecated.add_argument('-l','--list-unverified',action='store_true',
                            dest='list_unverified',default=False,
                            help="deprecated: does nothing (Fastqs with "
                            "missing QC outputs can no longer be listed)")
    deprecated.add_argument('--strand_stats',action='store_true',
                            dest='fastq_strand',default=False,
                            help="deprecated: does nothing (strand stats "
                            "are automatically included if present)")
    p.add_argument('dirs',metavar="DIR",nargs='+',
                   help="directory to report QC for")
    args = p.parse_args()

    # Report name and version
    print("%s version %s" % (os.path.basename(sys.argv[0]),__version__))

    # Check for MultiQC if required
    if args.multiqc:
        if find_program("multiqc") is None:
            logging.critical("MultiQC report requested but 'multiqc' "
                             "not available")
            sys.exit(1)

    # Examine projects i.e. supplied directories
    retval = 0
    for d in args.dirs:
        dir_path = os.path.abspath(d)
        p = AnalysisProject(dir_path)
        print("Project           : %s" % p.name)
        print("Primary Fastqs dir: %s" % p.fastq_dir)
        if args.qc_dir is None:
            qc_dir = p.qc_dir
        else:
            qc_dir = args.qc_dir
        if not os.path.isabs(qc_dir):
            qc_dir = os.path.join(p.dirn,qc_dir)
        qc_base = os.path.basename(qc_dir)
        # QC metadata
        qc_info = p.qc_info(qc_dir)
        print("QC output dir     : %s" % qc_dir)
        print("...stored protocol: %s" % qc_info.protocol)
        print("...stored Fastqdir: %s" % qc_info.fastq_dir)
        # Set QC protocol
        if args.qc_protocol is None:
            protocol = qc_info.protocol
            if protocol is None:
                protocol = determine_qc_protocol(p)
        else:
            protocol = args.qc_protocol
        print("QC protocol: %s" % protocol)
        # Fastq subdirectory
        if args.fastq_dir is None:
            fastq_dir = qc_info.fastq_dir
            if fastq_dir is None:
                fastq_dir = p.fastq_dir
        else:
            fastq_dir = args.fastq_dir
            if qc_info.fastq_dir is not None:
                if os.path.join(p.dirn,qc_info.fastq_dir) != fastq_dir:
                    logging.warning("Stored fastq dir mismatch "
                                    "(%s != %s)" % (fastq_dir,
                                                    qc_info.fastq_dir))
        p.use_fastq_dir(fastq_dir)
        print("Fastqs dir        : %s" % p.fastq_dir)
        print("-"*(len('Project: ')+len(p.name)))
        print("%d samples | %d fastqs" % (len(p.samples),len(p.fastqs)))
        # Setup reporter
        qc_reporter = QCReporter(p)
        # Verification step
        if len(p.fastqs) == 0:
            logging.critical("No Fastqs!")
            verified = False
        else:
            try:
                verified = qc_reporter.verify(qc_dir,protocol)
            except Exception as ex:
                logging.critical("Error: %s" % ex)
                verified = False
        if not verified:
            print("Verification: FAILED")
            if not args.force:
                retval = 1
                continue
            else:
                print("--force specified, ignoring previous errors")
        else:
            print("Verification: OK")
            if args.verify:
                continue
        # Filename and location for report
        if args.filename is None:
            out_file = '%s_report.html' % qc_base
        else:
            out_file = args.filename
        if not os.path.isabs(out_file):
            out_file = os.path.join(p.dirn,out_file)
        out_dir = os.path.dirname(out_file)
        # MultiQC report
        if args.multiqc:
            multiqc_report = os.path.join(out_dir,
                                          "multi%s_report.html" %
                                          qc_base)
            # Check if we need to rerun MultiQC
            if os.path.exists(multiqc_report) and not args.force:
                run_multiqc = False
                multiqc_mtime = os.path.getmtime(multiqc_report)
                for f in os.listdir(qc_dir):
                    if os.path.getmtime(os.path.join(qc_dir,f)) > \
                       multiqc_mtime:
                        # Input is newer than report
                        run_multiqc = True
                        break
            else:
                run_multiqc = True
            # (Re)run MultiQC
            if run_multiqc:
                multiqc_cmd = Command(
                    'multiqc',
                    '--title','%s' % args.title,
                    '--filename','%s' % multiqc_report,
                    '--force',
                    qc_dir)
                print("Running %s" % multiqc_cmd)
                multiqc_retval = multiqc_cmd.run_subprocess()
                if multiqc_retval == 0 and os.path.exists(multiqc_report):
                    print("MultiQC: %s" % multiqc_report)
                else:
                    print("MultiQC: FAILED")
                    retval += 1
            else:
                print("MultiQC: %s (already exists)" % multiqc_report)
        # Generate report
        report_html= qc_reporter.report(qc_dir=qc_dir,
                                        qc_protocol=protocol,
                                        title=args.title,
                                        filename=out_file,
                                        relative_links=True,
                                        make_zip=args.zip)
        print("Wrote QC report to %s" % out_file)
    # Finish with appropriate exit code
    print("%s completed: exit code %s (%s)" %
          (os.path.basename(sys.argv[0]),
           retval,
           ('ok' if retval == 0 else 'error')))
    sys.exit(retval)

if __name__ == '__main__':
    main()
