#!/usr/bin/env python
#
#     reportqc.py: generate report file for Illumina NGS qc runs
#     Copyright (C) University of Manchester 2015-2021 Peter Briggs
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
from auto_process_ngs.analysis import locate_project_info_file
from auto_process_ngs.command import Command
from auto_process_ngs.metadata import AnalysisProjectQCDirInfo
from auto_process_ngs.qc.constants import PROTOCOLS
from auto_process_ngs.qc.outputs import expected_outputs
from auto_process_ngs.qc.reporting import verify
from auto_process_ngs.qc.reporting import report
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
    data_dir_group = reporting.add_mutually_exclusive_group()
    data_dir_group.add_argument('--data-dir',action='store_true',
                                dest='use_data_dir',
                                help="create a data directory with copies "
                                "of QC artefacts needed for the HTML "
                                "report (NB data directory will always "
                                "be created for multi-project reports, "
                                "unless --no-data-dir is specified)")
    data_dir_group.add_argument('--no-data-dir',action='store_true',
                                dest='no_data_dir',
                                help="don't a data directory with copies "
                                "of QC artefacts (this is the default "
                                "except for multi-project reports)")
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
                   help="directory to report QC for; can be a project "
                   "directory (in which case the default QC directory "
                   "will be reported), or a QC directory within a "
                   "project")
    args = p.parse_args()

    # Report name and version
    print("%s version %s" % (os.path.basename(sys.argv[0]),__version__))

    # Report arguments
    if sys.argv[1:]:
        print("\n%s" % ' '.join(['"%s"' % arg if ' ' in arg else arg
                                 for arg in sys.argv[1:]]))

    # Report working directory
    print("\nCWD %s" % os.getcwd())

    # Check for MultiQC if required
    if args.multiqc:
        if find_program("multiqc") is None:
            logging.critical("MultiQC report requested but 'multiqc' "
                             "not available")
            sys.exit(1)

    # Get projects and QC dirs from supplied directories
    projects = []
    for d in args.dirs:
        print("\n**** Examining directory %s ****" % d)
        # Check if directory is a QC dir
        qc_dir = None
        # Look for 'qc.info' in current directory
        if os.path.exists(os.path.join(os.path.abspath(d),'qc.info')):
            print("...located 'qc.info', assuming this is QC dir")
            qc_dir = os.path.abspath(d)
            # Locate parent project dir
            metadata_file = locate_project_info_file(qc_dir)
            if metadata_file is not None:
                p = AnalysisProject(os.path.dirname(metadata_file))
                print("...located parent project: %s" % p.dirn)
            else:
                # Unable to locate project directory
                print("...failed to locate parent project metadata file")
                # Fall back to location of Fastq files
                qc_info = AnalysisProjectQCDirInfo(
                    os.path.join(qc_dir,'qc.info'))
                if qc_info.fastq_dir is not None:
                    project_dir = os.path.abspath(qc_info.fastq_dir)
                    if os.path.basename(project_dir).startswith('fastqs'):
                        # Use the next level up
                        project_dir = os.path.dirname(project_dir)
                    print("...putative parent project dir: %s (from "
                          " Fastq dir)" % project_dir)
                    p = AnalysisProject(project_dir)
                else:
                    # Failed to locate Fastqs
                    logging.fatal("Unable to locate parent project")
                    # Exit with an error
                    sys.exit(1)
            # Issue a warning if a QC dir was explicitly
            # specified on the command line
            if args.qc_dir is not None:
                logging.warning("--qc_dir has been ignored for this "
                                "directory")
        else:
            # Assume directory is a project
            p = AnalysisProject(os.path.abspath(d))
            print("...assuming this is a project dir")
            # Identify the QC directory
            if args.qc_dir is None:
                qc_dir = p.qc_dir
            else:
                qc_dir = args.qc_dir
            if not os.path.isabs(qc_dir):
                qc_dir = os.path.join(p.dirn,qc_dir)
            print("...QC directory: %s" % qc_dir)
        # Explicitly set the QC directory location)
        p.use_qc_dir(qc_dir)
        # Locate the Fastq dir
        qc_info = p.qc_info(qc_dir)
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
        print("...using Fastqs dir: %s" % p.fastq_dir)
        p.use_fastq_dir(fastq_dir,strict=False)
        projects.append(p)

    # Verify QC for projects
    print("\n**** Verifying QC ****")
    retval = 0
    report_projects = []
    for p in projects:
        print("\nProject: %s" % p.name)
        print("-"*(len('Project: ')+len(p.name)))
        print("%d sample%s | %d fastq%s" % (len(p.samples),
                                            's' if len(p.samples) != 1 else '',
                                            len(p.fastqs),
                                            's' if len(p.fastqs) != 1 else '',))
        # QC metadata
        qc_dir = p.qc_dir
        qc_info = p.qc_info(qc_dir)
        # Set QC protocol for verification
        if args.qc_protocol is None:
            protocol = qc_info.protocol
            if protocol is None:
                protocol = determine_qc_protocol(p)
        else:
            protocol = args.qc_protocol
        print("Verifying against QC protocol '%s'" % protocol)
        # Verification step
        if len(p.fastqs) == 0:
            logging.critical("No Fastqs!")
            verified = False
        else:
            try:
                verified = verify(p,qc_dir,protocol)
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
        report_projects.append(p)

    # Generate QC report
    if report_projects:
        # Set defaults from primary project
        p = report_projects[0]
        qc_base = os.path.basename(p.qc_dir)
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
                for p in report_projects:
                    multiqc_mtime = os.path.getmtime(multiqc_report)
                    for f in os.listdir(p.qc_dir):
                        if os.path.getmtime(os.path.join(p.qc_dir,f)) > \
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
                    '--force')
                for p in report_projects:
                    multiqc_cmd.add_args(p.qc_dir)
                print("\nRunning %s" % multiqc_cmd)
                multiqc_retval = multiqc_cmd.run_subprocess()
                if multiqc_retval == 0 and os.path.exists(multiqc_report):
                    print("MultiQC: %s\n" % multiqc_report)
                else:
                    print("MultiQC: FAILED")
                    retval += 1
            else:
                print("MultiQC: %s (already exists)\n" % multiqc_report)
        # Create data directory?
        use_data_dir = (len(projects) > 1)
        if args.use_data_dir:
            use_data_dir = True
        elif args.no_data_dir:
            use_data_dir = False
        # Generate report
        report_html = report(report_projects,
                             title=args.title,
                             filename=out_file,
                             relative_links=True,
                             use_data_dir=use_data_dir,
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
