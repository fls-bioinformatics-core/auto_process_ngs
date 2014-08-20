#!/bin/env python
#
#     manage_runs.py: utility for managing Illumina sequencing run directories
#     Copyright (C) University of Manchester 2014 Peter Briggs
#
#########################################################################
#
# manage_runs.py
#
#########################################################################

"""
Utility for managing Illumina sequencing run directories

"""

#######################################################################
# Imports
#######################################################################
import os
import sys
import logging
import shutil
import Md5sum
import bcf_utils
import auto_process_ngs.settings as settings
import auto_process_ngs.utils as utils
import auto_process_ngs.applications as applications
from auto_process_ngs.parser import CommandParser
from auto_process_ngs.parser import add_debug_option,add_dry_run_option

__version__ = settings.version

#######################################################################
# Functions
#######################################################################

def get_projects(seq_data,pattern=None):
    projects = []
    for p in seq_data.projects:
        if bcf_utils.name_matches(p.name,pattern):
            projects.append(p)
    return projects

def get_samples(seq_data,project_pattern=None,sample_pattern=None):
    samples = []
    for p in get_projects(seq_data,pattern=project_pattern):
        for s in p.samples:
            if bcf_utils.name_matches(s.name,sample_pattern):
                samples.append(s)
    return samples

def get_fastqs(seq_data,project_pattern=None,sample_pattern=None):
    fastqs = []
    for s in get_samples(seq_data,
                         project_pattern=project_pattern,
                         sample_pattern=sample_pattern):
        for fastq in s.fastq:
            fastqs.append(os.path.join(s.dirn,fastq))
    return fastqs

def copy_to_dest(f,dirn):
    """Copy a file to a local or remote destination

    Raises an exception if the copy operation fails.

    Arguments:
      f: file to copy (must be local)
      dirn: target directory, either local or of the form
        "[user@]host:dir"
    
    """
    if not os.path.exists(f):
        raise Exception("File %s doesn't exist" % f)
    user,host,dest = utils.split_user_host_dir(dirn)
    remote = (host is not None)
    if not remote:
        # Local copy
        shutil.copy(f,dirn)
    else:
        # Remote copy
        try:
            scp = applications.general.scp(user,host,f,dest)
            print "Running %s" % scp
            scp.run_subprocess()
        except Exception, ex:
            raise Exception("Failed to copy %s to %s: %s" % (f,dirn,ex))

#######################################################################
# Main program
#######################################################################

if __name__ == "__main__":
    # Set up command line parser
    p = CommandParser(description="Utility for managing processed and analysed Illumina "
                      "sequence data in ANALYSIS_DIR",
                      version="%prog "+__version__)
    # Add info command
    p.add_command('info',help="Get information about ANALYSIS_DIR",
                  usage="%prog info [OPTIONS] ANALYSIS_DIR",
                  description="Report information on processed Illumina "
                  "sequence data in ANALYSIS_DIR.")
    add_debug_option(p.parser_for('info'))
    # Add copy command
    p.add_command('copy',help="Copy fastqs from ANALYSIS_DIR",
                  usage="%prog copy [OPTIONS] ANALYSIS_DIR DEST_DIR",
                  description="Copy fastqs from ANALYSIS_DIR to DEST_DIR.")
    p.parser_for('copy').add_option('--projects',action='store',dest='projects',default=None,
                                    help="Restrict copying to projects matching the "
                                    "supplied pattern")
    add_dry_run_option(p.parser_for('copy'))
    add_debug_option(p.parser_for('copy'))
    # Process the command line
    cmd,options,args = p.parse_args()
    if len(args) < 1:
        p.error("Need to supply a directory to examine")
    if options.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    # Acquire data for putative analysis directory
    dirn = os.path.abspath(args[0])
    print "Examining %s" % dirn
    analysis_dir = utils.AnalysisDir(dirn)
    if not (analysis_dir.n_projects and analysis_dir.n_sequencing_data):
        logging.warning("Not an analysis directory")
        sys.exit()
    # Deal with commands
    if cmd == 'info':
        # Report projects from sequencing data
        if analysis_dir.n_sequencing_data:
            for data in analysis_dir.sequencing_data:
                print "Sequencing data projects found in '%s':\n" % \
                    os.path.basename(data.unaligned_dir)
                print "Project_name\t#smpl\tSample_names"
                for project in data.projects:
                    line = [project.name,
                            len(project.samples),
                            project.prettyPrintSamples()]
                    print "%s" % '\t'.join([str(x) for x in line])
                print ""
        # Report analysis project directories
        if analysis_dir.n_projects:
            print "Putative analysis project directories:\n"
            print "Project_name\t#smpl\tQC?\tSample_names"
            for project in analysis_dir.projects:
                line = [project.name,
                        len(project.samples),
                        ('ok' if project.verify_qc() else '?'),
                        project.prettyPrintSamples()]
                print "%s" % '\t'.join([str(x) for x in line])
            print ""
    elif cmd == 'copy':
        # Copy fastqs
        if len(args) != 2:
            p.error("Need to supply a destination directory")
        # Check destination
        user,host,dest_dir = utils.split_user_host_dir(args[1])
        if host is None:
            dest_dir = os.path.abspath(dest_dir)
            print "Copying to local directory: %s" % dest_dir
        else:
            dest_dir = args[1]
            print "Copying to remote directory: %s" % dest_dir
        # Fetch fastqs
        if options.projects is None:
            project_pattern = '*'
            sample_pattern = '*'
        else:
            project_pattern = options.projects.split('/')[0]
            try:
                sample_pattern = options.projects.split('/')[1]
            except IndexError:
                sample_pattern = '*'
        for data in analysis_dir.sequencing_data:
            fastqs = get_fastqs(data,project_pattern=project_pattern,
                                sample_pattern=sample_pattern)
        # File sizes
        total_size = 0
        for fq in fastqs:
            fsize = os.lstat(fq).st_size
            total_size += fsize
            print "%s\t%s" % (os.path.basename(fq),
                              bcf_utils.format_file_size(fsize))
        print "Total: %s" % bcf_utils.format_file_size(total_size)
        # Generate MD5 sums
        for fq in fastqs:
            chksum = Md5sum.md5sum(fq)
            print "%s  %s" % (chksum,os.path.basename(fq))
        # Copy them
        for fq in fastqs:
            if not options.dry_run:
                copy_to_dest(fq,args[1])
            else:
                print "%s" % fq

