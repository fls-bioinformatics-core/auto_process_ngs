#!/usr/bin/env python
#
#     manage_runs.py: utility for managing fastq files from auto_process
#     Copyright (C) University of Manchester 2014 Peter Briggs
#
#########################################################################
#
# manage_fastqs.py
#
#########################################################################

"""
Utility for managing Fastq files from auto_processing.

Functionality includes:

- report individal file sizes
- copy to another location
- create checksums
- make zip archive

"""

#######################################################################
# Imports
#######################################################################

import sys
import os
import argparse
import shutil
import tempfile
import zipfile
import fnmatch
import bcftbx.utils as bcf_utils
import bcftbx.Md5sum as md5sum
import auto_process_ngs.utils as utils
import auto_process_ngs.applications as applications
from auto_process_ngs.analysis import AnalysisDir
from auto_process_ngs import get_version

#######################################################################
# Functions
#######################################################################

def get_fastqs(project,pattern=None):
    """Return fastq files within an AnalysisProject

    Given an AnalysisProject, yields

    (sample_name,fastq,actual_fastq)

    tuples for each fastq file in all samples, where

    - 'fastq' = the original fastq name
    - 'actual_fastq' = the 'actual' fastq location (if
      'fastq' is a symbolic link)

    Arguments:
      project: AnalysisProject instance
      pattern: if supplied then use the supplied pattern
        to filter fastqs based on filename

    """
    for sample in project.samples:
        for fq in sample.fastq:
            # Filter on name
            if pattern is not None:
                if not fnmatch.fnmatch(os.path.basename(fq),pattern):
                    continue
            # Resolve links
            if os.path.islink(fq):
                target = bcf_utils.Symlink(fq).resolve_target()
            else:
                target = fq
            yield (sample.name,fq,target)

def write_checksums(project,pattern=None,filen=None,relative=True):
    """Write MD5 checksums for fastq files with an AnalysisProject

    Arguments:
      project: AnalysisProject instance
      pattern: if supplied then use the supplied pattern
        to filter fastqs based on filename
      filen: if supplied then checksums will be written
        to this file; otherwise they will be written to
        stdout (default)
      relative: if True (default) then fastq file names
        will be the basename; otherwise they will be the
        full paths.

    """
    if filen:
        fp = open(md5file,'w')
    else:
        fp = sys.stdout
    for sample_name,fastq,fq in get_fastqs(project,pattern=pattern):
        if relative:
            name = os.path.basename(fq)
        else:
            name = fq
        fp.write("%s  %s\n" % (md5sum.md5sum(fq),name))
    if filen:
        fp.close()

def copy_to_dest(f,dirn,chksum=None):
    """Copy a file to a local or remote destination

    Raises an exception if the copy operation fails.

    If 'chksum' argument is supplied then the MD5 sum of
    the copy is also verified against this and an
    exception is raised if this fails to match.

    Arguments:
      f: file to copy (must be local)
      dirn: target directory, either local or of the form
        "[user@]host:dir"
      chksum: (optional) MD5 sum of the original file
        to match against the copy
    
    """
    if not os.path.exists(f):
        raise Exception("File %s doesn't exist" % f)
    user,host,dest = utils.split_user_host_dir(dirn)
    remote = (host is not None)
    if not remote:
        # Local copy
        shutil.copy(f,dirn)
        if chksum is not None:
            if md5sum.md5sum(f) != chksum:
                raise Exception("MD5 checksum failed for copy of %s" % f)
    else:
        # Remote copy
        try:
            scp = applications.general.scp(user,host,f,dest)
            print("Running %s" % scp)
            scp.run_subprocess()
            # Run md5sum -c on the remote system
            if chksum is not None:
                md5sum_check = applications.general.ssh_command(
                    user,host,
                    ('echo',
                     '"%s  %s"' % (chksum,
                                   os.path.join(dest,os.path.basename(f))),
                    '|','md5sum','-c'))
                print("Running %s" % md5sum_check)
                md5sum_check.run_subprocess()
        except Exception, ex:
            raise Exception("Failed to copy %s to %s: %s" % (f,dirn,ex))

#######################################################################
# Main program
#######################################################################

if __name__ == "__main__":

    # Command line
    p = argparse.ArgumentParser(
        usage="\n\t%(prog)s DIR"
        "\n\t%(prog)s DIR PROJECT"
        "\n\t%(prog)s DIR PROJECT copy [[user@]host:]DEST"
        "\n\t%(prog)s DIR PROJECT md5"
        "\n\t%(prog)s DIR PROJECT zip",
        description="Fastq management utility. If only DIR is "
        "supplied then list the projects; if PROJECT is supplied "
        "then list the fastqs; 'copy' command copies fastqs for the "
        "specified PROJECT to DEST on a local or remote server; 'md5' "
        "command generates checksums for the fastqs; 'zip' command "
        "creates a zip file with the fastq files.",
        version="%(prog)s "+get_version())
    p.add_argument('--filter',action='store',dest='pattern',
                   default=None,
                   help="filter file names for reporting and copying "
                   "based on PATTERN")
    p.add_argument('--fastq_dir',action='store',dest='fastq_dir',
                   default=None,
                   help="explicitly specify subdirectory of DIR with "
                   "Fastq files to run the QC on.")
    options,args = p.parse_known_args()
    # Get analysis dir
    try:
        dirn = args[0]
        print("Loading data for analysis dir %s" % dirn)
    except IndexError:
        p.error("Need to supply the path to an analysis dir")
        sys.exit(1)
    analysis_dir = AnalysisDir(dirn)
    # Get specified project
    try:
        project_name = args[1]
    except IndexError:
        # List projects and exit
        print("Projects:")
        for project in analysis_dir.projects:
            print("%s" % project.name)
            if len(project.fastq_dirs) > 1:
                # List the fastq sets if there are more than one
                # and flag the primary set with an asterisk
                for d in project.fastq_dirs:
                    is_primary = (d == project.info.primary_fastq_dir)
                    print("- %s%s" % (d,
                                      (" *" if is_primary else "")))
        if analysis_dir.undetermined:
            print("_undetermined")
        sys.exit(0)
    sys.stdout.write("Checking for project '%s'..." % project_name)
    project = None
    for prj in analysis_dir.projects:
        if prj.name == project_name:
            project = prj
            break
    if project is None:
        if project_name == "_undetermined" and analysis_dir.undetermined:
            project = analysis_dir.undetermined
        else:
            print("not found")
            sys.stderr.write("FAILED cannot find project '%s'\n" % project_name)
            sys.exit(1)
    print("ok")

    # Switch to requested Fastq set
    if options.fastq_dir is not None:
        try:
            project.use_fastq_dir(fastq_dir=options.fastq_dir)
        except Exception as ex:
            sys.stderr.write("ERROR %s\n" % ex)
            sys.stderr.write("FAILED unable to switch to fastq set '%s'\n"
                             % options.fastq_dir)
            sys.exit(1)

    # Filter fastqs on pattern
    if options.pattern is not None:
        print("Filtering fastqs using pattern '%s'" % options.pattern)
    #    fastqs = [fq for fq in fastqs
    #              if fnmatch.fnmatch(os.path.basename(fq[2]),options.pattern)]
    # Check for a command
    try:
        cmd = args[2]
    except IndexError:
        # List fastqs and exit
        total_size = 0
        n_fastqs = 0
        sample_names = set()
        # Collect information
        fastq_set = os.path.relpath(project.fastq_dir,project.dirn)
        print("Fastq set: %s%s" % (
            ("default" if fastq_set == "fastqs" else fastq_set),
            (" (primary)"
             if fastq_set == project.info.primary_fastq_dir else "")))
        for sample_name,fastq,fq in get_fastqs(project,pattern=options.pattern):
            # File size
            fsize = os.lstat(fq).st_size
            print("%s\t%s%s\t%s" % (sample_name,
                                    os.path.basename(fq),
                                    ('*' if os.path.islink(fastq) else ''),
                                    bcf_utils.format_file_size(fsize)))
            sample_names.add(sample_name)
            total_size += fsize
            n_fastqs += 1
        # Summary
        print("Total:\t%s" % bcf_utils.format_file_size(total_size))
        print("%d %ssamples" % (len(sample_names),
                                ('paired-end '
                                 if project.info.paired_end else '')))
        print("%d fastqs" % n_fastqs)
        sys.exit(0)
    # Perform command
    if cmd not in ('copy','zip','md5'):
        p.error("Unrecognised command '%s'\n" % cmd)
        sys.exit(1)
    if cmd == 'copy':
        # Get the destination
        try:
            dest = args[3]
        except IndexError:
            p.error("Need to supply a destination for 'copy' command")
            sys.exit(1)
        # Make a temporary MD5 file
        tmp = tempfile.mkdtemp()
        try:
            md5file = os.path.join(tmp,"%s.chksums" % project.name)
            sys.stdout.write("Creating checksum file %s..." % md5file)
            write_checksums(project,pattern=options.pattern,filen=md5file)
            print("done")
            print("Copying to %s" % dest)
            copy_to_dest(md5file,dest)
            # Load checksums into dictionary
            chksums = dict()
            with open(md5file,'r') as fp:
                for line in fp:
                    chksum,filen = line.strip('\n').split()
                    chksums[filen] = chksum
        finally:
            shutil.rmtree(tmp)
        # Copy fastqs
        nfastqs = sum(1 for _ in get_fastqs(project,pattern=options.pattern))
        i = 0
        for sample_name,fastq,fq in get_fastqs(project,pattern=options.pattern):
            i += 1
            print("(% 2d/% 2d) %s" % (i,nfastqs,fq))
            copy_to_dest(fq,dest,chksums[os.path.basename(fq)])
    elif cmd == 'md5':
        # Generate MD5 checksums
        md5file = "%s.chksums" % project.name
        if os.path.exists(md5file):
            sys.stderr.write("ERROR checksum file '%s' already exists\n" % md5file)
            sys.exit(1)
        sys.stdout.write("Creating checksum file %s..." % md5file)
        write_checksums(project,pattern=options.pattern,filen=md5file)
        print("done")
    elif cmd == 'zip':
        # Create a zip file
        zip_file = "%s.zip" % project.name
        if os.path.exists(zip_file):
            sys.stderr.write("ERROR zip file '%s' already exists" % zip_file)
            sys.exit(1)
        print("Creating zip file %s" % zip_file)
        zz = zipfile.ZipFile(zip_file,'w',allowZip64=True)
        # Add fastqs
        for sample_name,fastq,fq in get_fastqs(project,pattern=options.pattern):
            zz.write(fq,arcname=os.path.basename(fq))
        # Make a temporary MD5 file
        tmp = tempfile.mkdtemp()
        try:
            md5file = os.path.join(tmp,"%s.chksums" % project.name)
            sys.stdout.write("Creating checksum file %s..." % md5file)
            write_checksums(project,filen=md5file)
            print("done")
            print("Adding to %s" % zip_file)
            zz.write(md5file,arcname=os.path.basename(md5file))
        finally:
            shutil.rmtree(tmp)
        zz.close()
        
