#!/usr/bin/env python
#
#     audit_projects.py: summarise disk usage for sequencing projects
#     Copyright (C) University of Manchester 2015 Peter Briggs
#
"""
Summarise the disk usage for runs processed using 'auto_process'

"""

########################################################################
# Imports
#######################################################################

import optparse
import fnmatch
import os
import sys
import bcftbx.utils as utils
from auto_process_ngs.analysis import AnalysisDir

#######################################################################
# Functions
#######################################################################

def get_size(f):
    """
    Return size (in bytes) for file or directory
    
    This wraps the 'get_blocks' function and returns the
    number of blocks * 512.

    """
    return get_blocks(f)*512

def get_blocks(f):
    """
    Return number of 512-byte blocks for file or directory
    
    For a file, returns the 'st_blocks' value from the
    os.lstat() function.

    For a directory, returns the sum of all 'st_blocks' values
    for the directory contents (recursing into subdirectories
    as required).

    """
    blocks = os.lstat(f).st_blocks
    if not os.path.isfile(f):
        for dirpath,dirnames,filenames in os.walk(f):
            for d in dirnames:
                blocks += os.lstat(os.path.join(dirpath,d)).st_blocks
            for f in filenames:
                blocks += os.lstat(os.path.join(dirpath,f)).st_blocks
    return blocks

#######################################################################
# Main program
#######################################################################

if __name__ == "__main__":
    p = optparse.OptionParser(usage="%prog DIR [DIR...]",
                              description="Summarise the disk usage for runs that have "
                              "been processed using auto_process. The supplied DIRs are "
                              "directories holding the top-level analysis directories "
                              "corresponding to different runs. The program reports "
                              "total disk usage for projects assigned to each PI across "
                              "all DIRs.")
    p.add_option("--pi",action='store',dest="pi_name",default=None,
                 help="List data for PI(s) matching PI_NAME (can use glob-style "
                 "patterns)")
    p.add_option("--unassigned",action='store_true',dest="unassigned",default=False,
                 help="List data for projects where PI is not assigned")
    opts,args = p.parse_args()
    # Collect data
    audit_data = {}
    unassigned = []
    undetermined = []
    for d in args:
        for dirn in utils.list_dirs(d):
            dirn = os.path.join(d,dirn)
            #print "Examining %s" % dirn
            try:
                run = AnalysisDir(dirn)
                for p in run.get_projects():
                    if p.name == "undetermined":
                        undetermined.append((p,get_size(p.dirn)))
                        continue
                    pi = p.info.PI
                    if pi is None:
                        # PI is not assigned
                        p.info['run'] = os.path.basename(dirn)
                        unassigned.append(p)
                        continue
                    elif opts.pi_name is not None:
                        if not fnmatch.fnmatch(pi,opts.pi_name):
                            # Skip non-matching name
                            continue
                    #print "\t%s: %s" % (p.name,pi)
                    if pi not in audit_data:
                        audit_data[pi] = []
                    # Acquire size of data
                    size = get_size(p.dirn)
                    if p.fastqs_are_symlinks:
                        # Actual fastq files are outside the project
                        # and need to be explicitly added
                        for fq in p.fastqs:
                            try:
                                size += get_size(utils.Symlink(fq).resolve_target())
                                #if os.path.islink(fq):
                                #    s = utils.Symlink(fq)
                                #    size += get_size(s.resolve_target())
                            except Exception,ex:
                                print "Failed to get size for fastq file: %s" % fq
                                print "%s" % ex
                                sys.exit(1)
                    # Add to list of projects
                    audit_data[pi].append((p,size))
            except Exception,ex:
                print "Failed to load as run: %s" % ex
                pass
    # Sort into order by disk usage
    pi_list = audit_data.keys()
    pi_list = sorted(pi_list,key=lambda x: sum([y[1] for y in audit_data[x]]),
                     reverse=True)
    # Report the unassigned projects, if requested
    if opts.unassigned:
        print "%d unassigned projects:" % len(unassigned)
        unassigned = sorted(unassigned,key=lambda x: str(x).split('_')[0])
        for project in unassigned:
            print "%s: %s" % (project.name,project.info.run)
        sys.exit(0)
    # Report if no PIs were found
    if len(pi_list) == 0:
        print "No projects assigned to PIs found"
        sys.exit(0)
    # Report PIs, projects etc
    print "Summary (PI, # of projects, total usage):"
    print "========================================="
    total_projects = 0
    total_size = 0
    for pi in pi_list:
        n_projects = len(audit_data[pi])
        size = sum([p[1] for p in audit_data[pi]])
        print "%s\t%d\t%s" % (pi,n_projects,
                              utils.format_file_size(size))
        total_projects += n_projects
        total_size += size
    print "Total usage\t%d\t%s" % (total_projects,
                                   utils.format_file_size(total_size))
    print "\nBreakdown by PI/project:"
    print "========================"
    for pi in pi_list:
        print "%s:" % pi
        for project,size in audit_data[pi]:
            print "\t%s:\t%s\t%s" % (project.info.run,project.name,
                                     utils.format_file_size(size))
    if undetermined:
        print "\nUsage for 'undetermined' reads:"
        print "==============================="
        total_size = 0
        for u,size in undetermined:
            print "%s\t%s" % (u.info.run,utils.format_file_size(size))
            total_size += size
        print "Total usage\t%s" % utils.format_file_size(total_size)
        
