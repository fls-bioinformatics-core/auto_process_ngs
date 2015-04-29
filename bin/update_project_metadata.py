#!/bin/env python
#
#     update_project_metadata.py: update data associated with a project
#     Copyright (C) University of Manchester 2015 Peter Briggs
#
"""
Update the metadata associated with project(s) in an Illumina run analysis

"""

########################################################################
# Imports
#######################################################################

import optparse
import fnmatch
import os
import sys
import bcftbx.utils as utils
import bcftbx.platforms as platforms
from auto_process_ngs.utils import AnalysisDir

#######################################################################
# Main program
#######################################################################

if __name__ == "__main__":
    p = optparse.OptionParser(usage="%prog DIR [PROJECT]")
    p.add_option("-i","--init",action="store_true",dest="init",default=False,
                 help="initialise metadata file for the selected project (nb can "
                 "only be applied to one project at a time)")
    p.add_option("-u","--update",action="append",dest="update",
                 help="update the metadata in the selected project by specifying "
                 "key=value pairs e.g. user='Peter Briggs' (nb can only be applied "
                 "to one project at a time)")
    opts,args = p.parse_args()
    if not args or len(args) > 2:
        p.error("Need to supply a run directory and an optional project name")
    # Try and load the data
    try:
        run = AnalysisDir(args[0])
        print "Loaded data from %s" % args[0]
        print "Found %d projects" % run.n_projects
        if run.undetermined is not None:
            print "Found 'undetermined' analysis"
        print "Found %d sequencing data directories" % run.n_sequencing_data
    except Exception,ex:
        sys.stderr.write("Failed to load data for %s: %s\n" % (args[0],ex))
        sys.exit(1)
    # Aquire the projects
    try:
        target_project = args[1]
    except IndexError:
        target_project = None
    projects = run.get_projects(pattern=target_project)
    if not projects:
        print "No projects found"
        sys.exit(1)
    print "%d projects selected" % len(projects)
    if (opts.update or opts.init) and len(projects) > 1:
        sys.stderr.write("-i/-u can't be used for multiple projects\n")
        sys.exit(1)
    # Check metadata for the projects
    print "Checking metadata for projects"
    for p in projects:
        # Check for metadata file
        has_metadata_file = os.path.exists(p.info_file)
        if not has_metadata_file:
            # No metadata file
            print "%s: missing metadata file '%s'" % (p.name,p.info_file)
            if opts.init:
                print "Initialising metadata file"
                run_name = os.path.basename(os.path.abspath(args[0]))
                if run_name.endswith('_analysis'):
                    run_name = run_name[:-9]
                platform = platforms.get_sequencer_platform(run_name)
                p.info['run'] = run_name
                p.info['platform'] = platform
                p.info['samples'] = p.prettyPrintSamples()
                p.info['paired_end'] = run.paired_end
                p.info.save(p.info_file)
        else:
            # Check for metadata that is not assigned
            missing_data = []
            for attr in ('run','user','PI','samples','library_type','organism','platform'):
                if p.info[attr] is None:
                    missing_data.append(attr)
            if missing_data:
                print "%s: missing values for %s" % (p.name,', '.join(["'%s'" % x 
                                                                       for x in missing_data]))
            else:
                print "%s: ok" % p.name
        # Do updates
        if opts.update:
            if not has_metadata_file:
                sys.stderr.write("No metadata file: use -i/--init to create one first\n")
                sys.exit(1)
            for pair in opts.update:
                try:
                    key = pair.split('=')[0]
                    value = '='.join(pair.split('=')[1:]).strip("'\"")
                except Exception,ex:
                    sys.stderr.write("Bad key=value pair: %s\n" % pair)
                    sys.exit(1)
                if key not in ('run','platform','user','PI','library_type','organism','comments'):
                    sys.stderr.write("Unable to update unrecognised item '%s'\n" % key)
                else:
                    print "Setting '%s' to '%s'" % (key,value)
                    p.info[key] = value
            # Save the metadata
            dry_run = False
            if not dry_run:
                p.info.save()


    
