#!/usr/bin/env python
import optparse
import difflib
import bcftbx.IlluminaData as IlluminaData

def predict_samplesheet_outputs(sample_sheet_file):
    """
    """
    sample_sheet = IlluminaData.SampleSheet(sample_sheet_file)
    projects = {}
    s_index = 0
    for line in sample_sheet.data:
        # Get project
        project = str(line[sample_sheet.sample_project_column])
        if project not in projects:
            projects[project] = {}
        # Get data on sample
        sample_id = str(line[sample_sheet.sample_id_column])
        try:
            sample_name = str(line[sample_sheet.sample_name_column])
        except KeyError:
            sample_name = ''
        if not sample_id:
            sample_id = sample_name
        if sample_id not in projects[project]:
            projects[project][sample_id] = {}
            projects[project][sample_id]['barcodes'] = {}
        # Barcodes
        indx = get_barcode_from_samplesheet_line(line)
        if not indx:
            indx = "NoIndex"
        # Lane
        if sample_sheet.has_lanes:
            lane = line['Lane']
        else:
            lane = None
        if indx not in projects[project][sample_id]['barcodes']:
            # Increment sample (S) index for bcl2fastq2
            s_index += 1
            # Store data
            projects[project][sample_id]['s_index'] = s_index
            projects[project][sample_id]['barcodes'][indx] = []
        if sample_sheet.has_lanes:
            projects[project][sample_id]['barcodes'][indx].append(lane)
    return projects

def get_barcode_from_samplesheet_line(line,delimiter='-'):
    """
    """
    indx = None
    try:
        # Try dual-indexed IEM4 format
        indx = "%s%s%s" %(line['index'].strip(),
                          delimiter,
                          line['index2'].strip())
    except KeyError:
        # Try single indexed IEM4 (no index2)
        try:
            indx = line['index'].strip()
        except KeyError:
            # Try CASAVA format
            indx = line['Index'].strip()
    return indx

def get_close_names(names):
    """
    """
    close_names = {}
    for i,name in enumerate(names[:-1]):
        matches = difflib.get_close_matches(name,names[i+1:])
        if len(matches):
            #print "*** %s: matches %s ***" % (name,matches)
            close_names[name] = matches
            for name1 in matches:
                if name1 not in close_names:
                    close_names[name1] = []
                close_names[name1].append(name)
    return close_names

if __name__ == "__main__":
    p = optparse.OptionParser()
    opts,args = p.parse_args()
    # Get prediction
    projects = predict_samplesheet_outputs(args[0])
    project_names = projects.keys()
    project_names.sort()
    # Check for close-matching project names
    close_names = get_close_names(project_names)
    # Report project names
    print
    print "Expecting %d projects:" % len(project_names)
    print
    warning_close_names = False
    for p in project_names:
        warning_flag = False
        if p in close_names:
            warning_flag = True
            warning_close_names = True
        print "%s %s" % (('*' if warning_flag else '-'),p)
    print
    # Report samples within each project
    warning_wrong_barcodes = False
    for p in project_names:
        samples = projects[p]
        print "%s\n%s" % (p,'-'*len(p))
        sample_names = samples.keys()
        sample_names.sort()
        for s in sample_names:
            warning_wrong_barcodes = warning_wrong_barcodes or \
                                     (len(samples[s]['barcodes']) > 1)
            warning_flag = (len(samples[s]['barcodes']) > 1)
            s_index = samples[s]['s_index']
            for barcode in samples[s]['barcodes']:
                #print samples[s][barcode]
                lanes = ','.join([str(x)
                                  for x in samples[s]['barcodes'][barcode]])
                print "%s%s\tS%s\t%s\t%s" % (('* ' if warning_flag else ''),
                                             s,
                                             s_index,
                                             barcode,
                                             (('L' + lanes) if lanes else ''))
        print
    # Issue warnings
    if warning_close_names:
        print "*** WARNING projects detected with closing-matching names ***"
    if warning_wrong_barcodes:
        print "*** WARNING samples detected with multiple barcodes ***"
