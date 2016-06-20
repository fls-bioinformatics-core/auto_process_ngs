#!/usr/bin/env python
import optparse
import difflib
import bcftbx.IlluminaData as IlluminaData

# Classes for sample sheet prediction and validation

# FIXME not tested on pre IEM sample sheets
# TODO implement fastq prediction
# TODO implement bcl2fastq output validation

class SampleSheetPredictor(object):
    """
    Sample sheet output predictor

    Processes a sample sheet file and builds a data
    structure which predicts the outputs.

    """
    def __init__(self,sample_sheet_file):
        self.projects = []
        sample_sheet = IlluminaData.SampleSheet(sample_sheet_file)
        s_index = 0
        for line in sample_sheet:
            # Get project and sample info
            project_name = str(line[sample_sheet.sample_project_column])
            sample_id = str(line[sample_sheet.sample_id_column])
            try:
                sample_name = str(line[sample_sheet.sample_name_column])
            except KeyError:
                sample_name = None
            project = self.add_project(project_name)
            sample = project.add_sample(sample_id,
                                        sample_name=sample_name)
            # Barcode and lane
            indx = get_barcode_from_samplesheet_line(line)
            if not indx:
                indx = "NoIndex"
            if sample_sheet.has_lanes:
                lane = line['Lane']
            else:
                lane = None
            sample.add_barcode(indx,lane=lane)
            # Sample index
            if sample.s_index is None:
                s_index += 1
                sample.s_index = s_index
    def get_project(self,project_name):
        # Fetch project by name
        for project in self.projects:
            if project.name == project_name:
                return project
        raise KeyError("%s: not found" % project_name)
    def add_project(self,project_name):
        # Add a new project or return existing one
        try:
            return self.get_project(project_name)
        except KeyError:
            project = SampleSheetProject(project_name)
            self.projects.append(project)
            return project
    @property
    def project_names(self):
        # Return sorted list of project names
        names = [str(p) for p in self.projects]
        names.sort()
        return names
    @property
    def nprojects(self):
        # Return number of projects
        return len(self.projects)

class SampleSheetProject(object):
    """
    """
    def __init__(self,project_name):
        self.name = project_name
        self.samples = []
    def get_sample(self,sample_id):
        # Fetch sample by name
        for sample in self.samples:
            if sample.sample_id == sample_id:
                return sample
        raise KeyError("%s: not found" % sample_id)
    def add_sample(self,sample_id,sample_name=None,s_index=None):
        # Add a new sample or return existing one
        try:
            return self.get_sample(sample_id)
        except KeyError:
            sample = SampleSheetSample(sample_id,
                                       sample_name=sample_name)
            self.samples.append(sample)
            return sample
    @property
    def sample_ids(self):
        # Return sorted list of sample ids
        ids = [s.sample_id for s in self.samples]
        ids.sort()
        return ids
    def __repr__(self):
        return str(self.name)

class SampleSheetSample(object):
    """
    """
    def __init__(self,sample_id,sample_name=None,s_index=None):
        self.sample_id = sample_id
        self.sample_name = sample_name
        self.s_index = s_index
        self.barcodes = {}
    def add_barcode(self,barcode,lane=None):
        # Add a new barcode or return existing one
        if barcode not in self.barcodes:
            self.barcodes[barcode] = []
        if lane and lane not in self.barcodes[barcode]:
            self.barcodes[barcode].append(lane)
        return self.barcodes[barcode]
    def __repr__(self):
        if not self.sample_id:
            return str(self.sample_name)
        if self.sample_name and self.sample_name != self.sample_id:
            return "%s/%s" % (self.sample_name,self.sample_id)
        return str(self.sample_id)

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

def main():
    p = optparse.OptionParser()
    opts,args = p.parse_args()
    # Get prediction
    predictor = SampleSheetPredictor(args[0])
    # Check for close-matching project names
    close_names = get_close_names(predictor.project_names)
    # Report project names
    print
    print "Expecting %d project%s:" % (predictor.nprojects,
                                       '' if predictor.nprojects == 1
                                       else 's')
    print
    warning_close_names = False
    for p in predictor.project_names:
        warning_flag = False
        if p in close_names:
            warning_flag = True
            warning_close_names = True
        print "%s %s" % (('*' if warning_flag else '-'),p)
    print
    # Report samples within each project
    warning_wrong_barcodes = False
    for pname in predictor.project_names:
        project = predictor.get_project(pname)
        print "%s\n%s" % (pname,'-'*len(pname))
        for sname in project.sample_ids:
            sample = project.get_sample(sname)
            warning_wrong_barcodes = warning_wrong_barcodes or \
                                     (len(sample.barcodes) > 1)
            warning_flag = (len(sample.barcodes) > 1)
            for barcode in sample.barcodes:
                #print samples[s][barcode]
                lanes = ','.join([str(x)
                                  for x in sample.barcodes[barcode]])
                print "%s%s\tS%s\t%s\t%s" % (('* ' if warning_flag else ''),
                                             sample,
                                             sample.s_index,
                                             barcode,
                                             (('L' + lanes) if lanes else ''))
        print
    # Issue warnings
    if warning_close_names:
        print "*** WARNING projects detected with closing-matching names ***"
    if warning_wrong_barcodes:
        print "*** WARNING samples detected with multiple barcodes ***"

if __name__ == "__main__":
    main()
