#!/usr/bin/env python
#
#     icell8_qc.py: run QC pipeline on iCell8 fastq data
#     Copyright (C) University of Manchester 2016 Peter Briggs
#
#########################################################################
#
# icell8_qc.py
#
#########################################################################

"""
Automated QC pipeline for iCell8 sequence data

The data is expected to be in a top-level directory containing:

 * Directory of Fastq files
 * Metadata file

The FASTQ files are named as e.g.

 * icell8_manchester_demo_fastqs.TGGCCAGCAGT.r1.fastq
 * icell8_manchester_demo_fastqs.TGGCCAGCAGT.r2.fastq

i.e. BASENAME.INDEX_SEQ.r[1|2].fastq

The metadata file consists of columns of information which map
onto the index sequences in the metadata file.

There are two stages:

 1. Create symlinks to FASTQ files in a new directory with names
    based on metadata
 2. Run the QC pipeline on the symlinks

"""

#######################################################################
# Imports
#######################################################################

import os
import optparse
from bcftbx.utils import AttributeDictionary
from bcftbx.utils import mklink
from bcftbx.utils import mkdir
from auto_process_ngs.utils import AnalysisProject
from auto_process_ngs.applications import Command
from auto_process_ngs.simple_scheduler import SimpleScheduler
from auto_process_ngs import envmod

import logging
logging.basicConfig(format='%(levelname) 8s: %(message)s')

import auto_process_ngs.settings
__settings = auto_process_ngs.settings.Settings()

try:
    __modulefiles = __settings.modulefiles['run_qc']
except KeyError:
    # No environment modules specified
    __modulefiles = None

#######################################################################
# Classes
#######################################################################

class ICell8Metadata(object):
    """
    Class representing an icell8 metadata file

    The metadata file consists of comment lines (preceeded by
    a hash) and tab-delimited columns:

    # Column 1  The Read 1 barcode sequence.
    # Column 2  The Read 2 barcode sequence (in this case, a dummy value)
    # Column 3  The sample type, e.g. Cell, RNA, Tube, etc.
    # Column 4  The source sample name.
    # Column 5  The barcode's dispense order, if applicable.
    # Column 6  The 384-well source plate location, if applicable.
    #               Range = A1 - H12
    # Column 7  The chip row ID, if applicable.  Range = 0 - 71
    # Column 8  The chip column ID, if applicable.  Range = 0 - 71
    # Column 9  The image ID in which the well is visualized, if applicable.

    """
    def __init__(self,metadata_file):
        """
        """
        self._attrs = ('barcode_r1',
                       'barcode_r2',
                       'sample_type',
                       'source_sample_name',
                       'barcode_dispense_order',
                       'source_plate_location',
                       'chip_row_id',
                       'chip_col_id',
                       'image_id')
        self._metadata = []
        # Read in the data
        with open(metadata_file,'r') as fp:
            for line in fp:
                if line.startswith('#'):
                    continue
                cols = line.rstrip('\n').split('\t')
                self.add_metadata(*cols)

    def add_metadata(self,barcode_r1,barcode_r2,sample_type,
                     source_sample_name,barcode_dispense_order,
                     source_plate_location,chip_row_id,chip_col_id,
                     image_id):
        """
        """
        entry = AttributeDictionary()
        entry['barcode_r1'] = barcode_r1
        entry['barcode_r2'] = barcode_r2
        entry['sample_type'] = sample_type
        entry['source_sample_name'] = source_sample_name
        entry['barcode_dispense_order'] = barcode_dispense_order
        entry['source_plate_location'] = source_plate_location
        entry['chip_row_id'] = chip_row_id
        entry['chip_col_id'] = chip_col_id
        entry['image_id'] = image_id
        for attr in self._attrs:
            if entry[attr] == 'N/A':
                entry[attr] = None
        self._metadata.append(entry)

    def lookup(self,**kws):
        """
        Look up entries based on one or more matching attributes
        """
        results = self._metadata
        for kw in kws:
            print "Filtering on %s..." % kw
            results=filter(lambda x: x[kw] == kws[kw],results)
        return results

#######################################################################
# Functions
#######################################################################

def announce(title):
    """
    """
    title = str(title)
    len_title = len(title)
    print "="*len_title
    print title
    print "="*len_title

#######################################################################
# Main program
#######################################################################

if __name__ == "__main__":
    # Make a command line parser
    p = optparse.OptionParser()

    # Parse the command line
    opt,args = p.parse_args()

    # Get the metadata file name
    try:
        metadata = ICell8Metadata(args[0])
    except IndexError:
        p.error("Need to supply a metadata file")
    
    # Directory with FASTQ files
    try:
        fastqs_dir = os.path.abspath(args[1])
    except IndexError:
        p.error("Need to supply a FASTQ directory")

    # Target directory for QC and analysis
    try:
        project_dir = os.path.abspath(args[2])
    except IndexError:
        p.error("Need to supply a target directory")

    # Get list of FASTQs and map to new names
    announce("Obtaining list of FASTQs")
    fastqs = sorted(os.listdir(fastqs_dir))
    mapping = {}
    for fq in fastqs:
        print fq
        name,barcode,read = fq.split('.')[0:3]
        print "-- %s" % name
        print "-- %s" % barcode
        print "-- %s" % read
        entries = metadata.lookup(barcode_r1=barcode)
        if not entries:
            raise KeyError("Can't find barcode '%s' in metadata" %
                           barcode)
        if len(entries) > 1:
            raise KeyError("Barcode '%s' matches multiple entries" %
                           barcode)
        print entries[0].sample_type    
        print entries[0].chip_row_id
        print entries[0].chip_col_id
        sample_type = entries[0].sample_type
        chip_row_id = entries[0].chip_row_id
        chip_col_id = entries[0].chip_col_id
        mapping[fq] = "%s_%s_%s.%s.fastq" % (sample_type,
                                             chip_col_id,
                                             chip_row_id,
                                             read)

    # Create symlinks with new names
    announce("Creating project dir with symlinks")
    if not os.path.isdir(project_dir):
        mkdir(project_dir)
        for fq in fastqs:
            print "%s -> %s" % (fq,mapping[fq])
            mklink(os.path.join(fastqs_dir,fq),
                   os.path.join(project_dir,mapping[fq]),
                   relative=True)
    else:
        logging.warning("'%s' already exists" % project_dir)

    # Set up environment
    if __modulefiles is not None:
        announce("Setting up environment")
        for modulefile in __modulefiles.split(','):
            envmod.load(modulefile)

    # Run the QC
    announce("Running QC")
    qc_runner = __settings.runners.qc
    max_jobs = __settings.general.max_concurrent_jobs
    sched = SimpleScheduler(runner=qc_runner,
                            max_concurrent=max_jobs)
    sched.start()
    project = AnalysisProject("test1",project_dir)
    qc_dir = os.path.join(project_dir,'qc')
    log_dir = os.path.join(project_dir,'qc','logs')
    if project.qc_dir is None:
        print "Making 'qc' subdirectory..."
        mkdir(qc_dir)
        mkdir(log_dir)
    if project.verify_qc():
        print "QC verified for these files"
    else:
        for sample in project.samples:
            print "Checking/setting up for sample '%s'" % sample.name
            for fq in sample.fastq:
                if sample.verify_qc(qc_dir,fq):
                    print "-- %s: QC ok"
                else:
                    print "-- %s: setting up QC" % fq
                    qc_cmd = Command('illumina_qc.sh',fq)
                    job = sched.submit(qc_cmd,
                                       wd=project_dir,
                                       name="qc.%s" % os.path.basename(fq),
                                       log_dir=log_dir)
                    print "Job: %s" % job
    # Wait for the scheduler to run all jobs
    sched.wait()
    sched.stop()

    # Verify the QC
    if not project.verify_qc():
        logging.error("QC failed")
    else:
        print "QC ok: generating report"
        project.qc_report()


