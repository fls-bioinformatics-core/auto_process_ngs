#!/bin/env python
#
#     auto_process.py: automated processing of Illumina sequence data
#     Copyright (C) University of Manchester 2013 Peter Briggs
#
#########################################################################
#
# auto_process.py
#
#########################################################################

"""
First attempt at an automated data processing & QC pipeline in Python

Implements a program for automating stages of a standard protocol for
processing and QC'ing Illumina sequencing data.

The stages are:

    setup
    make_fastqs
    setup_analysis_dirs
    run_qc

The 'setup' stage creates an analysis directory and acquires the basic
data about the sequencing run from a source directory. Subsequent stages
should be run in sequence to create fastq files, set up analysis
directories for each project, and run QC scripts for each sample in
each project.

"""

__version__ = "0.0.0"

#######################################################################
# Imports
#######################################################################

import sys
import os
import re
import subprocess
import logging
import optparse
import shutil
import time
import IlluminaData
import platforms
import TabFile
import FASTQFile
import JobRunner
import Pipeline
import bcf_utils
import analyse_illumina_run
import build_illumina_analysis_dir
import qcreporter
import bclToFastq
import applications
import auto_process_utils
# Import local settings
try:
    import auto_process_settings
except ImportError:
    logging.warning("Failed to import auto_process_settings.py")

#######################################################################
# Classes
#######################################################################

class AnalysisProject:
    # Class describing an analysis project
    # i.e. set of samples from a single sequencing experiment
    # THINKS can code be shared with IlluminaProject?
    def __init__(self,name,dirn,library_type=None,organism=None):
        # name = project name
        # dirn = directory path
        self.name = name
        self.dirn = dirn
        self.library_type = library_type
        self.organism = organism
        self.samples = []
        self.paired_end = False
        # Check for metadata file
        ## NOT IMPLEMENTED ##
        # Otherwise populate from files, using brute force
        # approach based on names
        print "Acquiring fastqs..."
        fastqs = Pipeline.GetFastqGzFiles(self.dirn)
        for fq in fastqs:
            print "%s" % fq
        print "Assigning fastqs to samples..."
        for fastq in fastqs:
            # GetFastqGzFile returns a list of tuples
            for fq in fastq:
                name = auto_process_utils.AnalysisFastq(fq).sample_name
                try:
                    sample = self.get_sample(name)
                except KeyError:
                    sample = AnalysisSample(name)
                    self.samples.append(sample)
                sample.add_fastq(fq)
        print "Listing samples and files:"
        for sample in self.samples:
            print "%s: %s" % (sample.name,sample.fastq)

    def get_sample(self,name):
        """Return sample that matches 'name'

        Arguments:
          name: name of a sample

        Returns:
          AnalysisSample object with the matching name; raises
          KeyError exception if no match is found.

        """
        for sample in self.samples:
            if sample.name == name: return sample
        raise KeyError, "No matching sample for '%s'" % name

    def prettyPrintSamples(self):
        """Return a nicely formatted string describing the sample names

        Wraps a call to 'pretty_print_names' function.
        """
        return bcf_utils.pretty_print_names(self.samples)

class AnalysisSample:
    # Class describing an analysis sample
    # i.e. set of fastqs from a single sample
    # Can be paired end and have multiple fastqs per
    # THINKS can code be shared with IlluminaSample?
    def __init__(self,name):
        self.name = name
        self.fastq = []
        self.paired_end = False

    def add_fastq(self,fastq):
        """Add a reference to a fastq file in the sample

        Arguments:
          fastq: name of the fastq file
        """
        self.fastq.append(fastq)
        # Sort fastq's into order
        self.fastq.sort()
        # Check paired-end status
        if not self.paired_end:
            fq = auto_process_utils.AnalysisFastq(fastq)
            if fq.read_number == 2:
                self.paired_end = True

    def fastq_subset(self,read_number=None,full_path=False):
        """Return a subset of fastq files from the sample

        Arguments:
          read_number: select subset based on read_number (1 or 2)
          full_path  : if True then fastq files will be returned
            with the full path, if False (default) then as file
            names only.

        Returns:
          List of fastq files matching the selection criteria.

        """
        # Build list of fastqs that match the selection criteria
        fastqs = []
        for fastq in self.fastq:
            fq = auto_process_utils.AnalysisFastq(fastq)
            if fq.read_number is None:
                raise IlluminaDataException, \
                    "Unable to determine read number for %s" % fastq
            if fq.read_number == read_number:
                if full_path:
                    fastqs.append(os.path.join(self.dirn,fastq))
                else:
                    fastqs.append(fastq)
        # Sort into dictionary order and return
        fastqs.sort()
        return fastqs

    def __repr__(self):
        """Implement __repr__ built-in

        Return string representation for the sample -
        i.e. the sample name."""
        return str(self.name)

class AutoProcess:
    # Class implementing an automatic fastq generation and QC
    # processing procedure

    def __init__(self,analysis_dir=None):
        # analysis_dir: name/path for existing analysis directory
        self.params = auto_process_utils.AttributeDict()
        self.analysis_dir = analysis_dir
        if self.analysis_dir is not None:
            self.analysis_dir = os.path.abspath(analysis_dir)
            self.load_parameters()
            self.params['analysis_dir'] = self.analysis_dir

    def add_directory(self,sub_dir):
        # Add a directory to the AutoProcess object
        dirn = os.path.join(self.analysis_dir,sub_dir)
        self.create_directory(dirn)
        return dirn

    def create_directory(self,dirn):
        # Make the specified directory, and any leading directories
        # that don't already exist
        if not os.path.exists(dirn):
            dir_path = os.sep
            for sub_dir in dirn.split(os.sep):
                dir_path = os.path.join(dir_path,sub_dir)
                if not os.path.exists(dir_path):
                    print "Making %s" % dir_path
                    bcf_utils.mkdir(dir_path)

    def load_parameters(self):
        # Get parameters from info file
        #
        # check for info file
        info_file_name = os.path.join(self.analysis_dir,'auto_process.info')
        if not os.path.isfile(info_file_name):
            raise Exception, "No info file %s" % info_file_name
        # Read contents of info file and assign values
        print "Loading settings from %s" % info_file_name
        info_file = TabFile.TabFile(info_file_name)
        for line in info_file:
            key = line[0]
            value = line[1]
            self.params[key] = value
            print "%s\t%s" % (key,self.params[key])

    def save_parameters(self):
        # Save parameters to info file
        #
        info_file_name = os.path.join(self.analysis_dir,'auto_process.info')
        info_file = TabFile.TabFile()
        for param in self.params:
            info_file.append(data=(param,self.params[param]))
        info_file.write(os.path.join(self.analysis_dir,'auto_process.info'))

    def log_path(self,*args):
        # Return path appended to log directory
        # Use for getting paths of files under the logs directory
        return os.path.join(self.log_dir,*args)

    def __del__(self):
        print "Saving parameters to file"
        self.save_parameters()

    @property
    def log_dir(self):
        # Generate and return full path to log directory
        return self.add_directory('logs')

    @property
    def tmp_dir(self):
        # Generate and return full path to tmp directory
        return self.add_directory('tmp')

    def setup(self,data_dir,analysis_dir=None):
        # Set up the initial analysis directory
        #
        # This does all the initialisation of the analysis directory
        # and processing parameters
        #
        # Arguments:
        # data_dir: source data directory
        # analysis_dir: corresponding analysis dir
        data_dir = data_dir.rstrip(os.sep)
        self.analysis_dir = analysis_dir
        if self.analysis_dir is None:
            self.analysis_dir = os.path.join(os.getcwd(),
                                             os.path.basename(data_dir))+'_analysis'
        # Create the analysis directory structure
        if not os.path.exists(self.analysis_dir):
            # Create directory structure
            self.create_directory(self.analysis_dir)
            # Identify platform
            platform = platforms.get_sequencer_platform(data_dir)
            print "Platform identified as '%s'" % platform
            # Fetch SampleSheet.csv file
            print "Acquiring sample sheet..."
            tmp_sample_sheet = os.path.join(self.tmp_dir,'SampleSheet.csv')
            rsync = applications.general.rsync(os.path.join(data_dir,
                                                            'Data/Intensities/Basecalls/SampleSheet.csv'),
                                               self.tmp_dir)
            status = rsync.run_subprocess(log=self.log_path('rsync.sample_sheet.log'))
            custom_sample_sheet = os.path.join(self.analysis_dir,'custom_SampleSheet.csv')
            sample_sheet = make_custom_sample_sheet(tmp_sample_sheet,
                                                    custom_sample_sheet)
            os.remove(tmp_sample_sheet)
            # Fetch RunInfo.xml file and get bases mask
            print "Acquiring RunInfo.xml..."
            tmp_run_info = os.path.join(self.tmp_dir,'RunInfo.xml')
            rsync = applications.general.rsync(os.path.join(data_dir,'RunInfo.xml'),
                                               self.tmp_dir)
            status = rsync.run_subprocess(log=self.log_path('rsync.run_info.log'))
            bases_mask = get_bases_mask(tmp_run_info,custom_sample_sheet)
            print "Corrected bases mask: %s" % bases_mask
            os.remove(tmp_run_info)
            # Print the predicted ouputs
            projects = sample_sheet.predict_output()
            print "Predicted output from sample sheet:"
            print "Project\tSample\tFastq"
            for project in projects:
                project_name = project[8:]
                for sample in projects[project]:
                    sample_name = sample[7:]
                    for fastq_base in projects[project][sample]:
                        print "%s\t%s\t%s" % (project_name,sample_name,fastq_base)
            # Store the parameters
            self.params['data_dir'] = data_dir
            self.params['analysis_dir'] = self.analysis_dir
            self.params['platform'] = platform
            self.params['sample_sheet'] = custom_sample_sheet
            self.params['bases_mask'] = bases_mask
        else:
            # Directory already exists
            logging.warning("Analysis directory already exists")
            # check for info file
            info_file_name = os.path.join(self.analysis_dir,'auto_process.info')
            if os.path.isfile(info_file_name):
                self.load_parameters()
            else:
                logging.error("No info file")

    def get_primary_data(self):
        # Copy the primary sequencing data (bcl files etc) to a local area
        # using rsync
        data_dir = self.params.data_dir
        self.params["primary_data_dir"] = self.add_directory('primary_data')
        try:
            rsync = applications.general.rsync(data_dir,self.params.primary_data_dir,
                                               prune_empty_dirs=True,
                                               extra_options=('--include=*/',
                                                              '--include=Data/**',
                                                              '--include=RunInfo.xml',
                                                              '--include=SampleSheet.csv',
                                                              '--exclude=*'))
            print "Running %s" % rsync
            status = rsync.run_subprocess(log=self.log_path('rsync.primary_data.log'))
        except Exception, ex:
            logging.error("Exception getting primary data: %s" % ex)
            status = -1
        if status != 0:
            logging.error("Failed to acquire primary data (status %s)" % status)
        return status
        
    def bcl_to_fastq(self):
        # Convert bcl files to fastq
        #
        # Directories
        analysis_dir = self.params.analysis_dir
        primary_data = os.path.join(self.params.primary_data_dir,
                                    os.path.basename(self.params.data_dir))
        self.params['unaligned_dir'] = 'bcl2fastq'
        bcl2fastq_dir = self.add_directory(self.params.unaligned_dir)
        # Get info about the run
        illumina_run = IlluminaData.IlluminaRun(primary_data)
        print "%s" % illumina_run.run_dir
        print "Platform  : %s" % illumina_run.platform
        print "Bcl format: %s" % illumina_run.bcl_extension
        # Get sample sheet and bases mask
        sample_sheet = self.params.sample_sheet
        bases_mask = self.params.bases_mask
        print "Sample sheet: %s" % sample_sheet
        print "Bases mask  : %s" % bases_mask
        # Get nmismatches
        nmismatches = bclToFastq.get_nmismatches(bases_mask)
        print "Nmismatches : %d (determined from bases mask)" % nmismatches
        # Set up runner
        runner = auto_process_settings.runners.bcl2fastq
        runner.log_dir(self.log_dir)
        # Run bcl2fastq
        bcl2fastq = applications.Command('bclToFastq.py',
                                         '--use-bases-mask',bases_mask,
                                         '--nmismatches',nmismatches,
                                         primary_data,
                                         bcl2fastq_dir,
                                         sample_sheet)
        print "Running %s" % bcl2fastq
        bcl2fastq_job = Pipeline.Job(runner,
                                     'bclToFastq',
                                     os.getcwd(),
                                     bcl2fastq.command,
                                     bcl2fastq.args)
        bcl2fastq_job.start()
        while bcl2fastq_job.isRunning():
            time.sleep(10)
        print "bcl2fastq completed"
        # Verify outputs
        illumina_data = IlluminaData.IlluminaData(self.analysis_dir,
                                                  unaligned_dir=self.params.unaligned_dir)
        if not analyse_illumina_run.verify_run_against_sample_sheet(illumina_data,
                                                                    sample_sheet):
            logging.error("Failed to verify bcl to fastq outputs against sample sheet")
            raise Exception, "Failed to verify bcl to fastq outputs against sample sheet"
        # Create a metadata file
        self.params['project_metadata'] = 'projects.info'
        metadata = get_project_metadata(illumina_data)
        metadata.write(os.path.join(analysis_dir,self.params.project_metadata),
                       include_header=True)
        print "Project metadata in %s" % self.params.project_metadata
        # Generate statistics
        self.params['stats_file'] = 'statistics.info'
        stats = illumina_data_statistics(illumina_data)
        stats.write(os.path.join(analysis_dir,self.params.stats_file),
                    include_header=True)
        print "Statistics in %s" % self.params.stats_file

    def remove_primary_data(self):
        # Remove primary data
        primary_data = os.path.join(self.params.primary_data_dir,
                                    os.path.basename(self.params.data_dir))
        if os.path.isdir(primary_data):
            print "Removing copy of primary data in %s" % primary_data
            shutil.rmtree(primary_data)

    def setup_analysis_dirs(self):
        # Construct and populate the analysis directories for each project
        illumina_data = IlluminaData.IlluminaData(self.analysis_dir,
                                                  unaligned_dir=self.params.unaligned_dir)
        project_metadata = TabFile.TabFile(os.path.join(analysis_dir,
                                                        self.params.project_metadata),
                                           first_line_is_header=True)
        for line in project_metadata:
            # Look up project data
            project_name = line['Project']
            project = illumina_data.get_project(project_name)
            # Make the project directory
            build_illumina_analysis_dir.create_analysis_dir(project,
                                                            top_dir=self.analysis_dir)
            # Add a logs subdirectory for each project
            log_dir = os.path.join(self.analysis_dir,project_name,'logs')
            bcf_utils.mkdir(log_dir)

    def run_qc(self):
        # Run QC pipeline for all projects
        #
        # Setup a pipeline runner
        qc_runner = auto_process_settings.runners.qc
        self.pipeline = Pipeline.PipelineRunner(qc_runner)
        # Get project data
        project_metadata = TabFile.TabFile(os.path.join(analysis_dir,
                                                        self.params.project_metadata),
                                           first_line_is_header=True)
        # Make list of project dirs
        projects = []
        for line in project_metadata:
            name = line['Project']
            print "Acquiring data for project %s" % name
            projects.append(AnalysisProject(name,
                                            os.path.join(self.analysis_dir,name)))
        # Populate pipeline with fastq.gz files
        for project in projects:
            print "*** Setting up QC for %s ***" % project.name
            qc_log_dir = os.path.join(self.analysis_dir,project.name,'logs')
            qc_runner.log_dir(qc_log_dir)
            for sample in project.samples:
                print "\tSetting up sample %s" % sample.name
                group = "%s.%s" % (project.name,sample.name)
                for fq in sample.fastq:
                    print "\t\tAdding %s" % fq
                    fastq = os.path.join(project.dirn,fq)
                    label = str(auto_process_utils.AnalysisFastq(fq))
                    self.pipeline.queueJob(project.dirn,'illumina_qc.sh',(fastq,),
                                           label=label,group=group)
        # Run the pipeline
        self.pipeline.run()
        # Verify the outputs
        for project in projects:
            qc_reporter = qcreporter.IlluminaQCReporter(project.dirn)
            if not qc_reporter.verify():
                logging.error("QC failed for one or more samples in %s" % project.name)
            else:
                print "QC okay, generating report for %s" % project.name
                qc_reporter.zip()

    def merge_fastqs(self):
        # Merge multiple fastq files for samples split across lanes
        raise NotImplementedError

    def copy_to_archive(self):
        # Copy the analysis directory and contents to an archive area
        raise NotImplementedError

#######################################################################
# Functions
#######################################################################

def get_project_metadata(illumina_data):
    # Make a TSV file with information on projects
    metadata = TabFile.TabFile(column_names=('Project',
                                             'Samples',
                                             'Paired_end'))
    for project in illumina_data.projects:
        metadata.append(data=(project.name,
                              project.prettyPrintSamples(),
                              'Y' if project.paired_end else 'N'))
    return metadata

def illumina_data_statistics(illumina_data):
    # Gather statistics for an Illumina run (bcl2fastq outputs)
    # Set up TabFile to hold the data collected
    stats = TabFile.TabFile(column_names=('Project',
                                          'Sample',
                                          'Fastq',
                                          'Size',
                                          'Nreads',
                                          'Paired_end'))
    # Get number of reads and file size for each file across all
    # projects and samples
    for project in illumina_data.projects:
        for sample in project.samples:
            for fastq in sample.fastq:
                fq = os.path.join(sample.dirn,fastq)
                nreads = FASTQFile.nreads(fq)
                fsize = os.path.getsize(fq)
                stats.append(data=(project.name,
                                   sample,fastq,
                                   bcf_utils.format_file_size(fsize),
                                   nreads,
                                   'Y' if sample.paired_end else 'N'))
    # Gather same information for undetermined reads (if present)
    if illumina_data.undetermined is not None:
        for sample in illumina_data.undetermined.samples:
            for fastq in sample.fastq:
                fq = os.path.join(sample.dirn,fastq)
                nreads = FASTQFile.nreads(fq)
                fsize = os.path.getsize(fq)
                stats.append(data=(illumina_data.undetermined.name,
                                   sample,fastq,
                                   bcf_utils.format_file_size(fsize),
                                   nreads,
                                   'Y' if sample.paired_end else 'N'))
    # Return the data
    return stats

def make_custom_sample_sheet(input_sample_sheet,output_sample_sheet=None):
    # Read sample sheet info from input_sample_sheet
    # Do clean up
    # Write to output_sample_sheet (if specified)
    # Return CasavaSampleSheet object
    sample_sheet = IlluminaData.get_casava_sample_sheet(input_sample_sheet)
    for line in sample_sheet:
        if not line['SampleProject']:
            line['SampleProject'] = line['SampleID']
    sample_sheet.fix_illegal_names()
    sample_sheet.fix_duplicated_names()
    if output_sample_sheet is not None:
        sample_sheet.write(output_sample_sheet)
    return sample_sheet

def get_bases_mask(run_info_xml,sample_sheet_file):
    # Return bases mask string generated from data in RunInfo.xml and
    # sample sheet files
    # Get initial bases mask
    bases_mask = IlluminaData.IlluminaRunInfo(run_info_xml).bases_mask
    print "Bases mask: %s (from RunInfo.xml)" % bases_mask
    # Update bases mask from sample sheet
    example_barcode = IlluminaData.get_casava_sample_sheet(sample_sheet_file)[0]['Index']
    bases_mask = IlluminaData.fix_bases_mask(bases_mask,example_barcode)
    print "Bases mask: %s (updated for barcode sequence '%s')" % (bases_mask,
                                                                  example_barcode)
    return bases_mask

#######################################################################
# Main program
#######################################################################

if __name__ == "__main__":

    # Set up option parser
    p = optparse.OptionParser(usage="%prog CMD [OPTIONS] [ARGS]",
                              version="%prog "+__version__,
                              description="Automatically process Illumina sequence data")
    p.add_option('--debug',action='store_true',dest='debug',default=False,
                 help="Turn on debugging output from Python libraries")
    if len(sys.argv) < 2:
        p.error("Takes at least one argument (command)")
        sys.exit()

    # Process command line
    cmd = sys.argv[1]
    if cmd not in ['setup','make_fastqs','setup_analysis_dirs','run_qc']:
        p.error("Unrecognised command '%s'" % cmd)
    options,args = p.parse_args(sys.argv[2:])

    # Debug output?
    if options.debug:
        logging.getLogger().setLevel(logging.DEBUG)

    # Report name and version
    print "%s %s" % (os.path.basename(sys.argv[0]),__version__)

    # Setup the processing object and run the requested command
    if cmd == 'setup':
        if len(args) != 1:
            p.error("Need to supply a source data location")
        d = AutoProcess()
        d.setup(args[0])
    else:
        # For other options check if an analysis
        # directory was specified
        if len(args) > 0:
            analysis_dir = args[0]
        else:
            analysis_dir = os.getcwd()
        d = AutoProcess(analysis_dir)
        # Run the specified stage
        if cmd == 'make_fastqs':
            d.get_primary_data()
            d.bcl_to_fastq()
            d.remove_primary_data()
        elif cmd == 'setup_analysis_dirs':
            d.setup_analysis_dirs()
        elif cmd == 'run_qc':
            d.run_qc()
