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
    archive
    publish_qc

The 'setup' stage creates an analysis directory and acquires the basic
data about the sequencing run from a source directory. Subsequent stages
should be run in sequence to create fastq files, set up analysis
directories for each project, and run QC scripts for each sample in
each project.

"""

__version__ = "0.0.14"

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
import qcreporter
import bclToFastq
import applications
import auto_process_utils
import auto_process_settings

#######################################################################
# Classes
#######################################################################

class AutoProcess:
    # Class implementing an automatic fastq generation and QC
    # processing procedure

    def __init__(self,analysis_dir=None):
        # analysis_dir: name/path for existing analysis directory
        self.params = auto_process_utils.AnalysisDirMetadata()
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
        self.params.load(info_file_name)

    def save_parameters(self):
        # Save parameters to info file
        #
        info_file_name = os.path.join(self.analysis_dir,'auto_process.info')
        self.params.save(info_file_name)

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
        else:
            # Directory already exists
            logging.warning("Analysis directory already exists")
            # check for info file
            info_file_name = os.path.join(self.analysis_dir,'auto_process.info')
            if os.path.isfile(info_file_name):
                self.load_parameters()
            else:
                logging.warning("No info file found in %s" % self.analysis_dir)
        # Identify missing data and attempt to acquire
        # Sequencing platform
        platform = self.params.platform
        if platform is None:
            print "Identifying platform from data directory name"
            platform = platforms.get_sequencer_platform(data_dir)
        print "Platform identified as '%s'" % platform
        # Custom SampleSheet.csv file
        custom_sample_sheet = self.params.sample_sheet
        if custom_sample_sheet is None:
            print "Acquiring sample sheet..."
            tmp_sample_sheet = os.path.join(self.tmp_dir,'SampleSheet.csv')
            rsync = applications.general.rsync(os.path.join(data_dir,
                                                            'Data/Intensities/BaseCalls/SampleSheet.csv'),
                                               self.tmp_dir)
            print "%s" % rsync
            status = rsync.run_subprocess(log=self.log_path('rsync.sample_sheet.log'))
            custom_sample_sheet = os.path.join(self.analysis_dir,'custom_SampleSheet.csv')
            sample_sheet = make_custom_sample_sheet(tmp_sample_sheet,
                                                    custom_sample_sheet)
            print "Keeping copy of original sample sheet"
            os.rename(tmp_sample_sheet,os.path.join(self.analysis_dir,'SampleSheet.orig.csv'))
        sample_sheet = IlluminaData.CasavaSampleSheet(custom_sample_sheet)
        print "Sample sheet '%s'" % custom_sample_sheet
        # Bases mask
        bases_mask = self.params.bases_mask
        if bases_mask is None:
            print "Acquiring RunInfo.xml to determine bases mask..."
            tmp_run_info = os.path.join(self.tmp_dir,'RunInfo.xml')
            rsync = applications.general.rsync(os.path.join(data_dir,'RunInfo.xml'),
                                               self.tmp_dir)
            status = rsync.run_subprocess(log=self.log_path('rsync.run_info.log'))
            bases_mask = get_bases_mask(tmp_run_info,custom_sample_sheet)
            os.remove(tmp_run_info)
        print "Corrected bases mask: %s" % bases_mask
        # Print the predicted ouputs
        projects = sample_sheet.predict_output()
        print "Predicted output from sample sheet:"
        print "Project\tSample\tFastq"
        for project in projects:
            project_name = project[8:]
            sample_names = []
            for sample in projects[project]:
                sample_name = sample[7:]
                for fastq_base in projects[project][sample]:
                    print "%s\t%s\t%s" % (project_name,sample_name,fastq_base)
                sample_names.append(sample_name)
        # Store the parameters
        self.params['data_dir'] = data_dir
        self.params['analysis_dir'] = self.analysis_dir
        self.params['platform'] = platform
        self.params['sample_sheet'] = custom_sample_sheet
        self.params['bases_mask'] = bases_mask

    def setup_from_fastq_dir(self,analysis_dir,fastq_dir):
        # Do setup for an existing directory containing fastq files
        # with the same structure as that produced by CASAVA and bcl2fastq
        #
        # Assumes that the files are in a subdirectory of the analysis
        # directory specified by the 'fastqe_dir' argument, and
        # that within that they are arranged in the structure
        # 'Project_<name>/Sample_<name>/<fastq>'
        self.analysis_dir = os.path.abspath(analysis_dir)
        # Create directory structure
        self.create_directory(self.analysis_dir)
        # Get information
        basespace = IlluminaData.IlluminaData(self.analysis_dir,
                                              unaligned_dir=fastq_dir)
        print "Identifying platform from data directory name"
        platform = platforms.get_sequencer_platform(analysis_dir)
        # Store the parameters
        self.params['analysis_dir'] = self.analysis_dir
        self.params['unaligned_dir'] = fastq_dir
        self.params['platform'] = platform
        # Generate statistics
        self.generate_stats()
        # Make a 'projects.info' metadata file
        self.make_project_metadata_file()

    def make_project_metadata_file(self,project_metadata_file='projects.info'):
        # Generate a project metadata file based on the fastq
        # files and directory structure
        illumina_data = IlluminaData.IlluminaData(self.analysis_dir,
                                                  unaligned_dir=self.params.unaligned_dir)
        project_metadata = ProjectMetadataFile()
        print "Project\tSample\tFastq"
        for project in illumina_data.projects:
            project_name = project.name
            sample_names = []
            for sample in project.samples:
                sample_name = sample.name
                for fastq in sample.fastq:
                    print "%s\t%s\t%s" % (project_name,sample_name,fastq)
                sample_names.append(sample_name)
            project_metadata.add_project(project_name,sample_names)
        project_metadata.save(os.path.join(self.params.analysis_dir,project_metadata_file))
        self.params['project_metadata'] = project_metadata_file
        print "Project metadata in %s" % self.params.project_metadata

    def get_analysis_projects(self,pattern=None):
        # Return the analysis projects in a list
        #
        # By default returns all projects within the analysis
        #
        # If the 'pattern' is not None then it should be a simple pattern
        # used to match against available names to select a subset of
        # projects (see bcf_utils.name_matches).
        project_metadata = ProjectMetadataFile(os.path.join(self.analysis_dir,
                                                            self.params.project_metadata))
        projects = []
        for line in project_metadata:
            name = line['Project']
            if pattern is not None:
                if not bcf_utils.name_matches(name,pattern):
                    # Name failed to match, ignore
                    continue
            print "Acquiring data for project %s" % name
            project_dir = os.path.join(self.analysis_dir,name)
            if os.path.isdir(project_dir):
                projects.append(
                    auto_process_utils.AnalysisProject(
                        name,
                        os.path.join(self.analysis_dir,name)))
            else:
                logging.warning("No directory found for project '%s'" % name)
        return projects

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
        
    def bcl_to_fastq(self,ignore_missing_bcl=False,ignore_missing_stats=False,
                     skip_rsync=False,keep_primary_data=False):
        # Convert bcl files to fastq
        #
        # Arguments:
        # ignore_missing_bcl  : if True then run bcl2fastq with --ignore-missing-bcl
        # ignore_missing_stats: if True then run bcl2fastq with --ignore-missing-stats
        # skip_rsync          : if True then don't rsync primary data at the start of
        #                       bcl2fastq conversion
        # keep_primary_data   : if True then don't remove primary data at
        #                       the end of bcl2fastq conversion
        #
        # Directories
        analysis_dir = self.params.analysis_dir
        self.params['unaligned_dir'] = 'bcl2fastq'
        bcl2fastq_dir = self.add_directory(self.params.unaligned_dir)
        sample_sheet = self.params.sample_sheet
        if self.verify_bcl_to_fastq():
            print "Bcl to fastq outputs already present, nothing to do"
            return
        # Fetch primary data
        if not skip_rsync:
            if self.get_primary_data() != 0:
                logging.error("Failed to acquire primary data")
                raise Exception, "Failed to acquire primary data"
        primary_data = os.path.join(self.params.primary_data_dir,
                                    os.path.basename(self.params.data_dir))
        # Get info about the run
        print "Primary data dir    : %s" % primary_data
        illumina_run = IlluminaData.IlluminaRun(primary_data)
        bases_mask = self.params.bases_mask
        nmismatches = bclToFastq.get_nmismatches(bases_mask)
        print "%s" % illumina_run.run_dir
        print "Platform            : %s" % illumina_run.platform
        print "Bcl format          : %s" % illumina_run.bcl_extension
        print "Sample sheet        : %s" % sample_sheet
        print "Bases mask          : %s" % bases_mask
        print "Nmismatches         : %d (determined from bases mask)" % nmismatches
        print "Ignore missing bcl  : %s" % ignore_missing_bcl
        print "Ignore missing stats: %s" % ignore_missing_bcl
        # Set up runner
        runner = auto_process_settings.runners.bcl2fastq
        runner.log_dir(self.log_dir)
        # Run bcl2fastq
        bcl2fastq = applications.Command('bclToFastq.py',
                                         '--use-bases-mask',bases_mask,
                                         '--nmismatches',nmismatches,
                                         '--ignore-missing-control')
        if ignore_missing_bcl:
            bcl2fastq.add_args('--ignore-missing-bcl')
        if ignore_missing_stats:
            bcl2fastq.add_args('--ignore-missing-stats')
        bcl2fastq.add_args(primary_data,
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
        # Remove primary data
        if not keep_primary_data:
            self.remove_primary_data()
        # Generate statistics
        self.generate_stats()
        # Make a 'projects.info' metadata file
        self.make_project_metadata_file()

    def generate_stats(self,stats_file='statistics.info'):
        # Generate statistics for initial fastq files from bcl2fastq
        print "Generating statistics:"
        self.params['stats_file'] = stats_file
        stats = illumina_data_statistics(
            IlluminaData.IlluminaData(self.params.analysis_dir,
                                      unaligned_dir=self.params.unaligned_dir))
        stats.write(os.path.join(self.params.analysis_dir,stats_file),
                    include_header=True)
        self.params['stats_file'] = stats_file
        print "Statistics in %s" % self.params.stats_file

    def remove_primary_data(self):
        # Remove primary data
        primary_data = os.path.join(self.params.primary_data_dir,
                                    os.path.basename(self.params.data_dir))
        if os.path.isdir(primary_data):
            print "Removing copy of primary data in %s" % primary_data
            shutil.rmtree(primary_data)

    def verify_bcl_to_fastq(self):
        # Check that bcl to fastq outputs match sample sheet predictions
        bcl_to_fastq_dir = os.path.join(self.analysis_dir,self.params.unaligned_dir)
        if not os.path.isdir(bcl_to_fastq_dir):
            # Directory doesn't exist
            return False
        # Try to create an IlluminaData object
        try:
            illumina_data = IlluminaData.IlluminaData(self.analysis_dir,
                                                      unaligned_dir=self.params.unaligned_dir)
        except IlluminaData.IlluminaDataError, ex:
            # Failed to initialise
            logging.debug("Failed to get information from %s: %s" % (bcl_to_fastq_dir,ex))
            return False
        # Do check
        return analyse_illumina_run.verify_run_against_sample_sheet(illumina_data,
                                                                    self.params.sample_sheet)

    def setup_analysis_dirs(self):
        # Construct and populate the analysis directories for each project
        illumina_data = IlluminaData.IlluminaData(self.analysis_dir,
                                                  unaligned_dir=self.params.unaligned_dir)
        project_metadata = TabFile.TabFile(os.path.join(analysis_dir,
                                                        self.params.project_metadata),
                                           first_line_is_header=True)
        for line in project_metadata:
            # Acquire the run name
            if self.params.data_dir is not None:
                run_name = os.path.basename(self.params.data_dir)
            else:
                run_name = os.path.basename(self.params.analysis_dir)
            # Look up project data
            project_name = line['Project']
            user = line['User']
            PI = line['PI']
            organism = line['Organism']
            library_type = line['Library']
            comments = line['Comments']
            # Create the project
            project = auto_process_utils.AnalysisProject(project_name,
                                                         os.path.join(self.analysis_dir,
                                                                      project_name),
                                                         user=user,
                                                         PI=PI,
                                                         organism=organism,
                                                         library_type=library_type,
                                                         run=run_name,
                                                         comments=comments,
                                                         platform=self.params.platform)
            project.create_directory(illumina_data.get_project(project_name))

    def run_qc(self,projects=None):
        # Run QC pipeline for all projects
        #
        # Tests whether QC outputs already exist and only runs
        # QC for those files where the outputs are not all present
        #
        # Setup a pipeline runner
        qc_runner = auto_process_settings.runners.qc
        pipeline = Pipeline.PipelineRunner(qc_runner)
        # Get project dir data
        projects = self.get_analysis_projects(projects)
        # Check we have projects
        if len(projects) == 0:
            raise Exception, "No projects found for QC analysis"
        # Look for samples with no/invalid QC outputs and populate
        # pipeline with the associated fastq.gz files
        for project in projects:
            print "*** Setting up QC for %s ***" % project.name
            # Make the qc directory if it doesn't exist
            qc_dir = os.path.join(project.dirn,'qc')
            if not os.path.exists(qc_dir):
                print "Making 'qc' subdirectory"
                bcf_utils.mkdir(qc_dir,mode=0775)
            # Loop over samples and queue up those where the QC
            # isn't validated
            for sample in project.samples:
                print "Examining files in sample %s" % sample.name
                for fq in sample.fastq:
                    if sample.verify_qc(qc_dir,fq):
                        logging.debug("\t%s: QC verified" % fq)
                    else:
                        print "\t%s: setting up QC run" % fq
                        group = "%s.%s" % (project.name,sample.name)
                        fastq = os.path.join(project.dirn,'fastqs',fq)
                        label = str(auto_process_utils.AnalysisFastq(fq))
                        pipeline.queueJob(project.dirn,'illumina_qc.sh',(fastq,),
                                          label=label,group=group)
        # Run the pipeline
        if pipeline.nWaiting() > 0:
            pipeline.run()
        # Verify the outputs
        for project in projects:
            if not project.verify_qc():
                logging.error("QC failed for one or more samples in %s" % project.name)
            else:
                print "QC okay, generating report for %s" % project.name
                project.qc_report

    def copy_to_archive(self):
        # Copy the analysis directory and contents to an archive area
        archive_dir = auto_process_settings.archive.dirn
        if archive_dir is None:
            raise Exception, "No archive directory specified in settings"
        # Construct subdirectory structure i.e. platform and year
        platform = self.params.platform
        year = time.strftime("%Y")
        archive_dir = os.path.join(archive_dir,year,platform)
        print "Copying to archive directory: %s" % archive_dir
        try:
            rsync = applications.general.rsync(self.analysis_dir,archive_dir,
                                               prune_empty_dirs=True)
            print "Running %s" % rsync
            status = rsync.run_subprocess(log=self.log_path('rsync.archive.log'))
        except Exception, ex:
            logging.error("Exception rsyncing to archive: %s" % ex)
            status = -1
        if status != 0:
            logging.error("Failed to rsync to archive (returned status %d)" % status)

    def log_analysis(self):
        # Add a record of the analysis to the logging file
        raise NotImplementedError

    def publish_qc(self):
        # Copy the QC reports to the webserver
        # Get project dir data
        if auto_process_settings.qc_web_server.server is None or \
           auto_process_settings.qc_web_server.webdir is None:
            raise Exception, "No server and/or directory specifed in settings"
        user = auto_process_settings.qc_web_server.user
        server = auto_process_settings.qc_web_server.server
        webdir = os.path.join(auto_process_settings.qc_web_server.webdir,
                              os.path.basename(self.analysis_dir))
        # Get project data
        projects = self.get_analysis_projects()
        # Make a remote directory for the QC reports
        try:
            mkdir_cmd = applications.general.ssh_command(user,server,('mkdir',webdir))
            mkdir_cmd.run_subprocess()
        except Exception, ex:
            raise Exception, "Exception making remote directory for QC reports: %s" % ex
        # Deal with QC for each project
        for project in projects:
            if project.verify_qc():
                try:
                    # scp the qc zip file
                    qc_zip = os.path.join(project.dirn,project.qc_report)
                    scp = applications.general.scp(user,server,qc_zip,webdir)
                    print "Running %s" % scp
                    scp.run_subprocess()
                    # Unpack at the other end
                    unzip_cmd = applications.general.ssh_command(
                        user,server,
                        ('unzip','-q','-o','-d',webdir,
                         os.path.join(webdir,project.qc_report)))
                    print "Running %s" % unzip_cmd
                    unzip_cmd.run_subprocess()
                except Exception, ex:
                    print "Failed to copy QC report: %s" % ex

class ProjectMetadataFile(TabFile.TabFile):
    def __init__(self,filen=None):
        self.__filen = filen
        TabFile.TabFile.__init__(self,filen=filen,
                                 column_names=('Project',
                                               'Samples',
                                               'User',
                                               'Library',
                                               'Organism',
                                               'PI',
                                               'Comments'),
                                 first_line_is_header=True)

    def add_project(self,project_name,sample_names,user=None,
                    library_type=None,organism=None,PI=None,
                    comments=None):
        # Add project info to the metadata file
        self.append(data=(project_name,
                          bcf_utils.pretty_print_names(sample_names),
                          '.' if user is None else user,
                          '.' if library_type is None else library_type,
                          '.' if organism is None else organism,
                          '.' if PI is None else PI,
                          '.' if comments is None else comments))

    def save(self,filen=None):
        # Save the data back to file
        if filen is not None:
            self.__filen = filen
        self.write(filen=self.__filen,include_header=True)

#######################################################################
# Functions
#######################################################################

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
                print "* %s" % fastq
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

def list_available_commands(cmds):
    # Pretty-print available commands
        print ""
        print "Available commands are:"
        for cmd in cmds:
            print "\t%s" % cmd
        print ""

def set_debug(flag):
    # Turn on debug output
    if flag:
        logging.getLogger().setLevel(logging.DEBUG)

# Command line parsers

def setup_parser():
    p = optparse.OptionParser(usage="%prog setup [OPTIONS] DIR",
                              version="%prog "+__version__,
                              description="Automatically process Illumina sequencing "
                              "data from DIR.")
    p.add_option('--fastq-dir',action='store',dest='fastq_dir',default=None,
                 help="Import fastq.gz files from FASTQ_DIR (which should be a "
                 "subdirectory of DIR with the same structure as that produced "
                 "by CASAVA/bcl2fastq).")
    p.add_option('--debug',action='store_true',dest='debug',default=False,
                 help="Turn on debugging output from Python libraries")
    return p

def make_fastqs_parser():
    p = optparse.OptionParser(usage="%prog setup [OPTIONS] [ANALYSIS_DIR]",
                              version="%prog "+__version__,
                              description="Automatically process Illumina sequence from "
                              "ANALYSIS_DIR.")
    p.add_option('--skip-rsync',action='store_true',
                 dest='skip_rsync',default=False,
                 help="don't rsync the primary data at the beginning of processing")
    p.add_option('--ignore-missing-bcl',action='store_true',
                 dest='ignore_missing_bcl',default=False,
                 help="Use the --ignore-missing-bcl option for bcl2fastq (treat "
                 "missing bcl files as no call)")
    p.add_option('--ignore-missing-stats',action='store_true',
                 dest='ignore_missing_stats',default=False,
                 help="Use the --ignore-missing-stats option for bcl2fastq (fill in "
                 "with zeroes when *.stats files are missing)")
    p.add_option('--keep-primary-data',action='store_true',
                 dest='keep_primary_data',default=False,
                 help="Don't delete the primary data at the end of processing")
    p.add_option('--debug',action='store_true',dest='debug',default=False,
                 help="Turn on debugging output from Python libraries")
    return p

def run_qc_parser():
    p = optparse.OptionParser(usage="%prog setup [OPTIONS] [ANALYSIS_DIR]",
                              version="%prog "+__version__,
                              description="Automatically process Illumina sequence from "
                              "ANALYSIS_DIR.")
    p.add_option('--projects',action='store',
                 dest='projects',default=None,
                 help="simple wildcard-based pattern specifying a subset of "
                 "projects to run the QC on")
    p.add_option('--debug',action='store_true',dest='debug',default=False,
                 help="Turn on debugging output from Python libraries")
    return p

def generic_parser():
    p  = optparse.OptionParser(usage="%prog setup [OPTIONS] [ANALYSIS_DIR]",
                              version="%prog "+__version__,
                              description="Automatically process Illumina sequence from "
                              "ANALYSIS_DIR.")
    p.add_option('--debug',action='store_true',dest='debug',default=False,
                 help="Turn on debugging output from Python libraries")
    return p

#######################################################################
# Main program
#######################################################################

if __name__ == "__main__":

    # Available commands and corresponding
    cmd_parsers = bcf_utils.OrderedDictionary()
    cmd_parsers['setup'] = setup_parser()
    cmd_parsers['make_fastqs'] = make_fastqs_parser()
    cmd_parsers['setup_analysis_dirs'] = generic_parser()
    cmd_parsers['run_qc'] = run_qc_parser()
    cmd_parsers['archive'] = generic_parser()
    cmd_parsers['publish_qc'] = generic_parser()

    # Process command line
    try:
        cmd = sys.argv[1]
    except IndexError:
        cmd = "help"
    if cmd == "help" or cmd == "--help" or cmd == "-h":
        list_available_commands(cmd_parsers)
        sys.exit(0)
    else:
        err = None
        if cmd not in cmd_parsers:
            err = "Unrecognised command '%s'" % cmd
        elif len(sys.argv) < 2:
            err = "Need to supply a command"
        if err is not None:
            print err
            list_available_commands(cmd_parsers)
            sys.stderr.write("%s\n" % err)
            sys.exit(1)
    options,args = cmd_parsers[cmd].parse_args(sys.argv[2:])

    # Report name and version
    print "%s version %s" % (os.path.basename(sys.argv[0]),__version__)

    # Turn on debugging?
    set_debug(options.debug)

    # Setup the processing object and run the requested command
    if cmd == 'setup':
        if len(args) != 1:
            sys.stderr.write("Need to supply a data source location\n")
            sys.exit(1)
        d = AutoProcess()
        if options.fastq_dir is None:
            d.setup(args[0])
        else:
            d.setup_from_fastq_dir(args[0],options.fastq_dir)
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
            d.bcl_to_fastq(skip_rsync=options.skip_rsync,
                           keep_primary_data=options.keep_primary_data,
                           ignore_missing_bcl=options.ignore_missing_bcl,
                           ignore_missing_stats=options.ignore_missing_stats)
        elif cmd == 'setup_analysis_dirs':
            d.setup_analysis_dirs()
        elif cmd == 'run_qc':
            d.run_qc(projects=options.projects)
        elif cmd == 'archive':
            d.copy_to_archive()
        elif cmd == 'publish_qc':
            d.publish_qc()
