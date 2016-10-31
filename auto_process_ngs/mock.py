#     mock.py: module providing mock Illumina data for testing
#     Copyright (C) University of Manchester 2012-2016 Peter Briggs
#
########################################################################

"""
mock.py

Provides classes for mocking up examples of inputs and outputs for
various parts of the process pipeline (including example directory
structures), to be used in testing.

These include:

- MockAnalysisDir
- MockAnalysisProject
- MockAnalysisDirFactory

"""

#######################################################################
# Import modules that this module depends on
#######################################################################

import os
import bcftbx.utils
from bcftbx.mock import MockIlluminaData

#######################################################################
# Class definitions
#######################################################################

class MockAnalysisDir(MockIlluminaData):
    """
    Utility class for creating mock auto-process analysis directories

    The MockAnalysisDir class allows artificial analysis directories
    to be defined, created and populated, and then destroyed.

    These artifical directories are intended to be used for testing
    purposes.

    Two styles of analysis directories can be produced: 'casava'-style
    aims to mimic that produced from the CASAVA and bcl2fastq 1.8
    processing software; 'bcl2fastq2' mimics that from the bcl2fastq
    2.* software.

    Basic example usage:

    >>> mockdir = MockAnalysisDir('130904_PJB_XXXXX','miseq',fmt='casava')
    >>> mockdir.add_fastq_batch('PJB','PJB1','PJB1_GCCAAT',lanes=[1,])
    >>> ...
    >>> mockdir.create()

    This will make a CASAVA-style directory structure like:

    1130904_PJB_XXXXX/
        metadata.info
        bcl2fastq/
            Project_PJB/
                Sample_PJB1/
                    PJB1_GCCAAT_L001_R1_001.fastq.gz
                    PJB1_GCCAAT_L001_R2_001.fastq.gz
        PJB/
            README.info
            fastqs/
                PJB1_GCCAAT_L001_R1_001.fastq.gz
                PJB1_GCCAAT_L001_R2_001.fastq.gz
        ...

    To delete the physical directory structure when finished:

    >>> mockdata.remove()

    """
    def __init__(self,
                 run_name,
                 platform,
                 unaligned_dir='bcl2fastq',
                 fmt='bcl2fastq2',
                 paired_end=True,
                 no_lane_splitting=False,
                 no_undetermined=False,
                 top_dir=None,
                 metadata=None,
                 readme=None):
        # Make a mock-up of an analysis dir
        self.run_name = os.path.basename(str(run_name))
        self.platform = str(platform).lower()
        self.readme = readme
        # Store metadata
        self.metadata = { 'run_name': self.run_name,
                          'platform': self.platform, }
        if metadata is not None:
            for item in metadata:
                self.metadata[item] = metadata[item]
        name = "%s_analysis" % run_name
        if top_dir is None:
            top_dir = os.getcwd()
        MockIlluminaData.__init__(self,name,fmt,
                                  unaligned_dir=unaligned_dir,
                                  paired_end=paired_end,
                                  no_lane_splitting=no_lane_splitting,
                                  top_dir=top_dir)
        # Add undetermined
        if not no_undetermined:
            self.add_undetermined()

    def create(self):
        """
        Build and populate the directory structure

        Creates the directory structure on disk which has been defined
        within the MockAnalysisDir object.

        Invoke the 'remove' method to delete the directory structure.

        The contents of the MockAnalysisDir object can be modified
        after the directory structure has been created, but changes will
        not be reflected on disk. Instead it is necessary to first
        remove the directory structure, and then re-invoke the create
        method.

        create raises an OSError exception if any part of the directory
        structure already exists.

        """
        MockIlluminaData.create(self)
        # Add (empty) metadata file
        with open(os.path.join(self.dirn,'metadata.info'),'w') as fp:
            if self.metadata is not None:
                for item in self.metadata:
                    fp.write("%s\t%s\n" % (item,self.metadata[item]))
            else:
                fp.write('')
        # Add auto_process.info file
        with open(os.path.join(self.dirn,'auto_process.info'),'w') as fp:
            fp.write("analysis_dir\t%s\n" % os.path.basename(self.dirn))
            fp.write("bases_mask\ty76,I8,I8,y76\n")
            fp.write("data_dir\t/mnt/data/%s\n" % self.run_name)
            fp.write("per_lane_stats_file\tper_lane_statistics.info\n")
            fp.write("primary_data_dir\t%s/primary_data/%s\n" % (self.dirn,
                                                                 self.run_name))
            fp.write("project_metadata\tprojects.info\n")
            fp.write("sample_sheet\t%s/custom_SampleSheet.csv\n" % self.dirn)
            fp.write("stats_file\tstatistics.info\n")
            fp.write("unaligned_dir\tbcl2fastq\n")
        # Add top-level README file
        if self.readme is not None:
            open(os.path.join(self.dirn,'README'),'w').write(self.readme)
        # Add empty original sample sheet
        open(os.path.join(self.dirn,'SampleSheet.orig.csv'),'w').write('')
        # Initialise a custom_SampleSheet.csv
        with open(os.path.join(self.dirn,'custom_SampleSheet.csv'),'w') as fp:
            fp.write('[Data]\n')
            fp.write('Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description\n')
        # Add top-level ScriptCode directory
        os.mkdir(os.path.join(self.dirn,'ScriptCode'))
        # Add top-level logs directory
        os.mkdir(os.path.join(self.dirn,'logs'))
        # Add project dirs
        projects_info = open(os.path.join(self.dirn,'projects.info'),'w')
        projects_info.write('#%s\n' % '\t'.join(('Project',
                                                 'Samples',
                                                 'User',
                                                 'Library',
                                                 'Organism',
                                                 'PI',
                                                 'Comments')))
        for project in self.projects:
            if project.startswith("Undetermined"):
                project_name = 'undetermined'
            else:
                project_name = project
            project_dir = MockAnalysisProject(project_name)
            sample_names = []
            for sample in self.samples_in_project(project):
                sample_names.append(sample)
                for fq in self.fastqs_in_sample(project,sample):
                    project_dir.add_fastq(fq)
            # Add line to projects.info
            if project_name != 'undetermined':
                projects_info.write('%s\n' % '\t'.join((project,
                                                        ','.join(sample_names),
                                                        '.',
                                                        '.',
                                                        '.',
                                                        '.',
                                                        '.')))
                # Add lines to custom_SampleSheet
                with open(os.path.join(self.dirn,'custom_SampleSheet.csv'),
                          'a') as fp:
                    for sample in self.samples_in_project(project):
                        fp.write('%s,,,,,,%s,\n' % (sample,
                                                    project_name))
            # Write the project directory to disk
            project_dir.create(top_dir=self.dirn)
        # Finished
        return self.dirn

class MockAnalysisProject(object):
    """
    Utility class for creating mock auto-process project directories

    Example usage:

    >>> m = MockAnalysisProject('PJB',('PJB1_S1_R1_001.fasta.gz,
    ...                                'PJB1_S1_R2_001.fasta.gz))
    >>> m.create()

    """
    def __init__(self,name,fastq_names=None):
        """
        Create a new MockAnalysisProject instance
        """
        self.name = name
        if fastq_names is None:
            fastq_names = []
        self.fastq_names = [fq for fq in fastq_names]

    def add_fastq(self,fq):
        """
        Add a Fastq file to the project
        """
        self.fastq_names.append(fq)

    def create(self,top_dir=None,readme=True,scriptcode=True):
        """
        Build and populate the directory structure

        Arguments:
          top_dir (str): path to directory to create project
            directory underneath (default is pwd)
          readme (boolean): if True then write a README file
          scriptcode (boolean): if True then write a ScriptCode
            subdirectory

        """
        # Create directory
        if top_dir is None:
            top_dir = os.getcwd()
        project_dir = os.path.join(top_dir,self.name)
        os.mkdir(project_dir)
        # Create fastqs subdirectory
        fqs_dir = os.path.join(project_dir,'fastqs')
        os.mkdir(fqs_dir)
        # Add Fastq files
        for fq in self.fastq_names:
            fq = os.path.basename(fq)
            with open(os.path.join(fqs_dir,fq),'w') as fp:
                fp.write('')
        # Add (empty) README.info
        if readme:
            open(os.path.join(project_dir,'README.info'),
                 'w').write('')
        # Add ScriptCode directory
        if scriptcode:
            os.mkdir(os.path.join(project_dir,'ScriptCode'))

class MockAnalysisDirFactory(object):
    """
    Collection of convenient pre-populated test cases
    """
    @classmethod
    def bcl2fastq2(self,
                   run_name,platform,
                   paired_end=True,
                   no_lane_splitting=True,
                   top_dir=None):
        """
        Basic analysis dir from bcl2fastq v2
        """
        mad = MockAnalysisDir(run_name,platform,
                              unaligned_dir='bcl2fastq',
                              fmt='bcl2fastq2',
                              paired_end=paired_end,
                              no_lane_splitting=no_lane_splitting,
                              top_dir=top_dir)
        if no_lane_splitting:
            lanes = []
        else:
            lanes = [1,2,3,4]
        mad.add_fastq_batch('AB','AB1','AB1_S1',lanes=lanes)
        mad.add_fastq_batch('AB','AB2','AB2_S2',lanes=lanes)
        mad.add_fastq_batch('CDE','CDE3','CDE3_S3',lanes=lanes)
        mad.add_fastq_batch('CDE','CDE4','CDE4_S4',lanes=lanes)
        return mad

    @classmethod
    def casava(self,
               run_name,platform,
               paired_end=True,
               top_dir=None):
        """
        Basic analysis dir from CASAVA/bcl2fastq v1.8
        """
        mad = MockAnalysisDir(run_name,platform,
                              unaligned_dir='bcl2fastq',
                              fmt='casava',
                              paired_end=paired_end,
                              top_dir=top_dir)
        lanes = [1,2,3,4]
        mad.add_fastq_batch('AB','AB1','AB1_GCCAAT',lanes=lanes)
        mad.add_fastq_batch('AB','AB2','AB2_AGTCAA',lanes=lanes)
        mad.add_fastq_batch('CDE','CDE3','CDE3_GCCAAT',lanes=lanes)
        mad.add_fastq_batch('CDE','CDE4','CDE4_AGTCAA',lanes=lanes)
        return mad
