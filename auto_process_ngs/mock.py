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
        # Add top-level README file
        if self.readme is not None:
            open(os.path.join(self.dirn,'README'),'w').write(self.readme)
        # Add top-level ScriptCode directory
        os.mkdir(os.path.join(self.dirn,'ScriptCode'))
        # Add top-level logs directory
        os.mkdir(os.path.join(self.dirn,'logs'))
        # Add project dirs
        for project in self.projects:
            if project.startswith("Undetermined"):
                project_name = 'undetermined'
            else:
                project_name = project
            os.mkdir(os.path.join(self.dirn,project_name))
            # Add fastqs
            fqs_dir = os.path.join(self.dirn,project_name,'fastqs')
            os.mkdir(fqs_dir)
            for sample in self.samples_in_project(project):
                for fq in self.fastqs_in_sample(project,sample):
                    with open(os.path.join(fqs_dir,fq),'w') as fp:
                        fp.write('')
            # Add (empty) README.info
            open(os.path.join(self.dirn,
                              project_name,
                              'README.info'),'w').write('')
            # Add ScriptCode directory
            os.mkdir(os.path.join(self.dirn,project_name,'ScriptCode'))
        # Finished
        return self.dirn

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

