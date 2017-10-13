#     mock.py: module providing mock Illumina data for testing
#     Copyright (C) University of Manchester 2012-2016 Peter Briggs
#
########################################################################

"""
mock.py

Provides classes for mocking up examples of inputs and outputs for
various parts of the process pipeline (including example directory
structures), to be used in testing.

The core classes are:

- MockAnalysisDir: create mock auto-process analysis directories
- MockAnalysisProject: create mock auto-process project
  directories

These can be used to configure and create mock directories mimicking
"minimal" versions of analysis directories and projects.

Additional mock artefacts (e.g. QC outputs, barcode analysis etc)
can be added to the mock directories once they have been created,
using the "updater" classes:

- UpdateAnalysisDir: add artefacts to an analysis directory
- UpdateAnalysisProject: add artefacts to an analysis project

Finally there is a factory class which provides methods to quickly
make "default" analysis directories for testing:

- MockAnalysisDirFactory

"""

#######################################################################
# Import modules that this module depends on
#######################################################################

import os
import bcftbx.utils
from bcftbx.mock import MockIlluminaData
from .utils import AnalysisProject
from .utils import ZipArchive
from .qc.illumina_qc import expected_qc_outputs


#######################################################################
# Classes for making mock directories
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
                 lanes=None,
                 no_lane_splitting=False,
                 no_undetermined=False,
                 top_dir=None,
                 metadata=None,
                 readme=None):
        """
        Create a mock-up of an analysis directory

        Arguments:
          run_name (str): name for the run e.g.
            '1130904_PJB_XXXXX'
          platform (str): name for the platform
            e.g. 'nextseq'
          unaligned_dir (str): for the bcl2fastq
            output directory (default: 'bcl2fastq')
          fmt (str): format of the outputs (can be
            'casava' or 'bcl2fastq2'; default is
            'bcl2fastq')
          paired_end (bool): whether run should be
            paired end (set True, default) or single
            end (set False)
          lanes (list): if not None then specify a
            list of lane numbers to include
          no_lane_splitting (bool): whether to
            mimic the '--no-lane-splitting' option
            of bcl2fastq2 in generating Fastq names
            (default: make separate Fastqs for each
            lane)
          no_undetermined (bool): whether to
            include 'undetermined' outputs (default:
            False, do include 'undetermined'
            outputs)
          top_dir (str): set parent directory to
            make the mock analysis directory in
            (default: current working directory)
          metadata (dict): if set then should be
            a dictionary of metadata items with
            corresponding values, which will be
            written to the metadata.info file
          readme (str): if set then will be
            written to a 'README' file in the mock
            analysis directory
        """
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
            self.add_undetermined(lanes=lanes)

    def create(self,no_project_dirs=False):
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

        'create' raises an OSError exception if any part of the
        directory structure already exists.

        Arguments:
          no_project_dirs (bool): if False then don't create
            analysis project subdirectories (these are created by
            default)
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
            if not no_project_dirs:
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
    def __init__(self,name,fastq_names=None,fastq_dir=None,metadata=dict()):
        """
        Create a new MockAnalysisProject instance
        """
        self.name = name
        if fastq_names is None:
            fastq_names = []
        self.fastq_names = [fq for fq in fastq_names]
        if fastq_dir is None:
            self.fastq_dir = 'fastqs'
        else:
            self.fastq_dir = fastq_dir
        self.metadata = metadata

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
        fqs_dir = os.path.join(project_dir,self.fastq_dir)
        os.mkdir(fqs_dir)
        # Add Fastq files
        for fq in self.fastq_names:
            fq = os.path.basename(fq)
            with open(os.path.join(fqs_dir,fq),'w') as fp:
                fp.write('')
        # Add README.info
        if readme:
            with open(os.path.join(project_dir,'README.info'),'w') as info:
                for key in self.metadata:
                    info.write("%s\t%s" % (key,self.metadata[key]))
        # Add ScriptCode directory
        if scriptcode:
            os.mkdir(os.path.join(project_dir,'ScriptCode'))

#######################################################################
# Classes for updating directories with mock artefacts
#######################################################################

class DirectoryUpdater(object):
    """
    Base class for updating mock directories

    Provides the following methods:

    - add_subdir: adds arbitrary subdirectory
    - add_file: adds arbitrary file
    """
    def __init__(self,base_dir):
        """
        Create DirectoryUpdater instance

        Must supply a "base directory" on instantiation;
        this is the directory within which subdirectories
        and files will be created.

        Arguments:
          base_dir (path): path to the base directory
        """
        self._base_dir = os.path.abspath(base_dir)

    def add_subdir(self,dirn):
        """
        Add an arbitrary directory to the base dir

        Arguments:
          dirn (str): path of directory to add
        """
        if not os.path.isabs(dirn):
            dirn = os.path.join(self._base_dir,
                                dirn)
        print "Adding directory: %s" % dirn
        os.mkdir(dirn)

    def add_file(self,filen,content=None):
        """
        Add an arbitrary file to the base dir

        Arguments:
          filen (str): path of file to add
          content (str): if supplied then will be
            written as content of new file
        """
        if not os.path.isabs(filen):
            filen = os.path.join(self._base_dir,
                                 filen)
        print "Adding file: %s" % filen
        with open(filen,'w') as fp:
            if content is not None:
                fp.write("%s" % content)
            else:
                fp.write("")

class UpdateAnalysisDir(DirectoryUpdater):
    """
    Utility class to add mock artefacts to an AnalysisDir

    Provides the following methods:

    - add_processing_report
    - add_barcode_analysis
    - add_cellranger_qc_output

    Example usage:

    >>> m = MockAnalysisDirFactory.bcl2fastq2(
    ...     '160621_M00879_0087_000000000-AGEW9',
    ...     'miseq')
    >>> m.create()
    >>> ap = AutoProcess(m.dirn)
    >>> UpdateAnalysisDir(ap).add_processing_report()

    """
    def __init__(self,ap):
        """
        Create a new UpdateAnalysisDir instance

        Arguments:
          ap (AutoProcess): AutoProcess or
            AnalysisDir instance to operate on
        """
        self._ap = ap
        DirectoryUpdater.__init__(self,ap.analysis_dir)

    def add_processing_report(self):
        """
        Add a 'processing_qc.html' file
        """
        self.add_file("processing_qc.html")

    def add_barcode_analysis(self,
                             barcode_analysis_dir="barcode_analysis"):
        """
        Add mock barcode analysis outputs

        Arguments:
          barcode_analysis_dir (str): name of
            barcode analysis subdirectory (default:
            'barcode_analysis'
        """
        self.add_subdir(barcode_analysis_dir)
        for f in ("barcodes.report",
                  "barcodes.xls",
                  "barcodes.html"):
            self.add_file(os.path.join(barcode_analysis_dir,f))

    def add_cellranger_qc_output(self,lanes=None):
        """
        Add mock cellranger QC report

        Arguments:
          lanes (str): optional, specify lane numbers for
            the report
        """
        self.add_file("cellranger_qc_summary%s.html"
                      % ("_%s" % lanes
                         if lanes is not None
                         else ""))

class UpdateAnalysisProject(DirectoryUpdater):
    """
    Utility class to add mock artefacts to an AnalysisDir

    Provides the following methods:

    - add_fastq_set
    - add_qc_outputs
    - add_icell8_outputs
    - add_cellranger_count_outputs

    Example usage:

    >>> m = MockAnalysisProject("PJB",('PJB1_S1_R1_001.fasta.gz,
    ...                                'PJB1_S1_R2_001.fasta.gz))
    >>> m.create()
    >>> p = AnalysisProject(m.name,m.name)
    >>> UpdateAnalysisProject(p).add_qc_outputs()

    """
    def __init__(self,project):
        """
        Create a new UpdateAnalysisProject instance

        Arguments:
          project (AnalysisProject): AnalysisProject
            instance to operate on
        """
        DirectoryUpdater.__init__(self,project.dirn)
        self._project = project
        self._reload_project()

    def _parent_dir(self):
        """
        Return parent directory of project directory
        """
        return os.path.dirname(self._project.dirn)

    def _reload_project(self):
        """
        Reload the project

        Creates new local AnalysisProject instance
        refreshed from the file system
        """
        project = AnalysisProject(
            self._project.name,
            self._project.dirn)
        self._project = project

    def add_fastq_set(self,fastq_set,fastqs):
        """
        Add an additional fastq set

        Arguments:
          fastq_set (str): name of the new Fastq
            set/subdirectory
          fastqs (list): list of Fastq filenames
            to create in the new set
        """
        self.add_subdir(fastq_set)
        for fq in fastqs:
            self.add_file(os.path.join(fastq_set,fq))
        self._reload_project()

    def add_qc_outputs(self,fastq_set=None,qc_dir=None,
                       include_multiqc=True):
        """
        Add mock QC outputs

        Arguments:
          fastq_set (str): specify non-default Fastq
            set to make QC outputs
        """
        print("Adding mock QC outputs to %s" %
              self._project.dirn)
        # Handle Fastq set
        if fastq_set is not None:
            self._project.use_fastq_dir(fastq_dir=fastq_set)
        print "- Fastq set: %s" % self._project.fastq_dir
        # Handle QC directory
        self._project.setup_qc_dir(qc_dir=qc_dir,
                                   fastq_dir=fastq_set)
        if qc_dir is not None:
            self._project.use_qc_dir(qc_dir)
        print "- QC dir: %s" % self._project.qc_dir
        # Generate base QC outputs (one set per fastq)
        for fq in self._project.fastqs:
            print "Adding outputs for %s" % fq
            for f in expected_qc_outputs(fq,self._project.qc_dir):
                self.add_file(f)
        # Make mock report
        fastq_set_name = os.path.basename(self._project.fastq_dir)[6:]
        qc_name = "qc%s_report" % fastq_set_name
        self.add_file("%s.html" % qc_name)
        # Make mock ZIP archive
        analysis_name = os.path.basename(self._parent_dir())
        report_zip = os.path.join(self._project.dirn,
                                  "%s.%s.%s.zip" %
                                  (qc_name,
                                   self._project.name,
                                   analysis_name))
        print "Building ZIP archive: %s" % report_zip
        zip_file = ZipArchive(report_zip,
                              relpath=self._project.dirn,
                              prefix="%s.%s.%s" %
                              (qc_name,
                               self._project.name,
                               analysis_name))
        zip_file.add_file(os.path.join(self._project.dirn,
                                       "%s.html" % qc_name))
        zip_file.add(self._project.qc_dir)
        zip_file.close()
        # MultiQC report
        if include_multiqc:
            multiqc_name = "multiqc%s_report" % fastq_set_name
            self.add_file("%s.html" % multiqc_name)
        # Reload the project
        self._reload_project()

    def add_icell8_outputs(self):
        """
        Add mock ICell8 outputs to the project
        """
        self.add_file("icell8_processing.html")
        self.add_subdir("stats")
        self.add_file(os.path.join("stats","icell8_stats.tsv"))
        self.add_file(os.path.join("stats","icell8_stats.xlsx"))
        self.add_subdir("icell8_processing_data")
        for png in ("poly_g_dist.png",
                    "read_dist.png",
                    "reads_per_stage.png",
                    "samples.png"):
            self.add_file(os.path.join("icell8_processing_data",
                                       png))
        # Build ZIP archive
        analysis_dir = os.path.basename(self._parent_dir())
        icell8_zip = os.path.join(self._project.dirn,
                                  "icell8_processing.%s.%s.zip" %
                                  (self._project.name,
                                   analysis_dir))
        zip_file = ZipArchive(icell8_zip,
                              relpath=self._project.dirn,
                              prefix="icell8_processing.%s.%s" %
                              (self._project.name,
                               analysis_dir))
        zip_file.add_file(os.path.join(self._project.dirn,
                                       "icell8_processing.html"))
        zip_file.add(os.path.join(self._project.dirn,"stats"))
        zip_file.add(os.path.join(self._project.dirn,
                                  "icell8_processing_data"))
        zip_file.close()
        self._reload_project()

    def add_cellranger_count_outputs(self):
        """
        Add mock 'cellranger count' outputs to project
        """
        print "Adding cellranger count outputs for %s" % self._project.dirn
        self.add_file("cellranger_count_report.html")
        self.add_subdir("cellranger_count")
        for sample in self._project.samples:
            sample_dir = os.path.join("cellranger_count",sample.name)
            self.add_subdir(sample_dir)
            self.add_subdir(os.path.join(sample_dir,"outs"))
            self.add_file(os.path.join(sample_dir,
                                       "outs",
                                       "web_summary.html"))
        # Build ZIP archive
        analysis_dir = os.path.basename(self._parent_dir())
        cellranger_zip = os.path.join(self._project.dirn,
                                      "cellranger_count_report.%s.%s.zip" %
                                      (self._project.name,
                                       analysis_dir))
        zip_file = ZipArchive(cellranger_zip,
                              relpath=self._project.dirn,
                              prefix="cellranger_count_report.%s.%s" %
                              (self._project.name,
                               analysis_dir))
        zip_file.add_file(os.path.join(self._project.dirn,
                                       "cellranger_count_report.html"))
        zip_file.add(os.path.join(self._project.dirn,
                                  "cellranger_count"))
        zip_file.close()
        self._reload_project()

#######################################################################
# Factory classes for quickly generating mock directories
#######################################################################

class MockAnalysisDirFactory(object):
    """
    Collection of convenient pre-populated test cases
    """
    @classmethod
    def bcl2fastq2(self,
                   run_name,platform,
                   paired_end=True,
                   no_lane_splitting=True,
                   top_dir=None,
                   metadata=None):
        """
        Basic analysis dir from bcl2fastq v2
        """
        if no_lane_splitting:
            lanes = []
        else:
            lanes = [1,2,3,4]
        mad = MockAnalysisDir(run_name,platform,
                              unaligned_dir='bcl2fastq',
                              fmt='bcl2fastq2',
                              paired_end=paired_end,
                              no_lane_splitting=no_lane_splitting,
                              lanes=lanes,
                              metadata=metadata,
                              top_dir=top_dir)
        mad.add_fastq_batch('AB','AB1','AB1_S1',lanes=lanes)
        mad.add_fastq_batch('AB','AB2','AB2_S2',lanes=lanes)
        mad.add_fastq_batch('CDE','CDE3','CDE3_S3',lanes=lanes)
        mad.add_fastq_batch('CDE','CDE4','CDE4_S4',lanes=lanes)
        return mad

    @classmethod
    def casava(self,
               run_name,platform,
               paired_end=True,
               metadata=None,
               top_dir=None):
        """
        Basic analysis dir from CASAVA/bcl2fastq v1.8
        """
        lanes = [1,2,3,4]
        mad = MockAnalysisDir(run_name,platform,
                              unaligned_dir='bcl2fastq',
                              fmt='casava',
                              paired_end=paired_end,
                              lanes=lanes,
                              metadata=metadata,
                              top_dir=top_dir)
        mad.add_fastq_batch('AB','AB1','AB1_GCCAAT',lanes=lanes)
        mad.add_fastq_batch('AB','AB2','AB2_AGTCAA',lanes=lanes)
        mad.add_fastq_batch('CDE','CDE3','CDE3_GCCAAT',lanes=lanes)
        mad.add_fastq_batch('CDE','CDE4','CDE4_AGTCAA',lanes=lanes)
        return mad
