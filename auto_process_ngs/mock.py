#     mock.py: module providing mock Illumina data for testing
#     Copyright (C) University of Manchester 2012-2022 Peter Briggs
#
########################################################################

"""
Provides classes for mocking up examples of inputs and outputs for
various parts of the process pipeline (including example directory
structures), as well as mock executables, to be used in testing.

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

There is also a convenience factory class which provides methods to
quickly make "default" analysis directories for testing:

- MockAnalysisDirFactory

It is also possible to make mock executables which mimick some of
the external software required for parts of the pipeline:

- MockBcl2fastq2Exe
- MockBclConvertExe
- Mock10xPackageExe
- MockFastqScreen
- MockFastQC
- MockFastqStrandPy
- MockMultiQC
- MockConda

There also is a wrapper for the 'Mock10xPackageExe' class which
is maintained for backwards compatibility:

- MockCellrangerExe

There are supporting standalone functions for mocking outputs:

- make_mock_bcl2fastq2_output: create mock output from `bcl2fastq`
- make_mock_analysis_project: create a mock analysis project directory
"""

#######################################################################
# Import modules that this module depends on
#######################################################################

import os
import sys
import argparse
import uuid
import shutil
import gzip
import bcftbx.utils
from bcftbx.mock import MockIlluminaData
from bcftbx.IlluminaData import IlluminaRun
from bcftbx.IlluminaData import IlluminaRunInfo
from bcftbx.IlluminaData import IlluminaData
from bcftbx.IlluminaData import IlluminaFastq
from bcftbx.IlluminaData import SampleSheet
from bcftbx.IlluminaData import SampleSheetPredictor
from bcftbx.qc.report import strip_ngs_extensions
from .analysis import AnalysisProject
from .analysis import AnalysisFastq
from .fastq_utils import pair_fastqs_by_name
from .tenx_genomics_utils import CellrangerMultiConfigCsv
from .tenx_genomics_utils import flow_cell_id
from .utils import ZipArchive
from .mockqc import MockQCOutputs
from .mockqc import make_mock_qc_dir
from . import mock10xdata

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
                 bases_mask='auto',
                 paired_end=True,
                 lanes=None,
                 no_lane_splitting=False,
                 no_undetermined=False,
                 top_dir=None,
                 params=None,
                 metadata=None,
                 readme=None,
                 project_metadata=None,
                 include_stats_files=False):
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
          bases_mask (str): bases mask string to
            put into 'auto_process.info'
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
          params (dict): if set then should be
            a dictionary of parameter items with
            corresponding values, which will be
            written to the auto_process.info file
          metadata (dict): if set then should be
            a dictionary of metadata items with
            corresponding values, which will be
            written to the metadata.info file
          readme (str): if set then will be
            written to a 'README' file in the mock
            analysis directory
          project_metadata (dict): if set then should
            be a dictionary where keys are names of
            projects and values are dictionaries of
            metadata items, which will be written to
            the README.info file for that project
          include_stats_files (bool): whether to
            include (empty) stats files (default:
            False, don't include stats files)
        """
        # Make a mock-up of an analysis dir
        self.run_name = os.path.basename(str(run_name))
        self.platform = str(platform).lower()
        self.bases_mask = bases_mask
        self.readme = readme
        self.include_stats_files = include_stats_files
        # Store params
        self.params = {}
        if params is not None:
            for item in params:
                self.params[item] = params[item]
        # Store metadata
        self.metadata = { 'run_name': self.run_name,
                          'platform': self.platform, }
        if metadata is not None:
            for item in metadata:
                self.metadata[item] = metadata[item]
        self.project_metadata = dict()
        if project_metadata is not None:
            for project in project_metadata:
                self.project_metadata[project] = \
                            project_metadata[project]
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
        default_params = {
            'analysis_dir': os.path.basename(self.dirn),
            'bases_mask': self.bases_mask,
            'data_dir': "/mnt/data/%s" % self.run_name,
            'primary_data_dir': "%s/primary_data" % self.dirn,
            'project_metadata': "projects.info",
            'sample_sheet': "%s/custom_SampleSheet.csv" % self.dirn,
            'unaligned_dir': "bcl2fastq",
        }
        if self.include_stats_files:
            default_params['per_lane_stats_files'] = "per_lane_statistics.info"
            default_params['stats_file'] = "statistics.info"
        else:
            default_params['per_lane_stats_files'] = "."
            default_params['stats_file'] = '.'
        for p in default_params:
            if p not in self.params:
                self.params[p] = default_params[p]
        with open(os.path.join(self.dirn,'auto_process.info'),'w') as fp:
            for item in self.params:
                fp.write("%s\t%s\n" % (item,self.params[item]))
        # Add (empty) stats files
        if self.include_stats_files:
            for stats_file in ('statistics.info',
                               'statistics_full.info',
                               'per_lane_statistics.info',
                               'per_lane_sample_stats.info',):
                with open(os.path.join(self.dirn,stats_file),'w') as fp:
                    fp.write('')
        # Add top-level README file
        if self.readme is not None:
            open(os.path.join(self.dirn,'README'),'w').write(self.readme)
        # Add empty original sample sheet
        with open(os.path.join(self.dirn,'SampleSheet.orig.csv'),'wt') as fp:
            fp.write('')
        # Initialise a custom_SampleSheet.csv
        with open(os.path.join(self.dirn,'custom_SampleSheet.csv'),'wt') as fp:
            fp.write('[Data]\n')
            fp.write('Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description\n')
        # Add top-level ScriptCode directory
        os.mkdir(os.path.join(self.dirn,'ScriptCode'))
        # Add top-level logs directory
        os.mkdir(os.path.join(self.dirn,'logs'))
        # Add project dirs
        with open(os.path.join(self.dirn,'projects.info'),
                  'wt') as projects_info:
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
                project_metadata = { 'Run': self.run_name,
                                     'Platform': self.platform }
                try:
                    for item in self.project_metadata[project_name]:
                        project_metadata[item] = \
                            self.project_metadata[project_name][item]
                except KeyError:
                    pass
                project_dir = MockAnalysisProject(project_name,
                                                  metadata=project_metadata)
                sample_names = []
                for sample in self.samples_in_project(project):
                    sample_names.append(sample)
                    for fq in self.fastqs_in_sample(project,sample):
                        project_dir.add_fastq(fq)
                # Add line to projects.info
                if project_name != 'undetermined':
                    projects_info.write('%s\n' % '\t'.join(
                        (project,
                         ','.join(sample_names),
                         '.',
                         '.',
                         '.',
                         '.',
                         '.')))
                    # Add lines to custom_SampleSheet
                    with open(os.path.join(self.dirn,
                                           'custom_SampleSheet.csv'),
                              'a') as fp:
                        for sample in self.samples_in_project(project):
                            fp.write('%s,,,,,,%s,\n' % (sample,
                                                        project_name))
                # Write the project directory to disk
                if not no_project_dirs:
                    project_dir.create(top_dir=self.dirn)
        # Finished
        return self.dirn

class MockAnalysisProject:
    """
    Utility class for creating mock auto-process project directories

    Example usage:

    >>> m = MockAnalysisProject('PJB',('PJB1_S1_R1_001.fastq.gz',
    ...                                'PJB1_S1_R2_001.fastq.gz'))
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

    def create(self,top_dir=None,readme=True,scriptcode=True,
               populate_fastqs=True):
        """
        Build and populate the directory structure

        Arguments:
          top_dir (str): path to directory to create project
            directory underneath (default is pwd)
          readme (boolean): if True then write a README file
          scriptcode (boolean): if True then write a ScriptCode
            subdirectory
          populate_fastqs (boolean): if True then write content
            to the Fastq files
        """
        # Create directory
        if top_dir is None:
            top_dir = os.getcwd()
        project_dir = os.path.join(top_dir,self.name)
        os.mkdir(project_dir)
        # Create fastqs subdirectory
        fqs_dir = os.path.join(project_dir,self.fastq_dir)
        if not os.path.exists(fqs_dir):
            os.mkdir(fqs_dir)
        # Add Fastq files
        for fq in self.fastq_names:
            fq = os.path.basename(fq)
            if populate_fastqs:
                read_number = AnalysisFastq(fq).read_number
                lane = AnalysisFastq(fq).lane_number
                read = """@ILLUMINA-545855:49:FC61RLR:%s:1:10979:1695 %s:N:0:TCCTGA
GCATACTCAGCTTTAGTAATAAGTGTGATTCTGGTA
+
IIIIIHIIIGHHIIDGHIIIIIIHIIIIIIIIIIIH\n""" % (lane,read_number)
            else:
                read = ""
            if fq.endswith('.gz'):
                open_func = gzip.open
            else:
                open_func = open
            with open_func(os.path.join(fqs_dir,fq),'wt') as fp:
                fp.write(read)
        # Add README.info
        if readme:
            with open(os.path.join(project_dir,'README.info'),'w') as info:
                for key in self.metadata:
                    info.write("%s\t%s\n" % (key,self.metadata[key]))
        # Add ScriptCode directory
        if scriptcode:
            os.mkdir(os.path.join(project_dir,'ScriptCode'))
        return project_dir

#######################################################################
# Classes for updating directories with mock artefacts
#######################################################################

class DirectoryUpdater:
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
        print("Adding directory: %s" % dirn)
        os.makedirs(dirn)

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
        print("Adding file: %s" % filen)
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
    - add_10x_mkfastq_qc_output
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

    def add_processing_report(self,name="processing_qc.html"):
        """
        Add a 'processing_qc.html' file

        Arguments:
          name (str): optionally, specify a non-standard
            report name (default: 'processing_qc.html')
        """
        self.add_file(name)

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

    def add_10x_mkfastq_qc_output(self,pkg,lanes=None):
        """
        Add mock 10xGenomics mkfastq QC report

        Arguments:
          pkg (str): 10xGenomics package (e.g. 'cellranger',
            'cellranger-atac' etc)
          lanes (str): optional, specify lane numbers for
            the report
        """
        self.add_file("%s_qc_summary%s.html"
                      % (pkg,
                         "_%s" % lanes
                         if lanes is not None
                         else ""))

    def add_cellranger_qc_output(self,lanes=None):
        """
        Add mock cellranger QC report

        Arguments:
          lanes (str): optional, specify lane numbers for
            the report
        """
        self.add_10x_mkfastq_qc_output(pkg="cellranger",
                                       lanes=lanes)

class UpdateAnalysisProject(DirectoryUpdater):
    """
    Utility class to add mock artefacts to an AnalysisDir

    Provides the following methods:

    - add_fastq_set
    - add_qc_outputs
    - add_icell8_outputs
    - add_cellranger_count_outputs
    - add_cellranger_multi_outputs

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
                       protocol="standardPE",
                       include_fastq_strand=True,
                       include_seqlens=True,
                       include_multiqc=True,
                       include_report=True,
                       include_zip_file=True,
                       legacy_screens=False,
                       legacy_zip_name=False):
        """
        Add mock QC outputs

        Arguments:
          fastq_set (str): specify non-default Fastq
            set to make QC outputs
          qc_dir (str): specify non-default QC output
            directory
          protocol (str): specify non-default QC
            protocol to use
          include_fastq_strand (bool): if True then
            add mock fastq_strand.py outputs
          include_seqlens (bool): if True then add
            mock sequence length outputs
          include_multiqc (bool): if True then add
            mock MultiQC outputs
          include_report (bool): if True then add
            mock QC report outputs
          include_zip_file (bool): if True then add
            mock ZIP archive for QC report
          legacy_screens (bool): if True then use
            old-style 'illumina_qc.sh' naming
            conventions for FastqScreen outputs
          legacy_zip_name (bool): if True then use
            old-style naming convention for ZIP file
            with QC outputs
        """
        print("Adding mock QC outputs to %s" %
              self._project.dirn)
        # Handle Fastq set
        if fastq_set is not None:
            self._project.use_fastq_dir(fastq_dir=fastq_set)
        print("- Fastq set: %s" % self._project.fastq_dir)
        # Handle QC directory
        self._project.setup_qc_dir(qc_dir=qc_dir,
                                   fastq_dir=fastq_set)
        if qc_dir is not None:
            self._project.use_qc_dir(qc_dir)
        print("- QC dir: %s" % self._project.qc_dir)
        # Generate base QC outputs (one set per fastq)
        for fq in self._project.fastqs:
            print("Adding outputs for %s" % fq)
            if include_seqlens:
                MockQCOutputs.seqlens(fq,self._project.qc_dir)
            MockQCOutputs.fastqc_v0_11_2(fq,self._project.qc_dir)
            if protocol in ('singlecell',
                            '10x_scRNAseq',
                            '10x_snRNAseq',
                            '10x_Multiome_GEX',) and \
                self._project.fastq_attrs(fq).read_number == 1:
                continue
            elif protocol in ('ParseEvercode',) and \
                 self._project.fastq_attrs(fq).read_number == 2:
                continue
            for screen in ("model_organisms",
                           "other_organisms",
                           "rRNA"):
                MockQCOutputs.fastq_screen_v0_9_2(fq,
                                                  self._project.qc_dir,
                                                  screen,
                                                  legacy=legacy_screens)
        # Handle fastq_strand, if requested
        if include_fastq_strand:
            fastq_strand_conf = os.path.join(self._project.dirn,
                                             "fastq_strand.conf")
            with open(fastq_strand_conf,'w') as fp:
                fp.write("")
            for fq_pair in pair_fastqs_by_name(self._project.fastqs):
                if protocol in ('singlecell',
                                '10x_scRNAseq',
                                '10x_snRNAseq',
                                '10x_Multiome_GEX',):
                    fq = fq_pair[1]
                else:
                    fq = fq_pair[0]
                MockQCOutputs.fastq_strand_v0_0_4(fq,
                                                  self._project.qc_dir)
        # Update protocol in qc.info
        qc_info = self._project.qc_info(self._project.qc_dir)
        qc_info['protocol'] = protocol
        qc_info['fastq_dir'] = self._project.fastq_dir
        if not legacy_screens:
            qc_info['fastq_screens'] = "model_organisms,other_organisms,rRNA"
        qc_info.save()
        # Make mock report
        fastq_set_name = os.path.basename(self._project.fastq_dir)[6:]
        if include_report:
            qc_name = "qc%s_report" % fastq_set_name
            self.add_file("%s.html" % qc_name)
        # MultiQC report
        if include_multiqc:
            MockQCOutputs.multiqc(self._project.dirn,
                                  multiqc_html="multi%s_report.html"
                                  % os.path.basename(self._project.qc_dir))
        # Make mock ZIP archive
        if include_zip_file:
            analysis_name = os.path.basename(self._parent_dir())
            if not legacy_zip_name:
                zip_prefix = "%s.%s%s" % (qc_name,
                                          self._project.name,
                                          ".%s" % self._project.info.run
                                          if self._project.info.run else '')
            else:
                zip_prefix = "%s.%s.%s" % (qc_name,
                                           self._project.name,
                                           analysis_name)
            report_zip = os.path.join(self._project.dirn,
                                      "%s.zip" % zip_prefix)
            print("Building ZIP archive: %s" % report_zip)
            zip_file = ZipArchive(report_zip,
                                  relpath=self._project.dirn,
                                  prefix=zip_prefix)
            zip_file.add_file(os.path.join(self._project.dirn,
                                           "%s.html" % qc_name))
            zip_file.add(self._project.qc_dir)
            zip_file.close()
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

    def add_cellranger_count_outputs(self,qc_dir=None,cellranger='cellranger',
                                     reference_data_path=
                                     "/data/refdata-cellranger-1.2.0",
                                     prefix="cellranger_count",
                                     legacy=False):
        """
        Add mock 'cellranger count' outputs to project

        Arguments:
          qc_dir (str): specify non-default QC output
            directory
          cellranger (str): specify the 10xGenomics software
            package to add outputs for (defaults to
            'cellranger'; alternatives are 'cellranger-atac')
          reference_data_path (str): optionally specify path
            for reference dataset (doesn't need to exist)
          prefix (str): leading subdirectory for cellranger
            count outputs (defaults to 'cellranger_count';
            ignored if 'legacy' style outputs are generated)
          legacy (bool): if True then generate 'legacy'
            style cellranger outputs
        """
        print("Adding cellranger count outputs for %s" % self._project.dirn)
        if legacy:
            print("- legacy mode")
            self.add_file("cellranger_count_report.html")
            self.add_subdir("cellranger_count")
        else:
            if qc_dir is not None:
                self._project.use_qc_dir(qc_dir)
            print("- QC dir: %s" % self._project.qc_dir)
            if not os.path.exists(self._project.qc_dir):
                self.add_subdir(self._project.qc_dir)
        for sample in self._project.samples:
            if legacy:
                MockQCOutputs.cellranger_count(
                    sample.name,
                    self._project.dirn,
                    cellranger=cellranger,
                    reference_data_path=
                    reference_data_path,
                    prefix="cellranger_count")
            else:
                MockQCOutputs.cellranger_count(
                    sample.name,
                    self._project.qc_dir,
                    cellranger=cellranger,
                    reference_data_path=
                    reference_data_path,
                    prefix=prefix)
        # Build ZIP archive
        if legacy:
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
        # Update cellranger reference data in qc.info
        if not legacy:
            qc_info = self._project.qc_info(self._project.qc_dir)
            qc_info['cellranger_refdata'] = reference_data_path
            qc_info.save()
        self._reload_project()

    def add_cellranger_multi_outputs(self,config_csv=None,sample_names=None,
                                     reference_data_path=None,qc_dir=None,
                                     prefix="cellranger_multi"):
        """
        Add mock 'cellranger multi' outputs to project

        If a 10x multiplexing config file is supplied then
        the mock outputs are generated using the data within
        that file; otherwise the sample names and reference
        dataset path should be explicitly supplied.

        Arguments:
          config_csv (str): path to a 10x multiplexing config
            file (if supplied then sample names and reference
            dataset path will be taken from this file)
          sample_names (list): optionally specify list of
            multiplexed sample names (ignored if 'config_csv'
            file is supplied)
          reference_data_path (str): optionally specify path
            to reference dataset (doesn't need to exist;
            ignored if 'config_csv' file is supplied)
          qc_dir (str): specify non-default QC output
            directory
          prefix (str): leading subdirectory for cellranger
            count outputs (defaults to 'cellranger_multi';
            ignored if 'legacy' style outputs are generated)
        """
        print("Adding cellranger multi outputs to %s" % self._project.dirn)
        if qc_dir is not None:
            self._project.use_qc_dir(qc_dir)
        print("- QC dir: %s" % self._project.qc_dir)
        if not os.path.exists(self._project.qc_dir):
            self.add_subdir(self._project.qc_dir)
        # Read in multiplexing config
        if config_csv:
            config = CellrangerMultiConfigCsv(config_csv)
            sample_names = [s for s in config.sample_names]
            reference_data_path = config.reference_data_path
        # Add outputs
        MockQCOutputs.cellranger_multi(sample_names,
                                       self._project.qc_dir,
                                       config_csv=config_csv,
                                       prefix=prefix)
        # Update cellranger reference data in qc.info
        qc_info = self._project.qc_info(self._project.qc_dir)
        qc_info['cellranger_refdata'] = reference_data_path
        qc_info.save()
        self._reload_project()

#######################################################################
# Factory classes for quickly generating mock directories
#######################################################################

class MockAnalysisDirFactory:
    """
    Collection of convenient pre-populated test cases
    """
    @classmethod
    def bcl2fastq2(self,
                   run_name,platform,
                   paired_end=True,
                   unaligned_dir='bcl2fastq',
                   no_lane_splitting=True,
                   reads=None,
                   top_dir=None,
                   params=None,
                   metadata=None,
                   project_metadata=None,
                   bases_mask='auto',
                   include_stats_files=False):
        """
        Basic analysis dir from bcl2fastq v2
        """
        if no_lane_splitting:
            lanes = []
        else:
            lanes = [1,2,3,4]
        mad = MockAnalysisDir(run_name,platform,
                              unaligned_dir=unaligned_dir,
                              fmt='bcl2fastq2',
                              bases_mask=bases_mask,
                              paired_end=paired_end,
                              no_lane_splitting=no_lane_splitting,
                              lanes=lanes,
                              params=params,
                              metadata=metadata,
                              project_metadata=project_metadata,
                              include_stats_files=include_stats_files,
                              top_dir=top_dir)
        mad.add_fastq_batch('AB','AB1','AB1_S1',lanes=lanes,reads=reads)
        mad.add_fastq_batch('AB','AB2','AB2_S2',lanes=lanes,reads=reads)
        mad.add_fastq_batch('CDE','CDE3','CDE3_S3',lanes=lanes,reads=reads)
        mad.add_fastq_batch('CDE','CDE4','CDE4_S4',lanes=lanes,reads=reads)
        return mad

    @classmethod
    def casava(self,
               run_name,platform,
               paired_end=True,
               unaligned_dir='bcl2fastq',
               params=None,
               metadata=None,
               top_dir=None):
        """
        Basic analysis dir from CASAVA/bcl2fastq v1.8
        """
        lanes = [1,2,3,4]
        mad = MockAnalysisDir(run_name,platform,
                              unaligned_dir=unaligned_dir,
                              fmt='casava',
                              paired_end=paired_end,
                              lanes=lanes,
                              params=params,
                              metadata=metadata,
                              top_dir=top_dir)
        mad.add_fastq_batch('AB','AB1','AB1_GCCAAT',lanes=lanes)
        mad.add_fastq_batch('AB','AB2','AB2_AGTCAA',lanes=lanes)
        mad.add_fastq_batch('CDE','CDE3','CDE3_GCCAAT',lanes=lanes)
        mad.add_fastq_batch('CDE','CDE4','CDE4_AGTCAA',lanes=lanes)
        return mad

#######################################################################
# Classes for making mock executables
#######################################################################

class MockBcl2fastq2Exe:
    """
    Create mock bcl2fastq2 executable

    This class can be used to create a mock bcl2fastq
    executable, which in turn can be used in place of
    the actual bcl2fastq software for testing purposes.

    To create a mock executable, use the 'create' static
    method, e.g.

    >>> MockBcl2fastq2Exe.create("/tmpbin/bcl2fastq")

    The resulting executable will generate mock outputs
    when run on actual or mock Illumina sequencer output
    directories (mock versions can be produced using the
    'mock.IlluminaRun' class in the genomics-bcftbx
    package).

    The executable can be configured on creation to
    produce different error conditions when run:

    - the exit code can be set to an arbitrary value
      via the `exit_code` argument
    - Fastqs can be removed from the output by
      specifying their names in the `missing_fastqs`
      argument

    The executable can also be configured to check
    supplied values:

    - the bases mask can be checked via the
      `assert_bases_mask` argument
    - lane splitting can be checked via the
      `assert_no_lane_splitting` argument
    - creation of Fastqs for index reads can
      be checked via the
      `assert_create_fastq_for_index_read`
      argument
    - adapter trimming and masking can be
      checked via the
      `assert_minimum_trimmed_read_length`
      and `assert_mask_short_adapter_reads`
      arguments
    - adapter sequences can be checked via
      the `assert_adapter` and
      `assert_adapter2` arguments
    - sliding window algorith for adapter
      trimming can be checked via
      `assert_find_adapters_with_sliding_window`
    """

    @staticmethod
    def create(path,exit_code=0,missing_fastqs=None,
               platform=None,assert_bases_mask=None,
               assert_no_lane_splitting=None,
               assert_create_fastq_for_index_read=None,
               assert_minimum_trimmed_read_length=None,
               assert_mask_short_adapter_reads=None,
               assert_adapter=None,assert_adapter2=None,
               assert_find_adapters_with_sliding_window=None,
               version='2.20.0.422'):
        """
        Create a "mock" bcl2fastq executable

        Arguments:
          path (str): path to the new executable
            to create. The final executable must
            not exist, however the directory it
            will be created in must.
          exit_code (int): exit code that the
            mock executable should complete
            with
          missing_fastqs (list): list of Fastq
            names that will not be created
          platform (str): platform for primary
            data (if it cannot be determined from
            the directory/instrument name)
          assert_bases_mask (str): if set then
            assert that bases mask matches the
            supplied string
          assert_lane_splitting (bool): if set
            then assert that --no-lane-splitting
            matches the supplied boolean value
          assert_create_fastq_for_index_read
            (bool): if set then assert that
            --create-fastq-for-index-read
            matches the supplied boolean value
          assert_minimum_trimmed_read_length (int):
            if set then assert that
            --minimum-trimmed-read-length matches
            the supplied value
          assert_mask_short_adapter_reads (int):
            if set then assert that
            --mask-short-adapter-reads matches
            the supplied value
          assert_adapter (str): if set then
            assert that the adapter sequence
            in the sample sheet matches the
            supplied value
          assert_adapter2 (str): if set then
            assert that the adapter sequence
            for read2 in the sample sheet matches
            the supplied value
          assert_find_adapters_with_sliding_window
            (bool): if set then assert that
            --find-adapters-with-sliding-window
            matches the supplied boolean value
          version (str): version of bcl2fastq2
            to imitate
        """
        path = os.path.abspath(path)
        print("Building mock executable: %s" % path)
        # Don't clobber an existing executable
        assert(os.path.exists(path) is False)
        with open(path,'w') as fp:
            fp.write("""#!/usr/bin/env python
import sys
from auto_process_ngs.mock import MockBcl2fastq2Exe
sys.exit(MockBcl2fastq2Exe(exit_code=%s,
                           missing_fastqs=%s,
                           platform=%s,
                           assert_bases_mask=%s,
                           assert_no_lane_splitting=%s,
                           assert_create_fastq_for_index_read=%s,
                           assert_minimum_trimmed_read_length=%s,
                           assert_mask_short_adapter_reads=%s,
                           assert_adapter=%s,
                           assert_adapter2=%s,
                           assert_find_adapters_with_sliding_window=%s,
                           version=%s).main(sys.argv[1:]))
            """ % (exit_code,
                   missing_fastqs,
                   ("\"%s\"" % platform
                    if platform is not None
                    else None),
                   ("\"%s\"" % assert_bases_mask
                    if assert_bases_mask is not None
                    else None),
                   assert_no_lane_splitting,
                   assert_create_fastq_for_index_read,
                   assert_minimum_trimmed_read_length,
                   assert_mask_short_adapter_reads,
                   ("\"%s\"" % assert_adapter
                    if assert_adapter is not None
                    else None),
                   ("\"%s\"" % assert_adapter2
                    if assert_adapter2 is not None
                    else None),
                   assert_find_adapters_with_sliding_window,
                   ("\"%s\"" % version
                    if version is not None
                    else None)))
            os.chmod(path,0o775)
        with open(path,'r') as fp:
            print("bcl2fastq:")
            print("%s" % fp.read())
        return path

    def __init__(self,exit_code=0,
                 missing_fastqs=None,
                 platform=None,
                 assert_bases_mask=None,
                 assert_no_lane_splitting=None,
                 assert_create_fastq_for_index_read=None,
                 assert_minimum_trimmed_read_length=None,
                 assert_mask_short_adapter_reads=None,
                 assert_adapter=None,
                 assert_adapter2=None,
                 assert_find_adapters_with_sliding_window=None,
                 version=None):
        """
        Internal: configure the mock bcl2fastq2
        """
        self._exit_code = exit_code
        self._missing_fastqs = missing_fastqs
        self._platform = platform
        self._assert_bases_mask = assert_bases_mask
        self._assert_no_lane_splitting = assert_no_lane_splitting
        self._assert_create_fastq_for_index_read = \
                                assert_create_fastq_for_index_read
        self._assert_minimum_trimmed_read_length = \
                                assert_minimum_trimmed_read_length
        self._assert_mask_short_adapter_reads = \
                                assert_mask_short_adapter_reads
        self._assert_adapter = assert_adapter
        self._assert_adapter2 = assert_adapter2
        self._assert_find_adapters_with_sliding_window = \
                                assert_find_adapters_with_sliding_window
        self._version = version

    def main(self,args):
        """
        Internal: provides mock bcl2fastq2 functionality
        """
        # Build generic header
        header = """BCL to FASTQ file converter
bcl2fastq v%s
""" % self._version
        if self._version.startswith("2.17."):
            header += \
            "Copyright (c) 2007-2015 Illumina, Inc.\n\n2015-12-17 14:08:00 [7fa113f3f780] Command-line invocation: bcl2fastq %s" % ' '.join(args)
        elif self._version.startswith("2.20."):
            header += "Copyright (c) 2007-2017 Illumina, Inc.\n"
        # Handle version request
        if "--version" in args:
            print(header)
            return self._exit_code
        # Deal with arguments
        p = argparse.ArgumentParser()
        p.add_argument("--runfolder-dir",action="store")
        p.add_argument("--output-dir",action="store")
        p.add_argument("--sample-sheet",action="store")
        p.add_argument("--use-bases-mask",action="store")
        p.add_argument("--barcode-mismatches",action="store")
        p.add_argument("--minimum-trimmed-read-length",action="store")
        p.add_argument("--mask-short-adapter-reads",action="store")
        p.add_argument("--ignore-missing-bcls",action="store_true")
        p.add_argument("--no-lane-splitting",action="store_true")
        p.add_argument("--create-fastq-for-index-read",action="store_true")
        p.add_argument("--find-adapters-with-sliding-window",
                       action="store_true")
        p.add_argument("-r",action="store")
        if self._version.startswith("2.17."):
            p.add_argument("-d",action="store")
        p.add_argument("-p",action="store")
        p.add_argument("-w",action="store")
        args = p.parse_args(args)
        # Check bases mask
        if self._assert_bases_mask:
            print("Checking bases mask: %s" % args.use_bases_mask)
            assert(args.use_bases_mask == self._assert_bases_mask)
        # Check --no-lane-splitting
        if self._assert_no_lane_splitting is not None:
            print("Checking --no-lane-splitting: %s" % args.no_lane_splitting)
            assert(args.no_lane_splitting == self._assert_no_lane_splitting)
        # Check --create-fastq-for-index-read
        if self._assert_create_fastq_for_index_read is not None:
            print("Checking --create-fastq-for-index-read: %s" %
                  args.create_fastq_for_index_read)
            assert(args.create_fastq_for_index_read ==
                   self._assert_create_fastq_for_index_read)
        # Check --minimum-trimmed-read-length
        if self._assert_minimum_trimmed_read_length is not None:
            print("Checking --minimum-trimmed-read-length: %s" %
                  args.minimum_trimmed_read_length)
            assert(int(args.minimum_trimmed_read_length) ==
                   int(self._assert_minimum_trimmed_read_length))
        # Check --mask-short-adapter-reads
        if self._assert_mask_short_adapter_reads is not None:
            print("Checking --mask-short-adapter-reads: %s" %
                  args.mask_short_adapter_reads)
            assert(int(args.mask_short_adapter_reads) ==
                   int(self._assert_mask_short_adapter_reads))
        # Check --find-adapters-with-sliding-window
        if self._assert_find_adapters_with_sliding_window is not None:
            print("Checking --find-adapters-with-sliding-window: %s" %
                  args.find_adapters_with_sliding_window)
            assert(args.find_adapters_with_sliding_window ==
                   self._assert_find_adapters_with_sliding_window)
        # Platform
        print("Platform (default): %s" % self._platform)
        # Run folder (input data)
        runfolder = args.runfolder_dir
        print("Runfolder dir: %s" % runfolder)
        if runfolder is None:
            return 1
        run_info_xml = os.path.join(runfolder,"RunInfo.xml")
        if not os.path.exists(run_info_xml):
            return 1
        # Modifiers
        no_lane_splitting = bool(args.no_lane_splitting)
        print("No lane splitting: %s" % no_lane_splitting)
        create_fastq_for_index_read = bool(args.create_fastq_for_index_read)
        print("Create fastq for index read: %s" %
              create_fastq_for_index_read)
        # Get reads and determine if run is paired end
        reads = []
        ii = 1
        ir = 1
        for r in IlluminaRunInfo(run_info_xml).reads:
            if r['is_indexed_read'] == 'N':
                reads.append("R%d" % ir)
                ir += 1
            elif create_fastq_for_index_read:
                reads.append("I%d" % ii)
                ii += 1
        if ir == 2:
            paired_end = True
        else:
            paired_end = False
        print("Paired-end: %s" % paired_end)
        # Lanes
        lanes = IlluminaRun(runfolder,platform=self._platform).lanes
        print("Lanes: %s" % lanes)
        # Output folder
        output_dir = args.output_dir
        if output_dir is None:
            output_dir = "bcl2fastq"
        print("Output dir: %s" % output_dir)
        # Sample sheet
        sample_sheet = args.sample_sheet
        if sample_sheet is None:
            for d in (runfolder,
                      os.path.join(runfolder,
                                   "Data",
                                   "Intensities",
                                   "BaseCalls")):
                sample_sheet = os.path.join(d,"SampleSheet.csv")
                if os.path.exists(sample_sheet):
                    break
                sample_sheet = None
        print("Sample sheet: %s" % sample_sheet)
        # Check the adapter sequences from the sample sheet
        if (self._assert_adapter is not None) or \
           (self._assert_adapter2 is not None):
            if sample_sheet is None:
                # No sample sheet
                adapter = None
                adapter2 = None
            else:
                # Extract adapter sequences
                s = SampleSheet(sample_sheet)
                try:
                    adapter = s.settings['Adapter']
                except KeyError:
                    adapter = ""
                try:
                    adapter2 = s.settings['AdapterRead2']
                except KeyError:
                    adapter2 = ""
            if self._assert_adapter is not None:
                assert(self._assert_adapter == adapter)
            if self._assert_adapter2 is not None:
                assert(self._assert_adapter2 == adapter2)
        # Generate mock output based on inputs
        tmpname = "tmp.%s" % uuid.uuid4()
        make_mock_bcl2fastq2_output(os.path.join(tmpname,"bcl2fastq"),
                                    lanes=lanes,
                                    reads=reads,
                                    sample_sheet=sample_sheet,
                                    no_lane_splitting=no_lane_splitting,
                                    exclude_fastqs=self._missing_fastqs)
        # Move to final location
        os.rename(os.path.join(tmpname,"bcl2fastq"),
                  output_dir)
        shutil.rmtree(tmpname)
        return self._exit_code

class MockBclConvertExe:
    """
    Create mock bcl-convert executable

    This class can be used to create a mock bcl-convert
    executable, which in turn can be used in place of
    the actual BCLConvert software for testing purposes.

    To create a mock executable, use the 'create' static
    method, e.g.

    >>> MockBclConvertExe.create("/tmpbin/bcl-convert")

    The resulting executable will generate mock outputs
    when run on actual or mock Illumina sequencer output
    directories (mock versions can be produced using the
    'mock.IlluminaRun' class in the genomics-bcftbx
    package).

    The executable can be configured on creation to
    produce different error conditions when run:

    - the exit code can be set to an arbitrary value
      via the `exit_code` argument
    - Fastqs can be removed from the output by
      specifying their names in the `missing_fastqs`
      argument

    The executable can also be configured to check
    supplied values:

    - the masking of cycles can be checked via the
      `assert_override_cycles` argument
    - lane splitting can be checked via the
      `assert_no_lane_splitting` argument
    - creation of Fastqs for index reads can
      be checked via the
      `assert_create_fastq_for_index_read`
      argument
    - adapter trimming and masking can be
      checked via the
      `assert_minimum_trimmed_read_length`
      and `assert_mask_short_reads`
      arguments
    - adapater sequences can be checked via
      the `assert_adapter1` and
      `assert_adapter2` arguments
    """

    @staticmethod
    def create(path,exit_code=0,missing_fastqs=None,
               platform=None,assert_override_cycles=None,
               assert_no_lane_splitting=None,
               assert_create_fastq_for_index_read=None,
               assert_minimum_trimmed_read_length=None,
               assert_mask_short_reads=None,
               assert_adapter1=None,assert_adapter2=None,
               version='3.7.5'):
        """
        Create a "mock" bcl-convert executable

        Arguments:
          path (str): path to the new executable
            to create. The final executable must
            not exist, however the directory it
            will be created in must.
          exit_code (int): exit code that the
            mock executable should complete
            with
          missing_fastqs (list): list of Fastq
            names that will not be created
          platform (str): platform for primary
            data (if it cannot be determined from
            the directory/instrument name)
          assert_override_cycles (str): if set then
            assert that the 'OverrideCycles' setting
            matches the supplied string
          assert_lane_splitting (bool): if set
            then assert that --no-lane-splitting
            matches the supplied boolean value
          assert_create_fastq_for_index_read
            (bool): if set then assert that
            --create-fastq-for-index-read
            matches the supplied boolean value
          assert_minimum_trimmed_read_length (int):
            if set then assert that
            --minimum-trimmed-read-length matches
            the supplied value
          assert_mask_short_reads (int):
            if set then assert that
            --mask-short-adapter-reads matches
            the supplied value
          assert_adapter1 (str): if set then
            assert that the adapter sequence
            in the sample sheet matches the
            supplied value
          assert_adapter2 (str): if set then
            assert that the adapter sequence
            for read2 in the sample sheet matches
            the supplied value
          version (str): version of BCLConvert
            to imitate
        """
        path = os.path.abspath(path)
        print("Building mock executable: %s" % path)
        # Don't clobber an existing executable
        assert(os.path.exists(path) is False)
        with open(path,'w') as fp:
            fp.write("""#!/usr/bin/env python
import sys
from auto_process_ngs.mock import MockBclConvertExe
sys.exit(MockBclConvertExe(exit_code=%s,
                           missing_fastqs=%s,
                           platform=%s,
                           assert_override_cycles=%s,
                           assert_no_lane_splitting=%s,
                           assert_create_fastq_for_index_read=%s,
                           assert_minimum_trimmed_read_length=%s,
                           assert_mask_short_reads=%s,
                           assert_adapter1=%s,
                           assert_adapter2=%s,
                           version=%s).main(sys.argv[1:]))
            """ % (exit_code,
                   missing_fastqs,
                   ("\"%s\"" % platform
                    if platform is not None
                    else None),
                   ("\"%s\"" % assert_override_cycles
                    if assert_override_cycles is not None
                    else None),
                   assert_no_lane_splitting,
                   assert_create_fastq_for_index_read,
                   assert_minimum_trimmed_read_length,
                   assert_mask_short_reads,
                   ("\"%s\"" % assert_adapter1
                    if assert_adapter1 is not None
                    else None),
                   ("\"%s\"" % assert_adapter2
                    if assert_adapter2 is not None
                    else None),
                   ("\"%s\"" % version
                    if version is not None
                    else None)))
            os.chmod(path,0o775)
        with open(path,'r') as fp:
            print("bcl2fastq:")
            print("%s" % fp.read())
        return path

    def __init__(self,exit_code=0,
                 missing_fastqs=None,
                 platform=None,
                 assert_override_cycles=None,
                 assert_no_lane_splitting=None,
                 assert_create_fastq_for_index_read=None,
                 assert_minimum_trimmed_read_length=None,
                 assert_mask_short_reads=None,
                 assert_adapter1=None,
                 assert_adapter2=None,
                 version=None):
        """
        Internal: configure the mock bcl-convert program
        """
        self._exit_code = exit_code
        self._missing_fastqs = missing_fastqs
        self._platform = platform
        self._assert_override_cycles = assert_override_cycles
        self._assert_no_lane_splitting = assert_no_lane_splitting
        self._assert_create_fastq_for_index_read = \
                                assert_create_fastq_for_index_read
        self._assert_minimum_trimmed_read_length = \
                                assert_minimum_trimmed_read_length
        self._assert_mask_short_reads = assert_mask_short_reads
        self._assert_adapter1 = assert_adapter1
        self._assert_adapter2 = assert_adapter2
        self._version = version

    def main(self,args):
        """
        Internal: provides mock bcl-convert functionality
        """
        # Build generic header
        header = """bcl-convert Version 00.000.000.%s
Copyright (c) 2014-2018 Illumina, Inc.
""" % self._version
        # Handle version request
        if "-V" in args:
            print(header)
            return self._exit_code
        # Deal with arguments
        p = argparse.ArgumentParser()
        p.add_argument("-r","--bcl-input-directory",action="store")
        p.add_argument("-o","--output-dir",action="store")
        p.add_argument("--sample-sheet",action="store")
        p.add_argument("--bcl-only-lane",action="store")
        p.add_argument("--mask-short-adapter-reads",action="store")
        p.add_argument("--no-lane-splitting",action="store")
        p.add_argument("--bcl-sampleproject-subdirectories",action="store")
        p.add_argument("--bcl-num-parallel-tiles",action="store")
        p.add_argument("--bcl-num-conversion-threads",action="store")
        p.add_argument("--bcl-num-compression-threads",action="store")
        p.add_argument("--bcl-num-decompression-threads",action="store")
        args = p.parse_args(args)
        # Check for input directory
        input_dir = args.bcl_input_directory
        print("BCL input directory: %s" % input_dir)
        if not input_dir:
            return 1
        # Output directory
        output_dir = args.output_dir
        print("Output directory: %s" % output_dir)
        if not output_dir:
            return 1
        elif os.path.exists(output_dir):
            return 1
        # Lane specified on command line
        if args.bcl_only_lane:
            try:
                lanes = [int(args.bcl_only_lane)]
            except ValueError:
                return 1
        else:
            lanes = []
        # Check for sample sheet file
        if not args.sample_sheet:
            # Default location is the root input folder
            sample_sheet_file = os.path.join(input_dir,"SampleSheet.csv")
        else:
            sample_sheet_file = args.sample_sheet
        print("Sample sheet: %s" % sample_sheet_file)
        if not os.path.exists(sample_sheet_file):
            return 1
        # Extract default settings from the sample sheet
        sample_sheet = SampleSheet(sample_sheet_file)
        create_fastq_for_index_read = 0
        barcode_mismatch_index1 = 1
        barcode_mismatch_index2 = 1
        adapter_sequence_read1 = None
        adapter_sequence_read2 = None
        minimum_trimmed_read_length = 35
        mask_short_reads = 22
        override_cycles = None
        for item in sample_sheet.settings_items:
            if item == "CreateFastqForIndexRead":
                # Create index read Fastqs
                create_fastq_for_index_read = \
                    sample_sheet.settings['CreateFastqForIndexRead']
            elif item.startswith("BarcodeMismatchIndex"):
                # Number of mismatches allowed
                indx = int(item[-1])
                if indx == 1:
                    barcode_mismatch_index1 = \
                        sample_sheet.settings['BarcodeMismatchIndex1']
                elif indx == 2:
                    barcode_mismatch_index2 = \
                        sample_sheet.settings['BarcodeMismatchIndex2']
            elif item.startswith("AdapterRead"):
                # Adapter sequences
                indx = int(item[-1])
                if indx == 1:
                    adapter_sequence_read1 = \
                        sample_sheet.settings['AdapterRead1']
                elif indx == 2:
                    adapter_sequence_read2 = \
                        sample_sheet.settings['AdapterRead2']
            elif item == "MinimumTrimmedReadLength":
                # Minimum sequence length after trimming
                minimum_trimmed_read_length = \
                    sample_sheet.settings['MinimumTrimmedReadLength']
            elif item == "MaskShortReads":
                # Masking short reads
                mask_short_reads = sample_sheet.settings['MaskShortReads']
            elif item == "OverrideCycles":
                override_cycles = sample_sheet.settings['OverrideCycles']
        if sample_sheet.has_lanes and lanes:
            # Check lane specifications match up between
            # command line and sample sheet
            sample_sheet_lanes = set()
            for line in sample_sheet:
                sample_sheet_lanes.add(int(line['Lane']))
            sample_sheet_lanes = sorted(list(lanes))
            if sample_sheet_lanes != lanes:
                return 1
        # Check the adapter sequences from the sample sheet
        if (self._assert_adapter1 is not None) or \
           (self._assert_adapter2 is not None):
            # Extract adapter sequences
            try:
                adapter1 = sample_sheet.settings['AdapterRead1']
            except KeyError:
                adapter1 = ""
            try:
                adapter2 = sample_sheet.settings['AdapterRead2']
            except KeyError:
                adapter2 = ""
            if self._assert_adapter1 is not None:
                print("Check AdapterRead1: should be '%s', got '%s'" %
                      (self._assert_adapter1,adapter1))
                assert(self._assert_adapter1 == adapter1)
            if self._assert_adapter2 is not None:
                print("Check AdapterRead2: should be '%s', got '%s'" %
                      (self._assert_adapter2,adapter2))
                assert(self._assert_adapter2 == adapter2)
        # No lane splitting
        no_lane_splitting = False
        if args.no_lane_splitting:
            if args.no_lane_splitting == "true":
                no_lane_splitting = True
        # SampleProject subdirectories
        sampleproject_subdirectories = False
        if args.bcl_sampleproject_subdirectories:
            if args.bcl_sampleproject_subdirectories == "true":
                sampleproject_subdirectories = True
        # Check override cycles
        if self._assert_override_cycles:
            print("Checking OverrideCycles: %s" % override_cycles)
            assert(override_cycles == self._assert_override_cycles)
        # Check CreateFastqForIndexRead
        if self._assert_create_fastq_for_index_read is not None:
            print("Checking CreateFastqForIndexRead: %s" %
                  create_fastq_for_index_read)
            assert(create_fastq_for_index_read ==
                   self._assert_create_fastq_for_index_read)
        # Check MinimumTrimmedReadLength
        if self._assert_minimum_trimmed_read_length is not None:
            print("Checking MinimumTrimmedReadLength: %s" %
                  minimum_trimmed_read_length)
            assert(int(minimum_trimmed_read_length) ==
                   int(self._assert_minimum_trimmed_read_length))
        # Check MaskShortReads
        if self._assert_mask_short_reads is not None:
            print("Checking MaskShortReads: %s" % mask_short_reads)
            assert(int(mask_short_reads) ==
                   int(self._assert_mask_short_reads))
        # Check --no-lane-splitting
        if self._assert_no_lane_splitting is not None:
            print("Checking --no-lane-splitting: %s" % no_lane_splitting)
            assert(no_lane_splitting == self._assert_no_lane_splitting)
        # Platform
        print("Platform (default): %s" % self._platform)
        # RunInfo.xml
        run_info_xml = os.path.join(input_dir,"RunInfo.xml")
        if not os.path.exists(run_info_xml):
            return 1
        # Modifiers
        no_lane_splitting = bool(args.no_lane_splitting)
        print("No lane splitting: %s" % no_lane_splitting)
        create_fastq_for_index_read = bool(create_fastq_for_index_read)
        print("Create fastq for index read: %s" %
              create_fastq_for_index_read)
        # Get reads and determine if run is paired end
        reads = []
        ii = 1
        ir = 1
        for r in IlluminaRunInfo(run_info_xml).reads:
            if r['is_indexed_read'] == 'N':
                reads.append("R%d" % ir)
                ir += 1
            elif create_fastq_for_index_read:
                reads.append("I%d" % ii)
                ii += 1
        if ir == 2:
            paired_end = True
        else:
            paired_end = False
        print("Paired-end: %s" % paired_end)
        # Lanes
        if not lanes:
            # Include all lanes
            lanes = IlluminaRun(input_dir,platform=self._platform).lanes
        print("Lanes: %s" % lanes)
        # Generate mock output based on inputs
        tmpname = "tmp.%s" % uuid.uuid4()
        make_mock_bcl2fastq2_output(os.path.join(tmpname,"bclconvert"),
                                    lanes=lanes,
                                    reads=reads,
                                    sample_sheet=sample_sheet_file,
                                    no_lane_splitting=no_lane_splitting,
                                    exclude_fastqs=self._missing_fastqs)
        # Move to final location
        os.rename(os.path.join(tmpname,"bclconvert"),
                  output_dir)
        shutil.rmtree(tmpname)
        return self._exit_code

class Mock10xPackageExe:
    """
    Create mock 10xGenomics pipeline executable

    This class can be used to create a mock 10xGenomics
    pipeline executable, which in turn can be used in
    place of the actual pipeline software (e.g. cellranger)
    for testing purposes.

    To create a mock executable, use the 'create' static
    method, e.g.

    >>> Mock10xPackageExe.create("/tmpbin/cellranger")

    The resulting executable will generate mock outputs
    when run on actual or mock Illumina sequencer output
    directories (mock versions can be produced using the
    'mock.IlluminaRun' class in the genomics-bcftbx
    package).

    The executable can be configured on creation to
    produce different error conditions when run:

    - the exit code can be set to an arbitrary value
      via the `exit_code` argument
    """

    @staticmethod
    def create(path,exit_code=0,missing_fastqs=None,
               platform=None,assert_bases_mask=None,
               assert_include_introns=None,
               assert_chemistry=None,
               assert_force_cells=None,
               assert_filter_single_index=None,
               assert_filter_dual_index=None,
               reads=None,multiome_data=None,
               version=None):
        """
        Create a "mock" 10xGenomics package executable

        Arguments:
          path (str): path to the new executable
            to create. The final executable must
            not exist, however the directory it
            will be created in must.
          exit_code (int): exit code that the
            mock executable should complete
            with
          missing_fastqs (list): list of Fastq
            names that will not be created
          platform (str): platform for primary
            data (if it cannot be determined from
            the directory/instrument name)
          assert_bases_mask (str): if set then
            check that the supplied bases mask
            matches this value
          assert_include_introns (bool): if set
            to True/False then check that the
            '--include-introns' option was/n't
            set (ignored if set to None)
          assert_chemistry (str): if set then
            check that the supplied chemistry
            specification matches this value
          assert_force_cells (int): if set then
            check that the '--force-cells' option
            was specified with this value
          assert_filter_single_index (bool): if
            set to True/False then check that
            the '--filter-single-index' option
            was/n't supplied (ignored if set to
            None)
          assert_filter_dual_index (bool): if
            set to True/False then check that
            the '--filter-dual-index' option
            was/n't supplied (ignored if set to
            None)
          reads (list): list of 'reads' that
            will be created
          multiome_data (str): either 'GEX' or
            'ATAC' (when mocking 'cellranger-arc')
          version (str): version of package to
            report
        """
        path = os.path.abspath(path)
        print("Building mock executable: %s" % path)
        # Don't clobber an existing executable
        assert(os.path.exists(path) is False)
        with open(path,'w') as fp:
            fp.write("""#!/usr/bin/env python
import sys
from auto_process_ngs.mock import Mock10xPackageExe
sys.exit(Mock10xPackageExe(path=sys.argv[0],
                           exit_code=%s,
                           platform=%s,
                           assert_bases_mask=%s,
                           assert_include_introns=%s,
                           assert_chemistry=%s,
                           assert_force_cells=%s,
                           assert_filter_single_index=%s,
                           assert_filter_dual_index=%s,
                           reads=%s,
                           multiome_data=%s,
                           version=%s).main(sys.argv[1:]))
            """ % (exit_code,
                   ("\"%s\"" % platform
                    if platform is not None
                    else None),
                   ("\"%s\"" % assert_bases_mask
                    if assert_bases_mask is not None
                    else None),
                   assert_include_introns,
                   ("\"%s\"" % assert_chemistry
                    if assert_chemistry is not None
                    else None),
                   assert_force_cells,
                   assert_filter_single_index,
                   assert_filter_dual_index,
                   reads,
                   ("\"%s\"" % multiome_data
                    if multiome_data is not None
                    else None),
                   ("\"%s\"" % version
                    if version is not None
                    else None)))
            os.chmod(path,0o775)
        with open(path,'r') as fp:
            print("%s:" % os.path.basename(path))
            print("%s" % fp.read())
        return path

    def __init__(self,path,
                 exit_code=0,
                 platform=None,
                 assert_bases_mask=None,
                 assert_include_introns=None,
                 assert_chemistry=None,
                 assert_force_cells=None,
                 assert_filter_single_index=None,
                 assert_filter_dual_index=None,
                 reads=None,
                 multiome_data=None,
                 version=None):
        """
        Internal: configure the mock 10xGenomics package
        """
        self._path = path
        self._package_name = os.path.basename(self._path)
        if version is None:
            if self._package_name == 'cellranger':
                self._version = '7.0.0'
            elif self._package_name == 'cellranger-atac':
                self._version = '2.0.0'
            elif self._package_name == 'cellranger-arc':
                self._version = '2.0.0'
            elif self._package_name == 'spaceranger':
                self._version = '1.3.1'
        else:
            self._version = version
        self._exit_code = exit_code
        self._platform = platform
        self._assert_bases_mask = assert_bases_mask
        self._assert_include_introns = assert_include_introns
        self._assert_chemistry = assert_chemistry
        self._assert_force_cells = assert_force_cells
        self._assert_filter_single_index = assert_filter_single_index
        self._assert_filter_dual_index = assert_filter_dual_index
        self._multiome_data = str(multiome_data).upper()
        if self._package_name == 'cellranger-arc':
            assert self._multiome_data is not None
        if reads is None:
            if self._package_name in ('cellranger',
                                      'spaceranger',):
                self._reads = ('R1','R2','I1',)
            elif self._package_name == 'cellranger-atac':
                self._reads = ('R1','R2','R3','I1',)
            elif self._package_name == 'cellranger-arc':
                if self._multiome_data == 'GEX':
                    self._reads = ('R1','R2','I1',)
                elif self._multiome_data == 'ATAC':
                    self._reads = ('R1','R2','R3','I1',)
        else:
            self._reads = tuple(reads)

    def _write_qc_summary_json(self,path,run_dir,out_dir,
                               sample_sheet,samples,lanes):
        """
        Internal: write qc_summary.json file
        """
        with open(path,'w') as fp:
            fp.write("""{
    "10x_software_version": "cellranger %s", 
    "PhiX_aligned": null, 
    "PhiX_error_worst_tile": null, 
""" % self._version)
            fp.write("""    "bcl2fastq_args": "bcl2fastq --minimum-trimmed-read-length 8 --mask-short-adapter-reads 8 --create-fastq-for-index-reads --ignore-missing-positions --ignore-missing-filter --ignore-missing-bcls --use-bases-mask=Y76,I8,Y76 -R %s --output-dir=%s --interop-dir=%s/HJCY7BBXX_56/MAKE_FASTQS_CS/MAKE_FASTQS/BCL2FASTQ_WITH_SAMPLESHEET/fork0/chnk0/files/interop_path --sample-sheet=%s/HJCY7BBXX_56/MAKE_FASTQS_CS/MAKE_FASTQS/PREPARE_SAMPLESHEET/fork0/files/samplesheet.csv --tiles s_[56] -p 6 -d 6 -r 6 -w 6", 
""" % (run_dir,
       out_dir,
       os.getcwd(),
       os.getcwd()))
            fp.write("""    "bcl2fastq_version": "2.17.1.14", 
    "experiment_name": "Example run", 
    "fwhm_A": null, 
    "fwhm_C": null, 
    "fwhm_G": null, 
    "fwhm_T": null, 
    "index1_PhiX_error_by_cycle": null, 
    "index1_mean_phasing": null, 
    "index1_mean_prephasing": null, 
    "index1_q20_fraction": null, 
    "index1_q20_fraction_by_cycle": null, 
    "index1_q30_fraction": null, 
    "index1_q30_fraction_by_cycle": null, 
    "index2_PhiX_error_by_cycle": null, 
    "index2_mean_phasing": null, 
    "index2_mean_prephasing": null, 
    "index2_q20_fraction": null, 
    "index2_q20_fraction_by_cycle": null, 
    "index2_q30_fraction": null, 
    "index2_q30_fraction_by_cycle": null, 
    "intensity_A": null, 
    "intensity_C": null, 
    "intensity_G": null, 
    "intensity_T": null, 
    "lanecount": 8, 
    "mean_cluster_density": null, 
    "mean_cluster_density_pf": null, 
    "num_clusters": null, 
    "percent_pf_clusters": null, 
    "read1_PhiX_error_by_cycle": null, 
    "read1_mean_phasing": null, 
    "read1_mean_prephasing": null, 
    "read1_q20_fraction": null, 
    "read1_q20_fraction_by_cycle": null, 
    "read1_q30_fraction": null, 
    "read1_q30_fraction_by_cycle": null, 
    "read2_PhiX_error_by_cycle": null, 
    "read2_mean_phasing": null, 
    "read2_mean_prephasing": null, 
    "read2_q20_fraction": null, 
    "read2_q20_fraction_by_cycle": null, 
    "read2_q30_fraction": null, 
    "read2_q30_fraction_by_cycle": null, 
    "rta_version": "2.7.6", 
    "run_id": "170426_K00311_0033_AHJCY7BBXX", 
    "sample_qc": {
""")
            for sample in samples:
                fp.write("""        "%s": {\n""" % sample)
                for lane in lanes:
                    fp.write("""            "%s": {
                "barcode_exact_match_ratio": 0.9716234673603158, 
                "barcode_q30_base_ratio": 0.9857481459248205, 
                "bc_on_whitelist": 0.9814401918935765, 
                "gem_count_estimate": 84073, 
                "mean_barcode_qscore": 38.592230237502555, 
                "number_reads": 22220650, 
                "read1_q30_base_ratio": 0.4826826169702818, 
                "read2_q30_base_ratio": 0.8117288414214294
            }, 
""" % lane)
                fp.write("""            "all": {
                "barcode_exact_match_ratio": 0.9693936169422027, 
                "barcode_q30_base_ratio": 0.9842256841986534, 
                "bc_on_whitelist": 0.9803349567432192, 
                "gem_count_estimate": 88278, 
                "mean_barcode_qscore": 38.547743820121745, 
                "number_reads": 42378905, 
                "read1_q30_base_ratio": 0.4825840484250696, 
                "read2_q30_base_ratio": 0.8140753810050518
            }
        }%s 
""" % ((',' if sample != samples[-1] else '',)))
            fp.write("""    },
    "signoise_ratio": null, 
    "start_datetime": "2017-04-26 00:00:00", 
    "surfacecount": 2, 
    "swathcount": 2, 
    "tilecount": 28, 
    "total_cluster_density": null, 
    "total_cluster_density_pf": null, 
    "yield": null, 
    "yield_pf": null, 
    "yield_pf_q30": null
}
""")

    def main(self,args):
        """
        Internal: provides mock 10xGenomics package functionality
        """
        # Store command line
        cmdline = "%s %s" % (self._path,' '.join(args))
        # Split version
        version = tuple([int(x) for x in self._version.split('.')])
        # Build generic header
        try:
            cmd = " %s" % args[0]
        except IndexError:
            cmd = ''
        # Construct header
        if (self._package_name == 'cellranger' and version[0] <= 3) or \
           self._package_name == 'cellranger-atac':
            header = """%s
%s%s (%s)
Copyright (c) 2018 10x Genomics, Inc.  All rights reserved.
-------------------------------------------------------------------------------
""" % (self._path,self._package_name,cmd,self._version)
        elif self._package_name == 'cellranger':
            header = "%s %s-%s" % (self._package_name,self._package_name,
                                   self._version)
        elif self._package_name == 'cellranger-arc':
            header = "%s %s-%s" % (self._package_name,self._package_name,
                                   self._version)
            if version[0] >= 2:
                header += '\n'
        elif self._package_name in ('spaceranger,'):
            if version[0] >= 2 or (version[0] == 1 and version[1] >= 3):
               header = "%s %s-%s\n" % (self._package_name,self._package_name,
                                        self._version)
            else:
                header = "%s %s" % (self._package_name,self._version)
        # Handle version request or no args
        sys.stdout.write(header)
        if cmd == " --version" or not cmd:
            return self._exit_code
        # Build top-level parser
        p = argparse.ArgumentParser()
        sp = p.add_subparsers(dest='command')
        # mkfastq subparser
        mkfastq = sp.add_parser("mkfastq")
        mkfastq.add_argument("--samplesheet",action="store")
        mkfastq.add_argument("--run",action="store")
        mkfastq.add_argument("--output-dir",action="store")
        mkfastq.add_argument("--lanes",action="store")
        mkfastq.add_argument("--use-bases-mask",action="store")
        mkfastq.add_argument("--minimum-trimmed-read-length",action="store")
        mkfastq.add_argument("--mask-short-adapter-reads",action="store")
        mkfastq.add_argument("--ignore-dual-index",action="store_true")
        mkfastq.add_argument("--filter-single-index",action="store_true")
        mkfastq.add_argument("--filter-dual-index",action="store_true")
        mkfastq.add_argument("--jobmode",action="store")
        mkfastq.add_argument("--localcores",action="store")
        mkfastq.add_argument("--localmem",action="store")
        mkfastq.add_argument("--mempercore",action="store")
        mkfastq.add_argument("--maxjobs",action="store")
        mkfastq.add_argument("--jobinterval",action="store")
        mkfastq.add_argument("--disable-ui",action="store_true")
        include_qc_arg = True
        if self._package_name == "cellranger" and version[0] >= 6:
            # --qc removed in cellranger 6.0.0
            include_qc_arg = False
        elif self._package_name == "cellranger-atac" and version[0] >= 2:
            # --qc removed in cellranger-atac 2.0.0
            include_qc_arg = False
        elif self._package_name == "cellranger-arc" and version[0] >= 2:
            # --qc removed in cellranger-arc 2.0.0
            include_qc_arg = False
        elif self._package_name == "spaceranger" and \
             (version[0] >= 2 or (version[0] == 1 and version[1] >= 3)):
            # --qc removed in spaceranger 1.3.0?
             include_qc_arg = False
        if include_qc_arg:
            mkfastq.add_argument("--qc",action="store_true")
        # count subparser
        count = sp.add_parser("count")
        count.add_argument("--id",action="store")
        count.add_argument("--fastqs",action="store")
        count.add_argument("--sample",action="store")
        if self._package_name == "cellranger":
            count.add_argument("--transcriptome",action="store")
            count.add_argument("--chemistry",action="store")
            count.add_argument("--force-cells",action="store",type=int)
            if version[0] == 7:
                # Cellranger 7: include introns on by default
                count.add_argument("--include-introns",
                                   choices=['true','false'],
                                   default='true')
            elif version[0] in (5,6):
                # Cellranger 5,6: include introns off by default
                count.add_argument("--include-introns",
                                   action="store_true")
            if version[0] >= 5:
                count.add_argument("--r1-length",action="store")
                count.add_argument("--r2-length",action="store")
        elif self._package_name == "cellranger-atac":
            count.add_argument("--reference",action="store")
            if version[0] >= 2:
                count.add_argument("--chemistry",choices=['ATAC-v1',
                                                          'ARC-v1'])
            count.add_argument("--force-cells",action="store",type=int)
        elif self._package_name == "cellranger-arc":
            count.add_argument("--reference",action="store")
            count.add_argument("--libraries",action="store")
        count.add_argument("--jobmode",action="store")
        count.add_argument("--localcores",action="store")
        count.add_argument("--localmem",action="store")
        count.add_argument("--mempercore",action="store")
        count.add_argument("--maxjobs",action="store")
        count.add_argument("--jobinterval",action="store")
        # multi subparser
        if self._package_name == "cellranger" and version[0] >= 6:
            # multi only implemented in cellranger 6.0.0
            multi = sp.add_parser("multi")
            multi.add_argument("--id",action="store")
            multi.add_argument("--csv",action="store")
            multi.add_argument("--jobmode",action="store")
            multi.add_argument("--localcores",action="store")
            multi.add_argument("--localmem",action="store")
            multi.add_argument("--mempercore",action="store")
            multi.add_argument("--maxjobs",action="store")
            multi.add_argument("--jobinterval",action="store")
        # Process command line
        args = p.parse_args()
        # Check bases mask
        if self._assert_bases_mask:
            print("Checking bases mask: %s" % args.use_bases_mask)
            assert(args.use_bases_mask == self._assert_bases_mask)
        # Check --include-introns
        if self._assert_include_introns is not None:
            print("Checking --include-introns")
            try:
                include_introns = args.include_introns
            except AttributeError:
                # Earlier versions don't support --include_introns
                include_introns = False
            if include_introns == 'true':
                include_introns = True
            elif include_introns == 'false':
                include_introns = False
            assert(include_introns == self._assert_include_introns)
        # Check --chemistry
        if self._assert_chemistry:
            print("Checking chemistry: %s" % args.chemistry)
            assert(args.chemistry == self._assert_chemistry)
        # Check --force-cells
        if self._assert_force_cells:
            print("Checking --force-cells: %s" % args.force_cells)
            assert(args.force_cells == self._assert_force_cells)
        # Check --filter-single-index
        if self._assert_filter_single_index is not None:
            print("Checking --filter-single-index: %s" %
                  args.filter_single_index)
            assert(args.filter_single_index ==
                   self._assert_filter_single_index)
        # Check --filter-dual-index
        if self._assert_filter_dual_index is not None:
            print("Checking --filter-dual-index: %s" %
                  args.filter_dual_index)
            assert(args.filter_dual_index ==
                   self._assert_filter_dual_index)
        # Handle commands
        if args.command == "mkfastq":
            ##################
            # mkfastq command
            ##################
            # Run folder
            run = args.run
            print("Run folder: %s" % run)
            # Sample sheet
            sample_sheet = args.samplesheet
            print("Sample sheet: %s" % sample_sheet)
            # Output folder
            output_dir = args.output_dir
            print("Output dir: %s" % output_dir)
            # Lanes
            s = SampleSheet(sample_sheet)
            if args.lanes:
                lanes = [int(l) for l in args.lanes.split(',')]
            elif s.has_lanes:
                lanes = [line['Lane'] for line in s.data]
            else:
                lanes = IlluminaRun(args.run).lanes
            print("Lanes: %s" % lanes)
            # Generate mock output based on inputs
            if self._package_name in ('cellranger',
                                      'cellranger-atac',):
                force_sample_dir = True
            elif self._package_name == 'cellranger-arc':
                if self._multiome_data == 'GEX':
                    force_sample_dir = False
                elif self._multiome_data == 'ATAC':
                    force_sample_dir = True
            elif self._package_name == 'spaceranger':
                force_sample_dir = False
            tmpname = "tmp.%s" % uuid.uuid4()
            output = MockIlluminaData(name=tmpname,
                                      package="bcl2fastq2",
                                      unaligned_dir="bcl2fastq")
            s = SampleSheetPredictor(sample_sheet_file=sample_sheet)
            s.set(paired_end=True,
                  lanes=lanes,
                  force_sample_dir=force_sample_dir)
            for project in s.projects:
                print("Adding project: %s" % project.name)
                for sample in project.samples:
                    for fq in sample.fastqs():
                        if sample.sample_name is None:
                            sample_name = sample.sample_id
                        else:
                            sample_name = sample.sample_name
                        output.add_fastq(project.name,
                                         sample_name,
                                         fq)
            # Add undetermined fastqs
            # NB Would like to use the 'add_undetermined'
            # method but this doesn't play well with using
            # the predictor-based approach above
            for r in self._reads:
                for lane in lanes:
                    output.add_fastq(
                        "Undetermined_indices",
                        "undetermined",
                        "Undetermined_S0_L%03d_%s_001.fastq.gz"
                        % (lane,r))
            # Build the output directory
            output.create(force_sample_dir=force_sample_dir)
            # Move to final location
            os.rename(os.path.join(tmpname,"bcl2fastq"),
                      output_dir)
            shutil.rmtree(tmpname)
            # Create cellranger-specific outputs
            if args.lanes:
                lanes_ext = "_%s" % ''.join([str(l) for l in lanes])
            else:
                lanes_ext = ''
            try:
                # Try to get the flow cell from RunInfo.xml
                run_info = IlluminaRunInfo(os.path.join(run,"RunInfo.xml"))
                flow_cell_dir = run_info.flowcell + lanes_ext
            except Exception:
                # Fallback to extracting from the name
                flow_cell_dir = flow_cell_id(run) + lanes_ext
            os.mkdir(flow_cell_dir)
            outs_dir = os.path.join(flow_cell_dir,"outs")
            os.mkdir(outs_dir)
            # Add qc metric files
            if include_qc_arg:
                if args.qc:
                    json_file = os.path.join(outs_dir,"qc_summary.json")
                    with open(json_file,'w') as fp:
                        fp.write(mock10xdata.QC_SUMMARY_JSON)
        elif args.command == "count":
            ###############
            # count command
            ###############
            missing_args = []
            if self._package_name == "cellranger-arc":
                if args.reference is None:
                    missing_args.append("--reference <PATH>")
                if args.libraries is None:
                    missing_args.append("--libraries <CSV>")
            if missing_args:
                sys.stderr.write("error: The following required arguments were not provided:\n")
                for missing_arg in missing_args:
                    sys.stderr.write("    %s\n" % missing_arg)
                sys.stderr.write("\nUSAGE:\n    cellranger-arc count --id <ID> --reference <PATH> --libraries <CSV> --jobmode <MODE>\n\nFor more information try --help\n")
                sys.exit(1)
            # Build outputs
            top_dir = str(args.id)
            os.mkdir(top_dir)
            # Outs
            outs_dir = os.path.join(top_dir,"outs")
            os.mkdir(outs_dir)
            if self._package_name == "cellranger":
                metrics_file = os.path.join(outs_dir,"metrics_summary.csv")
                with open(metrics_file,'w') as fp:
                    fp.write(mock10xdata.METRICS_SUMMARY)
            elif self._package_name in "cellranger-atac":
                summary_file = os.path.join(outs_dir,"summary.csv")
                with open(summary_file,'w') as fp:
                    if version[0] < 2:
                        # Format for cellranger-atac < 2.0.0
                        fp.write(mock10xdata.ATAC_SUMMARY)
                    else:
                        # Format for cellranger-atac >= 2.0.0
                        fp.write(mock10xdata.ATAC_SUMMARY_2_0_0)
            elif self._package_name == "cellranger-arc":
                summary_file = os.path.join(outs_dir,"summary.csv")
                with open(summary_file,'w') as fp:
                    if version[0] < 2:
                        # Format for cellranger-arc < 2.0.0
                        fp.write(mock10xdata.MULTIOME_SUMMARY)
                    else:
                        # Format for cellranger-arc >= 2.0.0
                        fp.write(mock10xdata.MULTIOME_SUMMARY_2_0_0)
            web_summary_file = os.path.join(outs_dir,"web_summary.html")
            with open(web_summary_file,'w') as fp:
                fp.write("PLACEHOLDER FOR WEB_SUMMARY.HTML")
            # _cmdline file
            cmdline_file = os.path.join(top_dir,"_cmdline")
            with open(cmdline_file,'w') as fp:
                fp.write("%s\n" % cmdline)
        elif args.command == "multi":
            ###############
            # multi command
            ###############
            missing_args = []
            if args.csv is None:
                missing_args.append("--csv <PATH>")
            if missing_args:
                sys.stderr.write("error: The following required arguments were not provided:\n")
                for missing_arg in missing_args:
                    sys.stderr.write("    %s\n" % missing_arg)
                sys.stderr.write("\nUSAGE:\n    cellranger multi --id <ID> --csv <CSV> --jobmode <MODE>\n\nFor more information try --help\n")
                sys.exit(1)
            # Build outputs
            top_dir = str(args.id)
            os.mkdir(top_dir)
            # Outs
            outs_dir = os.path.join(top_dir,"outs")
            os.mkdir(outs_dir)
            config = CellrangerMultiConfigCsv(args.csv)
            for d in ("per_sample_outs",
                      "multi",
                      "multi/multiplexing_analysis",):
                os.mkdir(os.path.join(outs_dir,d))
            # Per sample outs
            for sample in config.sample_names:
                sample_dir = os.path.join(outs_dir,
                                          "per_sample_outs",
                                          sample)
                os.mkdir(sample_dir)
                metrics_file = os.path.join(sample_dir,
                                            "metrics_summary.csv")
                with open(metrics_file,'wt') as fp:
                    fp.write(mock10xdata.CELLPLEX_METRICS_SUMMARY)
                web_summary_file = os.path.join(sample_dir,
                                                "web_summary.html")
                with open(web_summary_file,'wt') as fp:
                    fp.write("PLACEHOLDER FOR WEB_SUMMARY.HTML")
            # Multiplexing outs
            for f in ("assignment_confidence_table.csv",
                      "cells_per_tag.json",
                      "tag_calls_per_cell.csv",
                      "tag_calls_summary.csv"):
                multiplexing_output = os.path.join(outs_dir,
                                                   "multi",
                                                   "multiplexing_analysis",f)
                with open(multiplexing_output,'wt') as fp:
                    fp.write("PLACEHOLDER FOR %s"  %
                             os.path.basename(multiplexing_output).upper())
            # _cmdline file
            cmdline_file = os.path.join(top_dir,"_cmdline")
            with open(cmdline_file,'w') as fp:
                fp.write("%s\n" % cmdline)
        else:
            print("%s: not implemented" % args.command)
        print("Return exit code: %s" % self._exit_code)
        return self._exit_code

class MockFastqScreen:
    """
    Create mock fastq_screen

    This class can be used to create a mock
    fastq_screen executable, which in turn can be used
    in place of the actual fastq_screen program for
    testing purposes.

    To create a mock executable, use the 'create' static
    method, e.g.

    >>> MockFastqScreen.create("/tmpbin/fastq_screen")

    The resulting executable will generate mock outputs
    when run on a Fastq file (ignoring its content).

    The executable can be configured on creation to
    produce different error conditions when run:

    - the exit code can be set to an arbitrary value
      via the `exit_code` argument
    - outputs for specific stages can be removed by
      specifying their names in the `missing_fastqs`
      argument
    """

    @staticmethod
    def create(path,version=None,no_outputs=False,
               exit_code=0):
        """
        Create a "mock" fastq_screen executable

        Arguments:
          path (str): path to the new executable
            to create. The final executable must
            not exist, however the directory it
            will be created in must.
          version (str): explicit version string
          no_outputs (bool): if True then don't
            create outputs (default: False, do
            create outputs)
          exit_code (int): exit code that the
            mock executable should complete
            with
        """
        path = os.path.abspath(path)
        print("Building mock executable: %s" % path)
        # Don't clobber an existing executable
        assert(os.path.exists(path) is False)
        with open(path,'w') as fp:
            fp.write("""#!/usr/bin/env python
import sys
from auto_process_ngs.mock import MockFastqScreen
sys.exit(MockFastqScreen(version=%s,
                         no_outputs=%s,
                         exit_code=%s).main(sys.argv[1:]))
            """ % (("\"%s\"" % version
                    if version is not None
                    else None),
                   no_outputs,
                   exit_code))
            os.chmod(path,0o775)
        with open(path,'r') as fp:
            print("fastq_screen:")
            print("%s" % fp.read())
        return path

    def __init__(self,version=None,no_outputs=False,exit_code=0):
        """
        Internal: configure the fastq_screen
        """
        if version is None:
            version = "0.14.0"
        self._version = str(version)
        self._no_outputs = no_outputs
        self._exit_code = exit_code

    def main(self,args):
        """
        Internal: provides mock fastq_screen functionality
        """
        # No args
        if not args:
            return self._exit_code
        # Deal with arguments
        p = argparse.ArgumentParser()
        p.add_argument("--conf",action="store")
        p.add_argument("--outdir",action="store")
        p.add_argument("--threads",action="store")
        p.add_argument("--subset",action="store")
        p.add_argument("--force",action="store_true")
        p.add_argument("fastq")
        args = p.parse_args(args)
        # Check input file
        if not os.path.exists(args.fastq):
            print("%s: fastq file not found" % args.fastq)
            return 1
        # Output dir
        if args.outdir:
            outdir = args.outdir
        else:
            outdir = os.getcwd()
        outdir = os.path.abspath(outdir)
        # Fastq base name
        fastq_base = MockQCOutputs.fastq_basename(args.fastq)
        # Create screen outputs
        if not self._no_outputs:
            MockQCOutputs.fastq_screen_v0_9_2(args.fastq,outdir)
        return self._exit_code

class MockFastQC:
    """
    Create mock fastqc

    This class can be used to create a mock
    fastqc executable, which in turn can be used
    in place of the actual fastqc program for
    testing purposes.

    To create a mock script, use the 'create' static
    method, e.g.

    >>> MockFastQC.create("/tmpbin/fastqc")

    The resulting executable will generate mock outputs
    when run on Fastq files (ignoring their content).

    The executable can be configured on creation to
    produce different error conditions when run:

    - the exit code can be set to an arbitrary value
      via the `exit_code` argument
    """

    @staticmethod
    def create(path,version=None,no_outputs=False,
               exit_code=0):
        """
        Create a "mock" illumina.sh "script"

        Arguments:
          path (str): path to the new executable
            to create. The final executable must
            not exist, however the directory it
            will be created in must.
          version (str): explicit version string
          no_outputs (bool): if True then make
            don't create mock outputs for FastQC
          exit_code (int): exit code that the
            mock executable should complete
            with
        """
        path = os.path.abspath(path)
        print("Building mock executable: %s" % path)
        # Don't clobber an existing executable
        assert(os.path.exists(path) is False)
        with open(path,'w') as fp:
            fp.write("""#!/usr/bin/env python
import sys
from auto_process_ngs.mock import MockFastQC
sys.exit(MockFastQC(version=%s,
                    no_outputs=%s,
                    exit_code=%s).main(sys.argv[1:]))
            """ % (("\"%s\"" % version
                    if version is not None
                    else None),
                   no_outputs,
                   exit_code))
            os.chmod(path,0o775)
        with open(path,'r') as fp:
            print("fastqc:")
            print("%s" % fp.read())
        return path

    def __init__(self,version=None,no_outputs=False,
                 exit_code=0):
        """
        Internal: configure the mock fastqc
        """
        if version is None:
            version = "0.11.3"
        self._version = str(version)
        self._no_outputs = no_outputs
        self._exit_code = exit_code

    def main(self,args):
        """
        Internal: provides mock fastqc functionality
        """
        # No args
        if not args:
            return self._exit_code
        # Deal with arguments
        p = argparse.ArgumentParser()
        p.add_argument("--threads",action="store")
        p.add_argument("--outdir",action="store")
        p.add_argument("--nogroup",action="store_true")
        p.add_argument("--extract",action="store_true")
        p.add_argument("fastq",nargs='+')
        args = p.parse_args(args)
        # Output dir
        if args.outdir:
            outdir = args.outdir
        else:
            outdir = os.getcwd()
        outdir = os.path.abspath(outdir)
        # Create FastQC outputs
        if not self._no_outputs:
            for fastq in args.fastq:
                # Check input file
                if not os.path.exists(fastq):
                    print("%s: fastq file not found" % fastq)
                    return 1
                # Fastq base name
                MockQCOutputs.fastqc_v0_11_2(fastq,outdir)
        return self._exit_code

class MockMultiQC:
    """
    Create mock MultiQC executable

    This class can be used to create a mock multiqc
    executable, which in turn can be used in place of
    the actual multiqc executable for testing purposes.

    To create a mock executable, use the 'create' static
    method, e.g.

    >>> MockMultiQC.create("/tmpbin/multiqc")

    The resulting executable will generate mock outputs
    when run on a directory (ignoring its contents).

    The executable can be configured on creation to
    produce different error conditions when run:

    - the exit code can be set to an arbitrary value
      via the `exit_code` argument
    - the outputs can be suppressed by setting the
      `no_output` argument to `True`
    """

    @staticmethod
    def create(path,version=None,no_outputs=False,
               exit_code=0):
        """
        Create a "mock" multiqc executable

        Arguments:
          path (str): path to the new executable
            to create. The final executable must
            not exist, however the directory it
            will be created in must.
          version (str): explicit version string
          no_outputs (bool): if True then don't
            create any of the expected outputs
          exit_code (int): exit code that the
            mock executable should complete
            with
        """
        path = os.path.abspath(path)
        print("Building mock executable: %s" % path)
        # Don't clobber an existing executable
        assert(os.path.exists(path) is False)
        with open(path,'w') as fp:
            fp.write("""#!/usr/bin/env python
import sys
from auto_process_ngs.mock import MockMultiQC
sys.exit(MockMultiQC(version=%s,
                     no_outputs=%s,
                     exit_code=%s).main(sys.argv[1:]))
            """ % (("\"%s\"" % version if version else None),
                   no_outputs,
                   exit_code))
            os.chmod(path,0o775)
        with open(path,'r') as fp:
            print("multiqc:")
            print("%s" % fp.read())
        return path

    def __init__(self,version=None,no_outputs=False,exit_code=0):
        """
        Internal: configure the mock multiqc
        """
        if version:
            self._version = version
        else:
            self._version = "1.5"
        self._no_outputs = no_outputs
        self._exit_code = exit_code

    def main(self,args):
        """
        Internal: provides mock multiqc functionality
        """
        # No args
        if not args:
            return self._exit_code
        # Handle version request
        if args[0] == "--version":
            print("multiqc, version %s" % self._version)
            return self._exit_code
        # Deal with arguments
        p = argparse.ArgumentParser()
        p.add_argument("--title",action="store")
        p.add_argument("--filename",action="store")
        p.add_argument("--force",action="store_true")
        p.add_argument("analysis_directory",nargs='+')
        args = p.parse_args(args)
        # Check input directory
        for d in args.analysis_directory:
            if not os.path.exists(d):
                sys.stderr.write("""Usage: multiqc [OPTIONS] <analysis directory>

Error: Invalid value for "analysis_dir": Path "%s" does not exist.

This is MultiQC v%s

For more help, run 'multiqc --help' or visit http://multiqc.info
""" % (d,self._version))
                return 2
        # Outputs
        if args.filename is None:
            out_file = "multiqc_report.html"
            out_dir = "multiqc_data"
        else:
            out_file = args.filename
            out_dir = "%s_data" % os.path.splitext(out_file)[0]
        if not self._no_outputs:
            MockQCOutputs.multiqc(d,multiqc_html=out_file,
                                  version=self._version)
            os.mkdir(out_dir)
        # Exit
        return self._exit_code

class MockFastqStrandPy:
    """
    Create mock fastq_strand.py executable

    This class can be used to create a mock fastq_strand.py
    executable, which in turn can be used in place of
    the actual fastq_strand.py executable for testing
    purposes.

    To create a mock executable, use the 'create' static
    method, e.g.

    >>> MockFastqStrandPy.create("/tmpbin/fastq_strand.py")

    The resulting executable will generate mock outputs
    when run on a pair of Fastq files (ignoring their
    contents).

    The executable can be configured on creation to
    produce different error conditions when run:

    - the exit code can be set to an arbitrary value
      via the `exit_code` argument
    - the outputs can be suppressed by setting the
      `no_output` argument to `True`
    """

    @staticmethod
    def create(path,no_outputs=False,exit_code=0):
        """
        Create a "mock" fastq_strand.py executable

        Arguments:
          path (str): path to the new executable
            to create. The final executable must
            not exist, however the directory it
            will be created in must.
          no_outputs (bool): if True then don't
            create any of the expected outputs
          exit_code (int): exit code that the
            mock executable should complete
            with
        """
        path = os.path.abspath(path)
        print("Building mock executable: %s" % path)
        # Don't clobber an existing executable
        assert(os.path.exists(path) is False)
        with open(path,'w') as fp:
            fp.write("""#!/usr/bin/env python
import sys
from auto_process_ngs.mock import MockFastqStrandPy
sys.exit(MockFastqStrandPy(no_outputs=%s,
                     exit_code=%s).main(sys.argv[1:]))
            """ % (no_outputs,exit_code))
            os.chmod(path,0o775)
        with open(path,'r') as fp:
            print("fastq_strand.py:")
            print("%s" % fp.read())
        return path

    def __init__(self,no_outputs=False,exit_code=0):
        """
        Internal: configure the mock fastq_strand.py
        """
        self._version = "0.0.2"
        self._no_outputs = no_outputs
        self._exit_code = exit_code

    def main(self,args):
        """
        Internal: provides mock fastq_strand.py functionality
        """
        # No args
        if not args:
            return self._exit_code
        # Handle version request
        if args[0] == "--version":
            print("%s" % self._version)
            return self._exit_code
        # Deal with arguments
        p = argparse.ArgumentParser()
        p.add_argument("--subset",action="store")
        p.add_argument("-o","--outdir",action="store")
        p.add_argument("-c","--conf",action="store")
        p.add_argument("-n",action="store")
        p.add_argument("--counts",action="store_true")
        p.add_argument("--keep-star-output",action="store_true")
        p.add_argument("fastqr1")
        p.add_argument("fastqr2",nargs="?")
        args = p.parse_args(args)
        # Outputs
        if self._no_outputs:
            return self._exit_code
        # Create fake output file
        basename = strip_ngs_extensions(os.path.basename(args.fastqr1))
        if args.outdir is not None:
            outfile = os.path.join(args.outdir,"%s_fastq_strand.txt" % basename)
        else:
            outfile = os.path.join(os.getcwd(),"%s_fastq_strand.txt" % basename)
        with open(outfile,'w') as fp:
            fp.write("""#fastq_strand version: %s	#Aligner: STAR	#Reads in subset: 1000
#Genome	1st forward	2nd reverse
""" % self._version)
            with open(args.conf,'r') as conf:
                for line in conf:
                    if not line or line.startswith("#"):
                        continue
                    genome = line.split('\t')[0]
                    fp.write("%s	13.13	93.21\n" % genome)
        # Exit
        return self._exit_code

class MockConda:
    """
    Create mock conda installation

    This class can be used to create a mock conda
    installation consisting of:

    - ``bin`` subdirectory with mock ``conda`` executable
      and 'activate' script
    - ``envs`` subdirectory

    This can be used in place of an actual conda
    installation for testing purposes.

    To create a mock installation, use the 'create'
    static method, e.g.

    >>> MockCondaExe.create("/tmpbin/conda")

    The resulting ``conda`` executable supports
    ``--version`` and the ``create`` command, and
    will generate mock outputs for both.

    The executable can be configured on creation to
    produce different error conditions when run:

    - the exit code can be set to an arbitrary value
      via the `exit_code` argument
    - the 'create' command can be forced to fail for
      all inputs by setting the `create_fails` argument
    - the reported version can be set via the `version`
      argument
    """

    @staticmethod
    def create(path,version="4.10.3",create_fails=False,
               activate_fails=False,exit_code=0):
        """
        Create a "mock" fastq_strand.py executable

        Arguments:
          path (str): path to the top-level directory
            for the mock conda installation (which
            must not exist, however the directory it
            will be created in must be present).
          version (str): version that mock conda
            will claim to be
          create_fails (bool): if True then the
            'create' subcommand of the mock
            conda executable will fail.
          activate_fails (bool): if True then
            the 'activate' script in the mock conda
            installation will return with value 1
            (i.e. an error code).
          exit_code (int): exit code that the
            mock executable should complete with
        """
        path = os.path.abspath(path)
        print("Building mock installation: %s" % path)
        # Don't clobber an existing installation
        assert(os.path.exists(path) is False)
        # Set up directories
        os.mkdir(path)
        for d in ("bin","envs"):
            os.mkdir(os.path.join(path,d))
        # Create mock files
        bin_dir = os.path.join(path,"bin")
        conda_ = os.path.join(bin_dir,"conda")
        if not create_fails:
            with open(conda_,'wt') as fp:
                fp.write("""#!/bin/bash
if [ "$1" == "--version" ] ; then
   echo "conda %s"
   exit 0
elif [ "$1" == "create" ] ; then
   YES=
   PREFIX=
   PACKAGES=
   while [ ! -z "$2" ] ; do
     case "$2" in
       -n)
         shift
         PREFIX=$(dirname $(dirname $0))/envs/${2}
         ;;
       --prefix)
         shift
         PREFIX=$2
         ;;
       -y)
         YES=yes
         ;;
       -c)
         shift
         ;;
       --override-channels)
         ;;
       *)
         PACKAGES="$PACKAGES $2"
         ;;
     esac
     shift
   done
fi
if [ -z "$YES" ] ; then
   echo "Need to supply -y option"
   exit 1
fi
if [ -z "$PREFIX" ] ; then
   echo "Need to supply either -n or --prefix"
   exit 1
fi
# Make directory for new environment
mkdir -p $PREFIX
# Write package list to a 'packages.txt' file
echo $PACKAGES >${PREFIX}/packages.txt
# Make an executable script for each package name
for pkg in $PACKAGES ; do
   name=$(echo $pkg | cut -f1 -d=)
   cat >${PREFIX}/${name} <<EOF
#!/bin/bash
echo \\$1
exit %s
EOF
   chmod +x ${PREFIX}/${name}
done
""" % (version,exit_code))
        else:
            with open(conda_,'wt') as fp:
                fp.write("""#!/bin/bash
if [ "$1" == "--version" ] ; then
   echo "conda %s"
   exit 0
elif [ "$1" == "create" ] ; then
   echo "!!!! Failed to create environment !!!!"
   exit 1
fi
""" % version)
        os.chmod(conda_,0o775)
        # Make mock 'activate' script
        activate_ = os.path.join(bin_dir,"activate")
        with open(activate_,'wt') as fp:
            fp.write("""#!/bin/bash
export PATH=$PATH:${1}
""")
            if activate_fails:
                fp.write("""echo Activate failed
return 1
""")
            os.chmod(activate_,0o755)
        with open(conda_,'rt') as fp:
            print("conda:")
            print("%s" % fp.read())
        return path

class MockCellrangerExe(Mock10xPackageExe):
    """
    Wrapper for Mock10xPackageExe

    Maintained for backwards-compatibility
    """

#######################################################################
# Functions for creating mock data
#######################################################################

def make_mock_bcl2fastq2_output(out_dir,lanes,sample_sheet=None,
                                reads=None,no_lane_splitting=False,
                                exclude_fastqs=None,
                                create_fastq_for_index_read=False,
                                paired_end=False,
                                force_sample_dir=False):
    """
    Creates files & directories structure mimicking output from bcl2fastq2

    Arguments:
      out_dir (str): path to output directory
      lanes (iterable): list of lanes to create output
        for
      sample_sheet (str): path to sample sheet file
      reads (iterable): list of 'reads' to create (e.g.
        ('R1','R2'); defaults to ('R1') if not specified
      no_lane_splitting (bool): whether to produce mock
        Fastq files for each lane, or combine them
        across lanes (mimics the --no-lane-splitting
        option in `bcl2fastq`)
      exclude_fastqs (iterable): specifies a list of
        Fastq files to exclude from the outputs
      create_fastq_for_index_read (bool): whether to
        also include 'I1' etc Fastqs for index reads
        (ignored if 'reads' argument is set)
      paired_end (bool): whether to also include 'R2'
        and 'I2' Fastqs (ignored if 'reads' argument is
        set)
      force_sample_dir (bool): whether to force insertion
        of a 'sample name' directory for IEM4 sample
        sheets where sample name and ID are the same
    """
    # Setup up base mock output directory
    top_dir = os.path.dirname(os.path.abspath(out_dir))
    unaligned_dir = os.path.basename(out_dir)
    output = MockIlluminaData(name=top_dir,
                              package="bcl2fastq2",
                              unaligned_dir=unaligned_dir)
    # Sort out reads if not explicitly set
    if reads is None:
        reads = ['R1']
        if paired_end:
            reads.append('R2')
        if create_fastq_for_index_read:
            reads.append('I1')
            if paired_end:
                reads.append('I2')
    # Add outputs from sample sheet (if supplied)
    if sample_sheet is not None:
        # Cut down samplesheet if lanes were specified
        sample_sheet_ = SampleSheet(sample_sheet)
        i = 0
        if lanes and sample_sheet_.has_lanes:
            while i < len(sample_sheet_):
                line = sample_sheet_[i]
                if line['Lane'] in lanes:
                    i += 1
                else:
                    del(sample_sheet_[i])
        # Get lanes from modified sample sheet
        if sample_sheet_.has_lanes:
            lanes = sorted(list(set([int(line['Lane'])
                                     for line in sample_sheet_])))
            print("Updated lanes: %s" % lanes)
        # Predict outputs
        s = SampleSheetPredictor(sample_sheet=sample_sheet_)
        s.set(reads=reads,
              lanes=(lanes if not sample_sheet_.has_lanes else None),
              no_lane_splitting=no_lane_splitting,
              force_sample_dir=force_sample_dir)
        # Build and populate outputs
        for project in s.projects:
                print("Adding project: %s" % project.name)
                for sample in project.samples:
                    for fq in sample.fastqs():
                        print("- %s" % fq)
                        if exclude_fastqs and (fq in exclude_fastqs):
                            continue
                        if sample.sample_name is None:
                            sample_name = sample.sample_id
                        else:
                            sample_name = sample.sample_name
                        output.add_fastq(project.name,
                                         sample_name,
                                         fq)
    # Add undetermined fastqs
    # NB Would like to use the 'add_undetermined'
    # method but this doesn't play well with using
    # the predictor-based approach above
    print("Adding 'undetermined' Fastqs")
    for r in reads:
        if no_lane_splitting:
            output.add_fastq(
                "Undetermined_indices",
                "undetermined",
                "Undetermined_S0_%s_001.fastq.gz" % r)
        else:
            for lane in lanes:
                output.add_fastq(
                    "Undetermined_indices",
                    "undetermined",
                    "Undetermined_S0_L%03d_%s_001.fastq.gz"
                    % (lane,r))
    # Build the output directory
    output.create()
    # Populate the Fastqs with reads
    illuminadata = IlluminaData(top_dir,unaligned_dir=unaligned_dir)
    fastqs = []
    for p in illuminadata.projects:
        for s in p.samples:
            for fq in s.fastq:
                fastqs.append(os.path.join(s.dirn,fq))
    for s in illuminadata.undetermined.samples:
        for fq in s.fastq:
            fastqs.append(os.path.join(s.dirn,fq))
    for fastq in fastqs:
        read_number = IlluminaFastq(fastq).read_number
        with gzip.open(fastq,'wt') as fp:
            if no_lane_splitting:
                # Add one read per lane
                for lane in lanes:
                    read = """@ILLUMINA-545855:49:FC61RLR:%s:1:10979:1695 %s:N:0:TCCTGA
GCATACTCAGCTTTAGTAATAAGTGTGATTCTGGTA
+
IIIIIHIIIGHHIIDGHIIIIIIHIIIIIIIIIIIH\n""" % (lane,read_number)
                    fp.write(read)
            else:
                lane = IlluminaFastq(fastq).lane_number
                read = """@ILLUMINA-545855:49:FC61RLR:%s:1:10979:1695 %s:N:0:TCCTGA
GCATACTCAGCTTTAGTAATAAGTGTGATTCTGGTA
+
IIIIIHIIIGHHIIDGHIIIIIIHIIIIIIIIIIIH\n""" % (lane,read_number)
                fp.write(read)

def make_mock_analysis_project(name="PJB",top_dir=None,
                               protocol=None,paired_end=True,
                               qc_dir="qc",
                               fastq_dir='fastqs',
                               fastq_names=None,
                               sample_names=None,
                               screens=('model_organisms',
                                        'other_organisms',
                                        'rRNA',),
                               include_fastqc=True,
                               include_fastq_screen=True,
                               include_strandedness=True,
                               include_seqlens=True,
                               include_multiqc=True,
                               include_cellranger_count=False,
                               include_cellranger_multi=False,
                               cellranger_pipelines=('cellranger',),
                               cellranger_samples=None,
                               cellranger_multi_samples=None,
                               legacy_screens=False,
                               legacy_cellranger_outs=False):
    """
    Create a mock Analysis Project directory with QC artefacts

    Arguments:
      name (str): name for the mock project
      top_dir (str): path to the directory to create the
        mock project directory under
      protocol (str): QC protocol to emulate
      paired_end (bool): whether the mock project should
        be paired-end (the default)
      fastq_dir (str): optional, set a non-standard
        directory for the Fastq files
      fastq_names (list): optional, explicit list of Fastq
        names
      sample_names (list): optional, explicit list of sample
        names
      screens (list): optional, list of non-standard
        FastqScreen panel names
      include_fastqc (bool): include outputs from Fastqc
      include_fastq_screen (bool): include outputs from
        FastqScreen
      include_strandedness (bool): include outputs from
        strandedness
      include_seqlens (bool): include sequence length metrics
      include_multiqc (bool): include MultiQC outputs
      include_celllranger_count (bool): include 'cellranger
        count' outputs
      include_cellranger_multi (bool): include 'cellranger
        multi' outputs
      cellranger_pipelines (list): list of 10xGenomics pipelines
        to make mock outputs for (e.g. 'cellranger',
        'cellranger-atac' etc)
      cellranger_samples (list): list of sample names to
        produce 'cellranger count' outputs for
      cellranger_multi_samples (list): list of sample names to
        produce 'cellranger multi' outputs for
      legacy_screens (bool): if True then use legacy naming
        convention for FastqScreen outputs
      legacy_cellranger_outs (bool): if True then use legacy
        naming convention for 10xGenomics pipeline outputs

    Returns:
      String: path to the mock analysis project that was created.
    """
    # Read numbers
    if paired_end:
        reads = (1,2)
    else:
        reads = (1,)
    if not sample_names:
        sample_names = ('PJB1','PJB2')
    # Generate names for fastq files to add
    if fastq_names is None:
        fastq_names = []
        for i,sname in enumerate(sample_names,start=1):
            for read in reads:
                fq = "%s_S%d_R%d_001.fastq.gz" % (sname,i,read)
                fastq_names.append(fq)
    # Set up the analysis project
    analysis_project = MockAnalysisProject(name,
                                           fastq_names,
                                           fastq_dir=fastq_dir)
    project_dir = analysis_project.create(top_dir=top_dir)
    # Populate with fake QC products
    make_mock_qc_dir(os.path.join(project_dir,qc_dir),
                     fastq_names=fastq_names,
                     fastq_dir=os.path.join(project_dir,fastq_dir),
                     protocol=protocol,
                     screens=screens,
                     cellranger_pipelines=cellranger_pipelines,
                     cellranger_samples=cellranger_samples,
                     cellranger_multi_samples=cellranger_multi_samples,
                     include_fastqc=include_fastqc,
                     include_fastq_screen=include_fastq_screen,
                     include_strandedness=include_strandedness,
                     include_seqlens=include_seqlens,
                     include_multiqc=include_multiqc,
                     include_cellranger_count=include_cellranger_count,
                     include_cellranger_multi=include_cellranger_multi,
                     legacy_screens=legacy_screens,
                     legacy_cellranger_outs=legacy_cellranger_outs)
    return project_dir
