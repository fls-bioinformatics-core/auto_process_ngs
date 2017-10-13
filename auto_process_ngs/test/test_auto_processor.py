#######################################################################
# Tests for autoprocessor.py module
#######################################################################

import unittest
import tempfile
import shutil
import os
import gzip
from bcftbx.mock import MockIlluminaRun
from auto_process_ngs.auto_processor import AutoProcess
from auto_process_ngs.mock import MockAnalysisDirFactory
from auto_process_ngs.mock import MockAnalysisDir
from auto_process_ngs.mock import UpdateAnalysisDir
from auto_process_ngs.mock import UpdateAnalysisProject

# Unit tests

class TestAutoProcessSetup(unittest.TestCase):

    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestAutoProcessSetup')
        # Store original location so we can get back at the end
        self.pwd = os.getcwd()
        # Move to working dir
        os.chdir(self.dirn)
        # Placeholders for test objects
        self.mock_illumina_run = None
        self.analysis_dir = None

    def tearDown(self):
        # Return to original dir
        os.chdir(self.pwd)
        # Remove the temporary test directory
        shutil.rmtree(self.dirn)

    def test_autoprocess_setup(self):
        """AutoProcess.setup works for mock MISeq run
        """
        # Create mock Illumina run directory
        mock_illumina_run = MockIlluminaRun(
            '151125_M00879_0001_000000000-ABCDE1',
            'miseq',
            top_dir=self.dirn)
        mock_illumina_run.create()
        # Set up autoprocessor
        ap = AutoProcess()
        ap.setup(mock_illumina_run.dirn)
        analysis_dirn = "%s_analysis" % mock_illumina_run.name
        # Check parameters
        self.assertEqual(ap.analysis_dir,
                         os.path.join(self.dirn,analysis_dirn))
        self.assertEqual(ap.params.data_dir,mock_illumina_run.dirn)
        self.assertEqual(ap.params.sample_sheet,
                         os.path.join(self.dirn,analysis_dirn,
                                      'custom_SampleSheet.csv'))
        self.assertEqual(ap.params.bases_mask,
                         'y101,I8,I8,y101')
        # Delete to force write of data to disk
        del(ap)
        # Check directory exists
        self.assertTrue(os.path.isdir(analysis_dirn))
        # Check files exists
        for filen in ('SampleSheet.orig.csv',
                      'custom_SampleSheet.csv',
                      'auto_process.info',
                      'metadata.info',):
            self.assertTrue(os.path.exists(os.path.join(analysis_dirn,
                                                        filen)),
                            "Missing file: %s" % filen)
        # Check subdirs have been created
        for subdirn in ('ScriptCode',
                        'logs',):
            self.assertTrue(os.path.isdir(os.path.join(analysis_dirn,
                                                       subdirn)),
                            "Missing subdir: %s" % subdirn)

    def test_autoprocess_setup_existing_target_dir(self):
        """AutoProcess.setup works when target dir exists
        """
        # Create mock Illumina run directory
        mock_illumina_run = MockIlluminaRun(
            '160621_M00879_0087_000000000-AGEW9',
            'miseq',
            top_dir=self.dirn)
        mock_illumina_run.create()
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_M00879_0087_000000000-AGEW9',
            'miseq',
            top_dir=self.dirn)
        mockdir.create()
        # Do setup into existing analysis dir
        ap = AutoProcess()
        ap.setup(mock_illumina_run.dirn)
        self.assertTrue(os.path.isdir(
            '160621_M00879_0087_000000000-AGEW9'))

class TestAutoProcessPublishQc(unittest.TestCase):
    """
    Tests for AutoProcess.publish_qc
    """
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestAutoProcessPublishQc')
        # Store original location so we can get back at the end
        self.pwd = os.getcwd()
        # Move to working dir
        os.chdir(self.dirn)
        # Placeholders for test objects
        self.ap = None

    def tearDown(self):
        # Delete autoprocessor object
        if self.ap is not None:
            del(self.ap)
        # Return to original dir
        os.chdir(self.pwd)
        # Remove the temporary test directory
        shutil.rmtree(self.dirn)

    def test_publish_qc_metadata_missing(self):
        """publish_qc: raises exception if metadata not set
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_M00879_0087_000000000-AGEW9',
            'miseq',
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        ap = AutoProcess(mockdir.dirn)
        # Make a mock publication area
        publication_dir = os.path.join(self.dirn,'QC')
        os.mkdir(publication_dir)
        # Publish
        self.assertRaises(Exception,
                          ap.publish_qc,
                          location=publication_dir)

    def test_publish_qc_processing_qc(self):
        """publish_qc: processing QC report only
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        ap = AutoProcess(mockdir.dirn)
        # Add processing report
        UpdateAnalysisDir(ap).add_processing_report()
        # Make a mock publication area
        publication_dir = os.path.join(self.dirn,'QC')
        os.mkdir(publication_dir)
        # Publish QC
        ap.publish_qc(location=publication_dir)
        # Check outputs
        outputs = ("index.html",
                   "processing_qc.html")
        for item in outputs:
            f = os.path.join(publication_dir,
                             "160621_K00879_0087_000000000-AGEW9_analysis",
                             item)
            self.assertTrue(os.path.exists(f),"Missing %s" % f)

    def test_publish_qc_barcode_analysis(self):
        """publish_qc: barcode analysis outputs
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        ap = AutoProcess(mockdir.dirn)
        # Add processing and barcode analysis reports
        UpdateAnalysisDir(ap).add_processing_report()
        UpdateAnalysisDir(ap).add_barcode_analysis()
        # Make a mock publication area
        publication_dir = os.path.join(self.dirn,'QC')
        os.mkdir(publication_dir)
        # Publish QC
        ap.publish_qc(location=publication_dir)
        # Check outputs
        outputs = ("index.html",
                   "processing_qc.html",
                   os.path.join("barcodes","barcodes.report"),
                   os.path.join("barcodes","barcodes.xls"),
                   os.path.join("barcodes","barcodes.html"))
        for item in outputs:
            f = os.path.join(publication_dir,
                             "160621_K00879_0087_000000000-AGEW9_analysis",
                             item)
            self.assertTrue(os.path.exists(f),"Missing %s" % f)

    def test_publish_qc_with_projects_no_multiqc(self):
        """publish_qc: projects with QC outputs (no MultiQC)
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create()
        ap = AutoProcess(mockdir.dirn)
        # Add processing report and QC outputs
        UpdateAnalysisDir(ap).add_processing_report()
        for project in ap.get_analysis_projects():
            UpdateAnalysisProject(project).add_qc_outputs(
                include_multiqc=False)
        # Make a mock publication area
        publication_dir = os.path.join(self.dirn,'QC')
        os.mkdir(publication_dir)
        # Publish
        ap.publish_qc(location=publication_dir)
        self.assertTrue(
            os.path.exists(
                os.path.join(publication_dir,
                             "160621_K00879_0087_000000000-AGEW9_analysis",
                             "processing_qc.html")))
        # Check outputs
        outputs = ["index.html",
                   "processing_qc.html"]
        for project in ap.get_analysis_projects():
            project_qc = "qc_report.%s.%s" % (project.name,
                                              os.path.basename(
                                                  ap.analysis_dir))
            outputs.append(project_qc)
            outputs.append("%s.zip" % project_qc)
            outputs.append(os.path.join(project_qc,"qc_report.html"))
            outputs.append(os.path.join(project_qc,"qc"))
        for item in outputs:
            f = os.path.join(publication_dir,
                             "160621_K00879_0087_000000000-AGEW9_analysis",
                             item)
            self.assertTrue(os.path.exists(f),"Missing %s" % f)

    def test_publish_qc_with_projects(self):
        """publish_qc: projects with all QC outputs
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create()
        ap = AutoProcess(mockdir.dirn)
        # Add processing report and QC outputs
        UpdateAnalysisDir(ap).add_processing_report()
        for project in ap.get_analysis_projects():
            UpdateAnalysisProject(project).add_qc_outputs()
        # Make a mock publication area
        publication_dir = os.path.join(self.dirn,'QC')
        os.mkdir(publication_dir)
        # Publish
        ap.publish_qc(location=publication_dir)
        # Check outputs
        outputs = ["index.html",
                   "processing_qc.html"]
        for project in ap.get_analysis_projects():
            # Standard QC outputs
            project_qc = "qc_report.%s.%s" % (project.name,
                                              os.path.basename(
                                                  ap.analysis_dir))
            outputs.append(project_qc)
            outputs.append("%s.zip" % project_qc)
            outputs.append(os.path.join(project_qc,"qc_report.html"))
            outputs.append(os.path.join(project_qc,"qc"))
            # MultiQC output
            outputs.append("multiqc_report.%s.html" % project.name)
        for item in outputs:
            f = os.path.join(publication_dir,
                             "160621_K00879_0087_000000000-AGEW9_analysis",
                             item)
            self.assertTrue(os.path.exists(f),"Missing %s" % f)

    def test_publish_qc_with_projects_with_multiple_fastq_sets(self):
        """publish_qc: projects with multiple Fastq sets
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create()
        ap = AutoProcess(mockdir.dirn)
        # Add processing report and QC outputs
        UpdateAnalysisDir(ap).add_processing_report()
        for project in ap.get_analysis_projects():
            UpdateAnalysisProject(project).add_qc_outputs()
        # Add additional fastq set for first project
        multi_fastqs_project = ap.get_analysis_projects()[0]
        UpdateAnalysisProject(multi_fastqs_project).add_fastq_set(
            "fastqs.extra",
            ("Alt1.r1.fastq.gz","Alt2.r1.fastq.gz"))
        UpdateAnalysisProject(multi_fastqs_project).add_qc_outputs(
            fastq_set="fastqs.extra",
            qc_dir="qc.extra")
        # Make a mock publication area
        publication_dir = os.path.join(self.dirn,'QC')
        os.mkdir(publication_dir)
        # Publish
        ap.publish_qc(location=publication_dir)
        # Check outputs
        outputs = ["index.html",
                   "processing_qc.html"]
        for project in ap.get_analysis_projects():
            # Standard QC outputs
            project_qc = "qc_report.%s.%s" % (project.name,
                                              os.path.basename(
                                                  ap.analysis_dir))
            outputs.append(project_qc)
            outputs.append("%s.zip" % project_qc)
            outputs.append(os.path.join(project_qc,"qc_report.html"))
            outputs.append(os.path.join(project_qc,"qc"))
            # MultiQC output
            outputs.append("multiqc_report.%s.html" % project.name)
        # Additional QC for second fastq set in first project
        project_qc = "qc.extra_report.%s.%s" % (multi_fastqs_project.name,
                                                os.path.basename(
                                                    ap.analysis_dir))
        outputs.append(project_qc)
        outputs.append("%s.zip" % project_qc)
        outputs.append(os.path.join(project_qc,"qc.extra_report.html"))
        outputs.append(os.path.join(project_qc,"qc.extra"))
        # MultiQC output
        outputs.append("multiqc.extra_report.%s.html" %
                       multi_fastqs_project.name)
        for item in outputs:
            f = os.path.join(publication_dir,
                             "160621_K00879_0087_000000000-AGEW9_analysis",
                             item)
            self.assertTrue(os.path.exists(f),"Missing %s" % f)

    def test_publish_qc_with_project_missing_qc(self):
        """publish_qc: raises exception if project has missing QC
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create()
        ap = AutoProcess(mockdir.dirn)
        # Add processing report
        UpdateAnalysisDir(ap).add_processing_report()
        # Add QC outputs for subset of projects
        projects = ap.get_analysis_projects()[1:]
        for project in projects:
            UpdateAnalysisProject(project).add_qc_outputs()
        # Make a mock publication area
        publication_dir = os.path.join(self.dirn,'QC')
        os.mkdir(publication_dir)
        # Publish
        self.assertRaises(Exception,
                          ap.publish_qc,
                          location=publication_dir)

    def test_publish_qc_ignore_project_missing_qc(self):
        """publish_qc: skip project with missing QC
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create()
        ap = AutoProcess(mockdir.dirn)
        # Add processing report
        UpdateAnalysisDir(ap).add_processing_report()
        # Add QC outputs for subset of projects
        projects = ap.get_analysis_projects()
        missing_project = projects[0]
        projects = projects[1:]
        for project in projects:
            UpdateAnalysisProject(project).add_qc_outputs()
        # Make a mock publication area
        publication_dir = os.path.join(self.dirn,'QC')
        os.mkdir(publication_dir)
        # Publish
        ap.publish_qc(location=publication_dir,
                      ignore_missing_qc=True)
        # Check outputs
        outputs = ["index.html",
                   "processing_qc.html"]
        for project in projects:
            # Standard QC outputs
            project_qc = "qc_report.%s.%s" % (project.name,
                                              os.path.basename(
                                                  ap.analysis_dir))
            outputs.append(project_qc)
            outputs.append("%s.zip" % project_qc)
            outputs.append(os.path.join(project_qc,"qc_report.html"))
            outputs.append(os.path.join(project_qc,"qc"))
            # MultiQC output
            outputs.append("multiqc_report.%s.html" % project.name)
        for item in outputs:
            f = os.path.join(publication_dir,
                             "160621_K00879_0087_000000000-AGEW9_analysis",
                             item)
            self.assertTrue(os.path.exists(f),"Missing %s" % f)
        # Check that missing project wasn't copied
        self.assertFalse(os.path.exists(
            os.path.join(publication_dir,
                         "160621_K00879_0087_000000000-AGEW9_analysis",
                         "qc_report.%s.%s" % (missing_project.name,
                                              os.path.basename(
                                                  ap.analysis_dir)))))

    def test_publish_qc_subset_of_projects(self):
        """publish_qc: only publish subset of projects
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create()
        ap = AutoProcess(mockdir.dirn)
        # Add processing report
        UpdateAnalysisDir(ap).add_processing_report()
        # Add QC outputs for subset of projects
        projects = ap.get_analysis_projects()
        missing_projects = projects[1:]
        projects = projects[0:1]
        for project in ap.get_analysis_projects():
            UpdateAnalysisProject(project).add_qc_outputs()
        # Make a mock publication area
        publication_dir = os.path.join(self.dirn,'QC')
        os.mkdir(publication_dir)
        # Publish
        ap.publish_qc(location=publication_dir,
                      projects="AB*")
        # Check outputs
        outputs = ["index.html",
                   "processing_qc.html"]
        for project in projects:
            # Standard QC outputs
            project_qc = "qc_report.%s.%s" % (project.name,
                                              os.path.basename(
                                                  ap.analysis_dir))
            outputs.append(project_qc)
            outputs.append("%s.zip" % project_qc)
            outputs.append(os.path.join(project_qc,"qc_report.html"))
            outputs.append(os.path.join(project_qc,"qc"))
            # MultiQC output
            outputs.append("multiqc_report.%s.html" % project.name)
        for item in outputs:
            f = os.path.join(publication_dir,
                             "160621_K00879_0087_000000000-AGEW9_analysis",
                             item)
            self.assertTrue(os.path.exists(f),"Missing %s" % f)
        # Check that missing projects weren't copied
        for project in missing_projects:
            self.assertFalse(os.path.exists(
                os.path.join(publication_dir,
                             "160621_K00879_0087_000000000-AGEW9_analysis",
                             "qc_report.%s.%s" % (project.name,
                                                  os.path.basename(
                                                      ap.analysis_dir)))),
                             "%s exists, but shouldn't" % project.name)

    def test_publish_qc_with_icell8_outputs(self):
        """publish_qc: project with ICell8 QC outputs
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create()
        ap = AutoProcess(mockdir.dirn)
        # Add processing report and QC outputs
        UpdateAnalysisDir(ap).add_processing_report()
        projects = ap.get_analysis_projects()
        for project in projects:
            UpdateAnalysisProject(project).add_qc_outputs()
        # Add ICell8 report for one project
        icell8_project = projects[0]
        UpdateAnalysisProject(icell8_project).add_icell8_outputs()
        # Make a mock publication area
        publication_dir = os.path.join(self.dirn,'QC')
        os.mkdir(publication_dir)
        # Publish
        ap.publish_qc(location=publication_dir)
        # Check outputs
        outputs = ["index.html",
                   "processing_qc.html"]
        for project in ap.get_analysis_projects():
            # Standard QC outputs
            project_qc = "qc_report.%s.%s" % (project.name,
                                              os.path.basename(
                                                  ap.analysis_dir))
            outputs.append(project_qc)
            outputs.append("%s.zip" % project_qc)
            outputs.append(os.path.join(project_qc,"qc_report.html"))
            outputs.append(os.path.join(project_qc,"qc"))
            # MultiQC output
            outputs.append("multiqc_report.%s.html" % project.name)
        # ICell8 outputs
        icell8_dir = "icell8_processing.%s.%s" % (icell8_project.name,
                                                  os.path.basename(
                                                      ap.analysis_dir))
        outputs.append(icell8_dir)
        outputs.append("%s.zip" % icell8_dir)
        outputs.append(os.path.join(icell8_dir,"icell8_processing_data"))
        outputs.append(os.path.join(icell8_dir,"icell8_processing.html"))
        outputs.append(os.path.join(icell8_dir,"stats"))
        # Do checks
        for item in outputs:
            f = os.path.join(publication_dir,
                             "160621_K00879_0087_000000000-AGEW9_analysis",
                             item)
            self.assertTrue(os.path.exists(f),"Missing %s" % f)

    def test_publish_qc_with_cellranger_qc(self):
        """publish_qc: publish cellranger QC output (whole flowcell)
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        ap = AutoProcess(mockdir.dirn)
        # Add processing and cellranger QC reports
        UpdateAnalysisDir(ap).add_processing_report()
        UpdateAnalysisDir(ap).add_cellranger_qc_output()
        # Make a mock publication area
        publication_dir = os.path.join(self.dirn,'QC')
        os.mkdir(publication_dir)
        # Publish
        ap.publish_qc(location=publication_dir)
        # Check outputs
        outputs = ["index.html",
                   "processing_qc.html",
                   "cellranger_qc_summary.html"]
        # Do checks
        for item in outputs:
            f = os.path.join(publication_dir,
                             "160621_K00879_0087_000000000-AGEW9_analysis",
                             item)
            self.assertTrue(os.path.exists(f),"Missing %s" % f)

    def test_publish_qc_with_cellranger_qc_lanes_subset(self):
        """publish_qc: publish cellranger QC output (subset of lanes)
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        ap = AutoProcess(mockdir.dirn)
        # Add processing and cellranger QC reports
        UpdateAnalysisDir(ap).add_processing_report()
        UpdateAnalysisDir(ap).add_cellranger_qc_output(lanes="45")
        # Make a mock publication area
        publication_dir = os.path.join(self.dirn,'QC')
        os.mkdir(publication_dir)
        # Publish
        ap.publish_qc(location=publication_dir)
        # Check outputs
        outputs = ["index.html",
                   "processing_qc.html",
                   "cellranger_qc_summary_45.html"]
        # Do checks
        for item in outputs:
            f = os.path.join(publication_dir,
                             "160621_K00879_0087_000000000-AGEW9_analysis",
                             item)
            self.assertTrue(os.path.exists(f),"Missing %s" % f)

    def test_publish_qc_with_cellranger_qc_multiple_lanes_subsets(self):
        """publish_qc: publish cellranger QC output (multiple subsets of lanes)
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        ap = AutoProcess(mockdir.dirn)
        # Add processing and cellranger QC reports
        UpdateAnalysisDir(ap).add_processing_report()
        UpdateAnalysisDir(ap).add_cellranger_qc_output(lanes="45")
        UpdateAnalysisDir(ap).add_cellranger_qc_output(lanes="78")
        # Make a mock publication area
        publication_dir = os.path.join(self.dirn,'QC')
        os.mkdir(publication_dir)
        # Publish
        ap.publish_qc(location=publication_dir)
        # Check outputs
        outputs = ["index.html",
                   "processing_qc.html",
                   "cellranger_qc_summary_45.html",
                   "cellranger_qc_summary_78.html"]
        # Do checks
        for item in outputs:
            f = os.path.join(publication_dir,
                             "160621_K00879_0087_000000000-AGEW9_analysis",
                             item)
            self.assertTrue(os.path.exists(f),"Missing %s" % f)

    def test_publish_qc_with_cellranger_counts(self):
        """publish_qc: project with cellranger count output
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create()
        ap = AutoProcess(mockdir.dirn)
        # Add processing and cellranger QC reports
        UpdateAnalysisDir(ap).add_processing_report()
        UpdateAnalysisDir(ap).add_cellranger_qc_output()
        # Add QC outputs
        projects = ap.get_analysis_projects()
        for project in projects:
            UpdateAnalysisProject(project).add_qc_outputs()
        # Add cellranger count output for one project
        tenxgenomics_project = projects[0]
        UpdateAnalysisProject(tenxgenomics_project).add_cellranger_count_outputs()
        # Make a mock publication area
        publication_dir = os.path.join(self.dirn,'QC')
        os.mkdir(publication_dir)
        # Publish
        ap.publish_qc(location=publication_dir)
        # Check outputs
        outputs = ["index.html",
                   "processing_qc.html",
                   "cellranger_qc_summary.html"]
        for project in ap.get_analysis_projects():
            # Standard QC outputs
            project_qc = "qc_report.%s.%s" % (project.name,
                                              os.path.basename(
                                                  ap.analysis_dir))
            outputs.append(project_qc)
            outputs.append("%s.zip" % project_qc)
            outputs.append(os.path.join(project_qc,"qc_report.html"))
            outputs.append(os.path.join(project_qc,"qc"))
            # MultiQC output
            outputs.append("multiqc_report.%s.html" % project.name)
        # Cellranger count outputs
        cellranger_count_dir = "cellranger_count_report.%s.%s" % (
            tenxgenomics_project.name,
            os.path.basename(ap.analysis_dir))
        outputs.append(cellranger_count_dir)
        outputs.append("%s.zip" % cellranger_count_dir)
        outputs.append(os.path.join(cellranger_count_dir,
                                    "cellranger_count_report.html"))
        outputs.append(os.path.join(cellranger_count_dir,"cellranger_count"))
        for sample in tenxgenomics_project.samples:
            outputs.append(os.path.join(cellranger_count_dir,
                                        "cellranger_count",
                                        sample.name,
                                        "outs",
                                        "web_summary.html"))
        # Do checks
        for item in outputs:
            f = os.path.join(publication_dir,
                             "160621_K00879_0087_000000000-AGEW9_analysis",
                             item)
            self.assertTrue(os.path.exists(f),"Missing %s" % f)
        # Do checks
        for item in outputs:
            f = os.path.join(publication_dir,
                             "160621_K00879_0087_000000000-AGEW9_analysis",
                             item)
            self.assertTrue(os.path.exists(f),"Missing %s" % f)

fastq_reads_r1 = (
    "@HISEQ:1:000000000-A2Y1L:1:1101:19264:2433 1:N:0:AGATCGC",
    "AGATAGCCGA","+","?????BBB@B",
    "@HISEQ:1:000000000-A2Y1L:1:1101:18667:2435 1:N:0:AGATCGC",
    "ATATATTCAT","+","?????BBBDD",
    "@HISEQ:1:000000000-A2Y1L:1:1101:17523:2436 1:N:0:AGATCGC",
    "CATCACTACC","+","?<,<?BBBBB"
)
fastq_reads_r2 = (
    "@HISEQ:1:000000000-A2Y1L:1:1101:19264:2433 2:N:0:AGATCGC",
    "GCCGATATGC","+","??A??ABBDD",
    "@HISEQ:1:000000000-A2Y1L:1:1101:18667:2435 2:N:0:AGATCGC",
    "GATGACATCA","+","?????BBBDD",
    "@HISEQ:1:000000000-A2Y1L:1:1101:17523:2436 2:N:0:AGATCGC",
    "GAATATAGAA","+","??AAABBBDD"
)

class TestAutoProcessMergeFastqDirs(unittest.TestCase):
    """
    Tests for AutoProcess.merge_fastq_dirs
    """
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestAutoProcessMergeFastqDirs')
        # Store original location so we can get back at the end
        self.pwd = os.getcwd()
        # Move to working dir
        os.chdir(self.dirn)
        # Placeholders for test objects
        self.ap = None

    def tearDown(self):
        # Delete autoprocessor object
        if self.ap is not None:
            del(self.ap)
        # Return to original dir
        os.chdir(self.pwd)
        # Remove the temporary test directory
        shutil.rmtree(self.dirn)

    def _setup_bcl2fastq2_no_lane_splitting(self):
        # Create mock bcl2fastq2 dir structure with no lane splitting
        lanes = [1,2,3,4]
        mockdir1 = MockAnalysisDir("161209_K123_0001_BLAH",
                                   "hiseq",
                                   unaligned_dir='bcl2fastq.AB',
                                   fmt='bcl2fastq2',
                                   paired_end=True,
                                   no_lane_splitting=True,
                                   lanes=lanes,
                                   top_dir=self.dirn)
        mockdir1.add_fastq_batch('AB','AB1','AB1_S1',lanes=lanes)
        mockdir1.add_fastq_batch('AB','AB2','AB2_S2',lanes=lanes)
        m1 = mockdir1.create()
        print m1
        mockdir2 = MockAnalysisDir("161209_K123_0001_BLAH2",
                                   "hiseq",
                                   unaligned_dir='bcl2fastq.CDE',
                                   fmt='bcl2fastq2',
                                   paired_end=True,
                                   no_lane_splitting=True,
                                   lanes=lanes,
                                   top_dir=self.dirn)
        mockdir2.add_fastq_batch('CDE','CDE3','CDE3_S3',lanes=lanes)
        mockdir2.add_fastq_batch('CDE','CDE4','CDE4_S4',lanes=lanes)
        m2 = mockdir2.create()
        # Move the second 'unaligned' dir into the first analysis dir
        shutil.move(os.path.join(m2,'bcl2fastq.CDE'),m1)
        # Remove unwanted project dirs and files
        shutil.rmtree(os.path.join(m1,'AB'))
        shutil.rmtree(os.path.join(m1,'undetermined'))
        # Add content to the undetermined fastqs
        for n,dirn in enumerate(('bcl2fastq.AB','bcl2fastq.CDE')):
            undetermined_r1 = os.path.join(m1,dirn,'Undetermined_S0_R1_001.fastq.gz')
            undetermined_r2 = os.path.join(m1,dirn,'Undetermined_S0_R2_001.fastq.gz')
            with gzip.GzipFile(undetermined_r1,'wb') as fq:
                for i in xrange(4):
                    fq.write("%s\n" % fastq_reads_r1[n*4+i])
            with gzip.GzipFile(undetermined_r2,'wb') as fq:
                for i in xrange(4):
                    fq.write("%s\n" % fastq_reads_r2[n*4+i])
        print m2
        return m1

    def _setup_bcl2fastq2(self):
        # Create a mock bcl2fastq dir structure
        lanes = [1,2,3,4]
        mockdir1 = MockAnalysisDir("161209_K123_0002_BLAH",
                                   "hiseq",
                                   unaligned_dir='bcl2fastq.lanes1-2',
                                   fmt='bcl2fastq2',
                                   paired_end=True,
                                   no_lane_splitting=False,
                                   lanes=lanes,
                                   top_dir=self.dirn)
        mockdir1.add_fastq_batch('AB','AB1','AB1_S1',lanes=[1,2])
        mockdir1.add_fastq_batch('AB','AB2','AB2_S2',lanes=[1,2])
        m1 = mockdir1.create()
        print m1
        mockdir2 = MockAnalysisDir("161209_K123_0002_BLAH2",
                                   "hiseq",
                                   unaligned_dir='bcl2fastq.lanes3-4',
                                   fmt='bcl2fastq2',
                                   paired_end=True,
                                   no_lane_splitting=False,
                                   lanes=lanes,
                                   top_dir=self.dirn)
        mockdir2.add_fastq_batch('CDE','CDE3','CDE3_S3',lanes=[3,4])
        mockdir2.add_fastq_batch('CDE','CDE4','CDE4_S4',lanes=[3,4])
        m2 = mockdir2.create()
        # Move the second 'unaligned' dir into the first analysis dir
        shutil.move(os.path.join(m2,'bcl2fastq.lanes3-4'),m1)
        # Remove unwanted project dirs and files
        shutil.rmtree(os.path.join(m1,'AB'))
        shutil.rmtree(os.path.join(m1,'undetermined'))
        print m2
        return m1

    def _setup_casava(self):
        # Create a mock casava dir structure
        lanes = [1,2,3,4]
        mockdir1 = MockAnalysisDir("161209_K123_0003_BLAH",
                                   "hiseq",
                                   unaligned_dir='bcl2fastq.lanes1-2',
                                   fmt='casava',
                                   paired_end=True,
                                   no_lane_splitting=False,
                                   lanes=lanes,
                                   top_dir=self.dirn)
        mockdir1.add_fastq_batch('AB','AB1','AB1_GCCAAT',lanes=[1,2])
        mockdir1.add_fastq_batch('AB','AB2','AB2_AGTCAA',lanes=[1,2])
        m1 = mockdir1.create()
        print m1
        mockdir2 = MockAnalysisDir("161209_K123_0003_BLAH2",
                                   "hiseq",
                                   unaligned_dir='bcl2fastq.lanes3-4',
                                   fmt='casava',
                                   paired_end=True,
                                   no_lane_splitting=False,
                                   lanes=lanes,
                                   top_dir=self.dirn)
        mockdir2.add_fastq_batch('CDE','CDE3','CDE3_GCCAAT',lanes=[3,4])
        mockdir2.add_fastq_batch('CDE','CDE4','CDE4_AGTCAA',lanes=[3,4])
        m2 = mockdir2.create()
        # Move the second 'unaligned' dir into the first analysis dir
        shutil.move(os.path.join(m2,'bcl2fastq.lanes3-4'),m1)
        # Remove unwanted project dirs and files
        shutil.rmtree(os.path.join(m1,'AB'))
        shutil.rmtree(os.path.join(m1,'undetermined'))
        print m2
        return m1

    def _assert_dir_exists(self,path):
        self.assertTrue(os.path.isdir(path),
                        "Missing dir '%s'" % path)
        
    def _assert_dir_doesnt_exist(self,path):
        self.assertFalse(os.path.isdir(path),
                         "Dir '%s' shouldn't exist" % path)

    def _assert_file_exists(self,path):
        self.assertTrue(os.path.isfile(path),
                        "Missing file '%s'" % path)

    def _assert_file_doesnt_exist(self,path):
        self.assertFalse(os.path.isfile(path),
                        "File '%s' shouldn't exist" % path)

    def test_bcl2fastq2_no_lane_splitting_new_output_dir(self):
        """
        AutoProcess.merge_fastq_dirs: bcl2fastq v2 output with --no-lane-splitting, new output dir
        """
        analysis_dir = self._setup_bcl2fastq2_no_lane_splitting()
        # Merge the unaligned dirs
        self.ap = AutoProcess(analysis_dir)
        self.ap.merge_fastq_dirs("bcl2fastq.AB",output_dir="bcl2fastq")
        # Check outputs
        self._assert_dir_exists(os.path.join(analysis_dir,'save.bcl2fastq.AB'))
        self._assert_dir_exists(os.path.join(analysis_dir,'save.bcl2fastq.CDE'))
        self._assert_dir_exists(os.path.join(analysis_dir,'bcl2fastq'))
        self._assert_dir_doesnt_exist(os.path.join(analysis_dir,'bcl2fastq.AB'))
        self._assert_dir_doesnt_exist(os.path.join(analysis_dir,'bcl2fastq.CDE'))
        for f in ('AB/AB1_S1_R1_001.fastq.gz',
                  'AB/AB1_S1_R2_001.fastq.gz',
                  'AB/AB2_S2_R1_001.fastq.gz',
                  'AB/AB2_S2_R2_001.fastq.gz',
                  'CDE/CDE3_S3_R1_001.fastq.gz',
                  'CDE/CDE3_S3_R2_001.fastq.gz',
                  'CDE/CDE4_S4_R1_001.fastq.gz',
                  'CDE/CDE4_S4_R2_001.fastq.gz',
                  'Undetermined_S0_R1_001.fastq.gz',
                  'Undetermined_S0_R2_001.fastq.gz',):
            self._assert_file_exists(os.path.join(analysis_dir,'bcl2fastq',f))
        # Check merge of undetermined fastqs
        undetermined_r1 = gzip.GzipFile(
            os.path.join(analysis_dir,'bcl2fastq','Undetermined_S0_R1_001.fastq.gz'),
            'rb').read()
        expected_r1 = '\n'.join(fastq_reads_r1[:8])+'\n'
        self.assertEqual(undetermined_r1,expected_r1)
        undetermined_r2 = gzip.GzipFile(
            os.path.join(analysis_dir,'bcl2fastq','Undetermined_S0_R2_001.fastq.gz'),
            'rb').read()
        expected_r2 = '\n'.join(fastq_reads_r2[:8])+'\n'
        self.assertEqual(undetermined_r2,expected_r2)
        # Check projects.info files
        self._assert_file_exists(os.path.join(analysis_dir,'save.projects.info'))
        self._assert_file_exists(os.path.join(analysis_dir,'projects.info'))
        projects_info = open(os.path.join(analysis_dir,'projects.info'),'r').read()
        expected = """#Project	Samples	User	Library	Organism	PI	Comments
AB	AB1,AB2	.	.	.	.	.
CDE	CDE3,CDE4	.	.	.	.	.
"""
        self.assertEqual(projects_info,expected)

    def test_bcl2fastq2_new_output_dir(self):
        """
        AutoProcess.merge_fastq_dirs: bcl2fastq v2 output, new output dir
        """
        analysis_dir = self._setup_bcl2fastq2()
        # Merge the unaligned dirs
        self.ap = AutoProcess(analysis_dir)
        self.ap.merge_fastq_dirs("bcl2fastq.lanes1-2",output_dir="bcl2fastq")
        # Check outputs
        self._assert_dir_exists(os.path.join(analysis_dir,'save.bcl2fastq.lanes1-2'))
        self._assert_dir_exists(os.path.join(analysis_dir,'save.bcl2fastq.lanes3-4'))
        self._assert_dir_exists(os.path.join(analysis_dir,'bcl2fastq'))
        self._assert_dir_doesnt_exist(os.path.join(analysis_dir,'bcl2fastq.lanes1-2'))
        self._assert_dir_doesnt_exist(os.path.join(analysis_dir,'bcl2fastq.lanes3-4'))
        for f in ('AB/AB1_S1_L001_R1_001.fastq.gz',
                  'AB/AB1_S1_L001_R2_001.fastq.gz',
                  'AB/AB2_S2_L001_R1_001.fastq.gz',
                  'AB/AB2_S2_L001_R2_001.fastq.gz',
                  'AB/AB1_S1_L002_R1_001.fastq.gz',
                  'AB/AB1_S1_L002_R2_001.fastq.gz',
                  'AB/AB2_S2_L002_R1_001.fastq.gz',
                  'AB/AB2_S2_L002_R2_001.fastq.gz',
                  'CDE/CDE3_S3_L003_R1_001.fastq.gz',
                  'CDE/CDE3_S3_L003_R2_001.fastq.gz',
                  'CDE/CDE4_S4_L003_R1_001.fastq.gz',
                  'CDE/CDE4_S4_L003_R2_001.fastq.gz',
                  'CDE/CDE3_S3_L004_R1_001.fastq.gz',
                  'CDE/CDE3_S3_L004_R2_001.fastq.gz',
                  'CDE/CDE4_S4_L004_R1_001.fastq.gz',
                  'CDE/CDE4_S4_L004_R2_001.fastq.gz',
                  'Undetermined_S0_L001_R1_001.fastq.gz',
                  'Undetermined_S0_L001_R2_001.fastq.gz',
                  'Undetermined_S0_L002_R1_001.fastq.gz',
                  'Undetermined_S0_L002_R2_001.fastq.gz',
                  'Undetermined_S0_L003_R1_001.fastq.gz',
                  'Undetermined_S0_L003_R2_001.fastq.gz',
                  'Undetermined_S0_L004_R1_001.fastq.gz',
                  'Undetermined_S0_L004_R2_001.fastq.gz'):
            self._assert_file_exists(os.path.join(analysis_dir,'bcl2fastq',f))
        # Check projects.info files
        self._assert_file_exists(os.path.join(analysis_dir,'save.projects.info'))
        self._assert_file_exists(os.path.join(analysis_dir,'projects.info'))
        projects_info = open(os.path.join(analysis_dir,'projects.info'),'r').read()
        expected = """#Project	Samples	User	Library	Organism	PI	Comments
AB	AB1,AB2	.	.	.	.	.
CDE	CDE3,CDE4	.	.	.	.	.
"""
        self.assertEqual(projects_info,expected)

    def test_casava_new_output_dir(self):
        """
        AutoProcess.merge_fastq_dirs: casava/bcl2fastq v1.8.* output, new output dir
        """
        analysis_dir = self._setup_casava()
        # Merge the unaligned dirs
        self.ap = AutoProcess(analysis_dir)
        self.ap.merge_fastq_dirs("bcl2fastq.lanes1-2",output_dir="bcl2fastq")
        # Check outputs
        self._assert_dir_exists(os.path.join(analysis_dir,'save.bcl2fastq.lanes1-2'))
        self._assert_dir_exists(os.path.join(analysis_dir,'save.bcl2fastq.lanes3-4'))
        self._assert_dir_exists(os.path.join(analysis_dir,'bcl2fastq'))
        self._assert_dir_doesnt_exist(os.path.join(analysis_dir,'bcl2fastq.lanes1-2'))
        self._assert_dir_doesnt_exist(os.path.join(analysis_dir,'bcl2fastq.lanes3-4'))
        for f in ('Project_AB/Sample_AB1/AB1_GCCAAT_L001_R1_001.fastq.gz',
                  'Project_AB/Sample_AB1/AB1_GCCAAT_L001_R2_001.fastq.gz',
                  'Project_AB/Sample_AB2/AB2_AGTCAA_L001_R1_001.fastq.gz',
                  'Project_AB/Sample_AB2/AB2_AGTCAA_L001_R2_001.fastq.gz',
                  'Project_AB/Sample_AB1/AB1_GCCAAT_L002_R1_001.fastq.gz',
                  'Project_AB/Sample_AB1/AB1_GCCAAT_L002_R2_001.fastq.gz',
                  'Project_AB/Sample_AB2/AB2_AGTCAA_L002_R1_001.fastq.gz',
                  'Project_AB/Sample_AB2/AB2_AGTCAA_L002_R2_001.fastq.gz',
                  'Project_CDE/Sample_CDE3/CDE3_GCCAAT_L003_R1_001.fastq.gz',
                  'Project_CDE/Sample_CDE3/CDE3_GCCAAT_L003_R2_001.fastq.gz',
                  'Project_CDE/Sample_CDE4/CDE4_AGTCAA_L003_R1_001.fastq.gz',
                  'Project_CDE/Sample_CDE4/CDE4_AGTCAA_L003_R2_001.fastq.gz',
                  'Project_CDE/Sample_CDE3/CDE3_GCCAAT_L004_R1_001.fastq.gz',
                  'Project_CDE/Sample_CDE3/CDE3_GCCAAT_L004_R2_001.fastq.gz',
                  'Project_CDE/Sample_CDE4/CDE4_AGTCAA_L004_R1_001.fastq.gz',
                  'Project_CDE/Sample_CDE4/CDE4_AGTCAA_L004_R2_001.fastq.gz',
                  'Undetermined_indices/Sample_lane1/lane1_Undetermined_L001_R1_001.fastq.gz',
                  'Undetermined_indices/Sample_lane1/lane1_Undetermined_L001_R2_001.fastq.gz',
                  'Undetermined_indices/Sample_lane2/lane2_Undetermined_L002_R1_001.fastq.gz',
                  'Undetermined_indices/Sample_lane2/lane2_Undetermined_L002_R2_001.fastq.gz',
                  'Undetermined_indices/Sample_lane3/lane3_Undetermined_L003_R1_001.fastq.gz',
                  'Undetermined_indices/Sample_lane3/lane3_Undetermined_L003_R2_001.fastq.gz',
                  'Undetermined_indices/Sample_lane4/lane4_Undetermined_L004_R1_001.fastq.gz',
                  'Undetermined_indices/Sample_lane4/lane4_Undetermined_L004_R2_001.fastq.gz'):
            self._assert_file_exists(os.path.join(analysis_dir,'bcl2fastq',f))
        # Check projects.info files
        self._assert_file_exists(os.path.join(analysis_dir,'save.projects.info'))
        self._assert_file_exists(os.path.join(analysis_dir,'projects.info'))
        projects_info = open(os.path.join(analysis_dir,'projects.info'),'r').read()
        expected = """#Project	Samples	User	Library	Organism	PI	Comments
AB	AB1,AB2	.	.	.	.	.
CDE	CDE3,CDE4	.	.	.	.	.
"""
        self.assertEqual(projects_info,expected)

    def test_bcl2fastq2_no_lane_splitting(self):
        """
        AutoProcess.merge_fastq_dirs: bcl2fastq v2 output with --no-lane-splitting
        """
        analysis_dir = self._setup_bcl2fastq2_no_lane_splitting()
        # Merge the unaligned dirs
        self.ap = AutoProcess(analysis_dir)
        self.ap.merge_fastq_dirs("bcl2fastq.AB")
        # Check outputs
        self._assert_dir_exists(os.path.join(analysis_dir,'save.bcl2fastq.AB'))
        self._assert_dir_exists(os.path.join(analysis_dir,'save.bcl2fastq.CDE'))
        self._assert_dir_exists(os.path.join(analysis_dir,'bcl2fastq.AB'))
        self._assert_dir_doesnt_exist(os.path.join(analysis_dir,'bcl2fastq.CDE'))
        for f in ('AB/AB1_S1_R1_001.fastq.gz',
                  'AB/AB1_S1_R2_001.fastq.gz',
                  'AB/AB2_S2_R1_001.fastq.gz',
                  'AB/AB2_S2_R2_001.fastq.gz',
                  'CDE/CDE3_S3_R1_001.fastq.gz',
                  'CDE/CDE3_S3_R2_001.fastq.gz',
                  'CDE/CDE4_S4_R1_001.fastq.gz',
                  'CDE/CDE4_S4_R2_001.fastq.gz',
                  'Undetermined_S0_R1_001.fastq.gz',
                  'Undetermined_S0_R2_001.fastq.gz',):
            self._assert_file_exists(os.path.join(analysis_dir,'bcl2fastq.AB',f))
        # Check merge of undetermined fastqs
        undetermined_r1 = gzip.GzipFile(
            os.path.join(analysis_dir,'bcl2fastq.AB','Undetermined_S0_R1_001.fastq.gz'),
            'rb').read()
        expected_r1 = '\n'.join(fastq_reads_r1[:8])+'\n'
        self.assertEqual(undetermined_r1,expected_r1)
        undetermined_r2 = gzip.GzipFile(
            os.path.join(analysis_dir,'bcl2fastq.AB','Undetermined_S0_R2_001.fastq.gz'),
            'rb').read()
        expected_r2 = '\n'.join(fastq_reads_r2[:8])+'\n'
        self.assertEqual(undetermined_r2,expected_r2)
        # Check projects.info files
        self._assert_file_exists(os.path.join(analysis_dir,'save.projects.info'))
        self._assert_file_exists(os.path.join(analysis_dir,'projects.info'))
        projects_info = open(os.path.join(analysis_dir,'projects.info'),'r').read()
        expected = """#Project	Samples	User	Library	Organism	PI	Comments
AB	AB1,AB2	.	.	.	.	.
CDE	CDE3,CDE4	.	.	.	.	.
"""
        self.assertEqual(projects_info,expected)

    def test_bcl2fastq2(self):
        """
        AutoProcess.merge_fastq_dirs: bcl2fastq v2 output
        """
        analysis_dir = self._setup_bcl2fastq2()
        # Merge the unaligned dirs
        self.ap = AutoProcess(analysis_dir)
        self.ap.merge_fastq_dirs("bcl2fastq.lanes1-2")
        # Check outputs
        self._assert_dir_exists(os.path.join(analysis_dir,'save.bcl2fastq.lanes1-2'))
        self._assert_dir_exists(os.path.join(analysis_dir,'save.bcl2fastq.lanes3-4'))
        self._assert_dir_exists(os.path.join(analysis_dir,'bcl2fastq.lanes1-2'))
        self._assert_dir_doesnt_exist(os.path.join(analysis_dir,'bcl2fastq.lanes3-4'))
        for f in ('AB/AB1_S1_L001_R1_001.fastq.gz',
                  'AB/AB1_S1_L001_R2_001.fastq.gz',
                  'AB/AB2_S2_L001_R1_001.fastq.gz',
                  'AB/AB2_S2_L001_R2_001.fastq.gz',
                  'AB/AB1_S1_L002_R1_001.fastq.gz',
                  'AB/AB1_S1_L002_R2_001.fastq.gz',
                  'AB/AB2_S2_L002_R1_001.fastq.gz',
                  'AB/AB2_S2_L002_R2_001.fastq.gz',
                  'CDE/CDE3_S3_L003_R1_001.fastq.gz',
                  'CDE/CDE3_S3_L003_R2_001.fastq.gz',
                  'CDE/CDE4_S4_L003_R1_001.fastq.gz',
                  'CDE/CDE4_S4_L003_R2_001.fastq.gz',
                  'CDE/CDE3_S3_L004_R1_001.fastq.gz',
                  'CDE/CDE3_S3_L004_R2_001.fastq.gz',
                  'CDE/CDE4_S4_L004_R1_001.fastq.gz',
                  'CDE/CDE4_S4_L004_R2_001.fastq.gz',
                  'Undetermined_S0_L001_R1_001.fastq.gz',
                  'Undetermined_S0_L001_R2_001.fastq.gz',
                  'Undetermined_S0_L002_R1_001.fastq.gz',
                  'Undetermined_S0_L002_R2_001.fastq.gz',
                  'Undetermined_S0_L003_R1_001.fastq.gz',
                  'Undetermined_S0_L003_R2_001.fastq.gz',
                  'Undetermined_S0_L004_R1_001.fastq.gz',
                  'Undetermined_S0_L004_R2_001.fastq.gz'):
            self._assert_file_exists(os.path.join(analysis_dir,'bcl2fastq.lanes1-2',f))
        # Check projects.info files
        self._assert_file_exists(os.path.join(analysis_dir,'save.projects.info'))
        self._assert_file_exists(os.path.join(analysis_dir,'projects.info'))
        projects_info = open(os.path.join(analysis_dir,'projects.info'),'r').read()
        expected = """#Project	Samples	User	Library	Organism	PI	Comments
AB	AB1,AB2	.	.	.	.	.
CDE	CDE3,CDE4	.	.	.	.	.
"""
        self.assertEqual(projects_info,expected)

    def test_casava(self):
        """
        AutoProcess.merge_fastq_dirs: casava/bcl2fastq v1.8.* output
        """
        analysis_dir = self._setup_casava()
        # Merge the unaligned dirs
        self.ap = AutoProcess(analysis_dir)
        self.ap.merge_fastq_dirs("bcl2fastq.lanes1-2")
        # Check outputs
        self._assert_dir_exists(os.path.join(analysis_dir,'save.bcl2fastq.lanes1-2'))
        self._assert_dir_exists(os.path.join(analysis_dir,'save.bcl2fastq.lanes3-4'))
        self._assert_dir_exists(os.path.join(analysis_dir,'bcl2fastq.lanes1-2'))
        self._assert_dir_doesnt_exist(os.path.join(analysis_dir,'bcl2fastq.lanes3-4'))
        for f in ('Project_AB/Sample_AB1/AB1_GCCAAT_L001_R1_001.fastq.gz',
                  'Project_AB/Sample_AB1/AB1_GCCAAT_L001_R2_001.fastq.gz',
                  'Project_AB/Sample_AB2/AB2_AGTCAA_L001_R1_001.fastq.gz',
                  'Project_AB/Sample_AB2/AB2_AGTCAA_L001_R2_001.fastq.gz',
                  'Project_AB/Sample_AB1/AB1_GCCAAT_L002_R1_001.fastq.gz',
                  'Project_AB/Sample_AB1/AB1_GCCAAT_L002_R2_001.fastq.gz',
                  'Project_AB/Sample_AB2/AB2_AGTCAA_L002_R1_001.fastq.gz',
                  'Project_AB/Sample_AB2/AB2_AGTCAA_L002_R2_001.fastq.gz',
                  'Project_CDE/Sample_CDE3/CDE3_GCCAAT_L003_R1_001.fastq.gz',
                  'Project_CDE/Sample_CDE3/CDE3_GCCAAT_L003_R2_001.fastq.gz',
                  'Project_CDE/Sample_CDE4/CDE4_AGTCAA_L003_R1_001.fastq.gz',
                  'Project_CDE/Sample_CDE4/CDE4_AGTCAA_L003_R2_001.fastq.gz',
                  'Project_CDE/Sample_CDE3/CDE3_GCCAAT_L004_R1_001.fastq.gz',
                  'Project_CDE/Sample_CDE3/CDE3_GCCAAT_L004_R2_001.fastq.gz',
                  'Project_CDE/Sample_CDE4/CDE4_AGTCAA_L004_R1_001.fastq.gz',
                  'Project_CDE/Sample_CDE4/CDE4_AGTCAA_L004_R2_001.fastq.gz',
                  'Undetermined_indices/Sample_lane1/lane1_Undetermined_L001_R1_001.fastq.gz',
                  'Undetermined_indices/Sample_lane1/lane1_Undetermined_L001_R2_001.fastq.gz',
                  'Undetermined_indices/Sample_lane2/lane2_Undetermined_L002_R1_001.fastq.gz',
                  'Undetermined_indices/Sample_lane2/lane2_Undetermined_L002_R2_001.fastq.gz',
                  'Undetermined_indices/Sample_lane3/lane3_Undetermined_L003_R1_001.fastq.gz',
                  'Undetermined_indices/Sample_lane3/lane3_Undetermined_L003_R2_001.fastq.gz',
                  'Undetermined_indices/Sample_lane4/lane4_Undetermined_L004_R1_001.fastq.gz',
                  'Undetermined_indices/Sample_lane4/lane4_Undetermined_L004_R2_001.fastq.gz'):
            self._assert_file_exists(os.path.join(analysis_dir,'bcl2fastq.lanes1-2',f))
        # Check projects.info files
        self._assert_file_exists(os.path.join(analysis_dir,'save.projects.info'))
        self._assert_file_exists(os.path.join(analysis_dir,'projects.info'))
        projects_info = open(os.path.join(analysis_dir,'projects.info'),'r').read()
        expected = """#Project	Samples	User	Library	Organism	PI	Comments
AB	AB1,AB2	.	.	.	.	.
CDE	CDE3,CDE4	.	.	.	.	.
"""
        self.assertEqual(projects_info,expected)

    def test_bcl2fastq2_no_lane_splitting_dry_run(self):
        """
        AutoProcess.merge_fastq_dirs: dry run on bcl2fastq v2 output with --no-lane-splitting
        """
        analysis_dir = self._setup_bcl2fastq2_no_lane_splitting()
        # Merge the unaligned dirs
        self.ap = AutoProcess(analysis_dir)
        self.ap.merge_fastq_dirs("bcl2fastq.AB",dry_run=True)
        # Check outputs
        self._assert_dir_doesnt_exist(os.path.join(analysis_dir,'save.bcl2fastq.AB'))
        self._assert_dir_doesnt_exist(os.path.join(analysis_dir,'save.bcl2fastq.CDE'))
        self._assert_dir_exists(os.path.join(analysis_dir,'bcl2fastq.AB'))
        self._assert_dir_exists(os.path.join(analysis_dir,'bcl2fastq.CDE'))
        # Check projects.info files
        self._assert_file_doesnt_exist(os.path.join(analysis_dir,'save.projects.info'))
        self._assert_file_exists(os.path.join(analysis_dir,'projects.info'))

    def test_bcl2fastq2_dry_run(self):
        """
        AutoProcess.merge_fastq_dirs: dry run on bcl2fastq v2 output
        """
        analysis_dir = self._setup_bcl2fastq2()
        # Merge the unaligned dirs
        self.ap = AutoProcess(analysis_dir)
        self.ap.merge_fastq_dirs("bcl2fastq.lanes1-2",dry_run=True)
        # Check outputs
        self._assert_dir_doesnt_exist(os.path.join(analysis_dir,'save.bcl2fastq.lanes1-2'))
        self._assert_dir_doesnt_exist(os.path.join(analysis_dir,'save.bcl2fastq.lanes3-4'))
        self._assert_dir_exists(os.path.join(analysis_dir,'bcl2fastq.lanes1-2'))
        self._assert_dir_exists(os.path.join(analysis_dir,'bcl2fastq.lanes3-4'))
        # Check projects.info files
        self._assert_file_doesnt_exist(os.path.join(analysis_dir,'save.projects.info'))
        self._assert_file_exists(os.path.join(analysis_dir,'projects.info'))

    def test_casava_dry_run(self):
        """
        AutoProcess.merge_fastq_dirs: dry run on casava/bcl2fastq v1.8.* output
        """
        analysis_dir = self._setup_casava()
        # Merge the unaligned dirs
        self.ap = AutoProcess(analysis_dir)
        self.ap.merge_fastq_dirs("bcl2fastq.lanes1-2",dry_run=True)
        # Check outputs
        self._assert_dir_doesnt_exist(os.path.join(analysis_dir,'save.bcl2fastq.lanes1-2'))
        self._assert_dir_doesnt_exist(os.path.join(analysis_dir,'save.bcl2fastq.lanes3-4'))
        self._assert_dir_exists(os.path.join(analysis_dir,'bcl2fastq.lanes1-2'))
        self._assert_dir_exists(os.path.join(analysis_dir,'bcl2fastq.lanes3-4'))

class TestAutoProcessImportProject(unittest.TestCase):
    """Tests for AutoProcess.import_project

    """
    def setUp(self):
        self.dirn = tempfile.mkdtemp(suffix='TestAutoProcessImportProject')
        # Make a mock project
        project_dir = os.path.join(self.dirn,'NewProj')
        os.mkdir(project_dir)
        os.mkdir(os.path.join(project_dir,'fastqs'))
        for fq in ('NP01_S1_R1_001.fastq.gz','NP01_S1_R1_001.fastq.gz'):
            open(os.path.join(project_dir,'fastqs',fq),'w').write('')
        open(os.path.join(project_dir,'README.info'),'w').write(
            """Run\t160622_NB5001234_0011_ABCDE5AFXX
Platform\tnextseq
User\tPeter Briggs
PI\tAnne Cleaves
Organism\tHuman
Library type\tRNA-seq
Paired_end\tY
Samples\t1 sample (NP01)
Comments\t1% PhiX spike in
""")
        self.new_project_dir = project_dir

    def tearDown(self):
        # Remove the temporary test directory
        shutil.rmtree(self.dirn)

    def test_import_project(self):
        """Check AutoProcess.import_project imports a project
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_M00879_0087_000000000-AGEW9',
            'miseq',
            top_dir=self.dirn)
        mockdir.create()
        # Check that the project is not currently present
        ap = AutoProcess(mockdir.dirn)
        self.assertFalse('NewProj' in [p.name
                                       for p in ap.get_analysis_projects()])
        self.assertFalse('NewProj' in [p.name
                                       for p in ap.get_analysis_projects_from_dirs()])
        self.assertFalse(os.path.exists(os.path.join(ap.analysis_dir,'NewProj')))
        # Import the project
        ap.import_project(self.new_project_dir)
        self.assertTrue('NewProj' in [p.name
                                      for p in ap.get_analysis_projects()])
        self.assertTrue('NewProj' in [p.name
                                      for p in ap.get_analysis_projects_from_dirs()])
        self.assertTrue(os.path.exists(os.path.join(ap.analysis_dir,'NewProj')))
        # Verify via fresh AutoProcess object
        ap2 = AutoProcess(mockdir.dirn)
        self.assertTrue('NewProj' in [p.name
                                      for p in ap2.get_analysis_projects()])
        self.assertTrue('NewProj' in [p.name
                                      for p in ap2.get_analysis_projects_from_dirs()])
        self.assertTrue(os.path.exists(os.path.join(ap2.analysis_dir,'NewProj')))

    def test_import_project_already_in_metadata_file(self):
        """AutoProcess.import_project fails if project exists in projects.info
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_M00879_0087_000000000-AGEW9',
            'miseq',
            top_dir=self.dirn)
        mockdir.create()
        # Add the project to projects.info
        with open(os.path.join(mockdir.dirn,'projects.info'),'a') as fp:
            fp.write('%s\n' % '\t'.join(('NewProj','NP01',
                                         '.','.','.','.','.')))
        # Import the project
        ap = AutoProcess(mockdir.dirn)
        self.assertRaises(Exception,
                          ap.import_project,self.new_project_dir)


    def test_import_project_directory_already_exists(self):
        """AutoProcess.import_project fails if directory already exists
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_M00879_0087_000000000-AGEW9',
            'miseq',
            top_dir=self.dirn)
        mockdir.create()
        # Make an existing subdirectory with same name as target project
        os.mkdir(os.path.join(mockdir.dirn,'NewProj'))
        # Import the project
        ap = AutoProcess(mockdir.dirn)
        self.assertRaises(Exception,
                          ap.import_project,self.new_project_dir)
