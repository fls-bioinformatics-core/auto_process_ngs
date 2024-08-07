#######################################################################
# Tests for publish_qc_cmd.py module
#######################################################################

import unittest
import tempfile
import shutil
import os
from auto_process_ngs.settings import Settings
from auto_process_ngs.auto_processor import AutoProcess
from auto_process_ngs.mock import MockAnalysisDirFactory
from auto_process_ngs.mock import MockAnalysisProject
from auto_process_ngs.mock import UpdateAnalysisDir
from auto_process_ngs.mock import UpdateAnalysisProject
from auto_process_ngs.mock import MockMultiQC
from auto_process_ngs.commands.publish_qc_cmd import publish_qc

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

class TestAutoProcessPublishQc(unittest.TestCase):
    """
    Tests for AutoProcess.publish_qc
    """
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestAutoProcessPublishQc')
        # Create settings instance
        # This allows us to set the polling interval for the
        # unit tests
        settings_ini = os.path.join(self.dirn,"auto_process.ini")
        with open(settings_ini,'w') as s:
            s.write("""[general]
poll_interval = 0.1
""")
        self.settings = Settings(settings_ini)
        # Create a temp 'bin' dir
        self.bin = os.path.join(self.dirn,"bin")
        os.mkdir(self.bin)
        # Store original location so we can get back at the end
        self.pwd = os.getcwd()
        # Store original PATH
        self.path = os.environ['PATH']
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
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.dirn)

    def test_publish_qc_metadata_missing(self):
        """publish_qc: raise exception if metadata not set
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "160621" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        ap = AutoProcess(mockdir.dirn,
                         settings=self.settings)
        # Make a mock publication area
        publication_dir = os.path.join(self.dirn,'QC')
        os.mkdir(publication_dir)
        # Publish
        self.assertRaises(Exception,
                          publish_qc,
                          ap,
                          location=publication_dir)

    def test_publish_qc_processing_qc(self):
        """publish_qc: processing QC report only
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local",
                       "instrument_datestamp": "160621" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        ap = AutoProcess(mockdir.dirn,
                         settings=self.settings)
        # Add processing report
        UpdateAnalysisDir(ap).add_processing_report()
        # Make a mock publication area
        publication_dir = os.path.join(self.dirn,'QC')
        os.mkdir(publication_dir)
        # Publish QC
        publish_qc(ap,location=publication_dir)
        # Check outputs
        outputs = ("index.html",
                   "processing_qc.html")
        for item in outputs:
            f = os.path.join(publication_dir,
                             "160621_K00879_0087_000000000-AGEW9_analysis",
                             item)
            self.assertTrue(os.path.exists(f),"Missing %s" % f)

    def test_publish_qc_multiple_processing_qc(self):
        """publish_qc: handle multiple processing QC reports
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local",
                       "instrument_datestamp": "160621" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        ap = AutoProcess(mockdir.dirn,
                         settings=self.settings)
        # Add multiple processing reports with non-standard names
        UpdateAnalysisDir(ap).add_processing_report("processing_qc_v1.html")
        UpdateAnalysisDir(ap).add_processing_report("processing_qc_v2.html")
        # Make a mock publication area
        publication_dir = os.path.join(self.dirn,'QC')
        os.mkdir(publication_dir)
        # Publish QC
        publish_qc(ap,location=publication_dir)
        # Check outputs
        outputs = ("index.html",
                   "processing_qc_v1.html",
                   "processing_qc_v2.html",)
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
                       "source": "local",
                       "instrument_datestamp": "160621" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        ap = AutoProcess(mockdir.dirn,
                         settings=self.settings)
        # Add processing and barcode analysis reports
        UpdateAnalysisDir(ap).add_processing_report()
        UpdateAnalysisDir(ap).add_barcode_analysis()
        # Make a mock publication area
        publication_dir = os.path.join(self.dirn,'QC')
        os.mkdir(publication_dir)
        # Publish QC
        publish_qc(ap,location=publication_dir)
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

    def test_publish_qc_multiple_barcode_analyses(self):
        """publish_qc: handle multiple barcode analyses outputs
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local",
                       "instrument_datestamp": "160621" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        ap = AutoProcess(mockdir.dirn,
                         settings=self.settings)
        # Add processing and multiple barcode analysis reports
        UpdateAnalysisDir(ap).add_processing_report()
        UpdateAnalysisDir(ap).add_barcode_analysis("barcode_analysis_v1")
        UpdateAnalysisDir(ap).add_barcode_analysis("barcode_analysis_v2")
        # Make a mock publication area
        publication_dir = os.path.join(self.dirn,'QC')
        os.mkdir(publication_dir)
        # Publish QC
        publish_qc(ap,location=publication_dir)
        # Check outputs
        outputs = ("index.html",
                   "processing_qc.html",
                   os.path.join("barcode_analysis_v1","barcodes.report"),
                   os.path.join("barcode_analysis_v1","barcodes.xls"),
                   os.path.join("barcode_analysis_v1","barcodes.html"),
                   os.path.join("barcode_analysis_v2","barcodes.report"),
                   os.path.join("barcode_analysis_v2","barcodes.xls"),
                   os.path.join("barcode_analysis_v2","barcodes.html"))
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
                       "source": "local",
                       "instrument_datestamp": "160621" },
            top_dir=self.dirn)
        mockdir.create()
        ap = AutoProcess(mockdir.dirn,
                         settings=self.settings)
        # Add processing report and QC outputs
        UpdateAnalysisDir(ap).add_processing_report()
        for project in ap.get_analysis_projects():
            UpdateAnalysisProject(project).add_qc_outputs(
                include_multiqc=False)
        # Make a mock publication area
        publication_dir = os.path.join(self.dirn,'QC')
        os.mkdir(publication_dir)
        # Publish
        publish_qc(ap,location=publication_dir)
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
                                              project.info.run)
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
                       "source": "local",
                       "instrument_datestamp": "160621" },
            top_dir=self.dirn)
        mockdir.create()
        ap = AutoProcess(mockdir.dirn,
                         settings=self.settings)
        # Add processing report and QC outputs
        UpdateAnalysisDir(ap).add_processing_report()
        for project in ap.get_analysis_projects():
            UpdateAnalysisProject(project).add_qc_outputs()
        # Make a mock publication area
        publication_dir = os.path.join(self.dirn,'QC')
        os.mkdir(publication_dir)
        # Publish
        publish_qc(ap,location=publication_dir)
        # Check outputs
        outputs = ["index.html",
                   "processing_qc.html"]
        for project in ap.get_analysis_projects():
            # Standard QC outputs
            project_qc = "qc_report.%s.%s" % (project.name,
                                              project.info.run)
            outputs.append(project_qc)
            outputs.append("%s.zip" % project_qc)
            outputs.append(os.path.join(project_qc,"qc_report.html"))
            outputs.append(os.path.join(project_qc,"qc"))
        for item in outputs:
            f = os.path.join(publication_dir,
                             "160621_K00879_0087_000000000-AGEW9_analysis",
                             item)
            self.assertTrue(os.path.exists(f),"Missing %s" % f)

    def test_publish_qc_with_projects_old_zip_names(self):
        """publish_qc: projects with all QC outputs (old-style ZIP names)
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local",
                       "instrument_datestamp": "160621" },
            top_dir=self.dirn)
        mockdir.create()
        ap = AutoProcess(mockdir.dirn,
                         settings=self.settings)
        # Add processing report and QC outputs
        UpdateAnalysisDir(ap).add_processing_report()
        for project in ap.get_analysis_projects():
            UpdateAnalysisProject(project).add_qc_outputs(legacy_zip_name=True)
        # Make a mock publication area
        publication_dir = os.path.join(self.dirn,'QC')
        os.mkdir(publication_dir)
        # Publish
        publish_qc(ap,location=publication_dir)
        # Check outputs
        outputs = ["index.html",
                   "processing_qc.html"]
        for project in ap.get_analysis_projects():
            # Standard QC outputs with old-style ZIP file names
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

    def test_publish_qc_with_projects_no_reports(self):
        """publish_qc: projects with all QC outputs but no reports
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local",
                       "instrument_datestamp": "160621" },
            top_dir=self.dirn)
        mockdir.create()
        ap = AutoProcess(mockdir.dirn,
                         settings=self.settings)
        # Add processing report and QC outputs
        UpdateAnalysisDir(ap).add_processing_report()
        for project in ap.get_analysis_projects():
            UpdateAnalysisProject(project).add_qc_outputs()
        # Remove the QC reports
        for project in ap.get_analysis_projects():
            qc_reports = []
            qc_reports.append("qc_report.%s.%s.zip" %
                              (project.name,
                               project.info.run))
            qc_reports.append("qc_report.html")
            qc_reports.append("multiqc_report.html")
            for f in qc_reports:
                os.remove(os.path.join(project.dirn,f))
        # Make a mock multiqc
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make a mock publication area
        publication_dir = os.path.join(self.dirn,'QC')
        os.mkdir(publication_dir)
        # Publish
        publish_qc(ap,location=publication_dir)
        # Check outputs
        outputs = ["index.html",
                   "processing_qc.html"]
        for project in ap.get_analysis_projects():
            # Standard QC outputs
            project_qc = "qc_report.%s.%s" % (project.name,
                                              project.info.run)
            outputs.append(project_qc)
            outputs.append("%s.zip" % project_qc)
            outputs.append(os.path.join(project_qc,"qc_report.html"))
            outputs.append(os.path.join(project_qc,"qc"))
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
                       "source": "local",
                       "instrument_datestamp": "160621" },
            top_dir=self.dirn)
        mockdir.create()
        ap = AutoProcess(mockdir.dirn,
                         settings=self.settings)
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
        publish_qc(ap,location=publication_dir)
        # Check outputs
        outputs = ["index.html",
                   "processing_qc.html"]
        for project in ap.get_analysis_projects():
            # Standard QC outputs
            project_qc = "qc_report.%s.%s" % (project.name,
                                              project.info.run)
            outputs.append(project_qc)
            outputs.append("%s.zip" % project_qc)
            outputs.append(os.path.join(project_qc,"qc_report.html"))
            outputs.append(os.path.join(project_qc,"qc"))
        # Additional QC for second fastq set in first project
        project_qc = "qc.extra_report.%s.%s" % (multi_fastqs_project.name,
                                                project.info.run)
        outputs.append(project_qc)
        outputs.append("%s.zip" % project_qc)
        outputs.append(os.path.join(project_qc,"qc.extra_report.html"))
        outputs.append(os.path.join(project_qc,"qc.extra"))
        for item in outputs:
            f = os.path.join(publication_dir,
                             "160621_K00879_0087_000000000-AGEW9_analysis",
                             item)
            self.assertTrue(os.path.exists(f),"Missing %s" % f)

    def test_publish_qc_with_projects_custom_qc_protocols(self):
        """publish_qc: projects with custom QC protocols
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local",
                       "instrument_datestamp": "160621" },
            top_dir=self.dirn)
        mockdir.create()
        ap = AutoProcess(mockdir.dirn,
                         settings=self.settings)
        # Add processing report and QC outputs
        UpdateAnalysisDir(ap).add_processing_report()
        for project in ap.get_analysis_projects():
            UpdateAnalysisProject(project).add_qc_outputs()
        # Update QC metadata with custom protocol names and specifications
        for project in ap.get_analysis_projects():
            qc_info = project.qc_info("qc")
            qc_info['protocol'] = "Custom_QC"
            qc_info['protocol_specification'] = "Custom_QC:'Custom QC':seq_reads=[r1,r2]:index_reads=[]:qc_modules=[fastq_screen,fastqc,sequence_lengths]"
            qc_info.save()
        # Make a mock publication area
        publication_dir = os.path.join(self.dirn,'QC')
        os.mkdir(publication_dir)
        # Publish
        publish_qc(ap,location=publication_dir)
        # Check outputs
        outputs = ["index.html",
                   "processing_qc.html"]
        for project in ap.get_analysis_projects():
            # Standard QC outputs
            project_qc = "qc_report.%s.%s" % (project.name,
                                              project.info.run)
            outputs.append(project_qc)
            outputs.append("%s.zip" % project_qc)
            outputs.append(os.path.join(project_qc,"qc_report.html"))
            outputs.append(os.path.join(project_qc,"qc"))
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
                       "source": "local",
                       "instrument_datestamp": "160621" },
            top_dir=self.dirn)
        mockdir.create()
        ap = AutoProcess(mockdir.dirn,
                         settings=self.settings)
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
                          publish_qc,
                          ap,
                          location=publication_dir)

    def test_publish_qc_ignore_project_missing_qc(self):
        """publish_qc: skip project with missing QC
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local",
                       "instrument_datestamp": "160621" },
            top_dir=self.dirn)
        mockdir.create()
        ap = AutoProcess(mockdir.dirn,
                         settings=self.settings)
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
        publish_qc(ap,location=publication_dir,
                   ignore_missing_qc=True)
        # Check outputs
        outputs = ["index.html",
                   "processing_qc.html"]
        for project in projects:
            # Standard QC outputs
            project_qc = "qc_report.%s.%s" % (project.name,
                                              project.info.run)
            outputs.append(project_qc)
            outputs.append("%s.zip" % project_qc)
            outputs.append(os.path.join(project_qc,"qc_report.html"))
            outputs.append(os.path.join(project_qc,"qc"))
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

    def test_publish_qc_ignore_additional_project_dir(self):
        """publish_qc: ignore project-like directory not in projects.info
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local",
                       "instrument_datestamp": "160621" },
            top_dir=self.dirn)
        mockdir.create()
        ap = AutoProcess(mockdir.dirn,
                         settings=self.settings)
        # Add processing report
        UpdateAnalysisDir(ap).add_processing_report()
        # Add QC outputs for subset of projects
        projects = ap.get_analysis_projects()
        for project in ap.get_analysis_projects():
            UpdateAnalysisProject(project).add_qc_outputs()
        # Make an additional project-like directory
        project_like = MockAnalysisProject('additional_dir',
                                           ('random1_R1.fastq',
                                            'random2_R1.fastq',))
        project_like.create(top_dir=mockdir.dirn)
        # Make a mock publication area
        publication_dir = os.path.join(self.dirn,'QC')
        os.mkdir(publication_dir)
        # Publish
        publish_qc(ap,location=publication_dir)
        # Check outputs
        outputs = ["index.html",
                   "processing_qc.html"]
        for project in projects:
            # Standard QC outputs
            project_qc = "qc_report.%s.%s" % (project.name,
                                              project.info.run)
            outputs.append(project_qc)
            outputs.append("%s.zip" % project_qc)
            outputs.append(os.path.join(project_qc,"qc_report.html"))
            outputs.append(os.path.join(project_qc,"qc"))
        for item in outputs:
            f = os.path.join(publication_dir,
                             "160621_K00879_0087_000000000-AGEW9_analysis",
                             item)
            self.assertTrue(os.path.exists(f),"Missing %s" % f)

    def test_publish_qc_subset_of_projects(self):
        """publish_qc: only publish subset of projects
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local",
                       "instrument_datestamp": "160621" },
            top_dir=self.dirn)
        mockdir.create()
        ap = AutoProcess(mockdir.dirn,
                         settings=self.settings)
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
        publish_qc(ap,location=publication_dir,
                   projects="AB*")
        # Check outputs
        outputs = ["index.html",
                   "processing_qc.html"]
        for project in projects:
            # Standard QC outputs
            project_qc = "qc_report.%s.%s" % (project.name,
                                              project.info.run)
            outputs.append(project_qc)
            outputs.append("%s.zip" % project_qc)
            outputs.append(os.path.join(project_qc,"qc_report.html"))
            outputs.append(os.path.join(project_qc,"qc"))
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
                                                  project.info.run))),
                             "%s exists in final dir, but shouldn't" %
                             project.name)

    def test_publish_qc_with_icell8_outputs(self):
        """publish_qc: project with ICELL8 QC outputs
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local",
                       "instrument_datestamp": "160621" },
            top_dir=self.dirn)
        mockdir.create()
        ap = AutoProcess(mockdir.dirn,
                         settings=self.settings)
        # Add processing report and QC outputs
        UpdateAnalysisDir(ap).add_processing_report()
        projects = ap.get_analysis_projects()
        for project in projects:
            UpdateAnalysisProject(project).add_qc_outputs(
                protocol="singlecell")
        # Add ICELL8 report for one project
        icell8_project = projects[0]
        UpdateAnalysisProject(icell8_project).add_icell8_outputs()
        # Make a mock publication area
        publication_dir = os.path.join(self.dirn,'QC')
        os.mkdir(publication_dir)
        # Publish
        publish_qc(ap,location=publication_dir)
        # Check outputs
        outputs = ["index.html",
                   "processing_qc.html"]
        for project in ap.get_analysis_projects():
            # Standard QC outputs
            project_qc = "qc_report.%s.%s" % (project.name,
                                              project.info.run)
            outputs.append(project_qc)
            outputs.append("%s.zip" % project_qc)
            outputs.append(os.path.join(project_qc,"qc_report.html"))
            outputs.append(os.path.join(project_qc,"qc"))
        # Do checks
        for item in outputs:
            f = os.path.join(publication_dir,
                             "160621_K00879_0087_000000000-AGEW9_analysis",
                             item)
            self.assertTrue(os.path.exists(f),"Missing %s" % f)
        # ICELL8 outputs shouldn't be present
        icell8_dir = "icell8_processing.%s.%s" % (icell8_project.name,
                                                  os.path.basename(
                                                      ap.analysis_dir))
        outputs = []
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
            self.assertFalse(os.path.exists(f),"Found %s" % f)

    def test_publish_qc_with_icell8_outputs_legacy_mode(self):
        """publish_qc: project with ICELL8 QC outputs (legacy mode)
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local",
                       "instrument_datestamp": "160621" },
            top_dir=self.dirn)
        mockdir.create()
        ap = AutoProcess(mockdir.dirn,
                         settings=self.settings)
        # Add processing report and QC outputs
        UpdateAnalysisDir(ap).add_processing_report()
        projects = ap.get_analysis_projects()
        for project in projects:
            UpdateAnalysisProject(project).add_qc_outputs(
                protocol="singlecell")
        # Add ICELL8 report for one project
        icell8_project = projects[0]
        UpdateAnalysisProject(icell8_project).add_icell8_outputs()
        # Make a mock publication area
        publication_dir = os.path.join(self.dirn,'QC')
        os.mkdir(publication_dir)
        # Publish
        publish_qc(ap,location=publication_dir,legacy=True)
        # Check outputs
        outputs = ["index.html",
                   "processing_qc.html"]
        for project in ap.get_analysis_projects():
            # Standard QC outputs
            project_qc = "qc_report.%s.%s" % (project.name,
                                              project.info.run)
            outputs.append(project_qc)
            outputs.append("%s.zip" % project_qc)
            outputs.append(os.path.join(project_qc,"qc_report.html"))
            outputs.append(os.path.join(project_qc,"qc"))
        # ICELL8 outputs
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
                       "source": "local",
                       "instrument_datestamp": "160621" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        ap = AutoProcess(mockdir.dirn,
                         settings=self.settings)
        # Add processing and cellranger QC reports
        UpdateAnalysisDir(ap).add_processing_report()
        UpdateAnalysisDir(ap).add_cellranger_qc_output()
        # Make a mock publication area
        publication_dir = os.path.join(self.dirn,'QC')
        os.mkdir(publication_dir)
        # Publish
        publish_qc(ap,location=publication_dir)
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
                       "source": "local",
                       "instrument_datestamp": "160621" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        ap = AutoProcess(mockdir.dirn,
                         settings=self.settings)
        # Add processing and cellranger QC reports
        UpdateAnalysisDir(ap).add_processing_report()
        UpdateAnalysisDir(ap).add_cellranger_qc_output(lanes="45")
        # Make a mock publication area
        publication_dir = os.path.join(self.dirn,'QC')
        os.mkdir(publication_dir)
        # Publish
        publish_qc(ap,location=publication_dir)
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
                       "source": "local",
                       "instrument_datestamp": "160621" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        ap = AutoProcess(mockdir.dirn,
                         settings=self.settings)
        # Add processing and cellranger QC reports
        UpdateAnalysisDir(ap).add_processing_report()
        UpdateAnalysisDir(ap).add_cellranger_qc_output(lanes="45")
        UpdateAnalysisDir(ap).add_cellranger_qc_output(lanes="78")
        # Make a mock publication area
        publication_dir = os.path.join(self.dirn,'QC')
        os.mkdir(publication_dir)
        # Publish
        publish_qc(ap,location=publication_dir)
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

    def test_publish_qc_with_multiple_10x_mkfastq_qc(self):
        """publish_qc: publish 10xGenomics mkfastq QC output (multiple packages)
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local",
                       "instrument_datestamp": "160621" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        ap = AutoProcess(mockdir.dirn,
                         settings=self.settings)
        # Add processing and cellranger QC reports
        UpdateAnalysisDir(ap).add_processing_report()
        UpdateAnalysisDir(ap).add_10x_mkfastq_qc_output("cellranger-arc",
                                                        lanes="12")
        UpdateAnalysisDir(ap).add_10x_mkfastq_qc_output("spaceranger",
                                                        lanes="3")
        UpdateAnalysisDir(ap).add_10x_mkfastq_qc_output("cellranger-atac",
                                                        lanes="45")
        UpdateAnalysisDir(ap).add_10x_mkfastq_qc_output("cellranger",
                                                        lanes="78")
        # Make a mock publication area
        publication_dir = os.path.join(self.dirn,'QC')
        os.mkdir(publication_dir)
        # Publish
        publish_qc(ap,location=publication_dir)
        # Check outputs
        outputs = ["index.html",
                   "processing_qc.html",
                   "cellranger-arc_qc_summary_12.html",
                   "spaceranger_qc_summary_3.html",
                   "cellranger-atac_qc_summary_45.html",
                   "cellranger_qc_summary_78.html"]
        # Do checks
        for item in outputs:
            f = os.path.join(publication_dir,
                             "160621_K00879_0087_000000000-AGEW9_analysis",
                             item)
            self.assertTrue(os.path.exists(f),"Missing %s" % f)

    def test_publish_qc_with_cellranger_count(self):
        """publish_qc: project with cellranger count output
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local",
                       "instrument_datestamp": "160621" },
            top_dir=self.dirn)
        mockdir.create()
        ap = AutoProcess(mockdir.dirn,
                         settings=self.settings)
        # Add processing and cellranger QC reports
        UpdateAnalysisDir(ap).add_processing_report()
        UpdateAnalysisDir(ap).add_cellranger_qc_output()
        # Add QC outputs
        projects = ap.get_analysis_projects()
        for project in projects:
            UpdateAnalysisProject(project).add_qc_outputs(
                protocol="singlecell")
        # Add cellranger count output for one project
        tenxgenomics_project = projects[0]
        UpdateAnalysisProject(tenxgenomics_project).add_cellranger_count_outputs()
        # Make a mock publication area
        publication_dir = os.path.join(self.dirn,'QC')
        os.mkdir(publication_dir)
        # Publish
        publish_qc(ap,location=publication_dir)
        # Check outputs
        outputs = ["index.html",
                   "processing_qc.html",
                   "cellranger_qc_summary.html"]
        for project in ap.get_analysis_projects():
            # Standard QC outputs
            project_qc = "qc_report.%s.%s" % (project.name,
                                              project.info.run)
            outputs.append(project_qc)
            outputs.append("%s.zip" % project_qc)
            outputs.append(os.path.join(project_qc,"qc_report.html"))
            outputs.append(os.path.join(project_qc,"qc"))
        # NB cellranger count outputs shouldn't be present
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

    def test_publish_qc_with_cellranger_count_legacy_mode(self):
        """publish_qc: project with cellranger count output (legacy mode)
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local",
                       "instrument_datestamp": "160621" },
            top_dir=self.dirn)
        mockdir.create()
        ap = AutoProcess(mockdir.dirn,
                         settings=self.settings)
        # Add processing and cellranger QC reports
        UpdateAnalysisDir(ap).add_processing_report()
        UpdateAnalysisDir(ap).add_cellranger_qc_output()
        # Add QC outputs
        projects = ap.get_analysis_projects()
        for project in projects:
            UpdateAnalysisProject(project).add_qc_outputs(
                protocol="singlecell")
        # Add cellranger count output for one project
        tenxgenomics_project = projects[0]
        UpdateAnalysisProject(tenxgenomics_project).add_cellranger_count_outputs(legacy=True)
        # Make a mock publication area
        publication_dir = os.path.join(self.dirn,'QC')
        os.mkdir(publication_dir)
        # Publish
        publish_qc(ap,location=publication_dir,legacy=True)
        # Check outputs
        outputs = ["index.html",
                   "processing_qc.html",
                   "cellranger_qc_summary.html"]
        for project in ap.get_analysis_projects():
            # Standard QC outputs
            project_qc = "qc_report.%s.%s" % (project.name,
                                              project.info.run)
            outputs.append(project_qc)
            outputs.append("%s.zip" % project_qc)
            outputs.append(os.path.join(project_qc,"qc_report.html"))
            outputs.append(os.path.join(project_qc,"qc"))
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

    def test_publish_qc_use_hierarchy(self):
        """publish_qc: publish using YEAR/PLATFORM hierarchy
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local",
                       "instrument_datestamp": "160621" },
            top_dir=self.dirn)
        mockdir.create()
        ap = AutoProcess(mockdir.dirn,
                         settings=self.settings)
        # Add processing report and QC outputs
        UpdateAnalysisDir(ap).add_processing_report()
        for project in ap.get_analysis_projects():
            UpdateAnalysisProject(project).add_qc_outputs()
        # Make a mock publication area
        publication_dir = os.path.join(self.dirn,'QC')
        os.mkdir(publication_dir)
        # Publish
        publish_qc(ap,location=publication_dir,
                   use_hierarchy=True)
        # Check outputs
        final_dir = os.path.join(publication_dir,
                                 "2016",
                                 "hiseq")
        self.assertTrue(os.path.exists(final_dir))
        outputs = ["index.html",
                   "processing_qc.html"]
        for project in ap.get_analysis_projects():
            # Standard QC outputs
            project_qc = "qc_report.%s.%s" % (project.name,
                                              project.info.run)
            outputs.append(project_qc)
            outputs.append("%s.zip" % project_qc)
            outputs.append(os.path.join(project_qc,"qc_report.html"))
            outputs.append(os.path.join(project_qc,"qc"))
        for item in outputs:
            f = os.path.join(final_dir,
                             "160621_K00879_0087_000000000-AGEW9_analysis",
                             item)
            self.assertTrue(os.path.exists(f),"Missing %s" % f)

    def test_publish_qc_exclude_zip_files(self):
        """publish_qc: exclude ZIP files from publication
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local",
                       "instrument_datestamp": "160621" },
            top_dir=self.dirn)
        mockdir.create()
        ap = AutoProcess(mockdir.dirn,
                         settings=self.settings)
        # Add processing report and QC outputs
        UpdateAnalysisDir(ap).add_processing_report()
        projects = ap.get_analysis_projects()
        for project in projects:
            UpdateAnalysisProject(project).add_qc_outputs()
        # Add ICELL8 report for one project
        icell8_project = projects[0]
        UpdateAnalysisProject(icell8_project).add_icell8_outputs()
        # Add cellranger count output for one project
        tenxgenomics_project = projects[-1]
        UpdateAnalysisProject(tenxgenomics_project).add_cellranger_count_outputs()
        # Make a mock publication area
        publication_dir = os.path.join(self.dirn,'QC')
        os.mkdir(publication_dir)
        # Publish
        publish_qc(ap,location=publication_dir,
                   exclude_zip_files=True)
        # Check outputs
        outputs = ["index.html",
                   "processing_qc.html"]
        zip_files = []
        for project in ap.get_analysis_projects():
            # Standard QC outputs
            project_qc = "qc_report.%s.%s" % (project.name,
                                              project.info.run)
            outputs.append(project_qc)
            outputs.append(os.path.join(project_qc,"qc_report.html"))
            outputs.append(os.path.join(project_qc,"qc"))
            zip_files.append("%s.zip" % project_qc)
        # NB ICELL8 and cellranger count outputs shouldn't be present
        # Do checks
        for item in outputs:
            f = os.path.join(publication_dir,
                             "160621_K00879_0087_000000000-AGEW9_analysis",
                             item)
            self.assertTrue(os.path.exists(f),"Missing %s" % f)
        for item in outputs:
            f = os.path.join(publication_dir,
                             "160621_K00879_0087_000000000-AGEW9_analysis",
                             item)
            self.assertTrue(os.path.exists(f),"Missing %s" % f)
        # Check the ZIP files were excluded
        for zip_file in zip_files:
            self.assertFalse(os.path.exists(
                os.path.join(publication_dir,
                             "160621_K00879_0087_000000000-AGEW9_analysis",
                             zip_file)),
                             "ZIP file '%s' exists, but shouldn't" %
                             zip_file)

    def test_publish_qc_exclude_zip_files_legacy_mode(self):
        """publish_qc: exclude ZIP files from publication (legacy mode)
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local",
                       "instrument_datestamp": "160621" },
            top_dir=self.dirn)
        mockdir.create()
        ap = AutoProcess(mockdir.dirn,
                         settings=self.settings)
        # Add processing report and QC outputs
        UpdateAnalysisDir(ap).add_processing_report()
        projects = ap.get_analysis_projects()
        for project in projects:
            UpdateAnalysisProject(project).add_qc_outputs()
        # Add ICELL8 report for one project
        icell8_project = projects[0]
        UpdateAnalysisProject(icell8_project).add_icell8_outputs()
        # Add cellranger count output for one project
        tenxgenomics_project = projects[-1]
        UpdateAnalysisProject(tenxgenomics_project).add_cellranger_count_outputs(legacy=True)
        # Make a mock publication area
        publication_dir = os.path.join(self.dirn,'QC')
        os.mkdir(publication_dir)
        # Publish
        publish_qc(ap,location=publication_dir,
                   exclude_zip_files=True,
                   legacy=True)
        # Check outputs
        outputs = ["index.html",
                   "processing_qc.html"]
        zip_files = []
        for project in ap.get_analysis_projects():
            # Standard QC outputs
            project_qc = "qc_report.%s.%s" % (project.name,
                                              project.info.run)
            outputs.append(project_qc)
            outputs.append(os.path.join(project_qc,"qc_report.html"))
            outputs.append(os.path.join(project_qc,"qc"))
            zip_files.append("%s.zip" % project_qc)
        # ICELL8 outputs
        icell8_dir = "icell8_processing.%s.%s" % (icell8_project.name,
                                                  os.path.basename(
                                                      ap.analysis_dir))
        outputs.append(icell8_dir)
        outputs.append(os.path.join(icell8_dir,"icell8_processing_data"))
        outputs.append(os.path.join(icell8_dir,"icell8_processing.html"))
        outputs.append(os.path.join(icell8_dir,"stats"))
        zip_files.append("%s.zip" % icell8_dir)
        # Cellranger count outputs
        cellranger_count_dir = "cellranger_count_report.%s.%s" % (
            tenxgenomics_project.name,
            os.path.basename(ap.analysis_dir))
        outputs.append(cellranger_count_dir)
        outputs.append(os.path.join(cellranger_count_dir,
                                    "cellranger_count_report.html"))
        outputs.append(os.path.join(cellranger_count_dir,"cellranger_count"))
        for sample in tenxgenomics_project.samples:
            outputs.append(os.path.join(cellranger_count_dir,
                                        "cellranger_count",
                                        sample.name,
                                        "outs",
                                        "web_summary.html"))
        zip_files.append("%s.zip" % cellranger_count_dir)
        # Do checks
        for item in outputs:
            f = os.path.join(publication_dir,
                             "160621_K00879_0087_000000000-AGEW9_analysis",
                             item)
            self.assertTrue(os.path.exists(f),"Missing %s" % f)
        # Check the ZIP files were excluded
        for zip_file in zip_files:
            self.assertFalse(os.path.exists(
                os.path.join(publication_dir,
                             "160621_K00879_0087_000000000-AGEW9_analysis",
                             zip_file)),
                             "ZIP file '%s' exists, but shouldn't" %
                             zip_file)

    def test_publish_qc_missing_destination(self):
        """publish_qc: raise exception if destination doesn't exist
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local",
                       "instrument_datestamp": "160621" },
            top_dir=self.dirn)
        mockdir.create()
        ap = AutoProcess(mockdir.dirn,
                         settings=self.settings)
        # Add processing report and QC outputs
        UpdateAnalysisDir(ap).add_processing_report()
        for project in ap.get_analysis_projects():
            UpdateAnalysisProject(project).add_qc_outputs()
        # Reference publication area which doesn't exist
        publication_dir = os.path.join(self.dirn,'QC')
        self.assertFalse(os.path.exists(publication_dir))
        # Publish
        self.assertRaises(Exception,
                          publish_qc,
                          ap,
                          location=publication_dir)
        self.assertFalse(os.path.exists(publication_dir))

    def test_publish_qc_handles_four_digit_year_in_datestamp(self):
        """publish_qc: handle 4-digit year in datestamp using YEAR/PLATFORM hierarchy
        """
        # Make an auto-process directory with 4-digit year datestamp
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '20160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local",
                       "instrument_datestamp": "20160621" },
            top_dir=self.dirn)
        mockdir.create()
        ap = AutoProcess(mockdir.dirn,
                         settings=self.settings)
        # Add processing report and QC outputs
        UpdateAnalysisDir(ap).add_processing_report()
        for project in ap.get_analysis_projects():
            UpdateAnalysisProject(project).add_qc_outputs()
        # Make a mock publication area
        publication_dir = os.path.join(self.dirn,'QC')
        os.mkdir(publication_dir)
        # Publish
        publish_qc(ap,location=publication_dir,
                   use_hierarchy=True)
        # Check outputs
        final_dir = os.path.join(publication_dir,
                                 "2016",
                                 "hiseq")
        self.assertTrue(os.path.exists(final_dir))
        outputs = ["index.html",
                   "processing_qc.html"]
        for project in ap.get_analysis_projects():
            # Standard QC outputs
            project_qc = "qc_report.%s.%s" % (project.name,
                                              project.info.run)
            outputs.append(project_qc)
            outputs.append("%s.zip" % project_qc)
            outputs.append(os.path.join(project_qc,"qc_report.html"))
            outputs.append(os.path.join(project_qc,"qc"))
        for item in outputs:
            f = os.path.join(final_dir,
                             "20160621_K00879_0087_000000000-AGEW9_analysis",
                             item)
            self.assertTrue(os.path.exists(f),"Missing %s" % f)

    def test_publish_qc_with_projects_legacy_mode(self):
        """publish_qc: projects with all QC outputs (legacy mode)
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local",
                       "instrument_datestamp": "160621" },
            top_dir=self.dirn)
        mockdir.create()
        ap = AutoProcess(mockdir.dirn,
                         settings=self.settings)
        # Add processing report and QC outputs
        UpdateAnalysisDir(ap).add_processing_report()
        for project in ap.get_analysis_projects():
            UpdateAnalysisProject(project).add_qc_outputs()
        # Make a mock publication area
        publication_dir = os.path.join(self.dirn,'QC')
        os.mkdir(publication_dir)
        # Publish
        publish_qc(ap,location=publication_dir,legacy=True)
        # Check outputs
        outputs = ["index.html",
                   "processing_qc.html"]
        for project in ap.get_analysis_projects():
            # Standard QC outputs
            project_qc = "qc_report.%s.%s" % (project.name,
                                              project.info.run)
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

    def test_publish_qc_erases_previously_published_data(self):
        """publish_qc: erases previously published data
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local",
                       "instrument_datestamp": "160621" },
            top_dir=self.dirn)
        mockdir.create()
        ap = AutoProcess(mockdir.dirn,
                         settings=self.settings)
        # Add processing report and QC outputs
        UpdateAnalysisDir(ap).add_processing_report()
        for project in ap.get_analysis_projects():
            UpdateAnalysisProject(project).add_qc_outputs()
        # Make a mock publication area
        publication_dir = os.path.join(self.dirn,'QC')
        os.mkdir(publication_dir)
        # Publish
        publish_qc(ap,location=publication_dir)
        # Add additional content
        extra_file = os.path.join(publication_dir,
                                  "160621_K00879_0087_000000000-AGEW9_analysis",
                                  "extra_file.txt")
        with open(extra_file,'wt') as fp:
            fp.write("Placeholder")
        self.assertTrue(os.path.exists(extra_file))
        # Republish
        publish_qc(ap,location=publication_dir)
        # Check additional content was removed
        self.assertFalse(os.path.exists(extra_file))

    def test_publish_qc_doesnt_erase_non_qc_data(self):
        """publish_qc: doesn't erase non-QC report data
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local",
                       "instrument_datestamp": "160621" },
            top_dir=self.dirn)
        mockdir.create()
        ap = AutoProcess(mockdir.dirn,
                         settings=self.settings)
        # Add processing report and QC outputs
        UpdateAnalysisDir(ap).add_processing_report()
        for project in ap.get_analysis_projects():
            UpdateAnalysisProject(project).add_qc_outputs()
        # Try to publish to same location as analysis dir
        self.assertRaises(Exception,
                          publish_qc,
                          ap,
                          location=self.dirn)
        # Check analysis dir is unmolested
        self.assertTrue(os.path.exists(os.path.join(ap.analysis_dir,
                                                    "metadata.info")))
        # Make a directory with non-QC content
        target_dir = os.path.join(self.dirn,
                                  "data",
                                  os.path.basename(ap.analysis_dir))
        os.makedirs(target_dir)
        with open(os.path.join(target_dir,"stuff"),'wt') as fp:
            fp.write("Some arbitrary content")
        # Try to publish to location with non-QC content
        self.assertRaises(Exception,
                          publish_qc,
                          ap,
                          location=os.path.dirname(target_dir))
        # Check non-QC content is unmolested
        self.assertTrue(os.path.exists(os.path.join(target_dir,"stuff")))

    def test_publish_qc_with_backup_project(self):
        """publish_qc: handle 'backup' copy of project
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local",
                       "instrument_datestamp": "160621" },
            top_dir=self.dirn)
        mockdir.create()
        ap = AutoProcess(mockdir.dirn,
                         settings=self.settings)
        # Add processing report and QC outputs
        UpdateAnalysisDir(ap).add_processing_report()
        for project in ap.get_analysis_projects():
            UpdateAnalysisProject(project).add_qc_outputs()
        # Add extra 'backup' project by copying the
        # existing undetermine project
        shutil.copytree(
            os.path.join(self.dirn,
                         "160621_K00879_0087_000000000-AGEW9_analysis",
                         "undetermined"),
            os.path.join(self.dirn,
                         "undetermined.bak"))
        shutil.move(
            os.path.join(self.dirn,
                         "undetermined.bak"),
            os.path.join(self.dirn,
                         "160621_K00879_0087_000000000-AGEW9_analysis"))
        self.assertTrue(os.path.isdir(
            os.path.join(self.dirn,
                         "160621_K00879_0087_000000000-AGEW9_analysis",
                         "undetermined.bak")))
        # Make a mock publication area
        publication_dir = os.path.join(self.dirn,'QC')
        os.mkdir(publication_dir)
        # Publish
        publish_qc(ap,location=publication_dir)
        # Check outputs
        outputs = ["index.html",
                   "processing_qc.html"]
        for project in ap.get_analysis_projects():
            # Standard QC outputs
            project_qc = "qc_report.%s.%s" % (project.name,
                                              project.info.run)
            outputs.append(project_qc)
            outputs.append("%s.zip" % project_qc)
            outputs.append(os.path.join(project_qc,"qc_report.html"))
            outputs.append(os.path.join(project_qc,"qc"))
        for item in outputs:
            f = os.path.join(publication_dir,
                             "160621_K00879_0087_000000000-AGEW9_analysis",
                             item)
            self.assertTrue(os.path.exists(f),"Missing %s" % f)
        self.assertFalse(os.path.exists(
            os.path.join(publication_dir,
                         "160621_K00879_0087_000000000-AGEW9_analysis",
                         "undetermined.bak")))
