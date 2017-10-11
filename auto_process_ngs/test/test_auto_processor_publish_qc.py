#######################################################################
# Tests for 'auto_process publish_qc' command
#######################################################################

import unittest
import tempfile
import shutil
import os
import logging
from auto_process_ngs.auto_processor import AutoProcess
from auto_process_ngs.utils import AnalysisProject
from auto_process_ngs.utils import ZipArchive
from auto_process_ngs.qc.illumina_qc import expected_qc_outputs
from auto_process_ngs.mock import MockAnalysisDirFactory

# Module specific logger
logger = logging.getLogger(__name__)

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

    def _make_file(self,filen,content=None):
        # Make a file
        with open(filen,'w') as fp:
            if content is not None:
                fp.write("%s" % content)
            else:
                fp.write("")

    def _add_processing_report(self,dirn):
        # Add mock processing report
        self._make_file(os.path.join(dirn,"processing_qc.html"))

    def _add_qc_outputs(self,project_dir,include_multiqc=True):
        # Add mock QC outputs
        logger.debug("_add_qc_outputs: adding QC outputs for %s" %
                     project_dir)
        p = AnalysisProject(os.path.basename(project_dir),
                            project_dir)
        p.setup_qc_dir()
        logger.debug("_add_qc_outputs: QC dir is %s" % p.qc_dir)
        for fq in p.fastqs:
            logger.debug("_add_qc_outputs: adding outputs for %s" % fq)
            for f in expected_qc_outputs(fq,p.qc_dir):
                logger.debug("_add_qc_outputs: %s" % f)
                self._make_file(f)
        logger.debug("_add_qc_outputs: adding qc_report.html")
        self._make_file(os.path.join(p.dirn,"qc_report.html"))
        logger.debug("_add_qc_outputs: adding ZIP file")
        analysis_dir = os.path.basename(os.path.dirname(p.dirn))
        logger.debug("_add_qc_outputs: analysis dir %s" % analysis_dir)
        report_zip = os.path.join(p.dirn,
                                  "qc_report.%s.%s.zip" %
                                  (p.name,
                                   analysis_dir))
        logger.debug("_add_qc_outputs: ZIP file name %s" % report_zip)
        zip_file = ZipArchive(report_zip,
                              relpath=p.dirn,
                              prefix="qc_report.%s.%s" %
                              (p.name,
                               analysis_dir))
        zip_file.add_file(os.path.join(p.dirn,"qc_report.html"))
        zip_file.add(p.qc_dir)
        zip_file.close()
        if include_multiqc:
            logger.debug("_add_qc_outputs: adding MultiQC output")
            self._make_file(os.path.join(p.dirn,"multiqc_report.html"))

    def _add_icell8_outputs(self,project_dir):
        # Add mock ICell8 outputs
        logger.debug("_add_icell8_outputs: adding ICell8 outputs for %s" %
                     project_dir)
        p = AnalysisProject(os.path.basename(project_dir),
                            project_dir)
        self._make_file(os.path.join(p.dirn,"icell8_processing.html"))
        os.mkdir(os.path.join(p.dirn,"stats"))
        self._make_file(os.path.join(p.dirn,
                                     "stats",
                                     "icell8_stats.tsv"))
        self._make_file(os.path.join(p.dirn,
                                     "stats",
                                     "icell8_stats.xlsx"))
        os.mkdir(os.path.join(p.dirn,"icell8_processing_data"))
        for png in ("poly_g_dist.png",
                    "read_dist.png",
                    "reads_per_stage.png",
                    "samples.png"):
            self._make_file(os.path.join(p.dirn,
                                         "icell8_processing_data",
                                         png))
        # Build ZIP archive
        analysis_dir = os.path.basename(os.path.dirname(p.dirn))
        icell8_zip = os.path.join(p.dirn,
                                  "icell8_processing.%s.%s.zip" %
                                  (p.name,
                                   analysis_dir))
        zip_file = ZipArchive(icell8_zip,
                              relpath=p.dirn,
                              prefix="icell8_processing.%s.%s" %
                              (p.name,
                               analysis_dir))
        zip_file.add_file(os.path.join(p.dirn,"icell8_processing.html"))
        zip_file.add(os.path.join(p.dirn,"stats"))
        zip_file.add(os.path.join(p.dirn,"icell8_processing_data"))
        zip_file.close()

    def _add_cellranger_qc_output(self,dirn,lanes=None):
        # Add mock cellranger qc report
        cellranger_qc = "cellranger_qc_summary%s.html" % ("_%s" % lanes
                                                          if lanes is not None
                                                          else "")
        self._make_file(os.path.join(dirn,cellranger_qc))

    def _add_cellranger_count_outputs(self,project_dir):
        # Add mock cellranger count outputs
        logger.debug("_add_cellranger_count_outputs: adding cellranger "
                     "count outputs for %s" % project_dir)
        p = AnalysisProject(os.path.basename(project_dir),
                            project_dir)
        self._make_file(os.path.join(p.dirn,"cellranger_count_report.html"))
        os.mkdir(os.path.join(p.dirn,"cellranger_count"))
        for sample in p.samples:
            sample_dir = os.path.join(p.dirn,"cellranger_count",sample.name)
            os.mkdir(sample_dir)
            os.mkdir(os.path.join(sample_dir,"outs"))
            self._make_file(os.path.join(sample_dir,
                                         "outs",
                                         "web_summary.html"))
        # Build ZIP archive
        analysis_dir = os.path.basename(os.path.dirname(p.dirn))
        cellranger_zip = os.path.join(p.dirn,
                                      "cellranger_count_report.%s.%s.zip" %
                                      (p.name,
                                       analysis_dir))
        zip_file = ZipArchive(cellranger_zip,
                              relpath=p.dirn,
                              prefix="cellranger_count_report.%s.%s" %
                              (p.name,
                               analysis_dir))
        zip_file.add_file(os.path.join(p.dirn,"cellranger_count_report.html"))
        zip_file.add(os.path.join(p.dirn,"cellranger_count"))
        zip_file.close()

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
        self._add_processing_report(ap.analysis_dir)
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
        # Add processing report
        self._add_processing_report(ap.analysis_dir)
        # Add QC outputs
        for project in ap.get_analysis_projects():
            self._add_qc_outputs(project.dirn,
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
        # Add processing report
        self._add_processing_report(ap.analysis_dir)
        # Add QC outputs
        for project in ap.get_analysis_projects():
            self._add_qc_outputs(project.dirn)
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
        self._add_processing_report(ap.analysis_dir)
        # Add QC outputs for subset of projects
        projects = ap.get_analysis_projects()[1:]
        for project in projects:
            self._add_qc_outputs(project.dirn)
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
        self._add_processing_report(ap.analysis_dir)
        # Add QC outputs for subset of projects
        projects = ap.get_analysis_projects()
        missing_project = projects[0]
        projects = projects[1:]
        for project in projects:
            self._add_qc_outputs(project.dirn)
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
        self._add_processing_report(ap.analysis_dir)
        # Add QC outputs for subset of projects
        projects = ap.get_analysis_projects()
        missing_projects = projects[1:]
        projects = projects[0:1]
        for project in ap.get_analysis_projects():
            self._add_qc_outputs(project.dirn)
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
        # Add processing report
        self._add_processing_report(ap.analysis_dir)
        # Add QC outputs
        projects = ap.get_analysis_projects()
        for project in projects:
            self._add_qc_outputs(project.dirn)
        # Add ICell8 report for one project
        icell8_project = projects[0]
        self._add_icell8_outputs(icell8_project.dirn)
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
        # Add processing report
        self._add_processing_report(ap.analysis_dir)
        # Add cellranger QC report
        self._add_cellranger_qc_output(ap.analysis_dir)
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
        # Add processing report
        self._add_processing_report(ap.analysis_dir)
        # Add cellranger QC report
        self._add_cellranger_qc_output(ap.analysis_dir,lanes="45")
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
        # Add processing report
        self._add_processing_report(ap.analysis_dir)
        # Add cellranger QC reports
        self._add_cellranger_qc_output(ap.analysis_dir,lanes="45")
        self._add_cellranger_qc_output(ap.analysis_dir,lanes="78")
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
        # Add processing report
        self._add_processing_report(ap.analysis_dir)
        # Add cellranger QC report
        self._add_cellranger_qc_output(ap.analysis_dir)
        # Add QC outputs
        projects = ap.get_analysis_projects()
        for project in projects:
            self._add_qc_outputs(project.dirn)
        # Add cellranger count output for one project
        tenxgenomics_project = projects[0]
        self._add_cellranger_count_outputs(tenxgenomics_project.dirn)
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
