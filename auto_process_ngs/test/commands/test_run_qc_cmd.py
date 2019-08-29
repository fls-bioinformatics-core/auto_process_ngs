#######################################################################
# Tests for run_qc_cmd.py module
#######################################################################

import unittest
import tempfile
import shutil
import os
import zipfile
from bcftbx.JobRunner import SimpleJobRunner
from auto_process_ngs.auto_processor import AutoProcess
from auto_process_ngs.mock import MockAnalysisDirFactory
from auto_process_ngs.mock import MockIlluminaQcSh
from auto_process_ngs.mock import MockFastqStrandPy
from auto_process_ngs.mock import MockCellrangerExe
from auto_process_ngs.mock import MockMultiQC
from auto_process_ngs.settings import Settings
from auto_process_ngs.commands.run_qc_cmd import run_qc

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

class TestAutoProcessRunQc(unittest.TestCase):
    """
    Tests for AutoProcess.run_qc
    """
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestAutoProcessRunQc')
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
        # Return to original dir
        os.chdir(self.pwd)
        # Restore PATH
        os.environ['PATH'] = self.path
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.dirn)

    def test_run_qc(self):
        """run_qc: standard QC run
        """
        # Make mock illumina_qc.sh and multiqc
        MockIlluminaQcSh.create(os.path.join(self.bin,
                                             "illumina_qc.sh"))
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            top_dir=self.dirn)
        mockdir.create()
        # Settings file with polling interval
        settings_ini = os.path.join(self.dirn,"settings.ini")
        with open(settings_ini,'w') as s:
            s.write("""[general]
poll_interval = 0.5
""")
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn,
                         settings=Settings(settings_ini))
        # Run the QC
        status = run_qc(ap,
                        run_multiqc=True,
                        max_jobs=1)
        self.assertEqual(status,0)
        # Check output and reports
        for p in ("AB","CDE","undetermined"):
            for f in ("qc",
                      "qc_report.html",
                      "qc_report.%s.%s_analysis.zip" % (
                          p,
                          '170901_M00879_0087_000000000-AGEW9'),
                      "multiqc_report.html"):
                self.assertTrue(os.path.exists(os.path.join(mockdir.dirn,
                                                            p,f)),
                                "Missing %s in project '%s'" % (f,p))
            # Check zip file has MultiQC report
            zip_file = os.path.join(mockdir.dirn,p,
                                    "qc_report.%s.%s_analysis.zip" % (
                                        p,
                                        '170901_M00879_0087_000000000-AGEW9'))
            with zipfile.ZipFile(zip_file) as z:
                multiqc = os.path.join(
                    "qc_report.%s.%s_analysis" % (
                        p,'170901_M00879_0087_000000000-AGEW9'),
                    "multiqc_report.html")
                self.assertTrue(multiqc in z.namelist())

    def test_run_qc_with_strandedness(self):
        """run_qc: standard QC run with strandedness determination
        """
        # Make mock illumina_qc.sh and multiqc
        MockIlluminaQcSh.create(os.path.join(self.bin,
                                             "illumina_qc.sh"))
        MockFastqStrandPy.create(os.path.join(self.bin,
                                              "fastq_strand.py"))
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            project_metadata={ "AB": { "Organism": "human", },
                               "CDE": { "Organism": "mouse", } },
            top_dir=self.dirn)
        mockdir.create()
        # Settings file with fastq_strand indexes and
        # polling interval
        settings_ini = os.path.join(self.dirn,"settings.ini")
        with open(settings_ini,'w') as s:
            s.write("""[general]
poll_interval = 0.5

[fastq_strand_indexes]
human = /data/genomeIndexes/hg38/STAR
mouse = /data/genomeIndexes/mm10/STAR
""")
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn,
                         settings=Settings(settings_ini))
        # Run the QC
        status = run_qc(ap,
                        run_multiqc=True,
                        max_jobs=1)
        self.assertEqual(status,0)
        # Check the fastq_strand_conf files were created
        for p in ("AB","CDE"):
            self.assertTrue(os.path.exists(
                os.path.join(mockdir.dirn,p,"qc","fastq_strand.conf")))
        # Check fastq_strand outputs are present
        for p in ("AB","CDE"):
            fastq_strand_outputs = filter(lambda f:
                                          f.endswith("fastq_strand.txt"),
                                          os.listdir(os.path.join(
                                              mockdir.dirn,p,"qc")))
            self.assertTrue(len(fastq_strand_outputs) > 0)
        # Check output and reports
        for p in ("AB","CDE","undetermined"):
            for f in ("qc",
                      "qc_report.html",
                      "qc_report.%s.%s_analysis.zip" % (
                          p,
                          '170901_M00879_0087_000000000-AGEW9'),
                      "multiqc_report.html"):
                self.assertTrue(os.path.exists(os.path.join(mockdir.dirn,
                                                            p,f)),
                                "Missing %s in project '%s'" % (f,p))
            # Check zip file has MultiQC report
            zip_file = os.path.join(mockdir.dirn,p,
                                    "qc_report.%s.%s_analysis.zip" % (
                                        p,
                                        '170901_M00879_0087_000000000-AGEW9'))
            with zipfile.ZipFile(zip_file) as z:
                multiqc = os.path.join(
                    "qc_report.%s.%s_analysis" % (
                        p,'170901_M00879_0087_000000000-AGEW9'),
                    "multiqc_report.html")
                self.assertTrue(multiqc in z.namelist())

    def test_run_qc_single_end_with_strandedness(self):
        """run_qc: single-end QC run with strandedness determination
        """
        # Make mock illumina_qc.sh and multiqc
        MockIlluminaQcSh.create(os.path.join(self.bin,
                                             "illumina_qc.sh"))
        MockFastqStrandPy.create(os.path.join(self.bin,
                                              "fastq_strand.py"))
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            paired_end=False,
            metadata={ "instrument_datestamp": "170901" },
            project_metadata={ "AB": { "Organism": "human", },
                               "CDE": { "Organism": "mouse", } },
            top_dir=self.dirn)
        mockdir.create()
        # Settings file with fastq_strand indexes and
        # polling interval
        settings_ini = os.path.join(self.dirn,"settings.ini")
        with open(settings_ini,'w') as s:
            s.write("""[general]
poll_interval = 0.5

[fastq_strand_indexes]
human = /data/genomeIndexes/hg38/STAR
mouse = /data/genomeIndexes/mm10/STAR
""")
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn,
                         settings=Settings(settings_ini))
        # Run the QC
        status = run_qc(ap,
                        run_multiqc=True,
                        max_jobs=1)
        self.assertEqual(status,0)
        # Check the fastq_strand_conf files were created
        for p in ("AB","CDE"):
            self.assertTrue(os.path.exists(
                os.path.join(mockdir.dirn,p,"qc","fastq_strand.conf")))
        # Check fastq_strand outputs are present
        for p in ("AB","CDE"):
            fastq_strand_outputs = filter(lambda f:
                                          f.endswith("fastq_strand.txt"),
                                          os.listdir(os.path.join(
                                              mockdir.dirn,p,"qc")))
            self.assertTrue(len(fastq_strand_outputs) > 0)
        # Check output and reports
        for p in ("AB","CDE","undetermined"):
            for f in ("qc",
                      "qc_report.html",
                      "qc_report.%s.%s_analysis.zip" % (
                          p,
                          '170901_M00879_0087_000000000-AGEW9'),
                      "multiqc_report.html"):
                self.assertTrue(os.path.exists(os.path.join(mockdir.dirn,
                                                            p,f)),
                                "Missing %s in project '%s'" % (f,p))
            # Check zip file has MultiQC report
            zip_file = os.path.join(mockdir.dirn,p,
                                    "qc_report.%s.%s_analysis.zip" % (
                                        p,
                                        '170901_M00879_0087_000000000-AGEW9'))
            with zipfile.ZipFile(zip_file) as z:
                multiqc = os.path.join(
                    "qc_report.%s.%s_analysis" % (
                        p,'170901_M00879_0087_000000000-AGEW9'),
                    "multiqc_report.html")
                self.assertTrue(multiqc in z.namelist())

    def test_run_qc_single_cell_with_strandedness(self):
        """run_qc: ICELL8 scRNA-seq run with strandedness determination
        """
        # Make mock illumina_qc.sh and multiqc
        MockIlluminaQcSh.create(os.path.join(self.bin,
                                             "illumina_qc.sh"))
        MockFastqStrandPy.create(os.path.join(self.bin,
                                              "fastq_strand.py"))
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            project_metadata={ "AB": { "Organism": "human",
                                       "Single cell platform": "ICELL8", },
                               "CDE": { "Organism": "mouse",
                                        "Single cell platform": "ICELL8", } },
            top_dir=self.dirn)
        mockdir.create()
        # Settings file with fastq_strand indexes and
        # polling interval
        settings_ini = os.path.join(self.dirn,"settings.ini")
        with open(settings_ini,'w') as s:
            s.write("""[general]
poll_interval = 0.5

[fastq_strand_indexes]
human = /data/genomeIndexes/hg38/STAR
mouse = /data/genomeIndexes/mm10/STAR
""")
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn,
                         settings=Settings(settings_ini))
        # Run the QC
        status = run_qc(ap,
                        run_multiqc=True,
                        max_jobs=1)
        self.assertEqual(status,0)
        # Check the fastq_strand_conf files were created
        for p in ("AB","CDE"):
            self.assertTrue(os.path.exists(
                os.path.join(mockdir.dirn,p,"qc","fastq_strand.conf")))
        # Check fastq_strand outputs are present
        for p in ("AB","CDE"):
            fastq_strand_outputs = filter(lambda f:
                                          f.endswith("fastq_strand.txt"),
                                          os.listdir(os.path.join(
                                              mockdir.dirn,p,"qc")))
            self.assertTrue(len(fastq_strand_outputs) > 0)
        # Check output and reports
        for p in ("AB","CDE","undetermined"):
            for f in ("qc",
                      "qc_report.html",
                      "qc_report.%s.%s_analysis.zip" % (
                          p,
                          '170901_M00879_0087_000000000-AGEW9'),
                      "multiqc_report.html"):
                self.assertTrue(os.path.exists(os.path.join(mockdir.dirn,
                                                            p,f)),
                                "Missing %s in project '%s'" % (f,p))
            # Check zip file has MultiQC report
            zip_file = os.path.join(mockdir.dirn,p,
                                    "qc_report.%s.%s_analysis.zip" % (
                                        p,
                                        '170901_M00879_0087_000000000-AGEW9'))
            with zipfile.ZipFile(zip_file) as z:
                multiqc = os.path.join(
                    "qc_report.%s.%s_analysis" % (
                        p,'170901_M00879_0087_000000000-AGEW9'),
                    "multiqc_report.html")
                self.assertTrue(multiqc in z.namelist())

    def test_run_qc_10x_scRNAseq(self):
        """run_qc: 10x scRNA-seq with strandedness and single library analysis
        """
        # Make mock illumina_qc.sh and multiqc
        MockIlluminaQcSh.create(os.path.join(self.bin,
                                             "illumina_qc.sh"))
        MockFastqStrandPy.create(os.path.join(self.bin,
                                              "fastq_strand.py"))
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"))
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            project_metadata={ "AB": { "Organism": "human",
                                       "Single cell platform":
                                       "10xGenomics Chromium 3'v3",
                                       "Library type": "scRNA-seq", },
                               "CDE": { "Organism": "mouse",
                                        "Single cell platform":
                                        "10xGenomics Chromium 3'v3",
                                        "Library type": "scRNA-seq", } },
            top_dir=self.dirn)
        mockdir.create()
        # Settings file with fastq_strand indexes and
        # polling interval
        settings_ini = os.path.join(self.dirn,"settings.ini")
        with open(settings_ini,'w') as s:
            s.write("""[general]
poll_interval = 1.0

[fastq_strand_indexes]
human = /data/genomeIndexes/hg38/STAR
mouse = /data/genomeIndexes/mm10/STAR

[10xgenomics_transcriptomes]
human = /data/cellranger/transcriptomes/hg38
mouse = /data/cellranger/transcriptomes/mm10
""")
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn,
                         settings=Settings(settings_ini))
        # Run the QC
        status = run_qc(ap,
                        run_multiqc=True,
                        max_jobs=1)
        self.assertEqual(status,0)
        # Check the fastq_strand_conf files were created
        for p in ("AB","CDE"):
            self.assertTrue(os.path.exists(
                os.path.join(mockdir.dirn,p,"qc","fastq_strand.conf")))
        # Check fastq_strand outputs are present
        for p in ("AB","CDE"):
            fastq_strand_outputs = filter(lambda f:
                                          f.endswith("fastq_strand.txt"),
                                          os.listdir(os.path.join(
                                              mockdir.dirn,p,"qc")))
            self.assertTrue(len(fastq_strand_outputs) > 0)
        # Check cellranger count outputs are present
        for p in ("AB","CDE"):
            self.assertTrue(os.path.exists(os.path.join(mockdir.dirn,
                                                        p,
                                                        "cellranger_count")))
        # Check output and reports
        for p in ("AB","CDE","undetermined"):
            for f in ("qc",
                      "qc_report.html",
                      "qc_report.%s.%s_analysis.zip" % (
                          p,
                          '170901_M00879_0087_000000000-AGEW9'),
                      "multiqc_report.html"):
                self.assertTrue(os.path.exists(os.path.join(mockdir.dirn,
                                                            p,f)),
                                "Missing %s in project '%s'" % (f,p))
            # Check zip file has MultiQC report
            zip_file = os.path.join(mockdir.dirn,p,
                                    "qc_report.%s.%s_analysis.zip" % (
                                        p,
                                        '170901_M00879_0087_000000000-AGEW9'))
            with zipfile.ZipFile(zip_file) as z:
                multiqc = os.path.join(
                    "qc_report.%s.%s_analysis" % (
                        p,'170901_M00879_0087_000000000-AGEW9'),
                    "multiqc_report.html")
                self.assertTrue(multiqc in z.namelist())

    def test_run_qc_10x_scATACseq(self):
        """run_qc: 10x scATAC-seq with strandedness and single library analysis
        """
        # Make mock illumina_qc.sh and multiqc
        MockIlluminaQcSh.create(os.path.join(self.bin,
                                             "illumina_qc.sh"))
        MockFastqStrandPy.create(os.path.join(self.bin,
                                              "fastq_strand.py"))
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger-atac"))
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            project_metadata={ "AB": { "Organism": "human",
                                       "Single cell platform":
                                       "10xGenomics Single Cell ATAC",
                                       "Library type": "scATAC-seq", },
                               "CDE": { "Organism": "mouse",
                                        "Single cell platform":
                                        "10xGenomics Single Cell ATAC",
                                        "Library type": "scATAC-seq", } },
            reads=('R1','R2','R3','I1'),
            top_dir=self.dirn)
        mockdir.create()
        # Settings file with fastq_strand indexes and
        # polling interval
        settings_ini = os.path.join(self.dirn,"settings.ini")
        with open(settings_ini,'w') as s:
            s.write("""[general]
poll_interval = 1.0

[fastq_strand_indexes]
human = /data/genomeIndexes/hg38/STAR
mouse = /data/genomeIndexes/mm10/STAR

[10xgenomics_atac_genome_references]
human = /data/cellranger/atac_references/hg38
mouse = /data/cellranger/atac_references/mm10
""")
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn,
                         settings=Settings(settings_ini))
        # Run the QC
        status = run_qc(ap,
                        run_multiqc=True,
                        max_jobs=1)
        self.assertEqual(status,0)
        # Check the fastq_strand_conf files were created
        for p in ("AB","CDE"):
            self.assertTrue(os.path.exists(
                os.path.join(mockdir.dirn,p,"qc","fastq_strand.conf")))
        # Check fastq_strand outputs are present
        for p in ("AB","CDE"):
            fastq_strand_outputs = filter(lambda f:
                                          f.endswith("fastq_strand.txt"),
                                          os.listdir(os.path.join(
                                              mockdir.dirn,p,"qc")))
            self.assertTrue(len(fastq_strand_outputs) > 0)
        # Check cellranger count outputs are present
        for p in ("AB","CDE"):
            self.assertTrue(os.path.exists(os.path.join(mockdir.dirn,
                                                        p,
                                                        "cellranger_count")))
        # Check output and reports
        for p in ("AB","CDE","undetermined"):
            for f in ("qc",
                      "qc_report.html",
                      "qc_report.%s.%s_analysis.zip" % (
                          p,
                          '170901_M00879_0087_000000000-AGEW9'),
                      "multiqc_report.html"):
                self.assertTrue(os.path.exists(os.path.join(mockdir.dirn,
                                                            p,f)),
                                "Missing %s in project '%s'" % (f,p))
            # Check zip file has MultiQC report
            zip_file = os.path.join(mockdir.dirn,p,
                                    "qc_report.%s.%s_analysis.zip" % (
                                        p,
                                        '170901_M00879_0087_000000000-AGEW9'))
            with zipfile.ZipFile(zip_file) as z:
                multiqc = os.path.join(
                    "qc_report.%s.%s_analysis" % (
                        p,'170901_M00879_0087_000000000-AGEW9'),
                    "multiqc_report.html")
                self.assertTrue(multiqc in z.namelist())
