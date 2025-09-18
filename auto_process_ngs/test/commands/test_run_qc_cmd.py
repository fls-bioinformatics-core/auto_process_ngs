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
from auto_process_ngs.metadata import AnalysisProjectQCDirInfo
from auto_process_ngs.mock import MockAnalysisDirFactory
from auto_process_ngs.mock import MockFastqScreen
from auto_process_ngs.mock import MockFastQC
from auto_process_ngs.mock import MockGtf2bed
from auto_process_ngs.mock import MockSeqtk
from auto_process_ngs.mock import MockStar
from auto_process_ngs.mock import MockSamtools
from auto_process_ngs.mock import MockPicard
from auto_process_ngs.mock import MockRSeQC
from auto_process_ngs.mock import MockQualimap
from auto_process_ngs.mock import MockCellrangerExe
from auto_process_ngs.mock import MockMultiQC
from auto_process_ngs.qc.outputs import QCOutputs
from auto_process_ngs.settings import Settings
from auto_process_ngs.commands.run_qc_cmd import run_qc

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

# Polling interval for pipeline
POLL_INTERVAL = 0.1

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
        # Create a temp 'data' dir
        self.data = os.path.join(self.dirn,"data")
        os.mkdir(self.data)
        # Add (empty) FastqScreen conf files
        self.fastq_screens = {
            'model_organisms': "fastq_screen_model_organisms.conf",
            'other_organisms': "fastq_screen_other_organisms.conf",
            'rRNA': "fastq_screen_rRNA.conf",
        }
        for screen in self.fastq_screens:
            conf_file = os.path.join(self.data,
                                     self.fastq_screens[screen])
            with open(conf_file,'wt') as fp:
                fp.write("")
            self.fastq_screens[screen] = conf_file
        # Add (empty) reference data files
        self.ref_data = dict()
        for build in ('hg38','mm10',):
            self.ref_data[build] = {}
            build_dir = os.path.join(self.data,build)
            os.mkdir(build_dir)
            for ext in ('bed','gtf'):
                f = os.path.join(build_dir,"%s.%s" % (build,ext))
                with open(f,'wt') as fp:
                    fp.write("")
                self.ref_data[build][ext] = f
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
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        MockStar.create(os.path.join(self.bin,"STAR"))
        MockSamtools.create(os.path.join(self.bin,"samtools"))
        MockPicard.create(os.path.join(self.bin,"picard"))
        MockGtf2bed.create(os.path.join(self.bin,"gtf2bed"))
        MockRSeQC.create(os.path.join(self.bin,"infer_experiment.py"))
        MockRSeQC.create(os.path.join(self.bin,"geneBody_coverage.py"))
        MockQualimap.create(os.path.join(self.bin,"qualimap"))
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901", },
            project_metadata={ "AB": { "Library type": "RNA-seq",
                                       "Organism": "human", },
                               "CDE": { "Library type": "ChIP-seq",
                                        "Organism": "mouse", } },
            top_dir=self.dirn)
        mockdir.create()
        # Settings file with reference data and polling interval
        settings_ini = os.path.join(self.dirn,"auto_process.ini")
        with open(settings_ini,'w') as s:
            s.write("""[general]
poll_interval = {poll_interval}

[organism:human]
star_index = /data/hg38/star_index
annotation_bed = {hg38_bed}
annotation_gtf = {hg38_gtf}

[organism:mouse]
star_index = /data/mm10/star_index
annotation_bed = {mm10_bed}
annotation_gtf = {mm10_gtf}
""".format(hg38_bed=self.ref_data['hg38']['bed'],
           hg38_gtf=self.ref_data['hg38']['gtf'],
           mm10_bed=self.ref_data['mm10']['bed'],
           mm10_gtf=self.ref_data['mm10']['gtf'],
           poll_interval=POLL_INTERVAL))
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn,
                         settings=Settings(settings_ini))
        # Run the QC
        status = run_qc(ap,
                        fastq_screens=self.fastq_screens,
                        run_multiqc=True,
                        max_jobs=1)
        self.assertEqual(status,0)
        # Check detected outputs
        for p in ("AB","CDE"):
            qc_dir = os.path.join(mockdir.dirn,p,"qc")
            qcoutputs = QCOutputs(qc_dir)
            for qc_module in ("fastqc_r1",
                              "fastqc_r2",
                              "screens_r1",
                              "screens_r2",
                              "sequence_lengths",
                              "picard_insert_size_metrics",
                              "rseqc_genebody_coverage",
                              "rseqc_infer_experiment",
                              "qualimap_rnaseq",
                              "multiqc"):
                self.assertTrue(qc_module in qcoutputs.outputs,
                                "Project '%s': missing output '%s'" %
                                (p,qc_module))
        # Check output and reports
        for p in ("AB","CDE","undetermined"):
            for f in ("qc",
                      "qc_report.html",
                      "qc_report.%s.%s.zip" % (
                          p,
                          '170901_M00879_0087_000000000-AGEW9'),
                      "multiqc_report.html"):
                self.assertTrue(os.path.exists(os.path.join(mockdir.dirn,
                                                            p,f)),
                                "Missing %s in project '%s'" % (f,p))
            # Check zip file has MultiQC report
            zip_file = os.path.join(mockdir.dirn,p,
                                    "qc_report.%s.%s.zip" % (
                                        p,
                                        '170901_M00879_0087_000000000-AGEW9'))
            with zipfile.ZipFile(zip_file) as z:
                multiqc = os.path.join(
                    "qc_report.%s.%s" % (
                        p,'170901_M00879_0087_000000000-AGEW9'),
                    "multiqc_report.html")
                self.assertTrue(multiqc in z.namelist())
        # Check lane splitting for undetermined
        qc_info = AnalysisProjectQCDirInfo(os.path.join(mockdir.dirn,
                                                        "undetermined",
                                                        "qc",
                                                        "qc.info"))
        self.assertTrue(qc_info.fastqs_split_by_lane)

    def test_run_qc_single_end(self):
        """run_qc: single-end QC run
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        MockStar.create(os.path.join(self.bin,"STAR"))
        MockSamtools.create(os.path.join(self.bin,"samtools"))
        MockGtf2bed.create(os.path.join(self.bin,"gtf2bed"))
        MockRSeQC.create(os.path.join(self.bin,"infer_experiment.py"))
        MockRSeQC.create(os.path.join(self.bin,"geneBody_coverage.py"))
        MockQualimap.create(os.path.join(self.bin,"qualimap"))
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            paired_end=False,
            metadata={ "instrument_datestamp": "170901" },
            project_metadata={ "AB": { "Library type": "RNA-seq",
                                       "Organism": "human", },
                               "CDE": { "Library type": "ChIP-seq",
                                        "Organism": "mouse", } },
            top_dir=self.dirn)
        mockdir.create()
        # Settings file with reference data and polling interval
        settings_ini = os.path.join(self.dirn,"auto_process.ini")
        with open(settings_ini,'w') as s:
            s.write("""[general]
poll_interval = {poll_interval}

[organism:human]
star_index = /data/genomeIndexes/hg38/STAR
annotation_bed = {hg38_bed}
annotation_gtf = {hg38_gtf}

[organism:mouse]
star_index = /data/genomeIndexes/mm10/STAR
annotation_bed = {mm10_bed}
annotation_gtf = {mm10_gtf}
""".format(hg38_bed=self.ref_data['hg38']['bed'],
           hg38_gtf=self.ref_data['hg38']['gtf'],
           mm10_bed=self.ref_data['mm10']['bed'],
           mm10_gtf=self.ref_data['mm10']['gtf'],
           poll_interval=POLL_INTERVAL))
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn,
                         settings=Settings(settings_ini))
        # Run the QC
        status = run_qc(ap,
                        fastq_screens=self.fastq_screens,
                        run_multiqc=True,
                        max_jobs=1)
        self.assertEqual(status,0)
        # Check detected outputs
        for p in ("AB","CDE"):
            qc_dir = os.path.join(mockdir.dirn,p,"qc")
            qcoutputs = QCOutputs(qc_dir)
            for qc_module in ("fastqc_r1",
                              "screens_r1",
                              "sequence_lengths",
                              "rseqc_genebody_coverage",
                              "rseqc_infer_experiment",
                              "qualimap_rnaseq",
                              "multiqc"):
                self.assertTrue(qc_module in qcoutputs.outputs,
                                "Project '%s': missing output '%s'" %
                                (p,qc_module))
        # Check output and reports
        for p in ("AB","CDE","undetermined"):
            for f in ("qc",
                      "qc_report.html",
                      "qc_report.%s.%s.zip" % (
                          p,
                          '170901_M00879_0087_000000000-AGEW9'),
                      "multiqc_report.html"):
                self.assertTrue(os.path.exists(os.path.join(mockdir.dirn,
                                                            p,f)),
                                "Missing %s in project '%s'" % (f,p))
            # Check zip file has MultiQC report
            zip_file = os.path.join(mockdir.dirn,p,
                                    "qc_report.%s.%s.zip" % (
                                        p,
                                        '170901_M00879_0087_000000000-AGEW9'))
            with zipfile.ZipFile(zip_file) as z:
                multiqc = os.path.join(
                    "qc_report.%s.%s" % (
                        p,'170901_M00879_0087_000000000-AGEW9'),
                    "multiqc_report.html")
                self.assertTrue(multiqc in z.namelist())
        # Check lane splitting for undetermined
        qc_info = AnalysisProjectQCDirInfo(os.path.join(mockdir.dirn,
                                                        "undetermined",
                                                        "qc",
                                                        "qc.info"))
        self.assertTrue(qc_info.fastqs_split_by_lane)

    def test_run_qc_10x_scRNAseq(self):
        """run_qc: 10x scRNA-seq with single library analysis
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        MockStar.create(os.path.join(self.bin,"STAR"))
        MockSamtools.create(os.path.join(self.bin,"samtools"))
        MockGtf2bed.create(os.path.join(self.bin,"gtf2bed"))
        MockRSeQC.create(os.path.join(self.bin,"infer_experiment.py"))
        MockRSeQC.create(os.path.join(self.bin,"geneBody_coverage.py"))
        MockQualimap.create(os.path.join(self.bin,"qualimap"))
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
        # Settings file with reference data and polling interval
        settings_ini = os.path.join(self.dirn,"auto_process.ini")
        with open(settings_ini,'w') as s:
            s.write("""[general]
poll_interval = {poll_interval}

[organism:human]
star_index = /data/genomeIndexes/hg38/STAR
annotation_bed = {hg38_bed}
annotation_gtf = {hg38_gtf}
cellranger_reference = /data/cellranger/transcriptomes/hg38

[organism:mouse]
star_index = /data/genomeIndexes/mm10/STAR
annotation_bed = {mm10_bed}
annotation_gtf = {mm10_gtf}
cellranger_reference = /data/cellranger/transcriptomes/mm10
""".format(hg38_bed=self.ref_data['hg38']['bed'],
           hg38_gtf=self.ref_data['hg38']['gtf'],
           mm10_bed=self.ref_data['mm10']['bed'],
           mm10_gtf=self.ref_data['mm10']['gtf'],
           poll_interval=POLL_INTERVAL))
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn,
                         settings=Settings(settings_ini))
        # Run the QC
        status = run_qc(ap,
                        fastq_screens=self.fastq_screens,
                        run_multiqc=True,
                        max_jobs=1)
        self.assertEqual(status,0)
        # Check detected outputs
        for p in ("AB","CDE"):
            qc_dir = os.path.join(mockdir.dirn,p,"qc")
            qcoutputs = QCOutputs(qc_dir)
            print(qcoutputs.outputs)
            for qc_module in ("fastqc_r1",
                              "fastqc_r2",
                              "screens_r2",
                              "sequence_lengths",
                              "rseqc_genebody_coverage",
                              "rseqc_infer_experiment",
                              "qualimap_rnaseq",
                              "cellranger_count",
                              "multiqc"):
                self.assertTrue(qc_module in qcoutputs.outputs,
                                "Project '%s': missing output '%s'" %
                                (p,qc_module))
        # Check output and reports
        for p in ("AB","CDE","undetermined"):
            for f in ("qc",
                      "qc_report.html",
                      "qc_report.%s.%s.zip" % (
                          p,
                          '170901_M00879_0087_000000000-AGEW9'),
                      "multiqc_report.html"):
                self.assertTrue(os.path.exists(os.path.join(mockdir.dirn,
                                                            p,f)),
                                "Missing %s in project '%s'" % (f,p))
            # Check zip file has MultiQC report
            zip_file = os.path.join(mockdir.dirn,p,
                                    "qc_report.%s.%s.zip" % (
                                        p,
                                        '170901_M00879_0087_000000000-AGEW9'))
            with zipfile.ZipFile(zip_file) as z:
                multiqc = os.path.join(
                    "qc_report.%s.%s" % (
                        p,'170901_M00879_0087_000000000-AGEW9'),
                    "multiqc_report.html")
                self.assertTrue(multiqc in z.namelist())
        # Check lane splitting for undetermined
        qc_info = AnalysisProjectQCDirInfo(os.path.join(mockdir.dirn,
                                                        "undetermined",
                                                        "qc",
                                                        "qc.info"))
        self.assertTrue(qc_info.fastqs_split_by_lane)

    def test_run_qc_10x_snRNAseq_310(self):
        """run_qc: 10x snRNA-seq with single library analysis (cellranger 3.1.0)
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        MockStar.create(os.path.join(self.bin,"STAR"))
        MockSamtools.create(os.path.join(self.bin,"samtools"))
        MockGtf2bed.create(os.path.join(self.bin,"gtf2bed"))
        MockRSeQC.create(os.path.join(self.bin,"infer_experiment.py"))
        MockRSeQC.create(os.path.join(self.bin,"geneBody_coverage.py"))
        MockQualimap.create(os.path.join(self.bin,"qualimap"))
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 version="3.1.0",
                                 assert_include_introns=False)
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
                                       "Library type": "snRNA-seq", },
                               "CDE": { "Organism": "mouse",
                                        "Single cell platform":
                                        "10xGenomics Chromium 3'v3",
                                        "Library type": "snRNA-seq", } },
            top_dir=self.dirn)
        mockdir.create()
        # Settings file with reference data and polling interval
        settings_ini = os.path.join(self.dirn,"auto_process.ini")
        with open(settings_ini,'w') as s:
            s.write("""[general]
poll_interval = {poll_interval}

[organism:human]
star_index = /data/genomeIndexes/hg38/STAR
annotation_bed = {hg38_bed}
annotation_gtf = {hg38_gtf}
cellranger_premrna_reference = /data/cellranger/transcriptomes/hg38_pre_mrna

[organism:mouse]
star_index = /data/genomeIndexes/mm10/STAR
annotation_bed = {mm10_bed}
annotation_gtf = {mm10_gtf}
cellranger_premrna_reference = /data/cellranger/transcriptomes/mm10_pre_mrna
""".format(hg38_bed=self.ref_data['hg38']['bed'],
           hg38_gtf=self.ref_data['hg38']['gtf'],
           mm10_bed=self.ref_data['mm10']['bed'],
           mm10_gtf=self.ref_data['mm10']['gtf'],
           poll_interval=POLL_INTERVAL))
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn,
                         settings=Settings(settings_ini))
        # Run the QC
        status = run_qc(ap,
                        fastq_screens=self.fastq_screens,
                        run_multiqc=True,
                        max_jobs=1)
        self.assertEqual(status,0)
        # Check detected outputs
        for p in ("AB","CDE"):
            qc_dir = os.path.join(mockdir.dirn,p,"qc")
            qcoutputs = QCOutputs(qc_dir)
            print(qcoutputs.outputs)
            for qc_module in ("fastqc_r1",
                              "fastqc_r2",
                              "screens_r2",
                              "sequence_lengths",
                              "rseqc_genebody_coverage",
                              "rseqc_infer_experiment",
                              "qualimap_rnaseq",
                              "cellranger_count",
                              "multiqc"):
                self.assertTrue(qc_module in qcoutputs.outputs,
                                "Project '%s': missing output '%s'" %
                                (p,qc_module))
        # Check output and reports
        for p in ("AB","CDE","undetermined"):
            for f in ("qc",
                      "qc_report.html",
                      "qc_report.%s.%s.zip" % (
                          p,
                          '170901_M00879_0087_000000000-AGEW9'),
                      "multiqc_report.html"):
                self.assertTrue(os.path.exists(os.path.join(mockdir.dirn,
                                                            p,f)),
                                "Missing %s in project '%s'" % (f,p))
            # Check zip file has MultiQC report
            zip_file = os.path.join(mockdir.dirn,p,
                                    "qc_report.%s.%s.zip" % (
                                        p,
                                        '170901_M00879_0087_000000000-AGEW9'))
            with zipfile.ZipFile(zip_file) as z:
                multiqc = os.path.join(
                    "qc_report.%s.%s" % (
                        p,'170901_M00879_0087_000000000-AGEW9'),
                    "multiqc_report.html")
                self.assertTrue(multiqc in z.namelist())
        # Check lane splitting for undetermined
        qc_info = AnalysisProjectQCDirInfo(os.path.join(mockdir.dirn,
                                                        "undetermined",
                                                        "qc",
                                                        "qc.info"))
        self.assertTrue(qc_info.fastqs_split_by_lane)

    def test_run_qc_10x_snRNAseq_501(self):
        """run_qc: 10x snRNA-seq with single library analysis (cellranger 5.0.1)
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        MockStar.create(os.path.join(self.bin,"STAR"))
        MockSamtools.create(os.path.join(self.bin,"samtools"))
        MockGtf2bed.create(os.path.join(self.bin,"gtf2bed"))
        MockRSeQC.create(os.path.join(self.bin,"infer_experiment.py"))
        MockRSeQC.create(os.path.join(self.bin,"geneBody_coverage.py"))
        MockQualimap.create(os.path.join(self.bin,"qualimap"))
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 version="5.0.1",
                                 assert_include_introns=True)
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
                                       "Library type": "snRNA-seq", },
                               "CDE": { "Organism": "mouse",
                                        "Single cell platform":
                                        "10xGenomics Chromium 3'v3",
                                        "Library type": "snRNA-seq", } },
            top_dir=self.dirn)
        mockdir.create()
        # Settings file with fastq_strand indexes and
        # polling interval
        settings_ini = os.path.join(self.dirn,"auto_process.ini")
        with open(settings_ini,'w') as s:
            s.write("""[general]
poll_interval = {poll_interval}

[organism:human]
star_index = /data/genomeIndexes/hg38/STAR
annotation_bed = {hg38_bed}
annotation_gtf = {hg38_gtf}
cellranger_reference = /data/cellranger/transcriptomes/hg38

[organism:mouse]
star_index = /data/genomeIndexes/mm10/STAR
annotation_bed = {mm10_bed}
annotation_gtf = {mm10_gtf}
cellranger_reference = /data/cellranger/transcriptomes/mm10
""".format(hg38_bed=self.ref_data['hg38']['bed'],
           hg38_gtf=self.ref_data['hg38']['gtf'],
           mm10_bed=self.ref_data['mm10']['bed'],
           mm10_gtf=self.ref_data['mm10']['gtf'],
           poll_interval=POLL_INTERVAL))
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn,
                         settings=Settings(settings_ini))
        # Run the QC
        status = run_qc(ap,
                        fastq_screens=self.fastq_screens,
                        run_multiqc=True,
                        max_jobs=1)
        self.assertEqual(status,0)
        # Check detected outputs
        for p in ("AB","CDE"):
            qc_dir = os.path.join(mockdir.dirn,p,"qc")
            qcoutputs = QCOutputs(qc_dir)
            print(qcoutputs.outputs)
            for qc_module in ("fastqc_r1",
                              "fastqc_r2",
                              "screens_r2",
                              "sequence_lengths",
                              "rseqc_genebody_coverage",
                              "rseqc_infer_experiment",
                              "qualimap_rnaseq",
                              "cellranger_count",
                              "multiqc"):
                self.assertTrue(qc_module in qcoutputs.outputs,
                                "Project '%s': missing output '%s'" %
                                (p,qc_module))
        # Check output and reports
        for p in ("AB","CDE","undetermined"):
            for f in ("qc",
                      "qc_report.html",
                      "qc_report.%s.%s.zip" % (
                          p,
                          '170901_M00879_0087_000000000-AGEW9'),
                      "multiqc_report.html"):
                self.assertTrue(os.path.exists(os.path.join(mockdir.dirn,
                                                            p,f)),
                                "Missing %s in project '%s'" % (f,p))
            # Check zip file has MultiQC report
            zip_file = os.path.join(mockdir.dirn,p,
                                    "qc_report.%s.%s.zip" % (
                                        p,
                                        '170901_M00879_0087_000000000-AGEW9'))
            with zipfile.ZipFile(zip_file) as z:
                multiqc = os.path.join(
                    "qc_report.%s.%s" % (
                        p,'170901_M00879_0087_000000000-AGEW9'),
                    "multiqc_report.html")
                self.assertTrue(multiqc in z.namelist())
        # Check lane splitting for undetermined
        qc_info = AnalysisProjectQCDirInfo(os.path.join(mockdir.dirn,
                                                        "undetermined",
                                                        "qc",
                                                        "qc.info"))
        self.assertTrue(qc_info.fastqs_split_by_lane)

    def test_run_qc_10x_snRNAseq_600(self):
        """run_qc: 10x snRNA-seq with single library analysis (cellranger 6.0.0)
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        MockStar.create(os.path.join(self.bin,"STAR"))
        MockSamtools.create(os.path.join(self.bin,"samtools"))
        MockGtf2bed.create(os.path.join(self.bin,"gtf2bed"))
        MockRSeQC.create(os.path.join(self.bin,"infer_experiment.py"))
        MockRSeQC.create(os.path.join(self.bin,"geneBody_coverage.py"))
        MockQualimap.create(os.path.join(self.bin,"qualimap"))
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 version="6.0.0",
                                 assert_include_introns=True)
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
                                       "Library type": "snRNA-seq", },
                               "CDE": { "Organism": "mouse",
                                        "Single cell platform":
                                        "10xGenomics Chromium 3'v3",
                                        "Library type": "snRNA-seq", } },
            top_dir=self.dirn)
        mockdir.create()
        # Settings file with fastq_strand indexes and
        # polling interval
        settings_ini = os.path.join(self.dirn,"auto_process.ini")
        with open(settings_ini,'w') as s:
            s.write("""[general]
poll_interval = {poll_interval}

[organism:human]
star_index = /data/genomeIndexes/hg38/STAR
annotation_bed = {hg38_bed}
annotation_gtf = {hg38_gtf}
cellranger_reference = /data/cellranger/transcriptomes/hg38

[organism:mouse]
star_index = /data/genomeIndexes/mm10/STAR
annotation_bed = {mm10_bed}
annotation_gtf = {mm10_gtf}
cellranger_reference = /data/cellranger/transcriptomes/mm10
""".format(hg38_bed=self.ref_data['hg38']['bed'],
           hg38_gtf=self.ref_data['hg38']['gtf'],
           mm10_bed=self.ref_data['mm10']['bed'],
           mm10_gtf=self.ref_data['mm10']['gtf'],
           poll_interval=POLL_INTERVAL))
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn,
                         settings=Settings(settings_ini))
        # Run the QC
        status = run_qc(ap,
                        fastq_screens=self.fastq_screens,
                        run_multiqc=True,
                        max_jobs=1)
        self.assertEqual(status,0)
        # Check detected outputs
        for p in ("AB","CDE"):
            qc_dir = os.path.join(mockdir.dirn,p,"qc")
            qcoutputs = QCOutputs(qc_dir)
            print(qcoutputs.outputs)
            for qc_module in ("fastqc_r1",
                              "fastqc_r2",
                              "screens_r2",
                              "sequence_lengths",
                              "rseqc_genebody_coverage",
                              "rseqc_infer_experiment",
                              "qualimap_rnaseq",
                              "cellranger_count",
                              "multiqc"):
                self.assertTrue(qc_module in qcoutputs.outputs,
                                "Project '%s': missing output '%s'" %
                                (p,qc_module))
        # Check output and reports
        for p in ("AB","CDE","undetermined"):
            for f in ("qc",
                      "qc_report.html",
                      "qc_report.%s.%s.zip" % (
                          p,
                          '170901_M00879_0087_000000000-AGEW9'),
                      "multiqc_report.html"):
                self.assertTrue(os.path.exists(os.path.join(mockdir.dirn,
                                                            p,f)),
                                "Missing %s in project '%s'" % (f,p))
            # Check zip file has MultiQC report
            zip_file = os.path.join(mockdir.dirn,p,
                                    "qc_report.%s.%s.zip" % (
                                        p,
                                        '170901_M00879_0087_000000000-AGEW9'))
            with zipfile.ZipFile(zip_file) as z:
                multiqc = os.path.join(
                    "qc_report.%s.%s" % (
                        p,'170901_M00879_0087_000000000-AGEW9'),
                    "multiqc_report.html")
                self.assertTrue(multiqc in z.namelist())
        # Check lane splitting for undetermined
        qc_info = AnalysisProjectQCDirInfo(os.path.join(mockdir.dirn,
                                                        "undetermined",
                                                        "qc",
                                                        "qc.info"))
        self.assertTrue(qc_info.fastqs_split_by_lane)

    def test_run_qc_10x_scATACseq(self):
        """run_qc: 10x scATAC-seq with single library analysis
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        MockStar.create(os.path.join(self.bin,"STAR"))
        MockSamtools.create(os.path.join(self.bin,"samtools"))
        MockPicard.create(os.path.join(self.bin,"picard"))
        MockGtf2bed.create(os.path.join(self.bin,"gtf2bed"))
        MockRSeQC.create(os.path.join(self.bin,"infer_experiment.py"))
        MockRSeQC.create(os.path.join(self.bin,"geneBody_coverage.py"))
        MockQualimap.create(os.path.join(self.bin,"qualimap"))
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
        settings_ini = os.path.join(self.dirn,"auto_process.ini")
        with open(settings_ini,'w') as s:
            s.write("""[general]
poll_interval = {poll_interval}

[organism:human]
star_index = /data/genomeIndexes/hg38/STAR
annotation_bed = {hg38_bed}
annotation_gtf = {hg38_gtf}
cellranger_atac_reference = /data/cellranger/atac_references/hg38

[organism:mouse]
star_index = /data/genomeIndexes/mm10/STAR
annotation_bed = {mm10_bed}
annotation_gtf = {mm10_gtf}
cellranger_atac_reference = /data/cellranger/atac_references/mm10
""".format(hg38_bed=self.ref_data['hg38']['bed'],
           hg38_gtf=self.ref_data['hg38']['gtf'],
           mm10_bed=self.ref_data['mm10']['bed'],
           mm10_gtf=self.ref_data['mm10']['gtf'],
           poll_interval=POLL_INTERVAL))
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn,
                         settings=Settings(settings_ini))
        # Run the QC
        status = run_qc(ap,
                        fastq_screens=self.fastq_screens,
                        run_multiqc=True,
                        max_jobs=1)
        self.assertEqual(status,0)
        # Check detected outputs
        for p in ("AB","CDE"):
            qc_dir = os.path.join(mockdir.dirn,p,"qc")
            qcoutputs = QCOutputs(qc_dir)
            print(qcoutputs.outputs)
            for qc_module in ("fastqc_r1",
                              "fastqc_r3",
                              "screens_r1",
                              "screens_r3",
                              "sequence_lengths",
                              "picard_insert_size_metrics",
                              "rseqc_genebody_coverage",
                              "rseqc_infer_experiment",
                              "qualimap_rnaseq",
                              "cellranger-atac_count",
                              "multiqc"):
                self.assertTrue(qc_module in qcoutputs.outputs,
                                "Project '%s': missing output '%s'" %
                                (p,qc_module))
        # Check output and reports
        for p in ("AB","CDE","undetermined"):
            for f in ("qc",
                      "qc_report.html",
                      "qc_report.%s.%s.zip" % (
                          p,
                          '170901_M00879_0087_000000000-AGEW9'),
                      "multiqc_report.html"):
                self.assertTrue(os.path.exists(os.path.join(mockdir.dirn,
                                                            p,f)),
                                "Missing %s in project '%s'" % (f,p))
            # Check zip file has MultiQC report
            zip_file = os.path.join(mockdir.dirn,p,
                                    "qc_report.%s.%s.zip" % (
                                        p,
                                        '170901_M00879_0087_000000000-AGEW9'))
            with zipfile.ZipFile(zip_file) as z:
                multiqc = os.path.join(
                    "qc_report.%s.%s" % (
                        p,'170901_M00879_0087_000000000-AGEW9'),
                    "multiqc_report.html")
                self.assertTrue(multiqc in z.namelist())
        # Check lane splitting for undetermined
        qc_info = AnalysisProjectQCDirInfo(os.path.join(mockdir.dirn,
                                                        "undetermined",
                                                        "qc",
                                                        "qc.info"))
        self.assertTrue(qc_info.fastqs_split_by_lane)

    def test_run_qc_10x_visium_ffpe_spatial_gex(self):
        """run_qc: 10x Visium FFPE spatial GEX
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        MockStar.create(os.path.join(self.bin,"STAR"))
        MockSamtools.create(os.path.join(self.bin,"samtools"))
        MockGtf2bed.create(os.path.join(self.bin,"gtf2bed"))
        MockSeqtk.create(os.path.join(self.bin,"seqtk"))
        MockRSeQC.create(os.path.join(self.bin,"infer_experiment.py"))
        MockRSeQC.create(os.path.join(self.bin,"geneBody_coverage.py"))
        MockQualimap.create(os.path.join(self.bin,"qualimap"))
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
                                       "10xGenomics Visium",
                                       "Library type": "RNA-seq", },
                               "CDE": { "Organism": "mouse",
                                        "Single cell platform":
                                        "10xGenomics Visium",
                                        "Library type": "RNA-seq", } },
            top_dir=self.dirn)
        mockdir.create()
        # Settings file with fastq_strand indexes and
        # polling interval
        settings_ini = os.path.join(self.dirn,"auto_process.ini")
        with open(settings_ini,'w') as s:
            s.write("""[general]
poll_interval = {poll_interval}

[organism:human]
star_index = /data/genomeIndexes/hg38/STAR
annotation_bed = {hg38_bed}
annotation_gtf = {hg38_gtf}

[organism:mouse]
star_index = /data/genomeIndexes/mm10/STAR
annotation_bed = {mm10_bed}
annotation_gtf = {mm10_gtf}
""".format(hg38_bed=self.ref_data['hg38']['bed'],
           hg38_gtf=self.ref_data['hg38']['gtf'],
           mm10_bed=self.ref_data['mm10']['bed'],
           mm10_gtf=self.ref_data['mm10']['gtf'],
           poll_interval=POLL_INTERVAL))
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn,
                         settings=Settings(settings_ini))
        # Run the QC
        status = run_qc(ap,
                        fastq_screens=self.fastq_screens,
                        run_multiqc=True,
                        max_jobs=1)
        self.assertEqual(status,0)
        # Check detected outputs
        for p in ("AB","CDE"):
            qc_dir = os.path.join(mockdir.dirn,p,"qc")
            qcoutputs = QCOutputs(qc_dir)
            print(qcoutputs.outputs)
            for qc_module in ("fastqc_r1",
                              "fastqc_r2",
                              "screens_r2",
                              "sequence_lengths",
                              "rseqc_genebody_coverage",
                              "rseqc_infer_experiment",
                              "qualimap_rnaseq",
                              "multiqc"):
                self.assertTrue(qc_module in qcoutputs.outputs,
                                "Project '%s': missing output '%s'" %
                                (p,qc_module))
        # Check output and reports
        for p in ("AB","CDE","undetermined"):
            for f in ("qc",
                      "qc_report.html",
                      "qc_report.%s.%s.zip" % (
                          p,
                          '170901_M00879_0087_000000000-AGEW9'),
                      "multiqc_report.html"):
                self.assertTrue(os.path.exists(os.path.join(mockdir.dirn,
                                                            p,f)),
                                "Missing %s in project '%s'" % (f,p))
            # Check zip file has MultiQC report
            zip_file = os.path.join(mockdir.dirn,p,
                                    "qc_report.%s.%s.zip" % (
                                        p,
                                        '170901_M00879_0087_000000000-AGEW9'))
            with zipfile.ZipFile(zip_file) as z:
                multiqc = os.path.join(
                    "qc_report.%s.%s" % (
                        p,'170901_M00879_0087_000000000-AGEW9'),
                    "multiqc_report.html")
                self.assertTrue(multiqc in z.namelist())
        # Check lane splitting for undetermined
        qc_info = AnalysisProjectQCDirInfo(os.path.join(mockdir.dirn,
                                                        "undetermined",
                                                        "qc",
                                                        "qc.info"))
        self.assertTrue(qc_info.fastqs_split_by_lane)

    def test_run_qc_10x_multiome_atac(self):
        """run_qc: 10x Multiome ATAC
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        MockStar.create(os.path.join(self.bin,"STAR"))
        MockSamtools.create(os.path.join(self.bin,"samtools"))
        MockPicard.create(os.path.join(self.bin,"picard"))
        MockGtf2bed.create(os.path.join(self.bin,"gtf2bed"))
        MockRSeQC.create(os.path.join(self.bin,"infer_experiment.py"))
        MockRSeQC.create(os.path.join(self.bin,"geneBody_coverage.py"))
        MockQualimap.create(os.path.join(self.bin,"qualimap"))
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger-arc"))
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
                                       "10xGenomics Single Cell Multiome",
                                       "Library type": "ATAC", } },
            reads=('R1','R2','R3','I1'),
            top_dir=self.dirn)
        mockdir.create()
        # Settings file with fastq_strand indexes and
        # polling interval
        settings_ini = os.path.join(self.dirn,"auto_process.ini")
        with open(settings_ini,'w') as s:
            s.write("""[general]
poll_interval = {poll_interval}

[organism:human]
star_index = /data/genomeIndexes/hg38/STAR
annotation_bed = {hg38_bed}
annotation_gtf = {hg38_gtf}
cellranger_arc_reference = /data/cellranger/arc_references/hg38-arc
cellranger_atac_reference = /data/cellranger/atac_references/hg38-atac
""".format(hg38_bed=self.ref_data['hg38']['bed'],
           hg38_gtf=self.ref_data['hg38']['gtf'],
           poll_interval=POLL_INTERVAL))
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn,
                         settings=Settings(settings_ini))
        # Run the QC
        status = run_qc(ap,
                        fastq_screens=self.fastq_screens,
                        run_multiqc=True,
                        max_jobs=1)
        self.assertEqual(status,0)
        # Check detected outputs
        for p in ("AB",):
            qc_dir = os.path.join(mockdir.dirn,p,"qc")
            qcoutputs = QCOutputs(qc_dir)
            print(qcoutputs.outputs)
            for qc_module in ("fastqc_r1",
                              "fastqc_r3",
                              "screens_r1",
                              "screens_r3",
                              "sequence_lengths",
                              "picard_insert_size_metrics",
                              "rseqc_genebody_coverage",
                              "rseqc_infer_experiment",
                              "qualimap_rnaseq",
                              "cellranger-atac_count",
                              "multiqc"):
                self.assertTrue(qc_module in qcoutputs.outputs,
                                "Project '%s': missing output '%s'" %
                                (p,qc_module))
        # Check output and reports
        for p in ("AB","undetermined"):
            for f in ("qc",
                      "qc_report.html",
                      "qc_report.%s.%s.zip" % (
                          p,
                          '170901_M00879_0087_000000000-AGEW9'),
                      "multiqc_report.html"):
                self.assertTrue(os.path.exists(os.path.join(mockdir.dirn,
                                                            p,f)),
                                "Missing %s in project '%s'" % (f,p))
            # Check zip file has MultiQC report
            zip_file = os.path.join(mockdir.dirn,p,
                                    "qc_report.%s.%s.zip" % (
                                        p,
                                        '170901_M00879_0087_000000000-AGEW9'))
            with zipfile.ZipFile(zip_file) as z:
                multiqc = os.path.join(
                    "qc_report.%s.%s" % (
                        p,'170901_M00879_0087_000000000-AGEW9'),
                    "multiqc_report.html")
                self.assertTrue(multiqc in z.namelist())
        # Check lane splitting for undetermined
        qc_info = AnalysisProjectQCDirInfo(os.path.join(mockdir.dirn,
                                                        "undetermined",
                                                        "qc",
                                                        "qc.info"))
        self.assertTrue(qc_info.fastqs_split_by_lane)

    def test_run_qc_10x_multiome_gex(self):
        """run_qc: 10x Multiome GEX
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        MockStar.create(os.path.join(self.bin,"STAR"))
        MockSamtools.create(os.path.join(self.bin,"samtools"))
        MockPicard.create(os.path.join(self.bin,"picard"))
        MockGtf2bed.create(os.path.join(self.bin,"gtf2bed"))
        MockRSeQC.create(os.path.join(self.bin,"infer_experiment.py"))
        MockRSeQC.create(os.path.join(self.bin,"geneBody_coverage.py"))
        MockQualimap.create(os.path.join(self.bin,"qualimap"))
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"))
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger-arc"))
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            project_metadata={ "CDE": { "Organism": "mouse",
                                        "Single cell platform":
                                        "10xGenomics Single Cell Multiome",
                                        "Library type": "GEX", } },
            reads=('R1','R2','I1'),
            top_dir=self.dirn)
        mockdir.create()
        # Settings file with fastq_strand indexes and
        # polling interval
        settings_ini = os.path.join(self.dirn,"auto_process.ini")
        with open(settings_ini,'w') as s:
            s.write("""[general]
poll_interval = {poll_interval}

[organism:mouse]
star_index = /data/genomeIndexes/mm10/STAR
annotation_bed = {mm10_bed}
annotation_gtf = {mm10_gtf}
cellranger_reference = /data/cellranger/gex_references/mm10-gex
cellranger_arc_reference = /data/cellranger/arc_references/mm10-arc
""".format(mm10_bed=self.ref_data['mm10']['bed'],
           mm10_gtf=self.ref_data['mm10']['gtf'],
           poll_interval=POLL_INTERVAL))
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn,
                         settings=Settings(settings_ini))
        # Run the QC
        status = run_qc(ap,
                        fastq_screens=self.fastq_screens,
                        run_multiqc=True,
                        max_jobs=1)
        self.assertEqual(status,0)
        # Check detected outputs
        for p in ("CDE",):
            qc_dir = os.path.join(mockdir.dirn,p,"qc")
            qcoutputs = QCOutputs(qc_dir)
            print(qcoutputs.outputs)
            for qc_module in ("fastqc_r1",
                              "fastqc_r2",
                              "screens_r2",
                              "sequence_lengths",
                              "rseqc_genebody_coverage",
                              "rseqc_infer_experiment",
                              "qualimap_rnaseq",
                              "cellranger_count",
                              "multiqc"):
                self.assertTrue(qc_module in qcoutputs.outputs,
                                "Project '%s': missing output '%s'" %
                                (p,qc_module))
        # Check output and reports
        for p in ("CDE","undetermined"):
            for f in ("qc",
                      "qc_report.html",
                      "qc_report.%s.%s.zip" % (
                          p,
                          '170901_M00879_0087_000000000-AGEW9'),
                      "multiqc_report.html"):
                self.assertTrue(os.path.exists(os.path.join(mockdir.dirn,
                                                            p,f)),
                                "Missing %s in project '%s'" % (f,p))
            # Check zip file has MultiQC report
            zip_file = os.path.join(mockdir.dirn,p,
                                    "qc_report.%s.%s.zip" % (
                                        p,
                                        '170901_M00879_0087_000000000-AGEW9'))
            with zipfile.ZipFile(zip_file) as z:
                multiqc = os.path.join(
                    "qc_report.%s.%s" % (
                        p,'170901_M00879_0087_000000000-AGEW9'),
                    "multiqc_report.html")
                self.assertTrue(multiqc in z.namelist())
        # Check lane splitting for undetermined
        qc_info = AnalysisProjectQCDirInfo(os.path.join(mockdir.dirn,
                                                        "undetermined",
                                                        "qc",
                                                        "qc.info"))
        self.assertTrue(qc_info.fastqs_split_by_lane)

    def test_run_qc_parse_evercode_scRNAseq(self):
        """run_qc: Parse Evercode scRNA-seq
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        MockStar.create(os.path.join(self.bin,"STAR"))
        MockSamtools.create(os.path.join(self.bin,"samtools"))
        MockGtf2bed.create(os.path.join(self.bin,"gtf2bed"))
        MockRSeQC.create(os.path.join(self.bin,"infer_experiment.py"))
        MockRSeQC.create(os.path.join(self.bin,"geneBody_coverage.py"))
        MockQualimap.create(os.path.join(self.bin,"qualimap"))
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
                                       "Parse Evercode",
                                       "Library type": "scRNA-seq", } },
            reads=('R1','R2'),
            top_dir=self.dirn)
        mockdir.create()
        # Settings file with fastq_strand indexes and
        # polling interval
        settings_ini = os.path.join(self.dirn,"auto_process.ini")
        with open(settings_ini,'w') as s:
            s.write("""[general]
poll_interval = {poll_interval}

[organism:human]
star_index = /data/genomeIndexes/hg38/STAR
annotation_bed = {hg38_bed}
annotation_gtf = {hg38_gtf}
""".format(hg38_bed=self.ref_data['hg38']['bed'],
           hg38_gtf=self.ref_data['hg38']['gtf'],
           poll_interval=POLL_INTERVAL))
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn,
                         settings=Settings(settings_ini))
        # Run the QC
        status = run_qc(ap,
                        fastq_screens=self.fastq_screens,
                        run_multiqc=True,
                        max_jobs=1)
        self.assertEqual(status,0)
        # Check detected outputs
        for p in ("AB",):
            qc_dir = os.path.join(mockdir.dirn,p,"qc")
            qcoutputs = QCOutputs(qc_dir)
            print(qcoutputs.outputs)
            for qc_module in ("fastqc_r1",
                              "fastqc_r2",
                              "screens_r1",
                              "sequence_lengths",
                              "rseqc_genebody_coverage",
                              "rseqc_infer_experiment",
                              "qualimap_rnaseq",
                              "multiqc"):
                self.assertTrue(qc_module in qcoutputs.outputs,
                                "Project '%s': missing output '%s'" %
                                (p,qc_module))
        # Check output and reports
        for p in ("AB","undetermined"):
            for f in ("qc",
                      "qc_report.html",
                      "qc_report.%s.%s.zip" % (
                          p,
                          '170901_M00879_0087_000000000-AGEW9'),
                      "multiqc_report.html"):
                self.assertTrue(os.path.exists(os.path.join(mockdir.dirn,
                                                            p,f)),
                                "Missing %s in project '%s'" % (f,p))
            # Check zip file has MultiQC report
            zip_file = os.path.join(mockdir.dirn,p,
                                    "qc_report.%s.%s.zip" % (
                                        p,
                                        '170901_M00879_0087_000000000-AGEW9'))
            with zipfile.ZipFile(zip_file) as z:
                multiqc = os.path.join(
                    "qc_report.%s.%s" % (
                        p,'170901_M00879_0087_000000000-AGEW9'),
                    "multiqc_report.html")
                self.assertTrue(multiqc in z.namelist())
        # Check lane splitting for undetermined
        qc_info = AnalysisProjectQCDirInfo(os.path.join(mockdir.dirn,
                                                        "undetermined",
                                                        "qc",
                                                        "qc.info"))
        self.assertTrue(qc_info.fastqs_split_by_lane)

    def test_run_qc_fastqs_without_no_lane_splitting(self):
        """run_qc: Fastqs without --no-lane-splitting
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        MockStar.create(os.path.join(self.bin,"STAR"))
        MockSamtools.create(os.path.join(self.bin,"samtools"))
        MockGtf2bed.create(os.path.join(self.bin,"gtf2bed"))
        MockRSeQC.create(os.path.join(self.bin,"infer_experiment.py"))
        MockRSeQC.create(os.path.join(self.bin,"geneBody_coverage.py"))
        MockQualimap.create(os.path.join(self.bin,"qualimap"))
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            no_lane_splitting=False,
            paired_end=False,
            metadata={ "instrument_datestamp": "170901" },
            project_metadata={ "AB": { "Library type": "RNA-seq",
                                       "Organism": "human", },
                               "CDE": { "Library type": "ChIP-seq",
                                        "Organism": "mouse", } },
            top_dir=self.dirn)
        mockdir.create()
        # Settings file with reference data and polling interval
        settings_ini = os.path.join(self.dirn,"auto_process.ini")
        with open(settings_ini,'w') as s:
            s.write("""[general]
poll_interval = {poll_interval}

[organism:human]
star_index = /data/genomeIndexes/hg38/STAR
annotation_bed = {hg38_bed}
annotation_gtf = {hg38_gtf}

[organism:mouse]
star_index = /data/genomeIndexes/mm10/STAR
annotation_bed = {mm10_bed}
annotation_gtf = {mm10_gtf}
""".format(hg38_bed=self.ref_data['hg38']['bed'],
           hg38_gtf=self.ref_data['hg38']['gtf'],
           mm10_bed=self.ref_data['mm10']['bed'],
           mm10_gtf=self.ref_data['mm10']['gtf'],
           poll_interval=POLL_INTERVAL))
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn,
                         settings=Settings(settings_ini))
        # Run the QC
        status = run_qc(ap,
                        fastq_screens=self.fastq_screens,
                        run_multiqc=True,
                        max_jobs=1)
        self.assertEqual(status,0)
        # Check detected outputs
        for p in ("AB","CDE"):
            qc_dir = os.path.join(mockdir.dirn,p,"qc")
            qcoutputs = QCOutputs(qc_dir)
            for qc_module in ("fastqc_r1",
                              "screens_r1",
                              "sequence_lengths",
                              "rseqc_genebody_coverage",
                              "rseqc_infer_experiment",
                              "qualimap_rnaseq",
                              "multiqc"):
                self.assertTrue(qc_module in qcoutputs.outputs,
                                "Project '%s': missing output '%s'" %
                                (p,qc_module))
        # Check output and reports
        for p in ("AB","CDE","undetermined"):
            for f in ("qc",
                      "qc_report.html",
                      "qc_report.%s.%s.zip" % (
                          p,
                          '170901_M00879_0087_000000000-AGEW9'),
                      "multiqc_report.html"):
                self.assertTrue(os.path.exists(os.path.join(mockdir.dirn,
                                                            p,f)),
                                "Missing %s in project '%s'" % (f,p))
            # Check zip file has MultiQC report
            zip_file = os.path.join(mockdir.dirn,p,
                                    "qc_report.%s.%s.zip" % (
                                        p,
                                        '170901_M00879_0087_000000000-AGEW9'))
            with zipfile.ZipFile(zip_file) as z:
                multiqc = os.path.join(
                    "qc_report.%s.%s" % (
                        p,'170901_M00879_0087_000000000-AGEW9'),
                    "multiqc_report.html")
                self.assertTrue(multiqc in z.namelist())
        # Check lane splitting for undetermined
        qc_info = AnalysisProjectQCDirInfo(os.path.join(mockdir.dirn,
                                                        "undetermined",
                                                        "qc",
                                                        "qc.info"))
        self.assertFalse(qc_info.fastqs_split_by_lane)

    def test_run_qc_specify_qc_protocol(self):
        """run_qc: specify QC protocol for a project
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        MockStar.create(os.path.join(self.bin,"STAR"))
        MockSamtools.create(os.path.join(self.bin,"samtools"))
        MockPicard.create(os.path.join(self.bin,"picard"))
        MockGtf2bed.create(os.path.join(self.bin,"gtf2bed"))
        MockRSeQC.create(os.path.join(self.bin,"infer_experiment.py"))
        MockRSeQC.create(os.path.join(self.bin,"geneBody_coverage.py"))
        MockQualimap.create(os.path.join(self.bin,"qualimap"))
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901", },
            project_metadata={ "AB": { "Library type": "RNA-seq",
                                       "Organism": "human", },
                               "CDE": { "Library type": "ChIP-seq",
                                        "Organism": "mouse", } },
            top_dir=self.dirn)
        mockdir.create()
        # Settings file with reference data and polling interval
        settings_ini = os.path.join(self.dirn,"auto_process.ini")
        with open(settings_ini,'w') as s:
            s.write("""[general]
poll_interval = {poll_interval}

[organism:human]
star_index = /data/hg38/star_index
annotation_bed = {hg38_bed}
annotation_gtf = {hg38_gtf}

[organism:mouse]
star_index = /data/mm10/star_index
annotation_bed = {mm10_bed}
annotation_gtf = {mm10_gtf}
""".format(hg38_bed=self.ref_data['hg38']['bed'],
           hg38_gtf=self.ref_data['hg38']['gtf'],
           mm10_bed=self.ref_data['mm10']['bed'],
           mm10_gtf=self.ref_data['mm10']['gtf'],
           poll_interval=POLL_INTERVAL))
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn,
                         settings=Settings(settings_ini))
        # Run the QC
        status = run_qc(ap,
                        protocols={ "CDE": "minimal" },
                        fastq_screens=self.fastq_screens,
                        run_multiqc=True,
                        max_jobs=1)
        self.assertEqual(status,0)
        # Check output and reports
        for p in ("AB","CDE","undetermined"):
            for f in ("qc",
                      "qc_report.html",
                      "qc_report.%s.%s.zip" % (
                          p,
                          '170901_M00879_0087_000000000-AGEW9'),
                      "multiqc_report.html"):
                self.assertTrue(os.path.exists(os.path.join(mockdir.dirn,
                                                            p,f)),
                                "Missing %s in project '%s'" % (f,p))
            # Check zip file has MultiQC report
            zip_file = os.path.join(mockdir.dirn,p,
                                    "qc_report.%s.%s.zip" % (
                                        p,
                                        '170901_M00879_0087_000000000-AGEW9'))
            with zipfile.ZipFile(zip_file) as z:
                multiqc = os.path.join(
                    "qc_report.%s.%s" % (
                        p,'170901_M00879_0087_000000000-AGEW9'),
                    "multiqc_report.html")
                self.assertTrue(multiqc in z.namelist())
        # Check QC protocol for "CDE"
        qc_info = AnalysisProjectQCDirInfo(os.path.join(mockdir.dirn,
                                                        "CDE",
                                                        "qc",
                                                        "qc.info"))
        self.assertTrue(qc_info.protocol_specification.startswith("minimal:"))
