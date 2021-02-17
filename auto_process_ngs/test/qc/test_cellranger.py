#######################################################################
# Unit tests for qc/cellranger.py
#######################################################################

import unittest
import os
import shutil
import tempfile
from auto_process_ngs.mock import MockAnalysisProject
from auto_process_ngs.mock import UpdateAnalysisProject
from auto_process_ngs.analysis import AnalysisProject

from auto_process_ngs.qc.cellranger import CellrangerCount

class TestCellrangerCount(unittest.TestCase):
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestCellrangerCount')
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           "10xGenomics Chromium 3'v3" })
        p.create(top_dir=self.dirn)
        self.project = AnalysisProject("PJB",os.path.join(self.dirn,"PJB"))

    def tearDown(self):
        # Remove the temporary test directory
        shutil.rmtree(self.dirn)

    def test_cellrangercount_501(self):
        """
        CellrangerCount: check outputs from cellranger count (v5.0.1)
        """
        # Add cellranger count outputs
        UpdateAnalysisProject(self.project).add_cellranger_count_outputs()
        # Do tests
        count_dir = os.path.join(self.project.qc_dir,"cellranger_count","PJB1")
        cmdline = "/path/to/cellranger count --id PJB1 --fastqs /path/to/PJB/fastqs --sample PJB1 --transcriptome /data/refdata-gex-GRCh38-2020-A --chemistry auto --r1-length=26 --jobmode=local --localcores=16 --localmem=48 --maxjobs=1 --jobinterval=100"
        with open(os.path.join(count_dir,"_cmdline"),'wt') as fp:
            fp.write("%s\n" % cmdline)
        cellranger_count = CellrangerCount(count_dir)
        self.assertEqual(cellranger_count.dir,count_dir)
        self.assertEqual(cellranger_count.sample_name,"PJB1")
        self.assertEqual(cellranger_count.metrics_csv,
                         os.path.join(count_dir,"outs","metrics_summary.csv"))
        self.assertEqual(cellranger_count.web_summary,
                         os.path.join(count_dir,"outs","web_summary.html"))
        self.assertEqual(cellranger_count.cmdline_file,
                         os.path.join(count_dir,"_cmdline"))
        self.assertEqual(cellranger_count.cmdline,cmdline)
        self.assertEqual(cellranger_count.version,None)
        self.assertEqual(cellranger_count.reference_data,
                         "/data/refdata-gex-GRCh38-2020-A")
        self.assertEqual(cellranger_count.cellranger_exe,
                         "/path/to/cellranger")
        self.assertEqual(cellranger_count.pipeline_name,"cellranger")
    def test_cellrangercount_cellranger_310(self):
        """
        CellrangerCount: check outputs from cellranger count (v3.1.0)
        """
        # Add cellranger count outputs
        UpdateAnalysisProject(self.project).add_cellranger_count_outputs()
        # Do tests
        count_dir = os.path.join(self.project.qc_dir,"cellranger_count","PJB1")
        cmdline = "/path/to/cellranger-cs/3.1.0/bin/count --id PJB1 --fastqs /path/to/PJB/fastqs --sample PJB1 --transcriptome /data/refdata-cellranger-GRCh38-1.2.0 --chemistry auto --jobmode=local --localcores=16 --localmem=48 --maxjobs=1 --jobinterval=100"
        with open(os.path.join(count_dir,"_cmdline"),'wt') as fp:
            fp.write("%s\n" % cmdline)
        cellranger_count = CellrangerCount(count_dir)
        self.assertEqual(cellranger_count.dir,count_dir)
        self.assertEqual(cellranger_count.sample_name,"PJB1")
        self.assertEqual(cellranger_count.metrics_csv,
                         os.path.join(count_dir,"outs","metrics_summary.csv"))
        self.assertEqual(cellranger_count.web_summary,
                         os.path.join(count_dir,"outs","web_summary.html"))
        self.assertEqual(cellranger_count.cmdline_file,
                         os.path.join(count_dir,"_cmdline"))
        self.assertEqual(cellranger_count.cmdline,cmdline)
        self.assertEqual(cellranger_count.version,None)
        self.assertEqual(cellranger_count.reference_data,
                         "/data/refdata-cellranger-GRCh38-1.2.0")
        self.assertEqual(cellranger_count.cellranger_exe,
                         "/path/to/cellranger-cs/3.1.0/bin/count")
        self.assertEqual(cellranger_count.pipeline_name,"cellranger")

    def test_cellrangercount_cellranger_atac_120(self):
        """
        CellrangerCount: check outputs from cellranger-atac count (v1.2.0)
        """
        # Add cellranger count outputs
        UpdateAnalysisProject(self.project).add_cellranger_count_outputs(
            cellranger='cellranger-atac')
        # Do tests
        count_dir = os.path.join(self.project.qc_dir,"cellranger_count","PJB1")
        cmdline = "/path/to/cellranger-atac-cs/1.2.0/bin/count --id PJB1 --fastqs /path/to/PJB/fastqs --sample PJB1 --reference /data/refdata-cellranger-atac-GRCh38-1.2.0 --jobmode=local --localcores=16 --localmem=128 --maxjobs=48 --jobinterval=100"
        with open(os.path.join(count_dir,"_cmdline"),'wt') as fp:
            fp.write("%s\n" % cmdline)
        cellranger_count = CellrangerCount(count_dir)
        self.assertEqual(cellranger_count.dir,count_dir)
        self.assertEqual(cellranger_count.sample_name,"PJB1")
        self.assertEqual(cellranger_count.metrics_csv,
                         os.path.join(count_dir,"outs","summary.csv"))
        self.assertEqual(cellranger_count.web_summary,
                         os.path.join(count_dir,"outs","web_summary.html"))
        self.assertEqual(cellranger_count.cmdline_file,
                         os.path.join(count_dir,"_cmdline"))
        self.assertEqual(cellranger_count.cmdline,cmdline)
        self.assertEqual(cellranger_count.version,None)
        self.assertEqual(cellranger_count.reference_data,
                         "/data/refdata-cellranger-atac-GRCh38-1.2.0")
        self.assertEqual(cellranger_count.cellranger_exe,
                         "/path/to/cellranger-atac-cs/1.2.0/bin/count")
        self.assertEqual(cellranger_count.pipeline_name,"cellranger-atac")

    def test_cellrangercount_cellranger_arc_120(self):
        """
        CellrangerCount: check outputs from cellranger-arc count (v1.0.0)
        """
        # Add cellranger count outputs
        UpdateAnalysisProject(self.project).add_cellranger_count_outputs(
            cellranger='cellranger-atac')
        # Do tests
        count_dir = os.path.join(self.project.qc_dir,"cellranger_count","PJB1")
        cmdline = "/path/to/cellranger-arc count --id PJB1 --fastqs /path/to/PJB/fastqs --sample PJB1 --reference /data/refdata-cellranger-arc-GRCh38-2020-A --libraries /path/to/libraries.csv --jobmode=local --localcores=16 --localmem=128 --maxjobs=48 --jobinterval=100"
        with open(os.path.join(count_dir,"_cmdline"),'wt') as fp:
            fp.write("%s\n" % cmdline)
        cellranger_count = CellrangerCount(count_dir)
        self.assertEqual(cellranger_count.dir,count_dir)
        self.assertEqual(cellranger_count.sample_name,"PJB1")
        self.assertEqual(cellranger_count.metrics_csv,
                         os.path.join(count_dir,"outs","summary.csv"))
        self.assertEqual(cellranger_count.web_summary,
                         os.path.join(count_dir,"outs","web_summary.html"))
        self.assertEqual(cellranger_count.cmdline_file,
                         os.path.join(count_dir,"_cmdline"))
        self.assertEqual(cellranger_count.cmdline,cmdline)
        self.assertEqual(cellranger_count.version,None)
        self.assertEqual(cellranger_count.reference_data,
                         "/data/refdata-cellranger-arc-GRCh38-2020-A")
        self.assertEqual(cellranger_count.cellranger_exe,
                         "/path/to/cellranger-arc")
        self.assertEqual(cellranger_count.pipeline_name,"cellranger-arc")

    def test_cellrangercount_with_data(self):
        """
        CellrangerCount: check outputs when data are supplied
        """
        # Add cellranger count outputs
        UpdateAnalysisProject(self.project).add_cellranger_count_outputs()
        # Do tests
        count_dir = os.path.join(self.project.qc_dir,"cellranger_count","PJB1")
        cmdline = "/path/to/cellranger count --id PJB1 --fastqs /path/to/PJB/fastqs --sample PJB1 --transcriptome /data/refdata-gex-GRCh38-2020-A --chemistry auto --r1-length=26 --jobmode=local --localcores=16 --localmem=48 --maxjobs=1 --jobinterval=100"
        with open(os.path.join(count_dir,"_cmdline"),'wt') as fp:
            fp.write("%s\n" % cmdline)
        cellranger_count = CellrangerCount(
            count_dir,
            cellranger_exe="/alt/path/to/cellranger",
            version="5.0.1",
            reference_data="/alt/data/refdata-gex-GRCh38-2020-A")
        self.assertEqual(cellranger_count.dir,count_dir)
        self.assertEqual(cellranger_count.sample_name,"PJB1")
        self.assertEqual(cellranger_count.metrics_csv,
                         os.path.join(count_dir,"outs","metrics_summary.csv"))
        self.assertEqual(cellranger_count.web_summary,
                         os.path.join(count_dir,"outs","web_summary.html"))
        self.assertEqual(cellranger_count.cmdline_file,
                         os.path.join(count_dir,"_cmdline"))
        self.assertEqual(cellranger_count.cmdline,cmdline)
        self.assertEqual(cellranger_count.version,"5.0.1")
        self.assertEqual(cellranger_count.reference_data,
                         "/alt/data/refdata-gex-GRCh38-2020-A")
        self.assertEqual(cellranger_count.cellranger_exe,
                         "/alt/path/to/cellranger")
        self.assertEqual(cellranger_count.pipeline_name,"cellranger")

    def test_cellrangercount_missing_directory(self):
        """
        CellrangerCount: handle missing directory
        """
        # Do tests
        count_dir = os.path.join(self.project.qc_dir,"cellranger_count","PJB1")
        cellranger_count = CellrangerCount(count_dir)
        self.assertRaises(OSError,
                          getattr,cellranger_count,'dir')
        self.assertEqual(cellranger_count.sample_name,None)
        self.assertRaises(OSError,
                          getattr,cellranger_count,'metrics_csv')
        self.assertRaises(OSError,
                          getattr,cellranger_count,'web_summary')
        self.assertEqual(cellranger_count.cmdline_file,None)
        self.assertEqual(cellranger_count.cmdline,None)
        self.assertEqual(cellranger_count.version,None)
        self.assertEqual(cellranger_count.reference_data,None)
        self.assertEqual(cellranger_count.cellranger_exe,None)
        self.assertEqual(cellranger_count.pipeline_name,None)
