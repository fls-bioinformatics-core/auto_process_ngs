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
from auto_process_ngs.tenx.cellplex import CellrangerMultiConfigCsv
from auto_process_ngs.tenx.metrics import MultiplexSummary

from auto_process_ngs.qc.cellranger import CellrangerCount
from auto_process_ngs.qc.cellranger import CellrangerMulti
from auto_process_ngs.qc.cellranger import cellranger_count_output
from auto_process_ngs.qc.cellranger import cellranger_atac_count_output
from auto_process_ngs.qc.cellranger import cellranger_arc_count_output
from auto_process_ngs.qc.cellranger import cellranger_multi_output

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

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
        self.assertEqual(cellranger_count.mode,"count")
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
        self.assertEqual(cellranger_count.mode,"count")
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
        self.assertEqual(cellranger_count.mode,"count")
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
        self.assertEqual(cellranger_count.mode,"count")
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
        self.assertEqual(cellranger_count.mode,"count")
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
        self.assertEqual(cellranger_count.mode,"count")
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

class TestCellrangerMulti(unittest.TestCase):
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestCellrangerMulti')
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_GEX_S1_R1_001.fastq.gz",
                                       "PJB1_GEX_S1_R2_001.fastq.gz",
                                       "PJB2_MC_S2_R1_001.fastq.gz",
                                       "PJB2_MC_S2_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           "10xGenomics Chromium 3'v3" })
        p.create(top_dir=self.dirn)
        self.project = AnalysisProject("PJB",os.path.join(self.dirn,"PJB"))

    def tearDown(self):
        # Remove the temporary test directory
        shutil.rmtree(self.dirn)

    def test_cellrangermulti_cellplex(self):
        """
        CellrangerMulti: check outputs from cellranger multi for CellPlex
        """
        # Add config.csv file
        config_csv = os.path.join(self.project.dirn,
                                  "10x_multi_config.csv")
        with open(config_csv,'wt') as fp:
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB1_GEX,/data/runs/fastqs_gex,any,PJB1,gene expression,
PJB2_MC,/data/runs/fastqs_mc,any,PJB2,Multiplexing Capture,

[samples]
sample_id,cmo_ids,description
PBA,CMO301,PBA
PBB,CMO302,PBB
""")
        # Add cellranger multi outputs
        UpdateAnalysisProject(self.project).add_cellranger_multi_outputs(
            config_csv)
        # Do tests
        multi_dir = os.path.join(self.project.qc_dir,"cellranger_multi")
        cmdline = "/path/to/cellranger multi --id PJB --csv %s --jobmode=local --localcores=16 --localmem=48 --maxjobs=1 --jobinterval=100" % config_csv
        with open(os.path.join(multi_dir,"_cmdline"),'wt') as fp:
            fp.write("%s\n" % cmdline)
        cellranger_multi = CellrangerMulti(multi_dir)
        self.assertEqual(cellranger_multi.mode,"multi")
        self.assertEqual(cellranger_multi.dir,multi_dir)
        self.assertEqual(cellranger_multi.sample_names,["PBA","PBB"])
        self.assertEqual(cellranger_multi.metrics_csv('PBA'),
                         os.path.join(multi_dir,
                                      "outs",
                                      "per_sample_outs",
                                      "PBA",
                                      "metrics_summary.csv"))
        self.assertEqual(cellranger_multi.metrics_csv('PBB'),
                         os.path.join(multi_dir,
                                      "outs",
                                      "per_sample_outs",
                                      "PBB",
                                      "metrics_summary.csv"))
        self.assertTrue(isinstance(cellranger_multi.metrics('PBA'),
                                   MultiplexSummary))
        self.assertTrue(isinstance(cellranger_multi.metrics('PBB'),
                                   MultiplexSummary))
        self.assertEqual(cellranger_multi.web_summary('PBA'),
                         os.path.join(multi_dir,
                                      "outs",
                                      "per_sample_outs",
                                      "PBA",
                                      "web_summary.html"))
        self.assertEqual(cellranger_multi.web_summary('PBB'),
                         os.path.join(multi_dir,
                                      "outs",
                                      "per_sample_outs",
                                      "PBB",
                                      "web_summary.html"))
        self.assertEqual(cellranger_multi.cmdline_file,
                         os.path.join(multi_dir,"_cmdline"))
        self.assertEqual(cellranger_multi.cmdline,cmdline)
        self.assertEqual(cellranger_multi.version,None)
        self.assertEqual(cellranger_multi.reference_data,
                         "/data/refdata-cellranger-gex-GRCh38-2020-A")
        self.assertEqual(cellranger_multi.probe_set,None)
        self.assertEqual(cellranger_multi.cellranger_exe,
                         "/path/to/cellranger")
        self.assertEqual(cellranger_multi.pipeline_name,"cellranger")

    def test_cellrangermulti_flex(self):
        """
        CellrangerMulti: check outputs from cellranger multi for Flex
        """
        # Add config.csv file
        config_csv = os.path.join(self.project.dirn,
                                  "10x_multi_config.csv")
        with open(config_csv,'wt') as fp:
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A
probe-set,/data/probe_set_v1.0_GRCh38-2020-A.csv
no-bam,true

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB1_flex,/data/runs/fastqs_flex,any,PJB1,Gene Expression,

[samples]
sample_id,probe_barcode_ids,description
PBA,BC001,PBA
PBB,BC002,PBB
""")
        # Add cellranger multi outputs
        UpdateAnalysisProject(self.project).add_cellranger_multi_outputs(
            config_csv)
        # Do tests
        multi_dir = os.path.join(self.project.qc_dir,"cellranger_multi")
        cmdline = "/path/to/cellranger multi --id PJB --csv %s --jobmode=local --localcores=16 --localmem=48 --maxjobs=1 --jobinterval=100" % config_csv
        with open(os.path.join(multi_dir,"_cmdline"),'wt') as fp:
            fp.write("%s\n" % cmdline)
        cellranger_multi = CellrangerMulti(multi_dir)
        self.assertEqual(cellranger_multi.mode,"multi")
        self.assertEqual(cellranger_multi.dir,multi_dir)
        self.assertEqual(cellranger_multi.sample_names,["PBA","PBB"])
        self.assertEqual(cellranger_multi.metrics_csv('PBA'),
                         os.path.join(multi_dir,
                                      "outs",
                                      "per_sample_outs",
                                      "PBA",
                                      "metrics_summary.csv"))
        self.assertEqual(cellranger_multi.metrics_csv('PBB'),
                         os.path.join(multi_dir,
                                      "outs",
                                      "per_sample_outs",
                                      "PBB",
                                      "metrics_summary.csv"))
        self.assertTrue(isinstance(cellranger_multi.metrics('PBA'),
                                   MultiplexSummary))
        self.assertTrue(isinstance(cellranger_multi.metrics('PBB'),
                                   MultiplexSummary))
        self.assertEqual(cellranger_multi.web_summary('PBA'),
                         os.path.join(multi_dir,
                                      "outs",
                                      "per_sample_outs",
                                      "PBA",
                                      "web_summary.html"))
        self.assertEqual(cellranger_multi.web_summary('PBB'),
                         os.path.join(multi_dir,
                                      "outs",
                                      "per_sample_outs",
                                      "PBB",
                                      "web_summary.html"))
        self.assertEqual(cellranger_multi.cmdline_file,
                         os.path.join(multi_dir,"_cmdline"))
        self.assertEqual(cellranger_multi.cmdline,cmdline)
        self.assertEqual(cellranger_multi.version,None)
        self.assertEqual(cellranger_multi.reference_data,
                         "/data/refdata-cellranger-gex-GRCh38-2020-A")
        self.assertEqual(cellranger_multi.probe_set,
                         "/data/probe_set_v1.0_GRCh38-2020-A.csv")
        self.assertEqual(cellranger_multi.cellranger_exe,
                         "/path/to/cellranger")
        self.assertEqual(cellranger_multi.pipeline_name,"cellranger")

class TestCellrangerCountOutputFunction(unittest.TestCase):
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestCellrangerCountOutput')
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)

    def tearDown(self):
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    def test_cellranger_count_output(self):
        """cellranger_count_output: check for project
        """
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        self.assertEqual(cellranger_count_output(project),
                         ('cellranger_count/PJB1/outs/metrics_summary.csv',
                          'cellranger_count/PJB1/outs/web_summary.html',
                          'cellranger_count/PJB2/outs/metrics_summary.csv',
                          'cellranger_count/PJB2/outs/web_summary.html'))

    def test_cellranger_count_output_with_sample(self):
        """cellranger_count_output: check for project and sample
        """
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        self.assertEqual(cellranger_count_output(project,
                                                 sample_name="PJB2"),
                         ('cellranger_count/PJB2/outs/metrics_summary.csv',
                          'cellranger_count/PJB2/outs/web_summary.html'))

    def test_cellranger_count_output_with_prefix(self):
        """cellranger_count_output: check for project and prefix
        """
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        prefix = "cellranger_count/5.0.1/refdata-gex-mm10-2020-A"
        self.assertEqual(cellranger_count_output(project,
                                                 prefix=prefix),
                         ('%s/PJB1/outs/metrics_summary.csv' % prefix,
                          '%s/PJB1/outs/web_summary.html' % prefix,
                          '%s/PJB2/outs/metrics_summary.csv' % prefix,
                          '%s/PJB2/outs/web_summary.html' % prefix))

class TestCellrangerAtacCountOutputFunction(unittest.TestCase):
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestCellrangerAtacCountOutput')
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB1_S1_R3_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz",
                                       "PJB2_S2_R3_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)

    def tearDown(self):
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    def test_cellranger_atac_count_output(self):
        """cellranger_atac_count_output: check for project
        """
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        self.assertEqual(cellranger_atac_count_output(project),
                         ('cellranger_count/PJB1/outs/summary.csv',
                          'cellranger_count/PJB1/outs/web_summary.html',
                          'cellranger_count/PJB2/outs/summary.csv',
                          'cellranger_count/PJB2/outs/web_summary.html'))

    def test_cellranger_atac_count_output_with_sample(self):
        """cellranger_atac_count_output: check for project and sample
        """
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        self.assertEqual(cellranger_atac_count_output(project,
                                                 sample_name="PJB2"),
                         ('cellranger_count/PJB2/outs/summary.csv',
                          'cellranger_count/PJB2/outs/web_summary.html'))

    def test_cellranger_atac_count_output_with_prefix(self):
        """cellranger_atac_count_output: check for project and prefix
        """
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        prefix = "cellranger_count/1.2.0/refdata-cellranger-atac-GRCh38-1.2.0"
        self.assertEqual(cellranger_atac_count_output(project,
                                                      prefix=prefix),
                         ('%s/PJB1/outs/summary.csv' % prefix,
                          '%s/PJB1/outs/web_summary.html' % prefix,
                          '%s/PJB2/outs/summary.csv' % prefix,
                          '%s/PJB2/outs/web_summary.html' % prefix))

class TestCellrangerArcCountOutputFunction(unittest.TestCase):
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestCellrangerArcCountOutput')
        # Make mock analysis project
        p = MockAnalysisProject("PJB_ARC",("PJB1_S1_R1_001.fastq.gz",
                                           "PJB1_S1_R2_001.fastq.gz",
                                           "PJB1_S1_R3_001.fastq.gz",
                                           "PJB2_S2_R1_001.fastq.gz",
                                           "PJB2_S2_R2_001.fastq.gz",
                                           "PJB2_S2_R3_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)

    def test_cellranger_arc_count_output(self):
        """cellranger_arc_count_output: check for project
        """
        project = AnalysisProject(os.path.join(self.wd,"PJB_ARC"))
        self.assertEqual(cellranger_arc_count_output(project),
                         ('cellranger_count/PJB1/outs/summary.csv',
                          'cellranger_count/PJB1/outs/web_summary.html',
                          'cellranger_count/PJB2/outs/summary.csv',
                          'cellranger_count/PJB2/outs/web_summary.html'))

    def test_cellranger_arc_count_output_with_sample(self):
        """cellranger_arc_count_output: check for project and sample
        """
        project = AnalysisProject(os.path.join(self.wd,"PJB_ARC"))
        self.assertEqual(cellranger_arc_count_output(project,
                                                     sample_name="PJB2"),
                         ('cellranger_count/PJB2/outs/summary.csv',
                          'cellranger_count/PJB2/outs/web_summary.html'))

    def test_cellranger_arc_count_output_with_prefix(self):
        """cellranger_arc_count_output: check for project and prefix
        """
        project = AnalysisProject(os.path.join(self.wd,"PJB_ARC"))
        prefix = "cellranger_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A"
        self.assertEqual(cellranger_arc_count_output(project,
                                                     prefix=prefix),
                         ('%s/PJB1/outs/summary.csv' % prefix,
                          '%s/PJB1/outs/web_summary.html' % prefix,
                          '%s/PJB2/outs/summary.csv' % prefix,
                          '%s/PJB2/outs/web_summary.html' % prefix))

class TestCellrangerMultiOutputFunction(unittest.TestCase):
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestCellrangerMultiOutput')
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_GEX_S1_R1_001.fastq.gz",
                                       "PJB1_GEX_S1_R2_001.fastq.gz",
                                       "PJB2_MC_S2_R1_001.fastq.gz",
                                       "PJB2_MC_S2_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Add 10x_multi_config.csv file
        fastq_dir = os.path.join(self.wd,"PJB","fastqs")
        with open(os.path.join(self.wd,
                               "PJB",
                               "10x_multi_config.csv"),'wt') as fp:
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB1_GEX,%s,any,PJB1,gene expression,
PJB2_MC,%s,any,PJB2,Multiplexing Capture,

[samples]
sample_id,cmo_ids,description
PBA,CMO301,PBA
PBB,CMO302,PBB
""" % (fastq_dir,fastq_dir))

    def tearDown(self):
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    def test_cellranger_multi_output(self):
        """cellranger_multi_output: check for project
        """
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        config_csv = os.path.join(self.wd,"PJB","10x_multi_config.csv")
        self.assertEqual(cellranger_multi_output(project,
                                                 config_csv),
                         ('cellranger_multi/outs/per_sample_outs/PBA/metrics_summary.csv',
                          'cellranger_multi/outs/per_sample_outs/PBA/web_summary.html',
                          'cellranger_multi/outs/per_sample_outs/PBB/metrics_summary.csv',
                          'cellranger_multi/outs/per_sample_outs/PBB/web_summary.html',
                          'cellranger_multi/outs/multi/multiplexing_analysis/tag_calls_summary.csv',))

    def test_cellranger_multi_output_with_sample(self):
        """cellranger_multi_output: check for project and sample
        """
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        config_csv = os.path.join(self.wd,"PJB","10x_multi_config.csv")
        self.assertEqual(cellranger_multi_output(project,
                                                 config_csv,
                                                 sample_name="PBB"),
                         ('cellranger_multi/outs/per_sample_outs/PBB/metrics_summary.csv',
                          'cellranger_multi/outs/per_sample_outs/PBB/web_summary.html',
                          'cellranger_multi/outs/multi/multiplexing_analysis/tag_calls_summary.csv',))

    def test_cellranger_multi_output_with_prefix(self):
        """cellranger_multi_output: check for project and prefix
        """
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        config_csv = os.path.join(self.wd,"PJB","10x_multi_config.csv")
        prefix = "cellranger_multi/6.0.0/refdata-gex-mm10-2020-A"
        self.assertEqual(cellranger_multi_output(project,
                                                 config_csv,
                                                 prefix=prefix),
                         ('%s/outs/per_sample_outs/PBA/metrics_summary.csv' % prefix,
                          '%s/outs/per_sample_outs/PBA/web_summary.html' % prefix,
                          '%s/outs/per_sample_outs/PBB/metrics_summary.csv' % prefix,
                          '%s/outs/per_sample_outs/PBB/web_summary.html' % prefix,
                          '%s/outs/multi/multiplexing_analysis/tag_calls_summary.csv' % prefix,))

    def test_cellranger_multi_output_no_config_csv(self):
        """cellranger_multi_output: missing config.csv file
        """
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        config_csv = os.path.join(self.wd,"PJB","10x_multi_config.csv.missing")
        self.assertEqual(cellranger_multi_output(project,config_csv),[])
