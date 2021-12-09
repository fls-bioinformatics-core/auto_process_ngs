#######################################################################
# Unit tests for qc/outputs.py
#######################################################################

import unittest
import tempfile
import shutil
import os
from auto_process_ngs.mock import MockIlluminaQcSh
from auto_process_ngs.mock import MockAnalysisProject
from auto_process_ngs.mock import UpdateAnalysisProject
from auto_process_ngs.analysis import AnalysisProject
from auto_process_ngs.qc.outputs import fastq_screen_output
from auto_process_ngs.qc.outputs import fastqc_output
from auto_process_ngs.qc.outputs import fastq_strand_output
from auto_process_ngs.qc.outputs import cellranger_count_output
from auto_process_ngs.qc.outputs import cellranger_atac_count_output
from auto_process_ngs.qc.outputs import cellranger_arc_count_output
from auto_process_ngs.qc.outputs import cellranger_multi_output
from auto_process_ngs.qc.outputs import check_illumina_qc_outputs
from auto_process_ngs.qc.outputs import check_fastq_screen_outputs
from auto_process_ngs.qc.outputs import check_fastqc_outputs
from auto_process_ngs.qc.outputs import check_fastq_strand_outputs
from auto_process_ngs.qc.outputs import check_cellranger_count_outputs
from auto_process_ngs.qc.outputs import check_cellranger_atac_count_outputs
from auto_process_ngs.qc.outputs import check_cellranger_arc_count_outputs
from auto_process_ngs.qc.outputs import expected_outputs

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

class TestFastqScreenOutputFunction(unittest.TestCase):
    def test_fastq_screen_output(self):
        """fastq_screen_output: handles .fastq file
        """
        self.assertEqual(fastq_screen_output('/data/PB/PB1_ATTAGG_L001_R1_001.fastq',
                                             'model_organisms'),
                         ('PB1_ATTAGG_L001_R1_001_model_organisms_screen.png',
                          'PB1_ATTAGG_L001_R1_001_model_organisms_screen.txt'))
    def test_fastq_screen_output_fastqgz(self):
        """fastq_screen_output: handles fastq.gz file
        """
        self.assertEqual(fastq_screen_output('/data/PB/PB1_ATTAGG_L001_R1_001.fastq.gz',
                                             'model_organisms'),
                         ('PB1_ATTAGG_L001_R1_001_model_organisms_screen.png',
                          'PB1_ATTAGG_L001_R1_001_model_organisms_screen.txt'))

class TestFastqcOutputFunction(unittest.TestCase):
    def test_fastqc_output(self):
        """fastqc_output: handles .fastq file
        """
        self.assertEqual(fastqc_output('/data/PB/PB1_ATTAGG_L001_R1_001.fastq'),
                         ('PB1_ATTAGG_L001_R1_001_fastqc',
                          'PB1_ATTAGG_L001_R1_001_fastqc.html',
                          'PB1_ATTAGG_L001_R1_001_fastqc.zip'))
    def test_fastqc_output_fastqgz(self):
        """fastqc_output: handles fastq.gz file
        """
        self.assertEqual(fastqc_output('/data/PB/PB1_ATTAGG_L001_R1_001.fastq.gz'),
                         ('PB1_ATTAGG_L001_R1_001_fastqc',
                          'PB1_ATTAGG_L001_R1_001_fastqc.html',
                          'PB1_ATTAGG_L001_R1_001_fastqc.zip'))

class TestFastqStrandOutputFunction(unittest.TestCase):
    def test_fastq_strand_output(self):
        """fastq_strand_output: handles .fastq file
        """
        self.assertEqual(fastq_strand_output(
            '/data/PB/PB1_ATTAGG_L001_R1_001.fastq'),
                         'PB1_ATTAGG_L001_R1_001_fastq_strand.txt')
    def test_fastq_strand_output_fastqgz(self):
        """fastq_strand_output: handles fastq.gz file
        """
        self.assertEqual(fastq_strand_output(
            '/data/PB/PB1_ATTAGG_L001_R1_001.fastq.gz'),
                         'PB1_ATTAGG_L001_R1_001_fastq_strand.txt')

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

class TestCheckIlluminaQcOutputs(unittest.TestCase):
    """
    Tests for the 'check_illumina_qc_outputs' function
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestIlluminaQcOutputs')

    def tearDown(self):
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    def test_check_illumina_qc_outputs_standardPE_all_missing(self):
        """
        check_illumina_qc_outputs: all illumina_qc.sh outputs missing (standardPE)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Get the outputs
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        # Check
        self.assertEqual(check_illumina_qc_outputs(project,
                                                   qc_dir="qc",
                                                   qc_protocol="standardPE"),
                         project.fastqs)

    def test_check_illumina_qc_outputs_standardPE_all_present(self):
        """
        check_illumina_qc_outputs: all illumina_qc.sh outputs present (standardPE)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Add QC artefacts
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            include_fastq_strand=False,
            include_multiqc=False)
        # Check
        self.assertEqual(check_illumina_qc_outputs(project,
                                                   qc_dir="qc",
                                                   qc_protocol="standardPE"),
                         [])

    def test_check_illumina_qc_outputs_standardPE_some_missing(self):
        """
        check_illumina_qc_outputs: some illumina_qc.sh outputs missing (standardPE)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Add QC artefacts
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            include_fastq_strand=False,
            include_multiqc=False)
        # Remove some outputs
        for f in ("PJB1_S1_R1_001_fastqc.html",
                  "PJB1_S1_R1_001_model_organisms_screen.txt",):
            os.remove(os.path.join(project.qc_dir,f))
        # Check
        self.assertEqual(check_illumina_qc_outputs(project,
                                                   qc_dir="qc",
                                                   qc_protocol="standardPE"),
                         [os.path.join(project.fastq_dir,
                                       "PJB1_S1_R1_001.fastq.gz")])

    def test_check_illumina_qc_outputs_standardSE_all_missing(self):
        """
        check_illumina_qc_outputs: all illumina_qc.sh outputs missing (standardSE)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Get the outputs
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        # Check
        self.assertEqual(check_illumina_qc_outputs(project,
                                                   qc_dir="qc",
                                                   qc_protocol="standardSE"),
                         project.fastqs)

    def test_check_illumina_qc_outputs_standardSE_all_present(self):
        """
        check_illumina_qc_outputs: all illumina_qc.sh outputs present (standardSE)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Add QC artefacts
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            include_fastq_strand=False,
            include_multiqc=False)
        # Check
        self.assertEqual(check_illumina_qc_outputs(project,
                                                   qc_dir="qc",
                                                   qc_protocol="standardSE"),
                         [])

    def test_check_illumina_qc_outputs_standardSE_some_missing(self):
        """
        check_illumina_qc_outputs: some illumina_qc.sh outputs missing (standardSE)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Add QC artefacts
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            include_fastq_strand=False,
            include_multiqc=False)
        # Remove some outputs
        for f in ("PJB1_S1_R1_001_fastqc.html",
                  "PJB1_S1_R1_001_model_organisms_screen.txt",):
            os.remove(os.path.join(project.qc_dir,f))
        # Check
        self.assertEqual(check_illumina_qc_outputs(project,
                                                   qc_dir="qc",
                                                   qc_protocol="standardSE"),
                         [os.path.join(project.fastq_dir,
                                       "PJB1_S1_R1_001.fastq.gz")])

    def test_check_illumina_qc_outputs_singlecell_all_missing(self):
        """
        check_illumina_qc_outputs: all illumina_qc.sh outputs missing (singlecell)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Get the outputs
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        # Check
        self.assertEqual(check_illumina_qc_outputs(project,
                                                   qc_dir="qc",
                                                   qc_protocol="singlecell"),
                         project.fastqs)

    def test_check_illumina_qc_outputs_singlecell_all_present(self):
        """
        check_illumina_qc_outputs: all illumina_qc.sh outputs present (singlecell)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Add QC artefacts
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            include_fastq_strand=False,
            include_multiqc=False)
        # Check
        self.assertEqual(check_illumina_qc_outputs(project,
                                                   qc_dir="qc",
                                                   qc_protocol="singlecell"),
                         [])

    def test_check_illumina_qc_outputs_singlecell_some_missing(self):
        """
        check_illumina_qc_outputs: some illumina_qc.sh outputs missing (singlecell)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Add QC artefacts
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            include_fastq_strand=False,
            include_multiqc=False)
        # Remove some outputs
        for f in ("PJB1_S1_R2_001_fastqc.html",
                  "PJB1_S1_R2_001_model_organisms_screen.txt",):
            os.remove(os.path.join(project.qc_dir,f))
        # Check
        self.assertEqual(check_illumina_qc_outputs(project,
                                                   qc_dir="qc",
                                                   qc_protocol="singlecell"),
                         [os.path.join(project.fastq_dir,
                                       "PJB1_S1_R2_001.fastq.gz")])

class TestCheckFastqScreenOutputs(unittest.TestCase):
    """
    Tests for the 'check_fastq_screen_outputs' function
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestCheckFastqScreenOutputs')

    def tearDown(self):
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    def test_check_fastq_screen_outputs_standardPE_all_missing(self):
        """
        check_fastq_screen_outputs: all screen outputs missing (standardPE)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Get the outputs
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        # Check
        self.assertEqual(check_fastq_screen_outputs(project,
                                                    qc_dir="qc",
                                                    screen="model_organisms",
                                                    qc_protocol="standardPE"),
                         project.fastqs)

    def test_check_fastq_screen_outputs_standardPE_all_present(self):
        """
        check_fastq_screen_outputs: all screen outputs present (standardPE)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Add QC artefacts
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            include_fastq_strand=False,
            include_multiqc=False)
        # Check
        self.assertEqual(check_fastq_screen_outputs(project,
                                                    qc_dir="qc",
                                                    screen="model_organisms",
                                                    qc_protocol="standardPE"),
                         [])

    def test_check_fastq_screen_outputs_standardPE_some_missing(self):
        """
        check_fastq_screen_outputs: some screen outputs missing (standardPE)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Add QC artefacts
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            include_fastq_strand=False,
            include_multiqc=False)
        # Remove a screen output
        os.remove(os.path.join(
            project.qc_dir,
            "PJB1_S1_R1_001_model_organisms_screen.txt"))
        # Check
        self.assertEqual(check_fastq_screen_outputs(project,
                                                    qc_dir="qc",
                                                    screen="model_organisms",
                                                    qc_protocol="standardPE"),
                         [os.path.join(project.fastq_dir,
                                       "PJB1_S1_R1_001.fastq.gz")])

    def test_check_fastq_screen_outputs_standardSE_all_missing(self):
        """
        check_fastq_screen_outputs: all screen outputs missing (standardSE)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Get the outputs
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        # Check
        self.assertEqual(check_fastq_screen_outputs(project,
                                                    qc_dir="qc",
                                                    screen="model_organisms",
                                                    qc_protocol="standardSE"),
                         project.fastqs)

    def test_check_fastq_screen_outputs_standardSE_all_present(self):
        """
        check_fastq_screen_outputs: all screen outputs present (standardSE)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Add QC artefacts
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            include_fastq_strand=False,
            include_multiqc=False)
        # Check
        self.assertEqual(check_fastq_screen_outputs(project,
                                                    qc_dir="qc",
                                                    screen="model_organisms",
                                                    qc_protocol="standardSE"),
                         [])

    def test_check_fastq_screen_outputs_standardSE_some_missing(self):
        """
        check_fastq_screen_outputs: some screen outputs missing (standardSE)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Add QC artefacts
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            include_fastq_strand=False,
            include_multiqc=False)
        # Remove a screen output
        os.remove(os.path.join(
            project.qc_dir,
            "PJB1_S1_R1_001_model_organisms_screen.txt"))
        # Check
        self.assertEqual(check_fastq_screen_outputs(project,
                                                    qc_dir="qc",
                                                    screen="model_organisms",
                                                    qc_protocol="standardSE"),
                         [os.path.join(project.fastq_dir,
                                       "PJB1_S1_R1_001.fastq.gz")])

    def test_check_fastq_screen_outputs_singlecell_all_missing(self):
        """
        check_fastq_screen_outputs: all screen outputs missing (singlecell)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Get the outputs
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        # Check
        # NB no screens expected for R1 (only R2)
        self.assertEqual(check_fastq_screen_outputs(project,
                                                    qc_dir="qc",
                                                    screen="model_organisms",
                                                    qc_protocol="singlecell"),
                         [os.path.join(project.fastq_dir,
                                       "PJB1_S1_R2_001.fastq.gz")])

    def test_check_fastq_screen_outputs_singlecell_all_present(self):
        """
        check_fastq_screen_outputs: all screen outputs present (singlecell)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Add QC artefacts
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            include_fastq_strand=False,
            include_multiqc=False)
        # Check
        self.assertEqual(check_fastq_screen_outputs(project,
                                                    qc_dir="qc",
                                                    screen="model_organisms",
                                                    qc_protocol="singlecell"),
                         [])

    def test_check_fastq_screen_outputs_singlecell_some_missing(self):
        """
        check_fastq_screen_outputs: some screen outputs missing (singlecell)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Add QC artefacts
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            include_fastq_strand=False,
            include_multiqc=False)
        # Remove a screen output
        os.remove(os.path.join(
            project.qc_dir,
            "PJB1_S1_R2_001_model_organisms_screen.txt"))
        # Check
        self.assertEqual(check_fastq_screen_outputs(project,
                                                    qc_dir="qc",
                                                    screen="model_organisms",
                                                    qc_protocol="singlecell"),
                         [os.path.join(project.fastq_dir,
                                       "PJB1_S1_R2_001.fastq.gz")])

class TestCheckFastQCOutputs(unittest.TestCase):
    """
    Tests for the 'check_fastqc_outputs' function
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestFastQCOutputs')

    def tearDown(self):
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    def test_check_fastqc_outputs_standardPE_all_missing(self):
        """
        check_fastqc_outputs: all FastQC outputs missing (standardPE)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Get the outputs
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        # Check
        self.assertEqual(check_fastqc_outputs(project,
                                              qc_dir="qc",
                                              qc_protocol="standardPE"),
                         project.fastqs)

    def test_check_fastqc_outputs_standardPE_all_present(self):
        """
        check_fastqc_outputs: all FastQC outputs present (standardPE)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Add QC artefacts
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            include_fastq_strand=False,
            include_multiqc=False)
        # Check
        self.assertEqual(check_fastqc_outputs(project,
                                              qc_dir="qc",
                                              qc_protocol="standardPE"),
                         [])

    def test_check_fastqc_outputs_standardPE_some_missing(self):
        """
        check_fastqc_outputs: some FastQC outputs missing (standardPE)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Add QC artefacts
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            include_fastq_strand=False,
            include_multiqc=False)
        # Remove a FastQC output
        os.remove(os.path.join(project.qc_dir,
                               "PJB1_S1_R1_001_fastqc.html"))
        # Check
        self.assertEqual(check_fastqc_outputs(project,
                                              qc_dir="qc",
                                              qc_protocol="standardPE"),
                         [os.path.join(project.fastq_dir,
                                       "PJB1_S1_R1_001.fastq.gz")])

    def test_check_fastqc_outputs_standardSE_all_missing(self):
        """
        check_fastqc_outputs: all FastQC outputs missing (standardSE)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Get the outputs
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        # Check
        self.assertEqual(check_fastqc_outputs(project,
                                              qc_dir="qc",
                                              qc_protocol="standardSE"),
                         project.fastqs)

    def test_check_fastqc_outputs_standardSE_all_present(self):
        """
        check_fastqc_outputs: all FastQC outputs present (standardSE)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Add QC artefacts
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            include_fastq_strand=False,
            include_multiqc=False)
        # Check
        self.assertEqual(check_fastqc_outputs(project,
                                              qc_dir="qc",
                                              qc_protocol="standardSE"),
                         [])

    def test_check_fastqc_outputs_standardSE_some_missing(self):
        """
        check_fastqc_outputs: some FastQC outputs missing (standardSE)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Add QC artefacts
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            include_fastq_strand=False,
            include_multiqc=False)
        # Remove a FastQC output
        os.remove(os.path.join(project.qc_dir,
                               "PJB1_S1_R1_001_fastqc.html"))
        # Check
        self.assertEqual(check_fastqc_outputs(project,
                                              qc_dir="qc",
                                              qc_protocol="standardSE"),
                         [os.path.join(project.fastq_dir,
                                       "PJB1_S1_R1_001.fastq.gz")])

    def test_check_fastqc_outputs_singlecell_all_missing(self):
        """
        check_fastqc_outputs: all FastQC outputs missing (singlecell)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Get the outputs
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        # Check
        self.assertEqual(check_fastqc_outputs(project,
                                              qc_dir="qc",
                                              qc_protocol="singlecell"),
                         project.fastqs)

    def test_check_fastqc_outputs_singlecell_all_present(self):
        """
        check_fastqc_outputs: all FastQC outputs present (singlecell)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Add QC artefacts
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            include_fastq_strand=False,
            include_multiqc=False)
        # Check
        self.assertEqual(check_fastqc_outputs(project,
                                              qc_dir="qc",
                                              qc_protocol="singlecell"),
                         [])

    def test_check_fastqc_outputs_singlecell_some_missing(self):
        """
        check_fastqc_outputs: some FastQC outputs missing (singlecell)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Add QC artefacts
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            include_fastq_strand=False,
            include_multiqc=False)
        # Remove a FastQC output
        os.remove(os.path.join(project.qc_dir,
                               "PJB1_S1_R1_001_fastqc.html"))
        # Check
        self.assertEqual(check_fastqc_outputs(project,
                                              qc_dir="qc",
                                              qc_protocol="singlecell"),
                         [os.path.join(project.fastq_dir,
                                       "PJB1_S1_R1_001.fastq.gz")])

class TestCheckFastqStrandOutputs(unittest.TestCase):
    """
    Tests for the 'fastq_strand_outputs' function
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestFastqStrandOutputs')

    def tearDown(self):
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    def test_check_fastq_strand_outputs_standardPE_missing(self):
        """
        check_fastq_strand_outputs: fastq_strand.py output missing (standardPE)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        # Make fastq_strand.conf
        fastq_strand_conf = os.path.join(project.dirn,"fastq_strand.conf")
        with open(fastq_strand_conf,'w') as fp:
            fp.write("")
        # Check the outputs
        self.assertEqual(check_fastq_strand_outputs(project,
                                                    "qc",
                                                    fastq_strand_conf,
                                                    qc_protocol="standardPE"),
                         [(os.path.join(project.fastq_dir,
                                        "PJB1_S1_R1_001.fastq.gz"),
                           os.path.join(project.fastq_dir,
                                        "PJB1_S1_R2_001.fastq.gz")),])

    def test_check_fastq_strand_outputs_standardPE_present(self):
        """
        check_fastq_strand_outputs: fastq_strand.py output present (standardPE)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            protocol="standardPE",
            include_fastq_strand=True,
            include_multiqc=False)
        fastq_strand_conf = os.path.join(project.dirn,"fastq_strand.conf")
        # Check the outputs
        self.assertEqual(check_fastq_strand_outputs(project,
                                                    "qc",
                                                    fastq_strand_conf,
                                                    qc_protocol="standardPE"),
                         [])

    def test_check_fastq_strand_outputs_standardSE_missing(self):
        """
        check_fastq_strand_outputs: fastq_strand.py output missing (standardSE)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        # Make fastq_strand.conf
        fastq_strand_conf = os.path.join(project.dirn,"fastq_strand.conf")
        with open(fastq_strand_conf,'w') as fp:
            fp.write("")
        # Check the outputs
        self.assertEqual(check_fastq_strand_outputs(project,
                                                    "qc",
                                                    fastq_strand_conf,
                                                    qc_protocol="standardSE"),
                         [(os.path.join(project.fastq_dir,
                                        "PJB1_S1_R1_001.fastq.gz"),),])

    def test_check_fastq_strand_outputs_standardSE_present(self):
        """
        check_fastq_strand_outputs: fastq_strand.py output present (standardSE)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            protocol="standardSE",
            include_fastq_strand=True,
            include_multiqc=False)
        fastq_strand_conf = os.path.join(project.dirn,"fastq_strand.conf")
        # Check the outputs
        self.assertEqual(check_fastq_strand_outputs(project,
                                                    "qc",
                                                    fastq_strand_conf,
                                                    qc_protocol="standardSE"),
                         [])

    def test_check_fastq_strand_outputs_singlecell_missing(self):
        """
        check_fastq_strand_outputs: fastq_strand.py output missing (singlecell)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        # Make fastq_strand.conf
        fastq_strand_conf = os.path.join(project.dirn,"fastq_strand.conf")
        with open(fastq_strand_conf,'w') as fp:
            fp.write("")
        # Check the outputs
        self.assertEqual(check_fastq_strand_outputs(project,
                                                    "qc",
                                                    fastq_strand_conf,
                                                    qc_protocol="singlecell"),
                         [(os.path.join(project.fastq_dir,
                                        "PJB1_S1_R2_001.fastq.gz"),),])

    def test_check_fastq_strand_outputs_singlecell_present(self):
        """
        check_fastq_strand_outputs: fastq_strand.py output present (singlecell)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            protocol="singlecell",
            include_fastq_strand=True,
            include_multiqc=False)
        fastq_strand_conf = os.path.join(project.dirn,"fastq_strand.conf")
        # Check the outputs
        self.assertEqual(check_fastq_strand_outputs(project,
                                                    "qc",
                                                    fastq_strand_conf,
                                                    qc_protocol="singlecell"),
                         [])

class TestCheckCellrangerCountOutputs(unittest.TestCase):
    """
    Tests for the 'check_cellranger_count_outputs' function
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestCheckCellrangerCountOutputs')

    def tearDown(self):
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    def test_check_cellranger_count_outputs_singlecell_missing(self):
        """
        check_cellranger_count_outputs: cellranger count output missing (singlecell)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           "10xGenomics Chromium 3'v3" })
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            protocol="singlecell",
            include_fastq_strand=False,
            include_multiqc=False)
        # Check the outputs
        self.assertEqual(check_cellranger_count_outputs(project),["PJB1",])

    def test_check_cellranger_count_outputs_singlecell_present(self):
        """
        check_cellranger_count_outputs: cellranger count output present (singlecell)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           "10xGenomics Chromium 3'v3" })
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            protocol="singlecell",
            include_fastq_strand=False,
            include_multiqc=False)
        UpdateAnalysisProject(project).add_cellranger_count_outputs()
        # Check the outputs
        self.assertEqual(check_cellranger_count_outputs(project),[])

    def test_check_cellranger_count_outputs_singlecell_present_with_prefix(self):
        """
        check_cellranger_count_outputs: cellranger count output present with prefix (singlecell)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           "10xGenomics Chromium 3'v3" })
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        cellranger_count_prefix = "cellranger_count/5.0.1/refdata-gex-mm10-2020-A"
        UpdateAnalysisProject(project).add_qc_outputs(
            protocol="singlecell",
            include_fastq_strand=False,
            include_multiqc=False)
        UpdateAnalysisProject(project).add_cellranger_count_outputs(
            prefix=cellranger_count_prefix)
        # Check the outputs
        self.assertEqual(
            check_cellranger_count_outputs(project,
                                           prefix=cellranger_count_prefix),[])

    def test_check_cellranger_count_outputs_10x_scRNAseq_missing(self):
        """
        check_cellranger_count_outputs: cellranger count output missing (10x_scRNAseq)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           "10xGenomics Chromium 3'v3", })
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            protocol="10x_scRNAseq",
            include_fastq_strand=False,
            include_multiqc=False)
        # Check the outputs
        self.assertEqual(check_cellranger_count_outputs(project),["PJB1",])

    def test_check_cellranger_count_outputs_10x_scRNAseq_present(self):
        """
        check_cellranger_count_outputs: cellranger count output present (10x_scRNAseq)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           "10xGenomics Chromium 3'v3" })
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            protocol="10x_scRNAseq",
            include_fastq_strand=False,
            include_multiqc=False)
        UpdateAnalysisProject(project).add_cellranger_count_outputs()
        # Check the outputs
        self.assertEqual(check_cellranger_count_outputs(project),[])

    def test_check_cellranger_count_outputs_10x_scRNAseq_present_with_prefix(self):
        """
        check_cellranger_count_outputs: cellranger count output present with prefix (10x_scRNAseq)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           "10xGenomics Chromium 3'v3" })
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        cellranger_count_prefix = "cellranger_count/5.0.1/refdata-gex-mm10-2020-A"
        UpdateAnalysisProject(project).add_qc_outputs(
            protocol="10x_scRNAseq",
            include_fastq_strand=False,
            include_multiqc=False)
        UpdateAnalysisProject(project).add_cellranger_count_outputs(
            prefix=cellranger_count_prefix)
        # Check the outputs
        self.assertEqual(
            check_cellranger_count_outputs(project,
                                           prefix=cellranger_count_prefix),[])

    def test_check_cellranger_count_outputs_10x_snRNAseq_missing(self):
        """
        check_cellranger_count_outputs: cellranger count output missing (10x_snRNAseq)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           "10xGenomics Chromium 3'v3", })
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            protocol="10x_snRNAseq",
            include_fastq_strand=False,
            include_multiqc=False)
        # Check the outputs
        self.assertEqual(check_cellranger_count_outputs(project),["PJB1",])

    def test_check_cellranger_count_outputs_10x_snRNAseq_present(self):
        """
        check_cellranger_count_outputs: cellranger count output present (10x_snRNAseq)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           "10xGenomics Chromium 3'v3" })
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            protocol="10x_snRNAseq",
            include_fastq_strand=False,
            include_multiqc=False)
        UpdateAnalysisProject(project).add_cellranger_count_outputs()
        # Check the outputs
        self.assertEqual(check_cellranger_count_outputs(project),[])

    def test_check_cellranger_count_outputs_10x_snRNAseq_present_with_prefix(self):
        """
        check_cellranger_count_outputs: cellranger count output present with prefix (10x_snRNAseq)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           "10xGenomics Chromium 3'v3" })
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        cellranger_count_prefix = "cellranger_count/3.1.0/refdata-cellranger-GRCh38-3.0.0_premRNA_patch"
        UpdateAnalysisProject(project).add_qc_outputs(
            protocol="10x_snRNAseq",
            include_fastq_strand=False,
            include_multiqc=False)
        UpdateAnalysisProject(project).add_cellranger_count_outputs(
            prefix=cellranger_count_prefix)
        # Check the outputs
        self.assertEqual(
            check_cellranger_count_outputs(project,
                                           prefix=cellranger_count_prefix),[])

class TestCheckCellrangerAtacCountOutputs(unittest.TestCase):
    """
    Tests for the 'check_cellranger_atac_count_outputs' function
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestCheckCellrangerAtacCountOutputs')

    def tearDown(self):
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    def test_check_cellranger_atac_count_outputs_10x_scATAC_missing(self):
        """
        check_cellranger_atac_count_outputs: cellranger-atac count output missing (10x_scATAC)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           "10xGenomics Single Cell ATAC" })
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            protocol="10x_scATAC",
            include_fastq_strand=False,
            include_multiqc=False)
        # Check the outputs
        self.assertEqual(check_cellranger_atac_count_outputs(project),["PJB1",])

    def test_check_cellranger_atac_count_outputs_10x_scATAC_present(self):
        """
        check_cellranger_atac_count_outputs: cellranger-atac count output present (10x_scATAC)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           "10xGenomics Single Cell ATAC" })
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            protocol="10x_scATAC",
            include_fastq_strand=False,
            include_multiqc=False)
        UpdateAnalysisProject(project).add_cellranger_count_outputs(
            cellranger="cellranger-atac")
        # Check the outputs
        self.assertEqual(check_cellranger_atac_count_outputs(project),[])

    def test_check_cellranger_atac_count_outputs_10x_scATAC_present_with_prefix(self):
        """
        check_cellranger_atac_count_outputs: cellranger-atac count output present with prefix (10x_scATAC)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           "10xGenomics Single Cell ATAC" })
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        cellranger_count_prefix = "cellranger_count/1.2.0/refdata-cellranger-atac-GRCh38-1.2.0"
        UpdateAnalysisProject(project).add_qc_outputs(
            protocol="10x_scATAC",
            include_fastq_strand=False,
            include_multiqc=False)
        UpdateAnalysisProject(project).add_cellranger_count_outputs(
            cellranger="cellranger-atac",
            prefix=cellranger_count_prefix)
        # Check the outputs
        self.assertEqual(
            check_cellranger_atac_count_outputs(project,
                                                prefix=cellranger_count_prefix),
            [])

class TestCheckCellrangerArcCountOutputs(unittest.TestCase):
    """
    Tests for the 'check_cellranger_arc_count_outputs' function
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestCheckCellrangerArcCountOutputs')

    def tearDown(self):
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    def test_check_cellranger_arc_count_outputs_no_libraries_csv(self):
        """
        check_cellranger_arc_count_outputs: cellranger-arc count no libraries.csv (10x_Multiome_*)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           "10xGenomics Single Cell Multiome",
                                           'Library type': 'GEX' })
        p.create(top_dir=self.wd)
        project = AnalysisProject(os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            protocol="10x_Multiome_GEX",
            include_fastq_strand=False,
            include_multiqc=False)
        # Check the outputs
        self.assertEqual(check_cellranger_arc_count_outputs(project),[])

    def test_check_cellranger_arc_count_outputs_10x_multiome_missing(self):
        """
        check_cellranger_arc_count_outputs: cellranger-arc count output missing (10x_Multiome_*)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           "10xGenomics Single Cell Multiome",
                                           'Library type': 'GEX' })
        p.create(top_dir=self.wd)
        project = AnalysisProject(os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            protocol="10x_Multiome_GEX",
            include_fastq_strand=False,
            include_multiqc=False)
        # Add libraries.csv
        with open(os.path.join(project.qc_dir,"libraries.PJB1.csv"),'wt') as fp:
            fp.write("")
        # Check the outputs
        self.assertEqual(check_cellranger_arc_count_outputs(project),["PJB1",])

    def test_check_cellranger_arc_count_outputs_10x_multiome_present(self):
        """
        check_cellranger_arc_count_outputs: cellranger-arc count output present (10x_Multiome_*)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           "10xGenomics Single Cell Multiome",
                                           'Library type': 'GEX' })
        p.create(top_dir=self.wd)
        project = AnalysisProject(os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            protocol="10x_Multiome_GEX",
            include_fastq_strand=False,
            include_multiqc=False)
        UpdateAnalysisProject(project).add_cellranger_count_outputs(
            cellranger="cellranger-arc")
        # Add libraries.csv
        with open(os.path.join(project.qc_dir,"libraries.PJB1.csv"),'wt') as fp:
            fp.write("")
        # Check the outputs
        self.assertEqual(check_cellranger_arc_count_outputs(project),[])

    def test_check_cellranger_arc_count_outputs_10x_multiome_present_with_prefix(self):
        """
        check_cellranger_arc_count_outputs: cellranger-arc count output present with prefix (10x_Multiome_*)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           "10xGenomics Single Cell Multiome",
                                           'Library type': 'GEX' })
        p.create(top_dir=self.wd)
        project = AnalysisProject(os.path.join(self.wd,"PJB"))
        cellranger_count_prefix = "cellranger_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A"
        UpdateAnalysisProject(project).add_qc_outputs(
            protocol="10x_Multiome_GEX",
            include_fastq_strand=False,
            include_multiqc=False)
        UpdateAnalysisProject(project).add_cellranger_count_outputs(
            cellranger="cellranger-arc",
            prefix=cellranger_count_prefix)
        # Add libraries.csv
        with open(os.path.join(project.qc_dir,"libraries.PJB1.csv"),'wt') as fp:
            fp.write("")
        # Check the outputs
        self.assertEqual(
            check_cellranger_arc_count_outputs(project,
                                           prefix=cellranger_count_prefix),[])

class TestExpectedOutputs(unittest.TestCase):
    """
    Tests for the 'expected_outputs' function
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestExpectedOutputs')

    def tearDown(self):
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    def test_expected_outputs_standardPE(self):
        """
        expected_outputs: standard paired-end, no strandedness
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Reference outputs
        reference_outputs = ("PJB1_S1_R1_001_fastqc",
                             "PJB1_S1_R1_001_fastqc.html",
                             "PJB1_S1_R1_001_fastqc.zip",
                             "PJB1_S1_R1_001_model_organisms_screen.png",
                             "PJB1_S1_R1_001_model_organisms_screen.txt",
                             "PJB1_S1_R1_001_other_organisms_screen.png",
                             "PJB1_S1_R1_001_other_organisms_screen.txt",
                             "PJB1_S1_R1_001_rRNA_screen.png",
                             "PJB1_S1_R1_001_rRNA_screen.txt",
                             "PJB1_S1_R2_001_fastqc",
                             "PJB1_S1_R2_001_fastqc.html",
                             "PJB1_S1_R2_001_fastqc.zip",
                             "PJB1_S1_R2_001_model_organisms_screen.png",
                             "PJB1_S1_R2_001_model_organisms_screen.txt",
                             "PJB1_S1_R2_001_other_organisms_screen.png",
                             "PJB1_S1_R2_001_other_organisms_screen.txt",
                             "PJB1_S1_R2_001_rRNA_screen.png",
                             "PJB1_S1_R2_001_rRNA_screen.txt",)
        expected = expected_outputs(AnalysisProject(p.name,
                                                    os.path.join(self.wd,
                                                                 p.name)),
                                    "qc",
                                    qc_protocol="standardPE",
                                    fastq_screens=('model_organisms',
                                                   'other_organisms',
                                                   'rRNA'))
        for e in expected:
            self.assertEqual(os.path.dirname(e),os.path.join(self.wd,
                                                             p.name,
                                                             "qc"))
            self.assertTrue(os.path.basename(e) in reference_outputs)
        for r in reference_outputs:
            self.assertTrue(os.path.join(self.wd,p.name,"qc",r) in expected)

    def test_expected_outputs_standardPE_no_screens(self):
        """
        expected_outputs: standard paired-end, no screens
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Reference outputs
        reference_outputs = ("PJB1_S1_R1_001_fastqc",
                             "PJB1_S1_R1_001_fastqc.html",
                             "PJB1_S1_R1_001_fastqc.zip",
                             "PJB1_S1_R2_001_fastqc",
                             "PJB1_S1_R2_001_fastqc.html",
                             "PJB1_S1_R2_001_fastqc.zip",)
        expected = expected_outputs(AnalysisProject(p.name,
                                                    os.path.join(self.wd,
                                                                 p.name)),
                                    "qc",
                                    qc_protocol="standardPE")
        for e in expected:
            self.assertEqual(os.path.dirname(e),os.path.join(self.wd,
                                                             p.name,
                                                             "qc"))
            self.assertTrue(os.path.basename(e) in reference_outputs)
        for r in reference_outputs:
            self.assertTrue(os.path.join(self.wd,p.name,"qc",r) in expected)

    def test_expected_outputs_standardPE_with_strand(self):
        """
        expected_outputs: standard paired-end with strandedness
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Make mock fastq_strand
        mock_fastq_strand_conf = os.path.join(self.wd,
                                              p.name,
                                              "fastq_strand.conf")
        with open(mock_fastq_strand_conf,'w') as fp:
            fp.write("")
        # Reference outputs
        reference_outputs = ("PJB1_S1_R1_001_fastqc",
                             "PJB1_S1_R1_001_fastqc.html",
                             "PJB1_S1_R1_001_fastqc.zip",
                             "PJB1_S1_R1_001_model_organisms_screen.png",
                             "PJB1_S1_R1_001_model_organisms_screen.txt",
                             "PJB1_S1_R1_001_other_organisms_screen.png",
                             "PJB1_S1_R1_001_other_organisms_screen.txt",
                             "PJB1_S1_R1_001_rRNA_screen.png",
                             "PJB1_S1_R1_001_rRNA_screen.txt",
                             "PJB1_S1_R2_001_fastqc",
                             "PJB1_S1_R2_001_fastqc.html",
                             "PJB1_S1_R2_001_fastqc.zip",
                             "PJB1_S1_R2_001_model_organisms_screen.png",
                             "PJB1_S1_R2_001_model_organisms_screen.txt",
                             "PJB1_S1_R2_001_other_organisms_screen.png",
                             "PJB1_S1_R2_001_other_organisms_screen.txt",
                             "PJB1_S1_R2_001_rRNA_screen.png",
                             "PJB1_S1_R2_001_rRNA_screen.txt",
                             "PJB1_S1_R1_001_fastq_strand.txt",)
        expected = expected_outputs(AnalysisProject(p.name,
                                                    os.path.join(self.wd,
                                                                 p.name)),
                                    "qc",
                                    fastq_strand_conf=mock_fastq_strand_conf,
                                    qc_protocol="standardPE",
                                    fastq_screens=('model_organisms',
                                                   'other_organisms',
                                                   'rRNA'))
        for e in expected:
            self.assertEqual(os.path.dirname(e),os.path.join(self.wd,
                                                             p.name,
                                                             "qc"))
            self.assertTrue(os.path.basename(e) in reference_outputs)
        for r in reference_outputs:
            self.assertTrue(os.path.join(self.wd,p.name,"qc",r) in expected)

    def test_expected_outputs_standardSE(self):
        """
        expected_outputs: standard single-end, no strandedness
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Reference outputs
        reference_outputs = ("PJB1_S1_R1_001_fastqc",
                             "PJB1_S1_R1_001_fastqc.html",
                             "PJB1_S1_R1_001_fastqc.zip",
                             "PJB1_S1_R1_001_model_organisms_screen.png",
                             "PJB1_S1_R1_001_model_organisms_screen.txt",
                             "PJB1_S1_R1_001_other_organisms_screen.png",
                             "PJB1_S1_R1_001_other_organisms_screen.txt",
                             "PJB1_S1_R1_001_rRNA_screen.png",
                             "PJB1_S1_R1_001_rRNA_screen.txt",)
        expected = expected_outputs(AnalysisProject(p.name,
                                                    os.path.join(self.wd,
                                                                 p.name)),
                                    "qc",
                                    qc_protocol="standardSE",
                                    fastq_screens=('model_organisms',
                                                   'other_organisms',
                                                   'rRNA'))
        for e in expected:
            self.assertEqual(os.path.dirname(e),os.path.join(self.wd,
                                                             p.name,
                                                             "qc"))
            self.assertTrue(os.path.basename(e) in reference_outputs)
        for r in reference_outputs:
            self.assertTrue(os.path.join(self.wd,p.name,"qc",r) in expected)

    def test_expected_outputs_standardSE_with_strand(self):
        """
        expected_outputs: standard single-end with strandedness
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Make mock fastq_strand
        mock_fastq_strand_conf = os.path.join(self.wd,
                                              p.name,
                                              "fastq_strand.conf")
        with open(mock_fastq_strand_conf,'w') as fp:
            fp.write("")
        # Reference outputs
        reference_outputs = ("PJB1_S1_R1_001_fastqc",
                             "PJB1_S1_R1_001_fastqc.html",
                             "PJB1_S1_R1_001_fastqc.zip",
                             "PJB1_S1_R1_001_model_organisms_screen.png",
                             "PJB1_S1_R1_001_model_organisms_screen.txt",
                             "PJB1_S1_R1_001_other_organisms_screen.png",
                             "PJB1_S1_R1_001_other_organisms_screen.txt",
                             "PJB1_S1_R1_001_rRNA_screen.png",
                             "PJB1_S1_R1_001_rRNA_screen.txt",
                             "PJB1_S1_R1_001_fastq_strand.txt",)
        expected = expected_outputs(AnalysisProject(p.name,
                                                    os.path.join(self.wd,
                                                                 p.name)),
                                    "qc",
                                    fastq_strand_conf=mock_fastq_strand_conf,
                                    qc_protocol="standardSE",
                                    fastq_screens=('model_organisms',
                                                   'other_organisms',
                                                   'rRNA'))
        for e in expected:
            self.assertEqual(os.path.dirname(e),os.path.join(self.wd,
                                                             p.name,
                                                             "qc"))
            self.assertTrue(os.path.basename(e) in reference_outputs)
        for r in reference_outputs:
            self.assertTrue(os.path.join(self.wd,p.name,"qc",r) in expected)

    def test_expected_outputs_singlecell(self):
        """
        expected_outputs: single-cell, no strandedness
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Reference outputs
        reference_outputs = ("PJB1_S1_R1_001_fastqc",
                             "PJB1_S1_R1_001_fastqc.html",
                             "PJB1_S1_R1_001_fastqc.zip",
                             "PJB1_S1_R2_001_fastqc",
                             "PJB1_S1_R2_001_fastqc.html",
                             "PJB1_S1_R2_001_fastqc.zip",
                             "PJB1_S1_R2_001_model_organisms_screen.png",
                             "PJB1_S1_R2_001_model_organisms_screen.txt",
                             "PJB1_S1_R2_001_other_organisms_screen.png",
                             "PJB1_S1_R2_001_other_organisms_screen.txt",
                             "PJB1_S1_R2_001_rRNA_screen.png",
                             "PJB1_S1_R2_001_rRNA_screen.txt",)
        expected = expected_outputs(AnalysisProject(p.name,
                                                    os.path.join(self.wd,
                                                                 p.name)),
                                    "qc",
                                    qc_protocol="singlecell",
                                    fastq_screens=('model_organisms',
                                                   'other_organisms',
                                                   'rRNA'))
        for e in expected:
            self.assertEqual(os.path.dirname(e),os.path.join(self.wd,
                                                             p.name,
                                                             "qc"))
            self.assertTrue(os.path.basename(e) in reference_outputs)
        for r in reference_outputs:
            self.assertTrue(os.path.join(self.wd,p.name,"qc",r) in expected)

    def test_expected_outputs_singlecell_with_strand(self):
        """
        expected_outputs: single-cell with strandedness
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Make mock fastq_strand
        mock_fastq_strand_conf = os.path.join(self.wd,
                                              p.name,
                                              "fastq_strand.conf")
        with open(mock_fastq_strand_conf,'w') as fp:
            fp.write("")
        # Reference outputs
        reference_outputs = ("PJB1_S1_R1_001_fastqc",
                             "PJB1_S1_R1_001_fastqc.html",
                             "PJB1_S1_R1_001_fastqc.zip",
                             "PJB1_S1_R2_001_fastqc",
                             "PJB1_S1_R2_001_fastqc.html",
                             "PJB1_S1_R2_001_fastqc.zip",
                             "PJB1_S1_R2_001_model_organisms_screen.png",
                             "PJB1_S1_R2_001_model_organisms_screen.txt",
                             "PJB1_S1_R2_001_other_organisms_screen.png",
                             "PJB1_S1_R2_001_other_organisms_screen.txt",
                             "PJB1_S1_R2_001_rRNA_screen.png",
                             "PJB1_S1_R2_001_rRNA_screen.txt",
                             "PJB1_S1_R2_001_fastq_strand.txt",)
        expected = expected_outputs(AnalysisProject(p.name,
                                                    os.path.join(self.wd,
                                                                 p.name)),
                                    "qc",
                                    fastq_strand_conf=mock_fastq_strand_conf,
                                    qc_protocol="singlecell",
                                    fastq_screens=('model_organisms',
                                                   'other_organisms',
                                                   'rRNA'))
        for e in expected:
            self.assertEqual(os.path.dirname(e),os.path.join(self.wd,
                                                             p.name,
                                                             "qc"))
            self.assertTrue(os.path.basename(e) in reference_outputs)
        for r in reference_outputs:
            self.assertTrue(os.path.join(self.wd,p.name,"qc",r) in expected)

    def test_expected_outputs_10x_scRNA_seq_with_cellranger(self):
        """
        expected_outputs: 10xGenomics scRNA-seq with cellranger
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           "10xGenomics Chromium 3'v3" })
        p.create(top_dir=self.wd)
        # Make mock fastq_strand
        mock_fastq_strand_conf = os.path.join(self.wd,
                                              p.name,
                                              "fastq_strand.conf")
        with open(mock_fastq_strand_conf,'w') as fp:
            fp.write("")
        reference_outputs = ("qc/PJB1_S1_R1_001_fastqc",
                             "qc/PJB1_S1_R1_001_fastqc.html",
                             "qc/PJB1_S1_R1_001_fastqc.zip",
                             "qc/PJB1_S1_R2_001_fastqc",
                             "qc/PJB1_S1_R2_001_fastqc.html",
                             "qc/PJB1_S1_R2_001_fastqc.zip",
                             "qc/PJB1_S1_R2_001_model_organisms_screen.png",
                             "qc/PJB1_S1_R2_001_model_organisms_screen.txt",
                             "qc/PJB1_S1_R2_001_other_organisms_screen.png",
                             "qc/PJB1_S1_R2_001_other_organisms_screen.txt",
                             "qc/PJB1_S1_R2_001_rRNA_screen.png",
                             "qc/PJB1_S1_R2_001_rRNA_screen.txt",
                             "qc/PJB1_S1_R2_001_fastq_strand.txt",
                             "qc/cellranger_count/PJB1/outs/metrics_summary.csv",
                             "qc/cellranger_count/PJB1/outs/web_summary.html",)
        expected = expected_outputs(AnalysisProject(p.name,
                                                    os.path.join(self.wd,
                                                                 p.name)),
                                    "qc",
                                    fastq_strand_conf=mock_fastq_strand_conf,
                                    cellranger_refdata=
                                    "/data/refdata-cellranger-GRCh38-1.2.0",
                                    qc_protocol="10x_scRNAseq",
                                    fastq_screens=('model_organisms',
                                                   'other_organisms',
                                                   'rRNA'))
        for e in expected:
            print(e)
            self.assertTrue(e in [os.path.join(self.wd,p.name,r)
                                  for r in reference_outputs],
                            "%s not found in reference" % e)
        for r in reference_outputs:
            self.assertTrue(os.path.join(self.wd,p.name,r) in expected,
                            "%s not found in expected" % r)

    def test_expected_outputs_10x_scRNA_seq_with_cellranger_and_versioning(self):
        """
        expected_outputs: 10xGenomics scRNA-seq with cellranger (including versioning)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           "10xGenomics Chromium 3'v3" })
        p.create(top_dir=self.wd)
        # Make mock fastq_strand
        mock_fastq_strand_conf = os.path.join(self.wd,
                                              p.name,
                                              "fastq_strand.conf")
        with open(mock_fastq_strand_conf,'w') as fp:
            fp.write("")
        reference_outputs = ("qc/PJB1_S1_R1_001_fastqc",
                             "qc/PJB1_S1_R1_001_fastqc.html",
                             "qc/PJB1_S1_R1_001_fastqc.zip",
                             "qc/PJB1_S1_R2_001_fastqc",
                             "qc/PJB1_S1_R2_001_fastqc.html",
                             "qc/PJB1_S1_R2_001_fastqc.zip",
                             "qc/PJB1_S1_R2_001_model_organisms_screen.png",
                             "qc/PJB1_S1_R2_001_model_organisms_screen.txt",
                             "qc/PJB1_S1_R2_001_other_organisms_screen.png",
                             "qc/PJB1_S1_R2_001_other_organisms_screen.txt",
                             "qc/PJB1_S1_R2_001_rRNA_screen.png",
                             "qc/PJB1_S1_R2_001_rRNA_screen.txt",
                             "qc/PJB1_S1_R2_001_fastq_strand.txt",
                             "qc/cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB1/outs/metrics_summary.csv",
                             "qc/cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB1/outs/web_summary.html",)
        expected = expected_outputs(AnalysisProject(p.name,
                                                    os.path.join(self.wd,
                                                                 p.name)),
                                    "qc",
                                    fastq_strand_conf=mock_fastq_strand_conf,
                                    cellranger_version="5.0.1",
                                    cellranger_refdata=
                                    "/data/refdata-gex-GRCh38-2020-A",
                                    qc_protocol="10x_scRNAseq",
                                    fastq_screens=('model_organisms',
                                                   'other_organisms',
                                                   'rRNA'))
        for e in expected:
            print(e)
            self.assertTrue(e in [os.path.join(self.wd,p.name,r)
                                  for r in reference_outputs],
                            "%s not found in reference" % e)
        for r in reference_outputs:
            self.assertTrue(os.path.join(self.wd,p.name,r) in expected,
                            "%s not found in expected" % r)

    def test_expected_outputs_10x_scATAC_with_cellranger_atac(self):
        """
        expected_outputs: 10xGenomics scATAC-seq with cellranger-atac
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB1_S1_R3_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           "10xGenomics Single Cell ATAC" })
        p.create(top_dir=self.wd)
        # Make mock fastq_strand
        mock_fastq_strand_conf = os.path.join(self.wd,
                                              p.name,
                                              "fastq_strand.conf")
        with open(mock_fastq_strand_conf,'w') as fp:
            fp.write("")
        reference_outputs = ("qc/PJB1_S1_R1_001_fastqc",
                             "qc/PJB1_S1_R1_001_fastqc.html",
                             "qc/PJB1_S1_R1_001_fastqc.zip",
                             "qc/PJB1_S1_R1_001_fastqc",
                             "qc/PJB1_S1_R1_001_fastqc.html",
                             "qc/PJB1_S1_R1_001_fastqc.zip",
                             "qc/PJB1_S1_R1_001_model_organisms_screen.png",
                             "qc/PJB1_S1_R1_001_model_organisms_screen.txt",
                             "qc/PJB1_S1_R1_001_other_organisms_screen.png",
                             "qc/PJB1_S1_R1_001_other_organisms_screen.txt",
                             "qc/PJB1_S1_R1_001_rRNA_screen.png",
                             "qc/PJB1_S1_R1_001_rRNA_screen.txt",
                             "qc/PJB1_S1_R1_001_fastq_strand.txt",
                             "qc/PJB1_S1_R3_001_fastqc",
                             "qc/PJB1_S1_R3_001_fastqc.html",
                             "qc/PJB1_S1_R3_001_fastqc.zip",
                             "qc/PJB1_S1_R3_001_fastqc",
                             "qc/PJB1_S1_R3_001_fastqc.html",
                             "qc/PJB1_S1_R3_001_fastqc.zip",
                             "qc/PJB1_S1_R3_001_model_organisms_screen.png",
                             "qc/PJB1_S1_R3_001_model_organisms_screen.txt",
                             "qc/PJB1_S1_R3_001_other_organisms_screen.png",
                             "qc/PJB1_S1_R3_001_other_organisms_screen.txt",
                             "qc/PJB1_S1_R3_001_rRNA_screen.png",
                             "qc/PJB1_S1_R3_001_rRNA_screen.txt",
                             "qc/cellranger_count/PJB1/outs/summary.csv",
                             "qc/cellranger_count/PJB1/outs/web_summary.html",)
        expected = expected_outputs(AnalysisProject(p.name,
                                                    os.path.join(self.wd,
                                                                 p.name)),
                                    "qc",
                                    fastq_strand_conf=mock_fastq_strand_conf,
                                    cellranger_refdata=
                                    "/data/refdata-cellranger-atac-GRCh38-1.2.0",
                                    qc_protocol="10x_scATAC",
                                    fastq_screens=('model_organisms',
                                                   'other_organisms',
                                                   'rRNA'))
        for e in expected:
            print(e)
            self.assertTrue(e in [os.path.join(self.wd,p.name,r)
                                  for r in reference_outputs],
                            "%s not found in reference" % e)
        for r in reference_outputs:
            self.assertTrue(os.path.join(self.wd,p.name,r) in expected,
                            "%s not found in expected" % r)

    def test_expected_outputs_10x_scATAC_with_cellranger_atac_with_versioning(self):
        """
        expected_outputs: 10xGenomics scATAC-seq with cellranger-atac (including versioning)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB1_S1_R3_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           "10xGenomics Single Cell ATAC" })
        p.create(top_dir=self.wd)
        # Make mock fastq_strand
        mock_fastq_strand_conf = os.path.join(self.wd,
                                              p.name,
                                              "fastq_strand.conf")
        with open(mock_fastq_strand_conf,'w') as fp:
            fp.write("")
        reference_outputs = ("qc/PJB1_S1_R1_001_fastqc",
                             "qc/PJB1_S1_R1_001_fastqc.html",
                             "qc/PJB1_S1_R1_001_fastqc.zip",
                             "qc/PJB1_S1_R1_001_fastqc",
                             "qc/PJB1_S1_R1_001_fastqc.html",
                             "qc/PJB1_S1_R1_001_fastqc.zip",
                             "qc/PJB1_S1_R1_001_model_organisms_screen.png",
                             "qc/PJB1_S1_R1_001_model_organisms_screen.txt",
                             "qc/PJB1_S1_R1_001_other_organisms_screen.png",
                             "qc/PJB1_S1_R1_001_other_organisms_screen.txt",
                             "qc/PJB1_S1_R1_001_rRNA_screen.png",
                             "qc/PJB1_S1_R1_001_rRNA_screen.txt",
                             "qc/PJB1_S1_R1_001_fastq_strand.txt",
                             "qc/PJB1_S1_R3_001_fastqc",
                             "qc/PJB1_S1_R3_001_fastqc.html",
                             "qc/PJB1_S1_R3_001_fastqc.zip",
                             "qc/PJB1_S1_R3_001_fastqc",
                             "qc/PJB1_S1_R3_001_fastqc.html",
                             "qc/PJB1_S1_R3_001_fastqc.zip",
                             "qc/PJB1_S1_R3_001_model_organisms_screen.png",
                             "qc/PJB1_S1_R3_001_model_organisms_screen.txt",
                             "qc/PJB1_S1_R3_001_other_organisms_screen.png",
                             "qc/PJB1_S1_R3_001_other_organisms_screen.txt",
                             "qc/PJB1_S1_R3_001_rRNA_screen.png",
                             "qc/PJB1_S1_R3_001_rRNA_screen.txt",
                             "qc/cellranger_count/1.2.0/refdata-cellranger-atac-GRCh38-1.2.0/PJB1/outs/summary.csv",
                             "qc/cellranger_count/1.2.0/refdata-cellranger-atac-GRCh38-1.2.0/PJB1/outs/web_summary.html",)
        expected = expected_outputs(AnalysisProject(p.name,
                                                    os.path.join(self.wd,
                                                                 p.name)),
                                    "qc",
                                    fastq_strand_conf=mock_fastq_strand_conf,
                                    cellranger_version="1.2.0",
                                    cellranger_refdata=
                                    "/data/refdata-cellranger-atac-GRCh38-1.2.0",
                                    qc_protocol="10x_scATAC",
                                    fastq_screens=('model_organisms',
                                                   'other_organisms',
                                                   'rRNA'))
        for e in expected:
            print(e)
            self.assertTrue(e in [os.path.join(self.wd,p.name,r)
                                  for r in reference_outputs],
                            "%s not found in reference" % e)
        for r in reference_outputs:
            self.assertTrue(os.path.join(self.wd,p.name,r) in expected,
                            "%s not found in expected" % r)

    def test_expected_outputs_10x_scRNA_seq_with_no_cellranger_reference(self):
        """
        expected_outputs: 10xGenomics scRNA-seq with no cellranger reference data
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           "10xGenomics Chromium 3'v3" })
        p.create(top_dir=self.wd)
        # Make mock fastq_strand
        mock_fastq_strand_conf = os.path.join(self.wd,
                                              p.name,
                                              "fastq_strand.conf")
        with open(mock_fastq_strand_conf,'w') as fp:
            fp.write("")
        reference_outputs = ("qc/PJB1_S1_R1_001_fastqc",
                             "qc/PJB1_S1_R1_001_fastqc.html",
                             "qc/PJB1_S1_R1_001_fastqc.zip",
                             "qc/PJB1_S1_R2_001_fastqc",
                             "qc/PJB1_S1_R2_001_fastqc.html",
                             "qc/PJB1_S1_R2_001_fastqc.zip",
                             "qc/PJB1_S1_R2_001_model_organisms_screen.png",
                             "qc/PJB1_S1_R2_001_model_organisms_screen.txt",
                             "qc/PJB1_S1_R2_001_other_organisms_screen.png",
                             "qc/PJB1_S1_R2_001_other_organisms_screen.txt",
                             "qc/PJB1_S1_R2_001_rRNA_screen.png",
                             "qc/PJB1_S1_R2_001_rRNA_screen.txt",
                             "qc/PJB1_S1_R2_001_fastq_strand.txt",)
        expected = expected_outputs(AnalysisProject(p.name,
                                                    os.path.join(self.wd,
                                                                 p.name)),
                                    "qc",
                                    fastq_strand_conf=mock_fastq_strand_conf,
                                    cellranger_refdata=None,
                                    qc_protocol="10x_scRNAseq",
                                    fastq_screens=('model_organisms',
                                                   'other_organisms',
                                                   'rRNA'))
        for e in expected:
            print(e)
            self.assertTrue(e in [os.path.join(self.wd,p.name,r)
                                  for r in reference_outputs],
                            "%s not found in reference" % e)
        for r in reference_outputs:
            self.assertTrue(os.path.join(self.wd,p.name,r) in expected,
                            "%s not found in expected" % r)

    def test_expected_outputs_10x_visium(self):
        """
        expected_outputs: 10xGenomics Visium
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           "10xGenomics Visium" })
        p.create(top_dir=self.wd)
        # Make mock fastq_strand
        mock_fastq_strand_conf = os.path.join(self.wd,
                                              p.name,
                                              "fastq_strand.conf")
        with open(mock_fastq_strand_conf,'w') as fp:
            fp.write("")
        reference_outputs = ("qc/PJB1_S1_R1_001_fastqc",
                             "qc/PJB1_S1_R1_001_fastqc.html",
                             "qc/PJB1_S1_R1_001_fastqc.zip",
                             "qc/PJB1_S1_R2_001_fastqc",
                             "qc/PJB1_S1_R2_001_fastqc.html",
                             "qc/PJB1_S1_R2_001_fastqc.zip",
                             "qc/PJB1_S1_R2_001_model_organisms_screen.png",
                             "qc/PJB1_S1_R2_001_model_organisms_screen.txt",
                             "qc/PJB1_S1_R2_001_other_organisms_screen.png",
                             "qc/PJB1_S1_R2_001_other_organisms_screen.txt",
                             "qc/PJB1_S1_R2_001_rRNA_screen.png",
                             "qc/PJB1_S1_R2_001_rRNA_screen.txt",
                             "qc/PJB1_S1_R2_001_fastq_strand.txt",)
        expected = expected_outputs(AnalysisProject(p.name,
                                                    os.path.join(self.wd,
                                                                 p.name)),
                                    "qc",
                                    fastq_strand_conf=mock_fastq_strand_conf,
                                    cellranger_refdata=None,
                                    qc_protocol="10x_Visium",
                                    fastq_screens=('model_organisms',
                                                   'other_organisms',
                                                   'rRNA'))
        for e in expected:
            print(e)
            self.assertTrue(e in [os.path.join(self.wd,p.name,r)
                                  for r in reference_outputs],
                            "%s not found in reference" % e)
        for r in reference_outputs:
            self.assertTrue(os.path.join(self.wd,p.name,r) in expected,
                            "%s not found in expected" % r)

    def test_expected_outputs_10x_multiome_atac(self):
        """
        expected_outputs: 10xGenomics Multiome ATAC
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB1_S1_R3_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           '10xGenomics Single Cell Multiome',
                                           'Library type': 'ATAC' })
        p.create(top_dir=self.wd)
        # Make mock fastq_strand
        mock_fastq_strand_conf = os.path.join(self.wd,
                                              p.name,
                                              "fastq_strand.conf")
        with open(mock_fastq_strand_conf,'w') as fp:
            fp.write("")
        reference_outputs = ("qc/PJB1_S1_R1_001_fastqc",
                             "qc/PJB1_S1_R1_001_fastqc.html",
                             "qc/PJB1_S1_R1_001_fastqc.zip",
                             "qc/PJB1_S1_R1_001_fastqc",
                             "qc/PJB1_S1_R1_001_fastqc.html",
                             "qc/PJB1_S1_R1_001_fastqc.zip",
                             "qc/PJB1_S1_R1_001_model_organisms_screen.png",
                             "qc/PJB1_S1_R1_001_model_organisms_screen.txt",
                             "qc/PJB1_S1_R1_001_other_organisms_screen.png",
                             "qc/PJB1_S1_R1_001_other_organisms_screen.txt",
                             "qc/PJB1_S1_R1_001_rRNA_screen.png",
                             "qc/PJB1_S1_R1_001_rRNA_screen.txt",
                             "qc/PJB1_S1_R1_001_fastq_strand.txt",
                             "qc/PJB1_S1_R3_001_fastqc",
                             "qc/PJB1_S1_R3_001_fastqc.html",
                             "qc/PJB1_S1_R3_001_fastqc.zip",
                             "qc/PJB1_S1_R3_001_fastqc",
                             "qc/PJB1_S1_R3_001_fastqc.html",
                             "qc/PJB1_S1_R3_001_fastqc.zip",
                             "qc/PJB1_S1_R3_001_model_organisms_screen.png",
                             "qc/PJB1_S1_R3_001_model_organisms_screen.txt",
                             "qc/PJB1_S1_R3_001_other_organisms_screen.png",
                             "qc/PJB1_S1_R3_001_other_organisms_screen.txt",
                             "qc/PJB1_S1_R3_001_rRNA_screen.png",
                             "qc/PJB1_S1_R3_001_rRNA_screen.txt",)
        expected = expected_outputs(AnalysisProject(p.name,
                                                    os.path.join(self.wd,
                                                                 p.name)),
                                    "qc",
                                    fastq_strand_conf=mock_fastq_strand_conf,
                                    cellranger_refdata=None,
                                    qc_protocol="10x_Multiome_ATAC",
                                    fastq_screens=('model_organisms',
                                                   'other_organisms',
                                                   'rRNA'))
        for e in expected:
            print(e)
            self.assertTrue(e in [os.path.join(self.wd,p.name,r)
                                  for r in reference_outputs],
                            "%s not found in reference" % e)
        for r in reference_outputs:
            self.assertTrue(os.path.join(self.wd,p.name,r) in expected,
                            "%s not found in expected" % r)

    def test_expected_outputs_10x_multiome_gex(self):
        """
        expected_outputs: 10xGenomics Multiome GEX
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           '10xGenomics Single Cell Multiome',
                                           'Library type': 'GEX' })
        p.create(top_dir=self.wd)
        # Make mock fastq_strand
        mock_fastq_strand_conf = os.path.join(self.wd,
                                              p.name,
                                              "fastq_strand.conf")
        with open(mock_fastq_strand_conf,'w') as fp:
            fp.write("")
        reference_outputs = ("qc/PJB1_S1_R1_001_fastqc",
                             "qc/PJB1_S1_R1_001_fastqc.html",
                             "qc/PJB1_S1_R1_001_fastqc.zip",
                             "qc/PJB1_S1_R1_001_fastqc",
                             "qc/PJB1_S1_R1_001_fastqc.html",
                             "qc/PJB1_S1_R1_001_fastqc.zip",
                             "qc/PJB1_S1_R2_001_fastqc",
                             "qc/PJB1_S1_R2_001_fastqc.html",
                             "qc/PJB1_S1_R2_001_fastqc.zip",
                             "qc/PJB1_S1_R2_001_fastqc",
                             "qc/PJB1_S1_R2_001_fastqc.html",
                             "qc/PJB1_S1_R2_001_fastqc.zip",
                             "qc/PJB1_S1_R2_001_fastq_strand.txt",
                             "qc/PJB1_S1_R2_001_model_organisms_screen.png",
                             "qc/PJB1_S1_R2_001_model_organisms_screen.txt",
                             "qc/PJB1_S1_R2_001_other_organisms_screen.png",
                             "qc/PJB1_S1_R2_001_other_organisms_screen.txt",
                             "qc/PJB1_S1_R2_001_rRNA_screen.png",
                             "qc/PJB1_S1_R2_001_rRNA_screen.txt",)
        expected = expected_outputs(AnalysisProject(p.name,
                                                    os.path.join(self.wd,
                                                                 p.name)),
                                    "qc",
                                    fastq_strand_conf=mock_fastq_strand_conf,
                                    cellranger_refdata=None,
                                    qc_protocol="10x_Multiome_GEX",
                                    fastq_screens=('model_organisms',
                                                   'other_organisms',
                                                   'rRNA'))
        for e in expected:
            print(e)
            self.assertTrue(e in [os.path.join(self.wd,p.name,r)
                                  for r in reference_outputs],
                            "%s not found in reference" % e)
        for r in reference_outputs:
            self.assertTrue(os.path.join(self.wd,p.name,r) in expected,
                            "%s not found in expected" % r)

    def test_expected_outputs_10x_multiome_atac_with_libraries_csv(self):
        """
        expected_outputs: 10xGenomics Multiome ATAC with single library analysis
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB1_S1_R3_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           '10xGenomics Single Cell Multiome',
                                           'Library type': 'ATAC' })
        p.create(top_dir=self.wd)
        # Make mock fastq_strand
        mock_fastq_strand_conf = os.path.join(self.wd,
                                              p.name,
                                              "fastq_strand.conf")
        with open(mock_fastq_strand_conf,'w') as fp:
            fp.write("")
        # Make mock libraries.csv files
        os.mkdir(os.path.join(self.wd,
                              p.name,
                              "qc"))
        mock_libraries_csv = os.path.join(self.wd,
                                          p.name,
                                          "qc",
                                          "libraries.PJB1.csv")
        with open(mock_libraries_csv,'w') as fp:
            fp.write("")
        reference_outputs = ("qc/PJB1_S1_R1_001_fastqc",
                             "qc/PJB1_S1_R1_001_fastqc.html",
                             "qc/PJB1_S1_R1_001_fastqc.zip",
                             "qc/PJB1_S1_R1_001_fastqc",
                             "qc/PJB1_S1_R1_001_fastqc.html",
                             "qc/PJB1_S1_R1_001_fastqc.zip",
                             "qc/PJB1_S1_R1_001_model_organisms_screen.png",
                             "qc/PJB1_S1_R1_001_model_organisms_screen.txt",
                             "qc/PJB1_S1_R1_001_other_organisms_screen.png",
                             "qc/PJB1_S1_R1_001_other_organisms_screen.txt",
                             "qc/PJB1_S1_R1_001_rRNA_screen.png",
                             "qc/PJB1_S1_R1_001_rRNA_screen.txt",
                             "qc/PJB1_S1_R1_001_fastq_strand.txt",
                             "qc/PJB1_S1_R3_001_fastqc",
                             "qc/PJB1_S1_R3_001_fastqc.html",
                             "qc/PJB1_S1_R3_001_fastqc.zip",
                             "qc/PJB1_S1_R3_001_fastqc",
                             "qc/PJB1_S1_R3_001_fastqc.html",
                             "qc/PJB1_S1_R3_001_fastqc.zip",
                             "qc/PJB1_S1_R3_001_model_organisms_screen.png",
                             "qc/PJB1_S1_R3_001_model_organisms_screen.txt",
                             "qc/PJB1_S1_R3_001_other_organisms_screen.png",
                             "qc/PJB1_S1_R3_001_other_organisms_screen.txt",
                             "qc/PJB1_S1_R3_001_rRNA_screen.png",
                             "qc/PJB1_S1_R3_001_rRNA_screen.txt",
                             "qc/cellranger_count/PJB1/outs/summary.csv",
                             "qc/cellranger_count/PJB1/outs/web_summary.html",)
        expected = expected_outputs(AnalysisProject(os.path.join(self.wd,
                                                                 p.name)),
                                    "qc",
                                    fastq_strand_conf=mock_fastq_strand_conf,
                                    cellranger_refdata=
                                    "/data/refdata-cellranger-arc-GRCh38-2020-A",
                                    qc_protocol="10x_Multiome_ATAC",
                                    fastq_screens=('model_organisms',
                                                   'other_organisms',
                                                   'rRNA'))
        for e in expected:
            print(e)
            self.assertTrue(e in [os.path.join(self.wd,p.name,r)
                                  for r in reference_outputs],
                            "%s not found in reference" % e)
        for r in reference_outputs:
            self.assertTrue(os.path.join(self.wd,p.name,r) in expected,
                            "%s not found in expected" % r)

    def test_expected_outputs_10x_multiome_gex_with_libraries_csv(self):
        """
        expected_outputs: 10xGenomics Multiome GEX with single library analysis
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           '10xGenomics Single Cell Multiome',
                                           'Library type': 'GEX' })
        p.create(top_dir=self.wd)
        # Make mock fastq_strand
        mock_fastq_strand_conf = os.path.join(self.wd,
                                              p.name,
                                              "fastq_strand.conf")
        with open(mock_fastq_strand_conf,'w') as fp:
            fp.write("")
        # Make mock libraries.csv files
        os.mkdir(os.path.join(self.wd,
                              p.name,
                              "qc"))
        mock_libraries_csv = os.path.join(self.wd,
                                          p.name,
                                          "qc",
                                          "libraries.PJB1.csv")
        with open(mock_libraries_csv,'w') as fp:
            fp.write("")
        reference_outputs = ("qc/PJB1_S1_R1_001_fastqc",
                             "qc/PJB1_S1_R1_001_fastqc.html",
                             "qc/PJB1_S1_R1_001_fastqc.zip",
                             "qc/PJB1_S1_R1_001_fastqc",
                             "qc/PJB1_S1_R1_001_fastqc.html",
                             "qc/PJB1_S1_R1_001_fastqc.zip",
                             "qc/PJB1_S1_R2_001_fastqc",
                             "qc/PJB1_S1_R2_001_fastqc.html",
                             "qc/PJB1_S1_R2_001_fastqc.zip",
                             "qc/PJB1_S1_R2_001_fastqc",
                             "qc/PJB1_S1_R2_001_fastqc.html",
                             "qc/PJB1_S1_R2_001_fastqc.zip",
                             "qc/PJB1_S1_R2_001_fastq_strand.txt",
                             "qc/PJB1_S1_R2_001_model_organisms_screen.png",
                             "qc/PJB1_S1_R2_001_model_organisms_screen.txt",
                             "qc/PJB1_S1_R2_001_other_organisms_screen.png",
                             "qc/PJB1_S1_R2_001_other_organisms_screen.txt",
                             "qc/PJB1_S1_R2_001_rRNA_screen.png",
                             "qc/PJB1_S1_R2_001_rRNA_screen.txt",
                             "qc/cellranger_count/PJB1/outs/summary.csv",
                             "qc/cellranger_count/PJB1/outs/web_summary.html",)
        expected = expected_outputs(AnalysisProject(p.name,
                                                    os.path.join(self.wd,
                                                                 p.name)),
                                    "qc",
                                    fastq_strand_conf=mock_fastq_strand_conf,
                                    cellranger_refdata=
                                    "/data/refdata-cellranger-arc-GRCh38-2020-A",
                                    qc_protocol="10x_Multiome_GEX",
                                    fastq_screens=('model_organisms',
                                                   'other_organisms',
                                                   'rRNA'))
        for e in expected:
            print(e)
            self.assertTrue(e in [os.path.join(self.wd,p.name,r)
                                  for r in reference_outputs],
                            "%s not found in reference" % e)
        for r in reference_outputs:
            self.assertTrue(os.path.join(self.wd,p.name,r) in expected,
                            "%s not found in expected" % r)

    def test_expected_outputs_10x_cellplex_multiplexing(self):
        """
        expected_outputs: 10xGenomics CellPlex multiplexing with cellranger
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_GEX_S1_R1_001.fastq.gz",
                                       "PJB1_GEX_S1_R2_001.fastq.gz",
                                       "PJB2_MC_S2_R1_001.fastq.gz",
                                       "PJB2_MC_S2_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           "10xGenomics Chromium 3'v3" })
        p.create(top_dir=self.wd)
        # Make mock fastq_strand
        mock_fastq_strand_conf = os.path.join(self.wd,
                                              p.name,
                                              "fastq_strand.conf")
        with open(mock_fastq_strand_conf,'w') as fp:
            fp.write("")
        # Make mock cellranger multi config.csv file
        mock_config_csv = os.path.join(self.wd,
                               "PJB",
                               "10x_multi_config.csv")
        with open(mock_config_csv,'wt') as fp:
            fastq_dir = os.path.join(self.wd,
                                     "PJB",
                                     "fastqs")
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
        reference_outputs = ("qc/PJB1_GEX_S1_R1_001_fastqc",
                             "qc/PJB1_GEX_S1_R1_001_fastqc.html",
                             "qc/PJB1_GEX_S1_R1_001_fastqc.zip",
                             "qc/PJB1_GEX_S1_R2_001_fastqc",
                             "qc/PJB1_GEX_S1_R2_001_fastqc.html",
                             "qc/PJB1_GEX_S1_R2_001_fastqc.zip",
                             "qc/PJB1_GEX_S1_R2_001_model_organisms_screen.png",
                             "qc/PJB1_GEX_S1_R2_001_model_organisms_screen.txt",
                             "qc/PJB1_GEX_S1_R2_001_other_organisms_screen.png",
                             "qc/PJB1_GEX_S1_R2_001_other_organisms_screen.txt",
                             "qc/PJB1_GEX_S1_R2_001_rRNA_screen.png",
                             "qc/PJB1_GEX_S1_R2_001_rRNA_screen.txt",
                             "qc/PJB1_GEX_S1_R2_001_fastq_strand.txt",
                             "qc/PJB2_MC_S2_R1_001_fastqc",
                             "qc/PJB2_MC_S2_R1_001_fastqc.html",
                             "qc/PJB2_MC_S2_R1_001_fastqc.zip",
                             "qc/PJB2_MC_S2_R2_001_fastqc",
                             "qc/PJB2_MC_S2_R2_001_fastqc.html",
                             "qc/PJB2_MC_S2_R2_001_fastqc.zip",
                             "qc/PJB2_MC_S2_R2_001_model_organisms_screen.png",
                             "qc/PJB2_MC_S2_R2_001_model_organisms_screen.txt",
                             "qc/PJB2_MC_S2_R2_001_other_organisms_screen.png",
                             "qc/PJB2_MC_S2_R2_001_other_organisms_screen.txt",
                             "qc/PJB2_MC_S2_R2_001_rRNA_screen.png",
                             "qc/PJB2_MC_S2_R2_001_rRNA_screen.txt",
                             "qc/PJB2_MC_S2_R2_001_fastq_strand.txt",
                             "qc/cellranger_multi/6.0.0/refdata-cellranger-gex-GRCh38-2020-A/outs/per_sample_outs/PBA/metrics_summary.csv",
                             "qc/cellranger_multi/6.0.0/refdata-cellranger-gex-GRCh38-2020-A/outs/per_sample_outs/PBA/web_summary.html",
                             "qc/cellranger_multi/6.0.0/refdata-cellranger-gex-GRCh38-2020-A/outs/per_sample_outs/PBB/metrics_summary.csv",
                             "qc/cellranger_multi/6.0.0/refdata-cellranger-gex-GRCh38-2020-A/outs/per_sample_outs/PBB/web_summary.html",
                             "qc/cellranger_multi/6.0.0/refdata-cellranger-gex-GRCh38-2020-A/outs/multi/multiplexing_analysis/tag_calls_summary.csv",)
        expected = expected_outputs(AnalysisProject(p.name,
                                                    os.path.join(self.wd,
                                                                 p.name)),
                                    "qc",
                                    fastq_strand_conf=mock_fastq_strand_conf,
                                    cellranger_multi_config=mock_config_csv,
                                    cellranger_version="6.0.0",
                                    cellranger_refdata=
                                    "/data/refdata-cellranger-gex-GRCh38-2020-A",
                                    qc_protocol="10x_CellPlex",
                                    fastq_screens=('model_organisms',
                                                   'other_organisms',
                                                   'rRNA'))
        for e in expected:
            print(e)
            self.assertTrue(e in [os.path.join(self.wd,p.name,r)
                                  for r in reference_outputs],
                            "%s not found in reference" % e)
        for r in reference_outputs:
            self.assertTrue(os.path.join(self.wd,p.name,r) in expected,
                            "%s not found in expected" % r)

    def test_expected_outputs_10x_cellplex_multiplexing_no_config_csv(self):
        """
        expected_outputs: 10xGenomics CellPlex multiplexing with cellranger (no config.csv file)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_GEX_S1_R1_001.fastq.gz",
                                       "PJB1_GEX_S1_R2_001.fastq.gz",
                                       "PJB2_MC_S2_R1_001.fastq.gz",
                                       "PJB2_MC_S2_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           "10xGenomics Chromium 3'v3" })
        p.create(top_dir=self.wd)
        # Make mock fastq_strand
        mock_fastq_strand_conf = os.path.join(self.wd,
                                              p.name,
                                              "fastq_strand.conf")
        with open(mock_fastq_strand_conf,'w') as fp:
            fp.write("")
        reference_outputs = ("qc/PJB1_GEX_S1_R1_001_fastqc",
                             "qc/PJB1_GEX_S1_R1_001_fastqc.html",
                             "qc/PJB1_GEX_S1_R1_001_fastqc.zip",
                             "qc/PJB1_GEX_S1_R2_001_fastqc",
                             "qc/PJB1_GEX_S1_R2_001_fastqc.html",
                             "qc/PJB1_GEX_S1_R2_001_fastqc.zip",
                             "qc/PJB1_GEX_S1_R2_001_model_organisms_screen.png",
                             "qc/PJB1_GEX_S1_R2_001_model_organisms_screen.txt",
                             "qc/PJB1_GEX_S1_R2_001_other_organisms_screen.png",
                             "qc/PJB1_GEX_S1_R2_001_other_organisms_screen.txt",
                             "qc/PJB1_GEX_S1_R2_001_rRNA_screen.png",
                             "qc/PJB1_GEX_S1_R2_001_rRNA_screen.txt",
                             "qc/PJB1_GEX_S1_R2_001_fastq_strand.txt",
                             "qc/PJB2_MC_S2_R1_001_fastqc",
                             "qc/PJB2_MC_S2_R1_001_fastqc.html",
                             "qc/PJB2_MC_S2_R1_001_fastqc.zip",
                             "qc/PJB2_MC_S2_R2_001_fastqc",
                             "qc/PJB2_MC_S2_R2_001_fastqc.html",
                             "qc/PJB2_MC_S2_R2_001_fastqc.zip",
                             "qc/PJB2_MC_S2_R2_001_model_organisms_screen.png",
                             "qc/PJB2_MC_S2_R2_001_model_organisms_screen.txt",
                             "qc/PJB2_MC_S2_R2_001_other_organisms_screen.png",
                             "qc/PJB2_MC_S2_R2_001_other_organisms_screen.txt",
                             "qc/PJB2_MC_S2_R2_001_rRNA_screen.png",
                             "qc/PJB2_MC_S2_R2_001_rRNA_screen.txt",
                             "qc/PJB2_MC_S2_R2_001_fastq_strand.txt",)
        expected = expected_outputs(AnalysisProject(p.name,
                                                    os.path.join(self.wd,
                                                                 p.name)),
                                    "qc",
                                    fastq_strand_conf=mock_fastq_strand_conf,
                                    cellranger_multi_config=None,
                                    cellranger_version="6.0.0",
                                    cellranger_refdata=
                                    "/data/refdata-cellranger-gex-GRCh38-2020-A",
                                    qc_protocol="10x_CellPlex",
                                    fastq_screens=('model_organisms',
                                                   'other_organisms',
                                                   'rRNA'))
        for e in expected:
            print(e)
            self.assertTrue(e in [os.path.join(self.wd,p.name,r)
                                  for r in reference_outputs],
                            "%s not found in reference" % e)
        for r in reference_outputs:
            self.assertTrue(os.path.join(self.wd,p.name,r) in expected,
                            "%s not found in expected" % r)
