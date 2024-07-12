#######################################################################
# Unit tests for qc/modules/cellranger_count.py ('cellanger*_count' QC module)
#######################################################################

import unittest
import tempfile
import shutil
import os
from auto_process_ngs.mock import MockAnalysisProject
from auto_process_ngs.mock import UpdateAnalysisProject
from auto_process_ngs.analysis import AnalysisProject
from auto_process_ngs.qc.modules.cellranger_count import check_cellranger_count_outputs
from auto_process_ngs.qc.modules.cellranger_count import check_cellranger_atac_count_outputs
from auto_process_ngs.qc.modules.cellranger_count import check_cellranger_arc_count_outputs

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

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
