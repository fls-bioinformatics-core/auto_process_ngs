#######################################################################
# Unit tests for qc/protcols.py
#######################################################################

import unittest
import os
import tempfile
import shutil
from auto_process_ngs.analysis import AnalysisProject
from auto_process_ngs.mock import MockAnalysisProject
from auto_process_ngs.qc.protocols import QCProtocol
from auto_process_ngs.qc.protocols import determine_qc_protocol
from auto_process_ngs.qc.protocols import fetch_protocol_definition
from auto_process_ngs.qc.protocols import parse_protocol_spec

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

class TestQCProtocol(unittest.TestCase):
    """
    Tests for the QCProtocol class
    """
    def test_qcprotocol_null_protocol(self):
        """
        QCProtocol: check 'null' protocol
        """
        p = QCProtocol(name="null",
                       description=None,
                       seq_data_reads=None,
                       index_reads=None,
                       qc_modules=None)
        self.assertEqual(p.name,"null")
        self.assertEqual(p.description,"")
        self.assertEqual(p.reads.seq_data,())
        self.assertEqual(p.reads.index,())
        self.assertEqual(p.reads.qc,())
        self.assertEqual(p.read_numbers.seq_data,())
        self.assertEqual(p.read_numbers.index,())
        self.assertEqual(p.read_numbers.qc,())
        self.assertEqual(p.read_range,{})
        self.assertEqual(p.qc_modules,[])
        self.assertEqual(p.summarise(),"'null' protocol: no reads "
                         "explicitly assigned as biological data; no "
                         "reads explicitly assigned as index data")
        self.assertEqual(repr(p),
                         "null:'':seq_reads=[]:index_reads=[]:"
                         "qc_modules=[]")

    def test_qcprotocol_example_paired_end_protocol(self):
        """
        QCProtocol: check basic paired-end protocol
        """
        p = QCProtocol(name="basicPE",
                       description="Basic paired-end QC",
                       seq_data_reads=['r1','r2'],
                       index_reads=None,
                       qc_modules=("fastqc",
                                   "fastq_screen",
                                   "sequence_lengths"))
        self.assertEqual(p.name,"basicPE")
        self.assertEqual(p.description,"Basic paired-end QC")
        self.assertEqual(p.reads.seq_data,('r1','r2'))
        self.assertEqual(p.reads.index,())
        self.assertEqual(p.reads.qc,('r1','r2'))
        self.assertEqual(p.read_numbers.seq_data,(1,2))
        self.assertEqual(p.read_numbers.index,())
        self.assertEqual(p.read_numbers.qc,(1,2))
        self.assertEqual(p.read_range.r1,None)
        self.assertEqual(p.read_range.r2,None)
        self.assertEqual(p.qc_modules,
                         ["fastq_screen",
                          "fastqc",
                          "sequence_lengths"])
        self.assertEqual(p.summarise(),"'basicPE' protocol: biological "
                         "data in R1 and R2; no reads explicitly assigned "
                         "as index data; mapped metrics generated using "
                         "only biological data reads")
        self.assertEqual(repr(p),
                         "basicPE:'Basic paired-end QC':"
                         "seq_reads=[r1,r2]:index_reads=[]:"
                         "qc_modules=[fastq_screen,fastqc,sequence_lengths]")

    def test_qcprotocol_example_single_cell_protocol(self):
        """
        QCProtocol: check basic single cell protocol
        """
        p = QCProtocol(name="basicSC",
                       description="Basic single cell QC",
                       seq_data_reads=['r2'],
                       index_reads=['r1'],
                       qc_modules=("fastqc",
                                   "fastq_screen",
                                   "sequence_lengths"))
        self.assertEqual(p.name,"basicSC")
        self.assertEqual(p.description,"Basic single cell QC")
        self.assertEqual(p.reads.seq_data,('r2',))
        self.assertEqual(p.reads.index,('r1',))
        self.assertEqual(p.reads.qc,('r1','r2'))
        self.assertEqual(p.read_numbers.seq_data,(2,))
        self.assertEqual(p.read_numbers.index,(1,))
        self.assertEqual(p.read_numbers.qc,(1,2))
        self.assertEqual(p.read_range.r1,None)
        self.assertEqual(p.read_range.r2,None)
        self.assertEqual(p.qc_modules,
                         ["fastq_screen",
                          "fastqc",
                          "sequence_lengths"])
        self.assertEqual(p.summarise(),"'basicSC' protocol: biological "
                         "data in R2 only; index data in R1 only; mapped "
                         "metrics generated using only biological data "
                         "reads")
        self.assertEqual(repr(p),
                         "basicSC:'Basic single cell QC':"
                         "seq_reads=[r2]:index_reads=[r1]:"
                         "qc_modules=[fastq_screen,fastqc,sequence_lengths]")

    def test_qcprotocol_example_paired_end_protocol_with_ranges(self):
        """
        QCProtocol: check basic paired end protocol with ranges
        """
        p = QCProtocol(name="basicPE_with_ranges",
                       description="Basic paired-end QC with ranges",
                       seq_data_reads=['r1','r2:1-50'],
                       index_reads=None,
                       qc_modules=("fastqc",
                                   "fastq_screen",
                                   "sequence_lengths"))
        self.assertEqual(p.name,"basicPE_with_ranges")
        self.assertEqual(p.description,"Basic paired-end QC with ranges")
        self.assertEqual(p.reads.seq_data,('r1','r2',))
        self.assertEqual(p.reads.index,())
        self.assertEqual(p.reads.qc,('r1','r2'))
        self.assertEqual(p.read_numbers.seq_data,(1,2))
        self.assertEqual(p.read_numbers.index,())
        self.assertEqual(p.read_numbers.qc,(1,2))
        self.assertEqual(p.read_range.r1,None)
        self.assertEqual(p.read_range.r2,(1,50))
        self.assertEqual(p.qc_modules,
                         ["fastq_screen",
                          "fastqc",
                          "sequence_lengths"])
        self.assertEqual(p.summarise(),"'basicPE_with_ranges' protocol: "
                         "biological data in R1 and R2 (R2 bases 1 to 50); "
                         "no reads explicitly assigned as index data; mapped "
                         "metrics generated using only subsequences of "
                         "biological data reads")
        self.assertEqual(repr(p),
                         "basicPE_with_ranges:"
                         "'Basic paired-end QC with ranges':"
                         "seq_reads=[r1,r2:1-50]:index_reads=[]:"
                         "qc_modules=[fastq_screen,fastqc,sequence_lengths]")

class TestDetermineQCProtocolFunction(unittest.TestCase):
    """
    Tests for determine_qc_protocol function
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestDetermineQCProtocolFunction')

    def tearDown(self):
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    def test_determine_qc_protocol_standardPE(self):
        """determine_qc_protocol: standard paired-end run
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"))
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "standardPE")

    def test_determine_qc_protocol_standardSE(self):
        """determine_qc_protocol: standard single-end run
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",))
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "standardSE")

    def test_determine_qc_protocol_icell8(self):
        """determine_qc_protocol: single-cell run (ICELL8)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "ICELL8"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "singlecell")

    def test_determine_qc_protocol_10xchromium3(self):
        """determine_qc_protocol: single-cell run (10xGenomics Chromium 3')
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "10xGenomics Chromium 3'"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "singlecell")

    def test_determine_qc_protocol_10xchromium3v2(self):
        """determine_qc_protocol: single-cell run (10xGenomics Chromium 3'v2)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "10xGenomics Chromium 3'v2"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "singlecell")

    def test_determine_qc_protocol_10xchromium3v3(self):
        """determine_qc_protocol: single-cell run (10xGenomics Chromium 3'v3)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "10xGenomics Chromium 3'v3"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "singlecell")

    def test_determine_qc_protocol_10xchromium3_rna_seq(self):
        """determine_qc_protocol: single-cell RNA-seq (10xGenomics Chromium 3')
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "10xGenomics Chromium 3'",
                                          'Library type':
                                          "scRNA-seq"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "10x_scRNAseq")

    def test_determine_qc_protocol_10xchromium3v2_rna_seq(self):
        """determine_qc_protocol: single-cell RNA-seq (10xGenomics Chromium 3'v2)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "10xGenomics Chromium 3'v2",
                                          'Library type':
                                          "scRNA-seq"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "10x_scRNAseq")

    def test_determine_qc_protocol_10xchromium3v3_rna_seq(self):
        """determine_qc_protocol: single-cell RNA-seq (10xGenomics Chromium 3'v3)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "10xGenomics Chromium 3'v3",
                                          'Library type':
                                          "scRNA-seq"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "10x_scRNAseq")

    def test_determine_qc_protocol_10xchromium3_sn_rna_seq(self):
        """determine_qc_protocol: single-nuclei RNA-seq (10xGenomics Chromium 3')
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "10xGenomics Chromium 3'",
                                          'Library type':
                                          "snRNA-seq"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "10x_snRNAseq")

    def test_determine_qc_protocol_10xchromium3v3_sn_rna_seq(self):
        """determine_qc_protocol: single-nuclei RNA-seq (10xGenomics Chromium 3'v3)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "10xGenomics Chromium 3'v3",
                                          'Library type':
                                          "snRNA-seq"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "10x_snRNAseq")

    def test_determine_qc_protocol_10xchromium3_cellplex(self):
        """determine_qc_protocol: cell multiplexing CellPlex (10xGenomics Chromium 3')
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "10xGenomics Chromium 3'",
                                          'Library type':
                                          "CellPlex"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "10x_CellPlex")

    def test_determine_qc_protocol_10xchromium3_cellplex_scrnaseq(self):
        """determine_qc_protocol: cell multiplexing CellPlex scRNA-seq (10xGenomics Chromium 3')
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "10xGenomics Chromium 3'",
                                          'Library type':
                                          "CellPlex scRNA-seq"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "10x_CellPlex")

    def test_determine_qc_protocol_10xchromium3_cellplex_snrnaseq(self):
        """determine_qc_protocol: cell multiplexing CellPlex snRNA-seq (10xGenomics Chromium 3')
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "10xGenomics Chromium 3'",
                                          'Library type':
                                          "CellPlex snRNA-seq"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "10x_CellPlex")

    def test_determine_qc_protocol_10xchromium3v3_cellplex(self):
        """determine_qc_protocol: cell multiplexing CellPlex (10xGenomics Chromium 3'v3)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "10xGenomics Chromium 3'v3",
                                          'Library type':
                                          "CellPlex"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "10x_CellPlex")

    def test_determine_qc_protocol_10xchromium3_flex(self):
        """determine_qc_protocol: fixed RNA profiling (Flex) (10xGenomics Chromium 3')
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "10xGenomics Chromium 3'",
                                          'Library type':
                                          "Flex"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "10x_Flex")

    def test_determine_qc_protocol_10x_single_cell_atac_seq(self):
        """determine_qc_protocol: single-cell ATAC-seq (10xGenomics Single Cell ATAC)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "10xGenomics Single Cell ATAC",
                                          'Library type':
                                          "scATAC-seq"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "10x_scATAC")

    def test_determine_qc_protocol_10x_single_nuclei_atac_seq(self):
        """determine_qc_protocol: single-nuclei ATAC-seq (10xGenomics Single Cell ATAC)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "10xGenomics Single Cell ATAC",
                                          'Library type':
                                          "snATAC-seq"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "10x_scATAC")

    def test_determine_qc_protocol_icell8_atac_seq(self):
        """determine_qc_protocol: single-cell ATAC-seq (ICELL8)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "ICELL8",
                                          'Library type':
                                          "scATAC-seq"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "ICELL8_scATAC")

    def test_determine_qc_protocol_10x_visium(self):
        """determine_qc_protocol: spatial RNA-seq run (10xGenomics Visium)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "10xGenomics Visium",
                                          'Library type':
                                          "scATAC-seq"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "10x_Visium")

    def test_determine_qc_protocol_10x_multiome_atac(self):
        """determine_qc_protocol: single cell multiome ATAC run (10xGenomics Multiome ATAC)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"
                                       "PJB2_S2_R3_001.fastq.gz",
                                       "PJB2_S2_R3_001.fastq.gz",
                                       "PJB1_S1_I1_001.fastq.gz",
                                       "PJB1_S1_I2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "10xGenomics Single Cell Multiome",
                                          'Library type':
                                          "ATAC"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "10x_Multiome_ATAC")

    def test_determine_qc_protocol_10x_multiome_gex(self):
        """determine_qc_protocol: single cell multiome GEX run (10xGenomics Multiome GEX)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz",
                                       "PJB1_S1_I1_001.fastq.gz",
                                       "PJB1_S1_I2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "10xGenomics Single Cell Multiome",
                                          'Library type':
                                          "GEX"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "10x_Multiome_GEX")

    def test_determine_qc_protocol_parse_evercode_sc_rnaseq(self):
        """determine_qc_protocol: single-cell RNA-seq (Parse Evercode)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz",
                                       "PJB1_S1_I1_001.fastq.gz",
                                       "PJB1_S1_I2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "Parse Evercode",
                                          'Library type':
                                          "scRNA-seq"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "ParseEvercode")

    def test_determine_qc_protocol_parse_evercode_sn_rnaseq(self):
        """determine_qc_protocol: single-nuclei RNA-seq (Parse Evercode)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz",
                                       "PJB1_S1_I1_001.fastq.gz",
                                       "PJB1_S1_I2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "Parse Evercode",
                                          'Library type':
                                          "snRNA-seq"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "ParseEvercode")

class TestFetchProtocolDefinition(unittest.TestCase):

    def test_fetch_protocol_definition(self):
        """
        fetch_protocol_definition: check definition is returned
        """
        p = fetch_protocol_definition("standardPE")
        self.assertEqual(p.name,"standardPE")
        self.assertEqual(p.reads.seq_data,('r1','r2'))
        self.assertEqual(p.reads.index,())
        self.assertEqual(p.reads.qc,('r1','r2'))
        self.assertEqual(p.read_numbers.seq_data,(1,2))
        self.assertEqual(p.read_numbers.index,())
        self.assertEqual(p.read_numbers.qc,(1,2))
        self.assertEqual(p.qc_modules,['fastq_screen',
                                       'fastqc',
                                       'picard_insert_size_metrics',
                                       'qualimap_rnaseq',
                                       'rseqc_genebody_coverage',
                                       'sequence_lengths',
                                       'strandedness'])

    def test_fetch_protocol_definition_unknown_protocol(self):
        """
        fetch_protocol_definition: raises KeyError for unknown protocol
        """
        self.assertRaises(KeyError,
                          fetch_protocol_definition,
                          "whazzdis?")

class TestParseProtocolSpec(unittest.TestCase):
    def test_parse_protocol_spec(self):
        """
        parse_protocol_spec: simple protocol specification
        """
        s = "simple_pe_qc:'Simple PE QC protocol':" \
            "seq_reads=[r1,r2]:index_reads=[]:" \
            "qc_modules=[fastq_screen,fastqc]"
        p = parse_protocol_spec(s)
        self.assertEqual(p.name,"simple_pe_qc")
        self.assertEqual(p.description,"Simple PE QC protocol")
        self.assertEqual(p.reads.seq_data,('r1','r2'))
        self.assertEqual(p.reads.index,())
        self.assertEqual(p.reads.qc,('r1','r2'))
        self.assertEqual(p.read_numbers.seq_data,(1,2))
        self.assertEqual(p.read_numbers.index,())
        self.assertEqual(p.read_numbers.qc,(1,2))
        self.assertEqual(p.read_range.r1,None)
        self.assertEqual(p.read_range.r2,None)
        self.assertEqual(p.qc_modules,['fastq_screen','fastqc'])
        self.assertEqual(repr(p),s)

    def test_parse_protocol_spec_10x(self):
        """
        parse_protocol_spec: 10x single cell-style protocol specification
        """
        s = "simple_10x_qc:'Simple 10x single cell QC protocol':" \
            "seq_reads=[r2]:index_reads=[r1]:" \
            "qc_modules=[cellranger_count,fastqc]"
        p = parse_protocol_spec(s)
        self.assertEqual(p.name,"simple_10x_qc")
        self.assertEqual(p.description,"Simple 10x single cell QC protocol")
        self.assertEqual(p.reads.seq_data,('r2',))
        self.assertEqual(p.reads.index,('r1',))
        self.assertEqual(p.reads.qc,('r1','r2'))
        self.assertEqual(p.read_numbers.seq_data,(2,))
        self.assertEqual(p.read_numbers.index,(1,))
        self.assertEqual(p.read_numbers.qc,(1,2))
        self.assertEqual(p.read_range.r1,None)
        self.assertEqual(p.read_range.r2,None)
        self.assertEqual(p.qc_modules,['cellranger_count','fastqc'])
        self.assertEqual(repr(p),s)

    def test_parse_protocol_spec_simple_10x_flex(self):
        """
        parse_protocol_spec: 10x Flex-style protocol specification
        """
        s = "simple_flex_qc:'Simple 10x Flex QC protocol':" \
            "seq_reads=[r2:1-50]:index_reads=[r1]:" \
            "qc_modules=[cellranger_count,fastqc]"
        p = parse_protocol_spec(s)
        self.assertEqual(p.name,"simple_flex_qc")
        self.assertEqual(p.description,"Simple 10x Flex QC protocol")
        self.assertEqual(p.reads.seq_data,('r2',))
        self.assertEqual(p.reads.index,('r1',))
        self.assertEqual(p.reads.qc,('r1','r2'))
        self.assertEqual(p.read_numbers.seq_data,(2,))
        self.assertEqual(p.read_numbers.index,(1,))
        self.assertEqual(p.read_numbers.qc,(1,2))
        self.assertEqual(p.read_range.r1,None)
        self.assertEqual(p.read_range.r2,(1,50))
        self.assertEqual(p.qc_modules,['cellranger_count','fastqc'])
        self.assertEqual(repr(p),s)

    def test_parse_protocol_spec_simple_10x_gex_multiome(self):
        """
        parse_protocol_spec: 10x GEX multiome-style protocol specification
        """
        s = "simple_multiome_gex_qc:'Simple 10x multiome GEX QC protocol':" \
            "seq_reads=[r2]:index_reads=[r1]:" \
            "qc_modules=[cellranger-arc_count," \
            "cellranger_count(chemistry=ARC-v1;library=snRNA-seq;" \
            "cellranger_version=*;cellranger_refdata=*;" \
            "set_cell_count=false;set_metadata=False)]"
        p = parse_protocol_spec(s)
        self.assertEqual(p.name,"simple_multiome_gex_qc")
        self.assertEqual(p.description,"Simple 10x multiome GEX QC protocol")
        self.assertEqual(p.reads.seq_data,('r2',))
        self.assertEqual(p.reads.index,('r1',))
        self.assertEqual(p.reads.qc,('r1','r2'))
        self.assertEqual(p.read_numbers.seq_data,(2,))
        self.assertEqual(p.read_numbers.index,(1,))
        self.assertEqual(p.read_numbers.qc,(1,2))
        self.assertEqual(p.read_range.r1,None)
        self.assertEqual(p.read_range.r2,None)
        self.assertEqual(p.qc_modules,[
            'cellranger-arc_count',
            'cellranger_count(chemistry=ARC-v1;library=snRNA-seq;' \
            'cellranger_version=*;cellranger_refdata=*;' \
            'set_cell_count=false;set_metadata=False)'])
        self.assertEqual(repr(p),s)

    def test_parse_protocol_spec_no_index_reads(self):
        """
        parse_protocol_spec: handles missing index reads specification
        """
        s = "simple_pe_qc:'Simple PE QC protocol':" \
            "seq_reads=[r1,r2]:" \
            "qc_modules=[fastq_screen,fastqc]"
        p = parse_protocol_spec(s)
        self.assertEqual(p.name,"simple_pe_qc")
        self.assertEqual(p.description,"Simple PE QC protocol")
        self.assertEqual(p.reads.seq_data,('r1','r2'))
        self.assertEqual(p.reads.index,())
        self.assertEqual(p.reads.qc,('r1','r2'))
        self.assertEqual(p.read_numbers.seq_data,(1,2))
        self.assertEqual(p.read_numbers.index,())
        self.assertEqual(p.read_numbers.qc,(1,2))
        self.assertEqual(p.read_range.r1,None)
        self.assertEqual(p.read_range.r2,None)
        self.assertEqual(p.qc_modules,['fastq_screen','fastqc'])
        self.assertEqual(repr(p),
                         "simple_pe_qc:'Simple PE QC protocol':" \
                         "seq_reads=[r1,r2]:" \
                         "index_reads=[]:" \
                         "qc_modules=[fastq_screen,fastqc]")

    def test_parse_protocol_spec_no_seq_data_reads(self):
        """
        parse_protocol_spec: handles missing seq data reads specification
        """
        s = "simple_pe_qc:'Simple PE QC protocol':" \
            "index_reads=[r1,r2]:" \
            "qc_modules=[fastq_screen,fastqc]"
        p = parse_protocol_spec(s)
        self.assertEqual(p.name,"simple_pe_qc")
        self.assertEqual(p.description,"Simple PE QC protocol")
        self.assertEqual(p.reads.seq_data,())
        self.assertEqual(p.reads.index,('r1','r2'))
        self.assertEqual(p.reads.qc,('r1','r2'))
        self.assertEqual(p.read_numbers.seq_data,())
        self.assertEqual(p.read_numbers.index,(1,2))
        self.assertEqual(p.read_numbers.qc,(1,2))
        self.assertEqual(p.read_range.r1,None)
        self.assertEqual(p.read_range.r2,None)
        self.assertEqual(p.qc_modules,['fastq_screen','fastqc'])
        self.assertEqual(repr(p),
                         "simple_pe_qc:'Simple PE QC protocol':" \
                         "seq_reads=[]:" \
                         "index_reads=[r1,r2]:" \
                         "qc_modules=[fastq_screen,fastqc]")

    def test_parse_protocol_spec_no_qc_modules(self):
        """
        parse_protocol_spec: handles missing QC modules specification
        """
        s = "simple_pe_qc:'Simple PE QC protocol':" \
            "seq_reads=[r1,r2]:" \
            "index_reads=[]"
        p = parse_protocol_spec(s)
        self.assertEqual(p.name,"simple_pe_qc")
        self.assertEqual(p.description,"Simple PE QC protocol")
        self.assertEqual(p.reads.seq_data,('r1','r2'))
        self.assertEqual(p.reads.index,())
        self.assertEqual(p.reads.qc,('r1','r2'))
        self.assertEqual(p.read_numbers.seq_data,(1,2))
        self.assertEqual(p.read_numbers.index,())
        self.assertEqual(p.read_numbers.qc,(1,2))
        self.assertEqual(p.read_range.r1,None)
        self.assertEqual(p.read_range.r2,None)
        self.assertEqual(p.qc_modules,[])
        self.assertEqual(repr(p),
                         "simple_pe_qc:'Simple PE QC protocol':" \
                         "seq_reads=[r1,r2]:" \
                         "index_reads=[]:" \
                         "qc_modules=[]")

    def test_parse_protocol_spec_empty_qc_modules(self):
        """
        parse_protocol_spec: handles trailing colon in specification
        """
        s = "simple_pe_qc:'Simple PE QC protocol':" \
            "seq_reads=[r1,r2]:" \
            "index_reads=[]:" \
            "qc_modules=[fastqc]:"
        p = parse_protocol_spec(s)
        self.assertEqual(p.name,"simple_pe_qc")
        self.assertEqual(p.description,"Simple PE QC protocol")
        self.assertEqual(p.reads.seq_data,('r1','r2'))
        self.assertEqual(p.reads.index,())
        self.assertEqual(p.reads.qc,('r1','r2'))
        self.assertEqual(p.read_numbers.seq_data,(1,2))
        self.assertEqual(p.read_numbers.index,())
        self.assertEqual(p.read_numbers.qc,(1,2))
        self.assertEqual(p.read_range.r1,None)
        self.assertEqual(p.read_range.r2,None)
        self.assertEqual(p.qc_modules,['fastqc'])
        self.assertEqual(repr(p),s[:-1])

    def test_parse_protocol_spec_trailing_colon(self):
        """
        parse_protocol_spec: handles empty QC modules specification
        """
        s = "simple_pe_qc:'Simple PE QC protocol':" \
            "seq_reads=[r1,r2]:" \
            "index_reads=[]:" \
            "qc_modules=[]"
        p = parse_protocol_spec(s)
        self.assertEqual(p.name,"simple_pe_qc")
        self.assertEqual(p.description,"Simple PE QC protocol")
        self.assertEqual(p.reads.seq_data,('r1','r2'))
        self.assertEqual(p.reads.index,())
        self.assertEqual(p.reads.qc,('r1','r2'))
        self.assertEqual(p.read_numbers.seq_data,(1,2))
        self.assertEqual(p.read_numbers.index,())
        self.assertEqual(p.read_numbers.qc,(1,2))
        self.assertEqual(p.read_range.r1,None)
        self.assertEqual(p.read_range.r2,None)
        self.assertEqual(p.qc_modules,[])
        self.assertEqual(repr(p),s)
