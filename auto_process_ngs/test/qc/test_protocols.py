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
from auto_process_ngs.qc.protocols import determine_qc_protocol_from_metadata
from auto_process_ngs.qc.protocols import determine_qc_protocol
from auto_process_ngs.qc.protocols import fetch_protocol_definition
from auto_process_ngs.qc.protocols import parse_protocol_spec
from auto_process_ngs.qc.protocols import parse_qc_module_spec
from auto_process_ngs.qc.protocols import QCProtocolError
from auto_process_ngs.qc.protocols import QCProtocolParseSpecError

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

class TestQCProtocol(unittest.TestCase):
    """
    Tests for the QCProtocol class
    """
    def test_qcprotocol_eq(self):
        """
        QCProtocol: check equality operations
        """
        p = QCProtocol(name="example",
                       description="Example protocol",
                       seq_data_reads=['r1','r2'],
                       index_reads=None,
                       qc_modules=("fastqc",
                                    "fastq_screen",
                                    "sequence_lengths"))
        q = QCProtocol(name="example",
                       description="Example protocol",
                       seq_data_reads=['r1','r3'],
                       index_reads=None,
                       qc_modules=("fastqc",
                                   "fastq_screen",
                                   "sequence_lengths"))
        self.assertEqual(p,p)
        self.assertEqual(q,q)
        self.assertNotEqual(p,q)
        self.assertNotEqual(q,p)
        self.assertNotEqual(p,repr(p))

    def test_qcprotocol_repr(self):
        """
        QCProtocol: check 'repr' output
        """
        p = QCProtocol(name="example",
                       description="Example protocol",
                       seq_data_reads=['r1','r2'],
                       index_reads=None,
                       qc_modules=("fastqc",
                                   "fastq_screen",
                                   "sequence_lengths"))
        self.assertEqual(repr(p),
                         "example:'Example protocol':"
                         "seq_reads=[r1,r2]:"
                         "index_reads=[]:"
                         "qc_modules=[fastq_screen,"
                         "fastqc,sequence_lengths]")

    def test_qcprotocol_from_specification(self):
        """
        QCProtocol: check instantiating using 'from_specification'
        """
        p = QCProtocol(name="example",
                       description="Example protocol",
                       seq_data_reads=['r1','r2'],
                       index_reads=None,
                       qc_modules=("fastqc",
                                   "fastq_screen",
                                   "sequence_lengths"))
        q = QCProtocol.from_specification("example:'Example protocol':"
                                          "seq_reads=[r1,r2]:"
                                          "index_reads=[]:"
                                          "qc_modules=[fastq_screen,"
                                          "fastqc,sequence_lengths]")
        self.assertEqual(p,q)

    def test_qcprotocol_update(self):
        """
        QCProtocol: check protocol can be updated
        """
        p = QCProtocol(name="example",
                       description="Example protocol",
                       seq_data_reads=['r1','r2'],
                       index_reads=None,
                       qc_modules=("fastqc",
                                   "fastq_screen",
                                   "sequence_lengths"))
        self.assertEqual(p.seq_data_reads,['r1','r2'])
        self.assertEqual(p.index_reads,[])
        self.assertEqual(p.read_range,{ 'r1': None, 'r2': None })
        self.assertEqual(p.qc_modules,["fastq_screen",
                                       "fastqc",
                                       "sequence_lengths"])
        p.update(seq_data_reads=['r2'],
                 index_reads=['r1'],
                 qc_modules=["sequence_lengths","fastqc"])
        self.assertEqual(p.seq_data_reads,['r2'])
        self.assertEqual(p.index_reads,['r1'])
        self.assertEqual(p.read_range,{ 'r1': None, 'r2': None })
        self.assertEqual(p.qc_modules,["fastqc",
                                       "sequence_lengths"])

    def test_qcprotocol_handle_wildcard_reads(self):
        """
        QCProtocol: handle wildcard reads ('r*')
        """
        p = QCProtocol(name="example",
                       description="Example protocol",
                       seq_data_reads=['r*'],
                       index_reads=None,
                       qc_modules=("fastqc",
                                   "fastq_screen",
                                   "sequence_lengths"))
        self.assertEqual(p.seq_data_reads,['r*'])
        self.assertEqual(p.index_reads,[])
        self.assertEqual(p.read_range, { 'r*': None})
        self.assertEqual(p.qc_modules,["fastq_screen",
                                       "fastqc",
                                       "sequence_lengths"])

    def test_qcprotocol_unrecognised_module(self):
        """
        QCProtocol: raise exception for unrecognised QC module
        """
        self.assertRaises(QCProtocolError,
                          QCProtocol,
                          name="example",
                          description="Example protocol",
                          seq_data_reads=['r1','r2'],
                          index_reads=None,
                          qc_modules=("fastqc",
                                      "doesnt_exist",
                                      "fastq_screen",
                                      "sequence_lengths"))

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
        self.assertEqual(p.qc_module_names,[])
        self.assertEqual(p.expected_outputs,[])
        self.assertEqual(p.summarise(),"'null' protocol: no reads "
                         "explicitly assigned as biological data; no "
                         "reads explicitly assigned as index data")
        self.assertEqual(repr(p),
                         "null:'':seq_reads=[]:index_reads=[]:"
                         "qc_modules=[]")
        self.assertEqual(p,QCProtocol.from_specification(repr(p)))

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
        self.assertEqual(p.qc_module_names,
                         ["fastq_screen",
                          "fastqc",
                          "sequence_lengths"])
        self.assertEqual(p.expected_outputs,
                         ["fastqc_r1",
                          "fastqc_r2",
                          "screens_r1",
                          "screens_r2",
                          "sequence_lengths"])
        self.assertEqual(p.summarise(),"'basicPE' protocol: biological "
                         "data in R1 and R2; no reads explicitly assigned "
                         "as index data; mapped metrics generated using "
                         "only biological data reads")
        self.assertEqual(repr(p),
                         "basicPE:'Basic paired-end QC':"
                         "seq_reads=[r1,r2]:index_reads=[]:"
                         "qc_modules=[fastq_screen,fastqc,sequence_lengths]")
        self.assertEqual(p,QCProtocol.from_specification(repr(p)))

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
        self.assertEqual(p.qc_module_names,
                         ["fastq_screen",
                          "fastqc",
                          "sequence_lengths"])
        self.assertEqual(p.expected_outputs,
                         ["fastqc_r1",
                          "fastqc_r2",
                          "screens_r2",
                          "sequence_lengths"])
        self.assertEqual(p.summarise(),"'basicSC' protocol: biological "
                         "data in R2 only; index data in R1 only; mapped "
                         "metrics generated using only biological data "
                         "reads")
        self.assertEqual(repr(p),
                         "basicSC:'Basic single cell QC':"
                         "seq_reads=[r2]:index_reads=[r1]:"
                         "qc_modules=[fastq_screen,fastqc,sequence_lengths]")
        self.assertEqual(p,QCProtocol.from_specification(repr(p)))

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
        self.assertEqual(p.qc_module_names,
                         ["fastq_screen",
                          "fastqc",
                          "sequence_lengths"])
        self.assertEqual(p.expected_outputs,
                         ["fastqc_r1",
                          "fastqc_r2",
                          "screens_r1",
                          "screens_r2",
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
        self.assertEqual(p,QCProtocol.from_specification(repr(p)))

    def test_qcprotocol_example_paired_end_protocol_with_parameters(self):
        """
        QCProtocol: check basic paired end protocol with parameters
        """
        p = QCProtocol(name="basicPE_with_parameters",
                       description="Basic paired-end QC with parameters",
                       seq_data_reads=['r1','r2'],
                       index_reads=None,
                       qc_modules=("fastqc(dummy_var=true)",
                                   "fastq_screen(dummy_var=false)",
                                   "sequence_lengths(dummy_var=*)"))
        self.assertEqual(p.name,"basicPE_with_parameters")
        self.assertEqual(p.description,"Basic paired-end QC with parameters")
        self.assertEqual(p.reads.seq_data,('r1','r2',))
        self.assertEqual(p.reads.index,())
        self.assertEqual(p.reads.qc,('r1','r2'))
        self.assertEqual(p.read_numbers.seq_data,(1,2))
        self.assertEqual(p.read_numbers.index,())
        self.assertEqual(p.read_numbers.qc,(1,2))
        self.assertEqual(p.read_range.r1,None)
        self.assertEqual(p.read_range.r2,None)
        self.assertEqual(p.qc_modules,
                         ["fastq_screen(dummy_var=false)",
                          "fastqc(dummy_var=true)",
                          "sequence_lengths(dummy_var=*)"])
        self.assertEqual(p.qc_module_names,
                         ["fastq_screen",
                          "fastqc",
                          "sequence_lengths"])
        self.assertEqual(p.expected_outputs,
                         ["fastqc_r1",
                          "fastqc_r2",
                          "screens_r1",
                          "screens_r2",
                          "sequence_lengths"])
        self.assertEqual(p.summarise(),"'basicPE_with_parameters' protocol: "
                         "biological data in R1 and R2; "
                         "no reads explicitly assigned as index data; mapped "
                         "metrics generated using only biological data reads")
        self.assertEqual(repr(p),
                         "basicPE_with_parameters:"
                         "'Basic paired-end QC with parameters':"
                         "seq_reads=[r1,r2]:index_reads=[]:"
                         "qc_modules=[fastq_screen(dummy_var=false),"
                         "fastqc(dummy_var=true),"
                         "sequence_lengths(dummy_var=*)]")
        self.assertEqual(p,QCProtocol.from_specification(repr(p)))

class TestDetermineQCProtocolFromMetadataFunction(unittest.TestCase):
    """
    Tests for determine_qc_protocol_from_metadata function
    """
    def test_determine_qc_protocol_from_metadata_minimal(self):
        """
        determine_qc_protocol_from_metadata: no library type defined
        """
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type=None,
            single_cell_platform=None,
            paired_end=True),
                         "minimal")

    def test_determine_qc_protocol_from_metadata_standardPE(self):
        """
        determine_qc_protocol_from_metadata: standard paired-end data
        """
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="RNA-seq",
            single_cell_platform=None,
            paired_end=True),
                         "standardPE")

    def test_determine_qc_protocol_from_metadata_standardSE(self):
        """
        determine_qc_protocol_from_metadata: standard single-end data
        """
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="RNA-seq",
            single_cell_platform=None,
            paired_end=False),
                         "standardSE")

    def test_determine_qc_protocol_from_metadata_wgs(self):
        """
        determine_qc_protocol_from_metadata: WGS data
        """
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="WGS",
            single_cell_platform=None,
            paired_end=False),
                         "minimal")

    def test_determine_qc_protocol_from_metadata_dna_seq(self):
        """
        determine_qc_protocol_from_metadata: DNA-seq data
        """
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="DNA-seq",
            single_cell_platform=None,
            paired_end=False),
                         "minimal")

    def test_determine_qc_protocol_from_metadata_crispr(self):
        """
        determine_qc_protocol_from_metadata: CRISPR data
        """
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="CRISPR",
            single_cell_platform=None,
            paired_end=False),
                         "minimal")
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="CRISPR-Cas9",
            single_cell_platform=None,
            paired_end=False),
                         "minimal")

    def test_determine_qc_protocol_from_metadata_icell8(self):
        """
        determine_qc_protocol_from_metadata: ICELL8 data
        """
        # Non-ATAC ICELL8
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="any",
            single_cell_platform="ICELL8",
            paired_end=True),
                         "singlecell")
        # ATAC
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="scATAC-seq",
            single_cell_platform="ICELL8",
            paired_end=True),
                         "ICELL8_scATAC")

    def test_determine_qc_protocol_from_metadata_10xchromium3(self):
        """
        determine_qc_protocol_from_metadata: 10xGenomics Chromium 3' data
        """
        # Default
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="any",
            single_cell_platform="10xGenomics Chromium 3'",
            paired_end=True),
                         "singlecell")
        # scRNA-seq
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="scRNA-seq",
            single_cell_platform="10xGenomics Chromium 3'",
            paired_end=True),
                         "10x_scRNAseq")
        # snRNA-seq
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="snRNA-seq",
            single_cell_platform="10xGenomics Chromium 3'",
            paired_end=True),
                         "10x_snRNAseq")
        # CellPlex
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="CellPlex",
            single_cell_platform="10xGenomics Chromium 3'",
            paired_end=True),
                         "10x_CellPlex")
        # CellPlex scRNA-seq
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="CellPlex scRNA-seq",
            single_cell_platform="10xGenomics Chromium 3'",
            paired_end=True),
                         "10x_CellPlex")
        # CellPlex snRNA-seq
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="CellPlex snRNA-seq",
            single_cell_platform="10xGenomics Chromium 3'",
            paired_end=True),
                         "10x_CellPlex")
        # Flex
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="Flex",
            single_cell_platform="10xGenomics Chromium 3'",
            paired_end=True),
                         "10x_Flex")

    def test_determine_qc_protocol_from_metadata_10xchromium3v2(self):
        """
        determine_qc_protocol_from_metadata: 10xGenomics Chromium 3'v2 data
        """
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="any",
            single_cell_platform="10xGenomics Chromium 3'v2",
            paired_end=True),
                         "singlecell")
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="scRNA-seq",
            single_cell_platform="10xGenomics Chromium 3'v2",
            paired_end=True),
                         "10x_scRNAseq")
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="snRNA-seq",
            single_cell_platform="10xGenomics Chromium 3'v2",
            paired_end=True),
                         "10x_snRNAseq")

    def test_determine_qc_protocol_from_metadata_10xchromium3v3(self):
        """
        determine_qc_protocol_from_metadata: 10xGenomics Chromium 3'v3 data
        """
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="Any",
            single_cell_platform="10xGenomics Chromium 3'v3",
            paired_end=True),
                         "singlecell")
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="scRNA-seq",
            single_cell_platform="10xGenomics Chromium 3'v3",
            paired_end=True),
                         "10x_scRNAseq")
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="snRNA-seq",
            single_cell_platform="10xGenomics Chromium 3'v3",
            paired_end=True),
                         "10x_snRNAseq")
        # CellPlex
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="CellPlex",
            single_cell_platform="10xGenomics Chromium 3'v3",
            paired_end=True),
                         "10x_CellPlex")
        # CellPlex scRNA-seq
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="CellPlex scRNA-seq",
            single_cell_platform="10xGenomics Chromium 3'v3",
            paired_end=True),
                         "10x_CellPlex")
        # CellPlex snRNA-seq
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="CellPlex snRNA-seq",
            single_cell_platform="10xGenomics Chromium 3'v3",
            paired_end=True),
                         "10x_CellPlex")
        # Flex
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="Flex",
            single_cell_platform="10xGenomics Chromium 3'v3",
            paired_end=True),
                         "10x_Flex")

    def test_determine_qc_protocol_from_metadata_10xchromium_gemx3(self):
        """
        determine_qc_protocol_from_metadata: 10xGenomics Chromium GEM-X 3' data
        """
        # scRNA-seq
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="scRNA-seq",
            single_cell_platform="10xGenomics Chromium GEM-X 3'v4",
            paired_end=True),
                         "10x_scRNAseq")

    def test_determine_qc_protocol_from_metadata_10xchromium_gemx(self):
        """
        determine_qc_protocol_from_metadata: 10xGenomics Chromium GEM-X data
        """
        # Flex
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="Flex",
            single_cell_platform="10xGenomics Chromium GEM-X",
            paired_end=True),
                         "10x_Flex")

    def test_determine_qc_protocol_from_metadata_10xchromium_next_gem(self):
        """
        determine_qc_protocol_from_metadata: 10xGenomics Chromium Next GEM data
        """
        # scRNA-seq
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="scRNA-seq",
            single_cell_platform="10xGenomics Chromium Next GEM",
            paired_end=True),
                         "10x_scRNAseq")
        # CellPlex scRNA-seq
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="CellPlex scRNA-seq",
            single_cell_platform="10xGenomics Chromium Next GEM",
            paired_end=True),
                         "10x_CellPlex")
        # Flex
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="Flex",
            single_cell_platform="10xGenomics Chromium Next GEM",
            paired_end=True),
                         "10x_Flex")

    def test_determine_qc_protocol_from_metadata_10xchromium_next_gem3(self):
        """
        determine_qc_protocol_from_metadata: 10xGenomics Chromium Next GEM 3' data
        """
        # scRNA-seq
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="scRNA-seq",
            single_cell_platform="10xGenomics Chromium Next GEM 3'v3.1",
            paired_end=True),
                         "10x_scRNAseq")
        # CellPlex scRNA-seq
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="CellPlex scRNA-seq",
            single_cell_platform="10xGenomics Chromium Next GEM 3'v3.1",
            paired_end=True),
                         "10x_CellPlex")

    def test_determine_qc_protocol_from_metadata_10xchromium5(self):
        """
        determine_qc_protocol_from_metadata: 10xGenomics Chromium 5' data
        """
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="Single Cell Immune Profiling",
            single_cell_platform="10xGenomics Chromium 5'",
            paired_end=True),
                         "10x_ImmuneProfiling")

    def test_determine_qc_protocol_from_metadata_10x_atac(self):
        """
        determine_qc_protocol_from_metadata: 10xGenomics Single Cell ATAC
        """
        # Single cell ATAC
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="scATAC-seq",
            single_cell_platform="10xGenomics Single Cell ATAC",
            paired_end=True),
                         "10x_scATAC")
        # Single nuclei ATAC
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="snATAC-seq",
            single_cell_platform="10xGenomics Single Cell ATAC",
            paired_end=True),
                         "10x_scATAC")

    def test_determine_qc_protocol_from_metadata_10x_visium(self):
        """
        determine_qc_protocol_from_metadata: 10xGenomics Visium
        """
        # Visium Fresh Frozen Spatial Gene Expression
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="Fresh Frozen Spatial Gene Expression",
            single_cell_platform="10xGenomics Visium",
            paired_end=True),
                         "10x_Visium_GEX_90bp_insert")
        # Visium FFPE Spatial Gene Expression
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="FFPE Spatial Gene Expresion",
            single_cell_platform="10xGenomics Visium",
            paired_end=True),
                         "10x_Visium_GEX")
        # Visium (CytAssist) FFPE HD Spatial Gene Expression
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="FFPE spatial Gene Expression",
            single_cell_platform="10xGenomics Visium (CytAssist)",
            paired_end=True),
                         "10x_Visium_GEX")
        # Visium (CytAssist) FFPE Spatial Gene Expression
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="FFPE Spatial Gene Expression",
            single_cell_platform="10xGenomics Visium (CytAssist)",
            paired_end=True),
                         "10x_Visium_GEX")
        # Visium (CytAssist) Fixed Frozen spatial Gene Expression
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="Fixed Frozen Spatial Gene Expression",
            single_cell_platform="10xGenomics Visium (CytAssist)",
            paired_end=True),
                         "10x_Visium_GEX")
        # Visium (CytAssist) Fresh Frozen Spatial Gene Expression
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="Fresh Frozen Spatial Gene Expression",
            single_cell_platform="10xGenomics Visium (CytAssist)",
            paired_end=True),
                         "10x_Visium_GEX")
        # Visium (CytAssist) FFPE Spatial Protein Expression
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="FFPE Spatial Protein Expression",
            single_cell_platform="10xGenomics Visium (CytAssist)",
            paired_end=True),
                         "10x_Visium_PEX")

    def test_determine_qc_protocol_from_metadata_10x_visium_legacy(self):
        """
        determine_qc_protocol_from_metadata: 10xGenomics Visium (legacy)
        """
        # Spatial RNA-seq
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="spatial RNA-seq",
            single_cell_platform="10xGenomics Visium",
            paired_end=True),
                         "10x_Visium_legacy")
        # FFPE spatial RNA-seq
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="FFPE Spatial RNA-seq",
            single_cell_platform="10xGenomics Visium",
            paired_end=True),
                         "10x_Visium_GEX")
        # CytAssist spatial RNA-seq
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="spatial RNA-seq",
            single_cell_platform="10xGenomics CytAssist Visium",
            paired_end=True),
                         "10x_Visium_legacy")
        # CytAssist FFPE RNA-seq
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="FFPE Spatial RNA-seq",
            single_cell_platform="10xGenomics CytAssist Visium",
            paired_end=True),
                         "10x_Visium_GEX")
        # CytAssist FFPE spatial GEX
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="FFPE Spatial GEX",
            single_cell_platform="10xGenomics CytAssist Visium",
            paired_end=True),
                         "10x_Visium_GEX")
        # CytAssist HD spatial GEX
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="HD Spatial GEX",
            single_cell_platform="10xGenomics CytAssist Visium",
            paired_end=True),
                         "10x_Visium_GEX")
        # CytAssist FFPE spatial PEX
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="FFPE Spatial PEX",
            single_cell_platform="10xGenomics CytAssist Visium",
            paired_end=True),
                         "10x_Visium_PEX")

    def test_determine_qc_protocol_from_metadata_10x_multiome(self):
        """
        determine_qc_protocol_from_metadata: 10xGenomics single cell multiome data
        """
        # ATAC component
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="ATAC",
            single_cell_platform="10xGenomics Single Cell Multiome",
            paired_end=True),
                         "10x_Multiome_ATAC")
        # GEX component
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="GEX",
            single_cell_platform="10xGenomics Single Cell Multiome",
            paired_end=True),
                         "10x_Multiome_GEX")

    def test_determine_qc_protocol_from_metadata_parse_evercode(self):
        """
        determine_qc_protocol_from_metadata: Parse Evercode data
        """
        # scRNA-seq
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="scRNA-seq",
            single_cell_platform="Parse Evercode",
            paired_end=True),
                         "ParseEvercode")
        # snRNA-seq
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="snRNA-seq",
            single_cell_platform="Parse Evercode",
            paired_end=True),
                         "ParseEvercode")
        # TCR/TCR scRNA-seq
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="TCR",
            single_cell_platform="Parse Evercode",
            paired_end=True),
                         "ParseEvercode")
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="TCR scRNA-seq",
            single_cell_platform="Parse Evercode",
            paired_end=True),
                         "ParseEvercode")
        # WT/WT scRNA-seq
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="WT",
            single_cell_platform="Parse Evercode",
            paired_end=True),
                         "ParseEvercode")
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="WT scRNA-seq",
            single_cell_platform="Parse Evercode",
            paired_end=True),
                         "ParseEvercode")

    def test_determine_qc_protocol_from_metadata_biorad_ddseq_atac(self):
        """
        determine_qc_protocol_from_metadata: Bio-Rad ddSEQ Single Cell ATAC data
        """
        # scATAC-seq
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="scATAC-seq",
            single_cell_platform="Bio-Rad ddSEQ Single Cell ATAC",
            paired_end=True),
                         "BioRad_ddSEQ_ATAC")
        # snATAC-seq
        self.assertEqual(determine_qc_protocol_from_metadata(
            library_type="snATAC-seq",
            single_cell_platform="Bio-Rad ddSEQ Single Cell ATAC",
            paired_end=True),
                         "BioRad_ddSEQ_ATAC")

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

    def test_determine_qc_protocol_minimal(self):
        """determine_qc_protocol: no library type (minimal QC)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={'Library type':
                                          "RNA-seq"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "standardPE")

    def test_determine_qc_protocol_standardPE(self):
        """determine_qc_protocol: standard paired-end run
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={'Library type':
                                          "RNA-seq"})
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
                                       "PJB2_S2_R1_001.fastq.gz",),
                                metadata={'Library type':
                                          "RNA-seq"})
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

    def test_determine_qc_protocol_10xchromium5_immune_profiling(self):
        """determine_qc_protocol: single cell immune profiling (10xGenomics Chromium 5')
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "10xGenomics Chromium 5'",
                                          'Library type':
                                          "Single Cell Immune Profiling"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "10x_ImmuneProfiling")

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

    def test_determine_qc_protocol_10x_visium_spatial_rnaseq_legacy(self):
        """determine_qc_protocol: 10xGenomics Visium spatial RNA-seq (legacy)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "10xGenomics Visium",
                                          'Library type':
                                          "spatial RNA-seq"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "10x_Visium_legacy")

    def test_determine_qc_protocol_10x_visium_spatial_gex(self):
        """determine_qc_protocol: 10xGenomics Visium spatial GEX
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "10xGenomics Visium",
                                          'Library type':
                                          "spatial RNA-seq"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "10x_Visium_legacy")

    def test_determine_qc_protocol_10x_visium_ffpe_rnaseq(self):
        """determine_qc_protocol: FFPE spatial RNA-seq run (10xGenomics Visium)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "10xGenomics Visium",
                                          'Library type':
                                          "FFPE Spatial RNA-seq"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "10x_Visium_GEX")

    def test_determine_qc_protocol_10x_cytassist_visium(self):
        """determine_qc_protocol: spatial RNA-seq run (10xGenomics CytAssist Visium)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "10xGenomics CytAssist Visium",
                                          'Library type':
                                          "spatial RNA-seq"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "10x_Visium_legacy")

    def test_determine_qc_protocol_10x_visium_cytassist_ffpe_rnaseq(self):
        """determine_qc_protocol: FFPE spatial RNA-seq run (10xGenomics CytAssist Visium)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "10xGenomics CytAssist Visium",
                                          'Library type':
                                          "FFPE Spatial RNA-seq"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "10x_Visium_GEX")

    def test_determine_qc_protocol_10x_visium_cytassist_ffpe_gex(self):
        """determine_qc_protocol: FFPE spatial GEX run (10xGenomics CytAssist Visium)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "10xGenomics Visium (CytAssist)",
                                          'Library type':
                                          "FFPE Spatial GEX"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "10x_Visium_GEX")

    def test_determine_qc_protocol_10x_visium_cytassist_ffpe_pex(self):
        """determine_qc_protocol: FFPE spatial PEX run (10xGenomics CytAssist Visium)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "10xGenomics Visium (CytAssist)",
                                          'Library type':
                                          "FFPE Spatial PEX"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "10x_Visium_PEX")

    def test_determine_qc_protocol_10x_visium_cytassist_ffpe_hd_gex(self):
        """determine_qc_protocol: FFPE HD spatial GEX run (10xGenomics CytAssist Visium)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "10xGenomics Visium (CytAssist)",
                                          'Library type':
                                          "FFPE HD Spatial GEX"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "10x_Visium_GEX")

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

    def test_determine_qc_protocol_parse_evercode_tcr(self):
        """determine_qc_protocol: TCR (Parse Evercode)
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
                                          "TCR"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "ParseEvercode")

    def test_determine_qc_protocol_parse_evercode_tcr_sc_rnaseq(self):
        """determine_qc_protocol: TCR single-cell RNA-seq (Parse Evercode)
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
                                          "TCR scRNA-seq"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "ParseEvercode")

    def test_determine_qc_protocol_parse_evercode_wt(self):
        """determine_qc_protocol: WT (Parse Evercode)
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
                                          "WT"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "ParseEvercode")

    def test_determine_qc_protocol_parse_evercode_wt_sc_rnaseq(self):
        """determine_qc_protocol: WT single-cell RNA-seq (Parse Evercode)
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
                                          "WT scRNA-seq"})
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

    def test_determine_qc_protocol_biorad_ddseq_sc_atac(self):
        """
        determine_qc_protocol: single-cell ATAC-seq (Bio-Rad ddSEQ)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz",
                                       "PJB1_S1_I1_001.fastq.gz",
                                       "PJB1_S1_I2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "Bio-Rad ddSEQ Single Cell ATAC",
                                          'Library type':
                                          "scATAC-seq"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "BioRad_ddSEQ_ATAC")

    def test_determine_qc_protocol_biorad_ddseq_sn_atac(self):
        """
        determine_qc_protocol: single-nuclei ATAC-seq (Bio-Rad ddSEQ)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz",
                                       "PJB1_S1_I1_001.fastq.gz",
                                       "PJB1_S1_I2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "Bio-Rad ddSEQ Single Cell ATAC",
                                          'Library type':
                                          "snATAC-seq"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "BioRad_ddSEQ_ATAC")

    def test_determine_qc_protocol_unknown_single_cell_platform(self):
        """determine_qc_protocol: unknown single cell platform
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "NewSCPlatform",
                                          'Library type':
                                          "scRNA-seq"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "standardPE")

class TestFetchProtocolDefinition(unittest.TestCase):

    def test_fetch_protocol_definition_from_specification(self):
        """
        fetch_protocol_definition: get definition from specification
        """
        p = fetch_protocol_definition("custom:Custom protocol:"
                                      "seq_reads=[r1,r2]:"
                                      "index_reads=[]:"
                                      "qc_modules=[fastqc,"
                                      "fastq_screen]")
        self.assertEqual(p.name,"custom")
        self.assertEqual(p.reads.seq_data,('r1','r2'))
        self.assertEqual(p.reads.index,())
        self.assertEqual(p.reads.qc,('r1','r2'))
        self.assertEqual(p.read_numbers.seq_data,(1,2))
        self.assertEqual(p.read_numbers.index,())
        self.assertEqual(p.read_numbers.qc,(1,2))
        self.assertEqual(p.qc_modules,['fastq_screen',
                                       'fastqc'])
        self.assertEqual(p.qc_module_names,['fastq_screen',
                                            'fastqc'])
        self.assertEqual(p.expected_outputs,['fastqc_r1',
                                             'fastqc_r2',
                                             'screens_r1',
                                             'screens_r2'])

    def test_fetch_protocol_definition_from_specification_wildcard_reads(self):
        """
        fetch_protocol_definition: get definition from specification (wildcard reads)
        """
        p = fetch_protocol_definition("custom:Custom protocol:"
                                      "seq_reads=[r*]:"
                                      "index_reads=[]:"
                                      "qc_modules=[fastqc,"
                                      "fastq_screen]")
        self.assertEqual(p.name,"custom")
        self.assertEqual(p.reads.seq_data,("r*",))
        self.assertEqual(p.reads.index,())
        self.assertEqual(p.reads.qc,('r*',))
        self.assertEqual(p.read_numbers.seq_data,('*',))
        self.assertEqual(p.read_numbers.index,())
        self.assertEqual(p.read_numbers.qc,('*',))
        self.assertEqual(p.qc_modules,['fastq_screen',
                                       'fastqc'])
        self.assertEqual(p.qc_module_names,['fastq_screen',
                                            'fastqc'])
        self.assertEqual(p.expected_outputs,['fastqc_r*',
                                             'screens_r*'])

    def test_fetch_protocol_definition_from_name(self):
        """
        fetch_protocol_definition: get built-in definition from name
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
                                       'rseqc_infer_experiment',
                                       'sequence_lengths'])
        self.assertEqual(p.qc_module_names,['fastq_screen',
                                            'fastqc',
                                            'picard_insert_size_metrics',
                                            'qualimap_rnaseq',
                                            'rseqc_genebody_coverage',
                                            'rseqc_infer_experiment',
                                            'sequence_lengths'])
        self.assertEqual(p.expected_outputs,['fastqc_r1',
                                             'fastqc_r2',
                                             'picard_insert_size_metrics',
                                             'qualimap_rnaseq',
                                             'rseqc_genebody_coverage',
                                             'rseqc_infer_experiment',
                                             'screens_r1',
                                             'screens_r2',
                                             'sequence_lengths'])

    def test_fetch_protocol_definition_unknown_protocol_name(self):
        """
        fetch_protocol_definition: raises KeyError for unknown protocol name
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
        self.assertEqual(p.seq_data_reads,['r1','r2'])
        self.assertEqual(p.index_reads,[])
        self.assertEqual(p.qc_modules,['fastq_screen','fastqc'])

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
        self.assertEqual(p.seq_data_reads,['r2',])
        self.assertEqual(p.index_reads,['r1',])
        self.assertEqual(p.qc_modules,['cellranger_count','fastqc'])

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
        self.assertEqual(p.seq_data_reads,['r2:1-50',])
        self.assertEqual(p.index_reads,['r1',])
        self.assertEqual(p.qc_modules,['cellranger_count','fastqc'])

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
        self.assertEqual(p.seq_data_reads,['r2',])
        self.assertEqual(p.index_reads,['r1',])
        self.assertEqual(p.qc_modules,[
            'cellranger-arc_count',
            'cellranger_count(chemistry=ARC-v1;library=snRNA-seq;' \
            'cellranger_version=*;cellranger_refdata=*;' \
            'set_cell_count=false;set_metadata=False)'])

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
        self.assertEqual(p.seq_data_reads,['r1','r2'])
        self.assertEqual(p.index_reads,None)
        self.assertEqual(p.qc_modules,['fastq_screen','fastqc'])

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
        self.assertEqual(p.seq_data_reads,None)
        self.assertEqual(p.index_reads,['r1','r2'])
        self.assertEqual(p.qc_modules,['fastq_screen','fastqc'])

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
        self.assertEqual(p.seq_data_reads,['r1','r2'])
        self.assertEqual(p.index_reads,[])
        self.assertEqual(p.qc_modules,None)

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
        self.assertEqual(p.seq_data_reads,['r1','r2'])
        self.assertEqual(p.index_reads,[])
        self.assertEqual(p.qc_modules,['fastqc'])

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
        self.assertEqual(p.seq_data_reads,['r1','r2'])
        self.assertEqual(p.index_reads,[])
        self.assertEqual(p.qc_modules,[])

    def test_parse_protocol_spec_incomplete_specification(self):
        """
        parse_protocol_spec: raises exception for incomplete specification
        """
        self.assertRaises(QCProtocolParseSpecError,
                          parse_protocol_spec,
                          "test")
        self.assertRaises(QCProtocolParseSpecError,
                          parse_protocol_spec,
                          "test:description")

    def test_parse_protocol_spec_bad_description(self):
        """
        parse_protocol_spec: raises exception for incorrect description
        """
        self.assertRaises(QCProtocolParseSpecError,
                          parse_protocol_spec,
                          "test:'Continues'after quotes:"
                          "seq_reads=[r1,r2]:"
                          "index_reads=[]:"
                          "qc_modules=[]")
        self.assertRaises(QCProtocolParseSpecError,
                          parse_protocol_spec,
                          "test:'Mismatched quotes\":"
                          "seq_reads=[r1,r2]:"
                          "index_reads=[]:"
                          "qc_modules=[]")

    def test_parse_protocol_spec_bad_seq_data_reads(self):
        """
        parse_protocol_spec: raises exception for incorrect seq data reads
        """
        self.assertRaises(QCProtocolParseSpecError,
                          parse_protocol_spec,
                          "test:Reads beyond closing brace:"
                          "seq_reads=[r1,r2],r3:"
                          "index_reads=[]:"
                          "qc_modules=[]")
        self.assertRaises(QCProtocolParseSpecError,
                          parse_protocol_spec,
                          "test:Closing brace missing:"
                          "seq_reads=[r1,r2:"
                          "index_reads=[]:"
                          "qc_modules=[]")
        self.assertRaises(QCProtocolParseSpecError,
                          parse_protocol_spec,
                          "test:Opening brace missing:"
                          "seq_reads=r1,r2]:"
                          "index_reads=[]:"
                          "qc_modules=[]")
        self.assertRaises(QCProtocolParseSpecError,
                          parse_protocol_spec,
                          "test:No braces:"
                          "seq_reads=r1,r2:"
                          "index_reads=[]:"
                          "qc_modules=[]")

    def test_parse_protocol_spec_bad_index_reads(self):
        """
        parse_protocol_spec: raises exception for incorrect index reads
        """
        self.assertRaises(QCProtocolParseSpecError,
                          parse_protocol_spec,
                          "test:Reads beyond closing brace:"
                          "seq_reads=[]:"
                          "index_reads=[r1,r2],r3:"
                          "qc_modules=[]")
        self.assertRaises(QCProtocolParseSpecError,
                          parse_protocol_spec,
                          "test:Closing brace missing:"
                          "seq_reads=[]:"
                          "index_reads=[r1,r2:"
                          "qc_modules=[]")
        self.assertRaises(QCProtocolParseSpecError,
                          parse_protocol_spec,
                          "test:Opening brace missing:"
                          "seq_reads=[]:"
                          "index_reads=r1,r2]:"
                          "qc_modules=[]")
        self.assertRaises(QCProtocolParseSpecError,
                          parse_protocol_spec,
                          "test:No braces:"
                          "seq_reads=[]:"
                          "index_reads=r1,r2:"
                          "qc_modules=[]")

    def test_parse_protocol_spec_unrecognised_component(self):
        """
        parse_protocol_spec: raises exception for unrecognised component
        """
        self.assertRaises(QCProtocolParseSpecError,
                          parse_protocol_spec,
                          "test:Unrecognised component:" \
                          "seq_reads=[r1,r2]:" \
                          "index_reads=[]:" \
                          "magic_spell=True:" \
                          "qc_modules=[fastqc]")

class TestParseQCModuleSpec(unittest.TestCase):

    def test_parse_qc_module_spec(self):
        """
        parse_qc_module_spec: handle valid specifications
        """
        self.assertEqual(
            parse_qc_module_spec("fastqc"),
            ('fastqc',{}))
        self.assertEqual(
            parse_qc_module_spec("cellranger_count(cellranger_version=6.1.2)"),
            ('cellranger_count',{ 'cellranger_version': '6.1.2' }))
        self.assertEqual(
            parse_qc_module_spec("cellranger_count(cellranger_version=6.1.2;"
                                 "cellranger_refdata=*)"),
            ('cellranger_count',{ 'cellranger_version': '6.1.2',
                                  'cellranger_refdata': '*' }))

    def test_parse_qc_module_spec_quoted_string(self):
        """
        parse_qc_module_spec: handle quoted string values
        """
        self.assertEqual(
            parse_qc_module_spec("module(s=hello)"),
            ('module',{ 's': 'hello' }))
        self.assertEqual(
            parse_qc_module_spec("module(s='hello')"),
            ('module',{ 's': 'hello' }))
        self.assertEqual(
            parse_qc_module_spec("module(s=\"hello\")"),
            ('module',{ 's': 'hello' }))
        self.assertEqual(
            parse_qc_module_spec("module(s='\"hello\"')"),
            ('module',{ 's': '"hello"' }))
        self.assertEqual(
            parse_qc_module_spec("module(s='\"greeting=hello\"')"),
            ('module',{ 's': '"greeting=hello"' }))

    def test_parse_qc_module_spec_boolean(self):
        """
        parse_qc_module_spec: handle boolean values
        """
        self.assertEqual(
            parse_qc_module_spec("module(b=True)"),
            ('module',{ 'b': True }))
        self.assertEqual(
            parse_qc_module_spec("module(b=true)"),
            ('module',{ 'b': True }))
        self.assertEqual(
            parse_qc_module_spec("module(b='true')"),
            ('module',{ 'b': 'true' }))
        self.assertEqual(
            parse_qc_module_spec("module(b=False)"),
            ('module',{ 'b': False }))
        self.assertEqual(
            parse_qc_module_spec("module(b=false)"),
            ('module',{ 'b': False }))
        self.assertEqual(
            parse_qc_module_spec("module(b='false')"),
            ('module',{ 'b': 'false' }))

    def test_parse_qc_module_spec_quoted_special_characters(self):
        """
        parse_qc_module_spec: handle quoted special characters
        """
        self.assertEqual(
            parse_qc_module_spec("module(version='>=9')"),
            ('module',{ 'version': '>=9' }))
        self.assertEqual(
            parse_qc_module_spec("module(version=\"==9\")"),
            ('module',{ 'version': '==9' }))
        self.assertEqual(
            parse_qc_module_spec("module(s='semicolon;')"),
            ('module',{ 's': 'semicolon;' }))
        self.assertEqual(
            parse_qc_module_spec("module(s=\";semicolon\")"),
            ('module',{ 's': ';semicolon' }))
