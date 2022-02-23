#######################################################################
# Unit tests for qc/verification.py
#######################################################################

import unittest
import os
import tempfile
import shutil
from auto_process_ngs.mockqc import MockQCOutputs
from auto_process_ngs.metadata import AnalysisProjectQCDirInfo
from auto_process_ngs.qc.verification import QCVerifier
from auto_process_ngs.qc.verification import parse_qc_module_spec
from auto_process_ngs.qc.verification import filter_fastqs
from auto_process_ngs.qc.verification import filter_10x_pipelines

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

class TestQCVerifier(unittest.TestCase):
    def setUp(self):
        # Temporary working dir (if needed)
        self.wd = None
        # Keep working dir?
        self.remove_test_outputs = True

    def tearDown(self):
        # Remove temporary working dir
        if not (REMOVE_TEST_OUTPUTS and self.remove_test_outputs):
            return
        if self.wd is not None and os.path.isdir(self.wd):
            shutil.rmtree(self.wd)

    def _make_working_dir(self):
        # Create a temporary working directory
        if self.wd is None:
            self.wd = tempfile.mkdtemp(suffix='.test_QCVerifier')

    def _make_qc_dir(self,qc_dir,fastq_names,
                     protocol=None,
                     screens=('model_organisms','other_organisms','rRNA',),
                     cellranger_pipelines=('cellranger',),
                     cellranger_samples=None,
                     cellranger_multi_samples=None,
                     include_fastqc=True,
                     include_fastq_screen=True,
                     include_strandedness=True,
                     include_seqlens=True,
                     include_multiqc=True,
                     include_cellranger_count=False,
                     include_cellranger_multi=False,
                     legacy_screens=False,
                     legacy_cellranger_outs=False):
        # Create working directory and qc dir
        self._make_working_dir()
        qc_dir = os.path.join(self.wd,qc_dir)
        print("QC dir: %s" % qc_dir)
        os.mkdir(qc_dir)
        # QC metadata
        qc_info = AnalysisProjectQCDirInfo()
        qc_info['protocol'] = protocol
        # Populate with fake QC products
        for fq in fastq_names:
            # FastQC
            if include_fastqc:
                MockQCOutputs.fastqc_v0_11_2(fq,qc_dir)
            # Fastq_screen
            if include_fastq_screen:
                for screen in screens:
                    MockQCOutputs.fastq_screen_v0_9_2(
                        fq,qc_dir,screen,legacy=legacy_screens)
                qc_info['fastq_screens'] = ','.join(screens)
            # Strandedness
            if include_strandedness:
                MockQCOutputs.fastq_strand_v0_0_4(fq,qc_dir)
            # Sequence lengths
            if include_seqlens:
                MockQCOutputs.seqlens(fq,qc_dir)
        # Strandedness conf file
        if include_strandedness:
            with open(os.path.join(qc_dir,
                                   "fastq_strand.conf"),'wt') as fp:
                fp.write("Placeholder\n")
        # MultiQC
        if include_multiqc:
            out_file = "multi%s_report.html" % os.path.basename(qc_dir)
            MockQCOutputs.multiqc(self.wd,
                                  multiqc_html=out_file,
                                  version="1.8")
        # Cellranger count
        if include_cellranger_count:
            for cellranger in cellranger_pipelines:
                # Set defaults
                if cellranger == "cellranger":
                    version = "6.1.2"
                    refdata = "/data/refdata-cellranger-2020-A"
                elif cellranger == "cellranger-atac":
                    version = "2.0.0"
                    refdata = "/data/refdata-cellranger-atac-2020-A"
                elif cellranger == "cellranger-arc":
                    version = "2.0.0"
                    refdata = "/data/refdata-cellranger-arc-2020-A"
                # Set top-level output dir
                if not legacy_cellranger_outs:
                    count_dir = os.path.join("cellranger_count",
                                             version,
                                             os.path.basename(refdata))
                else:
                    count_dir = "cellranger_count"
                # Make pipeline outputs
                for sample in cellranger_samples:
                    MockQCOutputs.cellranger_count(
                        sample,
                        qc_dir,
                        cellranger=cellranger,
                        version=version,
                        reference_data_path=refdata,
                        prefix=count_dir)
                    if cellranger == "cellranger-arc":
                        multiome_config = os.path.join(qc_dir,
                                                       "libraries.%s.csv" %
                                                       sample)
                        with open(multiome_config,'wt') as fp:
                            fp.write("Placeholder\n")
                qc_info['cellranger_version'] = version
                qc_info['cellranger_refdata'] = refdata
        # Cellranger multi
        if include_cellranger_multi:
            # Make cellranger multi config.csv file
            multi_config = os.path.join(qc_dir,"10x_multi_config.csv")
            with open(multi_config,'wt') as fp:
                fastq_dir = os.path.join(self.wd,
                                         "PJB",
                                         "fastqs")
                fp.write("""[gene-expression]
reference,/data/refdata-cellranger-2020-A

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB1_GEX,%s,any,PJB1,gene expression,
PJB2_MC,%s,any,PJB2,Multiplexing Capture,

[samples]
sample_id,cmo_ids,description
PJB_CML1,CMO301,CML1
PJB_CML2,CMO302,CML2
""" % (fastq_dir,fastq_dir))
            # Cellranger version
            version = "6.1.2"
            # Set top-level output dir
            multi_dir = os.path.join("cellranger_multi",
                                     version,
                                     "refdata-cellranger-2020-A")
            # Make outputs
            MockQCOutputs.cellranger_multi(cellranger_multi_samples,
                                           qc_dir,
                                           config_csv=multi_config,
                                           prefix=multi_dir)
            qc_info['cellranger_version'] = version
        # Write out metadata file
        qc_info.save(os.path.join(qc_dir,"qc.info"))
        return qc_dir

    def test_qcverifier_verify_single_end(self):
        """
	QCVerifier: verify single-end data (standardSE)
        """
        ##self.remove_test_outputs = False
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',)
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=fastq_names)
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify(fastqs=fastq_names,
                                           qc_protocol="standardSE",
                                           fastq_screens=('model_organisms',
                                                          'other_organisms',
                                                          'rRNA')))

    def test_qcverifier_verify_single_end_no_fastqc(self):
        """
	QCVerifier: verify single-end data (standardSE, no FastQC)
        """
        ##self.remove_test_outputs = False
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',)
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=fastq_names,
                                   include_fastqc=False)
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(qc_verifier.verify(fastqs=fastq_names,
                                            qc_protocol="standardSE",
                                            fastq_screens=('model_organisms',
                                                           'other_organisms',
                                                           'rRNA')))

    def test_qcverifier_verify_single_end_no_screens(self):
        """
	QCVerifier: verify single-end data (standardSE, no screens)
        """
        ##self.remove_test_outputs = False
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',)
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=fastq_names,
                                   include_fastq_screen=False)
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(qc_verifier.verify(fastqs=fastq_names,
                                            qc_protocol="standardSE",
                                            fastq_screens=('model_organisms',
                                                           'other_organisms',
                                                           'rRNA')))

    def test_qcverifier_verify_single_end_no_seqlens(self):
        """
	QCVerifier: verify single-end data (standardSE, no sequence lengths)
        """
        ##self.remove_test_outputs = False
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',)
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=fastq_names,
                                   include_seqlens=False)
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(qc_verifier.verify(fastqs=fastq_names,
                                            qc_protocol="standardSE",
                                            fastq_screens=('model_organisms',
                                                           'other_organisms',
                                                           'rRNA')))

    def test_qcverifier_verify_single_end_no_standedness(self):
        """
	QCVerifier: verify single-end data (standardSE, no strandedness)
        """
        ##self.remove_test_outputs = False
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',)
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=fastq_names,
                                   include_strandedness=False)
        with open(os.path.join(qc_dir,"fastq_strand.conf"),'wt') as fp:
            fp.write("fastq_strand.conf\n")
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(qc_verifier.verify(fastqs=fastq_names,
                                            qc_protocol="standardSE",
                                            fastq_screens=('model_organisms',
                                                           'other_organisms',
                                                           'rRNA')))

    def test_qcverifier_verify_single_end_no_multiqc(self):
        """
	QCVerifier: verify single-end data (standardSE, no MultiQC)
        """
        ##self.remove_test_outputs = False
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',)
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=fastq_names,
                                   include_multiqc=False)
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(qc_verifier.verify(fastqs=fastq_names,
                                            qc_protocol="standardSE",
                                            fastq_screens=('model_organisms',
                                                           'other_organisms',
                                                           'rRNA')))

    def test_qcverifier_verify_single_end_legacy_screen_naming(self):
        """
	QCVerifier: verify single-end data (standardSE, legacy screen naming)
        """
        ##self.remove_test_outputs = False
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',)
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=fastq_names,
                                   legacy_screens=True)
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify(fastqs=fastq_names,
                                           qc_protocol="standardSE",
                                           fastq_screens=('model_organisms',
                                                          'other_organisms',
                                                          'rRNA')))

    def test_qcverifier_verify_paired_end(self):
        """
	QCVerifier: verify paired-end data (standardPE)
        """
        ##self.remove_test_outputs = False
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB1_S1_R2_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',
                     'PJB2_S2_R2_001.fastq.gz',)
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=fastq_names)
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify(fastqs=fastq_names,
                                           qc_protocol="standardPE",
                                           fastq_screens=('model_organisms',
                                                          'other_organisms',
                                                          'rRNA')))

    def test_qcverifier_verify_paired_end_no_fastqc(self):
        """
	QCVerifier: verify paired-end data (standardPE, no FastQC)
        """
        ##self.remove_test_outputs = False
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB1_S1_R2_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',
                     'PJB2_S2_R2_001.fastq.gz',)
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=fastq_names,
                                   include_fastqc=False)
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(qc_verifier.verify(fastqs=fastq_names,
                                            qc_protocol="standardPE",
                                            fastq_screens=('model_organisms',
                                                           'other_organisms',
                                                           'rRNA')))

    def test_qcverifier_verify_paired_end_no_screens(self):
        """
	QCVerifier: verify paired-end data (standardPE, no screens)
        """
        ##self.remove_test_outputs = False
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB1_S1_R2_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',
                     'PJB2_S2_R2_001.fastq.gz',)
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=fastq_names,
                                   include_fastq_screen=False)
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(qc_verifier.verify(fastqs=fastq_names,
                                            qc_protocol="standardPE",
                                            fastq_screens=('model_organisms',
                                                           'other_organisms',
                                                           'rRNA')))

    def test_qcverifier_verify_paired_end_no_seqlens(self):
        """
	QCVerifier: verify paired-end data (standardPE, no sequence lengths)
        """
        ##self.remove_test_outputs = False
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB1_S1_R2_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',
                     'PJB2_S2_R2_001.fastq.gz',)
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=fastq_names,
                                   include_seqlens=False)
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(qc_verifier.verify(fastqs=fastq_names,
                                            qc_protocol="standardPE",
                                            fastq_screens=('model_organisms',
                                                           'other_organisms',
                                                           'rRNA')))

    def test_qcverifier_verify_paired_end_no_strandedness(self):
        """
	QCVerifier: verify paired-end data (standardPE, no strandedness)
        """
        ##self.remove_test_outputs = False
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB1_S1_R2_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',
                     'PJB2_S2_R2_001.fastq.gz',)
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=fastq_names,
                                   include_strandedness=False)
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify(fastqs=fastq_names,
                                           qc_protocol="standardPE",
                                           fastq_screens=('model_organisms',
                                                          'other_organisms',
                                                          'rRNA')))

    def test_qcverifier_verify_paired_end_no_multiqc(self):
        """
	QCVerifier: verify paired-end data (standardPE, no MultiQC)
        """
        ##self.remove_test_outputs = False
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB1_S1_R2_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',
                     'PJB2_S2_R2_001.fastq.gz',)
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=fastq_names,
                                   include_multiqc=False)
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(qc_verifier.verify(fastqs=fastq_names,
                                            qc_protocol="standardPE",
                                            fastq_screens=('model_organisms',
                                                           'other_organisms',
                                                           'rRNA')))

    def test_qcverifier_verify_paired_end_legacy_screen_naming(self):
        """
	QCVerifier: verify paired-end data (standardPE, legacy screen naming)
        """
        ##self.remove_test_outputs = False
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB1_S1_R2_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',
                     'PJB2_S2_R2_001.fastq.gz',)
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=fastq_names,
                                   legacy_screens=True)
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify(fastqs=fastq_names,
                                           qc_protocol="standardPE",
                                           fastq_screens=('model_organisms',
                                                          'other_organisms',
                                                          'rRNA')))

    def test_qcverifier_verify_paired_end_non_canonical_fastq_names(self):
        """
	QCVerifier: verify paired-end data with non-canonical Fastq names
        """
        ##self.remove_test_outputs = False
        fastq_names=('PJB1_S1_R1_001_paired.fastq.gz',
                     'PJB1_S1_R2_001_paired.fastq.gz',
                     'PJB2_S2_R1_001_paired.fastq.gz',
                     'PJB2_S2_R2_001_paired.fastq.gz',)
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=fastq_names)
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify(fastqs=fastq_names,
                                           qc_protocol="standardPE",
                                           fastq_screens=('model_organisms',
                                                          'other_organisms',
                                                          'rRNA')))

    def test_qcverifier_verify_10x_cellranger_count(self):
        """
        QCVerifier: verify 10xGenomics scRNA-seq data (10x_scRNAseq)
        """
        ##self.remove_test_outputs = False
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB1_S1_R2_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',
                     'PJB2_S2_R2_001.fastq.gz',)
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=fastq_names,
                                   include_cellranger_count=True,
                                   cellranger_pipelines=('cellranger',),
                                   cellranger_samples=(
                                       'PJB1',
                                       'PJB2',
                                   ))
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify(fastqs=fastq_names,
                                           qc_protocol="10x_scRNAseq",
                                           fastq_screens=('model_organisms',
                                                          'other_organisms',
                                                          'rRNA'),
                                           cellranger_version="6.1.2",
                                           cellranger_refdata=\
                                           "/data/refdata-cellranger-2020-A"))

    def test_qcverifier_verify_10x_cellranger_count_different_version(self):
        """
        QCVerifier: verify 10xGenomics scRNA-seq data (10x_scRNAseq, different version)
        """
        ##self.remove_test_outputs = False
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB1_S1_R2_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',
                     'PJB2_S2_R2_001.fastq.gz',)
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=fastq_names,
                                   include_cellranger_count=True,
                                   cellranger_pipelines=('cellranger',),
                                   cellranger_samples=(
                                       'PJB1',
                                       'PJB2',
                                   ))
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(qc_verifier.verify(fastqs=fastq_names,
                                            qc_protocol="10x_scRNAseq",
                                            fastq_screens=('model_organisms',
                                                           'other_organisms',
                                                           'rRNA'),
                                            cellranger_version="5.0.0",
                                            cellranger_refdata=\
                                            "/data/refdata-cellranger-2020-A"))

    def test_qcverifier_verify_10x_cellranger_count_legacy(self):
        """
        QCVerifier: verify 10xGenomics scRNA-seq data (10x_scRNAseq, legacy format)
        """
        ##self.remove_test_outputs = False
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB1_S1_R2_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',
                     'PJB2_S2_R2_001.fastq.gz',)
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=fastq_names,
                                   include_cellranger_count=True,
                                   cellranger_pipelines=('cellranger',),
                                   cellranger_samples=(
                                       'PJB1',
                                       'PJB2',
                                   ),
                                   legacy_cellranger_outs=True)
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify(fastqs=fastq_names,
                                           qc_protocol="10x_scRNAseq",
                                           fastq_screens=('model_organisms',
                                                          'other_organisms',
                                                          'rRNA'),
                                           cellranger_version=None,
                                           cellranger_refdata=None))

    def test_qcverifier_10x_cellranger_atac_count(self):
        """
        QCVerifier: verify 10xGenomics scATAC-seq data (10x_scATAC)
        """
        ##self.remove_test_outputs = False
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB1_S1_R3_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',
                     'PJB2_S2_R3_001.fastq.gz',)
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=fastq_names,
                                   include_cellranger_count=True,
                                   cellranger_pipelines=(
                                       'cellranger-atac',
                                   ),
                                   cellranger_samples=(
                                       'PJB1',
                                       'PJB2',
                                   ))
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify(fastqs=fastq_names,
                                           qc_protocol="10x_scATAC",
                                           fastq_screens=('model_organisms',
                                                          'other_organisms',
                                                          'rRNA'),
                                           cellranger_version="2.0.0",
                                           cellranger_refdata=\
                                           "/data/refdata-cellranger-atac-2020-A"))

    def test_verify_qcverifier_10x_multiome_gex(self):
        """
        QCVerifier: verify 10xGenomics multiome GEX data (10x_Multiome_GEX)
        """
        ##self.remove_test_outputs = False
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB1_S1_R2_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',
                     'PJB2_S2_R2_001.fastq.gz',)
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=fastq_names,
                                   include_cellranger_count=True,
                                   cellranger_pipelines=(
                                       'cellranger',
                                       'cellranger-arc',
                                   ),
                                   cellranger_samples=(
                                       'PJB1',
                                       'PJB2',
                                   ))
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify(fastqs=fastq_names,
                                           qc_protocol="10x_Multiome_GEX",
                                           fastq_screens=('model_organisms',
                                                          'other_organisms',
                                                          'rRNA'),
                                           cellranger_version="2.0.0",
                                           cellranger_refdata=\
                                           "refdata-cellranger-arc-2020-A"))

    def test_verify_qcverifier_10x_multiome_gex_no_paired_samples(self):
        """
        QCVerifier: verify 10xGenomics multiome GEX data (10x_Multiome_GEX, no paired samples)
        """
        ##self.remove_test_outputs = False
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB1_S1_R2_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',
                     'PJB2_S2_R2_001.fastq.gz',)
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=fastq_names,
                                   include_cellranger_count=True,
                                   cellranger_pipelines=(
                                       'cellranger',
                                   ),
                                   cellranger_samples=(
                                       'PJB1',
                                       'PJB2',
                                   ))
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify(fastqs=fastq_names,
                                           qc_protocol="10x_Multiome_GEX",
                                           fastq_screens=('model_organisms',
                                                          'other_organisms',
                                                          'rRNA'),
                                           cellranger_version="2.0.0",
                                           cellranger_refdata=\
                                           "refdata-cellranger-arc-2020-A"))

    def test_verify_qcverifier_10x_multiome_gex_missing_paired_samples(self):
        """
        QCVerifier: verify 10xGenomics multiome GEX data (10x_Multiome_GEX, missing paired samples)
        """
        ##self.remove_test_outputs = False
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB1_S1_R2_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',
                     'PJB2_S2_R2_001.fastq.gz',)
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=fastq_names,
                                   include_cellranger_count=True,
                                   cellranger_pipelines=(
                                       'cellranger',
                                   ),
                                   cellranger_samples=(
                                       'PJB1',
                                       'PJB2',
                                   ))
        for s in ('PJB1','PJB2'):
            with open(os.path.join(qc_dir,"libraries.%s.csv" % s),'wt') as fp:
                fp.write("Placeholder\n")
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(qc_verifier.verify(fastqs=fastq_names,
                                            qc_protocol="10x_Multiome_GEX",
                                            fastq_screens=('model_organisms',
                                                           'other_organisms',
                                                           'rRNA'),
                                            cellranger_version="2.0.0",
                                            cellranger_refdata=\
                                            "refdata-cellranger-arc-2020-A"))

    def test_qcverifier_10x_multiome_atac(self):
        """
        QCVerifier: verify 10xGenomics multiome ATAC data (10x_Multiome_ATAC)
        """
        ##self.remove_test_outputs = False
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB1_S1_R3_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',
                     'PJB2_S2_R3_001.fastq.gz',)
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=fastq_names,
                                   include_cellranger_count=True,
                                   cellranger_pipelines=(
                                       'cellranger-arc',
                                       'cellranger-atac',
                                   ),
                                   cellranger_samples=(
                                       'PJB1',
                                       'PJB2',
                                   ))
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify(fastqs=fastq_names,
                                           qc_protocol="10x_Multiome_ATAC",
                                           fastq_screens=('model_organisms',
                                                          'other_organisms',
                                                          'rRNA'),
                                           cellranger_version="2.0.0",
                                           cellranger_refdata=\
                                           "refdata-cellranger-arc-2020-A"))

    def test_qcverifier_10x_multiome_atac_no_paired_samples(self):
        """
        QCVerifier: verify 10xGenomics multiome ATAC data (10x_Multiome_ATAC, no paired samples)
        """
        ##self.remove_test_outputs = False
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB1_S1_R3_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',
                     'PJB2_S2_R3_001.fastq.gz',)
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=fastq_names,
                                   include_cellranger_count=True,
                                   cellranger_pipelines=(
                                       'cellranger-atac',
                                   ),
                                   cellranger_samples=(
                                       'PJB1',
                                       'PJB2',
                                   ))
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify(fastqs=fastq_names,
                                           qc_protocol="10x_Multiome_ATAC",
                                           fastq_screens=('model_organisms',
                                                          'other_organisms',
                                                          'rRNA'),
                                           cellranger_version="2.0.0",
                                           cellranger_refdata=\
                                           "refdata-cellranger-arc-2020-A"))

    def test_qcverifier_10x_multiome_atac_missing_paired_samples(self):
        """
        QCVerifier: verify 10xGenomics multiome ATAC data (10x_Multiome_ATAC, missing paired samples)
        """
        ##self.remove_test_outputs = False
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB1_S1_R3_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',
                     'PJB2_S2_R3_001.fastq.gz',)
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=fastq_names,
                                   include_cellranger_count=True,
                                   cellranger_pipelines=(
                                       'cellranger-atac',
                                   ),
                                   cellranger_samples=(
                                       'PJB1',
                                       'PJB2',
                                   ))
        for s in ('PJB1','PJB2'):
            with open(os.path.join(qc_dir,"libraries.%s.csv" % s),'wt') as fp:
                fp.write("Placeholder\n")
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(qc_verifier.verify(fastqs=fastq_names,
                                            qc_protocol="10x_Multiome_ATAC",
                                            fastq_screens=('model_organisms',
                                                           'other_organisms',
                                                           'rRNA'),
                                            cellranger_version="2.0.0",
                                            cellranger_refdata=\
                                            "refdata-cellranger-arc-2020-A"))

    def test_qcverifier_verify_10x_visium(self):
        """
        QCVerifier: verify 10xGenomics Visium data (10x_scRNAseq)
        """
        ##self.remove_test_outputs = False
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB1_S1_R2_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',
                     'PJB2_S2_R2_001.fastq.gz',)
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=fastq_names)
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify(fastqs=fastq_names,
                                           qc_protocol="10x_Visium",
                                           fastq_screens=('model_organisms',
                                                          'other_organisms',
                                                          'rRNA')))

    def test_qcverifier_verify_10x_cellranger_multi(self):
        """
        QCVerifier: verify 10xGenomics CellPlex data (10x_CellPlex)
        """
        ##self.remove_test_outputs = False
        fastq_names=('PJB1_GEX_S1_R1_001.fastq.gz',
                     'PJB1_GEX_R2_001.fastq.gz',
                     'PJB2_MC_R1_001.fastq.gz',
                     'PJB2_MC_R2_001.fastq.gz',)
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=fastq_names,
                                   include_cellranger_count=True,
                                   include_cellranger_multi=True,
                                   cellranger_pipelines=('cellranger',),
                                   cellranger_samples=(
                                       'PJB1_GEX',
                                       'PJB2_MC',
                                   ),
                                   cellranger_multi_samples=(
                                       'PJB_CML1',
                                       'PJB_CML2',
                                   ))
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify(fastqs=fastq_names,
                                           qc_protocol="10x_CellPlex",
                                           fastq_screens=('model_organisms',
                                                          'other_organisms',
                                                          'rRNA'),
                                           cellranger_version="6.1.2",
                                           cellranger_refdata=\
                                           "/data/refdata-cellranger-2020-A"))

    def test_qcverifier_verify_10x_cellranger_multi_no_config(self):
        """
        QCVerifier: verify 10xGenomics CellPlex data (10x_CellPlex, no multi config file)
        """
        ##self.remove_test_outputs = False
        fastq_names=('PJB1_GEX_S1_R1_001.fastq.gz',
                     'PJB1_GEX_R2_001.fastq.gz',
                     'PJB2_MC_R1_001.fastq.gz',
                     'PJB2_MC_R2_001.fastq.gz',)
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=fastq_names,
                                   include_cellranger_count=True,
                                   cellranger_pipelines=('cellranger',),
                                   cellranger_samples=(
                                       'PJB1_GEX',
                                       'PJB2_MC',
                                   ))
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify(fastqs=fastq_names,
                                           qc_protocol="10x_CellPlex",
                                           fastq_screens=('model_organisms',
                                                          'other_organisms',
                                                          'rRNA'),
                                           cellranger_version="6.1.2",
                                           cellranger_refdata=\
                                           "/data/refdata-cellranger-2020-A"))

    def test_qcverifier_verify_10x_cellranger_multi(self):
        """
        QCVerifier: verify 10xGenomics CellPlex data (10x_CellPlex)
        """
        ##self.remove_test_outputs = False
        fastq_names=('PJB1_GEX_S1_R1_001.fastq.gz',
                     'PJB1_GEX_R2_001.fastq.gz',
                     'PJB2_MC_R1_001.fastq.gz',
                     'PJB2_MC_R2_001.fastq.gz',)
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=fastq_names,
                                   include_cellranger_count=True,
                                   include_cellranger_multi=True,
                                   cellranger_pipelines=('cellranger',),
                                   cellranger_samples=(
                                       'PJB1_GEX',
                                       'PJB2_MC',
                                   ),
                                   cellranger_multi_samples=(
                                       'PJB_CML1',
                                       'PJB_CML2',
                                   ))
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify(fastqs=fastq_names,
                                           qc_protocol="10x_CellPlex",
                                           fastq_screens=('model_organisms',
                                                          'other_organisms',
                                                          'rRNA'),
                                           cellranger_version="6.1.2",
                                           cellranger_refdata=\
                                           "/data/refdata-cellranger-2020-A"))

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

class TestFilterFastqs(unittest.TestCase):

    def test_filter_fastqs(self):
        """
        filter_fastqs: check Fastq names are correctly filtered
        """
        fastqs = ("PJB1_S1_R1_001.fastq.gz",
                  "PJB1_S1_R2_001.fastq.gz",
                  "PJB1_S1_R3_001.fastq.gz",
                  "PJB1_S1_I1_001.fastq.gz",
                  "PJB2_S2_R1_001.fastq.gz",
                  "PJB2_S2_R2_001.fastq.gz",
                  "PJB2_S2_R3_001.fastq.gz",
                  "PJB2_S2_I1_001.fastq.gz",)
        # Filter R1
        self.assertEqual(filter_fastqs(['r1'],fastqs),
                         ["PJB1_S1_R1_001",
                          "PJB2_S2_R1_001",])
        # Filter R1 & R3
        self.assertEqual(filter_fastqs(['r1','r3'],fastqs),
                         ["PJB1_S1_R1_001",
                          "PJB1_S1_R3_001",
                          "PJB2_S2_R1_001",
                          "PJB2_S2_R3_001",])
        # Filter I1
        self.assertEqual(filter_fastqs(['i1'],fastqs),
                         ["PJB1_S1_I1_001",
                          "PJB2_S2_I1_001",])
        # Filter R*
        self.assertEqual(filter_fastqs(['r*'],fastqs),
                         ["PJB1_S1_R1_001",
                          "PJB1_S1_R2_001",
                          "PJB1_S1_R3_001",
                          "PJB2_S2_R1_001",
                          "PJB2_S2_R2_001",
                          "PJB2_S2_R3_001",])
        # Filter *
        self.assertEqual(filter_fastqs(['*'],fastqs),
                         ["PJB1_S1_I1_001",
                          "PJB1_S1_R1_001",
                          "PJB1_S1_R2_001",
                          "PJB1_S1_R3_001",
                          "PJB2_S2_I1_001",
                          "PJB2_S2_R1_001",
                          "PJB2_S2_R2_001",
                          "PJB2_S2_R3_001",])
        # Filter everything
        self.assertEqual(filter_fastqs([],fastqs),[])

class TestFilter10xPipelines(unittest.TestCase):

    def test_filter_10x_pipelines(self):
        """
        filter_10x_pipelines: check pipelines are correctly filtered
        """
        pipelines = (("cellranger","5.0.0","refdata-gex-GRCh38-2.0.0"),
                     ("cellranger","6.0.1","refdata-gex-GRCh38-2020"),
                     ("cellranger","6.1.2","refdata-gex-GRCh38-2020"),
                     ("cellranger-arc","2.0.0","refdata-arc-GRCh38-2020"),
                     ("cellranger-arc","2.0.0","refdata-arc-mm10-2020"),)
        # Specific pipeline
        self.assertEqual(
            filter_10x_pipelines(
                ("cellranger","6.1.2","refdata-gex-GRCh38-2020"),
                pipelines),
            [("cellranger","6.1.2","refdata-gex-GRCh38-2020"),])
        # Specific package and reference data, any version
        self.assertEqual(
            filter_10x_pipelines(
                ("cellranger","*","refdata-gex-GRCh38-2020"),
                pipelines),
            [("cellranger","6.0.1","refdata-gex-GRCh38-2020"),
             ("cellranger","6.1.2","refdata-gex-GRCh38-2020"),])
        # Specific package and version, any reference
        self.assertEqual(
            filter_10x_pipelines(
                ("cellranger-arc","2.0.0","*"),
                pipelines),
            [("cellranger-arc","2.0.0","refdata-arc-GRCh38-2020"),
             ("cellranger-arc","2.0.0","refdata-arc-mm10-2020"),])
        # Specific package, any version, any reference
        self.assertEqual(
            filter_10x_pipelines(
                ("cellranger","*","*"),
                pipelines),
            [("cellranger","5.0.0","refdata-gex-GRCh38-2.0.0"),
             ("cellranger","6.0.1","refdata-gex-GRCh38-2020"),
             ("cellranger","6.1.2","refdata-gex-GRCh38-2020"),])
        # No matching package
        self.assertEqual(
            filter_10x_pipelines(
                ("cellranger-atac","*","*"),
                pipelines),[])
        # No matching version
        self.assertEqual(
            filter_10x_pipelines(
                ("cellranger","3.0.1","*"),
                pipelines),[])
        # No matching reference
        self.assertEqual(
            filter_10x_pipelines(
                ("cellranger","*","refdata-gex-mm10-2020"),
                pipelines),[])
