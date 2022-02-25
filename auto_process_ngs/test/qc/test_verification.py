#######################################################################
# Unit tests for qc/verification.py
#######################################################################

import unittest
import os
import tempfile
import shutil
from auto_process_ngs.analysis import AnalysisProject
from auto_process_ngs.mock import MockAnalysisProject
from auto_process_ngs.mock import UpdateAnalysisProject
from auto_process_ngs.mockqc import MockQCOutputs
from auto_process_ngs.metadata import AnalysisProjectQCDirInfo
from auto_process_ngs.qc.verification import QCVerifier
from auto_process_ngs.qc.verification import parse_qc_module_spec
from auto_process_ngs.qc.verification import filter_fastqs
from auto_process_ngs.qc.verification import filter_10x_pipelines
from auto_process_ngs.qc.verification import verify_project

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
                     include_multiqc=False,
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

    def test_qcverifier_verify_qc_module_fastqc(self):
        """
        QCVerifier: verify QC module 'fastqc'
        """
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB1_S1_R2_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',
                     'PJB2_S2_R2_001.fastq.gz',)
        # All outputs present
        qc_dir = self._make_qc_dir('qc.ok',
                                   fastq_names=fastq_names,
                                   include_fastqc=True)
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify_qc_module('fastqc',
                                                     fastqs=fastq_names,
                                                     qc_reads=('r1','r2')))
        # Some outputs missing
        qc_dir = self._make_qc_dir('qc.fail',
                                   fastq_names=fastq_names[:-1],
                                   include_fastqc=True)
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(qc_verifier.verify_qc_module('fastqc',
                                                      fastqs=fastq_names,
                                                      qc_reads=('r1','r2')))
        # Empty QC directory
        qc_dir = self._make_qc_dir('qc.empty',
                                   fastq_names=fastq_names,
                                   include_fastqc=False)
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(qc_verifier.verify_qc_module('fastqc',
                                                      fastqs=fastq_names,
                                                      qc_reads=('r1','r2')))

    def test_qcverifier_verify_qc_module_fastq_screen(self):
        """
        QCVerifier: verify QC module 'fastq_screen'
        """
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB1_S1_R2_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',
                     'PJB2_S2_R2_001.fastq.gz',)
        # All outputs present
        qc_dir = self._make_qc_dir('qc.ok',
                                   fastq_names=fastq_names,
                                   include_fastq_screen=True,
                                   screens=('model_organisms',
                                            'other_organisms',
                                            'rRNA',))
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify_qc_module('fastq_screen',
                                                     fastqs=fastq_names,
                                                     fastq_screens=(
                                                         'model_organisms',
                                                         'other_organisms',
                                                         'rRNA',),
                                                     data_reads=('r1','r2')))
        # Some outputs missing
        qc_dir = self._make_qc_dir('qc.fail',
                                   fastq_names=fastq_names[:-1],
                                   include_fastq_screen=True,
                                   screens=('model_organisms',
                                            'other_organisms',
                                            'rRNA',))
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(qc_verifier.verify_qc_module('fastq_screen',
                                                      fastqs=fastq_names,
                                                      fastq_screens=(
                                                          'model_organisms',
                                                          'other_organisms',
                                                          'rRNA',),
                                                      data_reads=('r1','r2')))
        # Some screens missing
        qc_dir = self._make_qc_dir('qc.missing_screens',
                                   fastq_names=fastq_names,
                                   include_fastq_screen=True,
                                   screens=('rRNA',))
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(qc_verifier.verify_qc_module('fastq_screen',
                                                      fastqs=fastq_names,
                                                      fastq_screens=(
                                                          'model_organisms',
                                                          'other_organisms',
                                                          'rRNA',),
                                                      data_reads=('r1','r2')))
        # Empty QC directory
        qc_dir = self._make_qc_dir('qc.empty',
                                   fastq_names=fastq_names,
                                   include_fastq_screen=False)
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(qc_verifier.verify_qc_module('fastq_screen',
                                                      fastqs=fastq_names,
                                                      fastq_screens=(
                                                          'model_organisms',
                                                          'other_organisms',
                                                          'rRNA',),
                                                      data_reads=('r1','r2')))

    def test_qcverifier_verify_qc_module_strandedness(self):
        """
        QCVerifier: verify QC module 'strandedness'
        """
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB1_S1_R2_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',
                     'PJB2_S2_R2_001.fastq.gz',)
        # All outputs present
        qc_dir = self._make_qc_dir('qc.ok',
                                   fastq_names=fastq_names,
                                   include_strandedness=True)
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify_qc_module('strandedness',
                                                     fastqs=fastq_names,
                                                     data_reads=('r1','r2')))
        # Some outputs missing
        qc_dir = self._make_qc_dir('qc.fail',
                                   fastq_names=fastq_names[:2],
                                   include_strandedness=True)
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(qc_verifier.verify_qc_module('strandedness',
                                                      fastqs=fastq_names,
                                                      data_reads=('r1','r2')))
        # Empty QC directory
        # NB this will verify as True because the fastq_strand.conf
        # file is missing (so no outputs are expected)
        qc_dir = self._make_qc_dir('qc.empty',
                                   fastq_names=fastq_names,
                                   include_strandedness=False)
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify_qc_module('strandedness',
                                                     fastqs=fastq_names,
                                                     data_reads=('r1','r2')))

    def test_qcverifier_verify_qc_module_sequence_lengths(self):
        """
        QCVerifier: verify QC module 'sequence_lengths'
        """
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB1_S1_R2_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',
                     'PJB2_S2_R2_001.fastq.gz',)
        # All outputs present
        qc_dir = self._make_qc_dir('qc.ok',
                                   fastq_names=fastq_names,
                                   include_seqlens=True)
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify_qc_module('sequence_lengths',
                                                     fastqs=fastq_names,
                                                     qc_reads=('r1','r2')))
        # Some outputs missing
        qc_dir = self._make_qc_dir('qc.fail',
                                   fastq_names=fastq_names[:-1],
                                   include_seqlens=True)
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(qc_verifier.verify_qc_module('sequence_lengths',
                                                      fastqs=fastq_names,
                                                      qc_reads=('r1','r2')))
        # Empty QC directory
        qc_dir = self._make_qc_dir('qc.empty',
                                   fastq_names=fastq_names,
                                   include_seqlens=False)
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(qc_verifier.verify_qc_module('sequence_lengths',
                                                      fastqs=fastq_names,
                                                      qc_reads=('r1','r2')))

    def test_qcverifier_verify_qc_module_multiqc(self):
        """
        QCVerifier: verify QC module 'multiqc'
        """
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB1_S1_R2_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',
                     'PJB2_S2_R2_001.fastq.gz',)
        # QC dir with MultiQC
        qc_dir = self._make_qc_dir('qc.multiqc',
                                   fastq_names=fastq_names,
                                   include_multiqc=True)
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify_qc_module('multiqc'))
        # QC dir without MultiQC
        qc_dir = self._make_qc_dir('qc.no_multiqc',
                                   fastq_names=fastq_names,
                                   include_multiqc=False)
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(qc_verifier.verify_qc_module('multiqc'))

    def test_qcverifier_verify_qc_module_cellranger_count(self):
        """
        QCVerifier: verify QC module 'cellranger_count'
        """
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB1_S1_R2_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',
                     'PJB2_S2_R2_001.fastq.gz',)
        # QC dir with cellranger count outputs
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=fastq_names,
                                   include_cellranger_count=True,
                                   cellranger_pipelines=('cellranger',),
                                   cellranger_samples=(
                                       'PJB1',
                                       'PJB2',
                                   ))
        qc_verifier = QCVerifier(qc_dir)
        # Implicitly match any version and reference
        self.assertTrue(qc_verifier.verify_qc_module(
            'cellranger_count',
            samples=('PJB1','PJB2')))
        # Explicitly match version with any reference
        self.assertTrue(qc_verifier.verify_qc_module('cellranger_count',
                                                     samples=('PJB1','PJB2'),
                                                     cellranger_version='6.1.2',
                                                     cellranger_refdata='*'))
        # Explicitly match reference with any version
        self.assertTrue(qc_verifier.verify_qc_module('cellranger_count',
                                                     samples=('PJB1','PJB2'),
                                                     cellranger_version='*',
                                                     cellranger_refdata=\
                                                     'refdata-cellranger-2020-A'))
        # Fail if version not found
        self.assertFalse(
            qc_verifier.verify_qc_module('cellranger_count',
                                         samples=('PJB1','PJB2'),
                                         cellranger_version='5.0.0',
                                         cellranger_refdata=\
                                         'refdata-cellranger-2020-A'))
        # Fail if reference not found
        self.assertFalse(
            qc_verifier.verify_qc_module('cellranger_count',
                                         samples=('PJB1','PJB2'),
                                         cellranger_version='6.1.2',
                                         cellranger_refdata=\
                                         'refdata-cellranger-2.0.0'))
        # Missing outputs for one sample
        qc_dir = self._make_qc_dir('qc.fail',
                                   fastq_names=fastq_names[:2],
                                   include_cellranger_count=True,
                                   cellranger_pipelines=('cellranger',),
                                   cellranger_samples=('PJB1',))
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(
            qc_verifier.verify_qc_module('cellranger_count',
                                         samples=('PJB1','PJB2'),
                                         cellranger_version='6.1.2',
                                         cellranger_refdata=\
                                         'refdata-cellranger-2020-A'))
        # Empty QC dir
        qc_dir = self._make_qc_dir('qc.empty',
                                   fastq_names=fastq_names,
                                   include_cellranger_count=False)
        qc_verifier = QCVerifier(qc_dir)
        # Okay if no reference data specified
        self.assertTrue(qc_verifier.verify_qc_module('cellranger_count',
                                                     samples=('PJB1','PJB2')))
        # Fail if reference data is specified
        self.assertFalse(
            qc_verifier.verify_qc_module('cellranger_count',
                                         samples=('PJB1','PJB2'),
                                         cellranger_refdata='*'))

    def test_qcverifier_verify_qc_module_cellranger_atac_count(self):
        """
        QCVerifier: verify QC module 'cellranger-atac_count'
        """
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB1_S1_R2_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',
                     'PJB2_S2_R2_001.fastq.gz',)
        # QC dir with cellranger count outputs
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=fastq_names,
                                   include_cellranger_count=True,
                                   cellranger_pipelines=('cellranger-atac',),
                                   cellranger_samples=(
                                       'PJB1',
                                       'PJB2',
                                   ))
        qc_verifier = QCVerifier(qc_dir)
        # Implicitly match any version and reference
        self.assertTrue(qc_verifier.verify_qc_module(
            'cellranger-atac_count',
            samples=('PJB1','PJB2')))
        # Explicitly match version with any reference
        self.assertTrue(qc_verifier.verify_qc_module('cellranger-atac_count',
                                                     samples=('PJB1','PJB2'),
                                                     cellranger_version='2.0.0',
                                                     cellranger_refdata='*'))
        # Explicitly match reference with any version
        self.assertTrue(
            qc_verifier.verify_qc_module('cellranger-atac_count',
                                         samples=('PJB1','PJB2'),
                                         cellranger_version='*',
                                         cellranger_refdata=\
                                         'refdata-cellranger-atac-2020-A'))
        # Fail if version not found
        self.assertFalse(
            qc_verifier.verify_qc_module('cellranger-atac_count',
                                         samples=('PJB1','PJB2'),
                                         cellranger_version='1.0.0',
                                         cellranger_refdata=\
                                         'refdata-cellranger-atac-2020-A'))
        # Fail if reference not found
        self.assertFalse(
            qc_verifier.verify_qc_module('cellranger-atac_count',
                                         samples=('PJB1','PJB2'),
                                         cellranger_version='2.0.0',
                                         cellranger_refdata=\
                                         'refdata-cellranger-atac-2.0.0'))
        # Missing outputs for one sample
        qc_dir = self._make_qc_dir('qc.fail',
                                   fastq_names=fastq_names[:2],
                                   include_cellranger_count=True,
                                   cellranger_pipelines=('cellranger-atac',),
                                   cellranger_samples=('PJB1',))
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(
            qc_verifier.verify_qc_module('cellranger-atac_count',
                                         samples=('PJB1','PJB2'),
                                         cellranger_version='2.0.0',
                                         cellranger_refdata=\
                                         'refdata-cellranger-atac-2020-A'))
        # Empty QC dir
        qc_dir = self._make_qc_dir('qc.empty',
                                   fastq_names=fastq_names,
                                   include_cellranger_count=False)
        qc_verifier = QCVerifier(qc_dir)
        # Okay if no reference data specified
        self.assertTrue(qc_verifier.verify_qc_module('cellranger-atac_count',
                                                     samples=('PJB1','PJB2')))
        # Fail if reference data is specified
        self.assertFalse(
            qc_verifier.verify_qc_module('cellranger-atac_count',
                                         samples=('PJB1','PJB2'),
                                         cellranger_refdata='*'))

    def test_qcverifier_verify_qc_module_cellranger_arc_count(self):
        """
        QCVerifier: verify QC module 'cellranger-arc_count'
        """
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB1_S1_R2_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',
                     'PJB2_S2_R2_001.fastq.gz',)
        # QC dir with cellranger count outputs
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=fastq_names,
                                   include_cellranger_count=True,
                                   cellranger_pipelines=('cellranger-arc',),
                                   cellranger_samples=(
                                       'PJB1',
                                       'PJB2',
                                   ))
        qc_verifier = QCVerifier(qc_dir)
        # Implicitly match any version and reference
        self.assertTrue(qc_verifier.verify_qc_module(
            'cellranger-arc_count',
            samples=('PJB1','PJB2')))
        # Explicitly match version with any reference
        self.assertTrue(qc_verifier.verify_qc_module('cellranger-arc_count',
                                                     samples=('PJB1','PJB2'),
                                                     cellranger_version='2.0.0',
                                                     cellranger_refdata='*'))
        # Explicitly match reference with any version
        self.assertTrue(
            qc_verifier.verify_qc_module('cellranger-arc_count',
                                         samples=('PJB1','PJB2'),
                                         cellranger_version='*',
                                         cellranger_refdata=\
                                         'refdata-cellranger-arc-2020-A'))
        # Fail if version not found
        self.assertFalse(
            qc_verifier.verify_qc_module('cellranger-arc_count',
                                         samples=('PJB1','PJB2'),
                                         cellranger_version='1.0.0',
                                         cellranger_refdata=\
                                         'refdata-cellranger-arc-2020-A'))
        # Fail if reference not found
        self.assertFalse(
            qc_verifier.verify_qc_module('cellranger-arc_count',
                                         samples=('PJB1','PJB2'),
                                         cellranger_version='2.0.0',
                                         cellranger_refdata=\
                                         'refdata-cellranger-arc-2.0.0'))
        # Missing outputs for one sample
        qc_dir = self._make_qc_dir('qc.fail',
                                   fastq_names=fastq_names[:2],
                                   include_cellranger_count=True,
                                   cellranger_pipelines=('cellranger-arc',),
                                   cellranger_samples=('PJB1',))
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(
            qc_verifier.verify_qc_module('cellranger-arc_count',
                                         samples=('PJB1','PJB2'),
                                         cellranger_version='2.0.0',
                                         cellranger_refdata=\
                                         'refdata-cellranger-arc-2020-A'))
        # Empty QC dir
        # NB this will verify as True because the multiome CSV config
        # files are missing (so no outputs are expected)
        qc_dir = self._make_qc_dir('qc.empty',
                                   fastq_names=fastq_names,
                                   include_cellranger_count=False)
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify_qc_module('cellranger-arc_count',
                                                     samples=('PJB1','PJB2')))

    def test_qcverifier_verify_qc_module_cellranger_multi(self):
        """
        QCVerifier: verify QC module 'cellranger_multi'
        """
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB1_S1_R2_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',
                     'PJB2_S2_R2_001.fastq.gz',)
        # QC dir with cellranger multi outputs
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=fastq_names,
                                   include_cellranger_multi=True,
                                   cellranger_pipelines=('cellranger',),
                                   cellranger_multi_samples=(
                                       'PJB_CML1',
                                       'PJB_CML2',
                                   ))
        qc_verifier = QCVerifier(qc_dir)
        # Implicitly match any version and reference
        self.assertTrue(qc_verifier.verify_qc_module('cellranger_multi'))
        # Explicitly match version with any reference
        self.assertTrue(qc_verifier.verify_qc_module('cellranger_multi',
                                                     cellranger_version='6.1.2',
                                                     cellranger_refdata='*'))
        # Explicitly match reference with any version
        self.assertTrue(qc_verifier.verify_qc_module('cellranger_multi',
                                                     cellranger_version='*',
                                                     cellranger_refdata=\
                                                     'refdata-cellranger-2020-A'))
        # Fail if version not found
        self.assertFalse(
            qc_verifier.verify_qc_module('cellranger_multi',
                                         cellranger_version='5.0.0',
                                         cellranger_refdata=\
                                         'refdata-cellranger-2020-A'))
        # Fail if reference not found
        self.assertFalse(
            qc_verifier.verify_qc_module('cellranger_multi',
                                         cellranger_version='6.1.2',
                                         cellranger_refdata=\
                                         'refdata-cellranger-2.0.0'))
        # Missing outputs for one sample
        qc_dir = self._make_qc_dir('qc.fail',
                                   fastq_names=fastq_names[:2],
                                   include_cellranger_multi=True,
                                   cellranger_pipelines=('cellranger',),
                                   cellranger_multi_samples=('PJB_CML1',))
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(
            qc_verifier.verify_qc_module('cellranger_multi',
                                         samples=('PJB1','PJB2'),
                                         cellranger_version='6.1.2',
                                         cellranger_refdata=\
                                         'refdata-cellranger-2020-A'))
        # Empty QC dir
        # NB this will verify as True because the 10x multi CSV config
        # file is missing (so no outputs are expected)
        qc_dir = self._make_qc_dir('qc.empty',
                                   fastq_names=fastq_names,
                                   include_cellranger_multi=False)
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify_qc_module('cellranger_multi'))

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
                                   # NB only GEX samples
                                   cellranger_samples=(
                                       'PJB1_GEX',
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
                                   include_cellranger_count=False,
                                   include_cellranger_multi=False)
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

class TestVerifyProject(unittest.TestCase):

    def setUp(self):
        # Temporary working dir (if needed)
        self.wd = None

    def tearDown(self):
        # Remove temporary working dir
        if not REMOVE_TEST_OUTPUTS:
            return
        if self.wd is not None and os.path.isdir(self.wd):
            shutil.rmtree(self.wd)

    def _make_working_dir(self):
        # Create a temporary working directory
        if self.wd is None:
            self.wd = tempfile.mkdtemp(suffix='.test_VerifyProject')

    def _make_analysis_project(self,protocol=None,paired_end=True,
                               fastq_dir=None,qc_dir="qc",
                               fastq_names=None,
                               include_multiqc=True,
                               legacy_screens=False):
        # Create a mock Analysis Project directory
        self._make_working_dir()
        # Generate names for fastq files to add
        if paired_end:
            reads = (1,2)
        else:
            reads = (1,)
        sample_names = ('PJB1','PJB2')
        if fastq_names is None:
            fastq_names = []
            for i,sname in enumerate(sample_names,start=1):
                for read in reads:
                    fq = "%s_S%d_R%d_001.fastq.gz" % (sname,i,read)
                    fastq_names.append(fq)
        self.analysis_dir = MockAnalysisProject('PJB',fastq_names)
        # Create the mock directory
        self.analysis_dir.create(top_dir=self.wd)
        # Populate with fake QC products
        project_dir = os.path.join(self.wd,self.analysis_dir.name)
        UpdateAnalysisProject(AnalysisProject(project_dir)).\
            add_qc_outputs(protocol=protocol,
                           qc_dir=qc_dir,
                           include_report=False,
                           include_zip_file=False,
                           include_multiqc=include_multiqc,
                           legacy_screens=legacy_screens)
        return project_dir

    def test_verify_project_single_end(self):
        """
        verify_project: single-end data
        """
        analysis_dir = self._make_analysis_project(protocol="standardSE")
        project = AnalysisProject(analysis_dir)
        self.assertTrue(verify_project(project))

    def test_verify_project_paired_end(self):
        """
        verify_project: paired-end data
        """
        analysis_dir = self._make_analysis_project(protocol="standardPE")
        project = AnalysisProject(analysis_dir)
        self.assertTrue(verify_project(project))

    def test_verify_project_paired_end_with_non_default_fastq_dir(self):
        """
        verify_project: paired-end data with non-default fastq dir
        """
        analysis_dir = self._make_analysis_project(protocol="standardPE",
                                                   fastq_dir=\
                                                   "fastqs.non_default")
        project = AnalysisProject(analysis_dir)
        self.assertTrue(verify_project(project))

    def test_verify_project_paired_end_with_no_fastq_dir(self):
        """
        verify_project: paired-end data with no fastq dir
        """
        analysis_dir = self._make_analysis_project(protocol="standardPE",
                                                   fastq_dir=".")
        project = AnalysisProject(analysis_dir)
        self.assertTrue(verify_project(project))

    def test_verify_project_paired_end_with_non_default_qc_dir(self):
        """
        verify_project: paired-end data with non-default QC dir
        """
        analysis_dir = self._make_analysis_project(protocol="standardPE",
                                                   qc_dir="qc.non_default")
        project = AnalysisProject(analysis_dir)
        self.assertTrue(verify_project(project,qc_dir="qc.non_default"))

    def test_verify_project_paired_end_with_non_canonical_fastq_names(self):
        """
        verify_project: paired-end data with non-canonical fastq names
        """
        analysis_dir = self._make_analysis_project(
            protocol="standardPE",
            fastq_names=
            ("PJB1_S1_R1_001_paired.fastq.gz",
             "PJB1_S1_R2_001_paired.fastq.gz",
             "PJB2_S2_R1_001_paired.fastq.gz",
             "PJB2_S2_R2_001_paired.fastq.gz",))
        project = AnalysisProject(analysis_dir)
        self.assertTrue(verify_project(project))

    def test_verify_project_paired_end_with_no_fastq_dir(self):
        """
        verify_project: paired-end data with legacy screen names
        """
        analysis_dir = self._make_analysis_project(protocol="standardPE",
                                                   legacy_screens=True)
        project = AnalysisProject(analysis_dir)
        self.assertTrue(verify_project(project))
