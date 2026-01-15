#######################################################################
# Unit tests for qc/verification.py
#######################################################################

import unittest
import os
import tempfile
import shutil
from bcftbx.utils import AttributeDictionary
from auto_process_ngs.analysis import AnalysisProject
from auto_process_ngs.mock import make_mock_analysis_project
from auto_process_ngs.mockqc import make_mock_qc_dir
from auto_process_ngs.metadata import AnalysisProjectQCDirInfo
from auto_process_ngs.qc.protocols import fetch_protocol_definition
from auto_process_ngs.qc.verification import QCVerifier
from auto_process_ngs.qc.verification import verify_project
from auto_process_ngs.tenx import DEFAULT_CELLRANGER_VERSION

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
                     organisms=('human',),
                     screens=('model_organisms','other_organisms','rRNA',),
                     cellranger_pipelines=('cellranger',),
                     cellranger_samples=None,
                     cellranger_multi_samples=None,
                     seq_data_samples=None,
                     include_fastqc=True,
                     include_fastq_screen=True,
                     include_strandedness=False,
                     include_seqlens=True,
                     include_picard_insert_size_metrics=False,
                     include_rseqc_genebody_coverage=False,
                     include_rseqc_infer_experiment=False,
                     include_qualimap_rnaseq=False,
                     include_multiqc=False,
                     include_cellranger_count=False,
                     include_cellranger_multi=False,
                     legacy_screens=False,
                     legacy_cellranger_outs=False,
                     legacy_cellranger_count_prefix=False):
        # Create working directory and qc dir
        self._make_working_dir()
        return make_mock_qc_dir(
            os.path.join(self.wd,qc_dir),
            fastq_names,
            protocol=protocol,
            organisms=organisms,
            screens=screens,
            cellranger_pipelines=cellranger_pipelines,
            cellranger_samples=cellranger_samples,
            cellranger_multi_samples=cellranger_multi_samples,
            seq_data_samples=seq_data_samples,
            include_fastqc=include_fastqc,
            include_fastq_screen=include_fastq_screen,
            include_strandedness=include_strandedness,
            include_seqlens=include_seqlens,
            include_picard_insert_size_metrics=\
            include_picard_insert_size_metrics,
            include_rseqc_genebody_coverage=\
            include_rseqc_genebody_coverage,
            include_rseqc_infer_experiment=\
            include_rseqc_infer_experiment,
            include_qualimap_rnaseq=include_qualimap_rnaseq,
            include_multiqc=include_multiqc,
            include_cellranger_count=include_cellranger_count,
            include_cellranger_multi=include_cellranger_multi,
            legacy_screens=legacy_screens,
            legacy_cellranger_outs=legacy_cellranger_outs,
            legacy_cellranger_count_prefix=legacy_cellranger_count_prefix)

    def _create_params_dict(self,**kws):
        # Return an AttributeDictionary with the
        # appropriate parameters set
        params = AttributeDictionary(
            qc_dir=None,
            fastqs=None,
            samples=None,
            seq_data_fastqs=None,
            seq_data_samples=None,
            seq_data_reads=None,
            qc_reads=None,
            organism=None,
            fastq_screens=None,
            star_index=None,
            annotation_bed=None,
            annotation_gtf=None,
            cellranger_version=None,
            cellranger_refdata=None,
            cellranger_use_multi_config=None
        )
        if kws:
            for k in kws:
                params[k] = kws[k]
        if not params.seq_data_fastqs:
            params['seq_data_fastqs'] = params.fastqs
        return params

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
        self.assertTrue(qc_verifier.verify_qc_module(
            'fastqc',
            self._create_params_dict(fastqs=fastq_names,
                                     qc_reads=('r1','r2'))))
        # Some outputs missing
        qc_dir = self._make_qc_dir('qc.fail',
                                   fastq_names=fastq_names[:-1],
                                   include_fastqc=True)
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(qc_verifier.verify_qc_module(
            'fastqc',
            self._create_params_dict(fastqs=fastq_names,
                                     qc_reads=('r1','r2'))))
        # Empty QC directory
        qc_dir = self._make_qc_dir('qc.empty',
                                   fastq_names=fastq_names,
                                   include_fastqc=False)
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(qc_verifier.verify_qc_module(
            'fastqc',
            self._create_params_dict(fastqs=fastq_names,
                                     qc_reads=('r1','r2'))))

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
        self.assertTrue(qc_verifier.verify_qc_module(
            'fastq_screen',
            self._create_params_dict(fastqs=fastq_names,
                                     fastq_screens=(
                                         'model_organisms',
                                         'other_organisms',
                                         'rRNA',),
                                     seq_data_reads=('r1','r2'))))
        # Some outputs missing
        qc_dir = self._make_qc_dir('qc.fail',
                                   fastq_names=fastq_names[:-1],
                                   include_fastq_screen=True,
                                   screens=('model_organisms',
                                            'other_organisms',
                                            'rRNA',))
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(qc_verifier.verify_qc_module(
            'fastq_screen',
            self._create_params_dict(fastqs=fastq_names,
                                     fastq_screens=(
                                         'model_organisms',
                                         'other_organisms',
                                         'rRNA',),
                                     seq_data_reads=('r1','r2'))))
        # Some screens missing
        qc_dir = self._make_qc_dir('qc.missing_screens',
                                   fastq_names=fastq_names,
                                   include_fastq_screen=True,
                                   screens=('rRNA',))
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(qc_verifier.verify_qc_module(
            'fastq_screen',
            self._create_params_dict(fastqs=fastq_names,
                                     fastq_screens=(
                                         'model_organisms',
                                         'other_organisms',
                                         'rRNA',),
                                     seq_data_reads=('r1','r2'))))
        # Empty QC directory
        qc_dir = self._make_qc_dir('qc.empty',
                                   fastq_names=fastq_names,
                                   include_fastq_screen=False)
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(qc_verifier.verify_qc_module(
            'fastq_screen',
            self._create_params_dict(fastqs=fastq_names,
                                     fastq_screens=(
                                         'model_organisms',
                                         'other_organisms',
                                         'rRNA',),
                                     seq_data_reads=('r1','r2'))))

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
        self.assertTrue(qc_verifier.verify_qc_module(
            'strandedness',
            self._create_params_dict(fastqs=fastq_names,
                                     seq_data_reads=('r1','r2'))))
        # Some outputs missing
        qc_dir = self._make_qc_dir('qc.fail',
                                   fastq_names=fastq_names[:2],
                                   include_strandedness=True)
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(qc_verifier.verify_qc_module(
            'strandedness',
            self._create_params_dict(fastqs=fastq_names,
                                     seq_data_reads=('r1','r2'))))
        # Empty QC directory
        # NB this will verify as None because the fastq_strand.conf
        # file is missing (so no outputs are expected)
        qc_dir = self._make_qc_dir('qc.empty',
                                   fastq_names=fastq_names,
                                   include_strandedness=False)
        qc_verifier = QCVerifier(qc_dir)
        self.assertEqual(None,
                         qc_verifier.verify_qc_module(
                             'strandedness',
                             self._create_params_dict(fastqs=fastq_names,
                                                      seq_data_reads=(
                                                          'r1','r2'))))

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
        self.assertTrue(qc_verifier.verify_qc_module(
            'sequence_lengths',
            self._create_params_dict(fastqs=fastq_names,
                                     qc_reads=('r1','r2'))))
        # Some outputs missing
        qc_dir = self._make_qc_dir('qc.fail',
                                   fastq_names=fastq_names[:-1],
                                   include_seqlens=True)
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(qc_verifier.verify_qc_module(
            'sequence_lengths',
            self._create_params_dict(fastqs=fastq_names,
                                     qc_reads=('r1','r2'))))
        # Empty QC directory
        qc_dir = self._make_qc_dir('qc.empty',
                                   fastq_names=fastq_names,
                                   include_seqlens=False)
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(qc_verifier.verify_qc_module(
            'sequence_lengths',
            self._create_params_dict(fastqs=fastq_names,
                                     qc_reads=('r1','r2'))))

    def test_qcverifier_verify_qc_module_rseqc_genebody_coverage(self):
        """
        QCVerifier: verify QC module 'rseqc_genebody_coverage'
        """
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB1_S1_R2_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',
                     'PJB2_S2_R2_001.fastq.gz',)
        # All outputs present
        qc_dir = self._make_qc_dir('qc.ok',
                                   fastq_names=fastq_names,
                                   organisms=('human',),
                                   include_rseqc_genebody_coverage=True)
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify_qc_module(
            'rseqc_genebody_coverage',
            self._create_params_dict(fastqs=fastq_names,
                                     organism="Human",
                                     star_index="/data/indexes/STAR",
                                     annotation_bed="/data/annot/human.bed")))
        # Returns None if organism, STAR index or annotation
        # is not supplied
        self.assertEqual(None,
                         qc_verifier.verify_qc_module(
                             'rseqc_genebody_coverage',
                             self._create_params_dict(
                                 fastqs=fastq_names,
                                 star_index="/data/indexes/STAR",
                                 annotation_bed="/data/annot/human.bed")))
        self.assertEqual(None,
                         qc_verifier.verify_qc_module(
                             'rseqc_genebody_coverage',
                             self._create_params_dict(
                                 fastqs=fastq_names,
                                 organism="Human",
                                 annotation_bed="/data/annot/human.bed")))
        self.assertEqual(None,
                         qc_verifier.verify_qc_module(
                             'rseqc_genebody_coverage',
                             self._create_params_dict(
                                 fastqs=fastq_names,
                                 organism="Human",
                                 star_index="/data/indexes/STAR")))
        # Organism name contains spaces
        qc_dir = self._make_qc_dir('qc.homo_sapiens',
                                   fastq_names=fastq_names[:-2],
                                   organisms=('Homo sapiens',),
                                   include_rseqc_genebody_coverage=True)
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify_qc_module(
            'rseqc_genebody_coverage',
            self._create_params_dict(fastqs=fastq_names,
                                     organism="Homo sapiens",
                                     star_index="/data/indexes/STAR",
                                     annotation_bed="/data/annot/human.bed")))
        # Empty QC directory
        qc_dir = self._make_qc_dir('qc.empty',
                                   fastq_names=fastq_names,
                                   organisms=('human',),
                                   include_rseqc_genebody_coverage=False)
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(qc_verifier.verify_qc_module(
            'rseqc_genebody_coverage',
            self._create_params_dict(fastqs=fastq_names,
                                     organism="Human",
                                     star_index="/data/indexes/STAR",
                                     annotation_bed="/data/annot/human.bed")))

    def test_qcverifier_verify_qc_module_rseqc_infer_experiment(self):
        """
        QCVerifier: verify QC module 'rseqc_infer_experiment'
        """
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB1_S1_R2_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',
                     'PJB2_S2_R2_001.fastq.gz',)
        # All outputs present
        qc_dir = self._make_qc_dir('qc.ok',
                                   fastq_names=fastq_names,
                                   organisms=('human',),
                                   include_rseqc_infer_experiment=True)
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify_qc_module(
            'rseqc_infer_experiment',
            self._create_params_dict(fastqs=fastq_names,
                                     seq_data_reads=('r1','r2'),
                                     organism="Human",
                                     star_index="/data/indexes/STAR",
                                     annotation_gtf="/data/annot/human.gtf",
                                     annotation_bed="/data/annot/human.bed")))
        # Returns None if organism, STAR index or annotation
        # is not supplied
        self.assertEqual(None,
                         qc_verifier.verify_qc_module(
                             'rseqc_infer_experiment',
                             self._create_params_dict(
                                 fastqs=fastq_names,
                                 seq_data_reads=('r1','r2'),
                                 star_index="/data/indexes/STAR",
                                 annotation_gtf="/data/annot/human.gtf",
                                 annotation_bed="/data/annot/human.bed")))
        self.assertEqual(None,
                         qc_verifier.verify_qc_module(
                             'rseqc_infer_experiment',
                             self._create_params_dict(
                                 fastqs=fastq_names,
                                 seq_data_reads=('r1','r2'),
                                 organism="Human",
                                 annotation_gtf="/data/annot/human.gtf",
                                 annotation_bed="/data/annot/human.bed")))
        self.assertEqual(None,
                         qc_verifier.verify_qc_module(
                             'rseqc_infer_experiment',
                             self._create_params_dict(
                                 fastqs=fastq_names,
                                 seq_data_reads=('r1','r2'),
                                 organism="Human",
                                 star_index="/data/indexes/STAR")))
        # Organism name contains spaces
        qc_dir = self._make_qc_dir('qc.homo_sapiens',
                                   fastq_names=fastq_names,
                                   organisms=('Homo sapiens',),
                                   include_rseqc_infer_experiment=True)
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify_qc_module(
            'rseqc_infer_experiment',
            self._create_params_dict(fastqs=fastq_names,
                                     seq_data_reads=('r1','r2'),
                                     organism="Homo sapiens",
                                     star_index="/data/indexes/STAR",
                                     annotation_gtf="/data/annot/human.gtf",
                                     annotation_bed="/data/annot/human.bed")))
        # Empty QC directory
        qc_dir = self._make_qc_dir('qc.empty',
                                   fastq_names=fastq_names,
                                   organisms=('human',),
                                   include_rseqc_infer_experiment=False)
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(qc_verifier.verify_qc_module(
            'rseqc_infer_experiment',
            self._create_params_dict(fastqs=fastq_names,
                                     seq_data_reads=('r1','r2'),
                                     organism="Human",
                                     star_index="/data/indexes/STAR",
                                     annotation_gtf="/data/annot/human.gtf",
                                     annotation_bed="/data/annot/human.bed")))

    def test_qcverifier_verify_qc_module_picard_insert_size_metrics(self):
        """
        QCVerifier: verify QC module 'picard_insert_size_metrics'
        """
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB1_S1_R2_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',
                     'PJB2_S2_R2_001.fastq.gz',)
        # All outputs present
        qc_dir = self._make_qc_dir('qc.ok',
                                   fastq_names=fastq_names,
                                   organisms=('human',),
                                   include_picard_insert_size_metrics=True)
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify_qc_module(
            'picard_insert_size_metrics',
            self._create_params_dict(fastqs=fastq_names,
                                     seq_data_reads=('r1','r2'),
                                     organism="Human",
                                     star_index="/data/indexes/STAR")))
        # Returns None if organism or STAR index is not supplied
        self.assertEqual(None,
                         qc_verifier.verify_qc_module(
                             'picard_insert_size_metrics',
                             self._create_params_dict(
                                 fastqs=fastq_names,
                                 seq_data_reads=('r1','r2'),
                                 star_index="/data/indexes/STAR")))
        self.assertEqual(None,
                         qc_verifier.verify_qc_module(
                             'picard_insert_size_metrics',
                             self._create_params_dict(
                                 fastqs=fastq_names,
                                 seq_data_reads=('r1','r2'),
                                 organism="Human")))
        # Some outputs missing
        qc_dir = self._make_qc_dir('qc.fail',
                                   fastq_names=fastq_names[:-2],
                                   organisms=('human',),
                                   include_picard_insert_size_metrics=True)
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(qc_verifier.verify_qc_module(
            'picard_insert_size_metrics',
            self._create_params_dict(fastqs=fastq_names,
                                     seq_data_reads=('r1','r2'),
                                     organism="Human",
                                     star_index="/data/indexes/STAR")))
        # Organism name contains spaces
        qc_dir = self._make_qc_dir('qc.homo_sapiens',
                                   fastq_names=fastq_names,
                                   organisms=('Homo sapiens',),
                                   include_picard_insert_size_metrics=True)
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify_qc_module(
            'picard_insert_size_metrics',
            self._create_params_dict(fastqs=fastq_names,
                                     seq_data_reads=('r1','r2'),
                                     organism="Homo sapiens",
                                     star_index="/data/indexes/STAR")))
        # Empty QC directory
        qc_dir = self._make_qc_dir('qc.empty',
                                   fastq_names=fastq_names,
                                   organisms=('human',),
                                   include_picard_insert_size_metrics=False)
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(qc_verifier.verify_qc_module(
            'picard_insert_size_metrics',
            self._create_params_dict(fastqs=fastq_names,
                                     seq_data_reads=('r1','r2'),
                                     organism="Human")))

    def test_qcverifier_verify_qc_module_qualimap_rnaseq(self):
        """
        QCVerifier: verify QC module 'qualimap_rnaseq'
        """
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB1_S1_R2_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',
                     'PJB2_S2_R2_001.fastq.gz',)
        # All outputs present
        qc_dir = self._make_qc_dir('qc.ok',
                                   fastq_names=fastq_names,
                                   organisms=('human',),
                                   include_qualimap_rnaseq=True)
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify_qc_module(
            'qualimap_rnaseq',
            self._create_params_dict(fastqs=fastq_names,
                                     seq_data_reads=('r1','r2'),
                                     organism="Human",
                                     star_index="/data/indexes/STAR",
                                     annotation_gtf="/data/annot/human.gtf",
                                     annotation_bed="/data/annot/human.bed")))
        # Returns None if organism, STAR index or annotation
        # is not supplied
        self.assertEqual(None,
                         qc_verifier.verify_qc_module(
                             'qualimap_rnaseq',
                             self._create_params_dict(
                                 fastqs=fastq_names,
                                 seq_data_reads=('r1','r2'),
                                 star_index="/data/indexes/STAR",
                                 annotation_gtf="/data/annot/human.gtf",
                                 annotation_bed="/data/annot/human.bed")))
        self.assertEqual(None,
                         qc_verifier.verify_qc_module(
                             'qualimap_rnaseq',
                             self._create_params_dict(
                                 fastqs=fastq_names,
                                 seq_data_reads=('r1','r2'),
                                 organism="Human",
                                 annotation_gtf="/data/annot/human.gtf",
                                 annotation_bed="/data/annot/human.bed")))
        self.assertEqual(None,
                         qc_verifier.verify_qc_module(
                             'qualimap_rnaseq',
                             self._create_params_dict(
                                 fastqs=fastq_names,
                                 seq_data_reads=('r1','r2'),
                                 organism="Human",
                                 star_index="/data/indexes/STAR")))
        # Organism name contains spaces
        qc_dir = self._make_qc_dir('qc.homo_sapiens',
                                   fastq_names=fastq_names,
                                   organisms=('Homo sapiens',),
                                   include_qualimap_rnaseq=True)
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify_qc_module(
            'qualimap_rnaseq',
            self._create_params_dict(fastqs=fastq_names,
                                     seq_data_reads=('r1','r2'),
                                     organism="Homo sapiens",
                                     star_index="/data/indexes/STAR",
                                     annotation_gtf="/data/annot/human.gtf",
                                     annotation_bed="/data/annot/human.bed")))
        # Empty QC directory
        qc_dir = self._make_qc_dir('qc.empty',
                                   fastq_names=fastq_names,
                                   organisms=('human',),
                                   include_qualimap_rnaseq=False)
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(qc_verifier.verify_qc_module(
            'qualimap_rnaseq',
            self._create_params_dict(fastqs=fastq_names,
                                     seq_data_reads=('r1','r2'),
                                     organism="Human",
                                     star_index="/data/indexes/STAR",
                                     annotation_gtf="/data/annot/human.gtf",
                                     annotation_bed="/data/annot/human.bed")))

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
        self.assertTrue(qc_verifier.verify_qc_module(
            'multiqc',
            self._create_params_dict()))
        # QC dir without MultiQC
        qc_dir = self._make_qc_dir('qc.no_multiqc',
                                   fastq_names=fastq_names,
                                   include_multiqc=False)
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(qc_verifier.verify_qc_module(
            'multiqc',
            self._create_params_dict()))

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
            self._create_params_dict(samples=('PJB1','PJB2'))))
        # Explicitly match version with any reference
        self.assertTrue(qc_verifier.verify_qc_module(
            'cellranger_count',
            self._create_params_dict(samples=('PJB1','PJB2'),
                                     cellranger_version=DEFAULT_CELLRANGER_VERSION,
                                     cellranger_refdata='*')))
        # Explicitly match reference with any version
        self.assertTrue(qc_verifier.verify_qc_module(
            'cellranger_count',
            self._create_params_dict(samples=('PJB1','PJB2'),
                                     cellranger_version='*',
                                     cellranger_refdata=\
                                     'refdata-cellranger-2020-A')))
        # Fail if version not found
        self.assertFalse(qc_verifier.verify_qc_module(
            'cellranger_count',
            self._create_params_dict(samples=('PJB1','PJB2'),
                                     cellranger_version='5.0.0',
                                     cellranger_refdata=\
                                     'refdata-cellranger-2020-A')))
        # Fail if reference not found
        self.assertFalse(qc_verifier.verify_qc_module(
            'cellranger_count',
            self._create_params_dict(samples=('PJB1','PJB2'),
                                     cellranger_version='7.1.0',
                                     cellranger_refdata=\
                                     'refdata-cellranger-2.0.0')))
        # Missing outputs for one sample
        qc_dir = self._make_qc_dir('qc.fail',
                                   fastq_names=fastq_names[:2],
                                   include_cellranger_count=True,
                                   cellranger_pipelines=('cellranger',),
                                   cellranger_samples=('PJB1',))
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(qc_verifier.verify_qc_module(
            'cellranger_count',
            self._create_params_dict(samples=('PJB1','PJB2'),
                                     cellranger_version='7.1.0',
                                     cellranger_refdata=\
                                     'refdata-cellranger-2020-A')))
        # Empty QC dir
        qc_dir = self._make_qc_dir('qc.empty',
                                   fastq_names=fastq_names,
                                   include_cellranger_count=False)
        qc_verifier = QCVerifier(qc_dir)
        # Okay if no reference data specified
        self.assertTrue(qc_verifier.verify_qc_module(
            'cellranger_count',
            self._create_params_dict(samples=('PJB1','PJB2'))))
        # Fail if reference data is specified
        self.assertFalse(qc_verifier.verify_qc_module(
            'cellranger_count',
            self._create_params_dict(samples=('PJB1','PJB2'),
                                     cellranger_refdata='*')))

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
            self._create_params_dict(samples=('PJB1','PJB2'))))
        # Explicitly match version with any reference
        self.assertTrue(qc_verifier.verify_qc_module(
            'cellranger-atac_count',
            self._create_params_dict(samples=('PJB1','PJB2'),
                                cellranger_version='2.0.0',
                                cellranger_refdata='*')))
        # Explicitly match reference with any version
        self.assertTrue(qc_verifier.verify_qc_module(
            'cellranger-atac_count',
            self._create_params_dict(samples=('PJB1','PJB2'),
                                cellranger_version='*',
                                cellranger_refdata=\
                                'refdata-cellranger-atac-2020-A')))
        # Fail if version not found
        self.assertFalse(qc_verifier.verify_qc_module(
            'cellranger-atac_count',
            self._create_params_dict(samples=('PJB1','PJB2'),
                                     cellranger_version='1.0.0',
                                     cellranger_refdata=\
                                     'refdata-cellranger-atac-2020-A')))
        # Fail if reference not found
        self.assertFalse(qc_verifier.verify_qc_module(
            'cellranger-atac_count',
            self._create_params_dict(samples=('PJB1','PJB2'),
                                     cellranger_version='2.0.0',
                                     cellranger_refdata=\
                                     'refdata-cellranger-atac-2.0.0')))
        # Missing outputs for one sample
        qc_dir = self._make_qc_dir('qc.fail',
                                   fastq_names=fastq_names[:2],
                                   include_cellranger_count=True,
                                   cellranger_pipelines=('cellranger-atac',),
                                   cellranger_samples=('PJB1',))
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(qc_verifier.verify_qc_module(
            'cellranger-atac_count',
            self._create_params_dict(samples=('PJB1','PJB2'),
                                     cellranger_version='2.0.0',
                                     cellranger_refdata=\
                                     'refdata-cellranger-atac-2020-A')))
        # Empty QC dir
        qc_dir = self._make_qc_dir('qc.empty',
                                   fastq_names=fastq_names,
                                   include_cellranger_count=False)
        qc_verifier = QCVerifier(qc_dir)
        # Okay if no reference data specified
        self.assertTrue(qc_verifier.verify_qc_module(
            'cellranger-atac_count',
            self._create_params_dict(samples=('PJB1','PJB2'))))
        # Fail if reference data is specified
        self.assertFalse(
            qc_verifier.verify_qc_module(
                'cellranger-atac_count',
                self._create_params_dict(samples=('PJB1','PJB2'),
                                         cellranger_refdata='*')))

    def test_qcverifier_verify_qc_module_cellranger_atac_count_legacy_dir(self):
        """
        QCVerifier: verify QC module 'cellranger-atac_count' (legacy directory name)
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
                                   ),
                                   legacy_cellranger_count_prefix=True)
        qc_verifier = QCVerifier(qc_dir)
        # Implicitly match any version and reference
        self.assertTrue(qc_verifier.verify_qc_module(
            'cellranger-atac_count',
            self._create_params_dict(samples=('PJB1','PJB2'))))
        # Explicitly match version with any reference
        self.assertTrue(qc_verifier.verify_qc_module(
            'cellranger-atac_count',
            self._create_params_dict(samples=('PJB1','PJB2'),
                                cellranger_version='2.0.0',
                                cellranger_refdata='*')))
        # Explicitly match reference with any version
        self.assertTrue(qc_verifier.verify_qc_module(
            'cellranger-atac_count',
            self._create_params_dict(samples=('PJB1','PJB2'),
                                cellranger_version='*',
                                cellranger_refdata=\
                                'refdata-cellranger-atac-2020-A')))
        # Fail if version not found
        self.assertFalse(qc_verifier.verify_qc_module(
            'cellranger-atac_count',
            self._create_params_dict(samples=('PJB1','PJB2'),
                                     cellranger_version='1.0.0',
                                     cellranger_refdata=\
                                     'refdata-cellranger-atac-2020-A')))
        # Fail if reference not found
        self.assertFalse(qc_verifier.verify_qc_module(
            'cellranger-atac_count',
            self._create_params_dict(samples=('PJB1','PJB2'),
                                     cellranger_version='2.0.0',
                                     cellranger_refdata=\
                                     'refdata-cellranger-atac-2.0.0')))
        # Missing outputs for one sample
        qc_dir = self._make_qc_dir('qc.fail',
                                   fastq_names=fastq_names[:2],
                                   include_cellranger_count=True,
                                   cellranger_pipelines=('cellranger-atac',),
                                   cellranger_samples=('PJB1',))
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(qc_verifier.verify_qc_module(
            'cellranger-atac_count',
            self._create_params_dict(samples=('PJB1','PJB2'),
                                     cellranger_version='2.0.0',
                                     cellranger_refdata=\
                                     'refdata-cellranger-atac-2020-A')))
        # Empty QC dir
        qc_dir = self._make_qc_dir('qc.empty',
                                   fastq_names=fastq_names,
                                   include_cellranger_count=False)
        qc_verifier = QCVerifier(qc_dir)
        # Okay if no reference data specified
        self.assertTrue(qc_verifier.verify_qc_module(
            'cellranger-atac_count',
            self._create_params_dict(samples=('PJB1','PJB2'))))
        # Fail if reference data is specified
        self.assertFalse(
            qc_verifier.verify_qc_module(
                'cellranger-atac_count',
                self._create_params_dict(samples=('PJB1','PJB2'),
                                         cellranger_refdata='*')))

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
            self._create_params_dict(samples=('PJB1','PJB2'))))
        # Explicitly match version with any reference
        self.assertTrue(qc_verifier.verify_qc_module(
            'cellranger-arc_count',
            self._create_params_dict(samples=('PJB1','PJB2'),
                                     cellranger_version='2.0.0',
                                     cellranger_refdata='*')))
        # Explicitly match reference with any version
        self.assertTrue(
            qc_verifier.verify_qc_module(
                'cellranger-arc_count',
                self._create_params_dict(samples=('PJB1','PJB2'),
                                         cellranger_version='*',
                                         cellranger_refdata=\
                                         'refdata-cellranger-arc-2020-A')))
        # Fail if version not found
        self.assertFalse(
            qc_verifier.verify_qc_module(
                'cellranger-arc_count',
                self._create_params_dict(samples=('PJB1','PJB2'),
                                         cellranger_version='1.0.0',
                                         cellranger_refdata=\
                                         'refdata-cellranger-arc-2020-A')))
        # Fail if reference not found
        self.assertFalse(
            qc_verifier.verify_qc_module(
                'cellranger-arc_count',
                self._create_params_dict(samples=('PJB1','PJB2'),
                                         cellranger_version='2.0.0',
                                         cellranger_refdata=\
                                         'refdata-cellranger-arc-2.0.0')))
        # Missing outputs for one sample
        qc_dir = self._make_qc_dir('qc.fail',
                                   fastq_names=fastq_names[:2],
                                   include_cellranger_count=True,
                                   cellranger_pipelines=('cellranger-arc',),
                                   cellranger_samples=('PJB1',))
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(
            qc_verifier.verify_qc_module(
                'cellranger-arc_count',
                self._create_params_dict(samples=('PJB1','PJB2'),
                                         cellranger_version='2.0.0',
                                         cellranger_refdata=\
                                         'refdata-cellranger-arc-2020-A')))
        # Empty QC dir
        # NB this will verify as True because the multiome CSV config
        # files are missing (so no outputs are expected)
        qc_dir = self._make_qc_dir('qc.empty',
                                   fastq_names=fastq_names,
                                   include_cellranger_count=False)
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify_qc_module(
            'cellranger-arc_count',
            self._create_params_dict(samples=('PJB1','PJB2'))))

    def test_qcverifier_verify_qc_module_cellranger_arc_count_legacy_dir(self):
        """
        QCVerifier: verify QC module 'cellranger-arc_count' (legacy directory name)
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
                                   ),
                                   legacy_cellranger_count_prefix=True)
        qc_verifier = QCVerifier(qc_dir)
        # Implicitly match any version and reference
        self.assertTrue(qc_verifier.verify_qc_module(
            'cellranger-arc_count',
            self._create_params_dict(samples=('PJB1','PJB2'))))
        # Explicitly match version with any reference
        self.assertTrue(qc_verifier.verify_qc_module(
            'cellranger-arc_count',
            self._create_params_dict(samples=('PJB1','PJB2'),
                                     cellranger_version='2.0.0',
                                     cellranger_refdata='*')))
        # Explicitly match reference with any version
        self.assertTrue(
            qc_verifier.verify_qc_module(
                'cellranger-arc_count',
                self._create_params_dict(samples=('PJB1','PJB2'),
                                         cellranger_version='*',
                                         cellranger_refdata=\
                                         'refdata-cellranger-arc-2020-A')))
        # Fail if version not found
        self.assertFalse(
            qc_verifier.verify_qc_module(
                'cellranger-arc_count',
                self._create_params_dict(samples=('PJB1','PJB2'),
                                         cellranger_version='1.0.0',
                                         cellranger_refdata=\
                                         'refdata-cellranger-arc-2020-A')))
        # Fail if reference not found
        self.assertFalse(
            qc_verifier.verify_qc_module(
                'cellranger-arc_count',
                self._create_params_dict(samples=('PJB1','PJB2'),
                                         cellranger_version='2.0.0',
                                         cellranger_refdata=\
                                         'refdata-cellranger-arc-2.0.0')))
        # Missing outputs for one sample
        qc_dir = self._make_qc_dir('qc.fail',
                                   fastq_names=fastq_names[:2],
                                   include_cellranger_count=True,
                                   cellranger_pipelines=('cellranger-arc',),
                                   cellranger_samples=('PJB1',))
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(
            qc_verifier.verify_qc_module(
                'cellranger-arc_count',
                self._create_params_dict(samples=('PJB1','PJB2'),
                                         cellranger_version='2.0.0',
                                         cellranger_refdata=\
                                         'refdata-cellranger-arc-2020-A')))
        # Empty QC dir
        # NB this will verify as True because the multiome CSV config
        # files are missing (so no outputs are expected)
        qc_dir = self._make_qc_dir('qc.empty',
                                   fastq_names=fastq_names,
                                   include_cellranger_count=False)
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify_qc_module(
            'cellranger-arc_count',
            self._create_params_dict(samples=('PJB1','PJB2'))))

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
        # Verification not possible if no version and no reference
        # explicitly specified
        self.assertEqual(qc_verifier.verify_qc_module(
            'cellranger_multi',
            self._create_params_dict(qc_dir=qc_dir,
                                     cellranger_version=None,
                                     cellranger_refdata=None)),
                         None)
        # Explicitly match version with any reference
        self.assertTrue(qc_verifier.verify_qc_module(
            'cellranger_multi',
            self._create_params_dict(qc_dir=qc_dir,
                                     cellranger_version=DEFAULT_CELLRANGER_VERSION,
                                     cellranger_refdata='*')))
        # Explicitly match reference with any version
        self.assertTrue(qc_verifier.verify_qc_module(
            'cellranger_multi',
            self._create_params_dict(qc_dir=qc_dir,
                                     cellranger_version='*',
                                     cellranger_refdata=\
                                     'refdata-cellranger-2020-A')))
        # Explicitly match reference with any version and any
        # reference
        self.assertTrue(qc_verifier.verify_qc_module(
            'cellranger_multi',
            self._create_params_dict(qc_dir=qc_dir,
                                     cellranger_version='*',
                                     cellranger_refdata='*')))
        # Fail if version not found
        self.assertFalse(qc_verifier.verify_qc_module(
            'cellranger_multi',
            self._create_params_dict(qc_dir=qc_dir,
                                     cellranger_version='5.0.0',
                                     cellranger_refdata=\
                                     'refdata-cellranger-2020-A')))
        # Fail if reference not found
        self.assertFalse(qc_verifier.verify_qc_module(
            'cellranger_multi',
            self._create_params_dict(qc_dir=qc_dir,
                                     cellranger_version='7.1.0',
                                     cellranger_refdata=\
                                     'refdata-cellranger-2.0.0')))
        # Verification when no multiplexed samples present
        qc_dir = self._make_qc_dir('qc.no_multiplexed',
                                   fastq_names=fastq_names[:2],
                                   include_cellranger_multi=True,
                                   cellranger_pipelines=('cellranger',))
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify_qc_module(
            'cellranger_multi',
            self._create_params_dict(qc_dir=qc_dir,
                                     cellranger_version=DEFAULT_CELLRANGER_VERSION,
                                     cellranger_refdata='*')))
        # Missing outputs for one sample
        qc_dir = self._make_qc_dir('qc.fail',
                                   fastq_names=fastq_names[:2],
                                   include_cellranger_multi=True,
                                   cellranger_pipelines=('cellranger',),
                                   cellranger_multi_samples=('PJB_CML1',))
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(qc_verifier.verify_qc_module(
            'cellranger_multi',
            self._create_params_dict(qc_dir=qc_dir,
                                     samples=('PJB1','PJB2'),
                                     cellranger_version='7.1.0',
                                     cellranger_refdata=\
                                     'refdata-cellranger-2020-A')))
        # Empty QC dir
        # NB this will verify as None because the 10x multi CSV config
        # file is missing (so no outputs are expected)
        qc_dir = self._make_qc_dir('qc.empty',
                                   fastq_names=fastq_names,
                                   include_cellranger_multi=False)
        qc_verifier = QCVerifier(qc_dir)
        self.assertEqual(qc_verifier.verify_qc_module(
            'cellranger_multi',
            self._create_params_dict(qc_dir=qc_dir)),
                         None)

    def test_qcverifier_verify_qc_module_cellranger_multi_multiple_samples(self):
        """
        QCVerifier: verify QC module 'cellranger_multi' (multiple samples)
        """
        fastq_names = ("PJB1_GEX_S1_R1_001.fastq.gz",
                       "PJB1_GEX_S1_R2_001.fastq.gz",
                       "PJB1_CML_S2_R1_001.fastq.gz",
                       "PJB1_CML_S2_R2_001.fastq.gz",
                       "PJB2_GEX_S3_R1_001.fastq.gz",
                       "PJB2_GEX_S3_R2_001.fastq.gz",
                       "PJB2_CML_S4_R1_001.fastq.gz",
                       "PJB2_CML_S4_R2_001.fastq.gz",)
        # QC dir with cellranger multi outputs
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=fastq_names,
                                   include_cellranger_multi=True,
                                   cellranger_pipelines=('cellranger',),
                                   cellranger_multi_samples={
                                       "PJB1": ("PB1", "PB2"),
                                       "PJB2": ("PB3", "PB4")
                                   })
        qc_verifier = QCVerifier(qc_dir)
        # Verification not possible if no version and no reference
        # were explicitly supplied
        self.assertEqual(qc_verifier.verify_qc_module(
            'cellranger_multi',
            self._create_params_dict(qc_dir=qc_dir)),
                         None)
        # Explicitly match version with any reference
        self.assertTrue(qc_verifier.verify_qc_module(
            'cellranger_multi',
            self._create_params_dict(qc_dir=qc_dir,
                                     cellranger_version=DEFAULT_CELLRANGER_VERSION,
                                     cellranger_refdata='*')))
        # Explicitly match reference with any version
        self.assertTrue(qc_verifier.verify_qc_module(
            'cellranger_multi',
            self._create_params_dict(qc_dir=qc_dir,
                                     cellranger_version='*',
                                     cellranger_refdata=\
                                     'refdata-cellranger-2020-A')))
        # Explicitly match reference with any version and any
        # reference
        self.assertTrue(qc_verifier.verify_qc_module(
            'cellranger_multi',
            self._create_params_dict(qc_dir=qc_dir,
                                     cellranger_version='*',
                                     cellranger_refdata='*')))
        # Fail if version not found
        self.assertFalse(qc_verifier.verify_qc_module(
            'cellranger_multi',
            self._create_params_dict(qc_dir=qc_dir,
                                     cellranger_version='5.0.0',
                                     cellranger_refdata=\
                                     'refdata-cellranger-2020-A')))
        # Fail if reference not found
        self.assertFalse(qc_verifier.verify_qc_module(
            'cellranger_multi',
            self._create_params_dict(qc_dir=qc_dir,
                                     cellranger_version='7.1.0',
                                     cellranger_refdata=\
                                     'refdata-cellranger-2.0.0')))
        # Verification when no multiplexed samples present
        qc_dir = self._make_qc_dir('qc.no_multiplexed',
                                   fastq_names=fastq_names[:2],
                                   include_cellranger_multi=True,
                                   cellranger_pipelines=('cellranger',),
                                   cellranger_multi_samples={
                                       "PJB1": (),
                                       "PJB2": ()
                                   })
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify_qc_module(
            'cellranger_multi',
            self._create_params_dict(qc_dir=qc_dir,
                                     cellranger_version=DEFAULT_CELLRANGER_VERSION,
                                     cellranger_refdata='*')))
        # Missing outputs for one sample
        qc_dir = self._make_qc_dir('qc.fail',
                                   fastq_names=fastq_names[:2],
                                   include_cellranger_multi=True,
                                   cellranger_pipelines=('cellranger',),
                                   cellranger_multi_samples=('PJB_CML1',))
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(qc_verifier.verify_qc_module(
            'cellranger_multi',
            self._create_params_dict(qc_dir=qc_dir,
                                     samples=('PJB1','PJB2'),
                                     cellranger_version='7.1.0',
                                     cellranger_refdata=\
                                     'refdata-cellranger-2020-A')))
        # Empty QC dir
        # NB this will verify as None because the 10x multi CSV config
        # file is missing (so no outputs are expected)
        qc_dir = self._make_qc_dir('qc.empty',
                                   fastq_names=fastq_names,
                                   include_cellranger_multi=False)
        qc_verifier = QCVerifier(qc_dir)
        self.assertEqual(qc_verifier.verify_qc_module(
            'cellranger_multi',
            self._create_params_dict(qc_dir=qc_dir)),
                         None)

    def test_qcverifier_verify_single_end(self):
        """
	QCVerifier: verify single-end data (standardSE)
        """
        ##self.remove_test_outputs = False
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',)
        qc_dir = self._make_qc_dir('qc',
                                   protocol="standardSE",
                                   fastq_names=fastq_names)
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify(
            fetch_protocol_definition("standardSE"),
            fastq_names,
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
                                   protocol="standardSE",
                                   fastq_names=fastq_names,
                                   include_fastqc=False)
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(qc_verifier.verify(
            fetch_protocol_definition("standardSE"),
            fastq_names,
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
                                   protocol="standardSE",
                                   fastq_names=fastq_names,
                                   include_fastq_screen=False)
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(qc_verifier.verify(
            fetch_protocol_definition("standardSE"),
            fastq_names,
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
                                   protocol="standardSE",
                                   fastq_names=fastq_names,
                                   include_seqlens=False)
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(qc_verifier.verify(
            fetch_protocol_definition("standardSE"),
            fastq_names,
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
                                   protocol="standardSE",
                                   fastq_names=fastq_names,
                                   legacy_screens=True)
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify(
            fetch_protocol_definition("standardSE"),
            fastq_names,
            fastq_screens=('model_organisms',
                           'other_organisms',
                           'rRNA')))

    def test_qcverifier_verify_single_end_biological_samples(self):
        """
	QCVerifier: verify single-end data (standardSE, biological samples defined)
        """
        ##self.remove_test_outputs = False
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',)
        qc_dir = self._make_qc_dir('qc',
                                   protocol="standardSE",
                                   fastq_names=fastq_names,
                                   seq_data_samples=["PJB1"])
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify(
            fetch_protocol_definition("standardSE"),
            fastq_names,
            seq_data_samples="PJB1",
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
                                   protocol="standardPE",
                                   fastq_names=fastq_names)
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify(
            fetch_protocol_definition("standardPE"),
            fastq_names,
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
                                   protocol="standardPE",
                                   fastq_names=fastq_names,
                                   include_fastqc=False)
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(qc_verifier.verify(
            fetch_protocol_definition("standardPE"),
            fastq_names,
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
                                   protocol="standardPE",
                                   fastq_names=fastq_names,
                                   include_fastq_screen=False)
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(qc_verifier.verify(
            fetch_protocol_definition("standardPE"),
            fastq_names,
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
                                   protocol="standardPE",
                                   fastq_names=fastq_names,
                                   include_seqlens=False)
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(qc_verifier.verify(
            fetch_protocol_definition("standardPE"),
            fastq_names,
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
                                   protocol="standardPE",
                                   fastq_names=fastq_names,
                                   legacy_screens=True)
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify(
            fetch_protocol_definition("standardPE"),
            fastq_names,
            fastq_screens=('model_organisms',
                           'other_organisms',
                           'rRNA')))

    def test_qcverifier_verify_paired_end_biological_samples(self):
        """
	QCVerifier: verify paired-end data (standardPE, biological samples defined)
        """
        ##self.remove_test_outputs = False
        fastq_names=('PJB1_S1_R1_001_paired.fastq.gz',
                     'PJB1_S1_R2_001_paired.fastq.gz',
                     'PJB2_S2_R1_001_paired.fastq.gz',
                     'PJB2_S2_R2_001_paired.fastq.gz',)
        qc_dir = self._make_qc_dir('qc',
                                   protocol="standardPE",
                                   fastq_names=fastq_names,
                                   seq_data_samples=["PJB1"])
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify(
            fetch_protocol_definition("standardPE"),
            fastq_names,
            seq_data_samples="PJB1",
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
                                   protocol="standardPE",
                                   fastq_names=fastq_names)
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify(
            fetch_protocol_definition("standardPE"),
            fastq_names,
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
                                   protocol="10x_scRNAseq",
                                   fastq_names=fastq_names,
                                   include_cellranger_count=True,
                                   cellranger_pipelines=('cellranger',),
                                   cellranger_samples=(
                                       'PJB1',
                                       'PJB2',
                                   ))
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify(
            fetch_protocol_definition("10x_scRNAseq"),
            fastq_names,
            fastq_screens=('model_organisms',
                           'other_organisms',
                           'rRNA'),
            cellranger_version=DEFAULT_CELLRANGER_VERSION,
            cellranger_refdata="/data/refdata-cellranger-2020-A"))

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
                                   protocol="10x_scRNAseq",
                                   fastq_names=fastq_names,
                                   include_cellranger_count=True,
                                   cellranger_pipelines=('cellranger',),
                                   cellranger_samples=(
                                       'PJB1',
                                       'PJB2',
                                   ))
        qc_verifier = QCVerifier(qc_dir)
        self.assertFalse(qc_verifier.verify(
            fetch_protocol_definition("10x_scRNAseq"),
            fastq_names,
            fastq_screens=('model_organisms',
                           'other_organisms',
                           'rRNA'),
            cellranger_version="5.0.0",
            cellranger_refdata="/data/refdata-cellranger-2020-A"))

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
                                   protocol="10x_scRNAseq",
                                   fastq_names=fastq_names,
                                   include_cellranger_count=True,
                                   cellranger_pipelines=('cellranger',),
                                   cellranger_samples=(
                                       'PJB1',
                                       'PJB2',
                                   ),
                                   legacy_cellranger_outs=True)
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify(
            fetch_protocol_definition("10x_scRNAseq"),
            fastq_names,
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
                                   protocol="10x_scATAC",
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
        self.assertTrue(qc_verifier.verify(
            fetch_protocol_definition("10x_scATAC"),
            fastq_names,
            fastq_screens=('model_organisms',
                           'other_organisms',
                           'rRNA'),
            cellranger_version="2.0.0",
            cellranger_refdata="/data/refdata-cellranger-atac-2020-A"))

    def test_qcverifier_10x_cellranger_atac_count_legacy_dir(self):
        """
        QCVerifier: verify 10xGenomics scATAC-seq data (10x_scATAC) (legacy directory name)
        """
        ##self.remove_test_outputs = False
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB1_S1_R3_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',
                     'PJB2_S2_R3_001.fastq.gz',)
        qc_dir = self._make_qc_dir('qc',
                                   protocol="10x_scATAC",
                                   fastq_names=fastq_names,
                                   include_cellranger_count=True,
                                   cellranger_pipelines=(
                                       'cellranger-atac',
                                   ),
                                   cellranger_samples=(
                                       'PJB1',
                                       'PJB2',
                                   ),
                                   legacy_cellranger_count_prefix=True)
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify(
            fetch_protocol_definition("10x_scATAC"),
            fastq_names,
            fastq_screens=('model_organisms',
                           'other_organisms',
                           'rRNA'),
            cellranger_version="2.0.0",
            cellranger_refdata="/data/refdata-cellranger-atac-2020-A"))

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
                                   protocol="10x_Multiome_GEX",
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
        self.assertTrue(qc_verifier.verify(
            fetch_protocol_definition("10x_Multiome_GEX"),
            fastq_names,
            fastq_screens=('model_organisms',
                           'other_organisms',
                           'rRNA'),
            cellranger_version="2.0.0",
            cellranger_refdata="refdata-cellranger-arc-2020-A"))

    def test_verify_qcverifier_10x_multiome_gex_legacy_dir(self):
        """
        QCVerifier: verify 10xGenomics multiome GEX data (10x_Multiome_GEX) (legacy directory name)
        """
        ##self.remove_test_outputs = False
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB1_S1_R2_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',
                     'PJB2_S2_R2_001.fastq.gz',)
        qc_dir = self._make_qc_dir('qc',
                                   protocol="10x_Multiome_GEX",
                                   fastq_names=fastq_names,
                                   include_cellranger_count=True,
                                   cellranger_pipelines=(
                                       'cellranger',
                                       'cellranger-arc',
                                   ),
                                   cellranger_samples=(
                                       'PJB1',
                                       'PJB2',
                                   ),
                                   legacy_cellranger_count_prefix=True)
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify(
            fetch_protocol_definition("10x_Multiome_GEX"),
            fastq_names,
            fastq_screens=('model_organisms',
                           'other_organisms',
                           'rRNA'),
            cellranger_version="2.0.0",
            cellranger_refdata="refdata-cellranger-arc-2020-A"))

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
                                   protocol="10x_Multiome_GEX",
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
        self.assertTrue(qc_verifier.verify(
            fetch_protocol_definition("10x_Multiome_GEX"),
            fastq_names,
            fastq_screens=('model_organisms',
                           'other_organisms',
                           'rRNA'),
            cellranger_version="2.0.0",
            cellranger_refdata="refdata-cellranger-arc-2020-A"))

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
                                   protocol="10x_Multiome_GEX",
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
        self.assertFalse(qc_verifier.verify(
            fetch_protocol_definition("10x_Multiome_GEX"),
            fastq_names,
            fastq_screens=('model_organisms',
                           'other_organisms',
                           'rRNA'),
            cellranger_version="2.0.0",
            cellranger_refdata="refdata-cellranger-arc-2020-A"))

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
                                   protocol="10x_Multiome_ATAC",
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
        self.assertTrue(qc_verifier.verify(
            fetch_protocol_definition("10x_Multiome_ATAC"),
            fastq_names,
            fastq_screens=('model_organisms',
                           'other_organisms',
                           'rRNA'),
            cellranger_version="2.0.0",
            cellranger_refdata="refdata-cellranger-arc-2020-A"))

    def test_qcverifier_10x_multiome_atac_legacy_dir(self):
        """
        QCVerifier: verify 10xGenomics multiome ATAC data (10x_Multiome_ATAC) (legacy directory name)
        """
        ##self.remove_test_outputs = False
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB1_S1_R3_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',
                     'PJB2_S2_R3_001.fastq.gz',)
        qc_dir = self._make_qc_dir('qc',
                                   protocol="10x_Multiome_ATAC",
                                   fastq_names=fastq_names,
                                   include_cellranger_count=True,
                                   cellranger_pipelines=(
                                       'cellranger-arc',
                                       'cellranger-atac',
                                   ),
                                   cellranger_samples=(
                                       'PJB1',
                                       'PJB2',
                                   ),
                                   legacy_cellranger_count_prefix=True)
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify(
            fetch_protocol_definition("10x_Multiome_ATAC"),
            fastq_names,
            fastq_screens=('model_organisms',
                           'other_organisms',
                           'rRNA'),
            cellranger_version="2.0.0",
            cellranger_refdata="refdata-cellranger-arc-2020-A"))

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
                                   protocol="10x_Multiome_ATAC",
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
        self.assertTrue(qc_verifier.verify(
            fetch_protocol_definition("10x_Multiome_ATAC"),
            fastq_names,
            fastq_screens=('model_organisms',
                           'other_organisms',
                           'rRNA'),
            cellranger_version="2.0.0",
            cellranger_refdata="refdata-cellranger-arc-2020-A"))

    def test_qcverifier_10x_multiome_atac_no_paired_samples_legacy_dir(self):
        """
        QCVerifier: verify 10xGenomics multiome ATAC data (10x_Multiome_ATAC, no paired samples) (legacy dir name)
        """
        ##self.remove_test_outputs = False
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB1_S1_R3_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',
                     'PJB2_S2_R3_001.fastq.gz',)
        qc_dir = self._make_qc_dir('qc',
                                   protocol="10x_Multiome_ATAC",
                                   fastq_names=fastq_names,
                                   include_cellranger_count=True,
                                   cellranger_pipelines=(
                                       'cellranger-atac',
                                   ),
                                   cellranger_samples=(
                                       'PJB1',
                                       'PJB2',
                                   ),
                                   legacy_cellranger_count_prefix=True)
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify(
            fetch_protocol_definition("10x_Multiome_ATAC"),
            fastq_names,
            fastq_screens=('model_organisms',
                           'other_organisms',
                           'rRNA'),
            cellranger_version="2.0.0",
            cellranger_refdata="refdata-cellranger-arc-2020-A"))

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
                                   protocol="10x_Multiome_ATAC",
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
        self.assertFalse(qc_verifier.verify(
            fetch_protocol_definition("10x_Multiome_ATAC"),
            fastq_names,
            fastq_screens=('model_organisms',
                           'other_organisms',
                           'rRNA'),
            cellranger_version="2.0.0",
            cellranger_refdata="refdata-cellranger-arc-2020-A"))

    def test_qcverifier_verify_10x_visium_gex(self):
        """
        QCVerifier: verify 10xGenomics Visium data (10x_Visium_GEX)
        """
        ##self.remove_test_outputs = False
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB1_S1_R2_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',
                     'PJB2_S2_R2_001.fastq.gz',)
        qc_dir = self._make_qc_dir('qc',
                                   protocol="10x_Visium_GEX",
                                   fastq_names=fastq_names)
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify(
            fetch_protocol_definition("10x_Visium_GEX"),
            fastq_names,
            fastq_screens=('model_organisms',
                           'other_organisms',
                           'rRNA')))

    def test_qcverifier_verify_10x_visium_gex_90bp_insert(self):
        """
        QCVerifier: verify 10xGenomics Visium data (10x_Visium_GEX_90bp_insert)
        """
        ##self.remove_test_outputs = False
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB1_S1_R2_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',
                     'PJB2_S2_R2_001.fastq.gz',)
        qc_dir = self._make_qc_dir('qc',
                                   protocol="10x_Visium",
                                   fastq_names=fastq_names)
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify(
            fetch_protocol_definition("10x_Visium_GEX_90bp_insert"),
            fastq_names,
            fastq_screens=('model_organisms',
                           'other_organisms',
                           'rRNA')))

    def test_qcverifier_verify_10x_visium_pex(self):
        """
        QCVerifier: verify 10xGenomics Visium data (10x_Visium_PEX)
        """
        ##self.remove_test_outputs = False
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB1_S1_R2_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',
                     'PJB2_S2_R2_001.fastq.gz',)
        qc_dir = self._make_qc_dir('qc',
                                   protocol="10x_Visium_PEX",
                                   fastq_names=fastq_names)
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify(
            fetch_protocol_definition("10x_Visium_PEX"),
            fastq_names,
            fastq_screens=('model_organisms',
                           'other_organisms',
                           'rRNA')))

    def test_qcverifier_verify_10x_visium_legacy(self):
        """
        QCVerifier: verify 10xGenomics Visium data (10x_Visium_legacy)
        """
        ##self.remove_test_outputs = False
        fastq_names=('PJB1_S1_R1_001.fastq.gz',
                     'PJB1_S1_R2_001.fastq.gz',
                     'PJB2_S2_R1_001.fastq.gz',
                     'PJB2_S2_R2_001.fastq.gz',)
        qc_dir = self._make_qc_dir('qc',
                                   protocol="10x_Visium_legacy",
                                   fastq_names=fastq_names)
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify(
            fetch_protocol_definition("10x_Visium_legacy"),
            fastq_names,
            fastq_screens=('model_organisms',
                           'other_organisms',
                           'rRNA')))

    def test_qcverifier_verify_10x_cellranger_multi(self):
        """
        QCVerifier: verify 10xGenomics CellPlex data (10x_CellPlex)
        """
        ##self.remove_test_outputs = False
        fastq_names=('PJB1_GEX_R1_001.fastq.gz',
                     'PJB1_GEX_R2_001.fastq.gz',
                     'PJB2_MC_R1_001.fastq.gz',
                     'PJB2_MC_R2_001.fastq.gz',)
        qc_dir = self._make_qc_dir('qc',
                                   protocol="10x_CellPlex",
                                   fastq_names=fastq_names,
                                   seq_data_samples=("PJB1_GEX",),
                                   include_rseqc_genebody_coverage=True,
                                   include_rseqc_infer_experiment=True,
                                   include_qualimap_rnaseq=True,
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
        verify_params = dict()
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify(
            fetch_protocol_definition("10x_CellPlex"),
            fastq_names,
            organism="Human",
            fastq_screens=('model_organisms',
                           'other_organisms',
                           'rRNA'),
            star_index="/data/star/hg38",
            annotation_bed="/data/annotation/hg38.bed",
            annotation_gtf="/data/annotation/hg38.gtf",
            cellranger_version=DEFAULT_CELLRANGER_VERSION,
            cellranger_refdata="/data/refdata-cellranger-2020-A"))

    def test_qcverifier_verify_10x_cellranger_multi_no_config(self):
        """
        QCVerifier: verify 10xGenomics CellPlex data (10x_CellPlex, no multi config file)
        """
        ##self.remove_test_outputs = False
        fastq_names=('PJB1_GEX_R1_001.fastq.gz',
                     'PJB1_GEX_R2_001.fastq.gz',
                     'PJB2_MC_R1_001.fastq.gz',
                     'PJB2_MC_R2_001.fastq.gz',)
        qc_dir = self._make_qc_dir('qc',
                                   protocol="10x_CellPlex",
                                   fastq_names=fastq_names,
                                   include_cellranger_count=False,
                                   include_cellranger_multi=False)
        qc_verifier = QCVerifier(qc_dir)
        self.assertTrue(qc_verifier.verify(
            fetch_protocol_definition("10x_CellPlex"),
            fastq_names,
            fastq_screens=('model_organisms',
                           'other_organisms',
                           'rRNA'),
            cellranger_version="7.1.0",
            cellranger_refdata="/data/refdata-cellranger-2020-A"))

    def test_qcverifier_identify_seq_data(self):
        """
        QCVerifier.identify_seq_data: no 10x multi config
        """
        fastq_names=('PJB1_R1_001.fastq.gz',
                     'PJB1_R2_001.fastq.gz',
                     'PJB2_R1_001.fastq.gz',
                     'PJB2_R2_001.fastq.gz',)
        qc_dir = self._make_qc_dir('qc',
                                   protocol="RNA-seq",
                                   fastq_names=fastq_names,
                                   include_cellranger_count=False,
                                   include_cellranger_multi=False)
        qc_verifier = QCVerifier(qc_dir)
        self.assertEqual(qc_verifier.identify_seq_data(["PJB1", "PJB2"]),
                         ["PJB1", "PJB2"])

    def test_qcverifier_identify_seq_data_using_10x_multi_config(self):
        """
        QCVerifier.identify_seq_data: using 10x multi config
        """
        fastq_names=('PJB1_GEX_R1_001.fastq.gz',
                     'PJB1_GEX_R2_001.fastq.gz',
                     'PJB1_CML_R1_001.fastq.gz',
                     'PJB1_CML_R2_001.fastq.gz',)
        qc_dir = self._make_qc_dir('qc',
                                   protocol="10x_CellPlex",
                                   fastq_names=fastq_names,
                                   include_cellranger_count=False,
                                   include_cellranger_multi=False)
        with open(os.path.join(qc_dir, "10x_multi_config.csv"), "wt") as fp:
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A

[vdj]
reference,/data/vdj_ref.csv

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB1_GEX,/data/runs/fastqs,any,PJB1,Gene Expression,
PJB2_CML,/data/runs/fastqs,any,PJB1,Multiplexing capture,
""")
        qc_verifier = QCVerifier(qc_dir)
        self.assertEqual(qc_verifier.identify_seq_data(["PJB1_GEX",
                                                        "PJB1_CML"]),
                         ["PJB1_GEX"])

    def test_qcverifier_identify_seq_data_using_multiple_10x_multi_configs(self):
        """
        QCVerifier.identify_seq_data: using multiple 10x multi configs
        """
        fastq_names=('PJB1_GEX_R1_001.fastq.gz',
                     'PJB1_GEX_R2_001.fastq.gz',
                     'PJB1_CML_R1_001.fastq.gz',
                     'PJB1_CML_R2_001.fastq.gz',
                     'PJB2_GEX_R1_001.fastq.gz',
                     'PJB2_GEX_R2_001.fastq.gz',
                     'PJB2_CML_R1_001.fastq.gz',
                     'PJB2_CML_R2_001.fastq.gz',)
        qc_dir = self._make_qc_dir('qc',
                                   protocol="10x_CellPlex",
                                   fastq_names=fastq_names,
                                   include_cellranger_count=False,
                                   include_cellranger_multi=False)
        with open(os.path.join(qc_dir, "10x_multi_config.PJB1.csv"),
                  "wt") as fp:
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A

[vdj]
reference,/data/vdj_ref.csv

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB1_GEX,/data/runs/fastqs,any,PJB1,Gene Expression,
PJB1_CML,/data/runs/fastqs,any,PJB1,Multiplexing capture,
""")
        with open(os.path.join(qc_dir, "10x_multi_config.PJB2.csv"),
                  "wt") as fp:
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A

[vdj]
reference,/data/vdj_ref.csv

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB2_GEX,/data/runs/fastqs,any,PJB2,Gene Expression,
PJB2_CML,/data/runs/fastqs,any,PJB2,Multiplexing capture,
""")
        qc_verifier = QCVerifier(qc_dir)
        self.assertEqual(qc_verifier.identify_seq_data(["PJB1_GEX",
                                                        "PJB1_CML",
                                                        "PJB2_GEX",
                                                        "PJB2_CML"]),
                         ["PJB1_GEX", "PJB2_GEX"])

    def test_qcverifier_identify_seq_data_handle_invalid_10x_multi_config(self):
        """
        QCVerifier.identify_seq_data: handle invalid 10x multi config
        """
        fastq_names=('PJB1_GEX_R1_001.fastq.gz',
                     'PJB1_GEX_R2_001.fastq.gz',
                     'PJB1_CML_R1_001.fastq.gz',
                     'PJB1_CML_R2_001.fastq.gz',
                     'PJB2_GEX_R1_001.fastq.gz',
                     'PJB2_GEX_R2_001.fastq.gz',
                     'PJB2_CML_R1_001.fastq.gz',
                     'PJB2_CML_R2_001.fastq.gz',)
        qc_dir = self._make_qc_dir('qc',
                                   protocol="10x_CellPlex",
                                   fastq_names=fastq_names,
                                   include_cellranger_count=False,
                                   include_cellranger_multi=False)
        with open(os.path.join(qc_dir, "10x_multi_config.csv"),
                  "wt") as fp:
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A

[vdj]
reference,/data/vdj_ref.csv

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB1_GEX,/data/runs/fastqs,any,PJB1,Gene Expression,
PJB1_CML,/data/runs/fastqs,any,PJB1,Multiplexing capture,
PJB2_GEX,/data/runs/fastqs,any,PJB2,Gene Expression,
PJB2_CML,/data/runs/fastqs,any,PJB2,Multiplexing capture,
""")
        qc_verifier = QCVerifier(qc_dir)
        self.assertEqual(qc_verifier.identify_seq_data(["PJB1_GEX",
                                                        "PJB1_CML",
                                                        "PJB2_GEX",
                                                        "PJB2_CML"]),
                         ["PJB1_GEX", "PJB2_GEX"])

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

    def _make_analysis_project(self,name="PJB",
                               protocol=None,paired_end=True,
                               fastq_dir="fastqs",qc_dir="qc",
                               fastq_names=None,
                               sample_names=None,
                               seq_data_samples=None,
                               screens=('model_organisms',
                                        'other_organisms',
                                        'rRNA',),
                               include_fastqc=True,
                               include_fastq_screen=True,
                               include_strandedness=True,
                               include_seqlens=True,
                               include_picard_insert_size_metrics=True,
                               include_rseqc_genebody_coverage=True,
                               include_rseqc_infer_experiment=True,
                               include_qualimap_rnaseq=True,
                               include_multiqc=True,
                               include_cellranger_count=False,
                               include_cellranger_multi=False,
                               cellranger_pipelines=('cellranger',),
                               cellranger_samples=None,
                               cellranger_multi_samples=None,
                               legacy_screens=False,
                               legacy_cellranger_outs=False):
        # Create a mock Analysis Project directory
        self._make_working_dir()
        return make_mock_analysis_project(
            name=name,
            top_dir=self.wd,
            protocol=protocol,
            paired_end=paired_end,
            fastq_dir=fastq_dir,
            qc_dir=qc_dir,
            fastq_names=fastq_names,
            sample_names=sample_names,
            seq_data_samples=seq_data_samples,
            screens=screens,
            include_fastqc=include_fastqc,
            include_fastq_screen=include_fastq_screen,
            include_strandedness=include_strandedness,
            include_seqlens=include_seqlens,
            include_picard_insert_size_metrics=\
            include_picard_insert_size_metrics,
            include_rseqc_genebody_coverage=\
            include_rseqc_genebody_coverage,
            include_rseqc_infer_experiment=\
            include_rseqc_infer_experiment,
            include_qualimap_rnaseq=include_qualimap_rnaseq,
            include_multiqc=include_multiqc,
            include_cellranger_count=include_cellranger_count,
            include_cellranger_multi=include_cellranger_multi,
            cellranger_pipelines=cellranger_pipelines,
            cellranger_samples=cellranger_samples,
            cellranger_multi_samples=cellranger_multi_samples,
            legacy_screens=legacy_screens,
            legacy_cellranger_outs=legacy_cellranger_outs)

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

    def test_verify_project_paired_end_with_biological_samples(self):
        """
        verify_project: paired-end data: biological samples defined
        """
        analysis_dir = self._make_analysis_project(
            protocol="standardPE",
            seq_data_samples=('PJB1',))
        project = AnalysisProject(analysis_dir)
        self.assertTrue(verify_project(project))

    def test_verify_project_paired_end_no_seq_lens(self):
        """
        verify_project: paired-end data: no sequence lengths
        """
        analysis_dir = self._make_analysis_project(protocol='standardPE',
                                                   include_seqlens=False)
        project = AnalysisProject(analysis_dir)
        self.assertFalse(verify_project(project))

    def test_verify_project_paired_end_cellranger_count(self):
        """
        verify_project: paired-end data with cellranger 'count'
        """
        analysis_dir = self._make_analysis_project(
            protocol='10x_scRNAseq',
            include_cellranger_count=True,
            cellranger_pipelines=('cellranger',),
            cellranger_samples=('PJB1','PJB2',))
        project = AnalysisProject(analysis_dir)
        self.assertTrue(verify_project(project))

    def test_verify_project_paired_end_cellranger_multi(self):
        """
        verify_project: paired-end data with cellranger 'multi'
        """
        analysis_dir = self._make_analysis_project(
            protocol="10x_CellPlex",
            sample_names=("PJB1_GEX", "PJB1_CML",),
            seq_data_samples=("PJB1_GEX",),
            include_cellranger_multi=True,
            cellranger_multi_samples=("PB1", "PB2"))
        project = AnalysisProject(analysis_dir)
        # Should fail to verify as no cellranger count outputs
        self.assertFalse(verify_project(project))

    def test_verify_project_paired_end_cellranger_count_and_multi(self):
        """
        verify_project: paired-end data with cellranger 'count' and 'multi'
        """
        analysis_dir = self._make_analysis_project(
            protocol="10x_CellPlex",
            sample_names=("PJB1_GEX", "PJB1_CML",),
            seq_data_samples=("PJB1_GEX",),
            include_cellranger_multi=True,
            include_cellranger_count=True,
            cellranger_pipelines=('cellranger',),
            # NB only GEX samples
            cellranger_samples=("PJB1_GEX",),
            cellranger_multi_samples=("PB1", "PB2"))
        project = AnalysisProject(analysis_dir)
        self.assertTrue(verify_project(project))

    def test_verify_project_paired_end_cellranger_multi_multiple_samples(self):
        """
        verify_project: paired-end data with cellranger 'multi' (multiple samples)
        """
        analysis_dir = self._make_analysis_project(
            protocol="10x_CellPlex",
            sample_names=("PJB1_GEX", "PJB1_CML", "PJB2_GEX", "PJB2_CML"),
            seq_data_samples=("PJB1_GEX", "PJB2_GEX"),
            include_cellranger_multi=True,
            cellranger_multi_samples={ "PJB1": ("PB1", "PB2"),
                                       "PJB2": ("PB3", "PB4") })
        project = AnalysisProject(analysis_dir)
        # Should fail to verify as no cellranger count outputs
        self.assertFalse(verify_project(project))

    def test_verify_project_paired_end_cellranger_count_and_multi_multiple_samples(self):
        """
        verify_project: paired-end data with cellranger 'count' and 'multi' (multiple samples)
        """
        analysis_dir = self._make_analysis_project(
            protocol="10x_CellPlex",
            sample_names=("PJB1_GEX", "PJB1_CML", "PJB2_GEX", "PJB2_CML"),
            seq_data_samples=("PJB1_GEX", "PJB2_GEX"),
            include_cellranger_multi=True,
            include_cellranger_count=True,
            cellranger_pipelines=('cellranger',),
            # NB only GEX samples
            cellranger_samples=("PJB1_GEX", "PJB2_GEX"),
            cellranger_multi_samples={ "PJB1": ("PB1", "PB2"),
                                       "PJB2": ("PB3", "PB4") })
        project = AnalysisProject(analysis_dir)
        self.assertTrue(verify_project(project))

    def test_verify_project_paired_end_cellranger_count_multiome(self):
        """
        verify_project: paired-end data with cellranger 'count' (multiome)
        """
        analysis_dir = self._make_analysis_project(
            protocol='10x_Multiome_GEX',
            include_cellranger_count=True,
            cellranger_pipelines=('cellranger-arc',),
            cellranger_samples=('PJB1','PJB2',))
        project = AnalysisProject(analysis_dir)
        self.assertFalse(verify_project(project))

    def test_verify_project_paired_end_cellranger_count_multiome_and_scrnaseq(self):
        """
        verify_project: paired-end data with cellranger 'count' (multiome+scRNAseq)
        """
        analysis_dir = self._make_analysis_project(
            protocol='10x_Multiome_GEX',
            include_cellranger_count=True,
            cellranger_pipelines=('cellranger',
                                  'cellranger-arc',),
            cellranger_samples=('PJB1','PJB2',))
        project = AnalysisProject(analysis_dir)
        self.assertTrue(verify_project(project))

    def test_verify_project_paired_end_legacy_cellranger_count(self):
        """
        verify_project: paired-end data with cellranger 'count' (legacy)
        """
        analysis_dir = self._make_analysis_project(
            protocol='10x_scRNAseq',
            paired_end=True,
            include_cellranger_count=True,
            cellranger_pipelines=('cellranger',),
            cellranger_samples=('PJB1','PJB2',),
            legacy_cellranger_outs=True)
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

    def test_verify_project_using_fastq_list_from_metadata(self):
        """
        verify_project: use implicit Fastq list from QC metadata
        """
        analysis_dir = self._make_analysis_project(
            protocol="standardPE",
            fastq_names=
            ("PJB1_S1_R1_001_paired.fastq.gz",
             "PJB1_S1_R2_001_paired.fastq.gz",
             "PJB2_S2_R1_001_paired.fastq.gz",
             "PJB2_S2_R2_001_paired.fastq.gz",))
        project = AnalysisProject(analysis_dir)
        # Remove some QC outputs from sample 'PJB1'
        qc_dir = os.path.join(analysis_dir,"qc")
        for f in os.listdir(qc_dir):
            if f.startswith("PJB1_"):
                ff = os.path.join(qc_dir,f)
                print("Removing %s" % ff)
                if os.path.isfile(ff):
                    os.remove(ff)
        # Rewrite the QC info to update the stored Fastqs list
        qc_info_file = os.path.join(analysis_dir,"qc","qc.info")
        with open(qc_info_file,'rt') as fp:
            qc_info = fp.read()
        with open(qc_info_file,'wt') as fp:
            for line in qc_info.split('\n'):
                # Update the "Fastqs" data item
                if line.startswith("Fastqs\t"):
                    line = "Fastqs\t"\
                           "PJB2_S2_R1_001_paired.fastq.gz,"\
                           "PJB2_S2_R2_001_paired.fastq.gz"
                if line:
                    fp.write("%s\n" % line)
        # Fails if full set of Fastqs is explcitly specified
        self.assertFalse(verify_project(project,
                                        fastqs=
                                        ["PJB1_S1_R1_001_paired.fastq.gz",
                                         "PJB1_S1_R2_001_paired.fastq.gz",
                                         "PJB2_S2_R1_001_paired.fastq.gz",
                                         "PJB2_S2_R2_001_paired.fastq.gz",]))
        # OK when using implicit Fastq list
        self.assertTrue(verify_project(project))

    def test_verify_project_using_fastq_list(self):
        """
        verify_project: use explicit Fastq list
        """
        analysis_dir = self._make_analysis_project(
            protocol="standardPE",
            fastq_names=
            ("PJB1_S1_R1_001_paired.fastq.gz",
             "PJB1_S1_R2_001_paired.fastq.gz",
             "PJB2_S2_R1_001_paired.fastq.gz",
             "PJB2_S2_R2_001_paired.fastq.gz",))
        project = AnalysisProject(analysis_dir)
        # Remove some QC outputs from sample 'PJB1'
        qc_dir = os.path.join(analysis_dir,"qc")
        for f in os.listdir(qc_dir):
            if f.startswith("PJB1_"):
                ff = os.path.join(qc_dir,f)
                print("Removing %s" % ff)
                if os.path.isfile(ff):
                    os.remove(ff)
        # Fails for default (all Fastqs in project)
        self.assertFalse(verify_project(project))
        # OK for subset specified explicitly
        self.assertTrue(verify_project(project,
                                       fastqs=
                                       ["PJB2_S2_R1_001_paired.fastq.gz",
                                        "PJB2_S2_R2_001_paired.fastq.gz",]))
