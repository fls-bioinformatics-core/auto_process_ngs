#######################################################################
# Unit tests for qc/outputs.py
#######################################################################

import unittest
import tempfile
import shutil
import os
from auto_process_ngs.mock import MockAnalysisProject
from auto_process_ngs.mock import UpdateAnalysisProject
from auto_process_ngs.mockqc import make_mock_qc_dir
from auto_process_ngs.mockqc import MockQCOutputs
from auto_process_ngs.analysis import AnalysisProject
from auto_process_ngs.qc.outputs import QCOutputs
from auto_process_ngs.qc.outputs import ExtraOutputs

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

class TestQCOutputs(unittest.TestCase):
    def setUp(self):
        # Temporary working dir (if needed)
        self.wd = None
    def tearDown(self):
        # Remove temporary working dir
        if not REMOVE_TEST_OUTPUTS:
            return
        if self.wd is not None and os.path.isdir(self.wd):
            shutil.rmtree(self.wd)
    def _make_qc_dir(self,qc_dir,fastq_names,
                     project_name='PJB',
                     organisms=('human',),
                     screens=('model_organisms','other_organisms','rRNA',),
                     cellranger_pipelines=('cellranger',),
                     cellranger_samples=None,
                     cellranger_multi_samples=None,
                     seq_data_samples=None,
                     include_fastqc=True,
                     include_fastq_screen=True,
                     include_strandedness=True,
                     include_seqlens=True,
                     include_rseqc_infer_experiment=False,
                     include_rseqc_genebody_coverage=False,
                     include_picard_insert_size_metrics=False,
                     include_qualimap_rnaseq=False,
                     include_multiqc=True,
                     include_cellranger_count=False,
                     include_cellranger_multi=False,
                     cellranger_version=None,
                     legacy_screens=False,
                     legacy_cellranger_outs=False,
                     protocol=None):
        # Create working directory and qc dir
        if self.wd is None:
            self.wd = tempfile.mkdtemp(suffix='.test_QCOutputs')
        return make_mock_qc_dir(
            os.path.join(self.wd,qc_dir),
            fastq_names,
            protocol=protocol,
            screens=screens,
            cellranger_pipelines=cellranger_pipelines,
            cellranger_samples=cellranger_samples,
            cellranger_multi_samples=cellranger_multi_samples,
            seq_data_samples=seq_data_samples,
            include_fastqc=include_fastqc,
            include_fastq_screen=include_fastq_screen,
            include_strandedness=include_strandedness,
            include_seqlens=include_seqlens,
            include_rseqc_infer_experiment=include_rseqc_infer_experiment,
            include_rseqc_genebody_coverage=include_rseqc_genebody_coverage,
            include_picard_insert_size_metrics=\
            include_picard_insert_size_metrics,
            include_qualimap_rnaseq=include_qualimap_rnaseq,
            include_multiqc=include_multiqc,
            include_cellranger_count=include_cellranger_count,
            include_cellranger_multi=include_cellranger_multi,
            cellranger_version=cellranger_version,
            legacy_screens=legacy_screens,
            legacy_cellranger_outs=legacy_cellranger_outs)

    def test_qcoutputs_single_end(self):
        """
	QCOutputs: single-end data
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB2_S2_R1_001',
                                   ))
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['fastqc_r1',
                          'multiqc',
                          'screens_r1',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB2_S2_R1_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.seq_data_samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.bams,[])
        self.assertEqual(qc_outputs.organisms,[])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.physical_samples,[])
        self.assertEqual(qc_outputs.reads,['r1'])
        self.assertEqual(qc_outputs.software,
                         { 'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.config_files,
                         ['fastq_strand.conf'])

    def test_qcoutputs_single_end_no_fastqc(self):
        """
	QCOutputs: single-end data (no FastQC)
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB2_S2_R1_001',
                                   ),
                                   include_fastqc=False)
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['multiqc',
                          'screens_r1',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB2_S2_R1_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.seq_data_samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.bams,[])
        self.assertEqual(qc_outputs.organisms,[])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.physical_samples,[])
        self.assertEqual(qc_outputs.reads,['r1'])
        self.assertEqual(qc_outputs.software,
                         { 'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.config_files,
                         ['fastq_strand.conf'])

    def test_qcoutputs_single_end_no_screens(self):
        """
	QCOutputs: single-end data (no screens)
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB2_S2_R1_001',
                                   ),
                                   include_fastq_screen=False)
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['fastqc_r1',
                          'multiqc',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB2_S2_R1_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.seq_data_samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.bams,[])
        self.assertEqual(qc_outputs.organisms,[])
        self.assertEqual(qc_outputs.fastq_screens,[])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.physical_samples,[])
        self.assertEqual(qc_outputs.reads,['r1'])
        self.assertEqual(qc_outputs.software,
                         { 'fastqc': [ '0.11.3' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.config_files,
                         ['fastq_strand.conf'])

    def test_qcoutputs_single_end_no_seqlens(self):
        """
	QCOutputs: single-end data (no sequence lengths)
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB2_S2_R1_001',
                                   ),
                                   include_seqlens=False)
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['fastqc_r1',
                          'multiqc',
                          'screens_r1',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB2_S2_R1_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.seq_data_samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.bams,[])
        self.assertEqual(qc_outputs.organisms,[])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.physical_samples,[])
        self.assertEqual(qc_outputs.reads,['r1'])
        self.assertEqual(qc_outputs.software,
                         { 'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,None)
        self.assertEqual(qc_outputs.stats.min_sequence_length,None)
        self.assertEqual(qc_outputs.stats.max_sequence_length,None)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         [])
        self.assertEqual(qc_outputs.config_files,
                         ['fastq_strand.conf'])

    def test_qcoutputs_single_end_no_standedness(self):
        """
	QCOutputs: single-end data (no strandedness)
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB2_S2_R1_001',
                                   ),
                                   include_strandedness=False)
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['fastqc_r1',
                          'multiqc',
                          'screens_r1',
                          'sequence_lengths'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB2_S2_R1_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.seq_data_samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.bams,[])
        self.assertEqual(qc_outputs.organisms,[])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.physical_samples,[])
        self.assertEqual(qc_outputs.reads,['r1'])
        self.assertEqual(qc_outputs.software,
                         { 'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.config_files,[])

    def test_qcoutputs_single_end_no_multiqc(self):
        """
	QCOutputs: single-end data (no MultiQC)
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB2_S2_R1_001',
                                   ),
                                   include_multiqc=False)
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['fastqc_r1',
                          'screens_r1',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB2_S2_R1_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.seq_data_samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.bams,[])
        self.assertEqual(qc_outputs.organisms,[])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.physical_samples,[])
        self.assertEqual(qc_outputs.reads,['r1'])
        self.assertEqual(qc_outputs.software,
                         { 'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.config_files,
                         ['fastq_strand.conf'])

    def test_qcoutputs_single_end_legacy_screen_naming(self):
        """
	QCOutputs: single-end data (legacy screen naming)
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB2_S2_R1_001',
                                   ),
                                   legacy_screens=True)
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['fastqc_r1',
                          'multiqc',
                          'screens_r1',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB2_S2_R1_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.seq_data_samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.bams,[])
        self.assertEqual(qc_outputs.organisms,[])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.physical_samples,[])
        self.assertEqual(qc_outputs.reads,['r1'])
        self.assertEqual(qc_outputs.software,
                         { 'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.config_files,
                         ['fastq_strand.conf'])

    def test_qcoutputs_paired_end(self):
        """
	QCOutputs: paired-end data
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB1_S1_R2_001',
                                       'PJB2_S2_R1_001',
                                       'PJB2_S2_R2_001',
                                   ))
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['fastqc_r1',
                          'fastqc_r2',
                          'multiqc',
                          'screens_r1',
                          'screens_r2',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB1_S1_R2_001',
                          'PJB2_S2_R1_001',
                          'PJB2_S2_R2_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.seq_data_samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.bams,[])
        self.assertEqual(qc_outputs.organisms,[])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.physical_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r2'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r2'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r2'],76)
        self.assertEqual(qc_outputs.config_files,
                         ['fastq_strand.conf'])

    def test_qcoutputs_paired_end_no_fastqc(self):
        """
	QCOutputs: paired-end data (no FastQC)
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB1_S1_R2_001',
                                       'PJB2_S2_R1_001',
                                       'PJB2_S2_R2_001',
                                   ),
                                   include_fastqc=False)
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['multiqc',
                          'screens_r1',
                          'screens_r2',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB1_S1_R2_001',
                          'PJB2_S2_R1_001',
                          'PJB2_S2_R2_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.seq_data_samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.bams,[])
        self.assertEqual(qc_outputs.organisms,[])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.physical_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r2'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r2'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r2'],76)
        self.assertEqual(qc_outputs.config_files,
                         ['fastq_strand.conf'])

    def test_qcoutputs_paired_end_rseqc_infer_experiment(self):
        """
	QCOutputs: paired-end data with RSeQC infer_experiment.py output
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB1_S1_R2_001',
                                       'PJB2_S2_R1_001',
                                       'PJB2_S2_R2_001',
                                   ),
                                   include_rseqc_infer_experiment=True)
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['fastqc_r1',
                          'fastqc_r2',
                          'multiqc',
                          'rseqc_infer_experiment',
                          'screens_r1',
                          'screens_r2',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB1_S1_R2_001',
                          'PJB2_S2_R1_001',
                          'PJB2_S2_R2_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.seq_data_samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.bams,['PJB1_S1_001',
                                          'PJB2_S2_001'])
        self.assertEqual(qc_outputs.organisms,['human'])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.physical_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                           'rseqc:infer_experiment': [ '4.0.0' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r2'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r2'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r2'],76)
        self.assertEqual(qc_outputs.config_files,
                         ['fastq_strand.conf'])

    def test_qcoutputs_paired_end_rseqc_genebody_coverage(self):
        """
	QCOutputs: paired-end data with RSeQC gene body coverage
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB1_S1_R2_001',
                                       'PJB2_S2_R1_001',
                                       'PJB2_S2_R2_001',
                                   ),
                                   include_rseqc_genebody_coverage=True)
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['fastqc_r1',
                          'fastqc_r2',
                          'multiqc',
                          'rseqc_genebody_coverage',
                          'screens_r1',
                          'screens_r2',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB1_S1_R2_001',
                          'PJB2_S2_R1_001',
                          'PJB2_S2_R2_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.seq_data_samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.bams,[])
        self.assertEqual(qc_outputs.organisms,['human'])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.physical_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                           'rseqc:genebody_coverage': [ '4.0.0' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r2'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r2'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r2'],76)
        self.assertEqual(qc_outputs.config_files,
                         ['fastq_strand.conf'])

    def test_qcoutputs_paired_end_picard_insert_size_metrics(self):
        """
	QCOutputs: paired-end data with Picard insert size metrics
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB1_S1_R2_001',
                                       'PJB2_S2_R1_001',
                                       'PJB2_S2_R2_001',
                                   ),
                                   include_picard_insert_size_metrics=True)
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['collated_insert_sizes',
                          'fastqc_r1',
                          'fastqc_r2',
                          'multiqc',
                          'picard_insert_size_metrics',
                          'screens_r1',
                          'screens_r2',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB1_S1_R2_001',
                          'PJB2_S2_R1_001',
                          'PJB2_S2_R2_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.seq_data_samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.bams,['PJB1_S1_001',
                                          'PJB2_S2_001'])
        self.assertEqual(qc_outputs.organisms,['human'])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.physical_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                           'picard': [ '2.27.1' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r2'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r2'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r2'],76)
        self.assertEqual(qc_outputs.config_files,
                         ['fastq_strand.conf'])

    def test_qcoutputs_paired_end_qualimap_rnaseq(self):
        """
	QCOutputs: paired-end data with Qualimap 'rnaseq' outputs
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB1_S1_R2_001',
                                       'PJB2_S2_R1_001',
                                       'PJB2_S2_R2_001',
                                   ),
                                   include_qualimap_rnaseq=True)
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['fastqc_r1',
                          'fastqc_r2',
                          'multiqc',
                          'qualimap_rnaseq',
                          'screens_r1',
                          'screens_r2',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB1_S1_R2_001',
                          'PJB2_S2_R1_001',
                          'PJB2_S2_R2_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.seq_data_samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.bams,['PJB1_S1_001',
                                          'PJB2_S2_001'])
        self.assertEqual(qc_outputs.organisms,['human'])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.physical_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                           'qualimap': [ 'v.2.2.2' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r2'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r2'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r2'],76)
        self.assertEqual(qc_outputs.config_files,
                         ['fastq_strand.conf'])

    def test_qcoutputs_paired_end_no_screens(self):
        """
	QCOutputs: paired-end data (no screens)
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB1_S1_R2_001',
                                       'PJB2_S2_R1_001',
                                       'PJB2_S2_R2_001',
                                   ),
                                   include_fastq_screen=False)
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['fastqc_r1',
                          'fastqc_r2',
                          'multiqc',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB1_S1_R2_001',
                          'PJB2_S2_R1_001',
                          'PJB2_S2_R2_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.seq_data_samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.bams,[])
        self.assertEqual(qc_outputs.organisms,[])
        self.assertEqual(qc_outputs.fastq_screens,[])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.physical_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'fastqc': [ '0.11.3' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r2'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r2'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r2'],76)
        self.assertEqual(qc_outputs.config_files,
                         ['fastq_strand.conf'])

    def test_qcoutputs_paired_end_no_seqlens(self):
        """
	QCOutputs: paired-end data (no sequence lengths)
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB1_S1_R2_001',
                                       'PJB2_S2_R1_001',
                                       'PJB2_S2_R2_001',
                                   ),
                                   include_seqlens=False)
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['fastqc_r1',
                          'fastqc_r2',
                          'multiqc',
                          'screens_r1',
                          'screens_r2',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB1_S1_R2_001',
                          'PJB2_S2_R1_001',
                          'PJB2_S2_R2_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.seq_data_samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.bams,[])
        self.assertEqual(qc_outputs.organisms,[])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.physical_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,None)
        self.assertEqual(qc_outputs.stats.min_sequence_length,None)
        self.assertEqual(qc_outputs.stats.max_sequence_length,None)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),[])
        self.assertEqual(qc_outputs.config_files,
                         ['fastq_strand.conf'])

    def test_qcoutputs_paired_end_no_strandedness(self):
        """
	QCOutputs: paired-end data (no strandedness)
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB1_S1_R2_001',
                                       'PJB2_S2_R1_001',
                                       'PJB2_S2_R2_001',
                                   ),
                                   include_strandedness=False)
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['fastqc_r1',
                          'fastqc_r2',
                          'multiqc',
                          'screens_r1',
                          'screens_r2',
                          'sequence_lengths'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB1_S1_R2_001',
                          'PJB2_S2_R1_001',
                          'PJB2_S2_R2_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.seq_data_samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.bams,[])
        self.assertEqual(qc_outputs.organisms,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.physical_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r2'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r2'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r2'],76)
        self.assertEqual(qc_outputs.config_files,[])

    def test_qcoutputs_paired_end_no_multiqc(self):
        """
	QCOutputs: paired-end data (no MultiQC)
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB1_S1_R2_001',
                                       'PJB2_S2_R1_001',
                                       'PJB2_S2_R2_001',
                                   ),
                                   include_multiqc=False)
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['fastqc_r1',
                          'fastqc_r2',
                          'screens_r1',
                          'screens_r2',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB1_S1_R2_001',
                          'PJB2_S2_R1_001',
                          'PJB2_S2_R2_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.seq_data_samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.bams,[])
        self.assertEqual(qc_outputs.organisms,[])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.physical_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r2'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r2'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r2'],76)
        self.assertEqual(qc_outputs.config_files,
                         ['fastq_strand.conf'])

    def test_qcoutputs_paired_end_legacy_screen_naming(self):
        """
	QCOutputs: paired-end data (legacy screen naming)
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB1_S1_R2_001',
                                       'PJB2_S2_R1_001',
                                       'PJB2_S2_R2_001',
                                   ),
                                   legacy_screens=True)
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['fastqc_r1',
                          'fastqc_r2',
                          'multiqc',
                          'screens_r1',
                          'screens_r2',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB1_S1_R2_001',
                          'PJB2_S2_R1_001',
                          'PJB2_S2_R2_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.seq_data_samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.bams,[])
        self.assertEqual(qc_outputs.organisms,[])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.physical_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r2'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r2'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r2'],76)
        self.assertEqual(qc_outputs.config_files,
                         ['fastq_strand.conf'])

    def test_qcoutputs_paired_end_explicit_biological_samples(self):
        """
	QCOutputs: paired-end data (biological samples explicitly specified)
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB1_S1_R2_001',
                                       'PJB2_S2_R1_001',
                                       'PJB2_S2_R2_001',
                                   ),
                                   seq_data_samples=('PJB2',))
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['fastqc_r1',
                          'fastqc_r2',
                          'multiqc',
                          'screens_r1',
                          'screens_r2',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB1_S1_R2_001',
                          'PJB2_S2_R1_001',
                          'PJB2_S2_R2_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.seq_data_samples,
                         ['PJB2'])
        self.assertEqual(qc_outputs.bams,[])
        self.assertEqual(qc_outputs.organisms,[])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.physical_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r2'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r2'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r2'],76)
        self.assertEqual(qc_outputs.config_files,
                         ['fastq_strand.conf'])

    def test_qcoutputs_paired_end_extra_outputs_tsv(self):
        """
        QCOutputs: paired-end data with 'extra_outputs.tsv'
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB1_S1_R2_001',
                                       'PJB2_S2_R1_001',
                                       'PJB2_S2_R2_001',
                                   ))
        extra_outputs_tsv = os.path.join(qc_dir,"extra_outputs.tsv")
        with open(extra_outputs_tsv,'wt') as fp:
            fp.write("# Extra files to include in QC reporting\nanalyser/index.html\tReport from 'analyser'\nfinal_result/main.html\tFinal result\tfinal_result/files")
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['extra_outputs',
                          'fastqc_r1',
                          'fastqc_r2',
                          'multiqc',
                          'screens_r1',
                          'screens_r2',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB1_S1_R2_001',
                          'PJB2_S2_R1_001',
                          'PJB2_S2_R2_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.seq_data_samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.bams,[])
        self.assertEqual(qc_outputs.organisms,[])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.physical_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r2'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r2'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r2'],76)
        self.assertEqual(qc_outputs.config_files,
                         ['fastq_strand.conf'])
        for f in ("analyser/index.html",
                  "final_result/main.html",
                  "final_result/files"):
            self.assertTrue(os.path.join(qc_dir,f)
                            in qc_outputs.output_files)

    def test_qcoutputs_10x_cellranger_count_612(self):
        """
        QCOutputs: 10xGenomics data with cellranger 'count' (6.1.2)
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB1_S1_R2_001',
                                       'PJB2_S2_R1_001',
                                       'PJB2_S2_R2_001',
                                   ),
                                   include_cellranger_count=True,
                                   cellranger_pipelines=('cellranger',),
                                   cellranger_samples=(
                                       'PJB1',
                                       'PJB2',
                                   ),
                                   cellranger_version='6.1.2')
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['cellranger_count',
                          'fastqc_r1',
                          'fastqc_r2',
                          'multiqc',
                          'screens_r1',
                          'screens_r2',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB1_S1_R2_001',
                          'PJB2_S2_R1_001',
                          'PJB2_S2_R2_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.seq_data_samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.bams,[])
        self.assertEqual(qc_outputs.organisms,[])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,
                         ['/data/refdata-cellranger-2020-A'])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.physical_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'cellranger': [ '6.1.2' ],
                           'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r2'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r2'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r2'],76)
        self.assertEqual(qc_outputs.config_files,
                         ['fastq_strand.conf'])

    def test_qcoutputs_10x_cellranger_count_710(self):
        """
        QCOutputs: 10xGenomics data with cellranger 'count' (7.1.2)
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB1_S1_R2_001',
                                       'PJB2_S2_R1_001',
                                       'PJB2_S2_R2_001',
                                   ),
                                   include_cellranger_count=True,
                                   cellranger_pipelines=('cellranger',),
                                   cellranger_samples=(
                                       'PJB1',
                                       'PJB2',
                                   ),
                                   cellranger_version='7.1.0')
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['cellranger_count',
                          'fastqc_r1',
                          'fastqc_r2',
                          'multiqc',
                          'screens_r1',
                          'screens_r2',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB1_S1_R2_001',
                          'PJB2_S2_R1_001',
                          'PJB2_S2_R2_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.seq_data_samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.bams,[])
        self.assertEqual(qc_outputs.organisms,[])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,
                         ['/data/refdata-cellranger-2020-A'])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.physical_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'cellranger': [ '7.1.0' ],
                           'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r2'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r2'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r2'],76)
        self.assertEqual(qc_outputs.config_files,
                         ['fastq_strand.conf'])

    def test_qcoutputs_10x_cellranger_count_800(self):
        """
        QCOutputs: 10xGenomics data with cellranger 'count' (8.0.0)
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB1_S1_R2_001',
                                       'PJB2_S2_R1_001',
                                       'PJB2_S2_R2_001',
                                   ),
                                   include_cellranger_count=True,
                                   cellranger_pipelines=('cellranger',),
                                   cellranger_samples=(
                                       'PJB1',
                                       'PJB2',
                                   ),
                                   cellranger_version='8.0.0')
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['cellranger_count',
                          'fastqc_r1',
                          'fastqc_r2',
                          'multiqc',
                          'screens_r1',
                          'screens_r2',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB1_S1_R2_001',
                          'PJB2_S2_R1_001',
                          'PJB2_S2_R2_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.seq_data_samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.bams,[])
        self.assertEqual(qc_outputs.organisms,[])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,
                         ['/data/refdata-cellranger-2020-A'])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.physical_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'cellranger': [ '8.0.0' ],
                           'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r2'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r2'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r2'],76)
        self.assertEqual(qc_outputs.config_files,
                         ['fastq_strand.conf'])

    def test_qcoutputs_10x_cellranger_count_900(self):
        """
        QCOutputs: 10xGenomics data with cellranger 'count' (9.0.0)
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB1_S1_R2_001',
                                       'PJB2_S2_R1_001',
                                       'PJB2_S2_R2_001',
                                   ),
                                   include_cellranger_count=True,
                                   cellranger_pipelines=('cellranger',),
                                   cellranger_samples=(
                                       'PJB1',
                                       'PJB2',
                                   ),
                                   cellranger_version='9.0.0')
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['cellranger_count',
                          'fastqc_r1',
                          'fastqc_r2',
                          'multiqc',
                          'screens_r1',
                          'screens_r2',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB1_S1_R2_001',
                          'PJB2_S2_R1_001',
                          'PJB2_S2_R2_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.seq_data_samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.bams,[])
        self.assertEqual(qc_outputs.organisms,[])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,
                         ['/data/refdata-cellranger-2020-A'])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.physical_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'cellranger': [ '9.0.0' ],
                           'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r2'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r2'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r2'],76)
        self.assertEqual(qc_outputs.config_files,
                         ['fastq_strand.conf'])

    def test_qcoutputs_10x_cellranger_count_legacy(self):
        """
        QCOutputs: 10xGenomics data with cellranger 'count' (legacy format)
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB1_S1_R2_001',
                                       'PJB2_S2_R1_001',
                                       'PJB2_S2_R2_001',
                                   ),
                                   include_cellranger_count=True,
                                   cellranger_pipelines=('cellranger',),
                                   cellranger_samples=(
                                       'PJB1',
                                       'PJB2',
                                   ),
                                   legacy_cellranger_outs=True)
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['cellranger_count',
                          'fastqc_r1',
                          'fastqc_r2',
                          'multiqc',
                          'screens_r1',
                          'screens_r2',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB1_S1_R2_001',
                          'PJB2_S2_R1_001',
                          'PJB2_S2_R2_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.seq_data_samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.bams,[])
        self.assertEqual(qc_outputs.organisms,[])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,
                         ['/data/refdata-cellranger-2020-A'])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.physical_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'cellranger': [ '?' ],
                           'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r2'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r2'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r2'],76)
        self.assertEqual(qc_outputs.config_files,
                         ['fastq_strand.conf'])

    def test_qcoutputs_10x_cellranger_atac_count(self):
        """
        QCOutputs: 10xGenomics data with cellranger-atac 'count'
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB1_S1_R3_001',
                                       'PJB2_S2_R1_001',
                                       'PJB2_S2_R3_001',
                                   ),
                                   include_cellranger_count=True,
                                   cellranger_pipelines=('cellranger-atac',),
                                   cellranger_samples=(
                                       'PJB1',
                                       'PJB2',
                                   ))
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['cellranger-atac_count',
                          'fastqc_r1',
                          'fastqc_r3',
                          'multiqc',
                          'screens_r1',
                          'screens_r3',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB1_S1_R3_001',
                          'PJB2_S2_R1_001',
                          'PJB2_S2_R3_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.seq_data_samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.bams,[])
        self.assertEqual(qc_outputs.organisms,[])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,
                         ['/data/refdata-cellranger-atac-2020-A'])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.physical_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r3'])
        self.assertEqual(qc_outputs.software,
                         { 'cellranger-atac': [ '2.0.0' ],
                           'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r3'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r3'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r3'],76)
        self.assertEqual(qc_outputs.config_files,
                         ['fastq_strand.conf'])

    def test_qcoutputs_10x_multiome_gex(self):
        """
        QCOutputs: 10xGenomics multiome data (GEX component)
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB1_S1_R2_001',
                                       'PJB2_S2_R1_001',
                                       'PJB2_S2_R2_001',
                                   ),
                                   include_cellranger_count=True,
                                   cellranger_pipelines=('cellranger',
                                                         'cellranger-arc',),
                                   cellranger_samples=(
                                       'PJB1',
                                       'PJB2',
                                   ))
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['cellranger-arc_count',
                          'cellranger_count',
                          'fastqc_r1',
                          'fastqc_r2',
                          'multiqc',
                          'screens_r1',
                          'screens_r2',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB1_S1_R2_001',
                          'PJB2_S2_R1_001',
                          'PJB2_S2_R2_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.seq_data_samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.bams,[])
        self.assertEqual(qc_outputs.organisms,[])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,
                         ['/data/refdata-cellranger-2020-A',
                          '/data/refdata-cellranger-arc-2020-A'])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.physical_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'cellranger': [ '8.0.0' ],
                           'cellranger-arc': [ '2.0.0' ],
                           'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r2'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r2'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r2'],76)
        self.assertEqual(qc_outputs.config_files,
                         ['fastq_strand.conf',
                          'libraries.PJB1.csv',
                          'libraries.PJB2.csv'])

    def test_qcoutputs_10x_multiome_atac(self):
        """
        QCOutputs: 10xGenomics multiome data (ATAC component)
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB1_S1_R3_001',
                                       'PJB2_S2_R1_001',
                                       'PJB2_S2_R3_001',
                                   ),
                                   include_cellranger_count=True,
                                   cellranger_pipelines=('cellranger-arc',
                                                         'cellranger-atac',),
                                   cellranger_samples=(
                                       'PJB1',
                                       'PJB2',
                                   ))
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['cellranger-arc_count',
                          'cellranger-atac_count',
                          'fastqc_r1',
                          'fastqc_r3',
                          'multiqc',
                          'screens_r1',
                          'screens_r3',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB1_S1_R3_001',
                          'PJB2_S2_R1_001',
                          'PJB2_S2_R3_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.seq_data_samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.bams,[])
        self.assertEqual(qc_outputs.organisms,[])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,
                         ['/data/refdata-cellranger-arc-2020-A',
                          '/data/refdata-cellranger-atac-2020-A'])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.physical_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r3'])
        self.assertEqual(qc_outputs.software,
                         { 'cellranger-arc': [ '2.0.0' ],
                           'cellranger-atac': [ '2.0.0' ],
                           'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r3'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r3'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r3'],76)
        self.assertEqual(qc_outputs.config_files,
                         ['fastq_strand.conf',
                          'libraries.PJB1.csv',
                          'libraries.PJB2.csv'])

    def test_qcoutputs_10x_cellranger_multi_cellplex_612(self):
        """
        QCOutputs: 10xGenomics CellPlex data with 'cellranger multi' (6.1.2)
        """
        qc_dir = self._make_qc_dir('qc',
                                   protocol="10x_CellPlex",
                                   fastq_names=(
                                       'PJB1_GEX_S1_R1_001',
                                       'PJB1_GEX_S1_R2_001',
                                       'PJB2_MC_S2_R1_001',
                                       'PJB2_MC_S2_R2_001',
                                   ),
                                   include_cellranger_multi=True,
                                   cellranger_pipelines=('cellranger',),
                                   cellranger_samples=(
                                       'PJB1_GEX',
                                       'PJB2_MC',
                                   ),
                                   cellranger_multi_samples=(
                                       'PJB_CML1',
                                       'PJB_CML2',
                                   ),
                                   cellranger_version="6.1.2")
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['cellranger_multi',
                          'fastqc_r1',
                          'fastqc_r2',
                          'multiqc',
                          'screens_r2',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_GEX_S1_R1_001',
                          'PJB1_GEX_S1_R2_001',
                          'PJB2_MC_S2_R1_001',
                          'PJB2_MC_S2_R2_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1_GEX','PJB2_MC'])
        self.assertEqual(qc_outputs.seq_data_samples,
                         ['PJB1_GEX'])
        self.assertEqual(qc_outputs.bams,[])
        self.assertEqual(qc_outputs.organisms,[])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,
                         ['/data/refdata-cellranger-2020-A'])
        self.assertEqual(qc_outputs.cellranger_probe_sets,[])
        self.assertEqual(qc_outputs.multiplexed_samples,
                         ['PJB_CML1','PJB_CML2'])
        self.assertEqual(qc_outputs.physical_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'cellranger': [ '6.1.2' ],
                           'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r2'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r2'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r2'],76)
        self.assertEqual(qc_outputs.config_files,
                         ['10x_multi_config.csv',
                          'fastq_strand.conf'])

    def test_qcoutputs_10x_cellranger_multi_cellplex_710(self):
        """
        QCOutputs: 10xGenomics CellPlex data with 'cellranger multi' (7.1.0)
        """
        qc_dir = self._make_qc_dir('qc',
                                   protocol="10x_CellPlex",
                                   fastq_names=(
                                       'PJB1_GEX_S1_R1_001',
                                       'PJB1_GEX_S1_R2_001',
                                       'PJB2_MC_S2_R1_001',
                                       'PJB2_MC_S2_R2_001',
                                   ),
                                   include_cellranger_multi=True,
                                   cellranger_pipelines=('cellranger',),
                                   cellranger_samples=(
                                       'PJB1_GEX',
                                       'PJB2_MC',
                                   ),
                                   cellranger_multi_samples=(
                                       'PJB_CML1',
                                       'PJB_CML2',
                                   ),
                                   cellranger_version="7.1.0")
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['cellranger_multi',
                          'fastqc_r1',
                          'fastqc_r2',
                          'multiqc',
                          'screens_r2',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_GEX_S1_R1_001',
                          'PJB1_GEX_S1_R2_001',
                          'PJB2_MC_S2_R1_001',
                          'PJB2_MC_S2_R2_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1_GEX','PJB2_MC'])
        self.assertEqual(qc_outputs.seq_data_samples,
                         ['PJB1_GEX'])
        self.assertEqual(qc_outputs.bams,[])
        self.assertEqual(qc_outputs.organisms,[])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,
                         ['/data/refdata-cellranger-2020-A'])
        self.assertEqual(qc_outputs.cellranger_probe_sets,[])
        self.assertEqual(qc_outputs.multiplexed_samples,
                         ['PJB_CML1','PJB_CML2'])
        self.assertEqual(qc_outputs.physical_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'cellranger': [ '7.1.0' ],
                           'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r2'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r2'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r2'],76)
        self.assertEqual(qc_outputs.config_files,
                         ['10x_multi_config.csv',
                          'fastq_strand.conf'])

    def test_qcoutputs_10x_cellranger_multi_cellplex_800(self):
        """
        QCOutputs: 10xGenomics CellPlex data with 'cellranger multi' (8.0.0)
        """
        qc_dir = self._make_qc_dir('qc',
                                   protocol="10x_CellPlex",
                                   fastq_names=(
                                       'PJB1_GEX_S1_R1_001',
                                       'PJB1_GEX_S1_R2_001',
                                       'PJB2_MC_S2_R1_001',
                                       'PJB2_MC_S2_R2_001',
                                   ),
                                   include_cellranger_multi=True,
                                   cellranger_pipelines=('cellranger',),
                                   cellranger_samples=(
                                       'PJB1_GEX',
                                       'PJB2_MC',
                                   ),
                                   cellranger_multi_samples=(
                                       'PJB_CML1',
                                       'PJB_CML2',
                                   ),
                                   cellranger_version="8.0.0")
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['cellranger_multi',
                          'fastqc_r1',
                          'fastqc_r2',
                          'multiqc',
                          'screens_r2',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_GEX_S1_R1_001',
                          'PJB1_GEX_S1_R2_001',
                          'PJB2_MC_S2_R1_001',
                          'PJB2_MC_S2_R2_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1_GEX','PJB2_MC'])
        self.assertEqual(qc_outputs.seq_data_samples,
                         ['PJB1_GEX'])
        self.assertEqual(qc_outputs.bams,[])
        self.assertEqual(qc_outputs.organisms,[])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,
                         ['/data/refdata-cellranger-2020-A'])
        self.assertEqual(qc_outputs.cellranger_probe_sets,[])
        self.assertEqual(qc_outputs.multiplexed_samples,
                         ['PJB_CML1','PJB_CML2'])
        self.assertEqual(qc_outputs.physical_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'cellranger': [ '8.0.0' ],
                           'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r2'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r2'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r2'],76)
        self.assertEqual(qc_outputs.config_files,
                         ['10x_multi_config.csv',
                          'fastq_strand.conf'])

    def test_qcoutputs_10x_cellranger_multi_cellplex_900(self):
        """
        QCOutputs: 10xGenomics CellPlex data with 'cellranger multi' (9.0.0)
        """
        qc_dir = self._make_qc_dir('qc',
                                   protocol="10x_CellPlex",
                                   fastq_names=(
                                       'PJB1_GEX_S1_R1_001',
                                       'PJB1_GEX_S1_R2_001',
                                       'PJB2_MC_S2_R1_001',
                                       'PJB2_MC_S2_R2_001',
                                   ),
                                   include_cellranger_multi=True,
                                   cellranger_pipelines=('cellranger',),
                                   cellranger_samples=(
                                       'PJB1_GEX',
                                       'PJB2_MC',
                                   ),
                                   cellranger_multi_samples=(
                                       'PJB_CML1',
                                       'PJB_CML2',
                                   ),
                                   cellranger_version="9.0.0")
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['cellranger_multi',
                          'fastqc_r1',
                          'fastqc_r2',
                          'multiqc',
                          'screens_r2',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_GEX_S1_R1_001',
                          'PJB1_GEX_S1_R2_001',
                          'PJB2_MC_S2_R1_001',
                          'PJB2_MC_S2_R2_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1_GEX','PJB2_MC'])
        self.assertEqual(qc_outputs.seq_data_samples,
                         ['PJB1_GEX'])
        self.assertEqual(qc_outputs.bams,[])
        self.assertEqual(qc_outputs.organisms,[])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,
                         ['/data/refdata-cellranger-2020-A'])
        self.assertEqual(qc_outputs.cellranger_probe_sets,[])
        self.assertEqual(qc_outputs.multiplexed_samples,
                         ['PJB_CML1','PJB_CML2'])
        self.assertEqual(qc_outputs.physical_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'cellranger': [ '9.0.0' ],
                           'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r2'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r2'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r2'],76)
        self.assertEqual(qc_outputs.config_files,
                         ['10x_multi_config.csv',
                          'fastq_strand.conf'])

    def test_qcoutputs_10x_cellranger_multi_and_count_cellplex_612(self):
        """
        QCOutputs: 10xGenomics CellPlex data with 'cellranger multi' and 'count' (6.1.2)
        """
        qc_dir = self._make_qc_dir('qc',
                                   protocol="10x_CellPlex",
                                   fastq_names=(
                                       'PJB1_GEX_S1_R1_001',
                                       'PJB1_GEX_S1_R2_001',
                                       'PJB2_MC_S2_R1_001',
                                       'PJB2_MC_S2_R2_001',
                                   ),
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
                                   ),
                                   cellranger_version="6.1.2")
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['cellranger_count',
                          'cellranger_multi',
                          'fastqc_r1',
                          'fastqc_r2',
                          'multiqc',
                          'screens_r2',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_GEX_S1_R1_001',
                          'PJB1_GEX_S1_R2_001',
                          'PJB2_MC_S2_R1_001',
                          'PJB2_MC_S2_R2_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1_GEX','PJB2_MC'])
        self.assertEqual(qc_outputs.seq_data_samples,
                         ['PJB1_GEX'])
        self.assertEqual(qc_outputs.bams,[])
        self.assertEqual(qc_outputs.organisms,[])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,
                         ['/data/refdata-cellranger-2020-A'])
        self.assertEqual(qc_outputs.cellranger_probe_sets,[])
        self.assertEqual(qc_outputs.multiplexed_samples,
                         ['PJB_CML1','PJB_CML2'])
        self.assertEqual(qc_outputs.physical_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'cellranger': [ '6.1.2' ],
                           'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r2'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r2'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r2'],76)
        self.assertEqual(qc_outputs.config_files,
                         ['10x_multi_config.csv',
                          'fastq_strand.conf'])

    def test_qcoutputs_10x_cellranger_multi_and_count_cellplex_710(self):
        """
        QCOutputs: 10xGenomics CellPlex data with 'cellranger multi' and 'count' (7.1.0)
        """
        qc_dir = self._make_qc_dir('qc',
                                   protocol="10x_CellPlex",
                                   fastq_names=(
                                       'PJB1_GEX_S1_R1_001',
                                       'PJB1_GEX_S1_R2_001',
                                       'PJB2_MC_S2_R1_001',
                                       'PJB2_MC_S2_R2_001',
                                   ),
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
                                   ),
                                   cellranger_version="7.1.0")
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['cellranger_count',
                          'cellranger_multi',
                          'fastqc_r1',
                          'fastqc_r2',
                          'multiqc',
                          'screens_r2',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_GEX_S1_R1_001',
                          'PJB1_GEX_S1_R2_001',
                          'PJB2_MC_S2_R1_001',
                          'PJB2_MC_S2_R2_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1_GEX','PJB2_MC'])
        self.assertEqual(qc_outputs.seq_data_samples,
                         ['PJB1_GEX'])
        self.assertEqual(qc_outputs.bams,[])
        self.assertEqual(qc_outputs.organisms,[])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,
                         ['/data/refdata-cellranger-2020-A'])
        self.assertEqual(qc_outputs.cellranger_probe_sets,[])
        self.assertEqual(qc_outputs.multiplexed_samples,
                         ['PJB_CML1','PJB_CML2'])
        self.assertEqual(qc_outputs.physical_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'cellranger': [ '7.1.0' ],
                           'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r2'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r2'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r2'],76)
        self.assertEqual(qc_outputs.config_files,
                         ['10x_multi_config.csv',
                          'fastq_strand.conf'])

    def test_qcoutputs_10x_cellranger_multi_and_count_cellplex_800(self):
        """
        QCOutputs: 10xGenomics CellPlex data with 'cellranger multi' and 'count' (8.0.0)
        """
        qc_dir = self._make_qc_dir('qc',
                                   protocol="10x_CellPlex",
                                   fastq_names=(
                                       'PJB1_GEX_S1_R1_001',
                                       'PJB1_GEX_S1_R2_001',
                                       'PJB2_MC_S2_R1_001',
                                       'PJB2_MC_S2_R2_001',
                                   ),
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
                                   ),
                                   cellranger_version="8.0.0")
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['cellranger_count',
                          'cellranger_multi',
                          'fastqc_r1',
                          'fastqc_r2',
                          'multiqc',
                          'screens_r2',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_GEX_S1_R1_001',
                          'PJB1_GEX_S1_R2_001',
                          'PJB2_MC_S2_R1_001',
                          'PJB2_MC_S2_R2_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1_GEX','PJB2_MC'])
        self.assertEqual(qc_outputs.seq_data_samples,
                         ['PJB1_GEX'])
        self.assertEqual(qc_outputs.bams,[])
        self.assertEqual(qc_outputs.organisms,[])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,
                         ['/data/refdata-cellranger-2020-A'])
        self.assertEqual(qc_outputs.cellranger_probe_sets,[])
        self.assertEqual(qc_outputs.multiplexed_samples,
                         ['PJB_CML1','PJB_CML2'])
        self.assertEqual(qc_outputs.physical_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'cellranger': [ '8.0.0' ],
                           'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r2'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r2'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r2'],76)
        self.assertEqual(qc_outputs.config_files,
                         ['10x_multi_config.csv',
                          'fastq_strand.conf'])

    def test_qcoutputs_10x_cellranger_multi_and_count_cellplex_900(self):
        """
        QCOutputs: 10xGenomics CellPlex data with 'cellranger multi' and 'count' (9.0.0)
        """
        qc_dir = self._make_qc_dir('qc',
                                   protocol="10x_CellPlex",
                                   fastq_names=(
                                       'PJB1_GEX_S1_R1_001',
                                       'PJB1_GEX_S1_R2_001',
                                       'PJB2_MC_S2_R1_001',
                                       'PJB2_MC_S2_R2_001',
                                   ),
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
                                   ),
                                   cellranger_version="9.0.0")
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['cellranger_count',
                          'cellranger_multi',
                          'fastqc_r1',
                          'fastqc_r2',
                          'multiqc',
                          'screens_r2',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_GEX_S1_R1_001',
                          'PJB1_GEX_S1_R2_001',
                          'PJB2_MC_S2_R1_001',
                          'PJB2_MC_S2_R2_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1_GEX','PJB2_MC'])
        self.assertEqual(qc_outputs.seq_data_samples,
                         ['PJB1_GEX'])
        self.assertEqual(qc_outputs.bams,[])
        self.assertEqual(qc_outputs.organisms,[])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,
                         ['/data/refdata-cellranger-2020-A'])
        self.assertEqual(qc_outputs.cellranger_probe_sets,[])
        self.assertEqual(qc_outputs.multiplexed_samples,
                         ['PJB_CML1','PJB_CML2'])
        self.assertEqual(qc_outputs.physical_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'cellranger': [ '9.0.0' ],
                           'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r2'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r2'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r2'],76)
        self.assertEqual(qc_outputs.config_files,
                         ['10x_multi_config.csv',
                          'fastq_strand.conf'])

    def test_qcoutputs_10x_cellranger_multi_flex_612(self):
        """
        QCOutputs: 10xGenomics Flex data with 'cellranger multi' (6.1.2)
        """
        qc_dir = self._make_qc_dir('qc',
                                   protocol="10x_Flex",
                                   fastq_names=(
                                       'PJB1_flex_S1_R1_001',
                                       'PJB1_flex_S1_R2_001',
                                   ),
                                   include_cellranger_multi=True,
                                   include_cellranger_count=True,
                                   cellranger_pipelines=('cellranger',),
                                   cellranger_samples=(
                                       'PJB1_flex',
                                   ),
                                   cellranger_multi_samples=(
                                       'PJB_BC1',
                                       'PJB_BC2',
                                   ),
                                   cellranger_version="6.1.2")
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['cellranger_count',
                          'cellranger_multi',
                          'fastqc_r1',
                          'fastqc_r2',
                          'multiqc',
                          'screens_r2',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_flex_S1_R1_001',
                          'PJB1_flex_S1_R2_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1_flex',])
        self.assertEqual(qc_outputs.seq_data_samples,
                         ['PJB1_flex'])
        self.assertEqual(qc_outputs.bams,[])
        self.assertEqual(qc_outputs.organisms,[])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,
                         ['/data/refdata-cellranger-2020-A'])
        self.assertEqual(qc_outputs.cellranger_probe_sets,
                         ['/data/probe-set-2020-A'])
        self.assertEqual(qc_outputs.multiplexed_samples,
                         ['PJB_BC1','PJB_BC2'])
        self.assertEqual(qc_outputs.physical_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'cellranger': [ '6.1.2' ],
                           'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r2'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r2'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r2'],76)
        self.assertEqual(qc_outputs.config_files,
                         ['10x_multi_config.csv',
                          'fastq_strand.conf'])

    def test_qcoutputs_10x_cellranger_single_cell_immune_profiling_710(self):
        """
        QCOutputs: 10xGenomics single cell immune profiling (7.1.0)
        """
        qc_dir = self._make_qc_dir('qc',
                                   protocol="10x_ImmuneProfiling",
                                   fastq_names=(
                                       'PJB1_GEX_S1_R1_001',
                                       'PJB1_GEX_S1_R2_001',
                                       'PJB1_TCR_S2_R1_001',
                                       'PJB1_TCR_S2_R2_001',
                                       'PJB2_GEX_S3_R1_001',
                                       'PJB2_GEX_S3_R2_001',
                                       'PJB2_TCR_S4_R1_001',
                                       'PJB2_TCR_S4_R2_001',
                                   ),
                                   include_cellranger_multi=False,
                                   include_cellranger_count=True,
                                   cellranger_pipelines=('cellranger',),
                                   cellranger_samples=(
                                       'PJB1_GEX',
                                       'PJB2_GEX',
                                   ),
                                   cellranger_version="7.1.0")
        # Make 10x_multi_config.csv files (one for each physical
        # sample)
        fastq_dir = os.path.join(os.path.dirname(qc_dir),'fastqs')
        with open(os.path.join(qc_dir,"10x_multi_config.PJB1.csv"),
                  'wt') as fp:
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A

[vdj]
reference,/data/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB1_GEX,{fastq_dir},any,PJB1,gene expression,
PJB1_TCR,{fastq_dir},any,PJB1,VDJ-T,
""".format(fastq_dir=fastq_dir))
        with open(os.path.join(qc_dir,"10x_multi_config.PJB2.csv"),
                  'wt') as fp:
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A

[vdj]
reference,/data/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB2_GEX,{fastq_dir},any,PJB2,gene expression,
PJB2_TCR,{fastq_dir},any,PJB2,VDJ-T,
""".format(fastq_dir=fastq_dir))
        # Collect and checkout outputs
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['cellranger_count',
                          'fastqc_r1',
                          'fastqc_r2',
                          'multiqc',
                          'screens_r2',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_GEX_S1_R1_001',
                          'PJB1_GEX_S1_R2_001',
                          'PJB1_TCR_S2_R1_001',
                          'PJB1_TCR_S2_R2_001',
                          'PJB2_GEX_S3_R1_001',
                          'PJB2_GEX_S3_R2_001',
                          'PJB2_TCR_S4_R1_001',
                          'PJB2_TCR_S4_R2_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1_GEX',
                          'PJB1_TCR',
                          'PJB2_GEX',
                          'PJB2_TCR'])
        self.assertEqual(qc_outputs.seq_data_samples,
                         ['PJB1_GEX',
                          'PJB2_GEX'])
        self.assertEqual(qc_outputs.bams,[])
        self.assertEqual(qc_outputs.organisms,[])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,
                         ['/data/refdata-cellranger-2020-A'])
        self.assertEqual(qc_outputs.cellranger_probe_sets,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.physical_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'cellranger': [ '7.1.0' ],
                           'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r2'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r2'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r2'],76)
        self.assertEqual(qc_outputs.config_files,
                         ['10x_multi_config.PJB1.csv',
                          '10x_multi_config.PJB2.csv',
                          'fastq_strand.conf'])

    def test_qcoutputs_10x_cellranger_single_cell_immune_profiling_800(self):
        """
        QCOutputs: 10xGenomics single cell immune profiling (8.0.0)
        """
        qc_dir = self._make_qc_dir('qc',
                                   protocol="10x_ImmuneProfiling",
                                   fastq_names=(
                                       'PJB1_GEX_S1_R1_001',
                                       'PJB1_GEX_S1_R2_001',
                                       'PJB1_TCR_S2_R1_001',
                                       'PJB1_TCR_S2_R2_001',
                                       'PJB2_GEX_S3_R1_001',
                                       'PJB2_GEX_S3_R2_001',
                                       'PJB2_TCR_S4_R1_001',
                                       'PJB2_TCR_S4_R2_001',
                                   ),
                                   include_cellranger_multi=False,
                                   include_cellranger_count=True,
                                   cellranger_pipelines=('cellranger',),
                                   cellranger_samples=(
                                       'PJB1_GEX',
                                       'PJB2_GEX',
                                   ),
                                   cellranger_version="8.0.0")
        # Make 10x_multi_config.csv files (one for each physical
        # sample)
        fastq_dir = os.path.join(os.path.dirname(qc_dir),'fastqs')
        with open(os.path.join(qc_dir,"10x_multi_config.PJB1.csv"),
                  'wt') as fp:
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A
create-bam,true

[vdj]
reference,/data/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB1_GEX,{fastq_dir},any,PJB1,gene expression,
PJB1_TCR,{fastq_dir},any,PJB1,VDJ-T,
""".format(fastq_dir=fastq_dir))
        with open(os.path.join(qc_dir,"10x_multi_config.PJB2.csv"),
                  'wt') as fp:
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A

[vdj]
reference,/data/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB2_GEX,{fastq_dir},any,PJB2,gene expression,
PJB2_TCR,{fastq_dir},any,PJB2,VDJ-T,
""".format(fastq_dir=fastq_dir))
        # Collect and checkout outputs
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['cellranger_count',
                          'fastqc_r1',
                          'fastqc_r2',
                          'multiqc',
                          'screens_r2',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_GEX_S1_R1_001',
                          'PJB1_GEX_S1_R2_001',
                          'PJB1_TCR_S2_R1_001',
                          'PJB1_TCR_S2_R2_001',
                          'PJB2_GEX_S3_R1_001',
                          'PJB2_GEX_S3_R2_001',
                          'PJB2_TCR_S4_R1_001',
                          'PJB2_TCR_S4_R2_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1_GEX',
                          'PJB1_TCR',
                          'PJB2_GEX',
                          'PJB2_TCR'])
        self.assertEqual(qc_outputs.seq_data_samples,
                         ['PJB1_GEX',
                          'PJB2_GEX'])
        self.assertEqual(qc_outputs.bams,[])
        self.assertEqual(qc_outputs.organisms,[])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,
                         ['/data/refdata-cellranger-2020-A'])
        self.assertEqual(qc_outputs.cellranger_probe_sets,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.physical_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'cellranger': [ '8.0.0' ],
                           'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r2'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r2'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r2'],76)
        self.assertEqual(qc_outputs.config_files,
                         ['10x_multi_config.PJB1.csv',
                          '10x_multi_config.PJB2.csv',
                          'fastq_strand.conf'])

    def test_qcoutputs_10x_cellranger_single_cell_immune_profiling_900(self):
        """
        QCOutputs: 10xGenomics single cell immune profiling (9.0.0)
        """
        qc_dir = self._make_qc_dir('qc',
                                   protocol="10x_ImmuneProfiling",
                                   fastq_names=(
                                       'PJB1_GEX_S1_R1_001',
                                       'PJB1_GEX_S1_R2_001',
                                       'PJB1_TCR_S2_R1_001',
                                       'PJB1_TCR_S2_R2_001',
                                       'PJB2_GEX_S3_R1_001',
                                       'PJB2_GEX_S3_R2_001',
                                       'PJB2_TCR_S4_R1_001',
                                       'PJB2_TCR_S4_R2_001',
                                   ),
                                   include_cellranger_multi=False,
                                   include_cellranger_count=True,
                                   cellranger_pipelines=('cellranger',),
                                   cellranger_samples=(
                                       'PJB1_GEX',
                                       'PJB2_GEX',
                                   ),
                                   cellranger_version="9.0.0")
        # Make 10x_multi_config.csv files (one for each physical
        # sample)
        fastq_dir = os.path.join(os.path.dirname(qc_dir),'fastqs')
        with open(os.path.join(qc_dir,"10x_multi_config.PJB1.csv"),
                  'wt') as fp:
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A
create-bam,true

[vdj]
reference,/data/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB1_GEX,{fastq_dir},any,PJB1,gene expression,
PJB1_TCR,{fastq_dir},any,PJB1,VDJ-T,
""".format(fastq_dir=fastq_dir))
        with open(os.path.join(qc_dir,"10x_multi_config.PJB2.csv"),
                  'wt') as fp:
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A

[vdj]
reference,/data/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB2_GEX,{fastq_dir},any,PJB2,gene expression,
PJB2_TCR,{fastq_dir},any,PJB2,VDJ-T,
""".format(fastq_dir=fastq_dir))
        # Collect and checkout outputs
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['cellranger_count',
                          'fastqc_r1',
                          'fastqc_r2',
                          'multiqc',
                          'screens_r2',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_GEX_S1_R1_001',
                          'PJB1_GEX_S1_R2_001',
                          'PJB1_TCR_S2_R1_001',
                          'PJB1_TCR_S2_R2_001',
                          'PJB2_GEX_S3_R1_001',
                          'PJB2_GEX_S3_R2_001',
                          'PJB2_TCR_S4_R1_001',
                          'PJB2_TCR_S4_R2_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1_GEX',
                          'PJB1_TCR',
                          'PJB2_GEX',
                          'PJB2_TCR'])
        self.assertEqual(qc_outputs.seq_data_samples,
                         ['PJB1_GEX',
                          'PJB2_GEX'])
        self.assertEqual(qc_outputs.bams,[])
        self.assertEqual(qc_outputs.organisms,[])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,
                         ['/data/refdata-cellranger-2020-A'])
        self.assertEqual(qc_outputs.cellranger_probe_sets,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.physical_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'cellranger': [ '9.0.0' ],
                           'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r2'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r2'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r2'],76)
        self.assertEqual(qc_outputs.config_files,
                         ['10x_multi_config.PJB1.csv',
                          '10x_multi_config.PJB2.csv',
                          'fastq_strand.conf'])

    def test_qcoutputs_10x_cellranger_multiple_multi_outputs(self):
        """
        QCOutputs: data with multiple 'cellranger multi' outputs
        """
        qc_dir = self._make_qc_dir('qc',
                                   protocol="10x_CellPlex",
                                   fastq_names=(
                                       'PJB1_GEX_S1_R1_001',
                                       'PJB1_GEX_S1_R2_001',
                                       'PJB1_CML_S2_R1_001',
                                       'PJB1_CML_S2_R2_001',
                                       'PJB2_GEX_S3_R1_001',
                                       'PJB2_GEX_S3_R2_001',
                                       'PJB2_CML_S4_R1_001',
                                       'PJB2_CML_S4_R2_001',
                                   ),
                                   include_cellranger_multi=False,
                                   include_cellranger_count=False,
                                   cellranger_version="9.0.0")
        # Add cellranger multi outputs for multiple physical samples
        for smpl in ("PJB1", "PJB2"):
            cellranger_multi_dir = os.path.join(
                "cellranger_multi",
                "9.0.0",
                "refdata-cellranger-gex-GRCh38-2020-A",
                smpl)
            if smpl == "PJB1":
                multiplexed_samples = ("PJB01", "PJB02")
            elif smpl == "PJB2":
                multiplexed_samples = ("PJB03", "PJB04")
            MockQCOutputs.cellranger_multi(multiplexed_samples,
                                           qc_dir,
                                           prefix=cellranger_multi_dir)
        # Check the outputs
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['cellranger_multi',
                          'fastqc_r1',
                          'fastqc_r2',
                          'multiqc',
                          'screens_r2',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_CML_S2_R1_001',
                          'PJB1_CML_S2_R2_001',
                          'PJB1_GEX_S1_R1_001',
                          'PJB1_GEX_S1_R2_001',
                          'PJB2_CML_S4_R1_001',
                          'PJB2_CML_S4_R2_001',
                          'PJB2_GEX_S3_R1_001',
                          'PJB2_GEX_S3_R2_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1_CML','PJB1_GEX','PJB2_CML','PJB2_GEX'])
        self.assertEqual(qc_outputs.seq_data_samples,
                         ['PJB1_CML','PJB1_GEX','PJB2_CML','PJB2_GEX'])
        self.assertEqual(qc_outputs.bams,[])
        self.assertEqual(qc_outputs.organisms,[])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.cellranger_probe_sets,[])
        self.assertEqual(qc_outputs.multiplexed_samples,
                         ["PJB01", "PJB02", "PJB03", "PJB04"])
        self.assertEqual(qc_outputs.physical_samples,
                         ["PJB1", "PJB2"])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'cellranger': [ '9.0.0' ],
                           'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r2'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r2'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r2'],76)
        self.assertEqual(qc_outputs.config_files,
                         ['fastq_strand.conf'])

    def test_qcoutputs_10x_cellranger_external_multi_outputs(self):
        """
        QCOutputs: data with external 'cellranger multi' outputs
        """
        qc_dir = self._make_qc_dir('qc',
                                   protocol="10x_CellPlex",
                                   fastq_names=(
                                       'PJB1_GEX_S1_R1_001',
                                       'PJB1_GEX_S1_R2_001',
                                       'PJB2_MC_S2_R1_001',
                                       'PJB2_MC_S2_R2_001',
                                   ),
                                   include_cellranger_multi=False,
                                   cellranger_pipelines=('cellranger',),
                                   cellranger_samples=(
                                       'PJB1_GEX',
                                       'PJB2_MC',
                                   ),
                                   cellranger_multi_samples=(
                                       'PJB_CML1',
                                       'PJB_CML2',
                                   ))
        external_multi_dir = os.path.join("cellranger_multi",
                                          "7.1.0",
                                          "external")
        MockQCOutputs.cellranger_multi(("EX1",),
                                       qc_dir,
                                       prefix=external_multi_dir)
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['cellranger_multi',
                          'fastqc_r1',
                          'fastqc_r2',
                          'multiqc',
                          'screens_r2',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_GEX_S1_R1_001',
                          'PJB1_GEX_S1_R2_001',
                          'PJB2_MC_S2_R1_001',
                          'PJB2_MC_S2_R2_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1_GEX','PJB2_MC'])
        self.assertEqual(qc_outputs.seq_data_samples,
                         ['PJB1_GEX','PJB2_MC'])
        self.assertEqual(qc_outputs.bams,[])
        self.assertEqual(qc_outputs.organisms,[])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.cellranger_probe_sets,[])
        self.assertEqual(qc_outputs.multiplexed_samples,
                         ['EX1'])
        self.assertEqual(qc_outputs.physical_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'cellranger': [ '7.1.0' ],
                           'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r2'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r2'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r2'],76)
        self.assertEqual(qc_outputs.config_files,
                         ['fastq_strand.conf'])

    def test_qcoutputs_with_non_standard_qc_dir(self):
        """
        QCOutputs: non-standard QC directory
        """
        qc_dir = self._make_qc_dir('my_qc_outs',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB1_S1_R2_001',
                                       'PJB2_S2_R1_001',
                                       'PJB2_S2_R2_001',
                                   ))
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['fastqc_r1',
                          'fastqc_r2',
                          'multiqc',
                          'screens_r1',
                          'screens_r2',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB1_S1_R2_001',
                          'PJB2_S2_R1_001',
                          'PJB2_S2_R2_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.seq_data_samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.bams,[])
        self.assertEqual(qc_outputs.organisms,[])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.physical_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r2'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r2'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r2'],76)
        self.assertEqual(qc_outputs.config_files,
                         ['fastq_strand.conf'])

    def test_qcoutputs_with_non_canonical_fastq_names(self):
        """
        QCOutputs: non-canonical fastq names
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001_paired',
                                       'PJB1_S1_R2_001_paired',
                                       'PJB2_S2_R1_001_paired',
                                       'PJB2_S2_R2_001_paired',
                                   ))
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['fastqc_r1',
                          'fastqc_r2',
                          'multiqc',
                          'screens_r1',
                          'screens_r2',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001_paired',
                          'PJB1_S1_R2_001_paired',
                          'PJB2_S2_R1_001_paired',
                          'PJB2_S2_R2_001_paired'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.seq_data_samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.bams,[])
        self.assertEqual(qc_outputs.organisms,[])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.physical_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r2'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r2'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r2'],76)
        self.assertEqual(qc_outputs.config_files,
                         ['fastq_strand.conf'])

class TestExtraOutputs(unittest.TestCase):
    def setUp(self):
        # Temporary working dir
        self.wd = tempfile.mkdtemp(suffix='.test_ExtraOutputs')
    def tearDown(self):
        # Remove temporary working dir
        if not REMOVE_TEST_OUTPUTS:
            return
        if self.wd is not None and os.path.isdir(self.wd):
            shutil.rmtree(self.wd)

    def test_extra_outputs(self):
        """
        ExtraOutputs: load data from TSV file
        """
        tsv_file = os.path.join(self.wd,"example.tsv")
        with open(tsv_file,'wt') as fp:
            fp.write("""# Example TSV file with external files
stuff.html\tSome random output
#more_stuff.html\tMore random output
external_stuff/index.html\tBunch of external stuff\texternal_stuff/results,external_stuff/more_results

""")
        extra_outputs = ExtraOutputs(tsv_file)
        self.assertEqual(len(extra_outputs.outputs),2)
        self.assertEqual(extra_outputs.outputs[0].file_path,
                         "stuff.html")
        self.assertEqual(extra_outputs.outputs[0].description,
                         "Some random output")
        self.assertEqual(extra_outputs.outputs[0].additional_files,
                         None)
        self.assertEqual(extra_outputs.outputs[1].file_path,
                         "external_stuff/index.html")
        self.assertEqual(extra_outputs.outputs[1].description,
                         "Bunch of external stuff")
        self.assertEqual(extra_outputs.outputs[1].additional_files,
                         ["external_stuff/results",
                          "external_stuff/more_results"])
