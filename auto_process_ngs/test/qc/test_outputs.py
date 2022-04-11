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
from auto_process_ngs.analysis import AnalysisProject
from auto_process_ngs.qc.outputs import QCOutputs
from auto_process_ngs.qc.outputs import fastq_screen_output
from auto_process_ngs.qc.outputs import fastqc_output
from auto_process_ngs.qc.outputs import fastq_strand_output
from auto_process_ngs.qc.outputs import cellranger_count_output
from auto_process_ngs.qc.outputs import cellranger_atac_count_output
from auto_process_ngs.qc.outputs import cellranger_arc_count_output
from auto_process_ngs.qc.outputs import cellranger_multi_output
from auto_process_ngs.qc.outputs import check_fastq_screen_outputs
from auto_process_ngs.qc.outputs import check_fastqc_outputs
from auto_process_ngs.qc.outputs import check_fastq_strand_outputs
from auto_process_ngs.qc.outputs import check_cellranger_count_outputs
from auto_process_ngs.qc.outputs import check_cellranger_atac_count_outputs
from auto_process_ngs.qc.outputs import check_cellranger_arc_count_outputs

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
        if self.wd is None:
            self.wd = tempfile.mkdtemp(suffix='.test_QCOutputs')
        return make_mock_qc_dir(
            os.path.join(self.wd,qc_dir),
            fastq_names,
            screens=screens,
            cellranger_pipelines=cellranger_pipelines,
            cellranger_samples=cellranger_samples,
            cellranger_multi_samples=cellranger_multi_samples,
            include_fastqc=include_fastqc,
            include_fastq_screen=include_fastq_screen,
            include_strandedness=include_strandedness,
            include_seqlens=include_seqlens,
            include_multiqc=include_multiqc,
            include_cellranger_count=include_cellranger_count,
            include_cellranger_multi=include_cellranger_multi,
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
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
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
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
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
        self.assertEqual(qc_outputs.fastq_screens,[])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
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
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
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
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
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
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
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
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
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
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
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
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
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
        self.assertEqual(qc_outputs.fastq_screens,[])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
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
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
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
        self.assertEqual(qc_outputs.multiplexed_samples,[])
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
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
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
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
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

    def test_qcoutputs_10x_cellranger_count(self):
        """
        QCOutputs: 10xGenomics data with cellranger 'count'
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
                                   ))
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
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,
                         ['/data/refdata-cellranger-2020-A'])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
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
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,
                         ['/data/refdata-cellranger-2020-A'])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
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
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,
                         ['/data/refdata-cellranger-atac-2020-A'])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
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
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,
                         ['/data/refdata-cellranger-2020-A',
                          '/data/refdata-cellranger-arc-2020-A'])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'cellranger': [ '6.1.2' ],
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
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,
                         ['/data/refdata-cellranger-arc-2020-A',
                          '/data/refdata-cellranger-atac-2020-A'])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
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

    def test_qcoutputs_10x_cellranger_multi(self):
        """
        QCOutputs: 10xGenomics CellPlex data with 'cellranger multi'
        """
        qc_dir = self._make_qc_dir('qc',
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
                                   ))
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['cellranger_multi',
                          'fastqc_r1',
                          'fastqc_r2',
                          'multiqc',
                          'screens_r1',
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
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,
                         ['/data/refdata-cellranger-2020-A'])
        self.assertEqual(qc_outputs.multiplexed_samples,
                         ['PJB_CML1','PJB_CML2'])
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

    def test_qcoutputs_10x_cellranger_multi_and_count(self):
        """
        QCOutputs: 10xGenomics CellPlex data with 'cellranger multi' and 'count'
        """
        qc_dir = self._make_qc_dir('qc',
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
                                   ))
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['cellranger_count',
                          'cellranger_multi',
                          'fastqc_r1',
                          'fastqc_r2',
                          'multiqc',
                          'screens_r1',
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
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,
                         ['/data/refdata-cellranger-2020-A'])
        self.assertEqual(qc_outputs.multiplexed_samples,
                         ['PJB_CML1','PJB_CML2'])
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
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
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
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
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

class TestFastqScreenOutputFunction(unittest.TestCase):
    def test_fastq_screen_output(self):
        """fastq_screen_output: handles .fastq file
        """
        self.assertEqual(fastq_screen_output('/data/PB/PB1_ATTAGG_L001_R1_001.fastq',
                                             'model_organisms'),
                         ('PB1_ATTAGG_L001_R1_001_screen_model_organisms.png',
                          'PB1_ATTAGG_L001_R1_001_screen_model_organisms.txt'))
    def test_fastq_screen_output_fastqgz(self):
        """fastq_screen_output: handles fastq.gz file
        """
        self.assertEqual(fastq_screen_output('/data/PB/PB1_ATTAGG_L001_R1_001.fastq.gz',
                                             'model_organisms'),
                         ('PB1_ATTAGG_L001_R1_001_screen_model_organisms.png',
                          'PB1_ATTAGG_L001_R1_001_screen_model_organisms.txt'))
    def test_fastq_screen_output_legacy_naming(self):
        """fastq_screen_output: handles legacy naming convention
        """
        self.assertEqual(fastq_screen_output('/data/PB/PB1_ATTAGG_L001_R1_001.fastq.gz',
                                             'model_organisms',
                                             legacy=True),
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
            "PJB1_S1_R1_001_screen_model_organisms.txt"))
        # Check
        self.assertEqual(check_fastq_screen_outputs(project,
                                                    qc_dir="qc",
                                                    screen="model_organisms",
                                                    qc_protocol="standardPE"),
                         [os.path.join(project.fastq_dir,
                                       "PJB1_S1_R1_001.fastq.gz")])

    def test_check_fastq_screen_outputs_standardPE_legacy(self):
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
            include_multiqc=False,
            legacy_screens=True)
        # Check
        self.assertEqual(check_fastq_screen_outputs(project,
                                                    qc_dir="qc",
                                                    screen="model_organisms",
                                                    qc_protocol="standardPE",
                                                    legacy=True),
                         [])

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
            "PJB1_S1_R1_001_screen_model_organisms.txt"))
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
            "PJB1_S1_R2_001_screen_model_organisms.txt"))
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
