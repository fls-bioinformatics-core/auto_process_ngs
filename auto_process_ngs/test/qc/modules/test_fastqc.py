#######################################################################
# Unit tests for qc/modules/fastqc.py ('fastqc' QC module)
#######################################################################

import unittest
import tempfile
import shutil
import os
from auto_process_ngs.mock import MockAnalysisProject
from auto_process_ngs.mock import UpdateAnalysisProject
from auto_process_ngs.analysis import AnalysisProject
from auto_process_ngs.qc.modules.fastqc import check_fastqc_outputs

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

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

    def test_check_fastqc_outputs_paired_end_all_missing(self):
        """
        check_fastqc_outputs: all FastQC outputs missing (paired end)
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
                                              read_numbers=(1,2)),
                         project.fastqs)

    def test_check_fastqc_outputs_paired_end_all_present(self):
        """
        check_fastqc_outputs: all FastQC outputs present (paired end)
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
                                              read_numbers=(1,2)),
                         [])

    def test_check_fastqc_outputs_paired_end_some_missing(self):
        """
        check_fastqc_outputs: some FastQC outputs missing (paired end)
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
                                              read_numbers=(1,2)),
                         [os.path.join(project.fastq_dir,
                                       "PJB1_S1_R1_001.fastq.gz")])

    def test_check_fastqc_outputs_single_end_all_missing(self):
        """
        check_fastqc_outputs: all FastQC outputs missing (single end)
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
                                              read_numbers=(1,)),
                         project.fastqs)

    def test_check_fastqc_outputs_single_end_all_present(self):
        """
        check_fastqc_outputs: all FastQC outputs present (single end)
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
                                              read_numbers=(1,)),
                         [])

    def test_check_fastqc_outputs_single_end_some_missing(self):
        """
        check_fastqc_outputs: some FastQC outputs missing (single end)
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
                                              read_numbers=(1,)),
                         [os.path.join(project.fastq_dir,
                                       "PJB1_S1_R1_001.fastq.gz")])

    def test_check_fastqc_outputs_single_cell_all_missing(self):
        """
        check_fastqc_outputs: all FastQC outputs missing (single cell)
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
                                              read_numbers=(1,2,)),
                         project.fastqs)

    def test_check_fastqc_outputs_single_cell_all_present(self):
        """
        check_fastqc_outputs: all FastQC outputs present (single cell)
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
                                              read_numbers=(1,2,)),
                         [])

    def test_check_fastqc_outputs_single_cell_some_missing(self):
        """
        check_fastqc_outputs: some FastQC outputs missing (single cell)
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
                                              read_numbers=(1,2,)),
                         [os.path.join(project.fastq_dir,
                                       "PJB1_S1_R1_001.fastq.gz")])

    def test_check_fastqc_outputs_from_fastq_list_all_missing(self):
        """
        check_fastqc_outputs: all FastQC outputs missing (Fastq list)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # List of Fastqs (subset of all Fastqs)
        fastqs = ["PJB2_S2_R1_001.fastq.gz",
                  "PJB2_S2_R2_001.fastq.gz"]
        # Get the outputs
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        # Check
        self.assertEqual(check_fastqc_outputs(project,
                                              qc_dir="qc",
                                              fastqs=fastqs,
                                              read_numbers=(1,2,)),
                         fastqs)

    def test_check_fastqc_outputs_from_fastq_list_all_present(self):
        """
        check_fastqc_outputs: all FastQC outputs present (Fastq list)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Add QC artefacts
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            include_fastq_strand=False,
            include_multiqc=False)
        # List of Fastqs (subset of all Fastqs)
        fastqs = ["PJB2_S2_R1_001.fastq.gz",
                  "PJB2_S2_R2_001.fastq.gz"]
        # Remove FastQC outputs for Fastqs not in the list
        for f in ("PJB1_S1_R1_001_fastqc.html",
                  "PJB1_S1_R2_001_fastqc.html"):
            os.remove(os.path.join(project.qc_dir,f))
        # Check
        self.assertEqual(check_fastqc_outputs(project,
                                              qc_dir="qc",
                                              fastqs=fastqs,
                                              read_numbers=(1,2,)),
                         [])

    def test_check_fastqc_outputs_from_fastq_list_some_missing(self):
        """
        check_fastqc_outputs: some FastQC outputs missing (Fastq list)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Add QC artefacts
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            include_fastq_strand=False,
            include_multiqc=False)
        # List of Fastqs (subset of all Fastqs)
        fastqs = ["PJB2_S2_R1_001.fastq.gz",
                  "PJB2_S2_R2_001.fastq.gz"]
        # Remove FastQC outputs for Fastqs not in the list
        for f in ("PJB1_S1_R1_001_fastqc.html",
                  "PJB1_S1_R2_001_fastqc.html"):
            os.remove(os.path.join(project.qc_dir,f))
        # Remove a FastQC output from the list
        os.remove(os.path.join(project.qc_dir,
                               "PJB2_S2_R1_001_fastqc.html"))
        # Check
        self.assertEqual(check_fastqc_outputs(project,
                                              qc_dir="qc",
                                              fastqs=fastqs,
                                              read_numbers=(1,2,)),
                         ["PJB2_S2_R1_001.fastq.gz"])
