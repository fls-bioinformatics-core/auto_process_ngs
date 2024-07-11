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
from auto_process_ngs.qc.modules.fastq_screen import check_fastq_screen_outputs

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

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

    def test_check_fastq_screen_outputs_paired_end_all_missing(self):
        """
        check_fastq_screen_outputs: all screen outputs missing (paired end)
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
                                                    read_numbers=(1,2)),
                         project.fastqs)

    def test_check_fastq_screen_outputs_paired_end_all_present(self):
        """
        check_fastq_screen_outputs: all screen outputs present (paired end)
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
                                                    read_numbers=(1,2)),
                         [])

    def test_check_fastq_screen_outputs_paired_end_some_missing(self):
        """
        check_fastq_screen_outputs: some screen outputs missing (paired end)
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
                                                    read_numbers=(1,2)),
                         [os.path.join(project.fastq_dir,
                                       "PJB1_S1_R1_001.fastq.gz")])

    def test_check_fastq_screen_outputs_paired_end_legacy(self):
        """
        check_fastq_screen_outputs: all screen outputs present (paired, legacy)
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
                                                    read_numbers=(1,2),
                                                    legacy=True),
                         [])

    def test_check_fastq_screen_outputs_single_end_all_missing(self):
        """
        check_fastq_screen_outputs: all screen outputs missing (single end)
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
                                                    read_numbers=(1,2)),
                         project.fastqs)

    def test_check_fastq_screen_outputs_single_end_all_present(self):
        """
        check_fastq_screen_outputs: all screen outputs present (single end)
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
                                                    read_numbers=(1,)),
                         [])

    def test_check_fastq_screen_outputs_single_end_some_missing(self):
        """
        check_fastq_screen_outputs: some screen outputs missing (single end)
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
                                                    read_numbers=(1,)),
                         [os.path.join(project.fastq_dir,
                                       "PJB1_S1_R1_001.fastq.gz")])

    def test_check_fastq_screen_outputs_single_cell_all_missing(self):
        """
        check_fastq_screen_outputs: all screen outputs missing (single cell)
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
                                                    read_numbers=(2,)),
                         [os.path.join(project.fastq_dir,
                                       "PJB1_S1_R2_001.fastq.gz")])

    def test_check_fastq_screen_outputs_single_cell_all_present(self):
        """
        check_fastq_screen_outputs: all screen outputs present (single cell)
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
                                                    read_numbers=(2,)),
                         [])

    def test_check_fastq_screen_outputs_single_cell_some_missing(self):
        """
        check_fastq_screen_outputs: some screen outputs missing (single cell)
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
                                                    read_numbers=(2,)),
                         [os.path.join(project.fastq_dir,
                                       "PJB1_S1_R2_001.fastq.gz")])

    def test_check_fastq_screen_outputs_parseevercode_all_missing(self):
        """
        check_fastq_screen_outputs: all screen outputs missing (ParseEvercode)
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
                                                    read_numbers=(1,)),
                         [os.path.join(project.fastq_dir,
                                       "PJB1_S1_R1_001.fastq.gz")])

    def test_check_fastq_screen_outputs_parseevercode_all_present(self):
        """
        check_fastq_screen_outputs: all screen outputs present (ParseEvercode)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Add QC artefacts
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            protocol="ParseEvercode",
            include_fastq_strand=False,
            include_multiqc=False)
        # Check
        self.assertEqual(check_fastq_screen_outputs(project,
                                                    qc_dir="qc",
                                                    screen="model_organisms",
                                                    read_numbers=(1,)),
                         [])

    def test_check_fastq_screen_outputs_parseevercode_some_missing(self):
        """
        check_fastq_screen_outputs: some screen outputs missing (ParseEvercode)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Add QC artefacts
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            protocol="ParseEvercode",
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
                                                    read_numbers=(1,)),
                         [os.path.join(project.fastq_dir,
                                       "PJB1_S1_R1_001.fastq.gz")])

    def test_check_fastq_screen_outputs_all_missing_specify_fastqs(self):
        """
        check_fastq_screen_outputs: all screen outputs missing (specify Fastqs)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Copy Fastqs to alternative location
        alt_fastqs_dir = os.path.join(self.wd,"alt_fastqs")
        os.mkdir(alt_fastqs_dir)
        project = AnalysisProject(os.path.join(self.wd,"PJB"))
        for fq in project.fastqs:
            shutil.copyfile(fq,os.path.join(alt_fastqs_dir,
                                            os.path.basename(fq)))
        alt_fastqs = [os.path.join(alt_fastqs_dir,os.path.basename(fq))
                      for fq in project.fastqs]
        # Check
        self.assertEqual(check_fastq_screen_outputs(project,
                                                    qc_dir="qc",
                                                    screen="model_organisms",
                                                    read_numbers=(1,2),
                                                    fastqs=alt_fastqs),
                         alt_fastqs)

    def test_check_fastq_screen_outputs_all_present_specify_fastqs(self):
        """
        check_fastq_screen_outputs: all screen outputs present (specify Fastqs)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Copy Fastqs to alternative location
        alt_fastqs_dir = os.path.join(self.wd,"alt_fastqs")
        os.mkdir(alt_fastqs_dir)
        project = AnalysisProject(os.path.join(self.wd,"PJB"))
        for fq in project.fastqs:
            shutil.copyfile(fq,os.path.join(alt_fastqs_dir,
                                            os.path.basename(fq)))
        alt_fastqs = [os.path.join(alt_fastqs_dir,os.path.basename(fq))
                      for fq in project.fastqs]
        # Add QC artefacts
        UpdateAnalysisProject(project).add_qc_outputs(
            include_fastq_strand=False,
            include_multiqc=False)
        # Check
        self.assertEqual(check_fastq_screen_outputs(project,
                                                    qc_dir="qc",
                                                    screen="model_organisms",
                                                    read_numbers=(1,2),
                                                    fastqs=alt_fastqs),
                         [])
