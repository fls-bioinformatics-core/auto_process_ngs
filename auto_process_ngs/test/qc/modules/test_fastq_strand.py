#######################################################################
# Unit tests for qc/modules/fastq_strand.py ('fastq_strand' QC module)
#######################################################################

import unittest
import tempfile
import shutil
import os
from auto_process_ngs.mock import MockAnalysisProject
from auto_process_ngs.mock import UpdateAnalysisProject
from auto_process_ngs.analysis import AnalysisProject
from auto_process_ngs.qc.modules.fastq_strand import check_fastq_strand_outputs

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

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

    def test_check_fastq_strand_outputs_paired_end_missing(self):
        """
        check_fastq_strand_outputs: fastq_strand.py output missing (paired end)
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
                                                    read_numbers=(1,2)),
                         [(os.path.join(project.fastq_dir,
                                        "PJB1_S1_R1_001.fastq.gz"),
                           os.path.join(project.fastq_dir,
                                        "PJB1_S1_R2_001.fastq.gz")),])

    def test_check_fastq_strand_outputs_paired_end_present(self):
        """
        check_fastq_strand_outputs: fastq_strand.py output present (paired end)
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
                                                    read_numbers=(1,2)),
                         [])

    def test_check_fastq_strand_outputs_single_end_missing(self):
        """
        check_fastq_strand_outputs: fastq_strand.py output missing (single end)
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
                                                    read_numbers=(1,)),
                         [(os.path.join(project.fastq_dir,
                                        "PJB1_S1_R1_001.fastq.gz"),),])

    def test_check_fastq_strand_outputs_single_end_present(self):
        """
        check_fastq_strand_outputs: fastq_strand.py output present (single end)
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
                                                    read_numbers=(1,)),
                         [])

    def test_check_fastq_strand_outputs_single_cell_missing(self):
        """
        check_fastq_strand_outputs: fastq_strand.py output missing (single cell)
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
                                                    read_numbers=(2,)),
                         [(os.path.join(project.fastq_dir,
                                        "PJB1_S1_R2_001.fastq.gz"),),])

    def test_check_fastq_strand_outputs_single_cell_present(self):
        """
        check_fastq_strand_outputs: fastq_strand.py output present (single cell)
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
                                                    read_numbers=(2,)),
                         [])

    def test_check_fastq_strand_outputs_parseevercode_missing(self):
        """
        check_fastq_strand_outputs: fastq_strand.py output missing (ParseEvercode)
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
                                                    read_numbers=(1,)),
                         [(os.path.join(project.fastq_dir,
                                        "PJB1_S1_R1_001.fastq.gz"),),])

    def test_check_fastq_strand_outputs_parseevercode_present(self):
        """
        check_fastq_strand_outputs: fastq_strand.py output present (ParseEvercode)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            protocol="ParseEvercode",
            include_fastq_strand=True,
            include_multiqc=False)
        fastq_strand_conf = os.path.join(project.dirn,"fastq_strand.conf")
        # Check the outputs
        self.assertEqual(check_fastq_strand_outputs(project,
                                                    "qc",
                                                    fastq_strand_conf,
                                                    read_numbers=(1,)),
                         [])

    def test_check_fastq_strand_outputs_all_missing_specify_fastqs(self):
        """
        check_fastq_strand_outputs: fastq_strand.py output missing (specify Fastqs)
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
        # Copy Fastqs to alternative location
        alt_fastqs_dir = os.path.join(self.wd,"alt_fastqs")
        os.mkdir(alt_fastqs_dir)
        for fq in project.fastqs:
            shutil.copyfile(fq,os.path.join(alt_fastqs_dir,
                                            os.path.basename(fq)))
        alt_fastqs = [os.path.join(alt_fastqs_dir,os.path.basename(fq))
                      for fq in project.fastqs]
        # Check the outputs
        self.assertEqual(check_fastq_strand_outputs(project,
                                                    "qc",
                                                    fastq_strand_conf,
                                                    read_numbers=(1,2),
                                                    fastqs=alt_fastqs),
                         [(os.path.join(alt_fastqs_dir,
                                        "PJB1_S1_R1_001.fastq.gz"),
                           os.path.join(alt_fastqs_dir,
                                        "PJB1_S1_R2_001.fastq.gz")),])

    def test_check_fastq_strand_outputs_all_present_specify_fastqs(self):
        """
        check_fastq_strand_outputs: fastq_strand.py output present (specify Fastqs)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        project = AnalysisProject(os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            protocol="standardPE",
            include_fastq_strand=True,
            include_multiqc=False)
        fastq_strand_conf = os.path.join(project.dirn,"fastq_strand.conf")
        with open(fastq_strand_conf,'w') as fp:
            fp.write("")
        # Copy Fastqs to alternative location
        alt_fastqs_dir = os.path.join(self.wd,"alt_fastqs")
        os.mkdir(alt_fastqs_dir)
        for fq in project.fastqs:
            shutil.copyfile(fq,os.path.join(alt_fastqs_dir,
                                            os.path.basename(fq)))
        alt_fastqs = [os.path.join(alt_fastqs_dir,os.path.basename(fq))
                      for fq in project.fastqs]
        # Check the outputs
        self.assertEqual(check_fastq_strand_outputs(project,
                                                    "qc",
                                                    fastq_strand_conf,
                                                    read_numbers=(1,2),
                                                    fastqs=alt_fastqs),
                         [])
