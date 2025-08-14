#######################################################################
# Tests for cli/transfer_data.py utility
#######################################################################

import unittest
import tempfile
import shutil
import os
from auto_process_ngs.cli.transfer_data import main as transfer_data
from auto_process_ngs.analysis import AnalysisProject
from auto_process_ngs.mock import MockAnalysisDirFactory
from auto_process_ngs.mock import UpdateAnalysisProject
from auto_process_ngs.settings import Settings

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

# Unit tests

class TestTransferData(unittest.TestCase):
    """
    Tests for the 'transfer_data' utility
    """
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestTransferData')
        # Create settings instance
        # This allows us to set the polling interval for the
        # unit tests
        settings_ini = os.path.join(self.dirn, "auto_process.ini")
        with open(settings_ini,'w') as s:
            s.write("""[general]
poll_interval = 0.5
""")
        self.settings = Settings(settings_ini)
        # Store original location so we can get back at the end
        self.pwd = os.getcwd()
        # Move to working dir
        os.chdir(self.dirn)

    def tearDown(self):
        # Return to original dir
        os.chdir(self.pwd)
        # Remove the temporary test directory
        def del_rw(action,name,excinfo):
            # Explicitly remove read only files/
            # dirs
            if os.path.isfile(name):
                os.chmod(os.path.dirname(name),0o755)
                os.chmod(name,0o655)
                os.remove(name)
            elif os.path.isdir(name):
                os.chmod(os.path.dirname(name),0o755)
                os.chmod(name,0o755)
                os.rmdir(name)
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.dirn,onerror=del_rw)

    def test_transfer_data_all_fastqs(self):
        """
        transfer_data: copy all Fastqs
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901",
                       "run_number": "89" },
            project_metadata={ "AB": { "Library type": "RNA-seq",
                                       "Organism": "Human" } },
            top_dir=self.dirn)
        mockdir.create()
        # Make a target directory
        target_dir = os.path.join(self.dirn, "shared")
        os.makedirs(target_dir)
        # Do data transfer
        self.assertEqual(transfer_data([target_dir,
                                        os.path.join(mockdir.dirn, "AB")]), 0)
        # Check transferred artefacts
        print(os.listdir(target_dir))
        expected_files = ("AB1_S1_R1_001.fastq.gz",
                          "AB1_S1_R2_001.fastq.gz",
                          "AB2_S2_R2_001.fastq.gz",
                          "AB2_S2_R1_001.fastq.gz",
                          "AB.chksums")
        for f in expected_files:
            self.assertTrue(os.path.exists(os.path.join(target_dir, f)),
                            f"'{f}': missing, should be present")
        for f in os.listdir(target_dir):
            self.assertTrue(f in expected_files,
                            f"'{f}': present, but not expected")

    def test_transfer_data_R1_fastqs(self):
        """
        transfer_data: copy R1 Fastqs only
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901",
                       "run_number": "89" },
            project_metadata={ "AB": { "Library type": "RNA-seq",
                                       "Organism": "Human" } },
            top_dir=self.dirn)
        mockdir.create()
        # Make a target directory
        target_dir = os.path.join(self.dirn, "shared")
        os.makedirs(target_dir)
        # Do data transfer (--filter *_R1_*)
        self.assertEqual(transfer_data(
            [target_dir,
             os.path.join(mockdir.dirn, "AB"),
             "--filter", "*_R1_*"]), 0)
        # Check transferred artefacts
        print(os.listdir(target_dir))
        expected_files = ("AB1_S1_R1_001.fastq.gz",
                          "AB2_S2_R1_001.fastq.gz",
                          "AB.chksums")
        for f in expected_files:
            self.assertTrue(os.path.exists(os.path.join(target_dir, f)),
                            f"'{f}': missing, should be present")
        for f in os.listdir(target_dir):
            self.assertTrue(f in expected_files,
                            f"'{f}': present, but not expected")

    def test_transfer_data_fastqs_for_single_sample(self):
        """
        transfer_data: copy Fastqs for single sample
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901",
                       "run_number": "89" },
            project_metadata={ "AB": { "Library type": "RNA-seq",
                                       "Organism": "Human" } },
            top_dir=self.dirn)
        mockdir.create()
        # Make a target directory
        target_dir = os.path.join(self.dirn, "shared")
        os.makedirs(target_dir)
        # Do data transfer (--samples AB2)
        self.assertEqual(transfer_data(
            [target_dir,
             os.path.join(mockdir.dirn, "AB"),
             "--samples", "AB2"]), 0)
        # Check transferred artefacts
        print(os.listdir(target_dir))
        expected_files = ("AB2_S2_R2_001.fastq.gz",
                          "AB2_S2_R1_001.fastq.gz",
                          "AB.chksums")
        for f in expected_files:
            self.assertTrue(os.path.exists(os.path.join(target_dir, f)),
                            f"'{f}': missing, should be present")
        for f in os.listdir(target_dir):
            self.assertTrue(f in expected_files,
                            f"'{f}': present, but not expected")

    def test_transfer_data_zip_fastqs(self):
        """
        transfer_data: put Fastqs into ZIP archive
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901",
                       "run_number": "89" },
            project_metadata={ "AB": { "Library type": "RNA-seq",
                                       "Organism": "Human" } },
            top_dir=self.dirn)
        mockdir.create()
        # Make a target directory
        target_dir = os.path.join(self.dirn, "shared")
        os.makedirs(target_dir)
        # Do data transfer (--zip_fastqs)
        self.assertEqual(transfer_data(
            [target_dir,
             os.path.join(mockdir.dirn, "AB"),
             "--zip_fastqs"]), 0)
        # Check transferred artefacts
        print(os.listdir(target_dir))
        expected_files = ("MISEQ_170901.89-AB-fastqs.zip",
                          "MISEQ_170901.89-AB-fastqs.checksums")
        for f in expected_files:
            self.assertTrue(os.path.exists(os.path.join(target_dir, f)),
                            f"'{f}': missing, should be present")
        for f in os.listdir(target_dir):
            self.assertTrue(f in expected_files,
                            f"'{f}': present, but not expected")

    def test_transfer_data_fastqs_with_qc_report(self):
        """
        transfer_data: copy Fastqs and QC report
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901",
                       "run_number": "89" },
            project_metadata={ "AB": { "Library type": "RNA-seq",
                                       "Organism": "Human" } },
            top_dir=self.dirn)
        mockdir.create()
        # Add QC outputs
        project = AnalysisProject(os.path.join(mockdir.dirn, "AB"))
        UpdateAnalysisProject(project).add_qc_outputs()
        # Make a target directory
        target_dir = os.path.join(self.dirn, "shared")
        os.makedirs(target_dir)
        # Do data transfer (--include_qc_report)
        self.assertEqual(transfer_data(
            [target_dir,
             os.path.join(mockdir.dirn, "AB"),
             "--include_qc_report"]), 0)
        # Check transferred artefacts
        print(os.listdir(target_dir))
        expected_files = ("AB1_S1_R1_001.fastq.gz",
                          "AB1_S1_R2_001.fastq.gz",
                          "AB2_S2_R2_001.fastq.gz",
                          "AB2_S2_R1_001.fastq.gz",
                          "AB.chksums",
                          "qc_report.AB.170901_M00879_0087_000000000-AGEW9.zip")
        for f in expected_files:
            self.assertTrue(os.path.exists(os.path.join(target_dir, f)),
                            f"'{f}': missing, should be present")
        for f in os.listdir(target_dir):
            self.assertTrue(f in expected_files,
                            f"'{f}': present, but not expected")

    def test_transfer_data_no_fastqs_with_qc_report(self):
        """
        transfer_data: copy QC report (no Fastqs)
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901",
                       "run_number": "89" },
            project_metadata={ "AB": { "Library type": "RNA-seq",
                                       "Organism": "Human" } },
            top_dir=self.dirn)
        mockdir.create()
        # Add QC outputs
        project = AnalysisProject(os.path.join(mockdir.dirn, "AB"))
        UpdateAnalysisProject(project).add_qc_outputs()
        # Make a target directory
        target_dir = os.path.join(self.dirn, "shared")
        os.makedirs(target_dir)
        # Do data transfer (--include_qc_report --no_fastqs)
        self.assertEqual(transfer_data(
            [target_dir,
             os.path.join(mockdir.dirn, "AB"),
             "--include_qc_report",
             "--no_fastqs"]), 0)
        # Check transferred artefacts
        print(os.listdir(target_dir))
        expected_files = ("qc_report.AB.170901_M00879_0087_000000000-AGEW9.zip",)
        for f in expected_files:
            self.assertTrue(os.path.exists(os.path.join(target_dir, f)),
                            f"'{f}': missing, should be present")
        for f in os.listdir(target_dir):
            self.assertTrue(f in expected_files,
                            f"'{f}': present, but not expected")

    def test_transfer_data_fastqs_with_10x_cellranger_count_outputs(self):
        """
        transfer_data: copy Fastqs and 10xGenomics 'cellranger count' outputs
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901",
                       "run_number": "89" },
            project_metadata={ "AB": { "Single cell platform":
                                       "10xGenomics Chromium 3'v3",
                                       "Library type": "RNA-seq",
                                       "Organism": "Human" } },
            top_dir=self.dirn)
        mockdir.create()
        # Add 10x Genomics pipeline outputs
        project = AnalysisProject(os.path.join(mockdir.dirn, "AB"))
        UpdateAnalysisProject(project).add_qc_outputs()
        UpdateAnalysisProject(project).add_cellranger_count_outputs()
        # Make a target directory
        target_dir = os.path.join(self.dirn, "shared")
        os.makedirs(target_dir)
        # Do data transfer (--include_10x_outputs)
        self.assertEqual(transfer_data(
            [target_dir,
             os.path.join(mockdir.dirn, "AB"),
             "--include_10x_outputs"]), 0)
        # Check transferred artefacts
        print(os.listdir(target_dir))
        expected_files = ("AB1_S1_R1_001.fastq.gz",
                          "AB1_S1_R2_001.fastq.gz",
                          "AB2_S2_R2_001.fastq.gz",
                          "AB2_S2_R1_001.fastq.gz",
                          "AB.chksums",
                          "cellranger_count.AB.170901_M00879_0087_000000000-AGEW9.tgz")
        for f in expected_files:
            self.assertTrue(os.path.exists(os.path.join(target_dir, f)),
                            f"'{f}': missing, should be present")
        for f in os.listdir(target_dir):
            self.assertTrue(f in expected_files,
                            f"'{f}': present, but not expected")

    def test_transfer_data_fastqs_with_10x_cellranger_multi_outputs(self):
        """
        transfer_data: copy Fastqs and 10xGenomics 'cellranger multi' outputs
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901",
                       "run_number": "89" },
            project_metadata={ "AB": { "Single cell platform":
                                       "10xGenomics Chromium 3'v3",
                                       "Library type": "RNA-seq",
                                       "Organism": "Human" } },
            top_dir=self.dirn)
        mockdir.create()
        # Add 10x Genomics pipeline outputs
        project = AnalysisProject(os.path.join(mockdir.dirn, "AB"))
        UpdateAnalysisProject(project).add_qc_outputs()
        UpdateAnalysisProject(project).add_cellranger_multi_outputs(
            sample_names=("AB1",))
        # Make a target directory
        target_dir = os.path.join(self.dirn, "shared")
        os.makedirs(target_dir)
        # Do data transfer (--include_10x_outputs)
        self.assertEqual(transfer_data(
            [target_dir,
             os.path.join(mockdir.dirn, "AB"),
             "--include_10x_outputs"]), 0)
        # Check transferred artefacts
        print(os.listdir(target_dir))
        expected_files = ("AB1_S1_R1_001.fastq.gz",
                          "AB1_S1_R2_001.fastq.gz",
                          "AB2_S2_R2_001.fastq.gz",
                          "AB2_S2_R1_001.fastq.gz",
                          "AB.chksums",
                          "cellranger_multi.AB.170901_M00879_0087_000000000-AGEW9.tgz")
        for f in expected_files:
            self.assertTrue(os.path.exists(os.path.join(target_dir, f)),
                            f"'{f}': missing, should be present")
        for f in os.listdir(target_dir):
            self.assertTrue(f in expected_files,
                            f"'{f}': present, but not expected")

    def test_transfer_data_fastqs_with_10x_cloupe_files(self):
        """
        transfer_data: copy Fastqs and 10xGenomics '.cloupe' files
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901",
                       "run_number": "89" },
            project_metadata={ "AB": { "Single cell platform":
                                       "10xGenomics Chromium 3'v3",
                                       "Library type": "RNA-seq",
                                       "Organism": "Human" } },
            top_dir=self.dirn)
        mockdir.create()
        # Add 10x Genomics pipeline outputs
        project = AnalysisProject(os.path.join(mockdir.dirn, "AB"))
        UpdateAnalysisProject(project).add_qc_outputs()
        UpdateAnalysisProject(project).add_cellranger_count_outputs()
        # Make a target directory
        target_dir = os.path.join(self.dirn, "shared")
        os.makedirs(target_dir)
        # Do data transfer (--include_cloupe_files)
        self.assertEqual(transfer_data(
            [target_dir,
             os.path.join(mockdir.dirn, "AB"),
             "--include_cloupe_files"]), 0)
        # Check transferred artefacts
        print(os.listdir(target_dir))
        expected_files = ("AB1_S1_R1_001.fastq.gz",
                          "AB1_S1_R2_001.fastq.gz",
                          "AB2_S2_R2_001.fastq.gz",
                          "AB2_S2_R1_001.fastq.gz",
                          "AB.chksums",
                          "10x_cloupe_files.zip")
        for f in expected_files:
            self.assertTrue(os.path.exists(os.path.join(target_dir, f)),
                            f"'{f}': missing, should be present")
        for f in os.listdir(target_dir):
            self.assertTrue(f in expected_files,
                            f"'{f}': present, but not expected")

    def test_transfer_data_fastqs_with_visium_images(self):
        """
        transfer_data: copy Fastqs and 10xGenomics Visium images
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901",
                       "run_number": "89" },
            project_metadata={ "AB": { "Single cell platform":
                                       "10xGenomics Chromium 3'v3",
                                       "Library type": "RNA-seq",
                                       "Organism": "Human" } },
            top_dir=self.dirn)
        mockdir.create()
        # Add mock 'Visium_images' dir to project
        visium_images = os.path.join(mockdir.dirn, "AB", "Visium_images")
        os.makedirs(visium_images)
        for img in ("Hires1.tif", "Hires.json"):
            with open(os.path.join(visium_images, img), "wt") as fp:
                fp.write("")
        # Make a target directory
        target_dir = os.path.join(self.dirn, "shared")
        os.makedirs(target_dir)
        # Do data transfer (--include_visium_images)
        self.assertEqual(transfer_data(
            [target_dir,
             os.path.join(mockdir.dirn, "AB"),
             "--include_visium_images"]), 0)
        # Check transferred artefacts
        print(os.listdir(target_dir))
        expected_files = ("AB1_S1_R1_001.fastq.gz",
                          "AB1_S1_R2_001.fastq.gz",
                          "AB2_S2_R2_001.fastq.gz",
                          "AB2_S2_R1_001.fastq.gz",
                          "AB.chksums",
                          "Visium_images.AB.170901_M00879_0087_000000000-AGEW9.tgz")
        for f in expected_files:
            self.assertTrue(os.path.exists(os.path.join(target_dir, f)),
                            f"'{f}': missing, should be present")
        for f in os.listdir(target_dir):
            self.assertTrue(f in expected_files,
                            f"'{f}': present, but not expected")

    def test_transfer_data_subdir_run_id(self):
        """
        transfer_data: copy to subdir with run ID
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901",
                       "run_number": "89" },
            project_metadata={ "AB": { "Library type": "RNA-seq",
                                       "Organism": "Human" } },
            top_dir=self.dirn)
        mockdir.create()
        # Make a target directory
        target_dir = os.path.join(self.dirn, "shared")
        os.makedirs(target_dir)
        # Do data transfer (--subdir run_id)
        self.assertEqual(transfer_data(
            [target_dir,
             os.path.join(mockdir.dirn, "AB"),
             "--subdir", "run_id"]), 0)
        # Update the target directory
        target_dir = os.path.join(target_dir, "MISEQ_170901.89-AB")
        # Check transferred artefacts
        print(os.listdir(target_dir))
        expected_files = ("AB1_S1_R1_001.fastq.gz",
                          "AB1_S1_R2_001.fastq.gz",
                          "AB2_S2_R2_001.fastq.gz",
                          "AB2_S2_R1_001.fastq.gz",
                          "AB.chksums")
        for f in expected_files:
            self.assertTrue(os.path.exists(os.path.join(target_dir, f)),
                            f"'{f}': missing, should be present")
        for f in os.listdir(target_dir):
            self.assertTrue(f in expected_files,
                            f"'{f}': present, but not expected")

    def test_transfer_data_subdir_random_bin(self):
        """
        transfer_data: copy to random bin subdir
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901",
                       "run_number": "89" },
            project_metadata={ "AB": { "Library type": "RNA-seq",
                                       "Organism": "Human" } },
            top_dir=self.dirn)
        mockdir.create()
        # Make a target directory
        target_dir = os.path.join(self.dirn, "shared")
        os.makedirs(target_dir)
        # Add an empty subdir
        subdir = "Ae2Eixor"
        os.makedirs(os.path.join(target_dir, subdir))
        # Do data transfer (--subdir random_bin)
        self.assertEqual(transfer_data(
            [target_dir,
             os.path.join(mockdir.dirn, "AB"),
             "--subdir", "random_bin"]), 0)
        # Update the target directory
        target_dir = os.path.join(target_dir, subdir)
        # Check transferred artefacts
        print(os.listdir(target_dir))
        expected_files = ("AB1_S1_R1_001.fastq.gz",
                          "AB1_S1_R2_001.fastq.gz",
                          "AB2_S2_R2_001.fastq.gz",
                          "AB2_S2_R1_001.fastq.gz",
                          "AB.chksums")
        for f in expected_files:
            self.assertTrue(os.path.exists(os.path.join(target_dir, f)),
                            f"'{f}': missing, should be present")
        for f in os.listdir(target_dir):
            self.assertTrue(f in expected_files,
                            f"'{f}': present, but not expected")
