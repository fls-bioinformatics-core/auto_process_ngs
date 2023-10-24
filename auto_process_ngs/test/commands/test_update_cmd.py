#######################################################################
# Tests for update_cmd.py module
#######################################################################

import unittest
import os
import tempfile
import shutil
from auto_process_ngs.analysis import AnalysisProject
from auto_process_ngs.auto_processor import AutoProcess
from auto_process_ngs.metadata import AnalysisProjectInfo
from auto_process_ngs.metadata import AnalysisProjectQCDirInfo
from auto_process_ngs.mock import MockAnalysisDirFactory
from auto_process_ngs.mock import UpdateAnalysisProject
from auto_process_ngs.commands.update_cmd import update

# Unit tests

class TestUpdate(unittest.TestCase):
    """
    Tests for the 'update' command
    """
    def setUp(self):
        # Project metadata items mapping
        self.metadata_map = {
            "User": "user",
            "PI": "PI",
            "Library": "library_type",
            "SC_Platform": "single_cell_platform",
            "Organism": "organism",
            "Comments": "comments"
        }
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestAutoProcess')
        # Store original location
        self.pwd = os.getcwd()
        # Move to working directory
        os.chdir(self.dirn)

    def tearDown(self):
        # Return to original dir
        os.chdir(self.pwd)
        # Remove the temporary test directory
        shutil.rmtree(self.dirn)

    def test_update_relocated_analysis_dir(self):
        """
        update: handle stored paths for relocated analysis dir
        """
        # Make an auto-process directory
        top_dir1 = os.path.join(self.dirn,"original")
        os.mkdir(top_dir1)
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '231021_A00879_0087_000000000-AGEW9',
            'novaseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=top_dir1)
        mockdir.create(no_project_dirs=True)
        original_path = mockdir.dirn
        # Relocate
        top_dir2 = os.path.join(self.dirn,"new")
        os.mkdir(top_dir2)
        new_path = os.path.join(top_dir2,os.path.basename(mockdir.dirn))
        os.rename(original_path,new_path)
        # Set up AutoProcess instance
        ap = AutoProcess(new_path)
        # Check metadata items pre-update
        self.assertEqual(ap.analysis_dir,new_path)
        self.assertEqual(ap.params.analysis_dir,original_path)
        self.assertEqual(ap.params.sample_sheet,
                         os.path.join(original_path,
                                      "custom_SampleSheet.csv"))
        self.assertEqual(ap.params.primary_data_dir,
                         os.path.join(original_path,"primary_data"))
        # Update the metadata
        update(ap)
        # Reload and check metadata items post-update
        ap = AutoProcess(new_path)
        self.assertEqual(ap.analysis_dir,new_path)
        self.assertEqual(ap.params.analysis_dir,new_path)
        self.assertEqual(ap.params.sample_sheet,
                         os.path.join(new_path,
                                      "custom_SampleSheet.csv"))
        self.assertEqual(ap.params.primary_data_dir,
                         os.path.join(new_path,"primary_data"))

    def test_update_relocated_analysis_dir_with_qc(self):
        """
        update: handle stored paths for QC in relocated analysis dir
        """
        # List of projects
        project_list = ("AB","CDE","undetermined")
        # Make an auto-process directory with projects
        top_dir1 = os.path.join(self.dirn,"original")
        os.mkdir(top_dir1)
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '231021_A00879_0087_000000000-AGEW9',
            'novaseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=top_dir1)
        mockdir.create()
        original_path = mockdir.dirn
        # Add mock QC outputs to projects
        for project_name in project_list:
            p = AnalysisProject(os.path.join(mockdir.dirn,project_name))
            UpdateAnalysisProject(p).add_qc_outputs()
        # Relocate
        top_dir2 = os.path.join(self.dirn,"new")
        os.mkdir(top_dir2)
        new_path = os.path.join(top_dir2,os.path.basename(mockdir.dirn))
        os.rename(original_path,new_path)
        # Set up AutoProcess instance
        ap = AutoProcess(new_path)
        # Check metadata items items pre-update
        self.assertEqual(ap.analysis_dir,new_path)
        self.assertEqual(ap.params.analysis_dir,original_path)
        self.assertEqual(ap.params.sample_sheet,
                         os.path.join(original_path,
                                      "custom_SampleSheet.csv"))
        self.assertEqual(ap.params.primary_data_dir,
                         os.path.join(original_path,"primary_data"))
        # Check QC paths in projects pre-update
        for proj in ap.get_analysis_projects():
            qc_dirs = proj.qc_dirs
            self.assertTrue(len(qc_dirs) == 1)
            qc_info = AnalysisProjectQCDirInfo(os.path.join(new_path,
                                                            proj.name,
                                                            qc_dirs[0],
                                                            "qc.info"))
            self.assertEqual(qc_info.fastq_dir,
                             os.path.join(original_path,
                                          proj.name,
                                          "fastqs"))
        # Update the metadata
        update(ap)
        # Reload and check metadata items post-update
        ap = AutoProcess(new_path)
        self.assertEqual(ap.analysis_dir,new_path)
        self.assertEqual(ap.params.analysis_dir,new_path)
        self.assertEqual(ap.params.sample_sheet,
                         os.path.join(new_path,
                                      "custom_SampleSheet.csv"))
        self.assertEqual(ap.params.primary_data_dir,
                         os.path.join(new_path,"primary_data"))
        # Check QC paths in projects post-update
        for proj in ap.get_analysis_projects():
            qc_dirs = proj.qc_dirs
            self.assertTrue(len(qc_dirs) == 1)
            qc_info = AnalysisProjectQCDirInfo(os.path.join(new_path,
                                                            proj.name,
                                                            qc_dirs[0],
                                                            "qc.info"))
            self.assertEqual(qc_info.fastq_dir,
                             os.path.join(new_path,
                                          proj.name,
                                          "fastqs"))

    def test_update_project_metadata_changed(self):
        """
        update: handle synchronising project metadata
        """
        # Metadata for projects
        project_metadata = {
            "AB": {
                "User": "Alan Bailey",
                "PI": "Archie Ballard",
                "Library": "RNA-seq",
                "Organism": "Human"
            },
            "CDE": {
                "User": "Charles Edwards",
                "PI": "Christian Eggars",
                "Library": "ChIP-seq",
                "Organism": "Mouse"
            }
        }
        # Make an auto-process directory with projects
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '231021_A00879_0087_000000000-AGEW9',
            'novaseq',
            project_metadata=project_metadata,
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create()
        # Set up AutoProcess instance
        ap = AutoProcess(mockdir.dirn)
        # Check metadata items in projects.info pre-update
        ap_project_metadata = ap.load_project_metadata()
        for pname in project_metadata:
            for item in project_metadata[pname]:
                self.assertEqual(ap_project_metadata.lookup(pname)[item],
                                 project_metadata[pname][item])
        # Check metadata items in projects pre-update
        for pname in project_metadata:
            p = ap.get_analysis_projects(pname)[0]
            for item in project_metadata[pname]:
                print(item)
                project_item = self.metadata_map[item]
                expected_value = project_metadata[pname][item]
                self.assertEqual(p.info[project_item],
                                 expected_value)
        # Overwrite projects.info and update metadata
        project_metadata["AB"]["Organism"] = "Mouse"
        project_metadata["AB"]["Comments"] = "1% PhiX spiked in"
        project_metadata["CDE"]["Comments"] = "1% PhiX spiked in"
        with open(os.path.join(mockdir.dirn,"projects.info"),'wt') as fp:
            fp.write("""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1,AB2\tAlan Bailey\tRNA-seq\t.\tMouse\tArchie Ballard\t1% PhiX spiked in
CDE\tCDE3,CDE4\tCharles Edwards\tChIP-seq\t.\tMouse\tChristian Eggars\t1% PhiX spiked in
""")
        update(ap)
        # Reload and confirm the updates in projects.info
        ap = AutoProcess(mockdir.dirn)
        ap_project_metadata = ap.load_project_metadata()
        for pname in project_metadata:
            for item in project_metadata[pname]:
                self.assertEqual(ap_project_metadata.lookup(pname)[item],
                                 project_metadata[pname][item])
        # Confirm the updates in the projects
        for pname in project_metadata:
            p = ap.get_analysis_projects(pname)[0]
            for item in project_metadata[pname]:
                print(item)
                project_item = self.metadata_map[item]
                expected_value = project_metadata[pname][item]
                self.assertEqual(p.info[project_item],
                                 expected_value)

    def test_update_sample_list_changed(self):
        """
        update: handle synchronising sample list metadata
        """
        # Metadata for projects
        project_metadata = {
            "AB": {
                "User": "Alan Bailey",
                "PI": "Archie Ballard",
                "Library": "RNA-seq",
                "Organism": "Human"
            },
            "CDE": {
                "User": "Charles Edwards",
                "PI": "Christian Eggars",
                "Library": "ChIP-seq",
                "Organism": "Mouse"
            }
        }
        # Make an auto-process directory with projects
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '231021_A00879_0087_000000000-AGEW9',
            'novaseq',
            project_metadata=project_metadata,
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create()
        # Set up AutoProcess instance
        ap = AutoProcess(mockdir.dirn)
        # Check sample lists in projects.info
        expected_samples = {
            "AB": ["AB1","AB2"],
            "CDE": ["CDE3","CDE4"]
        }
        ap_project_metadata = ap.load_project_metadata()
        for pname in expected_samples:
            self.assertEqual(ap_project_metadata.lookup(pname)["Samples"],
                             ','.join(expected_samples[pname]))
        # Check sample lists in project metadata
        for pname in expected_samples:
            p = ap.get_analysis_projects(pname)[0]
            self.assertEqual(p.info.samples,
                             "%s sample%s (%s)" %
                             (len(expected_samples[pname]),
                              's' if len(expected_samples[pname]) != 1 else '',
                              ', '.join(expected_samples[pname])))
        # Remove samples from projects
        for fq in ("AB/fastqs/AB1_S1_R1_001.fastq.gz",
                   "AB/fastqs/AB1_S1_R2_001.fastq.gz",
                   "CDE/fastqs/CDE4_S4_R1_001.fastq.gz",
                   "CDE/fastqs/CDE4_S4_R2_001.fastq.gz"):
            os.remove(os.path.join(mockdir.dirn,fq))
        expected_samples = {
            "AB": ["AB2"],
            "CDE": ["CDE3"]
        }
        # Update sample lists in metadata
        update(ap)
        # Reload and check updated sample lists in projects.info
        ap = AutoProcess(mockdir.dirn)
        ap_project_metadata = ap.load_project_metadata()
        for pname in expected_samples:
            self.assertEqual(ap_project_metadata.lookup(pname)["Samples"],
                             ','.join(expected_samples[pname]))
        # Check updated sample lists in project metadata
        for pname in expected_samples:
            p = ap.get_analysis_projects(pname)[0]
            self.assertEqual(p.info.samples,
                             "%s sample%s (%s)" %
                             (len(expected_samples[pname]),
                              's' if len(expected_samples[pname]) != 1 else '',
                              ', '.join(expected_samples[pname])))

    def test_update_paired_end_changed(self):
        """
        update: handle paired-end update in project metadata
        """
        # Make an auto-process directory with projects
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '231021_A00879_0087_000000000-AGEW9',
            'novaseq',
            #project_metadata=project_metadata,
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create()
        # Check paired-end metadata in projects
        # NB check the metadata directly (not via other objects)
        expected_paired_end = {
            "AB": True,
            "CDE": True
        }
        for pname in expected_paired_end:
            project_metadata = AnalysisProjectInfo(
                os.path.join(mockdir.dirn,pname,"README.info"))
            self.assertEqual(project_metadata.paired_end,
                             expected_paired_end[pname])
        # Remove R2 Fastqs from one of the projects
        for fq in ("AB/fastqs/AB1_S1_R2_001.fastq.gz",
                   "AB/fastqs/AB2_S2_R2_001.fastq.gz"):
            os.remove(os.path.join(mockdir.dirn,fq))
        # Set up AutoProcess instance and update metadata
        ap = AutoProcess(mockdir.dirn)
        update(ap)
        # Reload and check update paired-end metadata in projects
        ap = AutoProcess(mockdir.dirn)
        expected_paired_end = {
            "AB": False,
            "CDE": True
        }
        for pname in expected_paired_end:
            project_metadata = AnalysisProjectInfo(
                os.path.join(mockdir.dirn,pname,"README.info"))
            self.assertEqual(project_metadata.paired_end,
                             expected_paired_end[pname])
