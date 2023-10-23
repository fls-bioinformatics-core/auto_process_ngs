#######################################################################
# Tests for autoprocessor.py module
#######################################################################

import unittest
import tempfile
import shutil
import os
from bcftbx.mock import MockIlluminaRun
from auto_process_ngs.auto_processor import AutoProcess
from auto_process_ngs.analysis import AnalysisProject
from auto_process_ngs.metadata import AnalysisProjectInfo
from auto_process_ngs.mock import MockAnalysisDirFactory
from auto_process_ngs.mock import MockAnalysisProject

# Unit tests

class TestAutoProcess(unittest.TestCase):

    def setUp(self):
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

    def test_analysis_dir_path(self):
        """AutoProcess: analysis dir path is absolute and normalized
        """
        # Create mock Illumina run directory
        mock_illumina_run = MockIlluminaRun(
            '160621_M00879_0087_000000000-AGEW9',
            'miseq',
            top_dir=self.dirn)
        mock_illumina_run.create()
        # Set up new AutoProcess instance
        ap = AutoProcess()
        self.assertEqual(ap.analysis_dir,None)
        # Make a mock analysis dir
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_M00879_0087_000000000-AGEW9',
            'miseq',
            top_dir=self.dirn)
        mockdir.create()
        # Create Autoprocess instances from different
        # forms of path and check stored value
        rel_path = "160621_M00879_0087_000000000-AGEW9_analysis"
        abs_path = os.path.join(self.dirn,rel_path)
        rel_unnormalised = os.path.join("..",
                                        os.path.basename(self.dirn),
                                        rel_path)
        abs_unnormalised = os.path.join(self.dirn,
                                        rel_unnormalised)
        ap = AutoProcess(analysis_dir=abs_path)
        self.assertEqual(ap.analysis_dir,abs_path)
        ap = AutoProcess(analysis_dir=rel_path)
        self.assertEqual(ap.analysis_dir,abs_path)
        ap = AutoProcess(analysis_dir=abs_unnormalised)
        self.assertEqual(ap.analysis_dir,abs_path)
        ap = AutoProcess(analysis_dir=rel_unnormalised)
        self.assertEqual(ap.analysis_dir,abs_path)

    def test_set_log_dir(self):
        """AutoProcess.set_log_dir: sets correct log directory
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_M00879_0087_000000000-AGEW9',
            'miseq',
            top_dir=self.dirn)
        mockdir.create()
        # Load into AutoProcess object
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Check the initial log directory
        self.assertEqual(ap.log_dir,
                         os.path.join(mockdir.dirn,
                                      "logs"))
        # Set a log subdir
        subdir = ap.get_log_subdir("testing")
        self.assertEqual(subdir,"001_testing")
        ap.set_log_dir("001_testing")
        self.assertEqual(ap.log_dir,
                         os.path.join(mockdir.dirn,
                                      "logs",
                                      "001_testing"))
        self.assertEqual(ap.log_path("test.log"),
                         os.path.join(mockdir.dirn,
                                      "logs",
                                      "001_testing",
                                      "test.log"))
        self.assertEqual(ap.log_path("test_1",
                                     "test1.log"),
                         os.path.join(mockdir.dirn,
                                      "logs",
                                      "001_testing",
                                      "test_1",
                                      "test1.log"))
        # Set a second log subdir
        subdir = ap.get_log_subdir("testing")
        self.assertEqual(subdir,"002_testing")
        ap.set_log_dir("002_testing")
        self.assertEqual(ap.log_dir,
                         os.path.join(mockdir.dirn,
                                      "logs",
                                      "002_testing"))
        self.assertEqual(ap.log_path("test.log"),
                         os.path.join(mockdir.dirn,
                                      "logs",
                                      "002_testing",
                                      "test.log"))
        self.assertEqual(ap.log_path("test_2",
                                     "test2.log"),
                         os.path.join(mockdir.dirn,
                                      "logs",
                                      "002_testing",
                                      "test_2",
                                      "test2.log"))

class TestAutoProcessUpdateMetadata(unittest.TestCase):
    """
    Tests for the 'update_metadata' method
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

    def test_update_metadata_relocated_analysis_dir(self):
        """
        AutoProcess.update_metadata: handle relocated analysis dir
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
        ap.update_metadata()
        # Check metadata items post-update
        self.assertEqual(ap.analysis_dir,new_path)
        self.assertEqual(ap.params.analysis_dir,new_path)
        self.assertEqual(ap.params.sample_sheet,
                         os.path.join(new_path,
                                      "custom_SampleSheet.csv"))
        self.assertEqual(ap.params.primary_data_dir,
                         os.path.join(new_path,"primary_data"))

    def test_update_metadata_project_metadata_changed(self):
        """
        AutoProcess.update_metadata: handle project metadata update
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
        ap.update_metadata()
        # Confirm the updates in projects.info
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

    def test_update_metadata_sample_list_changed(self):
        """
        AutoProcess.update_metadata: handle sample list update
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
        ap.update_metadata()
        # Check updated sample lists in projects.info
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

    def test_update_metadata_paired_end_changed(self):
        """
        AutoProcess.update_metadata: handle paired-end update
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
        ap.update_metadata()
        # Check update paired-end metadata in projects
        expected_paired_end = {
            "AB": False,
            "CDE": True
        }
        for pname in expected_paired_end:
            project_metadata = AnalysisProjectInfo(
                os.path.join(mockdir.dirn,pname,"README.info"))
            self.assertEqual(project_metadata.paired_end,
                             expected_paired_end[pname])

class TestAutoProcessMakeProjectMetadataFile(unittest.TestCase):
    """
    Test for the 'make_project_metadata_file' method
    """
    def setUp(self):
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

    def test_make_project_metadata_file(self):
        """
        AutoProcess.make_project_metadata_file: new 'projects.info' from bcl2fastq output
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Remove the projects.info file
        os.remove(os.path.join(mockdir.dirn,"projects.info"))
        # Create a new projects.info file
        AutoProcess(mockdir.dirn).make_project_metadata_file()
        # Check outputs
        self.assertTrue(os.path.exists(
            os.path.join(mockdir.dirn,"projects.info")))
        with open(os.path.join(mockdir.dirn,"projects.info"),'rt') as fp:
            self.assertEqual(fp.read(),
                             """#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1,AB2\t.\t.\t.\t.\t.\t.
CDE\tCDE3,CDE4\t.\t.\t.\t.\t.\t.
""")

    def test_make_project_metadata_file_no_bcl2fastq_output(self):
        """
        AutoProcess.make_project_metadata_file: new 'projects.info' (no bcl2fastq output)
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Remove the projects.info file and the bcl2fastq output dir
        os.remove(os.path.join(mockdir.dirn,"projects.info"))
        shutil.rmtree(os.path.join(mockdir.dirn,"bcl2fastq"))
        # Create a new projects.info file
        AutoProcess(mockdir.dirn).make_project_metadata_file()
        # Check outputs
        self.assertTrue(os.path.exists(
            os.path.join(mockdir.dirn,"projects.info")))
        with open(os.path.join(mockdir.dirn,"projects.info"),'rt') as fp:
            self.assertEqual(fp.read(),
                             """#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
""")

    def test_make_project_metadata_file_exception_if_file_exists(self):
        """
        AutoProcess.make_project_metadata_file: raise exception if 'projects.info' already exists
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Raise exception if projects.info file exists already
        self.assertRaises(Exception,
                          AutoProcess(mockdir.dirn).make_project_metadata_file)

class TestAutoProcessUpdateProjectMetadataFile(unittest.TestCase):
    """
    Test for the 'update_project_metadata_file' method
    """
    def setUp(self):
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

    def test_update_project_metadata_file_empty_from_bcl2fastq_output(self):
        """
        AutoProcess.update_project_metadata_file: update 'empty' file from bcl2fastq output
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Create empty projects.info file
        with open(os.path.join(mockdir.dirn,"projects.info"),'wt') as fp:
            fp.write("#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments\n")
        # Update the projects.info file
        AutoProcess(mockdir.dirn).update_project_metadata_file()
        # Check output
        with open(os.path.join(mockdir.dirn,"projects.info"),'rt') as fp:
            self.assertEqual(fp.read(),
                             """#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1,AB2\t.\t.\t.\t.\t.\t.
CDE\tCDE3,CDE4\t.\t.\t.\t.\t.\t.
""")

    def test_update_project_metadata_file_partial_from_bcl2fastq_output(self):
        """
        AutoProcess.update_project_metadata_file: update partial file from bcl2fastq output
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Create projects.info file with one project already listed
        with open(os.path.join(mockdir.dirn,"projects.info"),'wt') as fp:
            fp.write("#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments\nCDE\tCDE3,CDE4\t.\t.\t.\t.\t.\tKeep me")
        # Update the projects.info file
        AutoProcess(mockdir.dirn).update_project_metadata_file()
        # Check output
        with open(os.path.join(mockdir.dirn,"projects.info"),'rt') as fp:
            self.assertEqual(fp.read(),
                             """#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1,AB2\t.\t.\t.\t.\t.\t.
CDE\tCDE3,CDE4\t.\t.\t.\t.\t.\tKeep me
""")

    def test_update_project_metadata_file_missing_from_bcl2fastq_output(self):
        """
        AutoProcess.update_project_metadata_file: make missing file and populate from bcl2fastq output
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Remove projects.info file
        os.remove(os.path.join(mockdir.dirn,"projects.info"))
        # Update the projects.info file
        AutoProcess(mockdir.dirn).update_project_metadata_file()
        # Check output
        with open(os.path.join(mockdir.dirn,"projects.info"),'rt') as fp:
            self.assertEqual(fp.read(),
                             """#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1,AB2\t.\t.\t.\t.\t.\t.
CDE\tCDE3,CDE4\t.\t.\t.\t.\t.\t.
""")

    def test_update_project_metadata_file_comment_out_missing_project(self):
        """
        AutoProcess.update_project_metadata_file: missing project is commented out
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Create projects.info file with one project already listed
        with open(os.path.join(mockdir.dirn,"projects.info"),'wt') as fp:
            fp.write("#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments\nFG\tFG5,FG6\t.\t.\t.\t.\t.\tKeep me")
        # Update the projects.info file
        AutoProcess(mockdir.dirn).update_project_metadata_file()
        # Check output - missing project kept but commented out
        with open(os.path.join(mockdir.dirn,"projects.info"),'rt') as fp:
            print(fp.read())
        with open(os.path.join(mockdir.dirn,"projects.info"),'rt') as fp:
            self.assertEqual(fp.read(),
                             """#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1,AB2\t.\t.\t.\t.\t.\t.
CDE\tCDE3,CDE4\t.\t.\t.\t.\t.\t.
#FG\tFG5,FG6\t.\t.\t.\t.\t.\tKeep me
""")

    def test_update_project_metadata_file_uncomment_existing_project(self):
        """
        AutoProcess.update_project_metadata_file: existing project is uncommented
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Create projects.info file with one project already listed
        with open(os.path.join(mockdir.dirn,"projects.info"),'wt') as fp:
            fp.write("#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments\n#CDE\tCDE3,CDE4\t.\t.\t.\t.\t.\tKeep me")
        # Update the projects.info file
        AutoProcess(mockdir.dirn).update_project_metadata_file()
        # Check output - missing project kept but commented out
        with open(os.path.join(mockdir.dirn,"projects.info"),'rt') as fp:
            self.assertEqual(fp.read(),
                             """#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1,AB2\t.\t.\t.\t.\t.\t.
CDE\tCDE3,CDE4\t.\t.\t.\t.\t.\tKeep me
""")

    def test_update_project_metadata_file_dont_comment_missing_project_when_dir_is_present(self):
        """
        AutoProcess.update_project_metadata_file: don't comment out 'missing' project when dir is present
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Create projects.info file with one project already listed
        with open(os.path.join(mockdir.dirn,"projects.info"),'wt') as fp:
            fp.write("#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments\nFG\tFG5,FG6\t.\t.\t.\t.\t.\tKeep me")
        # Create the corresponding project
        project = MockAnalysisProject('FG',('FG5_S1_R1_001.fastq.gz',
                                            'FG6_S1_R1_001.fastq.gz'))
        project.create(top_dir=mockdir.dirn)
        # Update the projects.info file
        AutoProcess(mockdir.dirn).update_project_metadata_file()
        # Check output - missing project kept but commented out
        with open(os.path.join(mockdir.dirn,"projects.info"),'rt') as fp:
            print(fp.read())
        with open(os.path.join(mockdir.dirn,"projects.info"),'rt') as fp:
            self.assertEqual(fp.read(),
                             """#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1,AB2\t.\t.\t.\t.\t.\t.
CDE\tCDE3,CDE4\t.\t.\t.\t.\t.\t.
FG\tFG5,FG6\t.\t.\t.\t.\t.\tKeep me
""")

    def test_update_project_metadata_file_dont_uncomment_missing_project_when_dir_is_present(self):
        """
        AutoProcess.update_project_metadata_file: don't uncomment 'missing' project when dir is present
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Create projects.info file with one project already listed
        with open(os.path.join(mockdir.dirn,"projects.info"),'wt') as fp:
            fp.write("#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments\n#FG\tFG5,FG6\t.\t.\t.\t.\t.\tKeep me")
        # Create the corresponding project
        project = MockAnalysisProject('FG',('FG5_S1_R1_001.fastq.gz',
                                            'FG6_S1_R1_001.fastq.gz'))
        project.create(top_dir=mockdir.dirn)
        # Update the projects.info file
        AutoProcess(mockdir.dirn).update_project_metadata_file()
        # Check output - missing project kept but commented out
        with open(os.path.join(mockdir.dirn,"projects.info"),'rt') as fp:
            print(fp.read())
        with open(os.path.join(mockdir.dirn,"projects.info"),'rt') as fp:
            self.assertEqual(fp.read(),
                             """#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1,AB2\t.\t.\t.\t.\t.\t.
CDE\tCDE3,CDE4\t.\t.\t.\t.\t.\t.
#FG\tFG5,FG6\t.\t.\t.\t.\t.\tKeep me
""")

class TestAutoProcessGetAnalysisProjectsMethod(unittest.TestCase):
    """
    Tests for the 'get_analysis_projects' method
    """
    def setUp(self):
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

    def test_no_projects_dot_info_no_project_dirs(self):
        """AutoProcess.get_analysis_projects: no project dirs (no projects.info)
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Remove the projects.info file
        os.remove(os.path.join(mockdir.dirn,"projects.info"))
        print(os.listdir(mockdir.dirn))
        # No projects should be listed
        projects = AutoProcess(mockdir.dirn).get_analysis_projects()
        self.assertEqual(projects,[])

    def test_projects_dot_info_no_project_dirs(self):
        """AutoProcess.get_analysis_projects: no project dirs
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # No projects should be listed
        projects = AutoProcess(mockdir.dirn).get_analysis_projects()
        self.assertEqual(projects,[])

    def test_with_project_dirs(self):
        """AutoProcess.get_analysis_projects: project dirs exist
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create()
        # List the projects
        projects = AutoProcess(mockdir.dirn).get_analysis_projects()
        expected = ('AB','CDE','undetermined')
        self.assertEqual(len(projects),len(expected))
        for p in projects:
            self.assertTrue(isinstance(p,AnalysisProject))
            self.assertTrue(p.name in expected)
        for p in expected:
            matched_projects = [x for x in projects if x.name == p]
            self.assertEqual(len(matched_projects),1)

    def test_with_project_dirs_no_projects_dot_info(self):
        """AutoProcess.get_analysis_projects: project dirs exist (no projects.info)
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create()
        # Remove the projects.info file
        os.remove(os.path.join(mockdir.dirn,"projects.info"))
        # List the projects
        projects = AutoProcess(mockdir.dirn).get_analysis_projects()
        expected = ('undetermined',)
        self.assertEqual(len(projects),1)
        for p in projects:
            self.assertTrue(isinstance(p,AnalysisProject))
            self.assertTrue(p.name in expected)
        for p in expected:
            matched_projects = [x for x in projects if x.name == p]
            self.assertEqual(len(matched_projects),1)

    def test_ignore_commented_projects(self):
        """AutoProcess.get_analysis_projects: ignore commented projects
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create()
        # Update the projects.info file
        projects_info = os.path.join(mockdir.dirn,"projects.info")
        with open(projects_info,"w") as fp:
            fp.write(
"""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
#AB\tAB1,AB2\tAlan Brown\tRNA-seq\t.\tHuman\tAudrey Benson\t1% PhiX
CDE\tCDE3,CDE4\tClive David Edwards\tChIP-seq\t.\tMouse\tClaudia Divine Eccleston\t1% PhiX
""")
        # List the projects
        projects = AutoProcess(mockdir.dirn).get_analysis_projects()
        expected = ('CDE','undetermined')
        self.assertEqual(len(projects),len(expected))
        for p in projects:
            self.assertTrue(isinstance(p,AnalysisProject))
            self.assertTrue(p.name in expected)
        for p in expected:
            matched_projects = [x for x in projects if x.name == p]
            self.assertEqual(len(matched_projects),1)

    def test_with_project_dirs_select_subset(self):
        """AutoProcess.get_analysis_projects: select subset of projects
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create()
        # List the projects
        projects = AutoProcess(mockdir.dirn).get_analysis_projects("C*")
        expected = ('CDE',)
        self.assertEqual(len(projects),len(expected))
        for p in projects:
            self.assertTrue(isinstance(p,AnalysisProject))
            self.assertTrue(p.name in expected)
        for p in expected:
            matched_projects = [x for x in projects if x.name == p]
            self.assertEqual(len(matched_projects),1)

class TestAutoProcessGetAnalysisProjectsFromDirsMethod(unittest.TestCase):
    """
    Tests for the 'get_analysis_projects_from_dirs' method
    """
    def setUp(self):
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

    def test_no_projects(self):
        """AutoProcess.get_analysis_projects_from_dirs: no project dirs
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # No projects should be listed
        projects = AutoProcess(mockdir.dirn).get_analysis_projects_from_dirs()
        self.assertEqual(projects,[])

    def test_with_projects(self):
        """AutoProcess.get_analysis_projects_from_dirs: fetches project dirs
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create()
        # List the projects
        projects = AutoProcess(mockdir.dirn).get_analysis_projects()
        expected = ('AB','CDE','undetermined')
        self.assertEqual(len(projects),len(expected))
        for p in projects:
            self.assertTrue(isinstance(p,AnalysisProject))
            self.assertTrue(p.name in expected)
        for p in expected:
            matched_projects = [x for x in projects if x.name == p]
            self.assertEqual(len(matched_projects),1)

    def test_with_projects_select_subset(self):
        """AutoProcess.get_analysis_projects_from_dirs: selects subset of projects
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create()
        # List the projects
        projects = AutoProcess(mockdir.dirn).get_analysis_projects_from_dirs("C*")
        expected = ('CDE',)
        self.assertEqual(len(projects),len(expected))
        for p in projects:
            self.assertTrue(isinstance(p,AnalysisProject))
            self.assertTrue(p.name in expected)
        for p in expected:
            matched_projects = [x for x in projects if x.name == p]
            self.assertEqual(len(matched_projects),1)

class TestAutoProcessPairedEndMethod(unittest.TestCase):
    """
    Tests for the 'paired_end' method
    """
    def setUp(self):
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

    def test_single_end(self):
        """AutoProcess.paired_end: check for single-ended data
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            paired_end=False,
            top_dir=self.dirn)
        mockdir.create()
        self.assertFalse(AutoProcess(mockdir.dirn).paired_end)

    def test_paired_end(self):
        """AutoProcess.paired_end: check for paired-ended data
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            paired_end=True,
            top_dir=self.dirn)
        mockdir.create()
        self.assertTrue(AutoProcess(mockdir.dirn).paired_end)

    def test_single_end_no_projects(self):
        """AutoProcess.paired_end: check for single-ended data (no project dirs)
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            paired_end=False,
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        self.assertFalse(AutoProcess(mockdir.dirn).paired_end)

    def test_paired_end_no_projects(self):
        """AutoProcess.paired_end: check for paired-ended data (no project dirs)
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            paired_end=True,
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        self.assertTrue(AutoProcess(mockdir.dirn).paired_end)

    def test_single_end_no_projects_no_aligned(self):
        """AutoProcess.paired_end: check for single-ended data (no project dirs, no unaligned)
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            paired_end=False,
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Remove the unaligned dir
        shutil.rmtree(os.path.join(mockdir.dirn,"bcl2fastq"))
        self.assertEqual(AutoProcess(mockdir.dirn).paired_end,None)

    def test_paired_end_no_projects_no_unaligned(self):
        """AutoProcess.paired_end: check for paired-ended data (no project dirs, no unaligned)
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            paired_end=True,
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Remove the unaligned dir
        shutil.rmtree(os.path.join(mockdir.dirn,"bcl2fastq"))
        self.assertEqual(AutoProcess(mockdir.dirn).paired_end,None)

    def test_single_end_extra_dir(self):
        """AutoProcess.paired_end: check for single-ended data (extra dir)
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            paired_end=False,
            top_dir=self.dirn)
        mockdir.create()
        # Add an extra project-like directory with Fastq pairs
        extra_dir = os.path.join(mockdir.dirn,'extra_project')
        os.mkdir(extra_dir)
        for fq in ('extra_R1.fastq','extra_R2.fastq,'):
            with open(os.path.join(extra_dir,fq),'wt') as fp:
                fp.write('')
        # Check that paired_end is false
        self.assertFalse(AutoProcess(mockdir.dirn).paired_end)

    def test_paired_end_extra_dir(self):
        """AutoProcess.paired_end: check for paired-ended data (extra dir)
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            paired_end=True,
            top_dir=self.dirn)
        mockdir.create()
        # Add an extra project-like directory with Fastqs
        extra_dir = os.path.join(mockdir.dirn,'extra_project')
        os.mkdir(extra_dir)
        for fq in ('extra1_R1.fastq','extra2_R1.fastq,'):
            with open(os.path.join(extra_dir,fq),'wt') as fp:
                fp.write('')
        # Check that paired_end is true
        self.assertTrue(AutoProcess(mockdir.dirn).paired_end)

class TestAutoProcessUndeterminedMethod(unittest.TestCase):
    """
    Tests for the 'undetermined' method
    """
    def setUp(self):
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

    def test_undetermined(self):
        """AutoProcess.undetermined: 'undetermined' exists
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            paired_end=True,
            top_dir=self.dirn)
        mockdir.create()
        # Check that undetermined is returned
        undetermined = AutoProcess(mockdir.dirn).undetermined()
        self.assertTrue(isinstance(undetermined,AnalysisProject))

    def test_undetermined_missing(self):
        """AutoProcess.undetermined: 'undetermined' missing
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            paired_end=True,
            top_dir=self.dirn)
        mockdir.create()
        # Remove undetermined
        shutil.rmtree(os.path.join(mockdir.dirn,'undetermined'))
        # Check that undetermined is returned
        undetermined = AutoProcess(mockdir.dirn).undetermined()
        self.assertEqual(undetermined,None)
