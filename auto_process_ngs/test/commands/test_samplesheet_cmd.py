#######################################################################
# Tests for samplesheet_cmd.py module
#######################################################################

import unittest
import tempfile
import shutil
import os
from bcftbx.IlluminaData import SampleSheet
from auto_process_ngs.auto_processor import AutoProcess
from auto_process_ngs.mock import MockAnalysisDir
from auto_process_ngs.commands.samplesheet_cmd import SampleSheetOperation
from auto_process_ngs.commands.samplesheet_cmd import samplesheet
from auto_process_ngs.commands.samplesheet_cmd import import_samplesheet

# Unit tests

class TestSampleSheetSetProject(unittest.TestCase):
    """
    Test the 'set_project' functionality of the 'samplesheet' command
    """
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestSampleSheetSetProject')
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
            os.chmod(os.path.dirname(name),0o755)
            os.chmod(name,0o655)
            os.remove(name)
        shutil.rmtree(self.dirn,onerror=del_rw)

    def test_samplesheet_set_project(self):
        """
        samplesheet: set project names
        """
        # Make an empty mock auto-process directory
        mockdir = MockAnalysisDir(
            '221219_SN00879_0087_000000000-AGEW9',
            'hiseq',
            no_undetermined=True,
            top_dir=self.dirn)
        mockdir.create()
        # Update content for sample sheet
        sample_sheet_data = u"""[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,AB1,Sample_AB1,,,N701,CGATGTAT,N501,TCTTTCCC,AB1,
2,AB2,Sample_AB2,,,N702,TGACCAAT,N502,TCTTTCCC,AB2,
3,CD1,Sample_CD1,,,N701,CGATGTAT,N501,TCTTTCCC,CD1,
4,FG2,Sample_FG2,,,N702,TGACCAAT,N502,TCTTTCCC,FG2,
"""
        sample_sheet_file = os.path.join(mockdir.dirn,
                                         "custom_SampleSheet.csv")
        with open(sample_sheet_file,'wt') as fp:
            fp.write(sample_sheet_data)
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Update sample sheet project for lane 3
        ap.samplesheet(SampleSheetOperation.SET_PROJECT,
                       "CarlDewey",
                       lanes=[3])
        for line in SampleSheet(sample_sheet_file):
            if line['Lane'] == 3:
                self.assertEqual(line["Sample_Project"],"CarlDewey")
        # Update sample sheet project matching on sample ID
        ap.samplesheet(SampleSheetOperation.SET_PROJECT,
                       "AndrewBloggs",
                       where=("SAMPLE_ID",
                              "AB*"))
        for line in SampleSheet(sample_sheet_file):
            if line['Lane'] in (1,2):
                self.assertEqual(line["Sample_Project"],"AndrewBloggs")

class TestSampleSheetSetSampleName(unittest.TestCase):
    """
    Test the 'set_sample_name' functionality of the 'samplesheet' command
    """
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestSampleSheetSetProject')
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
            os.chmod(os.path.dirname(name),0o755)
            os.chmod(name,0o655)
            os.remove(name)
        shutil.rmtree(self.dirn,onerror=del_rw)

    def test_samplesheet_set_sample_name(self):
        """
        samplesheet: set sample names
        """
        # Make an empty mock auto-process directory
        mockdir = MockAnalysisDir(
            '221219_SN00879_0087_000000000-AGEW9',
            'hiseq',
            no_undetermined=True,
            top_dir=self.dirn)
        mockdir.create()
        # Update content for sample sheet
        sample_sheet_data = u"""[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,AB1,Sample_AB1,,,N701,CGATGTAT,N501,TCTTTCCC,AndrewBloggs,
2,AB2,Sample_AB2,,,N702,TGACCAAT,N502,TCTTTCCC,AndrewBloggs,
3,CD1,Sample_CD1,,,N701,CGATGTAT,N501,TCTTTCCC,CarlDewey,
4,FG2,Sample_FG2,,,N702,TGACCAAT,N502,TCTTTCCC,FilipeGreer,
"""
        sample_sheet_file = os.path.join(mockdir.dirn,
                                         "custom_SampleSheet.csv")
        with open(sample_sheet_file,'wt') as fp:
            fp.write(sample_sheet_data)
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Update sample sheet name for lane 3
        ap.samplesheet(SampleSheetOperation.SET_SAMPLE_NAME,
                       "SAMPLE_ID",
                       lanes=[3])
        for line in SampleSheet(sample_sheet_file):
            if line['Lane'] == 3:
                self.assertEqual(line["Sample_Name"],"CD1")
        # Update sample sheet name matching on project name
        ap.samplesheet(SampleSheetOperation.SET_SAMPLE_NAME,
                       "SAMPLE_ID",
                       where=("SAMPLE_PROJECT",
                              "AndrewBloggs"))
        for line in SampleSheet(sample_sheet_file):
            if line['Lane'] in (1,2):
                self.assertEqual(line["Sample_Name"],
                                 line["Sample_ID"])

class TestSampleSheetSetSampleID(unittest.TestCase):
    """
    Test the 'set_sample_id' functionality of the 'samplesheet' command
    """
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestSampleSheetSetProject')
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
            os.chmod(os.path.dirname(name),0o755)
            os.chmod(name,0o655)
            os.remove(name)
        shutil.rmtree(self.dirn,onerror=del_rw)

    def test_samplesheet_set_sample_id(self):
        """
        samplesheet: set sample IDs
        """
        # Make an empty mock auto-process directory
        mockdir = MockAnalysisDir(
            '221219_SN00879_0087_000000000-AGEW9',
            'hiseq',
            no_undetermined=True,
            top_dir=self.dirn)
        mockdir.create()
        # Update content for sample sheet
        sample_sheet_data = u"""[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,Sample_AB1,AB1,,,N701,CGATGTAT,N501,TCTTTCCC,AndrewBloggs,
2,Sample_AB2,AB2,,,N702,TGACCAAT,N502,TCTTTCCC,AndrewBloggs,
3,Sample_CD1,CD1,,,N701,CGATGTAT,N501,TCTTTCCC,CarlDewey,
4,Sample_FG2,FG2,,,N702,TGACCAAT,N502,TCTTTCCC,FilipeGreer,
"""
        sample_sheet_file = os.path.join(mockdir.dirn,
                                         "custom_SampleSheet.csv")
        with open(sample_sheet_file,'wt') as fp:
            fp.write(sample_sheet_data)
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Update sample sheet name for lane 3
        ap.samplesheet(SampleSheetOperation.SET_SAMPLE_ID,
                       "SAMPLE_NAME",
                       lanes=[3])
        for line in SampleSheet(sample_sheet_file):
            if line['Lane'] == 3:
                self.assertEqual(line["Sample_ID"],"CD1")
        # Update sample sheet name matching on project name
        ap.samplesheet(SampleSheetOperation.SET_SAMPLE_ID,
                       "SAMPLE_NAME",
                       where=("SAMPLE_PROJECT",
                              "AndrewBloggs"))
        for line in SampleSheet(sample_sheet_file):
            if line['Lane'] in (1,2):
                self.assertEqual(line["Sample_Name"],
                                 line["Sample_ID"])

class TestSampleSheetImport(unittest.TestCase):
    """
    Test the 'import_samplesheet' functionality of the 'samplesheet' command
    """
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestSampleSheetImport')
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
            os.chmod(os.path.dirname(name),0o755)
            os.chmod(name,0o655)
            os.remove(name)
        shutil.rmtree(self.dirn,onerror=del_rw)

    def test_import_samplesheet(self):
        """
        samplesheet: import sample sheet data
        """
        # Make an empty mock auto-process directory
        mockdir = MockAnalysisDir(
            '221219_SN00879_0087_000000000-AGEW9',
            'hiseq',
            no_undetermined=True,
            top_dir=self.dirn)
        mockdir.create()
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Check initial status
        self.assertTrue(os.path.exists(
            os.path.join(mockdir.dirn,
                         "SampleSheet.orig.csv")))
        self.assertTrue(os.path.exists(
            os.path.join(mockdir.dirn,
                         "custom_SampleSheet.csv")))
        self.assertEqual(ap.params.sample_sheet,
                         os.path.join(mockdir.dirn,"custom_SampleSheet.csv"))
        self.assertFalse(os.path.exists(
            os.path.join(mockdir.dirn,
                         "SampleSheet.imported.csv")))
        self.assertFalse(os.path.exists(
            os.path.join(mockdir.dirn,
                         "custom_SampleSheet.imported.csv")))
        # Make new sample sheet to import
        sample_sheet_data = u"""[Header]

[Reads]

[Settings]

[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,AB1,Sample_AB1,,,N701,CGATGTAT,N501,TCTTTCCC,AndrewBloggs,
2,AB2,Sample_AB2,,,N702,TGACCAAT,N502,TCTTTCCC,AndrewBloggs,
3,CD1,Sample_CD1,,,N701,CGATGTAT,N501,TCTTTCCC,CarlDewey,
4,FG2,Sample_FG2,,,N702,TGACCAAT,N502,TCTTTCCC,FilipeGreer,
"""
        new_sample_sheet = os.path.join(self.dirn,"new_sample_sheet")
        with open(new_sample_sheet,'wt') as fp:
            fp.write(sample_sheet_data)
        # Import the new sample sheet
        import_samplesheet(ap,new_sample_sheet)
        # Check expected state after import
        self.assertTrue(os.path.exists(
            os.path.join(mockdir.dirn,
                         "SampleSheet.imported.csv")))
        self.assertTrue(os.path.exists(
            os.path.join(mockdir.dirn,
                         "custom_SampleSheet.imported.csv")))
        self.assertEqual(ap.params.sample_sheet,
                         os.path.join(mockdir.dirn,
                                      "custom_SampleSheet.imported.csv"))
        with open(ap.params.sample_sheet,'rt') as fp:
            self.assertEqual(fp.read(),
                             sample_sheet_data)
