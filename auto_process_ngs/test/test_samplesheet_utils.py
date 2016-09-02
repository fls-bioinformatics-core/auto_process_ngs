#######################################################################
# Tests for samplesheet_utils.py module
#######################################################################

import os
import unittest
import tempfile
import shutil
import codecs
import cStringIO
from bcftbx.IlluminaData import SampleSheet
from auto_process_ngs.samplesheet_utils import SampleSheetLinter
from auto_process_ngs.samplesheet_utils import has_invalid_characters
from auto_process_ngs.samplesheet_utils import get_close_names

sample_sheet_header = """[Header]
IEMFileVersion,4
Date,4/11/2014
Workflow,Metagenomics
Application,Metagenomics 16S rRNA
Assay,Nextera XT
Description,
Chemistry,Amplicon

[Reads]
150
150

[Settings]
Adapter,CTGTCTCTTATACACATCT

"""

class TestCloseProjectNames(unittest.TestCase):
    def test_sample_sheet_without_close_project_names(self):
        self.sample_sheet_no_close_project_names = """[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,AB1,AB1,,,N701,CGATGTAT,N501,TCTTTCCC,Andrew_Bloggs,
1,AB2,AB2,,,N702,TGACCAAT,N502,TCTTTCCC,Andrew_Bloggs,
2,CD1,CD1,,,N701,CGATGTAT,N501,TCTTTCCC,Carl_Dewey,
2,FG2,FG2,,,N702,TGACCAAT,N502,TCTTTCCC,Filipe_Greer,
"""
        linter = SampleSheetLinter(fp=cStringIO.StringIO(
            self.sample_sheet_no_close_project_names))
        self.assertEqual(linter.close_project_names(),{})
    def test_sample_sheet_with_close_project_names(self):
        self.sample_sheet_close_project_names = """[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,AB1,AB1,,,N701,CGATGTAT,N501,TCTTTCCC,Andrew_Blogs,
1,AB2,AB2,,,N702,TGACCAAT,N502,TCTTTCCC,Andrew_Bloggs,
2,CD1,CD1,,,N701,CGATGTAT,N501,TCTTTCCC,Carl_Dewey,
2,FG2,FG2,,,N702,TGACCAAT,N502,TCTTTCCC,Filipe_Greer,
"""
        linter = SampleSheetLinter(fp=cStringIO.StringIO(
            self.sample_sheet_close_project_names))
        self.assertEqual(linter.close_project_names(),
                         { 'Andrew_Bloggs': ['Andrew_Blogs'],
                           'Andrew_Blogs': ['Andrew_Bloggs'] })

class TestSamplesInMultipleProjects(unittest.TestCase):
    def test_sample_sheet_without_samples_in_multiple_projects(self):
        self.sample_sheet_no_samples_in_multiple_projects = """[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,AB1,AB1,,,N701,CGATGTAT,N501,TCTTTCCC,Andrew_Bloggs,
1,AB2,AB2,,,N702,TGACCAAT,N502,TCTTTCCC,Andrew_Bloggs,
2,AB1,AB1,,,N701,CGATGTAT,N501,TCTTTCCC,Andrew_Bloggs,
2,CD1,CD1,,,N701,CGATGTAT,N501,TCTTTCCC,Carl_Dewey,
2,FG2,FG2,,,N702,TGACCAAT,N502,TCTTTCCC,Filipe_Greer,
"""
        linter = SampleSheetLinter(fp=cStringIO.StringIO(
            self.sample_sheet_no_samples_in_multiple_projects))
        self.assertEqual(linter.samples_in_multiple_projects(),{})
    def test_sample_sheet_with_samples_in_multiple_projects(self):
        self.sample_sheet_samples_in_multiple_projects = """[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,AB1,AB1,,,N701,CGATGTAT,N501,TCTTTCCC,Andrew_Bloggs1,
1,AB2,AB2,,,N702,TGACCAAT,N502,TCTTTCCC,Andrew_Bloggs1,
2,AB1,AB1,,,N701,CGATGTAT,N501,TCTTTCCC,Andrew_Bloggs2,
2,CD1,CD1,,,N701,CGATGTAT,N501,TCTTTCCC,Carl_Dewey,
2,FG2,FG2,,,N702,TGACCAAT,N502,TCTTTCCC,Filipe_Greer,
"""
        linter = SampleSheetLinter(fp=cStringIO.StringIO(
            self.sample_sheet_samples_in_multiple_projects))
        self.assertEqual(linter.samples_in_multiple_projects(),
                         { 'AB1': ['Andrew_Bloggs1','Andrew_Bloggs2'] })

class TestSamplesWithMultipleBarcodes(unittest.TestCase):
    def test_sample_sheet_without_multiple_barcodes(self):
        self.sample_sheet_without_multiple_barcodes = """[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,AB1,AB1,,,N701,CGATGTAT,N501,TCTTTCCC,Andrew_Bloggs,
1,AB2,AB2,,,N702,TGACCAAT,N502,TCTTTCCC,Andrew_Bloggs,
2,AB1,AB1,,,N701,CGATGTAT,N501,TCTTTCCC,Andrew_Bloggs,
2,CD1,CD1,,,N701,CGATGTAT,N501,TCTTTCCC,Carl_Dewey,
2,FG2,FG2,,,N702,TGACCAAT,N502,TCTTTCCC,Filipe_Greer,
"""
        linter = SampleSheetLinter(fp=cStringIO.StringIO(
            self.sample_sheet_without_multiple_barcodes))
        self.assertEqual(linter.samples_with_multiple_barcodes(),{})
    def test_sample_sheet_with_multiple_barcodes(self):
        self.sample_sheet_with_multiple_barcodes = """[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,AB1,AB1,,,N701,CGATGTAT,N501,TCTTTCCC,Andrew_Bloggs,
1,AB2,AB2,,,N702,TGACCAAT,N502,TCTTTCCC,Andrew_Bloggs,
2,AB1,AB1,,,N701,CTGTAGTA,N501,TCTTTCCC,Andrew_Bloggs,
2,CD1,CD1,,,N701,CGATGTAT,N501,TCTTTCCC,Carl_Dewey,
2,FG2,FG2,,,N702,TGACCAAT,N502,TCTTTCCC,Filipe_Greer,
"""
        linter = SampleSheetLinter(fp=cStringIO.StringIO(
            self.sample_sheet_with_multiple_barcodes))
        self.assertEqual(linter.samples_with_multiple_barcodes(),
                         { 'AB1': ['CGATGTAT-TCTTTCCC','CTGTAGTA-TCTTTCCC'] })

class TestLinterHasInvalidLines(unittest.TestCase):
    def test_sample_sheet_without_invalid_lines(self):
        self.sample_sheet_without_invalid_lines = """[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,AB1,AB1,,,N701,CGATGTAT,N501,TCTTTCCC,Andrew_Bloggs,
1,AB2,AB2,,,N702,TGACCAAT,N502,TCTTTCCC,Andrew_Bloggs,
2,CD1,CD1,,,N701,CGATGTAT,N501,TCTTTCCC,Carl_Dewey,
2,FG2,FG2,,,N702,TGACCAAT,N502,TCTTTCCC,Filipe_Greer,
"""
        linter = SampleSheetLinter(fp=cStringIO.StringIO(
            self.sample_sheet_without_invalid_lines))
        self.assertFalse(linter.has_invalid_lines())
    def test_sample_sheet_with_invalid_lines_missing_lane(self):
        self.sample_sheet_with_invalid_lines = """[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,AB1,AB1,,,N701,CGATGTAT,N501,TCTTTCCC,Andrew_Bloggs,
1,AB2,AB2,,,N702,TGACCAAT,N502,TCTTTCCC,Andrew_Bloggs,
2,CD1,CD1,,,N701,CGATGTAT,N501,TCTTTCCC,Carl_Dewey,
2,FG2,FG2,,,N702,TGACCAAT,N502,TCTTTCCC,Filipe_Greer,
,,,,,,,,,Filipe_Greer,
,,,,,,,,,Filipe_Greer,
"""
        linter = SampleSheetLinter(fp=cStringIO.StringIO(
            self.sample_sheet_with_invalid_lines))
        self.assertTrue(linter.has_invalid_lines())
    def test_sample_sheet_with_invalid_lines_no_lane(self):
        self.sample_sheet_with_invalid_lines_no_lane = """[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
AB1,AB1,,,N701,CGATGTAT,N501,TCTTTCCC,Andrew_Bloggs,
AB2,AB2,,,N702,TGACCAAT,N502,TCTTTCCC,Andrew_Bloggs,
CD1,CD1,,,N701,CGATGTAT,N501,TCTTTCCC,Carl_Dewey,
FG2,FG2,,,N702,TGACCAAT,N502,TCTTTCCC,Filipe_Greer,
,,,,,,,,Filipe_Greer,
,,,,,,,,Filipe_Greer,
"""
        linter = SampleSheetLinter(fp=cStringIO.StringIO(
            self.sample_sheet_with_invalid_lines_no_lane))
        self.assertTrue(linter.has_invalid_lines())

class TestLinterHasInvalidCharacters(unittest.TestCase):
    def test_sample_sheet_with_non_printing_ascii_character(self):
        self.sample_sheet_with_non_printing_ascii_characters = """[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,AB1,AB1,,,N701,CGATGTAT,N501,TCTTTCCC,Andrew_Bloggs,\x13
1,AB2,AB2,,,N702,TGACCAAT,N502,TCTTTCCC,Andrew_Bloggs,
2,CD1,CD1,,,N701,CGATGTAT,N501,TCTTTCCC,Carl_Dewey,
2,FG2,FG2,,,N702,TGACCAAT,N502,TCTTTCCC,Filipe_Greer,
"""
        linter = SampleSheetLinter(fp=cStringIO.StringIO(
            self.sample_sheet_with_non_printing_ascii_characters))
        self.assertTrue(linter.has_invalid_characters())
    def test_sample_sheet_with_valid_characters(self):
        self.sample_sheet_with_valid_characters = """[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,AB1,AB1,,,N701,CGATGTAT,N501,TCTTTCCC,Andrew_Bloggs,
1,AB2,AB2,,,N702,TGACCAAT,N502,TCTTTCCC,Andrew_Bloggs,
2,CD1,CD1,,,N701,CGATGTAT,N501,TCTTTCCC,Carl_Dewey,
2,FG2,FG2,,,N702,TGACCAAT,N502,TCTTTCCC,Filipe_Greer,
"""
        linter = SampleSheetLinter(fp=cStringIO.StringIO(
            self.sample_sheet_with_valid_characters))
        self.assertFalse(linter.has_invalid_characters())

class TestHasInvalidCharactersFunction(unittest.TestCase):
    def setUp(self):
        self.wd = tempfile.mkdtemp()
    def tearDown(self):
        if os.path.isdir(self.wd):
            shutil.rmtree(self.wd)
    def test_file_with_non_printing_ascii_character(self):
        """has_invalid_characters: detects non-printing ASCII character in file
        """
        non_printing_ascii_file = os.path.join(self.wd,
                                               "test.nonprintingascii")
        with open(non_printing_ascii_file,'wb') as fp:
            fp.write(u"""This file contains:
a non-printing ASCII ctrl-S character here\x13
""")
        self.assertTrue(has_invalid_characters(non_printing_ascii_file))
    def test_text_with_non_printing_ascii_character(self):
        """has_invalid_characters: detects non-printing ASCII character in text
        """
        self.assertTrue(has_invalid_characters(text=u"""This file contains:
a non-printing ASCII ctrl-S character here\x13
"""))
    def test_file_with_non_ascii_character(self):
        """has_invalid_characters: detects non-ASCII character in file
        """
        non_ascii_file = os.path.join(self.wd,"test.nonascii")
        with codecs.open(non_ascii_file,'wb',encoding='utf-8') as fp:
            fp.write(u"""This file contains:
a non-ASCII character here\x80
""")
        self.assertTrue(has_invalid_characters(non_ascii_file))
    def test_text_with_non_ascii_character(self):
        """has_invalid_characters: detects non-ASCII character in text
        """
        self.assertTrue(has_invalid_characters(text=u"""This text contains:
a non-ASCII character here\x80
"""))
    def test_file_with_valid_characters(self):
        """has_invalid_characters: works for valid file
        """
        valid_file = os.path.join(self.wd,"test.valid")
        with open(valid_file,'wb') as fp:
            fp.write(u"""This file contains valid characters:
- ABCDEFGHIJKLMNOPQRSTUVWXYZ
- abcdefghijklmnopqrstuvwxyz
- 01234567890
- !"$%^&*()_-+=\|/:;'@#~`?
- {}[]
- \t\n
""")
        self.assertFalse(has_invalid_characters(valid_file))
    def test_text_with_valid_characters(self):
        """has_invalid_characters: works for valid text
        """
        self.assertFalse(has_invalid_characters(text=u"""This file contains valid characters:
- ABCDEFGHIJKLMNOPQRSTUVWXYZ
- abcdefghijklmnopqrstuvwxyz
- 01234567890
- !"$%^&*()_-+=\|/:;'@#~`?
- {}[]
- \t\n
"""))

class TestCloseNamesFunction(unittest.TestCase):
    def test_close_names(self):
        """get_close_names: check function returns expected results
        """
        self.assertEqual(get_close_names(("Andrew Bloggs",
                                          "Carl Dewey",
                                          "Filipe Greer",)),{})
        self.assertEqual(get_close_names(("Andrew Blogs",
                                          "Andrew Bloggs",
                                          "Carl Dewey",
                                          "Filipe Greer",)),
                         { "Andrew Bloggs": ["Andrew Blogs"],
                           "Andrew Blogs": ["Andrew Bloggs"] })
        
