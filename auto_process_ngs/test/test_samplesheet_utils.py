#######################################################################
# Tests for samplesheet_utils.py module
#######################################################################

import os
import unittest
import tempfile
import shutil
import codecs
from io import StringIO
from bcftbx.IlluminaData import SampleSheet
from auto_process_ngs.samplesheet_utils import SampleSheetLinter
from auto_process_ngs.samplesheet_utils import has_invalid_characters
from auto_process_ngs.samplesheet_utils import barcode_is_valid
from auto_process_ngs.samplesheet_utils import barcode_is_10xgenomics
from auto_process_ngs.samplesheet_utils import summarise_outputs
from auto_process_ngs.samplesheet_utils import get_close_names
from auto_process_ngs.samplesheet_utils import set_samplesheet_column

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
        self.sample_sheet_no_close_project_names = u"""[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,AB1,AB1,,,N701,CGATGTAT,N501,TCTTTCCC,Andrew_Bloggs,
1,AB2,AB2,,,N702,TGACCAAT,N502,TCTTTCCC,Andrew_Bloggs,
2,CD1,CD1,,,N701,CGATGTAT,N501,TCTTTCCC,Carl_Dewey,
2,FG2,FG2,,,N702,TGACCAAT,N502,TCTTTCCC,Filipe_Greer,
"""
        linter = SampleSheetLinter(fp=StringIO(
            self.sample_sheet_no_close_project_names))
        self.assertEqual(linter.close_project_names(),{})
    def test_sample_sheet_with_close_project_names(self):
        self.sample_sheet_close_project_names = u"""[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,AB1,AB1,,,N701,CGATGTAT,N501,TCTTTCCC,Andrew_Blogs,
1,AB2,AB2,,,N702,TGACCAAT,N502,TCTTTCCC,Andrew_Bloggs,
2,CD1,CD1,,,N701,CGATGTAT,N501,TCTTTCCC,Carl_Dewey,
2,FG2,FG2,,,N702,TGACCAAT,N502,TCTTTCCC,Filipe_Greer,
"""
        linter = SampleSheetLinter(fp=StringIO(
            self.sample_sheet_close_project_names))
        self.assertEqual(linter.close_project_names(),
                         { 'Andrew_Bloggs': ['Andrew_Blogs'],
                           'Andrew_Blogs': ['Andrew_Bloggs'] })

class TestSamplesInMultipleProjects(unittest.TestCase):
    def test_sample_sheet_without_samples_in_multiple_projects(self):
        self.sample_sheet_no_samples_in_multiple_projects = u"""[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,AB1,AB1,,,N701,CGATGTAT,N501,TCTTTCCC,Andrew_Bloggs,
1,AB2,AB2,,,N702,TGACCAAT,N502,TCTTTCCC,Andrew_Bloggs,
2,AB1,AB1,,,N701,CGATGTAT,N501,TCTTTCCC,Andrew_Bloggs,
2,CD1,CD1,,,N701,CGATGTAT,N501,TCTTTCCC,Carl_Dewey,
2,FG2,FG2,,,N702,TGACCAAT,N502,TCTTTCCC,Filipe_Greer,
"""
        linter = SampleSheetLinter(fp=StringIO(
            self.sample_sheet_no_samples_in_multiple_projects))
        self.assertEqual(linter.samples_in_multiple_projects(),{})
    def test_sample_sheet_with_samples_in_multiple_projects(self):
        self.sample_sheet_samples_in_multiple_projects = u"""[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,AB1,AB1,,,N701,CGATGTAT,N501,TCTTTCCC,Andrew_Bloggs1,
1,AB2,AB2,,,N702,TGACCAAT,N502,TCTTTCCC,Andrew_Bloggs1,
2,AB1,AB1,,,N701,CGATGTAT,N501,TCTTTCCC,Andrew_Bloggs2,
2,CD1,CD1,,,N701,CGATGTAT,N501,TCTTTCCC,Carl_Dewey,
2,FG2,FG2,,,N702,TGACCAAT,N502,TCTTTCCC,Filipe_Greer,
"""
        linter = SampleSheetLinter(fp=StringIO(
            self.sample_sheet_samples_in_multiple_projects))
        self.assertEqual(linter.samples_in_multiple_projects(),
                         { 'AB1': ['Andrew_Bloggs1','Andrew_Bloggs2'] })

class TestSamplesWithMultipleBarcodes(unittest.TestCase):
    def test_sample_sheet_without_multiple_barcodes(self):
        self.sample_sheet_without_multiple_barcodes = u"""[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,AB1,AB1,,,N701,CGATGTAT,N501,TCTTTCCC,Andrew_Bloggs,
1,AB2,AB2,,,N702,TGACCAAT,N502,TCTTTCCC,Andrew_Bloggs,
2,AB1,AB1,,,N701,CGATGTAT,N501,TCTTTCCC,Andrew_Bloggs,
2,CD1,CD1,,,N701,CGATGTAT,N501,TCTTTCCC,Carl_Dewey,
2,FG2,FG2,,,N702,TGACCAAT,N502,TCTTTCCC,Filipe_Greer,
"""
        linter = SampleSheetLinter(fp=StringIO(
            self.sample_sheet_without_multiple_barcodes))
        self.assertEqual(linter.samples_with_multiple_barcodes(),{})
    def test_sample_sheet_with_multiple_barcodes(self):
        self.sample_sheet_with_multiple_barcodes = u"""[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,AB1,AB1,,,N701,CGATGTAT,N501,TCTTTCCC,Andrew_Bloggs,
1,AB2,AB2,,,N702,TGACCAAT,N502,TCTTTCCC,Andrew_Bloggs,
2,AB1,AB1,,,N701,CTGTAGTA,N501,TCTTTCCC,Andrew_Bloggs,
2,CD1,CD1,,,N701,CGATGTAT,N501,TCTTTCCC,Carl_Dewey,
2,FG2,FG2,,,N702,TGACCAAT,N502,TCTTTCCC,Filipe_Greer,
"""
        linter = SampleSheetLinter(fp=StringIO(
            self.sample_sheet_with_multiple_barcodes))
        self.assertEqual(linter.samples_with_multiple_barcodes(),
                         { 'AB1': ['CGATGTAT-TCTTTCCC','CTGTAGTA-TCTTTCCC'] })

class TestLinterHasInvalidLines(unittest.TestCase):
    def test_sample_sheet_without_invalid_lines(self):
        self.sample_sheet_without_invalid_lines = u"""[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,AB1,AB1,,,N701,CGATGTAT,N501,TCTTTCCC,Andrew_Bloggs,
1,AB2,AB2,,,N702,TGACCAAT,N502,TCTTTCCC,Andrew_Bloggs,
2,CD1,CD1,,,N701,CGATGTAT,N501,TCTTTCCC,Carl_Dewey,
2,FG2,FG2,,,N702,TGACCAAT,N502,TCTTTCCC,Filipe_Greer,
"""
        linter = SampleSheetLinter(fp=StringIO(
            self.sample_sheet_without_invalid_lines))
        self.assertFalse(linter.has_invalid_lines())
    def test_sample_sheet_with_invalid_lines_missing_lane(self):
        self.sample_sheet_with_invalid_lines = u"""[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,AB1,AB1,,,N701,CGATGTAT,N501,TCTTTCCC,Andrew_Bloggs,
1,AB2,AB2,,,N702,TGACCAAT,N502,TCTTTCCC,Andrew_Bloggs,
2,CD1,CD1,,,N701,CGATGTAT,N501,TCTTTCCC,Carl_Dewey,
2,FG2,FG2,,,N702,TGACCAAT,N502,TCTTTCCC,Filipe_Greer,
,,,,,,,,,Filipe_Greer,
,,,,,,,,,Filipe_Greer,
"""
        linter = SampleSheetLinter(fp=StringIO(
            self.sample_sheet_with_invalid_lines))
        self.assertTrue(linter.has_invalid_lines())
    def test_sample_sheet_with_invalid_lines_no_lane(self):
        self.sample_sheet_with_invalid_lines_no_lane = u"""[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
AB1,AB1,,,N701,CGATGTAT,N501,TCTTTCCC,Andrew_Bloggs,
AB2,AB2,,,N702,TGACCAAT,N502,TCTTTCCC,Andrew_Bloggs,
CD1,CD1,,,N701,CGATGTAT,N501,TCTTTCCC,Carl_Dewey,
FG2,FG2,,,N702,TGACCAAT,N502,TCTTTCCC,Filipe_Greer,
,,,,,,,,Filipe_Greer,
,,,,,,,,Filipe_Greer,
"""
        linter = SampleSheetLinter(fp=StringIO(
            self.sample_sheet_with_invalid_lines_no_lane))
        self.assertTrue(linter.has_invalid_lines())

class TestLinterHasInvalidBarcodes(unittest.TestCase):
    def test_sample_sheet_with_valid_barcodes(self):
        self.sample_sheet_with_valid_barcodes = u"""[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,AB1,AB1,,,N701,CGATGTAT,N501,TCTTTCCC,Andrew_Bloggs,
1,AB2,AB2,,,N702,TGACCAAT,N502,TCTTTCCC,Andrew_Bloggs,
2,CD1,CD1,,,N701,CGATGTAT,N501,TCTTTCCC,Carl_Dewey,
2,FG2,FG2,,,N702,TGACCAAT,N502,TCTTTCCC,Filipe_Greer,
"""
        linter = SampleSheetLinter(fp=StringIO(
            self.sample_sheet_with_valid_barcodes))
        self.assertEqual(linter.has_invalid_barcodes(),list())
    def test_sample_sheet_with_no_barcodes(self):
        self.sample_sheet_with_no_barcodes = u"""[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,AB1,AB1,,,,,,,Andrew_Bloggs,
"""
        linter = SampleSheetLinter(fp=StringIO(
            self.sample_sheet_with_no_barcodes))
        self.assertEqual(linter.has_invalid_barcodes(),list())
    def test_sample_sheet_with_10xgenomics_barcodes(self):
        self.sample_sheet_with_10xgenomics_barcodes = u"""[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description
1,AB1,AB1,,,N701,SI-GA-B3,Philip_Crook,
2,AB2,AB2,,,N701,SI-GA-B3,Philip_Crook,
3,AB3,AB3,,,N701,SI-GA-G1,Philip_Crook,
4,AB4,AB4,,,N701,SI-GA-G1,Philip_Crook,
5,AB5,AB5,,,N701,SI-GA-H1,Philip_Crook,
6,AB6,AB6,,,N701,SI-GA-H1,Philip_Crook,
7,AB7,AB7,,,N701,SI-P03-C9,Philip_Crook,
8,AB8,AB8,,,N701,SI-P03-C9,Philip_Crook,
"""
        linter = SampleSheetLinter(fp=StringIO(
            self.sample_sheet_with_10xgenomics_barcodes))
        self.assertEqual(linter.has_invalid_barcodes(),list())
    def test_sample_sheet_with_invalid_barcodes(self):
        self.sample_sheet_with_invalid_barcodes = u"""[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,AB1,AB1,,,N701,CGATGTNN,N501,TCTTTCCC,Andrew_Bloggs,
1,AB2,AB2,,,N702,TGACCAAT,N502,TCTTNNNN,Andrew_Bloggs,
2,CD1,CD1,,,N701,CGATGTAT,N501,TCTTTCCC,Carl_Dewey,
2,FG2,FG2,,,N702,tgaccaat,N502,tctttccc,Filipe_Greer,
3,HJ1,HJ1,,,N702,AGHJATTC,N502,TCTTTCCC,Hugh_Jaxman,
3,HJ2,HJ2,,,N702,CGATGTAT,N502,TCTTKCCC,Hugh_Jaxman,
4,KL1,KL1,,,N702,TGACC&^&,N502,TCTTTCCC,Karl_Landseer,
4,KL2,KL2,,,N702,TGACCAAT,N502,TC%$TCCC,Karl_Landseer,
"""
        linter = SampleSheetLinter(fp=StringIO(
            self.sample_sheet_with_invalid_barcodes))
        self.assertTrue(linter.has_invalid_barcodes())

class TestLinterHasInvalidCharacters(unittest.TestCase):
    def test_sample_sheet_with_non_printing_ascii_character(self):
        self.sample_sheet_with_non_printing_ascii_characters = u"""[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,AB1,AB1,,,N701,CGATGTAT,N501,TCTTTCCC,Andrew_Bloggs,\x13
1,AB2,AB2,,,N702,TGACCAAT,N502,TCTTTCCC,Andrew_Bloggs,
2,CD1,CD1,,,N701,CGATGTAT,N501,TCTTTCCC,Carl_Dewey,
2,FG2,FG2,,,N702,TGACCAAT,N502,TCTTTCCC,Filipe_Greer,
"""
        linter = SampleSheetLinter(fp=StringIO(
            self.sample_sheet_with_non_printing_ascii_characters))
        self.assertTrue(linter.has_invalid_characters())
    def test_sample_sheet_with_valid_characters(self):
        self.sample_sheet_with_valid_characters = u"""[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,AB1,AB1,,,N701,CGATGTAT,N501,TCTTTCCC,Andrew_Bloggs,
1,AB2,AB2,,,N702,TGACCAAT,N502,TCTTTCCC,Andrew_Bloggs,
2,CD1,CD1,,,N701,CGATGTAT,N501,TCTTTCCC,Carl_Dewey,
2,FG2,FG2,,,N702,TGACCAAT,N502,TCTTTCCC,Filipe_Greer,
"""
        linter = SampleSheetLinter(fp=StringIO(
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
        with open(non_printing_ascii_file,'wt') as fp:
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
        with open(valid_file,'wt') as fp:
            fp.write(u"""This file contains valid characters:
- ABCDEFGHIJKLMNOPQRSTUVWXYZ
- abcdefghijklmnopqrstuvwxyz
- 01234567890
- !"$%^&*()_-+=\\|/:;'@#~`?
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
- !"$%^&*()_-+=\\|/:;'@#~`?
- {}[]
- \t\n
"""))

class TestBarcodeIs10xGenomics(unittest.TestCase):
    def test_barcode_is_10xgenomics(self):
        """barcode_is_10xgenomics: identifies 10xGenomics 'barcodes'
        """
        self.assertTrue(barcode_is_10xgenomics("SI-GA-B3"))
        self.assertTrue(barcode_is_10xgenomics("SI-GA-G1"))
        self.assertTrue(barcode_is_10xgenomics("SI-GA-H1"))
        self.assertTrue(barcode_is_10xgenomics("SI-P03-C9"))
    def test_barcode_is_not_10xgenomics(self):
        """barcode_is_10xgenomics: identifies non-10xGenomics barcodes
        """
        self.assertFalse(barcode_is_10xgenomics("TGACCAAT"))

class TestBarcodeIsValidFunction(unittest.TestCase):
    def test_standard_barcodes(self):
        """barcode_is_valid: standard barcode
        """
        self.assertTrue(barcode_is_valid("TGACCAAT"))
    def test_empty_barcode(self):
        """barcode_is_valid: 'empty' barcode
        """
        self.assertTrue(barcode_is_valid(""))
    def test_barcode_with_Ns(self):
        """barcode_is_valid: barcode with Ns
        """
        self.assertFalse(barcode_is_valid("TGACCANN"))
    def test_barcode_with_lower_case_characters(self):
        """barcode_is_valid: barcode with lower case characters
        """
        self.assertFalse(barcode_is_valid("tgaccaat"))
    def test_barcode_with_random_string(self):
        """barcode_is_valid: barcode with random string
        """
        self.assertFalse(barcode_is_valid("TADF%$12"))
    def test_10xgenomics_barcodes(self):
        """barcode_is_valid: 10xGenomics 'barcodes'
        """
        self.assertTrue(barcode_is_valid("SI-GA-B3"))
        self.assertTrue(barcode_is_valid("SI-GA-G1"))
        self.assertTrue(barcode_is_valid("SI-GA-H1"))
        self.assertTrue(barcode_is_valid("SI-P03-C9"))

class TestSummariseOutputs(unittest.TestCase):
    """
    Tests for the 'summarise_outputs' function
    """
    def setUp(self):
        self.wd = tempfile.mkdtemp()

    def tearDown(self):
        if os.path.isdir(self.wd):
            shutil.rmtree(self.wd)

    def write_sample_sheet(self,content):
        self.sample_sheet_file = os.path.join(self.wd,
                                              "SampleSheet.csv")
        with open(self.sample_sheet_file,'wt') as fp:
            fp.write(content)
        return self.sample_sheet_file

    def test_summarise_outputs_single_index_no_lanes(self):
        """
        summarise_outputs: sample sheet with single index, no lanes
        """
        self.write_sample_sheet("""[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description
SM1,SM1,,,A006,GCCAAT,SerenaMcCauley,
SM2,SM2,,,A012,CTTGTA,SerenaMcCauley,
""")
        self.assertEqual(summarise_outputs(
            sample_sheet_file=self.sample_sheet_file),"""Project: 'SerenaMcCauley':
- L1	SM1-2	2 samples	6bp indexes""")

    def test_summarise_outputs_single_index_multiple_lanes(self):
        """
        summarise_outputs: sample sheet with single index, multiple lanes
        """
        self.write_sample_sheet("""[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description
1,SM1,SM1,,,A006,GCCAAT,SerenaMcCauley,
1,SM2,SM2,,,A012,CTTGTA,SerenaMcCauley,
2,CD1,CD1,,,N701,CGATGT,CarlDewey,
2,CD2,CD2,,,N702,TCTTTC,CarlDewey,
""")
        self.assertEqual(summarise_outputs(
            sample_sheet_file=self.sample_sheet_file),"""Project: 'CarlDewey':
- L2	CD1-2	2 samples	6bp indexes
Project: 'SerenaMcCauley':
- L1	SM1-2	2 samples	6bp indexes""")

    def test_summarise_outputs_10x_barcodes(self):
        """
        summarise_outputs: sample sheet with 10x Genomics barcodes
        """
        self.write_sample_sheet("""[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description
SL1,SL1,,,N701,SI-GA-A2,SarahLee,
SL2,SL2,,,N702,SI-GA-B2,SarahLee,
SL3,SL3,,,N703,SI-GA-C2,SarahLee,
SL4,SL4,,,N704,SI-GA-D2,SarahLee,
""")
        self.assertEqual(summarise_outputs(
            sample_sheet_file=self.sample_sheet_file),"""Project: 'SarahLee':
- L1	SL1-4	4 samples	10x Genomics indexes""")

    def test_summarise_outputs_dual_index(self):
        """
        summarise_outputs: sample sheet with dual index barcodes
        """
        self.write_sample_sheet("""[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,AB1,AB1,,,N701,CGATGTAT,N501,TCTTTCCC,AndrewBloggs,
2,CD1,CD1,,,N701,CGATGTAT,N501,TCTTTCCC,CarlDewey,
""")
        self.assertEqual(summarise_outputs(
            sample_sheet_file=self.sample_sheet_file),"""Project: 'AndrewBloggs':
- L1	AB1	1 sample	8x8bp indexes
Project: 'CarlDewey':
- L2	CD1	1 sample	8x8bp indexes""")

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

class TestSetSampleSheetSetColumn(unittest.TestCase):
    def test_set_samplesheet_column_project_to_name(self):
        """
        set_samplesheet_column: set project to sample name
        """
        sample_sheet_data = u"""[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,AB1,Sample_AB1,,,N701,CGATGTAT,N501,TCTTTCCC,Andrew_Bloggs,
1,AB2,Sample_AB2,,,N702,TGACCAAT,N502,TCTTTCCC,Andrew_Bloggs,
2,CD1,Sample_CD1,,,N701,CGATGTAT,N501,TCTTTCCC,Carl_Dewey,
2,FG2,Sample_FG2,,,N702,TGACCAAT,N502,TCTTTCCC,Filipe_Greer,
"""
        s = SampleSheet(fp=StringIO(sample_sheet_data))
        s1 = set_samplesheet_column(s,"SAMPLE_PROJECT","SAMPLE_NAME")
        for line in s1:
            self.assertEqual(line['Sample_Project'],line['Sample_Name'])
    def test_set_samplesheet_column_project_to_id(self):
        """
        set_samplesheet_column: set project to sample ID
        """
        sample_sheet_data = u"""[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,AB1,Sample_AB1,,,N701,CGATGTAT,N501,TCTTTCCC,Andrew_Bloggs,
1,AB2,Sample_AB2,,,N702,TGACCAAT,N502,TCTTTCCC,Andrew_Bloggs,
2,CD1,Sample_CD1,,,N701,CGATGTAT,N501,TCTTTCCC,Carl_Dewey,
2,FG2,Sample_FG2,,,N702,TGACCAAT,N502,TCTTTCCC,Filipe_Greer,
"""
        s = SampleSheet(fp=StringIO(sample_sheet_data))
        s1 = set_samplesheet_column(s,"SAMPLE_PROJECT","SAMPLE_ID")
        for line in s1:
            self.assertEqual(line['Sample_Project'],line['Sample_ID'])
    def test_set_samplesheet_column_project_to_value(self):
        """
        set_samplesheet_column: set project to a specific value
        """
        sample_sheet_data = u"""[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,AB1,Sample_AB1,,,N701,CGATGTAT,N501,TCTTTCCC,Andrew_Bloggs,
1,AB2,Sample_AB2,,,N702,TGACCAAT,N502,TCTTTCCC,Andrew_Bloggs,
2,CD1,Sample_CD1,,,N701,CGATGTAT,N501,TCTTTCCC,Carl_Dewey,
2,FG2,Sample_FG2,,,N702,TGACCAAT,N502,TCTTTCCC,Filipe_Greer,
"""
        s = SampleSheet(fp=StringIO(sample_sheet_data))
        s1 = set_samplesheet_column(s,"SAMPLE_PROJECT","Project1")
        for line in s1:
            self.assertEqual(line['Sample_Project'],"Project1")
    def test_set_samplesheet_column_project_to_value_for_lane(self):
        """
        set_samplesheet_column: set project to specific value for a lane
        """
        sample_sheet_data = u"""[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,AB1,Sample_AB1,,,N701,CGATGTAT,N501,TCTTTCCC,Andrew_Bloggs,
1,AB2,Sample_AB2,,,N702,TGACCAAT,N502,TCTTTCCC,Andrew_Bloggs,
2,CD1,Sample_CD1,,,N701,CGATGTAT,N501,TCTTTCCC,Carl_Dewey,
2,FG2,Sample_FG2,,,N702,TGACCAAT,N502,TCTTTCCC,Filipe_Greer,
"""
        s = SampleSheet(fp=StringIO(sample_sheet_data))
        s1 = set_samplesheet_column(s,"SAMPLE_PROJECT","Project1",
                                    lanes=[1])
        for line in s1:
            if line['Lane'] == 1:
                self.assertEqual(line['Sample_Project'],"Project1")
            else:
                self.assertNotEqual(line['Sample_Project'],"Project1")
    def test_set_samplesheet_column_names_to_ids(self):
        """
        set_samplesheet_column: set sample names to sample IDs
        """
        sample_sheet_data = u"""[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,AB1,Sample_AB1,,,N701,CGATGTAT,N501,TCTTTCCC,Andrew_Bloggs,
1,AB2,Sample_AB2,,,N702,TGACCAAT,N502,TCTTTCCC,Andrew_Bloggs,
2,CD1,Sample_CD1,,,N701,CGATGTAT,N501,TCTTTCCC,Carl_Dewey,
2,FG2,Sample_FG2,,,N702,TGACCAAT,N502,TCTTTCCC,Filipe_Greer,
"""
        s = SampleSheet(fp=StringIO(sample_sheet_data))
        s1 = set_samplesheet_column(s,"SAMPLE_NAME","SAMPLE_ID")
        for line in s1:
            self.assertEqual(line['Sample_Name'],line['Sample_ID'])
    def test_set_samplesheet_column_ids_to_names(self):
        """
        set_samplesheet_column: set sample IDs to sample names
        """
        sample_sheet_data = u"""[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,AB1,Sample_AB1,,,N701,CGATGTAT,N501,TCTTTCCC,Andrew_Bloggs,
1,AB2,Sample_AB2,,,N702,TGACCAAT,N502,TCTTTCCC,Andrew_Bloggs,
2,CD1,Sample_CD1,,,N701,CGATGTAT,N501,TCTTTCCC,Carl_Dewey,
2,FG2,Sample_FG2,,,N702,TGACCAAT,N502,TCTTTCCC,Filipe_Greer,
"""
        s = SampleSheet(fp=StringIO(sample_sheet_data))
        s1 = set_samplesheet_column(s,"SAMPLE_ID","SAMPLE_NAME")
        for line in s1:
            self.assertEqual(line['Sample_Name'],line['Sample_ID'])
