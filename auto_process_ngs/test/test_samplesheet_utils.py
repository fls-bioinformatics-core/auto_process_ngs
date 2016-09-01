#######################################################################
# Tests for samplesheet_utils.py module
#######################################################################

import os
import unittest
import tempfile
import shutil
import codecs
from auto_process_ngs.samplesheet_utils import get_close_names
from auto_process_ngs.samplesheet_utils import has_invalid_characters

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

class TestHasInvalidCharacters(unittest.TestCase):
    def setUp(self):
        self.wd = tempfile.mkdtemp()
    def tearDown(self):
        if os.path.isdir(self.wd):
            shutil.rmtree(self.wd)
    def test_file_with_non_printing_ascii_character(self):
        """has_invalid_characters: detects non-printing ASCII character
        """
        non_printing_ascii_file = os.path.join(self.wd,
                                               "test.nonprintingascii")
        with open(non_printing_ascii_file,'wb') as fp:
            fp.write(u"""This file contains:
a non-printing ASCII ctrl-S character here\x13
""")
        self.assertTrue(has_invalid_characters(non_printing_ascii_file))
    def test_file_with_non_ascii_character(self):
        """has_invalid_characters: detects non-ASCII character
        """
        non_ascii_file = os.path.join(self.wd,"test.nonascii")
        with codecs.open(non_ascii_file,'wb',encoding='utf-8') as fp:
            fp.write(u"""This file contains:
a non-ASCII character here\x80
""")
        self.assertTrue(has_invalid_characters(non_ascii_file))
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
        
