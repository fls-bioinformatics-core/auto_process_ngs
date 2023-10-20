#######################################################################
# Tests for config.py module
#######################################################################

import unittest
from io import StringIO
from bcftbx.JobRunner import SimpleJobRunner,GEJobRunner
from auto_process_ngs.config import *

config = """
[basic]
greeting = hello
farewell = goodbye

[advanced]
bool = true
n = 12
m = 42.0

[none_types]
some_value = None
other_value = 

[runners]
my_runner = SimpleJobRunner

[names]
Name1 = Peter
name1 = John
"""

class TestConfig(unittest.TestCase):
    """Tests for the Config class
    """
    def setUp(self):
        fp = StringIO(config)
        self.config = Config()
        self.config.read_file(fp)

    def test_get(self):
        """
        Config: 'get' fetches correct value
        """
        self.assertEqual(self.config.get('basic','greeting'),'hello')

    def test_get_with_default(self):
        """
        Config: 'get' works with defaults
        """
        self.assertEqual(self.config.get('basic','salutation',
                                         default='bonjour'),
                         'bonjour')

    def test_get_with_null(self):
        """
        Config: 'get' returns NullValue for missing section/option
        """
        self.assertEqual(self.config.get('missing','value'),NullValue())
        self.assertEqual(self.config.get('basic','salutation'),NullValue())

    def test_getint(self):
        """
        Config: 'getint' fetches integer value
        """
        self.assertEqual(self.config.getint('advanced','n'),12)

    def test_getint_with_default(self):
        """
        Config: 'getint' works with defaults
        """
        self.assertEqual(self.config.getint('advanced','p',default=11),
                         11)

    def test_getint_with_null(self):
        """
        Config: 'getint' returns NullValue for missing section/option
        """
        self.assertEqual(self.config.getint('missing','value'),NullValue())
        self.assertEqual(self.config.getint('advanced','p'),NullValue())

    def test_getfloat(self):
        """
        Config: 'getfloat' fetches float value
        """
        self.assertEqual(self.config.getfloat('advanced','m'),42.0)

    def test_getfloat_with_default(self):
        """
        Config: 'getfloat' works with defaults
        """
        self.assertEqual(self.config.getfloat('advanced','p',default=5.0),
                         5.0)

    def test_getfloat_with_null(self):
        """
        Config: 'getfloat' returns NullValue for missing section/option
        """
        self.assertEqual(self.config.getint('missing','value'),NullValue())
        self.assertEqual(self.config.getfloat('advanced','p'),NullValue())

    def test_getboolean(self):
        """
        Config: 'getboolean' fetches boolean value
        """
        self.assertEqual(self.config.getboolean('advanced','bool'),True)

    def test_getboolean_with_default(self):
        """
        Config: 'getboolean' works with defaults
        """
        self.assertEqual(self.config.getboolean('advanced','p',default=True),
                         True)

    def test_getboolean_with_null(self):
        """
        Config: 'getboolean' returns NullValue for missing section/option
        """
        self.assertEqual(self.config.getboolean('missing','value'),NullValue())
        self.assertEqual(self.config.getboolean('advanced','p'),NullValue())

    def test_getrunner(self):
        """
        Config: 'getrunner' fetches correct value
        """
        self.assertTrue(
            isinstance(self.config.getrunner('runners',
                                             'my_runner'),
                       SimpleJobRunner))

    def test_getrunner_with_default(self):
        """
        Config: 'getrunner' works with defaults
        """
        self.assertTrue(
            isinstance(self.config.getrunner('runners',
                                             'your_runner',
                                             default='GEJobRunner'),
                       GEJobRunner))

    def test_getrunner_with_runner_instance_as_default(self):
        """
        Config: 'getrunner' works with runner instance as default
        """
        self.assertTrue(
            isinstance(self.config.getrunner('runners',
                                             'your_runner',
                                             default=GEJobRunner()),
                       GEJobRunner))

    def test_getrunner_with_null(self):
        """
        Config: 'getrunner' returns NullValue for missing section/option
        """
        self.assertEqual(self.config.getrunner('missing','value'),
                         NullValue())
        self.assertEqual(self.config.getrunner('runners','your_runner'),
                         NullValue())

    def test_get_with_None_value(self):
        """
        Config: 'get' works for parameters set to None
        """
        self.assertEqual(self.config.get('none_types','some_value'),None)
        self.assertEqual(self.config.get('none_types','some_value',
                                         default='something'),
                         'something')

    def test_get_with_empty_value(self):
        """
        Config: 'get' works for parameters with no value set
        """
        self.assertEqual(self.config.get('none_types','other_value'),None)
        self.assertEqual(self.config.get('none_types','other_value',
                                         default='something'),
                         'something')

    def test_option_case_is_preserved(self):
        """
        Config: 'get' preserves option case
        """
        self.assertEqual(self.config.get('names','Name1'),'Peter')
        self.assertEqual(self.config.get('names','name1'),'John')

class TestNullValue(unittest.TestCase):
    """
    Tests for the NullValue class
    """
    def test_nullvalue_repr(self):
        """
        NullValue: check __repr__
        """
        self.assertEqual(repr(NullValue()),"<NullValue>")

    def test_nullvalue_eq(self):
        """
        NullValue: check __eq__
        """
        nullvalue1 = NullValue()
        nullvalue2 = NullValue()
        self.assertTrue(nullvalue1 == nullvalue2)

    def test_nullvalue_bool(self):
        """
        NullValue: check __bool__
        """
        self.assertFalse(bool(NullValue()))
