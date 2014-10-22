#######################################################################
# Tests for config.py module
#######################################################################

import unittest
import cStringIO
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
"""

class TestConfig(unittest.TestCase):
    """Tests for the Config class
    """
    def setUp(self):
        fp = cStringIO.StringIO(config)
        self.config = Config()
        self.config.readfp(fp)

    def test_get(self):
        """Check Config.get fetches correct value
        """
        self.assertEqual(self.config.get('basic','greeting'),'hello')

    def test_get_with_default(self):
        """Check Config.get works with defaults
        """
        self.assertEqual(self.config.get('basic','salutation'),None)
        self.assertEqual(self.config.get('basic','salutation','bonjour'),'bonjour')

    def test_getint(self):
        """Check Config.getint fetches integer value
        """
        self.assertEqual(self.config.getint('advanced','n'),12)

    def test_getint_with_default(self):
        """Check Config.getint works with defaults
        """
        self.assertEqual(self.config.getint('advanced','p'),None)
        self.assertEqual(self.config.getint('advanced','p',11),11)

    def test_getrunner(self):
        """Check Config.getrunner fetches correct value
        """
        self.assertTrue(isinstance(self.config.getrunner('runners','my_runner'),
                                   SimpleJobRunner))

    def test_getrunner_with_default(self):
        """Check Config.getrunner works with defaults
        """
        self.assertTrue(isinstance(self.config.getrunner('runners','your_runner'),
                                   SimpleJobRunner))
        self.assertTrue(isinstance(self.config.getrunner('runners','your_runner','GEJobRunner'),
                                   GEJobRunner))

    def test_get_with_None_value(self):
        """Check Config.get works for parameters set to 'None'
        """
        self.assertEqual(self.config.get('none_types','some_value'),None)
        self.assertEqual(self.config.get('none_types','some_value','something'),'something')

    def test_get_with_empty_value(self):
        """Check Config.get works for parameters with no value set'
        """
        self.assertEqual(self.config.get('none_types','other_value'),None)
        self.assertEqual(self.config.get('none_types','other_value','something'),'something')
