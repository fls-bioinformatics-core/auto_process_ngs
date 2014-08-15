#######################################################################
# Tests for bcl2fastq_utils.py module
#######################################################################

import unittest
from auto_process_ngs.bcl2fastq_utils import get_nmismatches

class TestGetNmismatches(unittest.TestCase):
    """Tests for the get_nmismatches function

    """
    def test_n_mismatches(self):
        self.assertEqual(get_nmismatches('y50'),0)
        self.assertEqual(get_nmismatches('y50,I4'),0)
        self.assertEqual(get_nmismatches('y50,I6'),1)
        self.assertEqual(get_nmismatches('y101,I6,y101'),1)
        self.assertEqual(get_nmismatches('y250,I8,I8,y250'),1)
        self.assertEqual(get_nmismatches('y250,I6nn,I6nn,y250'),1)
        self.assertEqual(get_nmismatches('y250,I6n2,I6n2,y250'),1)
        self.assertEqual(get_nmismatches('y250,I16'),1)

