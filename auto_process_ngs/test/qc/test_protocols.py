#######################################################################
# Unit tests for qc/protcols.py
#######################################################################

import unittest
import os
import tempfile
import shutil
from auto_process_ngs.qc.protocols import fetch_protocol_definition

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

class TestFetchProtocolDefinition(unittest.TestCase):

    def test_fetch_protocol_definition(self):
        """
        fetch_protocol_definition: check definition is returned
        """
        reads,qc_modules = fetch_protocol_definition("standardPE")
        self.assertEqual(reads.data,['r1','r2'])
        self.assertEqual(reads.qc,['r1','r2'])
        self.assertEqual(qc_modules,['fastqc',
                                     'fastq_screen',
                                     'sequence_lengths',
                                     'strandedness'])

    def test_fetch_protocol_definition_unknown_protocol(self):
        """
        fetch_protocol_definition: raises KeyError for unknown protocol
        """
        self.assertRaises(KeyError,
                          fetch_protocol_definition,
                          "whazzdis?")
