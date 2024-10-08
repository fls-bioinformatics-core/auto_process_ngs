#######################################################################
# Unit tests for qc/qualimap.py
#######################################################################

import unittest
import os
import tempfile
import shutil

from auto_process_ngs.mockqcdata import QUALIMAP_RNASEQ_RESULTS
from auto_process_ngs.qc.apps.qualimap import QualimapRnaseq
from auto_process_ngs.qc.apps.qualimap import qualimap_rnaseq_output

class TestQualimapRnaseq(unittest.TestCase):
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestQualimapRnaseq')
        # Make an rnaseq_qc_results.txt file
        self.rnaseq_qc_results_txt = os.path.join(self.dirn,
                                                  "rnaseq_qc_results.txt")
        with open(self.rnaseq_qc_results_txt,'wt') as fp:
            fp.write(QUALIMAP_RNASEQ_RESULTS)

    def tearDown(self):
        # Remove the temporary test directory
        shutil.rmtree(self.dirn)

    def test_qualimap_rnaseq(self):
        """
        QualimapRnaseq: read in data from 'qualimap rnaseq'
        """
        rnaseq_results = QualimapRnaseq(self.dirn)
        self.assertEqual(rnaseq_results.html_report,
                         os.path.join(self.dirn,"qualimapReport.html"))
        self.assertEqual(rnaseq_results.results_txt,
                         self.rnaseq_qc_results_txt)
        self.assertEqual(rnaseq_results.input['bam file'],
                         "/path/to/SMP1_S1_R1_001.bam")
        self.assertEqual(rnaseq_results.input['gff file'],
                         "/path/to/annotation.gtf")
        self.assertEqual(rnaseq_results.input['counting algorithm'],
                         "uniquely-mapped-reads")
        self.assertEqual(rnaseq_results.input['protocol'],
                         "strand-specific-reverse")
        self.assertEqual(rnaseq_results.\
                         reads_genomic_origin['exonic'],
                         (63784,76.93))
        self.assertEqual(rnaseq_results.\
                         reads_genomic_origin['intronic'],
                         (4778,5.76))
        self.assertEqual(rnaseq_results.\
                         reads_genomic_origin['intergenic'],
                         (14352,17.31))
        self.assertEqual(rnaseq_results.\
                         reads_genomic_origin['overlapping exon'],
                         (6523,7.87))
        self.assertEqual(rnaseq_results.reads_genomic_origin['rRNA'],
                         (0,0.0))
        self.assertEqual(rnaseq_results.link_to_output("Reads Genomic Origin"),
                         os.path.join(self.dirn,
                                      "qualimapReport.html#Reads%20Genomic%20Origin"))

class TestQualimapRnaseqOutputFunction(unittest.TestCase):

    def test_qualimap_rnaseq_output(self):
        """
        qualimap_rnaseq_output: no prefix
        """
        self.assertEqual(qualimap_rnaseq_output(),
                         ('qualimapReport.html',
                          'rnaseq_qc_results.txt'))

    def test_qualimap_rnaseq_output_with_prefix(self):
        """
        qualimap_rnaseq_output: with prefix
        """
        self.assertEqual(
            qualimap_rnaseq_output(
                prefix="qualimap-rnaseq/human/PB1_ATTAGG_L001_R1_001"),
            ('qualimap-rnaseq/human/PB1_ATTAGG_L001_R1_001/qualimapReport.html',
             'qualimap-rnaseq/human/PB1_ATTAGG_L001_R1_001/rnaseq_qc_results.txt'))
