#######################################################################
# Unit tests for qc/rseqc.py
#######################################################################

import unittest
import os
import tempfile
import shutil

from auto_process_ngs.qc.rseqc import InferExperiment

class TestInferExperiment(unittest.TestCase):
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestRSeQC')

    def tearDown(self):
        # Remove the temporary test directory
        shutil.rmtree(self.dirn)

    def _make_file(self,f,data):
        # Write data to a file
        with open(os.path.join(self.dirn,f),'wt') as fp:
            fp.write(data)
        return os.path.join(self.dirn,f)

    def test_infer_experiment_paired_end(self):
        """
        InferExperiment: paired-end output
        """
        log = self._make_file("infer_experiment.log",
                              """This is PairEnd Data
Fraction of reads failed to determine: 0.0172
Fraction of reads explained by "1++,1--,2+-,2-+": 0.4903
Fraction of reads explained by "1+-,1-+,2++,2--": 0.4925
""")
        infer_expt = InferExperiment(log)
        self.assertEqual(infer_expt.log_file,log)
        self.assertTrue(infer_expt.paired_end)
        self.assertEqual(infer_expt.forward,0.4903)
        self.assertEqual(infer_expt.reverse,0.4925)
        self.assertEqual(infer_expt.unstranded,0.0172)

    def test_infer_experiment_single_end(self):
        """
        InferExperiment: single-end output
        """
        log = self._make_file("infer_experiment.log",
                              """This is SingleEnd Data
Fraction of reads failed to determine: 0.0170
Fraction of reads explained by "++,--": 0.9669
Fraction of reads explained by "+-,-+": 0.0161
""")
        infer_expt = InferExperiment(log)
        self.assertEqual(infer_expt.log_file,log)
        self.assertFalse(infer_expt.paired_end)
        self.assertEqual(infer_expt.forward,0.9669)
        self.assertEqual(infer_expt.reverse,0.0161)
        self.assertEqual(infer_expt.unstranded,0.0170)
