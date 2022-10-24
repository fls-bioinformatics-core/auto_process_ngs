#######################################################################
# Unit tests for indexes.py
#######################################################################

import unittest
import os
import tempfile
import shutil
from bcftbx.JobRunner import SimpleJobRunner
from auto_process_ngs.mock import MockStar
from auto_process_ngs.indexes import IndexBuilder

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

# Tests
class TestIndexBuilderSTAR(unittest.TestCase):

    def setUp(self):
        # Temporary working dir
        self.wd = tempfile.mkdtemp(suffix='TestIndexes')
        # Temporary 'bin' dir
        self.bin = os.path.join(self.wd,"bin")
        os.mkdir(self.bin)
        # Store original PATH
        self.path = os.environ['PATH']
        # Add 'bin' to PATH
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])

    def tearDown(self):
        # Restore PATH
        os.environ['PATH'] = self.path
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    def test_build_star_index(self):
        MockStar.create(os.path.join(self.bin,"STAR"))
        builder = IndexBuilder(SimpleJobRunner())
        retcode = builder.STAR("/data/example.fasta",
                               "/data/example.gtf",
                               os.path.join(self.wd,"star_index"))
        self.assertEqual(retcode,0)
        self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                    "star_index")))
        for f in ("chrLength.txt",
                  "chrName.txt",
                  "exonGeTrInfo.tab",
                  "geneInfo.tab",
                  "genomeParameters.txt",
                  "Log.out",
                  "SAindex",
                  "sjdbList.fromGTF.out.tab",
                  "transcriptInfo.tab",
                  "chrNameLength.txt",
                  "chrStart.txt",
                  "exonInfo.tab",
                  "Genome",
                  "SA",
                  "sjdbInfo.txt",
                  "sjdbList.out.tab"):
            self.assertTrue(os.path.exists(
                os.path.join(self.wd,
                             "star_index",
                             f)),
                            "Missing output file: %s" % f)
