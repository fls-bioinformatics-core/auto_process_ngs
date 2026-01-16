#######################################################################
# Unit tests for indexes.py
#######################################################################

import unittest
import os
import tempfile
import shutil
from bcftbx.JobRunner import SimpleJobRunner
from auto_process_ngs.mock import MockBowtieBuild
from auto_process_ngs.mock import MockBowtie2Build
from auto_process_ngs.mock import MockStar
from auto_process_ngs.indexes import IndexBuilder
from auto_process_ngs.indexes import bowtie_build_cmd
from auto_process_ngs.indexes import bowtie2_build_cmd
from auto_process_ngs.indexes import star_genome_generate_cmd

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

# Tests
class TestIndexBuilder(unittest.TestCase):

    def setUp(self):
        # Temporary working dir
        self.wd = tempfile.mkdtemp(suffix='TestIndexBuilder')
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
        """
        IndexBuilder: build STAR indexes
        """
        MockStar.create(os.path.join(self.bin,"STAR"))
        builder = IndexBuilder(SimpleJobRunner(),use_conda=False)
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

    def test_build_bowtie_index(self):
        """
        IndexBuilder: build Bowtie indexes
        """
        MockBowtieBuild.create(os.path.join(self.bin,"bowtie-build"))
        builder = IndexBuilder(SimpleJobRunner(),use_conda=False)
        retcode = builder.bowtie("/data/example.fasta",
                                 os.path.join(self.wd,"bowtie_index"))
        self.assertEqual(retcode,0)
        self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                    "bowtie_index")))
        for f in ("example.1.ebwt",
                  "example.2.ebwt",
                  "example.3.ebwt",
                  "example.4.ebwt",
                  "example.rev.1.ebwt",
                  "example.rev.2.ebwt"):
            self.assertTrue(os.path.exists(
                os.path.join(self.wd,
                             "bowtie_index",
                             f)),
                            "Missing output file: %s" % f)

    def test_build_bowtie2_index(self):
        """
        IndexBuilder: build Bowtie2 indexes
        """
        MockBowtie2Build.create(os.path.join(self.bin,"bowtie2-build"))
        builder = IndexBuilder(SimpleJobRunner(),use_conda=False)
        retcode = builder.bowtie2("/data/example.fasta",
                                  os.path.join(self.wd,"bowtie2_index"))
        self.assertEqual(retcode,0)
        self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                    "bowtie2_index")))
        for f in ("example.1.bt2",
                  "example.2.bt2",
                  "example.3.bt2",
                  "example.4.bt2",
                  "example.rev.1.bt2",
                  "example.rev.2.bt2"):
            self.assertTrue(os.path.exists(
                os.path.join(self.wd,
                             "bowtie2_index",
                             f)),
                            "Missing output file: %s" % f)

class TestBowtieBuildCmdFunction(unittest.TestCase):

    def test_bowtie_build_cmd(self):
        """
        bowtie_build_cmd: check index generation command
        """
        self.assertEqual(
            str(bowtie_build_cmd("/data/example.fasta","example")),
            "bowtie-build -f /data/example.fasta example")

class TestBowtie2BuildCmdFunction(unittest.TestCase):

    def test_bowtie2_build_cmd(self):
        """
        bowtie2_build_cmd: check index generation command
        """
        self.assertEqual(
            str(bowtie2_build_cmd("/data/example.fasta",
                                  "example")),
            "bowtie2-build -f /data/example.fasta example")

    def test_bowtie2_build_cmd_nthreads(self):
        """
        bowtie2_build_cmd: check index generation command with threads
        """
        self.assertEqual(
            str(bowtie2_build_cmd("/data/example.fasta",
                                  "example",
                                  nthreads=8)),
            "bowtie2-build --threads 8 -f /data/example.fasta example")

class TestStarGenomeGenerateCmdFunction(unittest.TestCase):

    def test_star_genome_generate_cmd_basic(self):
        """
        star_genome_generate_cmd: check basic index generation command
        """
        self.assertEqual(
            str(star_genome_generate_cmd("/data/example.fasta",
                                         "/data/example.gtf",
                                         "/indexes/example/star")),
            "STAR --runMode genomeGenerate "
            "--genomeFastaFiles /data/example.fasta "
            "--sjdbGTFfile /data/example.gtf "
            "--genomeDir /indexes/example/star")

    def test_star_genome_generate_cmd_full(self):
        """
        star_genome_generate_cmd: check full index generation command
        """
        self.assertEqual(
            str(star_genome_generate_cmd("/data/example.fasta",
                                         "/data/example.gtf",
                                         "/indexes/example/star",
                                         overhang=100,
                                         sa_index_nbases=14,
                                         nthreads=12,
                                         memory_limit=32)),
            "STAR --runMode genomeGenerate "
            "--genomeFastaFiles /data/example.fasta "
            "--sjdbGTFfile /data/example.gtf "
            "--genomeDir /indexes/example/star "
            "--sjdbOverhang 100 "
            "--genomeSAindexNbases 14 "
            "--runThreadN 12 "
            "--limitGenomeGenerateRAM 32")
