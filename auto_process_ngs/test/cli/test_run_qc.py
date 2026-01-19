#######################################################################
# Tests for cli/run_qc.py utility
#######################################################################

import unittest
import tempfile
import shutil
import os
import gzip
from auto_process_ngs.analysis import AnalysisFastq
from auto_process_ngs.cli.run_qc import main as run_qc
from auto_process_ngs.metadata import AnalysisProjectQCDirInfo
from auto_process_ngs.mock import MockAnalysisProject
from auto_process_ngs.mock import MockFastqScreen
from auto_process_ngs.mock import MockFastQC
from auto_process_ngs.mock import MockGtf2bed
from auto_process_ngs.mock import MockSeqtk
from auto_process_ngs.mock import MockStar
from auto_process_ngs.mock import MockSamtools
from auto_process_ngs.mock import MockPicard
from auto_process_ngs.mock import MockRSeQC
from auto_process_ngs.mock import MockQualimap
from auto_process_ngs.mock import MockMultiQC
from auto_process_ngs.mock import MockCellrangerExe
from auto_process_ngs.qc.protocols import fetch_protocol_definition

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

# Unit tests

class TestRunQc(unittest.TestCase):
    """
    Tests for the 'run_qc' utility
    """
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestRunQc')
        # Create a temp 'bin' dir
        self.bin = os.path.join(self.dirn, "bin")
        os.mkdir(self.bin)
        # Create a temp 'data' dir
        self.data = os.path.join(self.dirn, "data")
        os.mkdir(self.data)
        # Add (empty) FastqScreen conf files
        self.fastq_screens = {
            'model_organisms': "fastq_screen_model_organisms.conf",
            'other_organisms': "fastq_screen_other_organisms.conf",
            'rRNA': "fastq_screen_rRNA.conf",
        }
        for screen in self.fastq_screens:
            conf_file = os.path.join(self.data,
                                     self.fastq_screens[screen])
            with open(conf_file,'wt') as fp:
                fp.write("")
            self.fastq_screens[screen] = conf_file
        # Add (empty) reference data files
        self.ref_data = dict()
        for build in ("hg38", "mm10",):
            self.ref_data[build] = {}
            build_dir = os.path.join(self.data,build)
            os.mkdir(build_dir)
            for ext in ("bed", "gtf"):
                f = os.path.join(build_dir,"%s.%s" % (build, ext))
                with open(f,'wt') as fp:
                    fp.write("")
                self.ref_data[build][ext] = f
        # Create settings file and add reference data info
        settings_ini = os.path.join(self.dirn, "auto_process.ini")
        with open(settings_ini, "wt") as s:
            s.write(f"""[general]
poll_interval = 0.1
max_concurrent_jobs = 1

[qc]
fastq_screens = model_organisms,other_organisms,rRNA

[screen:model_organisms]
conf_file = {self.fastq_screens['model_organisms']}

[screen:other_organisms]
conf_file = {self.fastq_screens['other_organisms']}

[screen:rRNA]
conf_file = {self.fastq_screens['rRNA']}

[organism:human]
star_index = /data/hg38/star_index
annotation_bed = {self.ref_data['hg38']['bed']}
annotation_gtf = {self.ref_data['hg38']['gtf']}

[organism:mouse]
star_index = /data/mm10/star_index
annotation_bed = {self.ref_data['mm10']['bed']}
annotation_gtf = {self.ref_data['mm10']['gtf']}
""")
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        MockStar.create(os.path.join(self.bin,"STAR"))
        MockSamtools.create(os.path.join(self.bin,"samtools"))
        MockPicard.create(os.path.join(self.bin,"picard"))
        MockGtf2bed.create(os.path.join(self.bin,"gtf2bed"))
        MockRSeQC.create(os.path.join(self.bin,"infer_experiment.py"))
        MockRSeQC.create(os.path.join(self.bin,"geneBody_coverage.py"))
        MockQualimap.create(os.path.join(self.bin,"qualimap"))
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        # Store original PATH, then extend with "bin"
        self.path = os.environ['PATH']
        os.environ['PATH'] = "%s:%s" % (self.bin, os.environ['PATH'])
        # Store original location so we can get back at the end
        self.pwd = os.getcwd()
        # Move to working dir
        os.chdir(self.dirn)

    def tearDown(self):
        # Return to original dir
        os.chdir(self.pwd)
        # Restore PATH
        os.environ['PATH'] = self.path
        # Remove the temporary test directory
        def del_rw(action,name,excinfo):
            # Explicitly remove read only files/
            # dirs
            if os.path.isfile(name):
                os.chmod(os.path.dirname(name),0o755)
                os.chmod(name,0o655)
                os.remove(name)
            elif os.path.isdir(name):
                os.chmod(os.path.dirname(name),0o755)
                os.chmod(name,0o755)
                os.rmdir(name)
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.dirn,onerror=del_rw)

    def _mock_fastq(self, fastq):
        """
        Makes a mock Fastq file
        """
        fastq = os.path.abspath(fastq)
        fq = os.path.basename(fastq)
        read_number = AnalysisFastq(fq).read_number
        lane = AnalysisFastq(fq).lane_number
        if lane is None:
            lane = 1
            read = """@ILLUMINA-545855:49:FC61RLR:%s:1:10979:1695 %s:N:0:TCCTGA
GCATACTCAGCTTTAGTAATAAGTGTGATTCTGGTA
+
IIIIIHIIIGHHIIDGHIIIIIIHIIIIIIIIIIIH\n""" % (lane, read_number)
        if fq.endswith('.gz'):
            open_func = gzip.open
        else:
            open_func = open
        with open_func(fastq, "wt") as fp:
            fp.write(read)
        return fastq

    def test_run_qc_paired_end_project(self):
        """
        run_qc.py: paired-end project
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={ 'Organism': 'Human' })
        project_dir = p.create(top_dir=self.dirn)
        # Implicit output and QC directories
        out_dir = project_dir
        qc_dir = os.path.join(out_dir, "qc")
        # Run the QC
        self.assertEqual(run_qc([project_dir,]), 0)
        # Check output and reports
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.zip",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(out_dir, f)),
                            f"Missing '{f}' under {out_dir}")
        # Check collated Picard insert sizes
        collated_insert_sizes = os.path.join(qc_dir,
                                             "insert_sizes.human.tsv")
        self.assertTrue(os.path.exists(collated_insert_sizes),
                        "Missing collated insert sizes TSV")
        with open(collated_insert_sizes,'rt') as fp:
            self.assertEqual(fp.read(),
                             """#Bam file	Mean insert size	Standard deviation	Median insert size	Median absolute deviation
PJB1_S1_001.bam	153.754829	69.675347	139	37
PJB2_S2_001.bam	153.754829	69.675347	139	37
""")
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(os.path.join(qc_dir, "qc.info"))
        self.assertEqual(qc_info.protocol,"standardPE")
        self.assertEqual(qc_info.protocol_specification,
                         str(fetch_protocol_definition("standardPE")))
        self.assertEqual(qc_info.organism,"Human")
        self.assertEqual(qc_info.seq_data_samples,"PJB1,PJB2")
        self.assertEqual(qc_info.fastq_dir, None)
        self.assertEqual(qc_info.fastqs,
                         "PJB1_S1_R1_001.fastq.gz,"
                         "PJB1_S1_R2_001.fastq.gz,"
                         "PJB2_S2_R1_001.fastq.gz,"
                         "PJB2_S2_R2_001.fastq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane,False)
        self.assertEqual(qc_info.fastq_screens,
                         "model_organisms,other_organisms,rRNA")
        self.assertEqual(qc_info.star_index,"/data/hg38/star_index")
        self.assertEqual(qc_info.annotation_bed,self.ref_data['hg38']['bed'])
        self.assertEqual(qc_info.annotation_gtf,self.ref_data['hg38']['gtf'])
        self.assertEqual(qc_info.cellranger_version,None)
        self.assertEqual(qc_info.cellranger_refdata,None)
        self.assertEqual(qc_info.cellranger_probeset,None)

    def test_run_qc_single_end_project(self):
        """
        run_qc.py: single-end project
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz"),
                                metadata={ 'Organism': 'Human' })
        project_dir = p.create(top_dir=self.dirn)
        # Implicit output and QC directories
        out_dir = project_dir
        qc_dir = os.path.join(out_dir, "qc")
        # Run the QC
        self.assertEqual(run_qc([project_dir,]), 0)
        # Check output and reports
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.zip",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(out_dir, f)),
                            f"Missing '{f}' under {out_dir}")
        # Check collated Picard insert sizes
        collated_insert_sizes = os.path.join(qc_dir,
                                             "insert_sizes.human.tsv")
        self.assertFalse(os.path.exists(collated_insert_sizes),
                        "Collated insert sizes TSV shouldn't be present")
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(os.path.join(qc_dir, "qc.info"))
        self.assertEqual(qc_info.protocol,"standardSE")
        self.assertEqual(qc_info.protocol_specification,
                         str(fetch_protocol_definition("standardSE")))
        self.assertEqual(qc_info.organism,"Human")
        self.assertEqual(qc_info.seq_data_samples,"PJB1,PJB2")
        self.assertEqual(qc_info.fastq_dir, None)
        self.assertEqual(qc_info.fastqs,
                         "PJB1_S1_R1_001.fastq.gz,"
                         "PJB2_S2_R1_001.fastq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane,False)
        self.assertEqual(qc_info.fastq_screens,
                         "model_organisms,other_organisms,rRNA")
        self.assertEqual(qc_info.star_index,"/data/hg38/star_index")
        self.assertEqual(qc_info.annotation_bed,self.ref_data['hg38']['bed'])
        self.assertEqual(qc_info.annotation_gtf,self.ref_data['hg38']['gtf'])
        self.assertEqual(qc_info.cellranger_version,None)
        self.assertEqual(qc_info.cellranger_refdata,None)
        self.assertEqual(qc_info.cellranger_probeset,None)

    def test_run_qc_paired_end_fastqs_from_dir_no_metadata(self):
        """
        run_qc.py: paired-end Fastqs in directory (no metadata)
        """
        # Make directory with Fastqs
        fastq_dir = os.path.join(self.dirn, "fastqs")
        os.mkdir(fastq_dir)
        for fq in ("PJB1_S1_R1_001.fastq.gz",
                   "PJB1_S1_R2_001.fastq.gz",
                   "PJB2_S2_R1_001.fastq.gz",
                   "PJB2_S2_R2_001.fastq.gz"):
            self._mock_fastq(os.path.join(fastq_dir, fq))
        # Implicit output and QC directories
        out_dir = os.path.join(self.dirn, "fastqs")
        qc_dir = os.path.join(out_dir, "qc")
        # Run the QC
        self.assertEqual(run_qc([fastq_dir,]), 0)
        # Check output and reports
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.fastqs.zip",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(out_dir, f)),
                            f"Missing '{f}' under {out_dir}")
        # Check collated Picard insert sizes
        collated_insert_sizes = os.path.join(qc_dir,
                                             "insert_sizes.human.tsv")
        self.assertFalse(os.path.exists(collated_insert_sizes),
                        "Collated insert sizes TSV shouldn't be present")
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(os.path.join(qc_dir, "qc.info"))
        self.assertEqual(qc_info.protocol, "minimal")
        self.assertEqual(qc_info.organism, None)
        self.assertEqual(qc_info.seq_data_samples, "PJB1,PJB2")
        self.assertEqual(qc_info.fastq_dir, None)
        self.assertEqual(qc_info.fastqs,
                         "PJB1_S1_R1_001.fastq.gz,"
                         "PJB1_S1_R2_001.fastq.gz,"
                         "PJB2_S2_R1_001.fastq.gz,"
                         "PJB2_S2_R2_001.fastq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane, False)
        self.assertEqual(qc_info.fastq_screens,
                         "model_organisms,other_organisms,rRNA")
        self.assertEqual(qc_info.star_index, None)
        self.assertEqual(qc_info.annotation_bed, None)
        self.assertEqual(qc_info.annotation_gtf, None)
        self.assertEqual(qc_info.cellranger_version,None)
        self.assertEqual(qc_info.cellranger_refdata,None)
        self.assertEqual(qc_info.cellranger_probeset,None)

    def test_run_qc_paired_end_fastqs_from_dir_explicit_metadata(self):
        """
        run_qc.py: paired-end Fastqs in directory with metadata
        """
        # Make directory with Fastqs
        fastq_dir = os.path.join(self.dirn, "fastqs")
        os.mkdir(fastq_dir)
        for fq in ("PJB1_S1_R1_001.fastq.gz",
                   "PJB1_S1_R2_001.fastq.gz",
                   "PJB2_S2_R1_001.fastq.gz",
                   "PJB2_S2_R2_001.fastq.gz"):
            self._mock_fastq(os.path.join(fastq_dir, fq))
        # Implicit output and QC directories
        out_dir = os.path.join(self.dirn, "fastqs")
        qc_dir = os.path.join(out_dir, "qc")
        # Run the QC
        self.assertEqual(run_qc([fastq_dir,
                                 "--organism", "Human"]), 0)
        # Check output and reports
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.fastqs.zip",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(out_dir, f)),
                            f"Missing '{f}' under {out_dir}")
        # Check collated Picard insert sizes
        collated_insert_sizes = os.path.join(qc_dir,
                                             "insert_sizes.human.tsv")
        self.assertTrue(os.path.exists(collated_insert_sizes),
                        "Missing collated insert sizes TSV")
        with open(collated_insert_sizes,'rt') as fp:
            self.assertEqual(fp.read(),
                             """#Bam file	Mean insert size	Standard deviation	Median insert size	Median absolute deviation
PJB1_S1_001.bam	153.754829	69.675347	139	37
PJB2_S2_001.bam	153.754829	69.675347	139	37
""")
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(os.path.join(qc_dir, "qc.info"))
        self.assertEqual(qc_info.protocol,"standardPE")
        self.assertEqual(qc_info.protocol_specification,
                         str(fetch_protocol_definition("standardPE")))
        self.assertEqual(qc_info.organism,"Human")
        self.assertEqual(qc_info.seq_data_samples,"PJB1,PJB2")
        self.assertEqual(qc_info.fastq_dir, None)
        self.assertEqual(qc_info.fastqs,
                         "PJB1_S1_R1_001.fastq.gz,"
                         "PJB1_S1_R2_001.fastq.gz,"
                         "PJB2_S2_R1_001.fastq.gz,"
                         "PJB2_S2_R2_001.fastq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane,False)
        self.assertEqual(qc_info.fastq_screens,
                         "model_organisms,other_organisms,rRNA")
        self.assertEqual(qc_info.star_index,"/data/hg38/star_index")
        self.assertEqual(qc_info.annotation_bed,self.ref_data['hg38']['bed'])
        self.assertEqual(qc_info.annotation_gtf,self.ref_data['hg38']['gtf'])
        self.assertEqual(qc_info.cellranger_version,None)
        self.assertEqual(qc_info.cellranger_refdata,None)
        self.assertEqual(qc_info.cellranger_probeset,None)

    def test_run_qc_paired_end_fastqs_from_list_no_metadata(self):
        """
        run_qc.py: paired-end Fastqs as list (no metadata)
        """
        # Make directory with Fastqs
        fastq_dir = os.path.join(self.dirn, "fastqs")
        os.mkdir(fastq_dir)
        for fq in ("PJB1_S1_R1_001.fastq.gz",
                   "PJB1_S1_R2_001.fastq.gz",
                   "PJB2_S2_R1_001.fastq.gz",
                   "PJB2_S2_R2_001.fastq.gz"):
            self._mock_fastq(os.path.join(fastq_dir, fq))
        # Implicit output and QC directories
        out_dir = self.dirn
        qc_dir = os.path.join(out_dir, "qc")
        # Run the QC
        self.assertEqual(run_qc(
            [os.path.join(fastq_dir, "PJB1_S1_R1_001.fastq.gz"),
             os.path.join(fastq_dir, "PJB1_S1_R2_001.fastq.gz"),
             os.path.join(fastq_dir, "PJB2_S2_R1_001.fastq.gz"),
             os.path.join(fastq_dir, "PJB2_S2_R2_001.fastq.gz")]), 0)
        # Check output and reports
        for f in ("qc",
                  "qc_report.html",
                  f"qc_report.{os.path.basename(out_dir)}.zip",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(out_dir, f)),
                            f"Missing '{f}' under {out_dir}")
        # Check collated Picard insert sizes
        collated_insert_sizes = os.path.join(qc_dir,
                                             "insert_sizes.human.tsv")
        self.assertFalse(os.path.exists(collated_insert_sizes),
                        "Collated insert sizes TSV shouldn't be present")
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(os.path.join(qc_dir, "qc.info"))
        self.assertEqual(qc_info.protocol, "minimal")
        self.assertEqual(qc_info.organism, None)
        self.assertEqual(qc_info.seq_data_samples, "PJB1,PJB2")
        self.assertEqual(qc_info.fastq_dir, None)
        self.assertEqual(qc_info.fastqs,
                         "PJB1_S1_R1_001.fastq.gz,"
                         "PJB1_S1_R2_001.fastq.gz,"
                         "PJB2_S2_R1_001.fastq.gz,"
                         "PJB2_S2_R2_001.fastq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane, False)
        self.assertEqual(qc_info.fastq_screens,
                         "model_organisms,other_organisms,rRNA")
        self.assertEqual(qc_info.star_index, None)
        self.assertEqual(qc_info.annotation_bed, None)
        self.assertEqual(qc_info.annotation_gtf, None)
        self.assertEqual(qc_info.cellranger_version,None)
        self.assertEqual(qc_info.cellranger_refdata,None)
        self.assertEqual(qc_info.cellranger_probeset,None)

    def test_run_qc_paired_end_fastqs_from_list_explicit_metadata(self):
        """
        run_qc.py: paired-end Fastqs as list with metadata
        """
        # Make directory with Fastqs
        fastq_dir = os.path.join(self.dirn, "fastqs")
        os.mkdir(fastq_dir)
        for fq in ("PJB1_S1_R1_001.fastq.gz",
                   "PJB1_S1_R2_001.fastq.gz",
                   "PJB2_S2_R1_001.fastq.gz",
                   "PJB2_S2_R2_001.fastq.gz"):
            self._mock_fastq(os.path.join(fastq_dir, fq))
        # Implicit output and QC directories
        out_dir = self.dirn
        qc_dir = os.path.join(out_dir, "qc")
        # Run the QC
        self.assertEqual(run_qc(
            [os.path.join(fastq_dir, "PJB1_S1_R1_001.fastq.gz"),
             os.path.join(fastq_dir, "PJB1_S1_R2_001.fastq.gz"),
             os.path.join(fastq_dir, "PJB2_S2_R1_001.fastq.gz"),
             os.path.join(fastq_dir, "PJB2_S2_R2_001.fastq.gz"),
             "--organism", "Human"]), 0)
        # Check output and reports
        for f in ("qc",
                  "qc_report.html",
                  f"qc_report.{os.path.basename(out_dir)}.zip",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(out_dir, f)),
                            f"Missing '{f}' under {out_dir}")
        # Check collated Picard insert sizes
        collated_insert_sizes = os.path.join(qc_dir,
                                             "insert_sizes.human.tsv")
        self.assertTrue(os.path.exists(collated_insert_sizes),
                        "Missing collated insert sizes TSV")
        with open(collated_insert_sizes,'rt') as fp:
            self.assertEqual(fp.read(),
                             """#Bam file	Mean insert size	Standard deviation	Median insert size	Median absolute deviation
PJB1_S1_001.bam	153.754829	69.675347	139	37
PJB2_S2_001.bam	153.754829	69.675347	139	37
""")
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(os.path.join(qc_dir, "qc.info"))
        self.assertEqual(qc_info.protocol,"standardPE")
        self.assertEqual(qc_info.protocol_specification,
                         str(fetch_protocol_definition("standardPE")))
        self.assertEqual(qc_info.organism, "Human")
        self.assertEqual(qc_info.seq_data_samples, "PJB1,PJB2")
        self.assertEqual(qc_info.fastq_dir, None)
        self.assertEqual(qc_info.fastqs,
                         "PJB1_S1_R1_001.fastq.gz,"
                         "PJB1_S1_R2_001.fastq.gz,"
                         "PJB2_S2_R1_001.fastq.gz,"
                         "PJB2_S2_R2_001.fastq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane, False)
        self.assertEqual(qc_info.fastq_screens,
                         "model_organisms,other_organisms,rRNA")
        self.assertEqual(qc_info.star_index,"/data/hg38/star_index")
        self.assertEqual(qc_info.annotation_bed,self.ref_data['hg38']['bed'])
        self.assertEqual(qc_info.annotation_gtf,self.ref_data['hg38']['gtf'])
        self.assertEqual(qc_info.cellranger_version,None)
        self.assertEqual(qc_info.cellranger_refdata,None)
        self.assertEqual(qc_info.cellranger_probeset,None)

    def test_run_qc_paired_end_fastqs_from_wildcard_no_metadata(self):
        """
        run_qc.py: paired-end Fastqs from wildcard (no metadata)
        """
        # Make directory with Fastqs
        fastq_dir = os.path.join(self.dirn, "fastqs")
        os.mkdir(fastq_dir)
        for fq in ("PJB1_S1_R1_001.fastq.gz",
                   "PJB1_S1_R2_001.fastq.gz",
                   "PJB2_S2_R1_001.fastq.gz",
                   "PJB2_S2_R2_001.fastq.gz"):
            self._mock_fastq(os.path.join(fastq_dir, fq))
        # Implicit output and QC directories
        out_dir = self.dirn
        qc_dir = os.path.join(out_dir, "qc")
        # Run the QC
        self.assertEqual(run_qc([os.path.join(fastq_dir, "*"),]), 0)
        # Check output and reports
        for f in ("qc",
                  "qc_report.html",
                  f"qc_report.{os.path.basename(out_dir)}.zip",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(out_dir, f)),
                            f"Missing '{f}' under {out_dir}")
        # Check collated Picard insert sizes
        collated_insert_sizes = os.path.join(qc_dir,
                                             "insert_sizes.human.tsv")
        self.assertFalse(os.path.exists(collated_insert_sizes),
                        "Collated insert sizes TSV shouldn't be present")
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(os.path.join(qc_dir, "qc.info"))
        self.assertEqual(qc_info.protocol, "minimal")
        self.assertEqual(qc_info.organism, None)
        self.assertEqual(qc_info.seq_data_samples, "PJB1,PJB2")
        self.assertEqual(qc_info.fastq_dir, None)
        self.assertEqual(qc_info.fastqs,
                         "PJB1_S1_R1_001.fastq.gz,"
                         "PJB1_S1_R2_001.fastq.gz,"
                         "PJB2_S2_R1_001.fastq.gz,"
                         "PJB2_S2_R2_001.fastq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane, False)
        self.assertEqual(qc_info.fastq_screens,
                         "model_organisms,other_organisms,rRNA")
        self.assertEqual(qc_info.star_index, None)
        self.assertEqual(qc_info.annotation_bed, None)
        self.assertEqual(qc_info.annotation_gtf, None)
        self.assertEqual(qc_info.cellranger_version,None)
        self.assertEqual(qc_info.cellranger_refdata,None)
        self.assertEqual(qc_info.cellranger_probeset,None)

    def test_run_qc_paired_end_fastqs_from_wildcard_explicit_metadata(self):
        """
        run_qc.py: paired-end Fastqs from wildcard with metadata
        """
        # Make directory with Fastqs
        fastq_dir = os.path.join(self.dirn, "fastqs")
        os.mkdir(fastq_dir)
        for fq in ("PJB1_S1_R1_001.fastq.gz",
                   "PJB1_S1_R2_001.fastq.gz",
                   "PJB2_S2_R1_001.fastq.gz",
                   "PJB2_S2_R2_001.fastq.gz"):
            self._mock_fastq(os.path.join(fastq_dir, fq))
        # Implicit output and QC directories
        out_dir = self.dirn
        qc_dir = os.path.join(out_dir, "qc")
        # Run the QC
        self.assertEqual(run_qc([os.path.join(fastq_dir, "*"),
                                 "--organism", "Human"]), 0)
        # Check output and reports
        for f in ("qc",
                  "qc_report.html",
                  f"qc_report.{os.path.basename(out_dir)}.zip",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(out_dir, f)),
                            f"Missing '{f}' under {out_dir}")
        # Check collated Picard insert sizes
        collated_insert_sizes = os.path.join(qc_dir,
                                             "insert_sizes.human.tsv")
        self.assertTrue(os.path.exists(collated_insert_sizes),
                        "Missing collated insert sizes TSV")
        with open(collated_insert_sizes,'rt') as fp:
            self.assertEqual(fp.read(),
                             """#Bam file	Mean insert size	Standard deviation	Median insert size	Median absolute deviation
PJB1_S1_001.bam	153.754829	69.675347	139	37
PJB2_S2_001.bam	153.754829	69.675347	139	37
""")
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(os.path.join(qc_dir, "qc.info"))
        self.assertEqual(qc_info.protocol,"standardPE")
        self.assertEqual(qc_info.protocol_specification,
                         str(fetch_protocol_definition("standardPE")))
        self.assertEqual(qc_info.organism, "Human")
        self.assertEqual(qc_info.seq_data_samples, "PJB1,PJB2")
        self.assertEqual(qc_info.fastq_dir, None)
        self.assertEqual(qc_info.fastqs,
                         "PJB1_S1_R1_001.fastq.gz,"
                         "PJB1_S1_R2_001.fastq.gz,"
                         "PJB2_S2_R1_001.fastq.gz,"
                         "PJB2_S2_R2_001.fastq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane, False)
        self.assertEqual(qc_info.fastq_screens,
                         "model_organisms,other_organisms,rRNA")
        self.assertEqual(qc_info.star_index,"/data/hg38/star_index")
        self.assertEqual(qc_info.annotation_bed,self.ref_data['hg38']['bed'])
        self.assertEqual(qc_info.annotation_gtf,self.ref_data['hg38']['gtf'])
        self.assertEqual(qc_info.cellranger_version,None)
        self.assertEqual(qc_info.cellranger_refdata,None)
        self.assertEqual(qc_info.cellranger_probeset,None)

    def test_run_qc_pseudo_single_end_fastqs_from_wildcard_no_metadata(self):
        """
        run_qc.py: psuedo single-end Fastqs from wildcard (no metadata)
        """
        # Make directory with Fastqs
        fastq_dir = os.path.join(self.dirn, "fastqs")
        os.mkdir(fastq_dir)
        for fq in ("PJB1_S1_R1_001.fastq.gz",
                   "PJB1_S1_R2_001.fastq.gz",
                   "PJB2_S2_R1_001.fastq.gz",
                   "PJB2_S2_R2_001.fastq.gz"):
            self._mock_fastq(os.path.join(fastq_dir, fq))
        # Implicit output and QC directories
        out_dir = self.dirn
        qc_dir = os.path.join(out_dir, "qc")
        # Run the QC
        self.assertEqual(run_qc([os.path.join(fastq_dir, "*_R1_*"),]), 0)
        # Check output and reports
        for f in ("qc",
                  "qc_report.html",
                  f"qc_report.{os.path.basename(out_dir)}.zip",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(out_dir, f)),
                            f"Missing '{f}' under {out_dir}")
        # Check collated Picard insert sizes
        collated_insert_sizes = os.path.join(qc_dir,
                                             "insert_sizes.human.tsv")
        self.assertFalse(os.path.exists(collated_insert_sizes),
                        "Collated insert sizes TSV shouldn't be present")
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(os.path.join(qc_dir, "qc.info"))
        self.assertEqual(qc_info.protocol, "minimal")
        self.assertEqual(qc_info.organism, None)
        self.assertEqual(qc_info.seq_data_samples, "PJB1,PJB2")
        self.assertEqual(qc_info.fastq_dir, None)
        self.assertEqual(qc_info.fastqs,
                         "PJB1_S1_R1_001.fastq.gz,"
                         "PJB2_S2_R1_001.fastq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane, False)
        self.assertEqual(qc_info.fastq_screens,
                         "model_organisms,other_organisms,rRNA")
        self.assertEqual(qc_info.star_index, None)
        self.assertEqual(qc_info.annotation_bed, None)
        self.assertEqual(qc_info.annotation_gtf, None)
        self.assertEqual(qc_info.cellranger_version,None)
        self.assertEqual(qc_info.cellranger_refdata,None)
        self.assertEqual(qc_info.cellranger_probeset,None)

    def test_run_qc_pseudo_single_end_fastqs_from_wildcard_with_metadata(self):
        """
        run_qc.py: psuedo single-end Fastqs from wildcard with metadata
        """
        # Make directory with Fastqs
        fastq_dir = os.path.join(self.dirn, "fastqs")
        os.mkdir(fastq_dir)
        for fq in ("PJB1_S1_R1_001.fastq.gz",
                   "PJB1_S1_R2_001.fastq.gz",
                   "PJB2_S2_R1_001.fastq.gz",
                   "PJB2_S2_R2_001.fastq.gz"):
            self._mock_fastq(os.path.join(fastq_dir, fq))
        # Implicit output and QC directories
        out_dir = self.dirn
        qc_dir = os.path.join(out_dir, "qc")
        # Run the QC
        self.assertEqual(run_qc([os.path.join(fastq_dir, "*_R1_*"),
                                 "--organism", "Human"]), 0)
        # Check output and reports
        for f in ("qc",
                  "qc_report.html",
                  f"qc_report.{os.path.basename(out_dir)}.zip",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(out_dir, f)),
                            f"Missing '{f}' under {out_dir}")
        # Check collated Picard insert sizes
        collated_insert_sizes = os.path.join(qc_dir,
                                             "insert_sizes.human.tsv")
        self.assertFalse(os.path.exists(collated_insert_sizes),
                        "Collated insert sizes TSV shouldn't be present")
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(os.path.join(qc_dir, "qc.info"))
        self.assertEqual(qc_info.protocol, "standardSE")
        self.assertEqual(qc_info.protocol_specification,
                         str(fetch_protocol_definition("standardSE")))
        self.assertEqual(qc_info.organism, "Human")
        self.assertEqual(qc_info.seq_data_samples, "PJB1,PJB2")
        self.assertEqual(qc_info.fastq_dir, None)
        self.assertEqual(qc_info.fastqs,
                         "PJB1_S1_R1_001.fastq.gz,"
                         "PJB2_S2_R1_001.fastq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane, False)
        self.assertEqual(qc_info.fastq_screens,
                         "model_organisms,other_organisms,rRNA")
        self.assertEqual(qc_info.star_index,"/data/hg38/star_index")
        self.assertEqual(qc_info.annotation_bed,self.ref_data['hg38']['bed'])
        self.assertEqual(qc_info.annotation_gtf,self.ref_data['hg38']['gtf'])
        self.assertEqual(qc_info.cellranger_version,None)
        self.assertEqual(qc_info.cellranger_refdata,None)
        self.assertEqual(qc_info.cellranger_probeset,None)

    def test_run_qc_paired_end_fastqs_have_fq_extension(self):
        """
        run_qc.py: paired-end Fastqs with '.fq' extension
        """
        # Make directory with Fastqs
        fastq_dir = os.path.join(self.dirn, "fastqs")
        os.mkdir(fastq_dir)
        for fq in ("PJB1_S1_R1_001.fq.gz",
                   "PJB1_S1_R2_001.fq.gz",
                   "PJB2_S2_R1_001.fq.gz",
                   "PJB2_S2_R2_001.fq.gz"):
            self._mock_fastq(os.path.join(fastq_dir, fq))
        # Implicit output and QC directories
        out_dir = self.dirn
        qc_dir = os.path.join(out_dir, "qc")
        # Run the QC
        self.assertEqual(run_qc([os.path.join(fastq_dir, "*"),
                                 "--organism", "Human"]), 0)
        # Check output and reports
        for f in ("qc",
                  "qc_report.html",
                  f"qc_report.{os.path.basename(out_dir)}.zip",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(out_dir, f)),
                            f"Missing '{f}' under {out_dir}")
        # Check collated Picard insert sizes
        collated_insert_sizes = os.path.join(qc_dir,
                                             "insert_sizes.human.tsv")
        self.assertTrue(os.path.exists(collated_insert_sizes),
                        "Missing collated insert sizes TSV")
        with open(collated_insert_sizes,'rt') as fp:
            self.assertEqual(fp.read(),
                             """#Bam file	Mean insert size	Standard deviation	Median insert size	Median absolute deviation
PJB1_S1_001.bam	153.754829	69.675347	139	37
PJB2_S2_001.bam	153.754829	69.675347	139	37
""")
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(os.path.join(qc_dir, "qc.info"))
        self.assertEqual(qc_info.protocol,"standardPE")
        self.assertEqual(qc_info.protocol_specification,
                         str(fetch_protocol_definition("standardPE")))
        self.assertEqual(qc_info.organism,"Human")
        self.assertEqual(qc_info.seq_data_samples,"PJB1,PJB2")
        self.assertEqual(qc_info.fastq_dir, None)
        self.assertEqual(qc_info.fastqs,
                         "PJB1_S1_R1_001.fq.gz,"
                         "PJB1_S1_R2_001.fq.gz,"
                         "PJB2_S2_R1_001.fq.gz,"
                         "PJB2_S2_R2_001.fq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane,False)
        self.assertEqual(qc_info.fastq_screens,
                         "model_organisms,other_organisms,rRNA")
        self.assertEqual(qc_info.star_index,"/data/hg38/star_index")
        self.assertEqual(qc_info.annotation_bed,self.ref_data['hg38']['bed'])
        self.assertEqual(qc_info.annotation_gtf,self.ref_data['hg38']['gtf'])
        self.assertEqual(qc_info.cellranger_version,None)
        self.assertEqual(qc_info.cellranger_refdata,None)
        self.assertEqual(qc_info.cellranger_probeset,None)