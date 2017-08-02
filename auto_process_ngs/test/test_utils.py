#######################################################################
# Tests for utils.py module
#######################################################################

import unittest
import tempfile
import shutil
import zipfile
from bcftbx.JobRunner import SimpleJobRunner,GEJobRunner
from auto_process_ngs.mock import MockAnalysisDirFactory
from auto_process_ngs.mock import MockAnalysisProject
from auto_process_ngs.utils import *

class TestAnalysisFastq(unittest.TestCase):
    """Tests for the AnalysisFastq class

    """

    def test_full_name(self):
        """Handle full Illumina-style fastq name
        """
        fq = AnalysisFastq('NH1_ChIP-seq_Gli1_ACAGTG_L003_R2_001')
        self.assertEqual(fq.sample_name,'NH1_ChIP-seq_Gli1')
        self.assertEqual(fq.basename,
                         'NH1_ChIP-seq_Gli1_ACAGTG_L003_R2_001')
        self.assertEqual(fq.extension,'')
        self.assertEqual(fq.sample_number,None)
        self.assertEqual(fq.barcode_sequence,'ACAGTG')
        self.assertEqual(fq.lane_number,3)
        self.assertEqual(fq.read_number,2)
        self.assertEqual(fq.set_number,1)
        self.assertEqual(str(fq),'NH1_ChIP-seq_Gli1_ACAGTG_L003_R2_001')

    def test_full_name_dual_index(self):
        """Handle full Illumina-style fastq name with dual index
        """
        fq = AnalysisFastq('NH1_ChIP-seq_Gli1_ACAGTG-GTTCAC_L003_R2_001')
        self.assertEqual(fq.sample_name,'NH1_ChIP-seq_Gli1')
        self.assertEqual(fq.basename,
                         'NH1_ChIP-seq_Gli1_ACAGTG-GTTCAC_L003_R2_001')
        self.assertEqual(fq.extension,'')
        self.assertEqual(fq.sample_number,None)
        self.assertEqual(fq.barcode_sequence,'ACAGTG-GTTCAC')
        self.assertEqual(fq.lane_number,3)
        self.assertEqual(fq.read_number,2)
        self.assertEqual(fq.set_number,1)
        self.assertEqual(str(fq),'NH1_ChIP-seq_Gli1_ACAGTG-GTTCAC_L003_R2_001')

    def test_full_name_blc2fastq2(self):
        """Handle Illumina fastq name from bcl2fastq2
        """
        fq = AnalysisFastq('NH1_ChIP-seq_Gli1_S4_L003_R2_001')
        self.assertEqual(fq.sample_name,'NH1_ChIP-seq_Gli1')
        self.assertEqual(fq.basename,
                         'NH1_ChIP-seq_Gli1_S4_L003_R2_001')
        self.assertEqual(fq.extension,'')
        self.assertEqual(fq.sample_number,4)
        self.assertEqual(fq.barcode_sequence,None)
        self.assertEqual(fq.lane_number,3)
        self.assertEqual(fq.read_number,2)
        self.assertEqual(fq.set_number,1)
        self.assertEqual(str(fq),'NH1_ChIP-seq_Gli1_S4_L003_R2_001')

    def test_name_no_lane_blc2fastq2(self):
        """Handle Illumina fastq name from bcl2fastq2 (without lane)
        """
        fq = AnalysisFastq('NH1_ChIP-seq_Gli1_S4_R2_001')
        self.assertEqual(fq.sample_name,'NH1_ChIP-seq_Gli1')
        self.assertEqual(fq.basename,
                         'NH1_ChIP-seq_Gli1_S4_R2_001')
        self.assertEqual(fq.extension,'')
        self.assertEqual(fq.sample_number,4)
        self.assertEqual(fq.barcode_sequence,None)
        self.assertEqual(fq.lane_number,None)
        self.assertEqual(fq.read_number,2)
        self.assertEqual(fq.set_number,1)
        self.assertEqual(str(fq),'NH1_ChIP-seq_Gli1_S4_R2_001')

    def test_name_only(self):
        """Handle reduced fastq name (sample name only)
        """
        fq = AnalysisFastq('NH1_ChIP-seq_Gli1')
        self.assertEqual(fq.sample_name,'NH1_ChIP-seq_Gli1')
        self.assertEqual(fq.basename,'NH1_ChIP-seq_Gli1')
        self.assertEqual(fq.extension,'')
        self.assertEqual(fq.sample_number,None)
        self.assertEqual(fq.barcode_sequence,None)
        self.assertEqual(fq.lane_number,None)
        self.assertEqual(fq.read_number,None)
        self.assertEqual(fq.set_number,None)
        self.assertEqual(str(fq),'NH1_ChIP-seq_Gli1')

    def test_name_only_paired_end(self):
        """Handle reduced fastq name (sample name only, paired end)
        """
        fq = AnalysisFastq('NH1_ChIP-seq_Gli1_R2')
        self.assertEqual(fq.sample_name,'NH1_ChIP-seq_Gli1')
        self.assertEqual(fq.basename,'NH1_ChIP-seq_Gli1_R2')
        self.assertEqual(fq.extension,'')
        self.assertEqual(fq.sample_number,None)
        self.assertEqual(fq.barcode_sequence,None)
        self.assertEqual(fq.lane_number,None)
        self.assertEqual(fq.read_number,2)
        self.assertEqual(fq.set_number,None)
        self.assertEqual(str(fq),'NH1_ChIP-seq_Gli1_R2')

    def test_name_and_lane(self):
        """Handle reduced fastq name (sample name and lane)
        """
        fq = AnalysisFastq('NH1_ChIP-seq_Gli1_L001')
        self.assertEqual(fq.sample_name,'NH1_ChIP-seq_Gli1')
        self.assertEqual(fq.basename,'NH1_ChIP-seq_Gli1_L001')
        self.assertEqual(fq.extension,'')
        self.assertEqual(fq.sample_number,None)
        self.assertEqual(fq.barcode_sequence,None)
        self.assertEqual(fq.lane_number,1)
        self.assertEqual(fq.read_number,None)
        self.assertEqual(fq.set_number,None)
        self.assertEqual(str(fq),'NH1_ChIP-seq_Gli1_L001')

    def test_name_and_lane_paired_end(self):
        """Handle reduced fastq name (sample name and lane, paired end)
        """
        fq = AnalysisFastq('NH1_ChIP-seq_Gli1_L001_R2')
        self.assertEqual(fq.sample_name,'NH1_ChIP-seq_Gli1')
        self.assertEqual(fq.basename,'NH1_ChIP-seq_Gli1_L001_R2')
        self.assertEqual(fq.extension,'')
        self.assertEqual(fq.sample_number,None)
        self.assertEqual(fq.barcode_sequence,None)
        self.assertEqual(fq.lane_number,1)
        self.assertEqual(fq.read_number,2)
        self.assertEqual(fq.set_number,None)
        self.assertEqual(str(fq),'NH1_ChIP-seq_Gli1_L001_R2')

    def test_name_and_tag(self):
        """Handle reduced fastq name (sample name and barcode)
        """
        fq = AnalysisFastq('NH1_ChIP-seq_Gli1_ACAGTG')
        self.assertEqual(fq.sample_name,'NH1_ChIP-seq_Gli1')
        self.assertEqual(fq.basename,'NH1_ChIP-seq_Gli1_ACAGTG')
        self.assertEqual(fq.extension,'')
        self.assertEqual(fq.sample_number,None)
        self.assertEqual(fq.barcode_sequence,'ACAGTG')
        self.assertEqual(fq.lane_number,None)
        self.assertEqual(fq.read_number,None)
        self.assertEqual(fq.set_number,None)
        self.assertEqual(str(fq),'NH1_ChIP-seq_Gli1_ACAGTG')

    def test_name_and_tag_paired_end(self):
        """Handle reduced fastq name (sample name and barcode, paired end)
        """
        fq = AnalysisFastq('NH1_ChIP-seq_Gli1_ACAGTG_R2')
        self.assertEqual(fq.sample_name,'NH1_ChIP-seq_Gli1')
        self.assertEqual(fq.basename,'NH1_ChIP-seq_Gli1_ACAGTG_R2')
        self.assertEqual(fq.extension,'')
        self.assertEqual(fq.sample_number,None)
        self.assertEqual(fq.barcode_sequence,'ACAGTG')
        self.assertEqual(fq.lane_number,None)
        self.assertEqual(fq.read_number,2)
        self.assertEqual(fq.set_number,None)
        self.assertEqual(str(fq),'NH1_ChIP-seq_Gli1_ACAGTG_R2')

    def test_name_tag_and_lane(self):
        """Handle reduced fastq name (sample name, barcode and lane)
        """
        fq = AnalysisFastq('NH1_ChIP-seq_Gli1_ACAGTG_L001')
        self.assertEqual(fq.sample_name,'NH1_ChIP-seq_Gli1')
        self.assertEqual(fq.basename,'NH1_ChIP-seq_Gli1_ACAGTG_L001')
        self.assertEqual(fq.extension,'')
        self.assertEqual(fq.sample_number,None)
        self.assertEqual(fq.barcode_sequence,'ACAGTG')
        self.assertEqual(fq.lane_number,1)
        self.assertEqual(fq.read_number,None)
        self.assertEqual(fq.set_number,None)
        self.assertEqual(str(fq),'NH1_ChIP-seq_Gli1_ACAGTG_L001')

    def test_name_tag_and_lane_paired_end(self):
        """Handle reduced fastq name (sample name, barcode and lane, paired end)
        """
        fq = AnalysisFastq('NH1_ChIP-seq_Gli1_ACAGTG_L001_R2')
        self.assertEqual(fq.sample_name,'NH1_ChIP-seq_Gli1')
        self.assertEqual(fq.basename,'NH1_ChIP-seq_Gli1_ACAGTG_L001_R2')
        self.assertEqual(fq.extension,'')
        self.assertEqual(fq.sample_number,None)
        self.assertEqual(fq.barcode_sequence,'ACAGTG')
        self.assertEqual(fq.lane_number,1)
        self.assertEqual(fq.read_number,2)
        self.assertEqual(fq.set_number,None)
        self.assertEqual(str(fq),'NH1_ChIP-seq_Gli1_ACAGTG_L001_R2')

    def test_AGTC_sample_names(self):
        """Handle sample names consisting of letters 'A', 'G', 'T' and 'C'
        """
        for name in ('A','G','T','C','AGCT'):
            fq = AnalysisFastq('%s_R1' % name)
            self.assertEqual(fq.sample_name,name)
            self.assertEqual(fq.sample_number,None)
            self.assertEqual(fq.barcode_sequence,None)
            self.assertEqual(fq.lane_number,None)
            self.assertEqual(fq.read_number,1)
            self.assertEqual(fq.set_number,None)
            self.assertEqual(str(fq),'%s_R1' % name)

    def test_non_standard_sample_name(self):
        """Handle non-standard Fastq names with sample name only
        """
        fq = AnalysisFastq('NH1_ChIP-seq.r2')
        self.assertEqual(fq.sample_name,'NH1_ChIP-seq')
        self.assertEqual(fq.basename,'NH1_ChIP-seq.r2')
        self.assertEqual(fq.extension,'')
        self.assertEqual(fq.sample_number,None)
        self.assertEqual(fq.barcode_sequence,None)
        self.assertEqual(fq.lane_number,None)
        self.assertEqual(fq.read_number,2)
        self.assertEqual(fq.set_number,None)
        self.assertEqual(str(fq),'NH1_ChIP-seq.r2')

    def test_non_standard_sample_name_with_dots(self):
        """Handle non-standard Fastq names with sample name containing dots
        """
        fq = AnalysisFastq('NH1.2.r2')
        self.assertEqual(fq.sample_name,'NH1.2')
        self.assertEqual(fq.basename,'NH1.2.r2')
        self.assertEqual(fq.extension,'')
        self.assertEqual(fq.sample_number,None)
        self.assertEqual(fq.barcode_sequence,None)
        self.assertEqual(fq.lane_number,None)
        self.assertEqual(fq.read_number,2)
        self.assertEqual(fq.set_number,None)
        self.assertEqual(str(fq),'NH1.2.r2')

    def test_non_standard_sample_name_and_barcode(self):
        """Handle non-standard Fastq names with sample name and barcode
        """
        fq = AnalysisFastq('NH1_ChIP-seq.ACAGTG.r2')
        self.assertEqual(fq.sample_name,'NH1_ChIP-seq')
        self.assertEqual(fq.basename,'NH1_ChIP-seq.ACAGTG.r2')
        self.assertEqual(fq.extension,'')
        self.assertEqual(fq.sample_number,None)
        self.assertEqual(fq.barcode_sequence,'ACAGTG')
        self.assertEqual(fq.lane_number,None)
        self.assertEqual(fq.read_number,2)
        self.assertEqual(fq.set_number,None)
        self.assertEqual(str(fq),'NH1_ChIP-seq.ACAGTG.r2')

    def test_input_is_full_path(self):
        """Handle input as full path to Fastq file
        """
        fq = AnalysisFastq('/data/Project_NH/Sample_NH1/NH1_ChIP-seq_Gli1_ACAGTG_L003_R2_001.fastq.gz')
        self.assertEqual(fq.sample_name,'NH1_ChIP-seq_Gli1')
        self.assertEqual(fq.basename,'NH1_ChIP-seq_Gli1_ACAGTG_L003_R2_001')
        self.assertEqual(fq.extension,'.fastq.gz')
        self.assertEqual(fq.sample_number,None)
        self.assertEqual(fq.barcode_sequence,'ACAGTG')
        self.assertEqual(fq.lane_number,3)
        self.assertEqual(fq.read_number,2)
        self.assertEqual(fq.set_number,1)
        self.assertEqual(str(fq),'NH1_ChIP-seq_Gli1_ACAGTG_L003_R2_001')

class TestAnalysisDir(unittest.TestCase):
    """Tests for the AnalysisDir class

    """
    def setUp(self):
        self.dirn = tempfile.mkdtemp(suffix='TestAnalysisDir')

    def tearDown(self):
        # Remove the temporary test directory
        shutil.rmtree(self.dirn)

    def test_casava(self):
        """Check AnalysisDir against CASAVA-style output
        """
        mockdir = MockAnalysisDirFactory.casava(
            '160621_M00879_0087_000000000-AGEW9',
            'miseq',
            top_dir=self.dirn)
        mockdir.create()
        analysis_dir = AnalysisDir(mockdir.dirn)
        self.assertEqual(analysis_dir.run_name,mockdir.run_name)
        self.assertEqual(analysis_dir.n_sequencing_data,1)
        self.assertEqual(analysis_dir.n_projects,2)
        self.assertTrue(analysis_dir.paired_end)

    def test_bcl2fastq2(self):
        """Check AnalysisDir against bcl2fastq v2-style output
        """
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_M00879_0087_000000000-AGEW9',
            'miseq',
            top_dir=self.dirn)
        mockdir.create()
        analysis_dir = AnalysisDir(mockdir.dirn)
        self.assertEqual(analysis_dir.run_name,mockdir.run_name)
        self.assertEqual(analysis_dir.n_sequencing_data,1)
        self.assertEqual(analysis_dir.n_projects,2)
        self.assertTrue(analysis_dir.paired_end)

    def test_handle_non_project_dir(self):
        """Check AnalysisDir with non-project directory
        """
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_M00879_0087_000000000-AGEW9',
            'miseq',
            top_dir=self.dirn)
        mockdir.create()
        # Add non-project dir
        non_project_dir = os.path.join(self.dirn,
                                '160621_M00879_0087_000000000-AGEW9_analysis',
                                'extras')
        fqs_dir = os.path.join(non_project_dir,'fastqs')
        os.mkdir(non_project_dir)
        os.mkdir(fqs_dir)
        for fq in ('PB04_S4_R1_unpaired.fastq.gz',
                   'PB04_trimmoPE_bowtie2_notHg38.1.fastq.gz'):
            with open(os.path.join(fqs_dir,fq),'w') as fp:
                fp.write("")
        # Load and check the analysis dir
        analysis_dir = AnalysisDir(mockdir.dirn)
        self.assertEqual(analysis_dir.run_name,mockdir.run_name)
        self.assertEqual(analysis_dir.n_sequencing_data,1)
        self.assertEqual(analysis_dir.n_projects,2)
        self.assertTrue(analysis_dir.paired_end)

class TestAnalysisProject(unittest.TestCase):
    """Tests for the AnalysisProject class

    """
    def setUp(self):
        # Create a temporary directory for tests
        self.dirn = tempfile.mkdtemp(suffix='TestAnalysisProject')

    def make_data_dir(self,fastq_list):
        # Make a fake data source directory
        self.fastqs = []
        fake_fastqs_dir = os.path.join(self.dirn,'fake_fastqs')
        os.mkdir(fake_fastqs_dir)
        for fq in fastq_list:
            fastq = os.path.join(fake_fastqs_dir,fq)
            open(fastq,'w').close()
            self.fastqs.append(fastq)

    def make_mock_project_dir(self,name,fastq_list,fastq_dir='fastqs',
                              primary_fastq_dir=None):
        # Make a mock project directory
        if primary_fastq_dir is None:
            primary_fastq_dir = fastq_dir
        metadata = { 'Primary fastqs': primary_fastq_dir }
        MockAnalysisProject(name,
                            fastq_list,
                            fastq_dir=fastq_dir,
                            metadata=metadata).create(top_dir=self.dirn)

    def tearDown(self):
        # Remove the temporary test directory
        shutil.rmtree(self.dirn)

    def test_empty_analysis_project(self):
        """Check empty AnalysisProject class
        """
        dirn = os.path.join(self.dirn,'PJB')
        project = AnalysisProject('PJB',dirn)
        self.assertEqual(project.name,'PJB')
        self.assertEqual(project.dirn,dirn)
        self.assertEqual(project.samples,[])
        self.assertFalse(project.multiple_fastqs)
        self.assertEqual(project.fastq_dir,None)
        self.assertEqual(project.info.library_type,None)
        self.assertEqual(project.info.organism,None)
        self.assertFalse(project.info.paired_end)
        self.assertEqual(project.info.platform,None)
        self.assertEqual(project.info.primary_fastq_dir,None)
        self.assertEqual(project.fastq_dirs,[])

    def test_create_single_end_analysis_project(self):
        """Check creation of new single-end AnalysisProject directory
        """
        self.make_data_dir(('PJB1-A_ACAGTG_L001_R1_001.fastq.gz',
                            'PJB1-B_ACAGTG_L002_R1_001.fastq.gz',))
        dirn = os.path.join(self.dirn,'PJB')
        project = AnalysisProject('PJB',dirn)
        project.create_directory(fastqs=self.fastqs)
        self.assertEqual(project.name,'PJB')
        self.assertTrue(os.path.isdir(project.dirn))
        self.assertFalse(project.multiple_fastqs)
        self.assertFalse(project.info.paired_end)
        self.assertEqual(project.info.primary_fastq_dir,'fastqs')
        self.assertEqual(project.samples[0].name,'PJB1-A')
        self.assertEqual(project.samples[1].name,'PJB1-B')
        self.assertEqual(project.fastq_dir,
                         os.path.join(project.dirn,'fastqs'))
        self.assertEqual(project.fastq_dirs,['fastqs',])

    def test_create_single_end_analysis_project_multi_fastqs(self):
        """Check creation of new single-end AnalysisProject directory (multi-fastq/sample)
        """
        self.make_data_dir(('PJB1-A_ACAGTG_L001_R1_001.fastq.gz',
                            'PJB1-A_ACAGTG_L002_R1_001.fastq.gz',
                            'PJB1-B_ACAGTG_L001_R1_001.fastq.gz',
                            'PJB1-B_ACAGTG_L002_R1_001.fastq.gz',))
        dirn = os.path.join(self.dirn,'PJB')
        project = AnalysisProject('PJB',dirn)
        project.create_directory(fastqs=self.fastqs)
        self.assertEqual(project.name,'PJB')
        self.assertTrue(os.path.isdir(project.dirn))
        self.assertTrue(project.multiple_fastqs)
        self.assertFalse(project.info.paired_end)
        self.assertEqual(project.info.primary_fastq_dir,'fastqs')
        self.assertEqual(project.samples[0].name,'PJB1-A')
        self.assertEqual(project.samples[1].name,'PJB1-B')
        self.assertEqual(project.fastq_dir,
                         os.path.join(project.dirn,'fastqs'))
        self.assertEqual(project.fastq_dirs,['fastqs',])

    def test_create_paired_end_analysis_project(self):
        """Check creation of new paired-end AnalysisProject directory
        """
        self.make_data_dir(('PJB1-A_ACAGTG_L001_R1_001.fastq.gz',
                            'PJB1-B_ACAGTG_L002_R1_001.fastq.gz',
                            'PJB1-A_ACAGTG_L001_R2_001.fastq.gz',
                            'PJB1-B_ACAGTG_L002_R2_001.fastq.gz',))
        dirn = os.path.join(self.dirn,'PJB')
        project = AnalysisProject('PJB',dirn,fastq_dir='fastqs.test')
        project.create_directory(fastqs=self.fastqs)
        self.assertEqual(project.name,'PJB')
        self.assertTrue(os.path.isdir(project.dirn))
        self.assertFalse(project.multiple_fastqs)
        self.assertTrue(project.info.paired_end)
        self.assertEqual(project.info.primary_fastq_dir,'fastqs')
        self.assertEqual(project.samples[0].name,'PJB1-A')
        self.assertEqual(project.samples[1].name,'PJB1-B')
        self.assertEqual(project.samples[1].name,'PJB1-B')
        self.assertEqual(project.fastq_dir,
                         os.path.join(project.dirn,'fastqs'))
        self.assertEqual(project.fastq_dirs,['fastqs',])

    def test_create_paired_end_analysis_project_multi_fastqs(self):
        """Check creation of new paired-end AnalysisProject directory (multi-fastq/sample)
        """
        self.make_data_dir(('PJB1-B_ACAGTG_L001_R1_001.fastq.gz',
                            'PJB1-B_ACAGTG_L002_R1_001.fastq.gz',
                            'PJB1-B_ACAGTG_L001_R2_001.fastq.gz',
                            'PJB1-B_ACAGTG_L002_R2_001.fastq.gz',))
        dirn = os.path.join(self.dirn,'PJB')
        project = AnalysisProject('PJB',dirn)
        project.create_directory(fastqs=self.fastqs)
        self.assertEqual(project.name,'PJB')
        self.assertTrue(os.path.isdir(project.dirn))
        self.assertEqual(project.samples[0].name,'PJB1-B')
        self.assertTrue(project.multiple_fastqs)
        self.assertTrue(project.info.paired_end)
        self.assertEqual(project.fastq_dir,
                         os.path.join(project.dirn,'fastqs'))
        self.assertEqual(project.fastq_dirs,['fastqs',])
        self.assertEqual(project.info.primary_fastq_dir,'fastqs')

    def test_create_analysis_project_not_standard_fastq_dir(self):
        """Check creation of AnalysisProject directory with non-standard fastq dir
        """
        self.make_data_dir(('PJB1-B_ACAGTG_L001_R1_001.fastq.gz',
                            'PJB1-B_ACAGTG_L002_R1_001.fastq.gz',
                            'PJB1-B_ACAGTG_L001_R2_001.fastq.gz',
                            'PJB1-B_ACAGTG_L002_R2_001.fastq.gz',))
        dirn = os.path.join(self.dirn,'PJB')
        project = AnalysisProject('PJB',dirn)
        project.create_directory(fastqs=self.fastqs,fastq_dir='fastqs.test')
        self.assertEqual(project.name,'PJB')
        self.assertTrue(os.path.isdir(project.dirn))
        self.assertEqual(project.samples[0].name,'PJB1-B')
        self.assertTrue(project.multiple_fastqs)
        self.assertTrue(project.info.paired_end)
        self.assertEqual(project.fastq_dir,
                         os.path.join(project.dirn,'fastqs.test'))
        self.assertEqual(project.fastq_dirs,['fastqs.test',])
        self.assertEqual(project.info.primary_fastq_dir,'fastqs.test')

    def test_load_single_end_analysis_project(self):
        """Check loading of an existing single-end AnalysisProject directory
        """
        self.make_mock_project_dir(
            'PJB',
            ('PJB1-A_ACAGTG_L001_R1_001.fastq.gz',
             'PJB1-B_ACAGTG_L002_R1_001.fastq.gz',))
        dirn = os.path.join(self.dirn,'PJB')
        project = AnalysisProject('PJB',dirn)
        self.assertEqual(project.name,'PJB')
        self.assertTrue(os.path.isdir(project.dirn))
        self.assertFalse(project.multiple_fastqs)
        self.assertFalse(project.info.paired_end)
        self.assertEqual(project.samples[0].name,'PJB1-A')
        self.assertEqual(project.samples[1].name,'PJB1-B')
        self.assertEqual(project.fastq_dir,
                         os.path.join(project.dirn,'fastqs'))
        self.assertEqual(project.fastq_dirs,['fastqs',])
        self.assertEqual(project.info.primary_fastq_dir,'fastqs')

    def test_load_analysis_project_non_canonical_fastq_dir(self):
        """Check AnalysisProject loading for directory with non-canonical fastq directory
        """
        self.make_mock_project_dir(
            'PJB',
            ('PJB1-A_ACAGTG_L001_R1_001.fastq.gz',
             'PJB1-B_ACAGTG_L002_R1_001.fastq.gz',),
            fastq_dir='fastqs.test')
        dirn = os.path.join(self.dirn,'PJB')
        project = AnalysisProject('PJB',dirn,fastq_dir='fastqs.test')
        self.assertEqual(project.name,'PJB')
        self.assertTrue(os.path.isdir(project.dirn))
        self.assertFalse(project.multiple_fastqs)
        self.assertFalse(project.info.paired_end)
        self.assertEqual(project.samples[0].name,'PJB1-A')
        self.assertEqual(project.samples[1].name,'PJB1-B')
        self.assertEqual(project.fastq_dir,
                         os.path.join(project.dirn,'fastqs.test'))
        self.assertEqual(project.fastq_dirs,['fastqs.test',])
        self.assertEqual(project.info.primary_fastq_dir,'fastqs.test')

    def test_load_analysis_project_non_canonical_fastqs(self):
        """Check AnalysisProject loading fails for directory with non-canonical fastqs
        """
        self.make_mock_project_dir(
            'PJB',
            ('PB04_S4_R1_unpaired.fastq.gz',
             'PB04_trimmoPE_bowtie2_notHg38.1.fastq.gz',))
        dirn = os.path.join(self.dirn,'PJB')
        project = AnalysisProject('PJB',dirn)
        self.assertEqual(project.name,'PJB')
        self.assertTrue(os.path.isdir(project.dirn))
        self.assertFalse(project.multiple_fastqs)
        self.assertFalse(project.info.paired_end)
        self.assertEqual(project.fastq_dir,
                         os.path.join(project.dirn,'fastqs'))
        self.assertEqual(project.fastq_dirs,['fastqs',])
        self.assertEqual(project.info.primary_fastq_dir,'fastqs')

    def test_load_analysis_project_non_with_alternative_fastq_naming(self):
        """Check AnalysisProject loading for directory with alternative fastq naming
        """
        self.make_mock_project_dir(
            'PJB',
            ('PB02.trimmed.filtered.r1.fastq.gz',
             'PB02.trimmed.filtered.r2.fastq.gz',))
        # Create a class to handle the non-canonical Fastq names
        class NonCanonicalFastq(BaseFastqAttrs):
            def __init__(self,fastq):
                BaseFastqAttrs.__init__(self,fastq)
                self.sample_name = self.basename.split('.')[0]
                self.read_number = int(self.basename.split('.')[-1][1:])
        dirn = os.path.join(self.dirn,'PJB')
        project = AnalysisProject('PJB',dirn,fastq_attrs=NonCanonicalFastq)
        self.assertEqual(project.name,'PJB')
        self.assertTrue(os.path.isdir(project.dirn))
        self.assertFalse(project.multiple_fastqs)
        self.assertTrue(project.info.paired_end)
        self.assertEqual(project.samples[0].name,'PB02')
        self.assertEqual(project.fastq_dir,
                         os.path.join(project.dirn,'fastqs'))
        self.assertEqual(project.fastq_dirs,['fastqs',])
        self.assertEqual(project.info.primary_fastq_dir,'fastqs')

    def test_load_analysis_project_detect_multiple_fastq_dirs(self):
        """Check AnalysisProject detects multiple fastqs directories
        """
        # Construct test project with two fastq subdirectories
        self.make_mock_project_dir(
            'PJB',
            ('PJB1-A_ACAGTG_L001_R1_001.fastq.gz',
             'PJB1-B_ACAGTG_L002_R1_001.fastq.gz',))
        self.make_mock_project_dir(
            'PJB.untrimmed',
            ('PJB1-A-untrimmed_ACAGTG_L001_R1_001.fastq.gz',
             'PJB1-B-untrimmed_ACAGTG_L002_R1_001.fastq.gz',),
            fastq_dir='fastqs.untrimmed')
        shutil.move(os.path.join(self.dirn,
                                 'PJB.untrimmed',
                                 'fastqs.untrimmed'),
                    os.path.join(self.dirn,'PJB'))
        shutil.rmtree(os.path.join(self.dirn,'PJB.untrimmed'))
        # Load and check AnalysisProject: default fastqs dir
        dirn = os.path.join(self.dirn,'PJB')
        project = AnalysisProject('PJB',dirn)
        self.assertEqual(project.name,'PJB')
        self.assertTrue(os.path.isdir(project.dirn))
        self.assertFalse(project.multiple_fastqs)
        self.assertFalse(project.info.paired_end)
        self.assertEqual(project.samples[0].name,'PJB1-A')
        self.assertEqual(project.samples[1].name,'PJB1-B')
        self.assertEqual(project.fastq_dir,
                         os.path.join(project.dirn,'fastqs'))
        self.assertEqual(project.fastq_dirs,
                         ['fastqs','fastqs.untrimmed'])
        # Load and check AnalysisProject: default fastqs dir
        dirn = os.path.join(self.dirn,'PJB')
        project = AnalysisProject('PJB',dirn,
                                  fastq_dir='fastqs.untrimmed')
        self.assertEqual(project.name,'PJB')
        self.assertTrue(os.path.isdir(project.dirn))
        self.assertFalse(project.multiple_fastqs)
        self.assertFalse(project.info.paired_end)
        self.assertEqual(project.samples[0].name,'PJB1-A-untrimmed')
        self.assertEqual(project.samples[1].name,'PJB1-B-untrimmed')
        self.assertEqual(project.fastq_dir,
                         os.path.join(project.dirn,'fastqs.untrimmed'))
        self.assertEqual(project.fastq_dirs,
                         ['fastqs','fastqs.untrimmed'])
        self.assertEqual(project.info.primary_fastq_dir,'fastqs')

    def test_analysis_project_switch_fastq_dir(self):
        """Check AnalysisProject can switch between multiple fastqs directories
        """
        # Construct test project with two fastq subdirectories
        self.make_mock_project_dir(
            'PJB',
            ('PJB1-A_ACAGTG_L001_R1_001.fastq.gz',
             'PJB1-B_ACAGTG_L002_R1_001.fastq.gz',))
        self.make_mock_project_dir(
            'PJB.untrimmed',
            ('PJB1-A-untrimmed_ACAGTG_L001_R1_001.fastq.gz',
             'PJB1-B-untrimmed_ACAGTG_L002_R1_001.fastq.gz',),
            fastq_dir='fastqs.untrimmed')
        shutil.move(os.path.join(self.dirn,
                                 'PJB.untrimmed',
                                 'fastqs.untrimmed'),
                    os.path.join(self.dirn,'PJB'))
        shutil.rmtree(os.path.join(self.dirn,'PJB.untrimmed'))
        # Load and check AnalysisProject: default fastqs dir
        dirn = os.path.join(self.dirn,'PJB')
        project = AnalysisProject('PJB',dirn)
        self.assertEqual(project.name,'PJB')
        self.assertTrue(os.path.isdir(project.dirn))
        self.assertFalse(project.multiple_fastqs)
        self.assertFalse(project.info.paired_end)
        self.assertEqual(project.samples[0].name,'PJB1-A')
        self.assertEqual(project.samples[1].name,'PJB1-B')
        self.assertEqual(project.fastq_dir,
                         os.path.join(project.dirn,'fastqs'))
        self.assertEqual(project.fastq_dirs,
                         ['fastqs','fastqs.untrimmed'])
        # Switch to alternative fastqs dir
        project.use_fastq_dir('fastqs.untrimmed')
        self.assertEqual(project.name,'PJB')
        self.assertTrue(os.path.isdir(project.dirn))
        self.assertFalse(project.multiple_fastqs)
        self.assertFalse(project.info.paired_end)
        self.assertEqual(project.samples[0].name,'PJB1-A-untrimmed')
        self.assertEqual(project.samples[1].name,'PJB1-B-untrimmed')
        self.assertEqual(project.fastq_dir,
                         os.path.join(project.dirn,'fastqs.untrimmed'))
        self.assertEqual(project.fastq_dirs,
                         ['fastqs','fastqs.untrimmed'])
        self.assertEqual(project.info.primary_fastq_dir,'fastqs')

    def test_analysis_project_switch_to_default_fastq_dir(self):
        """Check AnalysisProject switches to default fastq set
        """
        # Construct test project with two fastq subdirectories
        self.make_mock_project_dir(
            'PJB',
            ('PJB1-A_ACAGTG_L001_R1_001.fastq.gz',
             'PJB1-B_ACAGTG_L002_R1_001.fastq.gz',))
        self.make_mock_project_dir(
            'PJB.untrimmed',
            ('PJB1-A-untrimmed_ACAGTG_L001_R1_001.fastq.gz',
             'PJB1-B-untrimmed_ACAGTG_L002_R1_001.fastq.gz',),
            fastq_dir='fastqs.untrimmed')
        shutil.move(os.path.join(self.dirn,
                                 'PJB.untrimmed',
                                 'fastqs.untrimmed'),
                    os.path.join(self.dirn,'PJB'))
        shutil.rmtree(os.path.join(self.dirn,'PJB.untrimmed'))
        # Load and check AnalysisProject with alternative fastq set
        dirn = os.path.join(self.dirn,'PJB')
        project = AnalysisProject('PJB',dirn,fastq_dir='fastqs.untrimmed')
        self.assertEqual(project.name,'PJB')
        self.assertTrue(os.path.isdir(project.dirn))
        self.assertFalse(project.multiple_fastqs)
        self.assertFalse(project.info.paired_end)
        self.assertEqual(project.samples[0].name,'PJB1-A-untrimmed')
        self.assertEqual(project.samples[1].name,'PJB1-B-untrimmed')
        self.assertEqual(project.fastq_dir,
                         os.path.join(project.dirn,'fastqs.untrimmed'))
        self.assertEqual(project.fastq_dirs,
                         ['fastqs','fastqs.untrimmed'])
        self.assertEqual(project.info.primary_fastq_dir,'fastqs')
        # Implicitly switch to primary fastq set
        project.use_fastq_dir()
        self.assertEqual(project.name,'PJB')
        self.assertTrue(os.path.isdir(project.dirn))
        self.assertFalse(project.multiple_fastqs)
        self.assertFalse(project.info.paired_end)
        self.assertEqual(project.samples[0].name,'PJB1-A')
        self.assertEqual(project.samples[1].name,'PJB1-B')
        self.assertEqual(project.fastq_dir,
                         os.path.join(project.dirn,'fastqs'))
        self.assertEqual(project.fastq_dirs,
                         ['fastqs','fastqs.untrimmed'])
        self.assertEqual(project.info.primary_fastq_dir,'fastqs')

    def test_analysis_project_switch_to_default_non_canonical_fastq_dir(self):
        """Check AnalysisProject switches to default (non-canonical) fastq set
        """
        # Construct test project with two fastq subdirectories
        # and make non-canonical named primary set
        self.make_mock_project_dir(
            'PJB',
            ('PJB1-A-untrimmed_ACAGTG_L001_R1_001.fastq.gz',
             'PJB1-B-untrimmed_ACAGTG_L002_R1_001.fastq.gz',),
            fastq_dir='fastqs.untrimmed')
        self.make_mock_project_dir(
            'PJB.trimmed',
            ('PJB1-A_ACAGTG_L001_R1_001.fastq.gz',
             'PJB1-B_ACAGTG_L002_R1_001.fastq.gz',))
        shutil.move(os.path.join(self.dirn,
                                 'PJB.trimmed',
                                 'fastqs'),
                    os.path.join(self.dirn,'PJB'))
        shutil.rmtree(os.path.join(self.dirn,'PJB.trimmed'))
        # Load and check AnalysisProject primary fastq set
        dirn = os.path.join(self.dirn,'PJB')
        project = AnalysisProject('PJB',dirn)
        self.assertEqual(project.info.primary_fastq_dir,'fastqs.untrimmed')
        self.assertEqual(project.fastq_dir,
                         os.path.join(project.dirn,'fastqs.untrimmed'))
        # Load again with alternative fastq set
        project = AnalysisProject('PJB',dirn,fastq_dir='fastqs')
        self.assertEqual(project.info.primary_fastq_dir,'fastqs.untrimmed')
        self.assertEqual(project.fastq_dir,
                         os.path.join(project.dirn,'fastqs'))
        # Implicitly switch to primary fastq set
        project.use_fastq_dir()
        self.assertEqual(project.info.primary_fastq_dir,'fastqs.untrimmed')
        self.assertEqual(project.fastq_dir,
                         os.path.join(project.dirn,'fastqs.untrimmed'))

    def test_analysis_project_update_primary_fastq_dir(self):
        """Check AnalysisProject primary fastq set can be updated
        """
        # Construct test project with two fastq subdirectories
        self.make_mock_project_dir(
            'PJB',
            ('PJB1-A_ACAGTG_L001_R1_001.fastq.gz',
             'PJB1-B_ACAGTG_L002_R1_001.fastq.gz',))
        self.make_mock_project_dir(
            'PJB.untrimmed',
            ('PJB1-A-untrimmed_ACAGTG_L001_R1_001.fastq.gz',
             'PJB1-B-untrimmed_ACAGTG_L002_R1_001.fastq.gz',),
            fastq_dir='fastqs.untrimmed')
        shutil.move(os.path.join(self.dirn,
                                 'PJB.untrimmed',
                                 'fastqs.untrimmed'),
                    os.path.join(self.dirn,'PJB'))
        shutil.rmtree(os.path.join(self.dirn,'PJB.untrimmed'))
        # Load and check AnalysisProject primary fastq set
        dirn = os.path.join(self.dirn,'PJB')
        project = AnalysisProject('PJB',dirn)
        self.assertEqual(project.info.primary_fastq_dir,'fastqs')
        self.assertEqual(project.fastq_dir,
                         os.path.join(project.dirn,'fastqs'))
        self.assertEqual(project.fastq_dirs,['fastqs','fastqs.untrimmed'])
        # Switch active fastq set
        project.use_fastq_dir('fastqs.untrimmed')
        self.assertEqual(project.info.primary_fastq_dir,'fastqs')
        self.assertEqual(project.fastq_dir,
                         os.path.join(project.dirn,'fastqs.untrimmed'))
        self.assertEqual(project.fastq_dirs,['fastqs','fastqs.untrimmed'])
        # Update the primary fastq set
        project.set_primary_fastq_dir('fastqs.untrimmed')
        self.assertEqual(project.info.primary_fastq_dir,'fastqs.untrimmed')
        self.assertEqual(project.fastq_dir,
                         os.path.join(project.dirn,'fastqs.untrimmed'))
        self.assertEqual(project.fastq_dirs,['fastqs','fastqs.untrimmed'])
        # Reload the project and check that the change has stuck
        project1 = AnalysisProject('PJB',dirn)
        self.assertEqual(project1.info.primary_fastq_dir,'fastqs.untrimmed')
        self.assertEqual(project1.fastq_dir,
                         os.path.join(project.dirn,'fastqs.untrimmed'))
        self.assertEqual(project.fastq_dirs,['fastqs','fastqs.untrimmed'])

    def test_analysis_project_switch_to_non_existant_fastq_dir(self):
        """Check AnalysisProject fails when switching to non-existant fastqs dir
        """
        self.make_mock_project_dir(
            'PJB',
            ('PJB1-A_ACAGTG_L001_R1_001.fastq.gz',
             'PJB1-B_ACAGTG_L002_R1_001.fastq.gz',))
        dirn = os.path.join(self.dirn,'PJB')
        project = AnalysisProject('PJB',dirn)
        self.assertRaises(Exception,
                          project.use_fastq_dir,'fastqs.non_existant')

class TestAnalysisSample(unittest.TestCase):
    """Tests for the AnalysisSample class

    """

    def test_empty_analysis_sample(self):
        """Check empty AnalysisSample class
        """
        sample = AnalysisSample('PJB1-A')
        self.assertEqual(sample.name,'PJB1-A')
        self.assertEqual(sample.fastq,[])
        self.assertFalse(sample.paired_end)
        self.assertEqual(str(sample),'PJB1-A')

    def test_single_end_analysis_sample(self):
        """Check AnalysisSample class with single-end sample
        """
        fq = '/run/sample1/PJB1-B_ACAGTG_L001_R1.fastq.gz'
        sample = AnalysisSample('PJB1-B')
        sample.add_fastq(fq)
        self.assertEqual(sample.name,'PJB1-B')
        self.assertEqual(sample.fastq,[fq])
        self.assertEqual(sample.fastq_subset(read_number=1),[fq])
        self.assertEqual(sample.fastq_subset(read_number=2),[])
        self.assertFalse(sample.paired_end)
        self.assertEqual(str(sample),'PJB1-B')

    def test_single_end_analysis_sample_multiple_fastqs(self):
        """Check AnalysisSample class with single-end sample (multiple fastqs)
        """
        sample = AnalysisSample('PJB1-B')
        fq_l1 = '/run/sample1/PJB1-B_ACAGTG_L001_R1.fastq.gz'
        fq_l2 = '/run/sample1/PJB1-B_ACAGTG_L002_R1.fastq.gz'
        sample.add_fastq(fq_l1)
        sample.add_fastq(fq_l2)
        self.assertEqual(sample.name,'PJB1-B')
        self.assertEqual(sample.fastq,[fq_l1,fq_l2])
        self.assertEqual(sample.fastq_subset(read_number=1),[fq_l1,fq_l2])
        self.assertEqual(sample.fastq_subset(read_number=2),[])
        self.assertFalse(sample.paired_end)
        self.assertEqual(str(sample),'PJB1-B')

    def test_paired_end_analysis_sample(self):
        """Check AnalysisSample class with paired-end sample
        """
        sample = AnalysisSample('PJB1-B')
        fq_r1 = '/run/sample1/PJB1-B_ACAGTG_L001_R1.fastq.gz'
        fq_r2 = '/run/sample1/PJB1-B_ACAGTG_L001_R2.fastq.gz'
        sample.add_fastq(fq_r1)
        sample.add_fastq(fq_r2)
        self.assertEqual(sample.name,'PJB1-B')
        self.assertEqual(sample.fastq,[fq_r1,fq_r2])
        self.assertEqual(sample.fastq_subset(read_number=1),[fq_r1])
        self.assertEqual(sample.fastq_subset(read_number=2),[fq_r2])
        self.assertTrue(sample.paired_end)
        self.assertEqual(str(sample),'PJB1-B')

    def test_paired_end_analysis_sample_multiple_fastqs(self):
        """Check AnalysisSample class with paired-end sample (multiple fastqs)
        """
        sample = AnalysisSample('PJB1-B')
        fq_l1_r1 = '/run/sample1/PJB1-B_ACAGTG_L001_R1.fastq.gz'
        fq_l2_r1 = '/run/sample1/PJB1-B_ACAGTG_L002_R1.fastq.gz'
        fq_l1_r2 = '/run/sample1/PJB1-B_ACAGTG_L001_R2.fastq.gz'
        fq_l2_r2 = '/run/sample1/PJB1-B_ACAGTG_L002_R2.fastq.gz'
        sample.add_fastq(fq_l1_r1)
        sample.add_fastq(fq_l2_r1)
        sample.add_fastq(fq_l1_r2)
        sample.add_fastq(fq_l2_r2)
        self.assertEqual(sample.name,'PJB1-B')
        self.assertEqual(sample.fastq,[fq_l1_r1,fq_l1_r2,
                                       fq_l2_r1,fq_l2_r2])
        self.assertEqual(sample.fastq_subset(read_number=1),[fq_l1_r1,fq_l2_r1])
        self.assertEqual(sample.fastq_subset(read_number=2),[fq_l1_r2,fq_l2_r2])
        self.assertTrue(sample.paired_end)
        self.assertEqual(str(sample),'PJB1-B')

    def test_analysis_sample_non_canonical_fastq_naming(self):
        """Check AnalysisSample class with non-canonical fastq naming
        """
        # Create a class to handle the non-canonical Fastq names
        class NonCanonicalFastq(BaseFastqAttrs):
            def __init__(self,fastq):
                BaseFastqAttrs.__init__(self,fastq)
                self.sample_name = self.basename.split('.')[0]
                self.read_number = int(self.basename.split('.')[-1][1:])
        sample = AnalysisSample('PJB1-B',fastq_attrs=NonCanonicalFastq)
        fq_r1 = '/run/sample1/PJB1-B.ACAGTG.L001.R1.fastq.gz'
        fq_r2 = '/run/sample1/PJB1-B.ACAGTG.L001.R2.fastq.gz'
        sample.add_fastq(fq_r1)
        sample.add_fastq(fq_r2)
        self.assertEqual(sample.name,'PJB1-B')
        self.assertEqual(sample.fastq,[fq_r1,fq_r2])
        self.assertEqual(sample.fastq_subset(read_number=1),[fq_r1])
        self.assertEqual(sample.fastq_subset(read_number=2),[fq_r2])
        self.assertTrue(sample.paired_end)
        self.assertEqual(str(sample),'PJB1-B')

class TestMetadataDict(unittest.TestCase):
    """Tests for the MetadataDict class

    """

    def setUp(self):
        self.metadata_file = None

    def tearDown(self):
        if self.metadata_file is not None:
            os.remove(self.metadata_file)

    def test_create_metadata_object(self):
        """Check creation of a metadata object
        """
        metadata = MetadataDict(attributes={'salutation':'Salutation',
                                            'valediction': 'Valediction'})
        self.assertEqual(metadata.salutation,None)
        self.assertEqual(metadata.valediction,None)

    def test_set_and_get(self):
        """Check metadata values can be stored and retrieved
        """
        metadata = MetadataDict(attributes={'salutation':'Salutation',
                                            'valediction': 'Valediction'})
        metadata['salutation'] = "hello"
        metadata['valediction'] = "goodbye"
        self.assertEqual(metadata.salutation,"hello")
        self.assertEqual(metadata.valediction,"goodbye")

    def test_save_and_load(self):
        """Check metadata can be saved to file and reloaded
        """
        self.metadata_file = tempfile.mkstemp()[1]
        metadata = MetadataDict(attributes={'salutation':'Salutation',
                                            'valediction': 'Valediction',
                                            'chat': 'Chit chat'})
        metadata['salutation'] = "hello"
        metadata['valediction'] = "goodbye"
        metadata.save(self.metadata_file)
        metadata2 = MetadataDict(attributes={'salutation':'Salutation',
                                             'valediction': 'Valediction',
                                             'chat': 'Chit chat'})
        metadata2.load(self.metadata_file)
        self.assertEqual(metadata2.salutation,"hello")
        self.assertEqual(metadata2.valediction,"goodbye")
        self.assertEqual(metadata2.chat,None)

    def test_get_non_existent_attribute(self):
        """Check that accessing non-existent attribute raises exception
        """
        metadata = MetadataDict(attributes={'salutation':'Salutation',
                                            'valediction': 'Valediction'})
        self.assertRaises(AttributeError,lambda: metadata.conversation)

    def test_set_non_existent_attribute(self):
        """Check that setting non-existent attribute raises exception
        """
        metadata = MetadataDict(attributes={'salutation':'Salutation',
                                            'valediction': 'Valediction'})
        try:
            metadata['conversation'] = 'hrm'
            self.fail('AttributeError not raised')
        except AttributeError,ex:
            pass

    def test_specify_key_order(self):
        """Check that specified key ordering is respected
        """
        self.metadata_file = tempfile.mkstemp()[1]
        expected_keys = ('Salutation',
                         'Chit chat',
                         'Valediction',)
        metadata = MetadataDict(attributes={'salutation':'Salutation',
                                            'valediction': 'Valediction',
                                            'chat': 'Chit chat'},
                                order=('salutation','chat','valediction'))
        metadata.save(self.metadata_file)
        fp = open(self.metadata_file,'rU')
        for line,expected_key in zip(fp,expected_keys):
            self.assertEqual(line.split('\t')[0],expected_key)

    def test_implicit_key_order(self):
        """Check that keys are implicitly ordered on output
        """
        self.metadata_file = tempfile.mkstemp()[1]
        metadata = MetadataDict(attributes={'salutation':'Salutation',
                                            'valediction': 'Valediction',
                                            'chat': 'Chit chat'})
        expected_keys = ('Chit chat',
                         'Salutation',
                         'Valediction',)
        metadata.save(self.metadata_file)
        fp = open(self.metadata_file,'rU')
        for line,expected_key in zip(fp,expected_keys):
            self.assertEqual(line.split('\t')[0],expected_key)

    def test_get_null_items(self):
        """Check fetching of items with null values
        """
        metadata = MetadataDict(attributes={'salutation':'.',
                                            'valediction': '.',
                                            'chat': '.'})
        self.assertEqual(metadata.null_items(),['chat',
                                                'salutation',
                                                'valediction'])
        #
        metadata['salutation'] = 'hello'
        metadata['valediction'] = 'ave'
        metadata['chat'] = 'awight'
        self.assertEqual(metadata.null_items(),[])
        #
        metadata['chat'] = None
        self.assertEqual(metadata.null_items(),['chat'])

    def test_undefined_items_in_file(self):
        """Check handling of additional undefined items in file
        """
        # Set up a metadata dictionary
        metadata = MetadataDict(attributes={'salutation':'salutation',
                                            'valediction': 'valediction'})
        # Create a file with an additional item
        self.metadata_file = tempfile.mkstemp()[1]
        contents = ('salutation\thello',
                    'valediction\tgoodbye',
                    'chit_chat\tstuff')
        with open(self.metadata_file,'w') as fp:
            for line in contents:
                fp.write("%s\n" % line)
        # Load into the dictionary and check that all
        # items are present
        metadata.load(self.metadata_file,strict=False)
        self.assertEqual(metadata.salutation,'hello')
        self.assertEqual(metadata.valediction,'goodbye')
        self.assertEqual(metadata.chit_chat,'stuff')

class TestAnalysisDirParameters(unittest.TestCase):
    """Tests for the AnalysisDirParameters class

    """

    def test_create_analysis_dir_parameters(self):
        """Check creation of an empty AnalysisDirParameters object
        """
        params = AnalysisDirParameters()
        self.assertEqual(params.analysis_dir,None)
        self.assertEqual(params.data_dir,None)
        self.assertEqual(params.sample_sheet,None)
        self.assertEqual(params.bases_mask,None)
        self.assertEqual(params.primary_data_dir,None)
        self.assertEqual(params.unaligned_dir,None)
        self.assertEqual(params.project_metadata,None)
        self.assertEqual(params.stats_file,None)

class TestAnalysisDirParameters(unittest.TestCase):
    """Tests for the AnalysisDirParameters class

    """

    def test_create_analysis_dir_parameters(self):
        """Check creation of an empty AnalysisDirParameters object
        """
        params = AnalysisDirParameters()
        self.assertEqual(params.analysis_dir,None)
        self.assertEqual(params.data_dir,None)
        self.assertEqual(params.sample_sheet,None)
        self.assertEqual(params.bases_mask,None)
        self.assertEqual(params.primary_data_dir,None)
        self.assertEqual(params.unaligned_dir,None)
        self.assertEqual(params.project_metadata,None)
        self.assertEqual(params.stats_file,None)

class TestAnalysisDirMetadata(unittest.TestCase):
    """Tests for the AnalysisDirMetadata class

    """

    def test_create_analysis_dir_metadata(self):
        """Check creation of an empty AnalysisDirMetadata object
        """
        metadata = AnalysisDirMetadata()
        self.assertEqual(metadata.run_number,None)
        self.assertEqual(metadata.platform,None)
        self.assertEqual(metadata.source,None)
        self.assertEqual(metadata.assay,None)
        self.assertEqual(metadata.bcl2fastq_software,None)

class TestProjectMetadataFile(unittest.TestCase):
    """Tests for the ProjectMetadataFile class

    """

    def setUp(self):
        self.metadata_file = tempfile.mkstemp()[1]
        self.projects = list()
        self.lines = list()

    def tearDown(self):
        if self.metadata_file is not None:
            os.remove(self.metadata_file)

    def test_empty_project_metadata_file(self):
        """Create and save empty ProjectMetadataFile
        """
        # Make an empty 'file'
        metadata = ProjectMetadataFile()
        contents = "#Project\tSamples\tUser\tLibrary\tOrganism\tPI\tComments\n"
        self.assertEqual(len(metadata),0)
        for project in metadata:
            self.fail()
        # Save to an actual file and check its contents
        metadata.save(self.metadata_file)
        self.assertEqual(open(self.metadata_file,'r').read(),contents)

    def test_create_new_project_metadata_file(self):
        """Create and save ProjectMetadataFile with content
        """
        # Make new 'file' and add projects
        metadata = ProjectMetadataFile()
        metadata.add_project('Charlie',['C1','C2'],
                             user="Charlie P",
                             library_type="RNA-seq",
                             organism="Yeast",
                             PI="Marley")
        metadata.add_project('Farley',['F3','F4'],
                             user="Farley G",
                             library_type="ChIP-seq",
                             organism="Mouse",
                             PI="Harley",
                             comments="Squeak!")
        contents = "#Project\tSamples\tUser\tLibrary\tOrganism\tPI\tComments\nCharlie\tC1,C2\tCharlie P\tRNA-seq\tYeast\tMarley\t.\nFarley\tF3,F4\tFarley G\tChIP-seq\tMouse\tHarley\tSqueak!\n"
        self.assertEqual(len(metadata),2)
        # Save to an actual file and check its contents
        metadata.save(self.metadata_file)
        self.assertEqual(open(self.metadata_file,'r').read(),contents)

    def test_read_existing_project_metadata_file(self):
        """Read contents from existing ProjectMetadataFile
        """
        # Create metadata file independently
        data = list()
        data.append(dict(Project="Charlie",
                         Samples="C1-2",
                         User="Charlie P",
                         Library="RNA-seq",
                         Organism="Yeast",
                         PI="Marley",
                         Comments="."))
        data.append(dict(Project="Farley",
                         Samples="F3-4",
                         User="Farley G",
                         Library="ChIP-seq",
                         Organism="Mouse",
                         PI="Harley",
                         Comments="Squeak!"))
        contents = "#Project\tSamples\tUser\tLibrary\tOrganism\tPI\tComments\nCharlie\tC1-2\tCharlie P\tRNA-seq\tYeast\tMarley\t.\nFarley\tF3-4\tFarley G\tChIP-seq\tMouse\tHarley\tSqueak!\n"
        open(self.metadata_file,'w').write(contents)
        # Load and check contents
        metadata = ProjectMetadataFile(self.metadata_file)
        self.assertEqual(len(metadata),2)
        for x,y in zip(data,metadata):
            for attr in ('Project','User','Library','Organism','PI','Comments'):
                self.assertEqual(x[attr],y[attr])

    def test_projects_with_numbers_for_names(self):
        """Handle projects with names that look like numbers
        """
        metadata = ProjectMetadataFile()
        metadata.add_project('1234',['C1','C2'],
                             user="Charlie P",
                             library_type="RNA-seq",
                             organism="Yeast",
                             PI="Marley")
        metadata.add_project('56.78',['F3','F4'],
                             user="Farley G",
                             library_type="ChIP-seq",
                             organism="Mouse",
                             PI="Harley",
                             comments="Squeak!")
        for project in metadata:
            print "%s" % project
            print "%s" % type(project['Project'])
            self.assertTrue(isinstance(project['Project'],str))

    def test_refuse_to_add_duplicate_projects(self):
        """Refuse to add duplicated project names
        """
        metadata = ProjectMetadataFile()
        # Make new 'file' and add project
        metadata = ProjectMetadataFile()
        metadata.add_project('Charlie',['C1','C2'],
                             user="Charlie P",
                             library_type="RNA-seq",
                             organism="Yeast",
                             PI="Marley")
        # Attempt to add same project name again
        self.assertRaises(Exception,
                          metadata.add_project,'Charlie',['C1','C2'])

class TestZipArchive(unittest.TestCase):
    """
    Tests for the ZipArchive class
    """
    def setUp(self):
        self.dirn = tempfile.mkdtemp(suffix='TestZipArchive')

    def tearDown(self):
        # Remove the temporary test directory
        shutil.rmtree(self.dirn)

    def test_create_zip_archive(self):
        """ZipArchive: create a new zip archive
        """
        # Make an example directory to zip up
        src_dir = os.path.join(self.dirn,'source')
        items = ('test1',
                 'test2',
                 'sub1/',
                 'sub1/test3',
                 'sub1/test4',
                 'sub1/sub2/',
                 'sub1/sub2/test5',
                 'sub3/',
                 'sub3/test6')
        os.mkdir(src_dir)
        for item in [os.path.join(self.dirn,x) for x in items]:
            if item.endswith('/'):
                try:
                    os.mkdir(item)
                except OSError as ex:
                    raise OSError("Failed to make %s" % item)
            else:
                open(item,'w').write('')
        # Create the zip archive
        zip_filename = os.path.join(self.dirn,'test.zip')
        self.assertFalse(os.path.exists(zip_filename))
        z = ZipArchive(zip_filename,relpath=self.dirn)
        for item in [os.path.join(self.dirn,x) for x in ('test1','sub1')]:
            z.add(item)
        z.close()
        # Check the zip archive exists
        self.assertTrue(os.path.exists(zip_filename),
                        "zip file not created")
        self.assertTrue(zipfile.is_zipfile(zip_filename),
                        "file exists but is not zip file")
        # Check the contents
        namelist = zipfile.ZipFile(zip_filename).namelist()
        print namelist
        for item in ('test1',
                     'sub1/test3',
                     'sub1/test4',
                     'sub1/sub2/test5'):
            self.assertTrue((item.rstrip('/') in namelist),
                            "Item '%s' not in zip file" % item)
        for name in namelist:
            self.assertTrue((name in ('test1',
                                      'sub1/test3',
                                      'sub1/test4',
                                      'sub1/sub2/test5')),
                            "'%s' in zip file but shouldn't be"
                            % name)

class TestOutputFiles(unittest.TestCase):
    """
    Tests for the OutputFiles class
    """
    def setUp(self):
        # Temporary working dir
        self.wd = tempfile.mkdtemp(suffix='.test_outputfile')
    def tearDown(self):
        # Remove temporary working dir
        if os.path.isdir(self.wd):
            shutil.rmtree(self.wd)
    def test_outputfiles(self):
        out = OutputFiles()
        out.open('test1',os.path.join(self.wd,'test1.txt'))
        out.open('test2',os.path.join(self.wd,'test2.txt'))
        self.assertEqual(
            out.file_name('test1'),os.path.join(self.wd,'test1.txt'))
        self.assertEqual(
            out.file_name('test2'),os.path.join(self.wd,'test2.txt'))
        out.write('test1','Some test text')
        out.write('test2','Some more\ntest text')
        out.close()
        self.assertEqual(open(out.file_name('test1'),'r').read(),
                         "Some test text\n")
        self.assertEqual(open(out.file_name('test2'),'r').read(),
                         "Some more\ntest text\n")
    def test_outputfiles_with_basedir(self):
        out = OutputFiles(self.wd)
        out.open('test1','test1.txt')
        out.open('test2','test2.txt')
        self.assertEqual(
            out.file_name('test1'),os.path.join(self.wd,'test1.txt'))
        self.assertEqual(
            out.file_name('test2'),os.path.join(self.wd,'test2.txt'))
    def test_outputfiles_contains(self):
        out = OutputFiles()
        out.open('test1',os.path.join(self.wd,'test1.txt'))
        self.assertTrue('test1' in out)
        self.assertFalse('test2' in out)
        out.close()
        self.assertFalse('test1' in out)
        self.assertEqual(out.file_name('test1'),
                         os.path.join(self.wd,'test1.txt'))
        self.assertFalse('test2' in out)
        self.assertRaises(KeyError,out.file_name,'test2')
    def test_outputfiles_append(self):
        out = OutputFiles()
        with open(os.path.join(self.wd,'test1.txt'),'w') as fp:
            fp.write("test1 already exists\n")
        with open(os.path.join(self.wd,'test2.txt'),'w') as fp:
            fp.write("test2 already exists\n")
        out.open('test1',os.path.join(self.wd,'test1.txt'),append=True)
        out.open('test2',os.path.join(self.wd,'test2.txt'))
        self.assertEqual(
            out.file_name('test1'),os.path.join(self.wd,'test1.txt'))
        self.assertEqual(
            out.file_name('test2'),os.path.join(self.wd,'test2.txt'))
        out.write('test1','Some test text')
        out.write('test2','Some more\ntest text')
        out.close()
        self.assertEqual(open(out.file_name('test1'),'r').read(),
                         "test1 already exists\nSome test text\n")
        self.assertEqual(open(out.file_name('test2'),'r').read(),
                         "Some more\ntest text\n")
    def test_outputfiles_reopen(self):
        out = OutputFiles(base_dir=self.wd)
        out.open('test1','test1.txt')
        out.open('test2','test2.txt')
        self.assertEqual(out.file_name('test1'),
                         os.path.join(self.wd,'test1.txt'))
        self.assertEqual(out.file_name('test2'),
                         os.path.join(self.wd,'test2.txt'))
        out.write('test1','Some test text')
        out.write('test2','Some more\ntest text')
        out.close()
        self.assertEqual(open(out.file_name('test1'),'r').read(),
                         "Some test text\n")
        self.assertEqual(open(out.file_name('test2'),'r').read(),
                         "Some more\ntest text\n")
        # Reopen files
        out.open('test1',append=True)
        out.open('test2')
        self.assertEqual(out.file_name('test1'),
                         os.path.join(self.wd,'test1.txt'))
        self.assertEqual(out.file_name('test2'),
                         os.path.join(self.wd,'test2.txt'))
        out.write('test1','Some extra test text')
        out.write('test2','Some different\ntest text')
        out.close()
        self.assertEqual(open(out.file_name('test1'),'r').read(),
                         "Some test text\nSome extra test text\n")
        self.assertEqual(open(out.file_name('test2'),'r').read(),
                         "Some different\ntest text\n")
    def test_outputfiles_len(self):
        out = OutputFiles(base_dir=self.wd)
        self.assertEqual(len(out),0)
        out.open('test1','test1.txt')
        self.assertEqual(len(out),1)
        out.close('test1')
        self.assertEqual(len(out),0)
        out.open('test2','test2.txt')
        out.open('test3','test3.txt')
        self.assertEqual(len(out),2)
        out.open('test1',append=True)
        self.assertEqual(len(out),3)
        out.close()
        self.assertEqual(len(out),0)
        out.open('test2',append=True)
        self.assertEqual(len(out),1)

class TestBasesMaskIsPairedEnd(unittest.TestCase):
    """Tests for the bases_mask_is_paired_end function

    """

    def test_single_end_bases_mask(self):
        """Check examples of single-ended bases masks
        """
        self.assertFalse(bases_mask_is_paired_end('y36'))
        self.assertFalse(bases_mask_is_paired_end('y50,I6'))
        self.assertFalse(bases_mask_is_paired_end('y50,I6n'))

    def test_paired_end_bases_mask(self):
        """Check examples of single-ended bases masks
        """
        self.assertTrue(bases_mask_is_paired_end('y101,I6n,y101'))
        self.assertTrue(bases_mask_is_paired_end('y101,I6,I6,y101'))

class TestSplitUserHostDir(unittest.TestCase):
    """Tests for the split_user_host_dir function

    """

    def test_user_host_dir(self):
        """Check we can split user@host.name:/path/to/somewhere
        """
        user,host,dirn = split_user_host_dir('user@host.name:/path/to/somewhere')
        self.assertEqual(user,'user')
        self.assertEqual(host,'host.name')
        self.assertEqual(dirn,'/path/to/somewhere')

    def test_host_dir(self):
        """Check we can split host.name:/path/to/somewhere
        """
        user,host,dirn = split_user_host_dir('host.name:/path/to/somewhere')
        self.assertEqual(user,None)
        self.assertEqual(host,'host.name')
        self.assertEqual(dirn,'/path/to/somewhere')

    def test_dir(self):
        """Check we can 'split' /path/to/somewhere
        """
        user,host,dirn = split_user_host_dir('/path/to/somewhere')
        self.assertEqual(user,None)
        self.assertEqual(host,None)
        self.assertEqual(dirn,'/path/to/somewhere')

    def test_extra_whitespace(self):
        """Check we can handle addition leading/trailing whitespace
        """
        user,host,dirn = split_user_host_dir('\tuser@host.name:/path/to/somewhere  \n')
        self.assertEqual(user,'user')
        self.assertEqual(host,'host.name')
        self.assertEqual(dirn,'/path/to/somewhere')

    def test_bad_input(self):
        """Check that empty string or None return sensible values
        """
        user,host,dirn = split_user_host_dir('')
        self.assertEqual(user,None)
        self.assertEqual(host,None)
        self.assertEqual(dirn,None)
        user,host,dirn = split_user_host_dir(None)
        self.assertEqual(user,None)
        self.assertEqual(host,None)
        self.assertEqual(dirn,None)

class TestPrettyPrintRows(unittest.TestCase):
    """Tests for the pretty_print_rows function

    """
    def test_pretty_print_rows_empty_inputs(self):
        """pretty_print_rows can handle empty list
        """
        self.assertEqual(pretty_print_rows([]),"")

    def test_pretty_print_rows_single_row(self):
        """pretty_print_rows prints a single row
        """
        rows = [['hello:','A salutation']]
        self.assertEqual(pretty_print_rows(rows),
                         "hello: A salutation")

    def test_pretty_print_rows_multiple_rows(self):
        """pretty_print_rows prints multiple rows
        """
        rows = [['-','hello:','A salutation'],
                ['-','goodbye:','The End']]
        self.assertEqual(pretty_print_rows(rows),
                         "- hello:   A salutation\n- goodbye: The End     ")

    def test_pretty_print_rows_prepend(self):
        """pretty_print_rows prepend right-justifies values
        """
        rows = [['-','hello:','A salutation'],
                ['-','goodbye:','The End']]
        self.assertEqual(pretty_print_rows(rows,prepend=True),
                         "-   hello: A salutation\n- goodbye:      The End")

class TestWriteScriptFile(unittest.TestCase):
    """Tests for the write_script_file function

    """
    def setUp(self):
        self.dirn = tempfile.mkdtemp(suffix='TestWriteScripFile')

    def tearDown(self):
        # Remove the temporary test directory
        shutil.rmtree(self.dirn)

    def test_write_script_file(self):
        """write_script_file creates new script file
        """
        script_file = os.path.join(self.dirn,'test.sh')
        self.assertFalse(os.path.exists(script_file))
        write_script_file(script_file,"sleep 50\n#")
        self.assertTrue(os.path.exists(script_file))
        self.assertEqual(open(script_file).read(),
                         "sleep 50\n#\n")

    def test_write_script_file_with_shell(self):
        """write_script_file creates new script file with shell
        """
        script_file = os.path.join(self.dirn,'test.sh')
        self.assertFalse(os.path.exists(script_file))
        write_script_file(script_file,"sleep 50\n#",shell='/bin/sh')
        self.assertTrue(os.path.exists(script_file))
        self.assertEqual(open(script_file).read(),
                         "#!/bin/sh\nsleep 50\n#\n")

    def test_write_script_file_append(self):
        """write_script_file appends to an existing file
        """
        script_file = os.path.join(self.dirn,'test.sh')
        with open(script_file,'w') as fp:
            fp.write("#\necho Going to sleep\n")
        write_script_file(script_file,"sleep 50\n#",append=True)
        self.assertTrue(os.path.exists(script_file))
        self.assertEqual(open(script_file).read(),
                         "#\necho Going to sleep\nsleep 50\n#\n")
