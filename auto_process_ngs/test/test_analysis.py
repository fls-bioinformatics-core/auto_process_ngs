#######################################################################
# Tests for analysis.py module
#######################################################################

import unittest
import tempfile
import shutil
import zipfile
import pickle
import cloudpickle
from bcftbx.JobRunner import SimpleJobRunner,GEJobRunner
from bcftbx.utils import find_program
from auto_process_ngs.mock import MockAnalysisDirFactory
from auto_process_ngs.mock import MockAnalysisProject
from auto_process_ngs.applications import Command
from auto_process_ngs.fastq_utils import BaseFastqAttrs
from auto_process_ngs.analysis import *

class TestAnalysisFastq(unittest.TestCase):
    """
    Tests for the AnalysisFastq class
    """
    def test_illumina_style_with_extras(self):
        """AnalysisFastq: Illumina-style fastq name with extra elements appended
        """
        fq = AnalysisFastq('M_19_0040_S13-A40h-R_D709-D505_L007_R1_001_repeated_1')
        self.assertEqual(fq.sample_name,'M_19_0040_S13-A40h-R_D709-D505')
        self.assertEqual(fq.basename,
                         'M_19_0040_S13-A40h-R_D709-D505_L007_R1_001_repeated_1')
        self.assertEqual(fq.extension,'')
        self.assertEqual(fq.sample_number,None)
        self.assertEqual(fq.barcode_sequence,None)
        self.assertEqual(fq.lane_number,7)
        self.assertEqual(fq.read_number,1)
        self.assertEqual(fq.set_number,1)
        self.assertFalse(fq.is_index_read)
        self.assertEqual(fq.canonical_name,'M_19_0040_S13-A40h-R_D709-D505_L007_R1_001')
        self.assertEqual(fq.extras,'_repeated_1')
        self.assertEqual(str(fq),
                         'M_19_0040_S13-A40h-R_D709-D505_L007_R1_001_repeated_1')
    def test_non_canonical_fastq_name(self):
        """AnalysisFastq: Illumina-style fastq name with extra elements appended
        """
        fq = AnalysisFastq('PB04_trimmoPE_bowtie2_notHg38.1')
        self.assertEqual(fq.sample_name,'PB04_trimmoPE_bowtie2_notHg38.1')
        self.assertEqual(fq.basename,'PB04_trimmoPE_bowtie2_notHg38.1')
        self.assertEqual(fq.extension,'')
        self.assertEqual(fq.sample_number,None)
        self.assertEqual(fq.barcode_sequence,None)
        self.assertEqual(fq.lane_number,None)
        self.assertEqual(fq.read_number,None)
        self.assertEqual(fq.set_number,None)
        self.assertFalse(fq.is_index_read)
        self.assertEqual(fq.canonical_name,None)
        self.assertEqual(fq.extras,None)
        self.assertEqual(str(fq),'PB04_trimmoPE_bowtie2_notHg38.1')
    def test_reproduces_illumina_fastq_attrs(self):
        """AnalysisFastq: reproduces IlluminaFastqAttr behaviour
        """
        for name in ('NH1_ChIP-seq_Gli1_ACAGTG_L003_R2_001',
                     'NH1_ChIP-seq_Gli1_ACAGTG-GTTCAC_L003_R2_001',
                     'NH1_ChIP-seq_Gli1_S4_L003_R2_001',
                     'NH1_ChIP-seq_Gli1_S4_L003_I1_001',
                     'NH1_ChIP-seq_Gli1_S4_R2_001',
                     'NH1_ChIP-seq_Gli1',
                     'NH1_ChIP-seq_Gli1_R2',
                     'NH1_ChIP-seq_Gli1_L001',
                     'NH1_ChIP-seq_Gli1_L001_R2',
                     'NH1_ChIP-seq_Gli1_ACAGTG',
                     'NH1_ChIP-seq_Gli1_ACAGTG_R2',
                     'NH1_ChIP-seq_Gli1_ACAGTG_L001',
                     'NH1_ChIP-seq_Gli1_ACAGTG_L001_R2',):
            illumina_fastq_attrs = IlluminaFastqAttrs(name)
            fq = AnalysisFastq(name)
            self.assertEqual(fq.sample_name,illumina_fastq_attrs.sample_name)
            self.assertEqual(fq.basename,illumina_fastq_attrs.basename)
            self.assertEqual(fq.extension,illumina_fastq_attrs.extension)
            self.assertEqual(fq.sample_number,illumina_fastq_attrs.sample_number)
            self.assertEqual(fq.barcode_sequence,illumina_fastq_attrs.barcode_sequence)
            self.assertEqual(fq.lane_number,illumina_fastq_attrs.lane_number)
            self.assertEqual(fq.read_number,illumina_fastq_attrs.read_number)
            self.assertEqual(fq.is_index_read,illumina_fastq_attrs.is_index_read)
            self.assertEqual(str(fq),str(illumina_fastq_attrs))
    def test_changing_attribute_updates_canonical_name(self):
        """AnalysisFastq: changing attributes is reflected in canonical name
        """
        fq = AnalysisFastq('M_19_0040_S13-A40h-R_D709-D505_L007_R1_001_repeated_1')
        fq.read_number = 2
        self.assertEqual(fq.canonical_name,
                         'M_19_0040_S13-A40h-R_D709-D505_L007_R2_001')
    def test_changing_attribute_updates_repr(self):
        """AnalysisFastq: changing attributes is reflected in repr
        """
        fq = AnalysisFastq('M_19_0040_S13-A40h-R_D709-D505_L007_R1_001_repeated_1')
        fq.read_number = 2
        self.assertEqual(str(fq),
                         'M_19_0040_S13-A40h-R_D709-D505_L007_R2_001_repeated_1')

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
        sample_names = list(
            set([AnalysisFastq(fq).sample_name for fq in fastq_list]))
        sample_names.sort()
        n_samples = len(sample_names)
        sample_names = "%d sample%s (%s)" % \
                       (n_samples,
                        '' if n_samples == 1 else 's',
                        ', '.join(sample_names))
        metadata = { 'Primary fastqs': primary_fastq_dir,
                     'Samples': sample_names }
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
        self.assertEqual(project.info.single_cell_platform,None)
        self.assertEqual(project.info.organism,None)
        self.assertEqual(project.info.number_of_cells,None)
        self.assertEqual(project.info.icell8_well_list,None)
        self.assertFalse(project.info.paired_end)
        self.assertEqual(project.info.platform,None)
        self.assertEqual(project.info.primary_fastq_dir,None)
        self.assertEqual(project.info.samples,None)
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
        self.assertEqual(project.info.samples,'2 samples (PJB1-A, PJB1-B)')
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
        self.assertEqual(project.info.samples,
                         '2 samples (PJB1-A, PJB1-B, multiple fastqs per sample)')
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
        project = AnalysisProject('PJB',dirn)
        project.create_directory(fastqs=self.fastqs)
        self.assertEqual(project.name,'PJB')
        self.assertTrue(os.path.isdir(project.dirn))
        self.assertFalse(project.multiple_fastqs)
        self.assertTrue(project.info.paired_end)
        self.assertEqual(project.info.primary_fastq_dir,'fastqs')
        self.assertEqual(project.info.samples,'2 samples (PJB1-A, PJB1-B)')
        self.assertEqual(project.samples[0].name,'PJB1-A')
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
        self.assertEqual(project.info.samples,
                         '1 sample (PJB1-B, multiple fastqs per sample)')

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
        self.assertEqual(project.info.samples,
                         '1 sample (PJB1-B, multiple fastqs per sample)')
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
        self.assertEqual(project.info.samples,'2 samples (PJB1-A, PJB1-B)')
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
        self.assertEqual(project.info.samples,'2 samples (PJB1-A, PJB1-B)')
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
        self.assertEqual(project.info.samples,
                         '2 samples (PB04, PB04_trimmoPE_bowtie2_notHg38.1)')
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
        self.assertEqual(project.info.samples,
                         '1 sample (PB02.trimmed.filtered)')
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
        self.assertEqual(project.info.samples,'2 samples (PJB1-A, PJB1-B)')
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
        self.assertEqual(project.info.samples,'2 samples (PJB1-A, PJB1-B)')
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
        self.assertEqual(project.info.samples,'2 samples (PJB1-A, PJB1-B)')
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
        self.assertEqual(project.info.samples,'2 samples (PJB1-A, PJB1-B)')
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
        self.assertEqual(project.info.samples,'2 samples (PJB1-A, PJB1-B)')
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
        self.assertEqual(project.info.samples,
                         '2 samples (PJB1-A-untrimmed, PJB1-B-untrimmed)')
        # Load again with alternative fastq set
        project = AnalysisProject('PJB',dirn,fastq_dir='fastqs')
        self.assertEqual(project.info.primary_fastq_dir,'fastqs.untrimmed')
        self.assertEqual(project.fastq_dir,
                         os.path.join(project.dirn,'fastqs'))
        self.assertEqual(project.info.samples,
                         '2 samples (PJB1-A-untrimmed, PJB1-B-untrimmed)')
        # Implicitly switch to primary fastq set
        project.use_fastq_dir()
        self.assertEqual(project.info.primary_fastq_dir,'fastqs.untrimmed')
        self.assertEqual(project.fastq_dir,
                         os.path.join(project.dirn,'fastqs.untrimmed'))
        self.assertEqual(project.info.samples,
                         '2 samples (PJB1-A-untrimmed, PJB1-B-untrimmed)')

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
        self.assertEqual(project.info.samples,'2 samples (PJB1-A, PJB1-B)')
        self.assertEqual(project.fastq_dir,
                         os.path.join(project.dirn,'fastqs'))
        self.assertEqual(project.fastq_dirs,['fastqs','fastqs.untrimmed'])
        # Update the primary fastq set
        project.set_primary_fastq_dir('fastqs.untrimmed')
        self.assertEqual(project.info.primary_fastq_dir,'fastqs.untrimmed')
        self.assertEqual(project.fastq_dir,
                         os.path.join(project.dirn,'fastqs'))
        self.assertEqual(project.info.samples,
                         '2 samples (PJB1-A-untrimmed, PJB1-B-untrimmed)')
        self.assertEqual(project.fastq_dirs,['fastqs','fastqs.untrimmed'])
        # Reload the project and check that the change has stuck
        project1 = AnalysisProject('PJB',dirn)
        self.assertEqual(project1.info.primary_fastq_dir,'fastqs.untrimmed')
        self.assertEqual(project1.fastq_dir,
                         os.path.join(project.dirn,'fastqs.untrimmed'))
        self.assertEqual(project1.fastq_dirs,['fastqs','fastqs.untrimmed'])
        self.assertEqual(project1.info.samples,
                         '2 samples (PJB1-A-untrimmed, PJB1-B-untrimmed)')

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

    def test_analysis_project_switch_fastq_dir_preserves_qc_dir(self):
        """Check AnalysisProject switch fastqs dirs preserves QC dir
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
        self.assertEqual(project.qc_dir,
                         os.path.join(project.dirn,'qc'))
        # Set new QC dir
        project.use_qc_dir('qc.new')
        self.assertEqual(project.qc_dir,
                         os.path.join(project.dirn,'qc.new'))
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
        self.assertEqual(project.info.samples,'2 samples (PJB1-A, PJB1-B)')
        self.assertEqual(project.fastq_dirs,
                         ['fastqs','fastqs.untrimmed'])
        self.assertEqual(project.info.primary_fastq_dir,'fastqs')
        self.assertEqual(project.qc_dir,
                         os.path.join(project.dirn,'qc.new'))

    def test_sample_summary_single_ended(self):
        """AnalysisProject: sample_summary works for SE data
        """
        self.make_mock_project_dir(
            'PJB',
            ('PJB1-A_ACAGTG_L001_R1_001.fastq.gz',
             'PJB1-B_ACAGTG_L002_R1_001.fastq.gz',))
        project = AnalysisProject('PJB',os.path.join(self.dirn,'PJB'))
        self.assertEqual(project.sample_summary(),
                         "2 samples (PJB1-A, PJB1-B)")

    def test_sample_summary_paired_ended(self):
        """AnalysisProject: sample_summary works for PE data
        """
        self.make_mock_project_dir(
            'PJB',
            ('PJB1-A_ACAGTG_L001_R1_001.fastq.gz',
             'PJB1-A_ACAGTG_L001_R2_001.fastq.gz',
             'PJB1-B_ACAGTG_L002_R1_001.fastq.gz',
             'PJB1-B_ACAGTG_L002_R2_001.fastq.gz',))
        project = AnalysisProject('PJB',os.path.join(self.dirn,'PJB'))
        self.assertEqual(project.sample_summary(),
                         "2 samples (PJB1-A, PJB1-B)")

    def test_sample_summary_single_ended_multiple_fastqs(self):
        """AnalysisProject: sample_summary works for SE data, multiple fastqs
        """
        self.make_mock_project_dir(
            'PJB',
            ('PJB1-A_ACAGTG_L001_R1_001.fastq.gz',
             'PJB1-A_ACAGTG_L002_R1_001.fastq.gz',
             'PJB1-B_ACAGTG_L001_R1_001.fastq.gz',
             'PJB1-B_ACAGTG_L002_R1_001.fastq.gz',))
        project = AnalysisProject('PJB',os.path.join(self.dirn,'PJB'))
        self.assertEqual(project.sample_summary(),
                         "2 samples (PJB1-A, PJB1-B, multiple fastqs per sample)")

    def test_sample_summary_paired_ended_multiple_fastqs(self):
        """AnalysisProject: sample_summary works for PE data, multiple fastqs
        """
        self.make_mock_project_dir(
            'PJB',
            ('PJB1-A_ACAGTG_L001_R1_001.fastq.gz',
             'PJB1-A_ACAGTG_L001_R2_001.fastq.gz',
             'PJB1-A_ACAGTG_L002_R1_001.fastq.gz',
             'PJB1-A_ACAGTG_L002_R2_001.fastq.gz',
             'PJB1-B_ACAGTG_L001_R1_001.fastq.gz',
             'PJB1-B_ACAGTG_L002_R2_001.fastq.gz',
             'PJB1-B_ACAGTG_L001_R1_001.fastq.gz',
             'PJB1-B_ACAGTG_L002_R2_001.fastq.gz',))
        project = AnalysisProject('PJB',os.path.join(self.dirn,'PJB'))
        self.assertEqual(project.sample_summary(),
                         "2 samples (PJB1-A, PJB1-B, multiple fastqs per sample)")

    def test_sample_summary_paired_ended_ignore_index_reads(self):
        """AnalysisProject: sample_summary works for PE data with index reads
        """
        self.make_mock_project_dir(
            'PJB',
            ('PJB1-A_ACAGTG_L001_R1_001.fastq.gz',
             'PJB1-A_ACAGTG_L001_R2_001.fastq.gz',
             'PJB1-A_ACAGTG_L001_I1_001.fastq.gz',
             'PJB1-B_ACAGTG_L002_R1_001.fastq.gz',
             'PJB1-B_ACAGTG_L002_R2_001.fastq.gz',
             'PJB1-B_ACAGTG_L002_I1_001.fastq.gz',))
        project = AnalysisProject('PJB',os.path.join(self.dirn,'PJB'))
        self.assertEqual(project.sample_summary(),
                         "2 samples (PJB1-A, PJB1-B)")

    def test_sample_summary_no_samples(self):
        """AnalysisProject: sample_summary works when there are no samples
        """
        self.make_mock_project_dir('PJB',())
        project = AnalysisProject('PJB',os.path.join(self.dirn,'PJB'))
        self.assertEqual(project.sample_summary(),"No samples")

    def test_order_samples_by_name(self):
        """AnalysisProject: sample_summary works for SE data
        """
        self.make_mock_project_dir(
            'PJB',
            ('PJB1_ACAGTG_L001_R1_001.fastq.gz',
             'PJB2_ACAGTG_L001_R1_001.fastq.gz',
             'PJB3_ACAGTG_L001_R1_001.fastq.gz',
             'PJB10_ACAGTG_L001_R1_001.fastq.gz',
             'PJB20_ACAGTG_L001_R1_001.fastq.gz',
             'PJB21_ACAGTG_L001_R1_001.fastq.gz',))
        project = AnalysisProject('PJB',os.path.join(self.dirn,'PJB'))
        sample_names = [s.name for s in project.samples]
        self.assertEqual(sample_names,
                         ['PJB1','PJB2','PJB3','PJB10','PJB20','PJB21'])

    def test_pickle_analysis_project(self):
        """AnalysisProject: check serialisation with 'pickle'
        """
        self.make_mock_project_dir(
            'PJB',
            ('PJB1-A_ACAGTG_L001_R1_001.fastq.gz',
             'PJB1-B_ACAGTG_L002_R1_001.fastq.gz',))
        dirn = os.path.join(self.dirn,'PJB')
        project = AnalysisProject('PJB',dirn)
        # Pickle project
        pickled = pickle.dumps(project)
        # Unpickle it
        unpickled = pickle.loads(pickled)
        # Check the unpickled data
        self.assertEqual(unpickled.name,'PJB')
        self.assertTrue(os.path.isdir(unpickled.dirn))
        self.assertFalse(unpickled.multiple_fastqs)
        self.assertFalse(unpickled.info.paired_end)
        self.assertEqual(unpickled.info.primary_fastq_dir,'fastqs')
        self.assertEqual(unpickled.info.samples,'2 samples (PJB1-A, PJB1-B)')
        self.assertEqual(unpickled.samples[0].name,'PJB1-A')
        self.assertEqual(unpickled.samples[1].name,'PJB1-B')
        self.assertEqual(unpickled.fastq_dir,
                         os.path.join(unpickled.dirn,'fastqs'))
        self.assertEqual(unpickled.fastq_dirs,['fastqs',])

    def test_cloudpickle_analysis_project(self):
        """AnalysisProject: check serialisation with 'cloudpickle'
        """
        self.make_mock_project_dir(
            'PJB',
            ('PJB1-A_ACAGTG_L001_R1_001.fastq.gz',
             'PJB1-B_ACAGTG_L002_R1_001.fastq.gz',))
        dirn = os.path.join(self.dirn,'PJB')
        project = AnalysisProject('PJB',dirn)
        # Pickle project
        pickled = cloudpickle.dumps(project)
        # Unpickle it
        unpickled = cloudpickle.loads(pickled)
        # Check the unpickled data
        self.assertEqual(unpickled.name,'PJB')
        self.assertTrue(os.path.isdir(unpickled.dirn))
        self.assertFalse(unpickled.multiple_fastqs)
        self.assertFalse(unpickled.info.paired_end)
        self.assertEqual(unpickled.info.primary_fastq_dir,'fastqs')
        self.assertEqual(unpickled.info.samples,'2 samples (PJB1-A, PJB1-B)')
        self.assertEqual(unpickled.samples[0].name,'PJB1-A')
        self.assertEqual(unpickled.samples[1].name,'PJB1-B')
        self.assertEqual(unpickled.fastq_dir,
                         os.path.join(unpickled.dirn,'fastqs'))
        self.assertEqual(unpickled.fastq_dirs,['fastqs',])

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

    def test_paired_end_analysis_sample_index_read_fastq(self):
        """Check AnalysisSample class with index read fastqs
        """
        sample = AnalysisSample('PJB1-B')
        fq_l1_r1 = '/run/sample1/PJB1-B_S1_L001_R1.fastq.gz'
        fq_l1_r2 = '/run/sample1/PJB1-B_S1_L001_R2.fastq.gz'
        fq_l1_i1 = '/run/sample1/PJB1-B_S1_L001_I1.fastq.gz'
        sample.add_fastq(fq_l1_r1)
        sample.add_fastq(fq_l1_r2)
        sample.add_fastq(fq_l1_i1)
        self.assertEqual(sample.name,'PJB1-B')
        self.assertEqual(sample.fastq,[fq_l1_i1,fq_l1_r1,fq_l1_r2])
        self.assertEqual(sample.fastq_subset(read_number=1),[fq_l1_r1,])
        self.assertEqual(sample.fastq_subset(read_number=2),[fq_l1_r2,])
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

class TestRunReferenceIdFunction(unittest.TestCase):
    """
    Tests for the 'run_reference_id' function
    """
    def test_run_reference_id(self):
        """run_reference_id: run name, platform and facility run number
        """
        self.assertEqual(run_reference_id("160621_M00879_0087_000000000-AGEW9",
                                          platform="miseq",
                                          facility_run_number=87),
                         "MISEQ_160621#87")
        self.assertEqual(run_reference_id("/data/160621_M00879_0087_000000000-AGEW9/",
                                          platform="miseq",
                                          facility_run_number=87),
                         "MISEQ_160621#87")

    def test_run_reference_id_no_platform(self):
        """run_reference_id: run name and facility run number (no platform)
        """
        self.assertEqual(run_reference_id("160621_M00879_0087_000000000-AGEW9",
                                          platform=None,
                                          facility_run_number=87),
                         "M00879_160621#87")
        self.assertEqual(run_reference_id("160621_M00879_0087_000000000-AGEW9",
                                          platform=None,
                                          facility_run_number=88),
                         "M00879_160621/87#88")
        self.assertEqual(run_reference_id("/data/160621_M00879_0087_000000000-AGEW9/",
                                          platform=None,
                                          facility_run_number=87),
                         "M00879_160621#87")

    def test_run_reference_id_no_facility_run_number(self):
        """run_reference_id: run name and platform (no facility run number)
        """
        self.assertEqual(run_reference_id("160621_M00879_0087_000000000-AGEW9",
                                          platform="miseq",
                                          facility_run_number=None),
                         "MISEQ_160621/87")

    def test_run_reference_id_facility_run_number_differs(self):
        """run_reference_id: instrument and facility run numbers differ
        """
        self.assertEqual(run_reference_id("160621_M00879_0087_000000000-AGEW9",
                                          platform="miseq",
                                          facility_run_number=90),
                         "MISEQ_160621/87#90")
        self.assertEqual(run_reference_id("/data/160621_M00879_0087_000000000-AGEW9/",
                                          platform="miseq",
                                          facility_run_number=90),
                         "MISEQ_160621/87#90")

    def test_run_reference_id_bad_run_name(self):
        """run_reference_id: handle 'bad' run name (cannot be split)
        """
        self.assertEqual(run_reference_id("rag_05_2017",
                                          platform=None,
                                          facility_run_number=None),
                         "rag_05_2017")
        self.assertEqual(run_reference_id("rag_05_2017",
                                          platform="miseq",
                                          facility_run_number=None),
                         "MISEQ_rag_05_2017")
        self.assertEqual(run_reference_id("rag_05_2017",
                                          platform=None,
                                          facility_run_number=90),
                         "rag_05_2017#90")
        self.assertEqual(run_reference_id("rag_05_2017",
                                          platform="miseq",
                                          facility_run_number=90),
                         "MISEQ_rag_05_2017#90")

class TestSplitSampleNameFunction(unittest.TestCase):
    """
    Tests for the 'split_sample_name' function
    """
    def test_split_sample_name(self):
        """split_sample_name: check names are split correctly
        """
        self.assertEqual(split_sample_name("PJB"),["PJB"])
        self.assertEqual(split_sample_name("PJB1"),["PJB",1])
        self.assertEqual(split_sample_name("PJB0001"),["PJB",1])
        self.assertEqual(split_sample_name("PJB_1-10"),["PJB_",1,"-",10])

class TestCopyAnalysisProject(unittest.TestCase):
    """
    Tests for the 'copy_analysis_project' function
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestCopyProject')

    def tearDown(self):
        # Remove the temporary test directory
        shutil.rmtree(self.wd)

    def test_copy_project(self):
        """
        copy_project: copies project instance
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Make initial project
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        # Make a copy
        project2 = copy_analysis_project(project)
        # Check copy
        self.assertEqual(project.name,project2.name)
        self.assertEqual(project.dirn,project2.dirn)
        self.assertEqual(project.fastq_dir,project2.fastq_dir)
        self.assertEqual(project.fastq_dirs,project2.fastq_dirs)
        self.assertEqual(project.fastqs,project2.fastqs)
        self.assertEqual(project.info.organism,project2.info.organism)
