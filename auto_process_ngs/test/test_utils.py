#######################################################################
# Tests for utils.py module
#######################################################################

import unittest
import tempfile
import shutil
from bcftbx.JobRunner import SimpleJobRunner,GEJobRunner
from auto_process_ngs.utils import *

class TestAnalysisFastq(unittest.TestCase):
    """Tests for the AnalysisFastq class

    """

    def test_full_name(self):
        """Handle full Illumina-style fastq name
        """
        fq = AnalysisFastq('NH1_ChIP-seq_Gli1_ACAGTG_L003_R2_001')
        self.assertEqual(fq.sample_name,'NH1_ChIP-seq_Gli1')
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

    def test_input_is_full_path(self):
        """Handle input as full path to Fastq file
        """
        fq = AnalysisFastq('/data/Project_NH/Sample_NH1/NH1_ChIP-seq_Gli1_ACAGTG_L003_R2_001.fastq.gz')
        self.assertEqual(fq.sample_name,'NH1_ChIP-seq_Gli1')
        self.assertEqual(fq.sample_number,None)
        self.assertEqual(fq.barcode_sequence,'ACAGTG')
        self.assertEqual(fq.lane_number,3)
        self.assertEqual(fq.read_number,2)
        self.assertEqual(fq.set_number,1)
        self.assertEqual(str(fq),'NH1_ChIP-seq_Gli1_ACAGTG_L003_R2_001')

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
        self.assertEqual(project.samples[0].name,'PJB1-A')
        self.assertEqual(project.samples[1].name,'PJB1-B')

    def test_create_single_end_analysis_project_multi_fastqs(self):
        """Check creation of new single-end AnalysisProject directory (multi-fastq/sample)
        """
        self.make_data_dir(('PJB1-B_ACAGTG_L001_R1_001.fastq.gz',
                            'PJB1-B_ACAGTG_L002_R1_001.fastq.gz',))
        dirn = os.path.join(self.dirn,'PJB')
        project = AnalysisProject('PJB',dirn)
        project.create_directory(fastqs=self.fastqs)
        self.assertEqual(project.name,'PJB')
        self.assertTrue(os.path.isdir(project.dirn))
        self.assertTrue(project.multiple_fastqs)
        self.assertFalse(project.info.paired_end)
        self.assertEqual(project.samples[0].name,'PJB1-B')

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
        self.assertEqual(project.samples[0].name,'PJB1-A')
        self.assertEqual(project.samples[1].name,'PJB1-B')

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

        
class TestBclToFastqInfo(unittest.TestCase):
    """
    Tests for the bcl_to_fastq_info function

    """
    def setUp(self):
        # Make some fake directories for different
        # software versions
        self.dirn = tempfile.mkdtemp(suffix='TestBclToFastqInfo')
        # No version
        os.makedirs(os.path.join(self.dirn,'unknown'))
        open(os.path.join(self.dirn,'unknown','configureBclToFastq.pl'),'w').write("")
        os.chmod(os.path.join(self.dirn,'unknown','configureBclToFastq.pl'),0775)
        # CASAVA 1.8.2
        os.makedirs(os.path.join(self.dirn,'casava','1.8.2','bin'))
        os.makedirs(os.path.join(self.dirn,'casava','1.8.2','etc','CASAVA-1.8.2'))
        open(os.path.join(self.dirn,'casava','1.8.2','bin',
                          'configureBclToFastq.pl'),'w').write("")
        os.chmod(os.path.join(self.dirn,'casava','1.8.2','bin',
                              'configureBclToFastq.pl'),0775)
        # bcl2fastq 1.8.3
        os.makedirs(os.path.join(self.dirn,'bcl2fastq','1.8.3','bin'))
        os.makedirs(os.path.join(self.dirn,'bcl2fastq','1.8.3','etc','bcl2fastq-1.8.3'))
        open(os.path.join(self.dirn,'bcl2fastq','1.8.3','bin',
                          'configureBclToFastq.pl'),'w').write("")
        os.chmod(os.path.join(self.dirn,'bcl2fastq','1.8.3','bin',
                              'configureBclToFastq.pl'),0775)
        # bcl2fastq 1.8.4
        os.makedirs(os.path.join(self.dirn,'bcl2fastq','1.8.4','bin'))
        os.makedirs(os.path.join(self.dirn,'bcl2fastq','1.8.4','etc','bcl2fastq-1.8.4'))
        open(os.path.join(self.dirn,'bcl2fastq','1.8.4','bin',
                          'configureBclToFastq.pl'),'w').write("")
        os.chmod(os.path.join(self.dirn,'bcl2fastq','1.8.4','bin',
                              'configureBclToFastq.pl'),0775)

    def tearDown(self):
        # Remove the temporary test directory
        shutil.rmtree(self.dirn)

    def test_casava_1_8_2(self):
        """
        Collect info for CASAVA 1.8.2

        """
        exe = os.path.join(self.dirn,'casava','1.8.2','bin','configureBclToFastq.pl')
        os.environ['PATH'] = os.path.dirname(exe)+':'+os.environ['PATH']
        path,name,version = bcl_to_fastq_info()
        self.assertEqual(path,exe)
        self.assertEqual(name,'CASAVA')
        self.assertEqual(version,'1.8.2')

    def test_bcl2fastq_1_8_3(self):
        """
        Collect info for bcl2fastq 1.8.3

        """
        exe = os.path.join(self.dirn,'bcl2fastq','1.8.3','bin','configureBclToFastq.pl')
        os.environ['PATH'] = os.path.dirname(exe)+':'+os.environ['PATH']
        path,name,version = bcl_to_fastq_info()
        self.assertEqual(path,exe)
        self.assertEqual(name,'bcl2fastq')
        self.assertEqual(version,'1.8.3')

    def test_bcl2fastq_1_8_4(self):
        """
        Collect info for bcl2fastq 1.8.4

        """
        exe = os.path.join(self.dirn,'bcl2fastq','1.8.4','bin','configureBclToFastq.pl')
        os.environ['PATH'] = os.path.dirname(exe)+':'+os.environ['PATH']
        path,name,version = bcl_to_fastq_info()
        self.assertEqual(path,exe)
        self.assertEqual(name,'bcl2fastq')
        self.assertEqual(version,'1.8.4')

    def test_configurebcltofastq_not_found(self):
        """
        Collect info when configureBclToFastq.pl is not found

        """
        path,name,version = bcl_to_fastq_info()
        self.assertEqual(path,None)
        self.assertEqual(name,None)
        self.assertEqual(version,None)

    def test_configurebcltofastq_no_package(self):
        """
        Collect info when there is no package info for configureBclToFastq.pl

        """
        exe = os.path.join(self.dirn,'unknown','configureBclToFastq.pl')
        os.environ['PATH'] = os.path.dirname(exe)+':'+os.environ['PATH']
        path,name,version = bcl_to_fastq_info()
        self.assertEqual(path,exe)
        self.assertEqual(name,None)
        self.assertEqual(version,None)
