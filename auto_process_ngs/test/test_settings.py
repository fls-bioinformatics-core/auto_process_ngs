#######################################################################
# Tests for settings.py module
#######################################################################

import unittest
import tempfile
import shutil
from bcftbx.JobRunner import SimpleJobRunner,GEJobRunner
from auto_process_ngs.settings import *

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

class TestSettings(unittest.TestCase):
    """Tests for the Settings class
    """
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestSettings')

    def tearDown(self):
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.dirn)

    def test_sample_settings_file(self):
        """Settings: load from sample file
        """
        sample_settings_file = os.path.join(get_config_dir(),
                                            'auto_process.ini.sample')
        self.assertTrue(os.path.isfile(sample_settings_file),
                        "Missing sample file %s" % sample_settings_file)
        s = Settings(sample_settings_file)
        # General settings
        self.assertTrue(isinstance(s.general.default_runner,SimpleJobRunner))
        self.assertEqual(s.general.max_concurrent_jobs,12)
        self.assertEqual(s.general.max_cores,None)
        self.assertEqual(s.general.max_batches,None)
        self.assertEqual(s.general.poll_interval,5.0)
        # Modulefiles
        self.assertEqual(s.modulefiles.make_fastqs,None)
        self.assertEqual(s.modulefiles.run_qc,None)
        self.assertEqual(s.modulefiles.publish_qc,None)
        self.assertEqual(s.modulefiles.process_icell8,None)
        self.assertEqual(s.modulefiles.bcl2fastq,None)
        self.assertEqual(s.modulefiles.bcl_convert,None)
        self.assertEqual(s.modulefiles.cellranger_mkfastq,None)
        self.assertEqual(s.modulefiles.cellranger_atac_mkfastq,None)
        self.assertEqual(s.modulefiles.cellranger_arc_mkfastq,None)
        self.assertEqual(s.modulefiles.spaceranger_mkfastq,None)
        self.assertEqual(s.modulefiles.fastqc,None)
        self.assertEqual(s.modulefiles.fastq_screen,None)
        self.assertEqual(s.modulefiles.cellranger,None)
        self.assertEqual(s.modulefiles.report_qc,None)
        # Conda
        self.assertEqual(s.conda.enable_conda,False)
        self.assertEqual(s.conda.env_dir,None)
        # BCL conversion software
        self.assertEqual(s.bcl_conversion.bcl_converter,
                         'bcl2fastq>=1.8.4')
        self.assertEqual(s.bcl_conversion.nprocessors,None)
        self.assertEqual(s.bcl_conversion.no_lane_splitting,False)
        self.assertEqual(s.bcl_conversion.create_empty_fastqs,False)
        # NextSeq-specific
        self.assertEqual(s.platform.nextseq.bcl_converter,'bcl2fastq>=2.0')
        self.assertEqual(s.platform.nextseq.no_lane_splitting,True)
        # Fastq_stats
        self.assertEqual(s.fastq_stats.nprocessors,1)
        # Job-specific runners
        self.assertTrue(isinstance(s.runners.bcl2fastq,SimpleJobRunner))
        self.assertTrue(isinstance(s.runners.qc,SimpleJobRunner))
        self.assertTrue(isinstance(s.runners.star,SimpleJobRunner))
        self.assertTrue(isinstance(s.runners.stats,SimpleJobRunner))
        # Archiving
        self.assertEqual(s.archive.dirn,None)
        self.assertEqual(s.archive.log,None)
        self.assertEqual(s.archive.group,None)
        self.assertEqual(s.archive.chmod,None)
        # QC reporting
        self.assertEqual(s.qc_web_server.dirn,None)
        self.assertEqual(s.qc_web_server.url,None)
        # Metadata
        self.assertEqual(s.metadata.default_data_source,None)

    def test_partial_settings_file(self):
        """Settings: load a partial auto_process.ini file
        """
        # Partial file
        partial_settings_file = os.path.join(self.dirn,
                                             "auto_process.ini")
        with open(partial_settings_file,'w') as s:
            s.write("""[fastq_stats]
nprocessors = 8
""")
        # Load settings
        s = Settings(partial_settings_file)
        # General settings
        self.assertTrue(isinstance(s.general.default_runner,SimpleJobRunner))
        self.assertEqual(s.general.max_concurrent_jobs,12)
        self.assertEqual(s.general.max_cores,None)
        self.assertEqual(s.general.max_batches,None)
        self.assertEqual(s.general.poll_interval,5.0)
        # Modulefiles
        self.assertEqual(s.modulefiles.make_fastqs,None)
        self.assertEqual(s.modulefiles.bcl2fastq,None)
        self.assertEqual(s.modulefiles.bcl_convert,None)
        self.assertEqual(s.modulefiles.cellranger_mkfastq,None)
        self.assertEqual(s.modulefiles.cellranger_atac_mkfastq,None)
        self.assertEqual(s.modulefiles.cellranger_arc_mkfastq,None)
        self.assertEqual(s.modulefiles.spaceranger_mkfastq,None)
        self.assertEqual(s.modulefiles.run_qc,None)
        # Conda
        self.assertEqual(s.conda.enable_conda,False)
        self.assertEqual(s.conda.env_dir,None)
        # Fastq_stats
        self.assertEqual(s.fastq_stats.nprocessors,8)

    def test_get_item(self):
        """Settings: get_item fetches a value
        """
        sample_settings_file = os.path.join(get_config_dir(),
                                            'auto_process.ini.sample')
        s = Settings(sample_settings_file)
        max_concurrent_jobs = s.general.max_concurrent_jobs
        self.assertEqual(s['general'].max_concurrent_jobs,
                         max_concurrent_jobs)

    def test_preserve_option_case(self):
        """Settings: case of option names is preserved
        """
        # Settings file
        settings_file = os.path.join(self.dirn,"auto_process.ini")
        with open(settings_file,'w') as s:
            s.write("""[sequencers]
SN7001251 = hiseq
""")
        # Load settings
        s = Settings(settings_file)
        # Check case of option
        self.assertTrue('SN7001251' in s.sequencers)
        self.assertEqual(s.sequencers['SN7001251']['platform'],'hiseq')
        self.assertEqual(s.sequencers.SN7001251.platform,'hiseq')

    def test_conda_settings(self):
        """Settings: check conda options are set correctly
        """
        # Settings file
        settings_file = os.path.join(self.dirn,"auto_process.ini")
        with open(settings_file,'w') as s:
            s.write("""[conda]
enable_conda = true
env_dir = /scratch/conda_envs
""")
        # Load settings
        s = Settings(settings_file)
        # Check conda settings
        self.assertTrue(s.conda.enable_conda)
        self.assertEqual(s.conda.env_dir,"/scratch/conda_envs")

    def test_conda_env_dir(self):
        """Settings: check conda env dir expands env variable
        """
        # Settings file
        settings_file = os.path.join(self.dirn,"auto_process.ini")
        with open(settings_file,'w') as s:
            s.write("""[conda]
enable_conda = true
env_dir = /scratch/$USER/conda_envs
""")
        # Load settings
        s = Settings(settings_file)
        # Check conda settings
        self.assertTrue(s.conda.enable_conda)
        self.assertEqual(s.conda.env_dir,os.path.join("/scratch",
                                                      os.environ['USER'],
                                                      "conda_envs"))

    def test_destination_definitions(self):
        """Settings: handle 'destination:...' sections
        """
        # Settings file
        settings_file = os.path.join(self.dirn,"auto_process.ini")
        with open(settings_file,'w') as s:
            s.write("""[destination:webserver]
directory = /mnt/www/data
subdir = random_bin
readme_template = README.webserver.template
url = https://our.data.com/data
include_downloader = true
include_qc_report = true
hard_links = true

[destination:local]
directory = /mnt/shared
subdir = run_id
""")
        # Load settings
        s = Settings(settings_file)
        # Check destination settings
        self.assertTrue('webserver' in s.destination)
        self.assertEqual(s.destination['webserver']['directory'],
                         '/mnt/www/data')
        self.assertEqual(s.destination['webserver']['subdir'],'random_bin')
        self.assertEqual(s.destination['webserver']['readme_template'],
                         'README.webserver.template')
        self.assertEqual(s.destination['webserver']['url'],
                         'https://our.data.com/data')
        self.assertEqual(s.destination['webserver']['include_downloader'],
                         True)
        self.assertEqual(s.destination['webserver']['include_qc_report'],
                         True)
        self.assertEqual(s.destination['webserver']['hard_links'],True)
        self.assertTrue('local' in s.destination)
        self.assertEqual(s.destination['local']['directory'],'/mnt/shared')
        self.assertEqual(s.destination['local']['subdir'],'run_id')
        self.assertEqual(s.destination['local']['readme_template'],None)
        self.assertEqual(s.destination['local']['url'],None)
        self.assertEqual(s.destination['local']['include_downloader'],False)
        self.assertEqual(s.destination['local']['include_qc_report'],False)
        self.assertEqual(s.destination['local']['hard_links'],False)

    def test_sequencer_definitions(self):
        """Settings: handle 'sequencer:...' sections
        """
        # Settings file
        settings_file = os.path.join(self.dirn,"auto_process.ini")
        with open(settings_file,'w') as s:
            s.write("""[sequencer:SN7001250]
platform = hiseq
model = "HiSeq 2500"

[sequencer:NB110920]
platform = nextseq
model = "NextSeq 500"
""")
        # Load settings
        s = Settings(settings_file)
        # Check sequencer settings
        self.assertTrue('SN7001250' in s.sequencers)
        self.assertEqual(s.sequencers['SN7001250']['platform'],'hiseq')
        self.assertEqual(s.sequencers['SN7001250']['model'],"HiSeq 2500")
        self.assertTrue('NB110920' in s.sequencers)
        self.assertEqual(s.sequencers['NB110920']['platform'],'nextseq')
        self.assertEqual(s.sequencers['NB110920']['model'],"NextSeq 500")

    def test_legacy_sequencer_definitions(self):
        """Settings: handle 'sequencers' section (no 'sequencer:...'s)
        """
        # Settings file
        settings_file = os.path.join(self.dirn,"auto_process.ini")
        with open(settings_file,'w') as s:
            s.write("""[sequencers]
SN7001250 = hiseq
NB110920 = nextseq
""")
        # Load settings
        s = Settings(settings_file)
        # Check sequencer settings
        self.assertTrue('SN7001250' in s.sequencers)
        self.assertEqual(s.sequencers['SN7001250']['platform'],'hiseq')
        self.assertEqual(s.sequencers['SN7001250']['model'],None)
        self.assertTrue('NB110920' in s.sequencers)
        self.assertEqual(s.sequencers['NB110920']['platform'],'nextseq')
        self.assertEqual(s.sequencers['NB110920']['model'],None)

    def test_mixed_sequencer_definitions(self):
        """Settings: handle mixture of 'sequencer:...' & 'sequencers' sections
        """
        # Settings file
        settings_file = os.path.join(self.dirn,"auto_process.ini")
        with open(settings_file,'w') as s:
            s.write("""[sequencers]
SN7001250 = hiseq
NB110920 = nextseq

[sequencer:K00129]
platform = hiseq4000
model = "HiSeq 4000"
""")
        # Load settings
        s = Settings(settings_file)
        # Check sequencer settings
        self.assertTrue('K00129' in s.sequencers)
        self.assertEqual(s.sequencers['K00129']['platform'],'hiseq4000')
        self.assertEqual(s.sequencers['K00129']['model'],"HiSeq 4000")
        self.assertTrue('SN7001250' in s.sequencers)
        self.assertEqual(s.sequencers['SN7001250']['platform'],'hiseq')
        self.assertEqual(s.sequencers['SN7001250']['model'],None)
        self.assertTrue('NB110920' in s.sequencers)
        self.assertEqual(s.sequencers['NB110920']['platform'],'nextseq')
        self.assertEqual(s.sequencers['NB110920']['model'],None)

    def test_sequencer_definitions_fails_if_platform_not_set(self):
        """Settings: fail to load if 'sequencer:...' section missing 'platform'
        """
        # Settings file
        settings_file = os.path.join(self.dirn,"auto_process.ini")
        with open(settings_file,'w') as s:
            s.write("""[sequencer:SN7001250]
model = "HiSeq 2500"
""")
        # Load settings
        self.assertRaises(Exception,
                          Settings,
                          settings_file)

    def test_organism_definitions(self):
        """Settings: handle 'organism:...' sections
        """
        # Settings file
        settings_file = os.path.join(self.dirn,"auto_process.ini")
        with open(settings_file,'w') as s:
            s.write("""[organism:human]
star_index = /data/hg38/star
bowtie_index = /data/hg38/bowtie
cellranger_reference = /data/10x/refdata-gex-GRCh38-2020-A
cellranger_premrna_reference = /data/10x/refdata-cellranger-GRCh38-1.0.1-pre_mrna
cellranger_atac_reference = /data/10x/refdata-cellranger-atac-GRCh38-2020-A-2.0.0
cellranger_arc_reference = /data/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0

[organism:mouse]
star_index = /data/mm10/star
bowtie_index = /data/mm10/bowtie
cellranger_reference = /data/10x/refdata-gex-mm10-2020-A
cellranger_atac_reference = /data/10x/refdata-cellranger-atac-mm10-2020-A-2.0.0
cellranger_arc_reference = /data/10x/refdata-cellranger-arc-mm10-2020-A-2.0.0
""")
        # Load settings
        s = Settings(settings_file)
        # Check organism settings
        self.assertTrue('human' in s.organisms)
        self.assertEqual(s.organisms['human']['star_index'],
                         '/data/hg38/star')
        self.assertEqual(s.organisms['human']['bowtie_index'],
                         '/data/hg38/bowtie')
        self.assertEqual(s.organisms['human']['cellranger_reference'],
                         '/data/10x/refdata-gex-GRCh38-2020-A')
        self.assertEqual(s.organisms['human']['cellranger_premrna_reference'],
                         '/data/10x/refdata-cellranger-GRCh38-1.0.1-pre_mrna')
        self.assertEqual(s.organisms['human']['cellranger_atac_reference'],
                         '/data/10x/refdata-cellranger-atac-GRCh38-2020-A-2.0.0')
        self.assertEqual(s.organisms['human']['cellranger_arc_reference'],
                         '/data/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0')
        self.assertTrue('mouse' in s.organisms)
        self.assertEqual(s.organisms['mouse']['star_index'],
                         '/data/mm10/star')
        self.assertEqual(s.organisms['mouse']['bowtie_index'],
                         '/data/mm10/bowtie')
        self.assertEqual(s.organisms['mouse']['cellranger_reference'],
                         '/data/10x/refdata-gex-mm10-2020-A')
        self.assertEqual(s.organisms['mouse']['cellranger_premrna_reference'],
                         None)
        self.assertEqual(s.organisms['mouse']['cellranger_atac_reference'],
                         '/data/10x/refdata-cellranger-atac-mm10-2020-A-2.0.0')
        self.assertEqual(s.organisms['mouse']['cellranger_arc_reference'],
                         '/data/10x/refdata-cellranger-arc-mm10-2020-A-2.0.0')

    def test_legacy_organism_definitions(self):
        """Settings: handle sections for specific indices (no 'organism:...' sections)
        """
        # Settings file
        settings_file = os.path.join(self.dirn,"auto_process.ini")
        with open(settings_file,'w') as s:
            s.write("""[fastq_strand_indexes]
human = /data/hg38/star
mouse = /data/mm10/star

[10xgenomics_transcriptomes]
human = /data/10x/refdata-gex-GRCh38-2020-A
mouse = /data/10x/refdata-gex-mm10-2020-A

[10xgenomics_premrna_references]
human = /data/10x/refdata-cellranger-GRCh38-1.0.1-pre_mrna

[10xgenomics_atac_genome_references]
human = /data/10x/refdata-cellranger-atac-GRCh38-2020-A-2.0.0
mouse = /data/10x/refdata-cellranger-atac-mm10-2020-A-2.0.0

[10xgenomics_multiome_references]
human = /data/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0
mouse = /data/10x/refdata-cellranger-arc-mm10-2020-A-2.0.0
""")
        # Load settings
        s = Settings(settings_file)
        # Check organism settings
        self.assertTrue('human' in s.organisms)
        self.assertEqual(s.organisms['human']['star_index'],
                         '/data/hg38/star')
        self.assertEqual(s.organisms['human']['bowtie_index'],None)
        self.assertEqual(s.organisms['human']['cellranger_reference'],
                         '/data/10x/refdata-gex-GRCh38-2020-A')
        self.assertEqual(s.organisms['human']['cellranger_premrna_reference'],
                         '/data/10x/refdata-cellranger-GRCh38-1.0.1-pre_mrna')
        self.assertEqual(s.organisms['human']['cellranger_atac_reference'],
                         '/data/10x/refdata-cellranger-atac-GRCh38-2020-A-2.0.0')
        self.assertEqual(s.organisms['human']['cellranger_arc_reference'],
                         '/data/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0')
        self.assertTrue('mouse' in s.organisms)
        self.assertEqual(s.organisms['mouse']['star_index'],
                         '/data/mm10/star')
        self.assertEqual(s.organisms['mouse']['bowtie_index'],None)
        self.assertEqual(s.organisms['mouse']['cellranger_reference'],
                         '/data/10x/refdata-gex-mm10-2020-A')
        self.assertEqual(s.organisms['mouse']['cellranger_premrna_reference'],
                         None)
        self.assertEqual(s.organisms['mouse']['cellranger_atac_reference'],
                         '/data/10x/refdata-cellranger-atac-mm10-2020-A-2.0.0')
        self.assertEqual(s.organisms['mouse']['cellranger_arc_reference'],
                         '/data/10x/refdata-cellranger-arc-mm10-2020-A-2.0.0')

    def test_mixed_organism_definitions(self):
        """Settings: handle mixture of 'organism:...' and index sections
        """
        # Settings file
        settings_file = os.path.join(self.dirn,"auto_process.ini")
        with open(settings_file,'w') as s:
            s.write("""[organism:human]
bowtie_index = /data/hg38/bowtie

[organism:mouse]
star_index = /data/mm10/star
bowtie_index = /data/mm10/bowtie
cellranger_reference = /data/10x/refdata-gex-mm10-2020-A
cellranger_atac_reference = /data/10x/refdata-cellranger-atac-mm10-2020-A-2.0.0
cellranger_arc_reference = /data/10x/refdata-cellranger-arc-mm10-2020-A-2.0.0

[fastq_strand_indexes]
human = /data/hg38/star
mouse = /data/mm10/star

[10xgenomics_transcriptomes]
human = /data/10x/refdata-gex-GRCh38-2020-A

[10xgenomics_premrna_references]
human = /data/10x/refdata-cellranger-GRCh38-1.0.1-pre_mrna

[10xgenomics_atac_genome_references]
human = /data/10x/refdata-cellranger-atac-GRCh38-2020-A-2.0.0

[10xgenomics_multiome_references]
human = /data/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0
""")
        # Load settings
        s = Settings(settings_file)
        # Check organism settings
        self.assertTrue('human' in s.organisms)
        self.assertEqual(s.organisms['human']['star_index'],
                         '/data/hg38/star')
        self.assertEqual(s.organisms['human']['bowtie_index'],
                         '/data/hg38/bowtie')
        self.assertEqual(s.organisms['human']['cellranger_reference'],
                         '/data/10x/refdata-gex-GRCh38-2020-A')
        self.assertEqual(s.organisms['human']['cellranger_premrna_reference'],
                         '/data/10x/refdata-cellranger-GRCh38-1.0.1-pre_mrna')
        self.assertEqual(s.organisms['human']['cellranger_atac_reference'],
                         '/data/10x/refdata-cellranger-atac-GRCh38-2020-A-2.0.0')
        self.assertEqual(s.organisms['human']['cellranger_arc_reference'],
                         '/data/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0')
        self.assertTrue('mouse' in s.organisms)
        self.assertEqual(s.organisms['mouse']['star_index'],
                         '/data/mm10/star')
        self.assertEqual(s.organisms['mouse']['bowtie_index'],
                         '/data/mm10/bowtie')
        self.assertEqual(s.organisms['mouse']['cellranger_reference'],
                         '/data/10x/refdata-gex-mm10-2020-A')
        self.assertEqual(s.organisms['mouse']['cellranger_premrna_reference'],
                         None)
        self.assertEqual(s.organisms['mouse']['cellranger_atac_reference'],
                         '/data/10x/refdata-cellranger-atac-mm10-2020-A-2.0.0')
        self.assertEqual(s.organisms['mouse']['cellranger_arc_reference'],
                         '/data/10x/refdata-cellranger-arc-mm10-2020-A-2.0.0')

    def test_screen_definitions(self):
        """Settings: handle 'screen:...' sections
        """
        # Settings file
        settings_file = os.path.join(self.dirn,"auto_process.ini")
        with open(settings_file,'w') as s:
            s.write("""[qc]
fastq_screens = model_organisms,other_organisms

[screen:other_organisms]
conf_file = /data/fastq_screen/other_organisms.conf

[screen:model_organisms]
conf_file = /data/fastq_screen/model_organisms.conf
""")
        # Load settings
        s = Settings(settings_file)
        # Check screen settings
        self.assertTrue('model_organisms' in s.screens)
        self.assertEqual(s.screens['model_organisms']['conf_file'],
                         '/data/fastq_screen/model_organisms.conf')
        self.assertEqual(s.screens['other_organisms']['conf_file'],
                         '/data/fastq_screen/other_organisms.conf')
        self.assertEqual(s.qc.fastq_screens,
                         "model_organisms,other_organisms")

    def test_legacy_fastq_screen_naming_setting(self):
        """Settings: handle 'use_legacy_screen_names' setting
        """
        # Settings file without 'use_legacy_screen_names'
        settings_file = os.path.join(self.dirn,"auto_process.ini")
        with open(settings_file,'w') as s:
            s.write("""[qc]
#use_legacy_screen_names = False
""")
        # Load and check settings
        s = Settings(settings_file)
        self.assertEqual(s.qc.use_legacy_screen_names,False)
        # Settings file with 'use_legacy_screen_names' turned off
        settings_file = os.path.join(self.dirn,"auto_process.ini")
        with open(settings_file,'w') as s:
            s.write("""[qc]
use_legacy_screen_names = False
""")
        # Load and check settings
        s = Settings(settings_file)
        self.assertEqual(s.qc.use_legacy_screen_names,False)
        # Settings file with 'use_legacy_screen_names' turned on
        settings_file = os.path.join(self.dirn,"auto_process.ini")
        with open(settings_file,'w') as s:
            s.write("""[qc]
use_legacy_screen_names = True
""")
        # Load settings
        s = Settings(settings_file)
        # Load and check settings
        s = Settings(settings_file)
        self.assertEqual(s.qc.use_legacy_screen_names,True)

    def test_legacy_illumina_qc_modulefile_setting(self):
        """Settings: handle legacy 'illumina_qc' modulefile setting
        """
        # Settings file with 'illumina_qc' only
        settings_file = os.path.join(self.dirn,"auto_process.ini")
        with open(settings_file,'w') as s:
            s.write("""[modulefiles]
illumina_qc = SimpleJobRunner(nslots=12)
""")
        # Load settings
        s = Settings(settings_file)
        # Check fastqc and fastq_screen modulefile settings
        self.assertEqual(str(s.modulefiles.fastqc),
                         'SimpleJobRunner(nslots=12)')
        self.assertEqual(str(s.modulefiles.fastq_screen),
                         'SimpleJobRunner(nslots=12)')
        # Settings file with 'illumina_qc' and 'fastqc' and
        # 'fastq_screen'
        with open(settings_file,'w') as s:
            s.write("""[modulefiles]
illumina_qc = SimpleJobRunner(nslots=12)
fastqc = SimpleJobRunner(nslots=1)
fastq_screen = SimpleJobRunner(nslots=8)
""")
        # Load settings
        s = Settings(settings_file)
        # Check fastqc and fastq_screen modulefile settings
        self.assertEqual(str(s.modulefiles.fastqc),
                         'SimpleJobRunner(nslots=1)')
        self.assertEqual(str(s.modulefiles.fastq_screen),
                         'SimpleJobRunner(nslots=8)')

    def test_legacy_bcl2fastq_settings_no_bcl_conversion(self):
        """Settings: handle legacy 'bcl2fastq' section (no 'bcl_conversion' settings)
        """
        # Settings file
        settings_file = os.path.join(self.dirn,"auto_process.ini")
        with open(settings_file,'w') as s:
            s.write("""[bcl2fastq]
default_version = >=2.20
nprocessors = 8
no_lane_splitting = true
create_empty_fastqs = false
""")
        # Load settings
        s = Settings(settings_file)
        # Check bcl_conversion settings
        self.assertEqual(s.bcl_conversion.bcl_converter,
                         'bcl2fastq>=2.20')
        self.assertEqual(s.bcl_conversion.nprocessors,8)
        self.assertEqual(s.bcl_conversion.no_lane_splitting,True)
        self.assertEqual(s.bcl_conversion.create_empty_fastqs,False)

    def test_legacy_platform_settings_no_bcl_conversion(self):
        """Settings: handle legacy 'platform:...' section (no 'bcl_conversion' settings)
        """
        # Settings file
        settings_file = os.path.join(self.dirn,"auto_process.ini")
        with open(settings_file,'w') as s:
            s.write("""[platform:nextseq]
bcl2fastq = >=2.20
nprocessors = 1
""")
        # Load settings
        s = Settings(settings_file)
        # Check bcl_conversion settings (should be defaults)
        self.assertEqual(s.bcl_conversion.bcl_converter,None)
        self.assertEqual(s.bcl_conversion.nprocessors,None)
        self.assertEqual(s.bcl_conversion.no_lane_splitting,None)
        self.assertEqual(s.bcl_conversion.create_empty_fastqs,None)
        # Check platform-specific options
        self.assertEqual(s.platform['nextseq'].bcl_converter,
                         'bcl2fastq>=2.20')
        self.assertEqual(s.platform['nextseq'].nprocessors,1)

    def test_platform_settings_override_bcl_conversion_section(self):
        """Settings: 'platform:...' section overrides 'bcl_conversion' settings
        """
        # Settings file
        settings_file = os.path.join(self.dirn,"auto_process.ini")
        with open(settings_file,'w') as s:
            s.write("""[bcl_conversion]
bcl_converter = bcl-convert=3.7.5
nprocessors = 16

[platform:nextseq]
bcl_converter = bcl2fastq>=2.20
nprocessors = 8
no_lane_splitting = true
""")
        # Load settings
        s = Settings(settings_file)
        # Check bcl_conversion settings
        self.assertEqual(s.bcl_conversion.bcl_converter,
                         'bcl-convert=3.7.5')
        self.assertEqual(s.bcl_conversion.nprocessors,16)
        self.assertEqual(s.bcl_conversion.no_lane_splitting,None)
        self.assertEqual(s.bcl_conversion.create_empty_fastqs,None)
        # Check platform-specific options
        self.assertEqual(s.platform['nextseq'].bcl_converter,
                         'bcl2fastq>=2.20')
        self.assertEqual(s.platform['nextseq'].nprocessors,8)
        self.assertEqual(s.platform['nextseq'].no_lane_splitting,True)
        self.assertEqual(s.platform['nextseq'].create_empty_fastqs,None)
