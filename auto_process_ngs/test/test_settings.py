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
                         'bcl2fastq>=2.20')
        self.assertEqual(s.bcl_conversion.nprocessors,None)
        self.assertEqual(s.bcl_conversion.no_lane_splitting,False)
        self.assertEqual(s.bcl_conversion.create_empty_fastqs,False)
        # Fastq_stats
        self.assertEqual(s.fastq_stats.nprocessors,None)
        # Job-specific runners
        self.assertTrue(isinstance(s.runners.bcl2fastq,SimpleJobRunner))
        self.assertTrue(isinstance(s.runners.bcl_convert,SimpleJobRunner))
        self.assertTrue(isinstance(s.runners.barcode_analysis,SimpleJobRunner))
        self.assertTrue(isinstance(s.runners.fastqc,SimpleJobRunner))
        self.assertTrue(isinstance(s.runners.fastq_screen,SimpleJobRunner))
        self.assertTrue(isinstance(s.runners.merge_fastqs,SimpleJobRunner))
        self.assertTrue(isinstance(s.runners.picard,SimpleJobRunner))
        self.assertTrue(isinstance(s.runners.qualimap,SimpleJobRunner))
        self.assertTrue(isinstance(s.runners.rseqc,SimpleJobRunner))
        self.assertTrue(isinstance(s.runners.rsync,SimpleJobRunner))
        self.assertTrue(isinstance(s.runners.star,SimpleJobRunner))
        self.assertTrue(isinstance(s.runners.stats,SimpleJobRunner))
        self.assertTrue(isinstance(s.runners.cellranger_count,SimpleJobRunner))
        self.assertTrue(isinstance(s.runners.cellranger_mkfastq,
                                   SimpleJobRunner))
        self.assertTrue(isinstance(s.runners.cellranger_multi,SimpleJobRunner))
        # Legacy runners no longer in config file
        self.assertTrue(isinstance(s.runners.cellranger,SimpleJobRunner))
        self.assertTrue(isinstance(s.runners.qc,SimpleJobRunner))
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

    def test_contains(self):
        """Settings: 'in' works for checking if section is defined
        """
        # Partial file
        partial_settings_file = os.path.join(self.dirn,
                                             "auto_process.ini")
        with open(partial_settings_file,'w') as s:
            s.write("""[general]
max_concurrent_jobs = 12

[organism:human]
star_index = /data/hg38/star_index
""")
        # Load settings
        s = Settings(partial_settings_file)
        # Check that sections are located
        self.assertTrue("general" in s)
        self.assertTrue("organisms" in s)
        # Check that section with subsection is located
        self.assertTrue("organism:human" in s)
        # Check that parameters can be located
        self.assertTrue("general.max_concurrent_jobs" in s)
        self.assertTrue("organism:human.star_index" in s)
        # Check that missing sections, subsections and parameters
        # are not located
        self.assertFalse("missing" in s)
        self.assertFalse("organism:missing" in s)
        self.assertFalse("organism:human.missing" in s)

    def test_add_section(self):
        """Settings: add_section creates new section
        """
        # Empty file
        empty_settings_file = os.path.join(self.dirn,
                                           "auto_process.ini")
        with open(empty_settings_file,'wt') as s:
            s.write("")
        s = Settings(empty_settings_file)
        # Arbitrary section name
        self.assertFalse("new_section" in s)
        s.add_section("new_section")
        self.assertTrue("new_section" in s)
        # Arbitrary section and subsection
        self.assertFalse("new_section2" in s)
        self.assertFalse("new_section2:subsection" in s)
        s.add_section("new_section2:subsection")
        self.assertTrue("new_section2" in s)
        self.assertTrue("new_section2:subsection" in s)
        # Special case: organism
        self.assertFalse("organism:human" in s)
        s.add_section("organism:human")
        self.assertTrue("organism:human" in s)

    def test_set_item(self):
        """Settings: 'set' updates a value
        """
        # Partial file
        partial_settings_file = os.path.join(self.dirn,
                                             "auto_process.ini")
        with open(partial_settings_file,'wt') as s:
            s.write("""[general]
max_concurrent_jobs = None

[organism:human]
""")
        # Load settings
        s = Settings(partial_settings_file)
        s.set("general.max_concurrent_jobs",8)
        self.assertEqual(s['general'].max_concurrent_jobs,8)
        s.set("organism:human.star_index","/data/hg38/star_index")
        self.assertEqual(s['organisms']['human']['star_index'],
                         "/data/hg38/star_index")

    def test_save_raw(self):
        """
        Settings: test saving "raw" settings file
        """
        # Settings content
        content = """[general]
default_runner = SimpleJobRunner(join_logs=True)
max_concurrent_jobs = 12
poll_interval = 5.0

[modulefiles]
bcl2fastq = apps/bcl2fastq/2.20.0.422
cellranger = apps/cellranger/6.1.2

[conda]
enable_conda = True
env_dir = $HOME/qc_conda_envs

[organism:human]
star_index = /data/indexes/hg38

[organism:mouse]
star_index = /data/indexes/mm10

[sequencer:A01234]
platform = novaseq6000
model = NovaSeq 6000

[sequencer:NB543201]
platform = nextseq
model = NextSeq 500

[runners]
bcl2fastq = SimpleJobRunner(nslots=8 join_logs=True)
star = SimpleJobRunner(nslots=18 join_logs=True)
"""
        # Settings file
        settings_file = os.path.join(self.dirn,"auto_process.ini")
        with open(settings_file,'wt') as s:
            s.write(content)
        # Load settings
        s = Settings(settings_file,resolve_undefined=False)
        # Save to new file
        saved_settings_file = os.path.join(self.dirn,"auto_process.ini.bak")
        s.save(out_file=saved_settings_file)
        # Compare with original
        self.assertTrue(os.path.exists(saved_settings_file))
        self.assertEqual(open(settings_file,'rt').read().rstrip(),
                         open(saved_settings_file,'rt').read().rstrip())

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

[destination:zips]
directory = /mnt/www/zipped
subdir = run_id
zip_fastqs = true
max_zip_size = 5G
""")
        # Load settings
        s = Settings(settings_file)
        # Check destination settings
        self.assertTrue('webserver' in s.destination)
        self.assertEqual(s.destination['webserver']['directory'],
                         '/mnt/www/data')
        self.assertEqual(s.destination['webserver']['subdir'],'random_bin')
        self.assertEqual(s.destination['webserver']['zip_fastqs'],False)
        self.assertEqual(s.destination['webserver']['max_zip_size'],None)
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
        self.assertEqual(s.destination['local']['zip_fastqs'],False)
        self.assertEqual(s.destination['local']['max_zip_size'],None)
        self.assertEqual(s.destination['local']['readme_template'],None)
        self.assertEqual(s.destination['local']['url'],None)
        self.assertEqual(s.destination['local']['include_downloader'],False)
        self.assertEqual(s.destination['local']['include_qc_report'],False)
        self.assertEqual(s.destination['local']['hard_links'],False)
        self.assertTrue('zips' in s.destination)
        self.assertEqual(s.destination['zips']['directory'],'/mnt/www/zipped')
        self.assertEqual(s.destination['zips']['subdir'],'run_id')
        self.assertEqual(s.destination['zips']['zip_fastqs'],True)
        self.assertEqual(s.destination['zips']['max_zip_size'],'5G')
        self.assertEqual(s.destination['zips']['readme_template'],None)
        self.assertEqual(s.destination['zips']['url'],None)
        self.assertEqual(s.destination['zips']['include_downloader'],False)
        self.assertEqual(s.destination['zips']['include_qc_report'],False)
        self.assertEqual(s.destination['zips']['hard_links'],False)

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
annotation_bed = /data/hg38/annotation/hg38.bed
annotation_gtf = /data/hg38/annotation/hg38.gtf
cellranger_reference = /data/10x/refdata-gex-GRCh38-2020-A
cellranger_premrna_reference = /data/10x/refdata-cellranger-GRCh38-1.0.1-pre_mrna
cellranger_atac_reference = /data/10x/refdata-cellranger-atac-GRCh38-2020-A-2.0.0
cellranger_arc_reference = /data/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0
cellranger_probe_set = /data/10x/Probe_Set_v1.0_GRCh38-2020-A.csv

[organism:mouse]
star_index = /data/mm10/star
bowtie_index = /data/mm10/bowtie
annotation_bed = /data/mm10/annotation/mm10.bed
annotation_gtf = /data/mm10/annotation/mm10.gtf
cellranger_reference = /data/10x/refdata-gex-mm10-2020-A
cellranger_atac_reference = /data/10x/refdata-cellranger-atac-mm10-2020-A-2.0.0
cellranger_arc_reference = /data/10x/refdata-cellranger-arc-mm10-2020-A-2.0.0
cellranger_probe_set = /data/10x/Probe_Set_v1.0_mm10-2020-A.csv
""")
        # Load settings
        s = Settings(settings_file)
        # Check organism settings
        self.assertTrue('human' in s.organisms)
        self.assertEqual(s.organisms['human']['star_index'],
                         '/data/hg38/star')
        self.assertEqual(s.organisms['human']['bowtie_index'],
                         '/data/hg38/bowtie')
        self.assertEqual(s.organisms['human']['annotation_bed'],
                         '/data/hg38/annotation/hg38.bed')
        self.assertEqual(s.organisms['human']['annotation_gtf'],
                         '/data/hg38/annotation/hg38.gtf')
        self.assertEqual(s.organisms['human']['cellranger_reference'],
                         '/data/10x/refdata-gex-GRCh38-2020-A')
        self.assertEqual(s.organisms['human']['cellranger_premrna_reference'],
                         '/data/10x/refdata-cellranger-GRCh38-1.0.1-pre_mrna')
        self.assertEqual(s.organisms['human']['cellranger_atac_reference'],
                         '/data/10x/refdata-cellranger-atac-GRCh38-2020-A-2.0.0')
        self.assertEqual(s.organisms['human']['cellranger_arc_reference'],
                         '/data/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0')
        self.assertEqual(s.organisms['human']['cellranger_probe_set'],
                         '/data/10x/Probe_Set_v1.0_GRCh38-2020-A.csv')
        self.assertTrue('mouse' in s.organisms)
        self.assertEqual(s.organisms['mouse']['star_index'],
                         '/data/mm10/star')
        self.assertEqual(s.organisms['mouse']['bowtie_index'],
                         '/data/mm10/bowtie')
        self.assertEqual(s.organisms['mouse']['annotation_bed'],
                         '/data/mm10/annotation/mm10.bed')
        self.assertEqual(s.organisms['mouse']['annotation_gtf'],
                         '/data/mm10/annotation/mm10.gtf')
        self.assertEqual(s.organisms['mouse']['cellranger_reference'],
                         '/data/10x/refdata-gex-mm10-2020-A')
        self.assertEqual(s.organisms['mouse']['cellranger_premrna_reference'],
                         None)
        self.assertEqual(s.organisms['mouse']['cellranger_atac_reference'],
                         '/data/10x/refdata-cellranger-atac-mm10-2020-A-2.0.0')
        self.assertEqual(s.organisms['mouse']['cellranger_arc_reference'],
                         '/data/10x/refdata-cellranger-arc-mm10-2020-A-2.0.0')
        self.assertEqual(s.organisms['mouse']['cellranger_probe_set'],
                         '/data/10x/Probe_Set_v1.0_mm10-2020-A.csv')

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
        self.assertEqual(s.organisms['human']['annotation_bed'],None)
        self.assertEqual(s.organisms['human']['annotation_gtf'],None)
        self.assertEqual(s.organisms['human']['cellranger_reference'],
                         '/data/10x/refdata-gex-GRCh38-2020-A')
        self.assertEqual(s.organisms['human']['cellranger_premrna_reference'],
                         '/data/10x/refdata-cellranger-GRCh38-1.0.1-pre_mrna')
        self.assertEqual(s.organisms['human']['cellranger_atac_reference'],
                         '/data/10x/refdata-cellranger-atac-GRCh38-2020-A-2.0.0')
        self.assertEqual(s.organisms['human']['cellranger_arc_reference'],
                         '/data/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0')
        self.assertEqual(s.organisms['human']['cellranger_probe_set'],None)
        self.assertTrue('mouse' in s.organisms)
        self.assertEqual(s.organisms['mouse']['star_index'],
                         '/data/mm10/star')
        self.assertEqual(s.organisms['mouse']['bowtie_index'],None)
        self.assertEqual(s.organisms['mouse']['annotation_bed'],None)
        self.assertEqual(s.organisms['mouse']['annotation_gtf'],None)
        self.assertEqual(s.organisms['mouse']['cellranger_reference'],
                         '/data/10x/refdata-gex-mm10-2020-A')
        self.assertEqual(s.organisms['mouse']['cellranger_premrna_reference'],
                         None)
        self.assertEqual(s.organisms['mouse']['cellranger_atac_reference'],
                         '/data/10x/refdata-cellranger-atac-mm10-2020-A-2.0.0')
        self.assertEqual(s.organisms['mouse']['cellranger_arc_reference'],
                         '/data/10x/refdata-cellranger-arc-mm10-2020-A-2.0.0')
        self.assertEqual(s.organisms['mouse']['cellranger_probe_set'],None)

    def test_mixed_organism_definitions(self):
        """Settings: handle mixture of 'organism:...' and index sections
        """
        # Settings file
        settings_file = os.path.join(self.dirn,"auto_process.ini")
        with open(settings_file,'w') as s:
            s.write("""[organism:human]
bowtie_index = /data/hg38/bowtie
annotation_bed = /data/hg38/annotation/hg38.bed
annotation_gtf = /data/hg38/annotation/hg38.gtf

[organism:mouse]
star_index = /data/mm10/star
bowtie_index = /data/mm10/bowtie
annotation_bed = /data/mm10/annotation/mm10.bed
annotation_gtf = /data/mm10/annotation/mm10.gtf
cellranger_reference = /data/10x/refdata-gex-mm10-2020-A
cellranger_atac_reference = /data/10x/refdata-cellranger-atac-mm10-2020-A-2.0.0
cellranger_arc_reference = /data/10x/refdata-cellranger-arc-mm10-2020-A-2.0.0
cellranger_probe_set = /data/10x/Probe_Set_v1.0_mm10-2020-A.csv

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
        self.assertEqual(s.organisms['human']['annotation_bed'],
                         '/data/hg38/annotation/hg38.bed')
        self.assertEqual(s.organisms['human']['annotation_gtf'],
                         '/data/hg38/annotation/hg38.gtf')
        self.assertEqual(s.organisms['human']['cellranger_reference'],
                         '/data/10x/refdata-gex-GRCh38-2020-A')
        self.assertEqual(s.organisms['human']['cellranger_premrna_reference'],
                         '/data/10x/refdata-cellranger-GRCh38-1.0.1-pre_mrna')
        self.assertEqual(s.organisms['human']['cellranger_atac_reference'],
                         '/data/10x/refdata-cellranger-atac-GRCh38-2020-A-2.0.0')
        self.assertEqual(s.organisms['human']['cellranger_arc_reference'],
                         '/data/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0')
        self.assertEqual(s.organisms['human']['cellranger_probe_set'],None)
        self.assertTrue('mouse' in s.organisms)
        self.assertEqual(s.organisms['mouse']['star_index'],
                         '/data/mm10/star')
        self.assertEqual(s.organisms['mouse']['bowtie_index'],
                         '/data/mm10/bowtie')
        self.assertEqual(s.organisms['mouse']['annotation_bed'],
                         '/data/mm10/annotation/mm10.bed')
        self.assertEqual(s.organisms['mouse']['annotation_gtf'],
                         '/data/mm10/annotation/mm10.gtf')
        self.assertEqual(s.organisms['mouse']['cellranger_reference'],
                         '/data/10x/refdata-gex-mm10-2020-A')
        self.assertEqual(s.organisms['mouse']['cellranger_premrna_reference'],
                         None)
        self.assertEqual(s.organisms['mouse']['cellranger_atac_reference'],
                         '/data/10x/refdata-cellranger-atac-mm10-2020-A-2.0.0')
        self.assertEqual(s.organisms['mouse']['cellranger_arc_reference'],
                         '/data/10x/refdata-cellranger-arc-mm10-2020-A-2.0.0')
        self.assertEqual(s.organisms['mouse']['cellranger_probe_set'],
                         '/data/10x/Probe_Set_v1.0_mm10-2020-A.csv')

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
            illumina_qc = apps/fastqc/0.11.3,apps/fastq_screen/0.14.0
""")
        # Load settings
        s = Settings(settings_file)
        # Check fastqc and fastq_screen modulefile settings
        self.assertEqual(str(s.modulefiles.fastqc),
                         'apps/fastqc/0.11.3,apps/fastq_screen/0.14.0')
        self.assertEqual(str(s.modulefiles.fastq_screen),
                         'apps/fastqc/0.11.3,apps/fastq_screen/0.14.0')
        # Settings file with 'illumina_qc' and 'fastqc' and
        # 'fastq_screen'
        with open(settings_file,'w') as s:
            s.write("""[modulefiles]
illumina_qc = apps/fastqc/0.11.3,apps/fastq_screen/0.14.0
fastqc = apps/fastqc/0.11.3
fastq_screen = apps/fastq_screen/0.14.0
""")
        # Load settings
        s = Settings(settings_file)
        # Check fastqc and fastq_screen modulefile settings
        self.assertEqual(str(s.modulefiles.fastqc),'apps/fastqc/0.11.3')
        self.assertEqual(str(s.modulefiles.fastq_screen),
                         'apps/fastq_screen/0.14.0')

    def test_legacy_runner_settings(self):
        """Settings: handle legacy runner settings
        """
        # Settings file with 'qc' and 'cellranger' runners only
        settings_file = os.path.join(self.dirn,"auto_process.ini")
        with open(settings_file,'w') as s:
            s.write("""[general]
default_runner = SimpleJobRunner

[runners]
qc = SimpleJobRunner(nslots=4)
cellranger = SimpleJobRunner(nslots=12)
""")
        # Load settings
        s = Settings(settings_file)
        # Check fastqc and fastq_screen runner settings
        self.assertEqual(str(s.runners.qc),
                         'SimpleJobRunner(nslots=4 join_logs=True)')
        self.assertEqual(str(s.runners.fastqc),
                         'SimpleJobRunner(nslots=4 join_logs=True)')
        self.assertEqual(str(s.runners.fastq_screen),
                         'SimpleJobRunner(nslots=4 join_logs=True)')
        self.assertEqual(str(s.runners.cellranger_mkfastq),
                         'SimpleJobRunner(nslots=12 join_logs=True)')
        self.assertEqual(str(s.runners.cellranger_count),
                         'SimpleJobRunner(nslots=12 join_logs=True)')
        self.assertEqual(str(s.runners.cellranger_multi),
                         'SimpleJobRunner(nslots=12 join_logs=True)')
        # Settings file with 'qc' and 'fastqc' and
        # 'fastq_screen'
        with open(settings_file,'w') as s:
            s.write("""[general]
default_runner = SimpleJobRunner

[runners]
qc = SimpleJobRunner(nslots=4)
fastqc = SimpleJobRunner(nslots=1)
fastq_screen = SimpleJobRunner(nslots=8)
cellranger_mkfastq = SimpleJobRunner(nslots=8)
cellranger_count = SimpleJobRunner(nslots=16)
cellranger_multi = SimpleJobRunner(nslots=12)
""")
        # Load settings
        s = Settings(settings_file)
        # Check fastqc and fastq_screen runner settings
        self.assertEqual(str(s.runners.fastqc),
                         'SimpleJobRunner(join_logs=True)')
        self.assertEqual(str(s.runners.fastq_screen),
                         'SimpleJobRunner(nslots=8 join_logs=True)')
        self.assertEqual(str(s.runners.cellranger_mkfastq),
                         'SimpleJobRunner(nslots=8 join_logs=True)')
        self.assertEqual(str(s.runners.cellranger_count),
                         'SimpleJobRunner(nslots=16 join_logs=True)')
        self.assertEqual(str(s.runners.cellranger_multi),
                         'SimpleJobRunner(nslots=12 join_logs=True)')

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
nprocessors = 8
""")
        # Load settings
        s = Settings(settings_file)
        # Check bcl_conversion settings (should be defaults)
        self.assertEqual(s.bcl_conversion.bcl_converter,
                         'bcl2fastq>=2.20')
        self.assertEqual(s.bcl_conversion.nprocessors,None)
        self.assertEqual(s.bcl_conversion.no_lane_splitting,False)
        self.assertEqual(s.bcl_conversion.create_empty_fastqs,False)
        # Check platform-specific options
        self.assertEqual(s.platform['nextseq'].bcl_converter,
                         'bcl2fastq>=2.20')
        self.assertEqual(s.platform['nextseq'].nprocessors,8)

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
        self.assertEqual(s.bcl_conversion.no_lane_splitting,False)
        self.assertEqual(s.bcl_conversion.create_empty_fastqs,False)
        # Check platform-specific options
        self.assertEqual(s.platform['nextseq'].bcl_converter,
                         'bcl2fastq>=2.20')
        self.assertEqual(s.platform['nextseq'].nprocessors,8)
        self.assertEqual(s.platform['nextseq'].no_lane_splitting,True)
        self.assertEqual(s.platform['nextseq'].create_empty_fastqs,None)

class TestFetchReferenceData(unittest.TestCase):
    """
    Tests for the 'fetch_reference_data' function
    """
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestFetchRefData')

    def tearDown(self):
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.dirn)

    def test_fetch_reference_data(self):
        """
        fetch_reference_data: check expected data are returned
        """
        # Settings file
        settings_file = os.path.join(self.dirn,"auto_process.ini")
        with open(settings_file,'w') as s:
            s.write("""[organism:human]
star_index = /data/hg38/star
bowtie_index = /data/hg38/bowtie

[organism:mouse]
star_index = /data/mm10/star
bowtie_index = /data/mm10/bowtie

[organism:fly]
star_index = /data/dm6/star
""")
        # Load settings
        s = Settings(settings_file)
        # Check reference data can be fetched
        self.assertEqual(fetch_reference_data(s,'star_index'),
                         {
                             'human': '/data/hg38/star',
                             'mouse': '/data/mm10/star',
                             'fly': '/data/dm6/star'
                         })
        self.assertEqual(fetch_reference_data(s,'bowtie_index'),
                         {
                             'human': '/data/hg38/bowtie',
                             'mouse': '/data/mm10/bowtie'
                         })
