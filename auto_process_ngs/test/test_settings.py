#######################################################################
# Tests for settings.py module
#######################################################################

import unittest
import cStringIO
from bcftbx.JobRunner import SimpleJobRunner,GEJobRunner
from auto_process_ngs.settings import *

class TestSettings(unittest.TestCase):
    """Tests for the Settings class
    """
    def test_sample_settings_file(self):
        sample_settings_file = os.path.join(get_config_dir(),
                                            'settings.ini.sample')
        self.assertTrue(os.path.isfile(sample_settings_file),
                        "Missing sample file %s" % sample_settings_file)
        s = Settings(sample_settings_file)
        # General settings
        self.assertTrue(isinstance(s.general.default_runner,SimpleJobRunner))
        self.assertEqual(s.general.max_concurrent_jobs,12)
        self.assertEqual(s.modulefiles.make_fastqs,None)
        self.assertEqual(s.modulefiles.run_qc,None)
        # Bcl2fastq
        self.assertEqual(s.bcl2fastq.nprocessors,1)
        self.assertEqual(s.bcl2fastq.default_version,'<=1.8.4')
        self.assertEqual(s.bcl2fastq.no_lane_splitting,False)
        # NextSeq-specific
        self.assertEqual(s.platform.nextseq.bcl2fastq,'>=2.0')
        self.assertEqual(s.platform.nextseq.no_lane_splitting,True)
        # Fastq_stats
        self.assertEqual(s.fastq_stats.nprocessors,1)
        # Job-specific runners
        self.assertTrue(isinstance(s.runners.bcl2fastq,SimpleJobRunner))
        self.assertTrue(isinstance(s.runners.qc,SimpleJobRunner))
        self.assertTrue(isinstance(s.runners.stats,SimpleJobRunner))
        # Archiving
        self.assertEqual(s.archive.dirn,None)
        self.assertEqual(s.archive.log,None)
        self.assertEqual(s.archive.group,None)
        self.assertEqual(s.archive.chmod,None)
        # QC reporting
        self.assertEqual(s.qc_web_server.dirn,None)
        self.assertEqual(s.qc_web_server.url,None)
