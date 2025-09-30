#######################################################################
# Tests for report_cmd.py module
#######################################################################

import unittest
import tempfile
import shutil
import os
from auto_process_ngs.auto_processor import AutoProcess
from auto_process_ngs.analysis import AnalysisProject
from auto_process_ngs.mock import MockAnalysisDirFactory
from auto_process_ngs.commands.report_cmd import ReportingMode
from auto_process_ngs.commands.report_cmd import report
from auto_process_ngs.commands.report_cmd import report_info
from auto_process_ngs.commands.report_cmd import report_concise
from auto_process_ngs.commands.report_cmd import report_summary
from auto_process_ngs.commands.report_cmd import report_projects
from auto_process_ngs.commands.report_cmd import fetch_value
from auto_process_ngs.commands.report_cmd import default_value
from auto_process_ngs.commands.report_cmd import get_multiplexed_samples

# Unit tests

class TestReportInfo(unittest.TestCase):
    """
    Tests for the 'report' command ('info' mode)
    """
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestReportInfo')
        # Store original location so we can get back at the end
        self.pwd = os.getcwd()
        # Move to working dir
        os.chdir(self.dirn)

    def tearDown(self):
        # Return to original dir
        os.chdir(self.pwd)
        # Remove the temporary test directory
        def del_rw(action,name,excinfo):
            # Explicitly remove read only files/
            # dirs
            os.chmod(os.path.dirname(name),0o755)
            os.chmod(name,0o655)
            os.remove(name)
        shutil.rmtree(self.dirn,onerror=del_rw)

    def test_report_info(self):
        """report: report run in 'info' mode
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": 87, },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "RNA-seq",
                        "Organism": "Human",
                        "PI": "Audrey Bower" },
                "CDE": { "User": "Charles David Edwards",
                         "Library type": "ChIP-seq",
                         "Organism": "Mouse",
                         "PI": "Colin Delaney Eccleston" }
            },
            top_dir=self.dirn)
        mockdir.create()
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Generate full report
        expected = """Run ID       : MISEQ_170901#87
Directory    : %s
Platform     : miseq
Unaligned dir: bcl2fastq

Summary of data in 'bcl2fastq' dir:

- AB: AB1-2 (2 paired end samples)
- CDE: CDE3-4 (2 paired end samples)

3 analysis projects:

- AB
  --
  User    : Alison Bell
  PI      : Audrey Bower
  Library : RNA-seq
  SC Plat.: None
  Organism: Human
  Dir     : AB
  #samples: 2
  #cells  : 
  Samples : AB1-2
  QC      : not verified
  Comments: None

- CDE
  ---
  User    : Charles David Edwards
  PI      : Colin Delaney Eccleston
  Library : ChIP-seq
  SC Plat.: None
  Organism: Mouse
  Dir     : CDE
  #samples: 2
  #cells  : 
  Samples : CDE3-4
  QC      : not verified
  Comments: None

- undetermined
  ------------
  User    : None
  PI      : None
  Library : None
  SC Plat.: None
  Organism: None
  Dir     : undetermined
  #samples: 1
  #cells  : 
  Samples : Undetermined
  QC      : not verified
  Comments: None""" % mockdir.dirn
        for o,e in zip(report_info(ap).split('\n'),
                       expected.split('\n')):
            self.assertEqual(o,e)

    def test_report_info_single_cell(self):
        """report: report single-cell run in 'info' mode
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": 87, },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "scRNA-seq",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Single cell platform": "10xGenomics Chromium 3'",
                        "Number of cells": 1311
                        },
                "CDE": { "User": "Charles David Edwards",
                         "Library type": "ChIP-seq",
                         "Organism": "Mouse",
                         "PI": "Colin Delaney Eccleston" }
            },
            top_dir=self.dirn)
        mockdir.create()
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Generate full report
        expected = """Run ID       : MISEQ_170901#87
Directory    : %s
Platform     : miseq
Unaligned dir: bcl2fastq

Summary of data in 'bcl2fastq' dir:

- AB: AB1-2 (2 paired end samples)
- CDE: CDE3-4 (2 paired end samples)

3 analysis projects:

- AB
  --
  User    : Alison Bell
  PI      : Audrey Bower
  Library : scRNA-seq
  SC Plat.: 10xGenomics Chromium 3'
  Organism: Human
  Dir     : AB
  #samples: 2
  #cells  : 1311
  Samples : AB1-2
  QC      : not verified
  Comments: None

- CDE
  ---
  User    : Charles David Edwards
  PI      : Colin Delaney Eccleston
  Library : ChIP-seq
  SC Plat.: None
  Organism: Mouse
  Dir     : CDE
  #samples: 2
  #cells  : 
  Samples : CDE3-4
  QC      : not verified
  Comments: None

- undetermined
  ------------
  User    : None
  PI      : None
  Library : None
  SC Plat.: None
  Organism: None
  Dir     : undetermined
  #samples: 1
  #cells  : 
  Samples : Undetermined
  QC      : not verified
  Comments: None""" % mockdir.dirn
        for o,e in zip(report_info(ap).split('\n'),
                       expected.split('\n')):
            self.assertEqual(o,e)

    def test_report_info_no_projects(self):
        """report: report run with no projects in 'info' mode
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": 87, },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Generate full report
        expected = """Run ID       : MISEQ_170901#87
Directory    : %s
Platform     : miseq
Unaligned dir: bcl2fastq

Summary of data in 'bcl2fastq' dir:

- AB: AB1-2 (2 paired end samples)
- CDE: CDE3-4 (2 paired end samples)

No analysis projects found""" % mockdir.dirn
        for o,e in zip(report_info(ap).split('\n'),
                       expected.split('\n')):
            self.assertEqual(o,e)

    def test_report_info_10x_cellplex(self):
        """report: report 10xGenomics CellPlex run in 'info' mode
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": 87, },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "CellPlex",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Single cell platform": "10xGenomics Chromium 3'v3",
                        "Number of cells": 1311
                        },
                "CDE": { "User": "Charles David Edwards",
                         "Library type": "ChIP-seq",
                         "Organism": "Mouse",
                         "PI": "Colin Delaney Eccleston" }
            },
            top_dir=self.dirn)
        mockdir.create()
        # Add a cellranger multi config.csv file
        with open(os.path.join(mockdir.dirn,
                               "AB",
                               "10x_multi_config.csv"),'wt') as fp:
            fastq_dir = os.path.join(mockdir.dirn,
                                     "AB",
                                     "fastqs")
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
AB1,%s,any,AB1,gene expression,
AB2,%s,any,AB2,Multiplexing Capture,

[samples]
sample_id,cmo_ids,description
ABM1,CMO301,ABM1
ABM2,CMO302,ABM2
ABM3,CMO303,ABM3
ABM4,CMO304,ABM4
""" % (fastq_dir,fastq_dir))
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Generate full report
        expected = """Run ID       : MISEQ_170901#87
Directory    : %s
Platform     : miseq
Unaligned dir: bcl2fastq

Summary of data in 'bcl2fastq' dir:

- AB: AB1-2 (2 paired end samples)
- CDE: CDE3-4 (2 paired end samples)

3 analysis projects:

- AB
  --
  User    : Alison Bell
  PI      : Audrey Bower
  Library : CellPlex
  SC Plat.: 10xGenomics Chromium 3'v3
  Organism: Human
  Dir     : AB
  #samples: 4 multiplexed (2 physical)
  #cells  : 1311
  Samples : ABM1-4 (AB1-2)
  QC      : not verified
  Comments: None

- CDE
  ---
  User    : Charles David Edwards
  PI      : Colin Delaney Eccleston
  Library : ChIP-seq
  SC Plat.: None
  Organism: Mouse
  Dir     : CDE
  #samples: 2
  #cells  : 
  Samples : CDE3-4
  QC      : not verified
  Comments: None

- undetermined
  ------------
  User    : None
  PI      : None
  Library : None
  SC Plat.: None
  Organism: None
  Dir     : undetermined
  #samples: 1
  #cells  : 
  Samples : Undetermined
  QC      : not verified
  Comments: None""" % mockdir.dirn
        for o,e in zip(report_info(ap).split('\n'),
                       expected.split('\n')):
            self.assertEqual(o,e)

    def test_report_info_10x_cellplex_scrnaseq(self):
        """report: report 10xGenomics CellPlex scRNA-seq run in 'info' mode
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": 87, },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "CellPlex scRNA-seq",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Single cell platform": "10xGenomics Chromium 3'v3",
                        "Number of cells": 1311
                        },
                "CDE": { "User": "Charles David Edwards",
                         "Library type": "ChIP-seq",
                         "Organism": "Mouse",
                         "PI": "Colin Delaney Eccleston" }
            },
            top_dir=self.dirn)
        mockdir.create()
        # Add a cellranger multi config.csv file
        with open(os.path.join(mockdir.dirn,
                               "AB",
                               "10x_multi_config.csv"),'wt') as fp:
            fastq_dir = os.path.join(mockdir.dirn,
                                     "AB",
                                     "fastqs")
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
AB1,%s,any,AB1,gene expression,
AB2,%s,any,AB2,Multiplexing Capture,

[samples]
sample_id,cmo_ids,description
ABM1,CMO301,ABM1
ABM2,CMO302,ABM2
ABM3,CMO303,ABM3
ABM4,CMO304,ABM4
""" % (fastq_dir,fastq_dir))
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Generate full report
        expected = """Run ID       : MISEQ_170901#87
Directory    : %s
Platform     : miseq
Unaligned dir: bcl2fastq

Summary of data in 'bcl2fastq' dir:

- AB: AB1-2 (2 paired end samples)
- CDE: CDE3-4 (2 paired end samples)

3 analysis projects:

- AB
  --
  User    : Alison Bell
  PI      : Audrey Bower
  Library : CellPlex scRNA-seq
  SC Plat.: 10xGenomics Chromium 3'v3
  Organism: Human
  Dir     : AB
  #samples: 4 multiplexed (2 physical)
  #cells  : 1311
  Samples : ABM1-4 (AB1-2)
  QC      : not verified
  Comments: None

- CDE
  ---
  User    : Charles David Edwards
  PI      : Colin Delaney Eccleston
  Library : ChIP-seq
  SC Plat.: None
  Organism: Mouse
  Dir     : CDE
  #samples: 2
  #cells  : 
  Samples : CDE3-4
  QC      : not verified
  Comments: None

- undetermined
  ------------
  User    : None
  PI      : None
  Library : None
  SC Plat.: None
  Organism: None
  Dir     : undetermined
  #samples: 1
  #cells  : 
  Samples : Undetermined
  QC      : not verified
  Comments: None""" % mockdir.dirn
        for o,e in zip(report_info(ap).split('\n'),
                       expected.split('\n')):
            self.assertEqual(o,e)

    def test_report_info_10x_flex(self):
        """report: report 10xGenomics Flex run in 'info' mode
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '240213_A00879_0087_000000000-AGEW9',
            'novaseq',
            metadata={ "source": "testing",
                       "run_number": 87, },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "Flex",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Single cell platform": "10xGenomics Chromium 3'v3",
                        "Number of cells": 1311
                        },
                "CDE": { "User": "Charles David Edwards",
                         "Library type": "ChIP-seq",
                         "Organism": "Mouse",
                         "PI": "Colin Delaney Eccleston" }
            },
            top_dir=self.dirn)
        mockdir.create()
        # Add a cellranger multi config.csv file
        with open(os.path.join(mockdir.dirn,
                               "AB",
                               "10x_multi_config.csv"),'wt') as fp:
            fastq_dir = os.path.join(mockdir.dirn,
                                     "AB",
                                     "fastqs")
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A
probe-set,/data/Probe_Set_v1.0_GRCh38-2020-A.csv
no-bam,true

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB1_Flex,{fastq_dir},any,PJB1,Gene Expression,

[samples]
sample_id,probe_barcode_ids,description
ABF1,BC001,ABF1
ABF2,BC002,ABF2
ABF3,BC003,ABF3
ABF4,BC004,ABF4
""".format(fastq_dir=fastq_dir))
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Generate full report
        expected = """Run ID       : NOVASEQ_240213#87
Directory    : %s
Platform     : novaseq
Unaligned dir: bcl2fastq

Summary of data in 'bcl2fastq' dir:

- AB: AB1-2 (2 paired end samples)
- CDE: CDE3-4 (2 paired end samples)

3 analysis projects:

- AB
  --
  User    : Alison Bell
  PI      : Audrey Bower
  Library : Flex
  SC Plat.: 10xGenomics Chromium 3'v3
  Organism: Human
  Dir     : AB
  #samples: 4 multiplexed (2 physical)
  #cells  : 1311
  Samples : ABF1-4 (AB1-2)
  QC      : not verified
  Comments: None

- CDE
  ---
  User    : Charles David Edwards
  PI      : Colin Delaney Eccleston
  Library : ChIP-seq
  SC Plat.: None
  Organism: Mouse
  Dir     : CDE
  #samples: 2
  #cells  : 
  Samples : CDE3-4
  QC      : not verified
  Comments: None

- undetermined
  ------------
  User    : None
  PI      : None
  Library : None
  SC Plat.: None
  Organism: None
  Dir     : undetermined
  #samples: 1
  #cells  : 
  Samples : Undetermined
  QC      : not verified
  Comments: None""" % mockdir.dirn
        for o,e in zip(report_info(ap).split('\n'),
                       expected.split('\n')):
            self.assertEqual(o,e)

    def test_report_info_explicit_multiplexed_samples(self):
        """report: report run with multiplexed samples in 'info' mode
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": 87, },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "RNA-seq",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Multiplexed samples": "M1,M2,M3,M4" },
                "CDE": { "User": "Charles David Edwards",
                         "Library type": "ChIP-seq",
                         "Organism": "Mouse",
                         "PI": "Colin Delaney Eccleston" }
            },
            top_dir=self.dirn)
        mockdir.create()
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Generate full report
        expected = """Run ID       : MISEQ_170901#87
Directory    : %s
Platform     : miseq
Unaligned dir: bcl2fastq

Summary of data in 'bcl2fastq' dir:

- AB: AB1-2 (2 paired end samples)
- CDE: CDE3-4 (2 paired end samples)

3 analysis projects:

- AB
  --
  User    : Alison Bell
  PI      : Audrey Bower
  Library : RNA-seq
  SC Plat.: None
  Organism: Human
  Dir     : AB
  #samples: 4 multiplexed (2 physical)
  #cells  : 
  Samples : M1-4 (AB1-2)
  QC      : not verified
  Comments: None

- CDE
  ---
  User    : Charles David Edwards
  PI      : Colin Delaney Eccleston
  Library : ChIP-seq
  SC Plat.: None
  Organism: Mouse
  Dir     : CDE
  #samples: 2
  #cells  : 
  Samples : CDE3-4
  QC      : not verified
  Comments: None

- undetermined
  ------------
  User    : None
  PI      : None
  Library : None
  SC Plat.: None
  Organism: None
  Dir     : undetermined
  #samples: 1
  #cells  : 
  Samples : Undetermined
  QC      : not verified
  Comments: None""" % mockdir.dirn
        for o,e in zip(report_info(ap).split('\n'),
                       expected.split('\n')):
            self.assertEqual(o,e)

    def test_report_info_no_projects(self):
        """report: report run with no projects in 'info' mode
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": 87, },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Generate full report
        expected = """Run ID       : MISEQ_170901#87
Directory    : %s
Platform     : miseq
Unaligned dir: bcl2fastq

Summary of data in 'bcl2fastq' dir:

- AB: AB1-2 (2 paired end samples)
- CDE: CDE3-4 (2 paired end samples)

No analysis projects found""" % mockdir.dirn
        for o,e in zip(report_info(ap).split('\n'),
                       expected.split('\n')):
            self.assertEqual(o,e)

    def test_report_info_explicit_unknown_multiplexed_samples(self):
        """report: report run with "unknown" multiplexed samples in 'info' mode
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": 87, },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "RNA-seq",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Multiplexed samples": "?" },
                "CDE": { "User": "Charles David Edwards",
                         "Library type": "ChIP-seq",
                         "Organism": "Mouse",
                         "PI": "Colin Delaney Eccleston" }
            },
            top_dir=self.dirn)
        mockdir.create()
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Generate full report
        expected = """Run ID       : MISEQ_170901#87
Directory    : %s
Platform     : miseq
Unaligned dir: bcl2fastq

Summary of data in 'bcl2fastq' dir:

- AB: AB1-2 (2 paired end samples)
- CDE: CDE3-4 (2 paired end samples)

3 analysis projects:

- AB
  --
  User    : Alison Bell
  PI      : Audrey Bower
  Library : RNA-seq
  SC Plat.: None
  Organism: Human
  Dir     : AB
  #samples: ? multiplexed (2 physical)
  #cells  : 
  Samples : ? (AB1-2)
  QC      : not verified
  Comments: None

- CDE
  ---
  User    : Charles David Edwards
  PI      : Colin Delaney Eccleston
  Library : ChIP-seq
  SC Plat.: None
  Organism: Mouse
  Dir     : CDE
  #samples: 2
  #cells  : 
  Samples : CDE3-4
  QC      : not verified
  Comments: None

- undetermined
  ------------
  User    : None
  PI      : None
  Library : None
  SC Plat.: None
  Organism: None
  Dir     : undetermined
  #samples: 1
  #cells  : 
  Samples : Undetermined
  QC      : not verified
  Comments: None""" % mockdir.dirn
        for o,e in zip(report_info(ap).split('\n'),
                       expected.split('\n')):
            self.assertEqual(o,e)

    def test_report_info_no_projects(self):
        """report: report run with no projects in 'info' mode
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": 87, },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Generate full report
        expected = """Run ID       : MISEQ_170901#87
Directory    : %s
Platform     : miseq
Unaligned dir: bcl2fastq

Summary of data in 'bcl2fastq' dir:

- AB: AB1-2 (2 paired end samples)
- CDE: CDE3-4 (2 paired end samples)

No analysis projects found""" % mockdir.dirn
        for o,e in zip(report_info(ap).split('\n'),
                       expected.split('\n')):
            self.assertEqual(o,e)

    def test_report_info_legacy_metadata_items(self):
        """report: report run in 'info' mode with legacy metadata items
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": 87, },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "RNA-seq",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "ICELL8 well list": None },
                "CDE": { "User": "Charles David Edwards",
                         "Library type": "ChIP-seq",
                         "Organism": "Mouse",
                         "PI": "Colin Delaney Eccleston",
                         "ICELL8 well list": None }
            },
            top_dir=self.dirn)
        mockdir.create()
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Generate full report
        expected = """Run ID       : MISEQ_170901#87
Directory    : %s
Platform     : miseq
Unaligned dir: bcl2fastq

Summary of data in 'bcl2fastq' dir:

- AB: AB1-2 (2 paired end samples)
- CDE: CDE3-4 (2 paired end samples)

3 analysis projects:

- AB
  --
  User    : Alison Bell
  PI      : Audrey Bower
  Library : RNA-seq
  SC Plat.: None
  Organism: Human
  Dir     : AB
  #samples: 2
  #cells  : 
  Samples : AB1-2
  QC      : not verified
  Comments: None

- CDE
  ---
  User    : Charles David Edwards
  PI      : Colin Delaney Eccleston
  Library : ChIP-seq
  SC Plat.: None
  Organism: Mouse
  Dir     : CDE
  #samples: 2
  #cells  : 
  Samples : CDE3-4
  QC      : not verified
  Comments: None

- undetermined
  ------------
  User    : None
  PI      : None
  Library : None
  SC Plat.: None
  Organism: None
  Dir     : undetermined
  #samples: 1
  #cells  : 
  Samples : Undetermined
  QC      : not verified
  Comments: None""" % mockdir.dirn
        for o,e in zip(report_info(ap).split('\n'),
                       expected.split('\n')):
            self.assertEqual(o,e)

class TestReportConcise(unittest.TestCase):
    """
    Tests for the 'report' command ('concise' mode)
    """
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestReportConcise')
        # Store original location so we can get back at the end
        self.pwd = os.getcwd()
        # Move to working dir
        os.chdir(self.dirn)

    def tearDown(self):
        # Return to original dir
        os.chdir(self.pwd)
        # Remove the temporary test directory
        def del_rw(action,name,excinfo):
            # Explicitly remove read only files/
            # dirs
            os.chmod(os.path.dirname(name),0o755)
            os.chmod(name,0o655)
            os.remove(name)
        shutil.rmtree(self.dirn,onerror=del_rw)

    def test_report_concise(self):
        """report: report run in 'concise' mode
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": 87,
                       "sequencer_model": "MiSeq" },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "RNA-seq",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Sequencer model": "MiSeq" },
                "CDE": { "User": "Charles David Edwards",
                         "Library type": "ChIP-seq",
                         "Organism": "Mouse",
                         "PI": "Colin Delaney Eccleston",
                         "Sequencer model": "MiSeq" }
            },
            top_dir=self.dirn)
        mockdir.create()
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Generate concise report
        self.assertEqual(report_concise(ap),
                         "Paired end: 'AB': Alison Bell, Human RNA-seq (PI: Audrey Bower) (2 samples); 'CDE': Charles David Edwards, Mouse ChIP-seq (PI: Colin Delaney Eccleston) (2 samples)")

    def test_report_concise_single_cell(self):
        """report: report single-cell run in 'concise' mode
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": 87, },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "scRNA-seq",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Single cell platform": "10xGenomics Chromium 3'",
                        "Number of cells": 1311 },
                "CDE": { "User": "Charles David Edwards",
                         "Library type": "ChIP-seq",
                         "Organism": "Mouse",
                         "PI": "Colin Delaney Eccleston" }
            },
            top_dir=self.dirn)
        mockdir.create()
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Generate concise report
        self.assertEqual(report_concise(ap),
                         "Paired end: 'AB': Alison Bell, Human 10xGenomics Chromium 3' scRNA-seq (PI: Audrey Bower) (2 samples/1311 cells); 'CDE': Charles David Edwards, Mouse ChIP-seq (PI: Colin Delaney Eccleston) (2 samples)")

    def test_report_concise_10x_cellplex(self):
        """report: report 10xGenomics CellPlex run in 'concise' mode
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": 87, },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "CellPlex",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Single cell platform": "10xGenomics Chromium 3'v3",
                        "Number of cells": 1311 },
                "CDE": { "User": "Charles David Edwards",
                         "Library type": "ChIP-seq",
                         "Organism": "Mouse",
                         "PI": "Colin Delaney Eccleston" }
            },
            top_dir=self.dirn)
        mockdir.create()
        # Add a cellranger multi config.csv file
        with open(os.path.join(mockdir.dirn,
                               "AB",
                               "10x_multi_config.csv"),'wt') as fp:
            fastq_dir = os.path.join(mockdir.dirn,
                                     "AB",
                                     "fastqs")
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
AB1,%s,any,AB1,gene expression,
AB2,%s,any,AB2,Multiplexing Capture,

[samples]
sample_id,cmo_ids,description
ABM1,CMO301,ABM1
ABM2,CMO302,ABM2
ABM3,CMO303,ABM3
ABM4,CMO304,ABM4
""" % (fastq_dir,fastq_dir))
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Generate concise report
        self.assertEqual(report_concise(ap),
                         "Paired end: 'AB': Alison Bell, Human 10xGenomics Chromium 3'v3 CellPlex (PI: Audrey Bower) (4 multiplexed samples/1311 cells); 'CDE': Charles David Edwards, Mouse ChIP-seq (PI: Colin Delaney Eccleston) (2 samples)")

    def test_report_concise_10x_cellplex_scrnaseq(self):
        """report: report 10xGenomics CellPlex scRNA-seq run in 'concise' mode
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": 87, },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "CellPlex scRNA-seq",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Single cell platform": "10xGenomics Chromium 3'v3",
                        "Number of cells": 1311 },
                "CDE": { "User": "Charles David Edwards",
                         "Library type": "ChIP-seq",
                         "Organism": "Mouse",
                         "PI": "Colin Delaney Eccleston" }
            },
            top_dir=self.dirn)
        mockdir.create()
        # Add a cellranger multi config.csv file
        with open(os.path.join(mockdir.dirn,
                               "AB",
                               "10x_multi_config.csv"),'wt') as fp:
            fastq_dir = os.path.join(mockdir.dirn,
                                     "AB",
                                     "fastqs")
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
AB1,%s,any,AB1,gene expression,
AB2,%s,any,AB2,Multiplexing Capture,

[samples]
sample_id,cmo_ids,description
ABM1,CMO301,ABM1
ABM2,CMO302,ABM2
ABM3,CMO303,ABM3
ABM4,CMO304,ABM4
""" % (fastq_dir,fastq_dir))
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Generate concise report
        self.assertEqual(report_concise(ap),
                         "Paired end: 'AB': Alison Bell, Human 10xGenomics Chromium 3'v3 CellPlex scRNA-seq (PI: Audrey Bower) (4 multiplexed samples/1311 cells); 'CDE': Charles David Edwards, Mouse ChIP-seq (PI: Colin Delaney Eccleston) (2 samples)")

    def test_report_concise_10x_flex(self):
        """report: report 10xGenomics Flex run in 'concise' mode
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '240213_A00879_0087_000000000-AGEW9',
            'novaseq',
            metadata={ "source": "testing",
                       "run_number": 87, },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "Flex",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Single cell platform": "10xGenomics Chromium 3'v3",
                        "Number of cells": 1311 },
                "CDE": { "User": "Charles David Edwards",
                         "Library type": "ChIP-seq",
                         "Organism": "Mouse",
                         "PI": "Colin Delaney Eccleston" }
            },
            top_dir=self.dirn)
        mockdir.create()
        # Add a cellranger multi config.csv file
        with open(os.path.join(mockdir.dirn,
                               "AB",
                               "10x_multi_config.csv"),'wt') as fp:
            fastq_dir = os.path.join(mockdir.dirn,
                                     "AB",
                                     "fastqs")
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A
probe-set,/data/Probe_Set_v1.0_GRCh38-2020-A.csv
no-bam,true

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB1_Flex,{fastq_dir},any,PJB1,Gene Expression,

[samples]
sample_id,probe_barcode_ids,description
ABF1,BC001,ABF1
ABF2,BC002,ABF2
ABF3,BC003,ABF3
ABF4,BC004,ABF4
""".format(fastq_dir=fastq_dir))
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Generate concise report
        self.assertEqual(report_concise(ap),
                         "Paired end: 'AB': Alison Bell, Human 10xGenomics Chromium 3'v3 Flex (PI: Audrey Bower) (4 multiplexed samples/1311 cells); 'CDE': Charles David Edwards, Mouse ChIP-seq (PI: Colin Delaney Eccleston) (2 samples)")

    def test_report_concise_multiplexed_samples(self):
        """report: report run with multiplexed samples in 'concise' mode
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": 87,
                       "sequencer_model": "MiSeq" },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "RNA-seq",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Sequencer model": "MiSeq",
                        "Multiplexed samples": "M1,M2,M3,M4" },
                "CDE": { "User": "Charles David Edwards",
                         "Library type": "ChIP-seq",
                         "Organism": "Mouse",
                         "PI": "Colin Delaney Eccleston",
                         "Sequencer model": "MiSeq" }
            },
            top_dir=self.dirn)
        mockdir.create()
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Generate concise report
        self.assertEqual(report_concise(ap),
                         "Paired end: 'AB': Alison Bell, Human RNA-seq (PI: Audrey Bower) (4 multiplexed samples); 'CDE': Charles David Edwards, Mouse ChIP-seq (PI: Colin Delaney Eccleston) (2 samples)")

    def test_report_concise_unknown_multiplexed_samples(self):
        """report: report run with "unknown" multiplexed samples in 'concise' mode
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": 87,
                       "sequencer_model": "MiSeq" },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "RNA-seq",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Sequencer model": "MiSeq",
                        "Multiplexed samples": "?" },
                "CDE": { "User": "Charles David Edwards",
                         "Library type": "ChIP-seq",
                         "Organism": "Mouse",
                         "PI": "Colin Delaney Eccleston",
                         "Sequencer model": "MiSeq" }
            },
            top_dir=self.dirn)
        mockdir.create()
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Generate concise report
        self.assertEqual(report_concise(ap),
                         "Paired end: 'AB': Alison Bell, Human RNA-seq (PI: Audrey Bower) (? multiplexed samples); 'CDE': Charles David Edwards, Mouse ChIP-seq (PI: Colin Delaney Eccleston) (2 samples)")

    def test_report_concise_legacy_metadata_items(self):
        """report: report run in 'concise' mode with legacy metadata items
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": 87,
                       "sequencer_model": "MiSeq" },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "RNA-seq",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Sequencer model": "MiSeq",
                        "ICELL8 well list": None },
                "CDE": { "User": "Charles David Edwards",
                         "Library type": "ChIP-seq",
                         "Organism": "Mouse",
                         "PI": "Colin Delaney Eccleston",
                         "Sequencer model": "MiSeq",
                         "ICELL8 well list": None}
            },
            top_dir=self.dirn)
        mockdir.create()
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Generate concise report
        self.assertEqual(report_concise(ap),
                         "Paired end: 'AB': Alison Bell, Human RNA-seq (PI: Audrey Bower) (2 samples); 'CDE': Charles David Edwards, Mouse ChIP-seq (PI: Colin Delaney Eccleston) (2 samples)")

    def test_report_concise_no_projects(self):
        """report: report run with no projects in 'concise' mode
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": 87, },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Generate concise report
        self.assertEqual(report_concise(ap),
                         "Paired end: no projects found; contents of 'bcl2fastq' are: 'AB' (2 samples), 'CDE' (2 samples)")

class TestReportSummary(unittest.TestCase):
    """
    Tests for the 'report' command ('summary' mode)
    """
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestReportSummary')
        # Store original location so we can get back at the end
        self.pwd = os.getcwd()
        # Move to working dir
        os.chdir(self.dirn)

    def tearDown(self):
        # Return to original dir
        os.chdir(self.pwd)
        # Remove the temporary test directory
        def del_rw(action,name,excinfo):
            # Explicitly remove read only files/
            # dirs
            os.chmod(os.path.dirname(name),0o755)
            os.chmod(name,0o655)
            os.remove(name)
        shutil.rmtree(self.dirn,onerror=del_rw)

    def test_report_summary(self):
        """report: report run in 'summary' mode
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": 87,
                       "sequencer_model": "MiSeq" },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "RNA-seq",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Sequencer model": "MiSeq" },
                "CDE": { "User": "Charles David Edwards",
                         "Library type": "ChIP-seq",
                         "Organism": "Mouse",
                         "PI": "Colin Delaney Eccleston",
                         "Sequencer model": "MiSeq",
                         "Comments": "Repeat of previous run" }
            },
            top_dir=self.dirn)
        mockdir.create()
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Generate summary report
        expected = """MISEQ run #87 datestamped 170901
================================
Run name : 170901_M00879_0087_000000000-AGEW9
Reference: MISEQ_170901#87
Platform : MISEQ
Sequencer: MiSeq
Directory: %s
Endedness: Paired end
Bcl2fastq: Unknown

2 projects:
- 'AB':  Alison Bell           Human RNA-seq  2 samples (PI Audrey Bower)           
- 'CDE': Charles David Edwards Mouse ChIP-seq 2 samples (PI Colin Delaney Eccleston)

Additional notes/comments:
- CDE: Repeat of previous run
""" % mockdir.dirn
        for o,e in zip(report_summary(ap).split('\n'),
                       expected.split('\n')):
            self.assertEqual(o,e)

    def test_report_summary_with_analysis_number(self):
        """report: report run with analysis number in 'summary' mode
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": 87,
                       "sequencer_model": "MiSeq",
                       "analysis_number": 2 },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "RNA-seq",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Sequencer model": "MiSeq" },
                "CDE": { "User": "Charles David Edwards",
                         "Library type": "ChIP-seq",
                         "Organism": "Mouse",
                         "PI": "Colin Delaney Eccleston",
                         "Sequencer model": "MiSeq",
                         "Comments": "Repeat of previous run" }
            },
            top_dir=self.dirn)
        mockdir.create()
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Generate summary report
        expected = """MISEQ run #87 datestamped 170901 [analysis #2]
==============================================
Run name : 170901_M00879_0087_000000000-AGEW9
Reference: MISEQ_170901#87.2
Platform : MISEQ
Sequencer: MiSeq
Directory: %s
Endedness: Paired end
Bcl2fastq: Unknown

2 projects:
- 'AB':  Alison Bell           Human RNA-seq  2 samples (PI Audrey Bower)           
- 'CDE': Charles David Edwards Mouse ChIP-seq 2 samples (PI Colin Delaney Eccleston)

Additional notes/comments:
- CDE: Repeat of previous run
""" % mockdir.dirn
        for o,e in zip(report_summary(ap).split('\n'),
                       expected.split('\n')):
            self.assertEqual(o,e)

    def test_report_summary_single_cell(self):
        """report: report single-cell run in 'summary' mode
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": 87,
                       "bcl2fastq_software":
                       "('/usr/bin/bcl2fastq', 'bcl2fastq', '2.17.1.14')",
                       "cellranger_software":
                       "('/usr/bin/cellranger', 'cellranger', '3.0.1')",
                       "sequencer_model": "MiSeq" },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "scRNA-seq",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Single cell platform": "10xGenomics Chromium 3'",
                        "Number of cells": 1311 },
                "CDE": { "User": "Charles David Edwards",
                         "Library type": "ChIP-seq",
                         "Organism": "Mouse",
                         "PI": "Colin Delaney Eccleston",
                         "Comments": "Repeat of previous run" }
            },
            top_dir=self.dirn)
        mockdir.create()
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Generate summary report
        expected = """MISEQ run #87 datestamped 170901
================================
Run name  : 170901_M00879_0087_000000000-AGEW9
Reference : MISEQ_170901#87
Platform  : MISEQ
Sequencer : MiSeq
Directory : %s
Endedness : Paired end
Bcl2fastq : bcl2fastq 2.17.1.14
Cellranger: cellranger 3.0.1

2 projects:
- 'AB':  Alison Bell           Human scRNA-seq (10xGenomics Chromium 3') 2 samples/1311 cells (PI Audrey Bower)           
- 'CDE': Charles David Edwards Mouse ChIP-seq                            2 samples            (PI Colin Delaney Eccleston)

Additional notes/comments:
- CDE: Repeat of previous run
""" % mockdir.dirn
        for o,e in zip(report_summary(ap).split('\n'),
                       expected.split('\n')):
            self.assertEqual(o,e)

    def test_report_summary_10x_cellplex(self):
        """report: report 10xGenomics CellPlex run in 'summary' mode
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": 87,
                       "bcl2fastq_software":
                       "('/usr/bin/bcl2fastq', 'bcl2fastq', '2.17.1.14')",
                       "cellranger_software":
                       "('/usr/bin/cellranger', 'cellranger', '3.0.1')",
                       "sequencer_model": "MiSeq" },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "CellPlex",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Single cell platform": "10xGenomics Chromium 3'v3",
                        "Number of cells": 1311 },
                "CDE": { "User": "Charles David Edwards",
                         "Library type": "ChIP-seq",
                         "Organism": "Mouse",
                         "PI": "Colin Delaney Eccleston",
                         "Comments": "Repeat of previous run" }
            },
            top_dir=self.dirn)
        mockdir.create()
        # Add a cellranger multi config.csv file
        with open(os.path.join(mockdir.dirn,
                               "AB",
                               "10x_multi_config.csv"),'wt') as fp:
            fastq_dir = os.path.join(mockdir.dirn,
                                     "AB",
                                     "fastqs")
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
AB1,%s,any,AB1,gene expression,
AB2,%s,any,AB2,Multiplexing Capture,

[samples]
sample_id,cmo_ids,description
ABM1,CMO301,ABM1
ABM2,CMO302,ABM2
ABM3,CMO303,ABM3
ABM4,CMO304,ABM4
""" % (fastq_dir,fastq_dir))
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Generate summary report
        expected = """MISEQ run #87 datestamped 170901
================================
Run name  : 170901_M00879_0087_000000000-AGEW9
Reference : MISEQ_170901#87
Platform  : MISEQ
Sequencer : MiSeq
Directory : %s
Endedness : Paired end
Bcl2fastq : bcl2fastq 2.17.1.14
Cellranger: cellranger 3.0.1

2 projects:
- 'AB':  Alison Bell           Human CellPlex (10xGenomics Chromium 3'v3) 4 multiplexed samples/1311 cells (PI Audrey Bower)           
- 'CDE': Charles David Edwards Mouse ChIP-seq                             2 samples                        (PI Colin Delaney Eccleston)

Additional notes/comments:
- CDE: Repeat of previous run
""" % mockdir.dirn
        for o,e in zip(report_summary(ap).split('\n'),
                       expected.split('\n')):
            self.assertEqual(o,e)

    def test_report_summary_10x_cellplex_scrnaseq(self):
        """report: report 10xGenomics CellPlex scRNA-seq run in 'summary' mode
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": 87,
                       "bcl2fastq_software":
                       "('/usr/bin/bcl2fastq', 'bcl2fastq', '2.17.1.14')",
                       "cellranger_software":
                       "('/usr/bin/cellranger', 'cellranger', '3.0.1')",
                       "sequencer_model": "MiSeq" },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "CellPlex scRNA-seq",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Single cell platform": "10xGenomics Chromium 3'v3",
                        "Number of cells": 1311 },
                "CDE": { "User": "Charles David Edwards",
                         "Library type": "ChIP-seq",
                         "Organism": "Mouse",
                         "PI": "Colin Delaney Eccleston",
                         "Comments": "Repeat of previous run" }
            },
            top_dir=self.dirn)
        mockdir.create()
        # Add a cellranger multi config.csv file
        with open(os.path.join(mockdir.dirn,
                               "AB",
                               "10x_multi_config.csv"),'wt') as fp:
            fastq_dir = os.path.join(mockdir.dirn,
                                     "AB",
                                     "fastqs")
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
AB1,%s,any,AB1,gene expression,
AB2,%s,any,AB2,Multiplexing Capture,

[samples]
sample_id,cmo_ids,description
ABM1,CMO301,ABM1
ABM2,CMO302,ABM2
ABM3,CMO303,ABM3
ABM4,CMO304,ABM4
""" % (fastq_dir,fastq_dir))
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Generate summary report
        expected = """MISEQ run #87 datestamped 170901
================================
Run name  : 170901_M00879_0087_000000000-AGEW9
Reference : MISEQ_170901#87
Platform  : MISEQ
Sequencer : MiSeq
Directory : %s
Endedness : Paired end
Bcl2fastq : bcl2fastq 2.17.1.14
Cellranger: cellranger 3.0.1

2 projects:
- 'AB':  Alison Bell           Human CellPlex scRNA-seq (10xGenomics Chromium 3'v3) 4 multiplexed samples/1311 cells (PI Audrey Bower)           
- 'CDE': Charles David Edwards Mouse ChIP-seq                                       2 samples                        (PI Colin Delaney Eccleston)

Additional notes/comments:
- CDE: Repeat of previous run
""" % mockdir.dirn
        for o,e in zip(report_summary(ap).split('\n'),
                       expected.split('\n')):
            self.assertEqual(o,e)

    def test_report_summary_10x_flex(self):
        """report: report 10xGenomics Flex run in 'summary' mode
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '240213_A00879_0087_000000000-AGEW9',
            'novaseq',
            metadata={ "source": "testing",
                       "run_number": 87,
                       "bcl2fastq_software":
                       "('/usr/bin/bcl2fastq', 'bcl2fastq', '2.17.1.14')",
                       "cellranger_software":
                       "('/usr/bin/cellranger', 'cellranger', '3.0.1')",
                       "sequencer_model": "NovaSeq" },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "Flex",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Single cell platform": "10xGenomics Chromium 3'v3",
                        "Number of cells": 1311 },
                "CDE": { "User": "Charles David Edwards",
                         "Library type": "ChIP-seq",
                         "Organism": "Mouse",
                         "PI": "Colin Delaney Eccleston",
                         "Comments": "Repeat of previous run" }
            },
            top_dir=self.dirn)
        mockdir.create()
        # Add a cellranger multi config.csv file
        with open(os.path.join(mockdir.dirn,
                               "AB",
                               "10x_multi_config.csv"),'wt') as fp:
            fastq_dir = os.path.join(mockdir.dirn,
                                     "AB",
                                     "fastqs")
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A
probe-set,/data/Probe_Set_v1.0_GRCh38-2020-A.csv
no-bam,true

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB1_Flex,{fastq_dir},any,PJB1,Gene Expression,

[samples]
sample_id,probe_barcode_ids,description
ABF1,BC001,ABF1
ABF2,BC002,ABF2
ABF3,BC003,ABF3
ABF4,BC004,ABF4
""".format(fastq_dir=fastq_dir))
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Generate summary report
        expected = """NOVASEQ run #87 datestamped 240213
==================================
Run name  : 240213_A00879_0087_000000000-AGEW9
Reference : NOVASEQ_240213#87
Platform  : NOVASEQ
Sequencer : NovaSeq
Directory : %s
Endedness : Paired end
Bcl2fastq : bcl2fastq 2.17.1.14
Cellranger: cellranger 3.0.1

2 projects:
- 'AB':  Alison Bell           Human Flex (10xGenomics Chromium 3'v3) 4 multiplexed samples/1311 cells (PI Audrey Bower)           
- 'CDE': Charles David Edwards Mouse ChIP-seq                         2 samples                        (PI Colin Delaney Eccleston)

Additional notes/comments:
- CDE: Repeat of previous run
""" % mockdir.dirn
        for o,e in zip(report_summary(ap).split('\n'),
                       expected.split('\n')):
            self.assertEqual(o,e)

    def test_report_summary_multiplexed_samples(self):
        """report: report run with multiplexed samples in 'summary' mode
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": 87,
                       "sequencer_model": "MiSeq" },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "RNA-seq",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Sequencer model": "MiSeq",
                        "Multiplexed samples": "M1,M2,M3,M4" },
                "CDE": { "User": "Charles David Edwards",
                         "Library type": "ChIP-seq",
                         "Organism": "Mouse",
                         "PI": "Colin Delaney Eccleston",
                         "Sequencer model": "MiSeq",
                         "Comments": "Repeat of previous run" }
            },
            top_dir=self.dirn)
        mockdir.create()
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Generate summary report
        expected = """MISEQ run #87 datestamped 170901
================================
Run name : 170901_M00879_0087_000000000-AGEW9
Reference: MISEQ_170901#87
Platform : MISEQ
Sequencer: MiSeq
Directory: %s
Endedness: Paired end
Bcl2fastq: Unknown

2 projects:
- 'AB':  Alison Bell           Human RNA-seq  4 multiplexed samples (PI Audrey Bower)           
- 'CDE': Charles David Edwards Mouse ChIP-seq 2 samples             (PI Colin Delaney Eccleston)

Additional notes/comments:
- CDE: Repeat of previous run
""" % mockdir.dirn
        for o,e in zip(report_summary(ap).split('\n'),
                       expected.split('\n')):
            self.assertEqual(o,e)

    def test_report_summary_unknown_multiplexed_samples(self):
        """report: report run with "unknown" multiplexed samples in 'summary' mode
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": 87,
                       "sequencer_model": "MiSeq" },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "RNA-seq",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Sequencer model": "MiSeq",
                        "Multiplexed samples": "?" },
                "CDE": { "User": "Charles David Edwards",
                         "Library type": "ChIP-seq",
                         "Organism": "Mouse",
                         "PI": "Colin Delaney Eccleston",
                         "Sequencer model": "MiSeq",
                         "Comments": "Repeat of previous run" }
            },
            top_dir=self.dirn)
        mockdir.create()
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Generate summary report
        expected = """MISEQ run #87 datestamped 170901
================================
Run name : 170901_M00879_0087_000000000-AGEW9
Reference: MISEQ_170901#87
Platform : MISEQ
Sequencer: MiSeq
Directory: %s
Endedness: Paired end
Bcl2fastq: Unknown

2 projects:
- 'AB':  Alison Bell           Human RNA-seq  ? multiplexed samples (PI Audrey Bower)           
- 'CDE': Charles David Edwards Mouse ChIP-seq 2 samples             (PI Colin Delaney Eccleston)

Additional notes/comments:
- CDE: Repeat of previous run
""" % mockdir.dirn
        for o,e in zip(report_summary(ap).split('\n'),
                       expected.split('\n')):
            self.assertEqual(o,e)

    def test_report_summary_no_projects(self):
        """report: report run with no projects in 'summary' mode
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": 87,
                       "processing_software":
                       "{ 'bcl2fastq': ('/usr/bin/bcl2fastq', 'bcl2fastq', '2.17.1.14'), 'cellranger': ('/usr/bin/cellranger', 'cellranger', '3.0.1') }",
                       "bcl2fastq_software":
                       "('/usr/bin/bcl2fastq', 'bcl2fastq', '2.17.1.14')",
                       "cellranger_software":
                       "('/usr/bin/cellranger', 'cellranger', '3.0.1')",
                       "sequencer_model": "MiSeq" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Generate summary report
        expected = """MISEQ run #87 datestamped 170901
================================
Run name  : 170901_M00879_0087_000000000-AGEW9
Reference : MISEQ_170901#87
Platform  : MISEQ
Sequencer : MiSeq
Directory : %s
Endedness : Paired end
Bcl2fastq : bcl2fastq 2.17.1.14
Cellranger: cellranger 3.0.1

No projects found; 'bcl2fastq' directory contains the following data:

- 'AB':  2 samples
- 'CDE': 2 samples""" % mockdir.dirn
        report = report_summary(ap)
        self.assertEqual(len(report.split('\n')),
                         len(expected.split('\n')))
        for o,e in zip(report.split('\n'),
                       expected.split('\n')):
            self.assertEqual(o,e)

    def test_report_summary_stored_path_differs_from_actual_path(self):
        """report: report run in 'summary' mode (stored path is wrong)
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": 87,
                       "sequencer_model": "MiSeq" },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "RNA-seq",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Sequencer model": "MiSeq" },
                "CDE": { "User": "Charles David Edwards",
                         "Library type": "ChIP-seq",
                         "Organism": "Mouse",
                         "PI": "Colin Delaney Eccleston",
                         "Sequencer model": "MiSeq",
                         "Comments": "Repeat of previous run" }
            },
            top_dir=self.dirn)
        mockdir.create()
        # Relocate (so stored paths will be wrong)
        new_dir = os.path.join(self.dirn,"elsewhere")
        os.mkdir(new_dir)
        analysis_dir_path = os.path.join(new_dir,
                                         os.path.basename(mockdir.dirn))
        os.rename(mockdir.dirn,analysis_dir_path)
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=analysis_dir_path)
        # Generate summary report
        expected = """MISEQ run #87 datestamped 170901
================================
Run name : 170901_M00879_0087_000000000-AGEW9
Reference: MISEQ_170901#87
Platform : MISEQ
Sequencer: MiSeq
Directory: %s
Endedness: Paired end
Bcl2fastq: Unknown

2 projects:
- 'AB':  Alison Bell           Human RNA-seq  2 samples (PI Audrey Bower)           
- 'CDE': Charles David Edwards Mouse ChIP-seq 2 samples (PI Colin Delaney Eccleston)

Additional notes/comments:
- CDE: Repeat of previous run
""" % analysis_dir_path
        for o,e in zip(report_summary(ap).split('\n'),
                       expected.split('\n')):
            self.assertEqual(o,e)

    def test_report_summary_with_legacy_metadata_items(self):
        """report: report run in 'summary' mode with legacy metadata items
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": 87,
                       "sequencer_model": "MiSeq" },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "RNA-seq",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Sequencer model": "MiSeq",
                        "ICELL8 well list": None },
                "CDE": { "User": "Charles David Edwards",
                         "Library type": "ChIP-seq",
                         "Organism": "Mouse",
                         "PI": "Colin Delaney Eccleston",
                         "Sequencer model": "MiSeq",
                         "Comments": "Repeat of previous run",
                         "ICELL8 well list": None }
            },
            top_dir=self.dirn)
        mockdir.create()
        # Make autoprocess instance
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Generate summary report
        expected = """MISEQ run #87 datestamped 170901
================================
Run name : 170901_M00879_0087_000000000-AGEW9
Reference: MISEQ_170901#87
Platform : MISEQ
Sequencer: MiSeq
Directory: %s
Endedness: Paired end
Bcl2fastq: Unknown

2 projects:
- 'AB':  Alison Bell           Human RNA-seq  2 samples (PI Audrey Bower)           
- 'CDE': Charles David Edwards Mouse ChIP-seq 2 samples (PI Colin Delaney Eccleston)

Additional notes/comments:
- CDE: Repeat of previous run
""" % mockdir.dirn
        for o,e in zip(report_summary(ap).split('\n'),
                       expected.split('\n')):
            self.assertEqual(o,e)

class TestReportProjects(unittest.TestCase):
    """
    Tests for the 'report' command ('projects' mode)
    """
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestProjectsSummary')
        # Store original location so we can get back at the end
        self.pwd = os.getcwd()
        # Move to working dir
        os.chdir(self.dirn)

    def tearDown(self):
        # Return to original dir
        os.chdir(self.pwd)
        # Remove the temporary test directory
        def del_rw(action,name,excinfo):
            # Explicitly remove read only files/
            # dirs
            os.chmod(os.path.dirname(name),0o755)
            os.chmod(name,0o655)
            os.remove(name)
        shutil.rmtree(self.dirn,onerror=del_rw)

    def test_report_projects(self):
        """report: report run in 'projects' mode
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": 87,
                       "sequencer_model": "MiSeq" },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "RNA-seq",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Sequencer model": "MiSeq" },
                "CDE": { "User": "Charles David Edwards",
                         "Library type": "ChIP-seq",
                         "Organism": "Mouse",
                         "PI": "Colin Delaney Eccleston",
                         "Sequencer model": "MiSeq" }
            },
            top_dir=self.dirn)
        mockdir.create()
        # Make autoprocess instance and set required metadata
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Generate projects report
        expected = """MISEQ_170901#87\t87\ttesting\t\tAlison Bell\tAudrey Bower\tRNA-seq\t\tHuman\tMISEQ\t2\t\tyes\tAB1-2
MISEQ_170901#87\t87\ttesting\t\tCharles David Edwards\tColin Delaney Eccleston\tChIP-seq\t\tMouse\tMISEQ\t2\t\tyes\tCDE3-4
"""
        for o,e in zip(report_projects(ap).split('\n'),
                       expected.split('\n')):
            self.assertEqual(o,e)

    def test_report_projects_custom_fields(self):
        """report: report run in 'projects' mode with custom fields
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": 87,
                       "sequencer_model": "MiSeq" },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "RNA-seq",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Sequencer model": "MiSeq" },
                "CDE": { "User": "Charles David Edwards",
                         "Library type": "ChIP-seq",
                         "Organism": "Mouse",
                         "PI": "Colin Delaney Eccleston",
                         "Sequencer model": "MiSeq" }
            },
            top_dir=self.dirn)
        mockdir.create()
        # Make autoprocess instance and set required metadata
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Generate projects report
        expected = """170901\tMISEQ_170901#87\t87\ttesting\t\tAB\tAlison Bell\tAudrey Bower\tRNA-seq\t\tHuman\tMISEQ\tMiSeq\t2\t\tyes\tAB1-2\t%s
170901\tMISEQ_170901#87\t87\ttesting\t\tCDE\tCharles David Edwards\tColin Delaney Eccleston\tChIP-seq\t\tMouse\tMISEQ\tMiSeq\t2\t\tyes\tCDE3-4\t%s
""" % (ap.params.analysis_dir,ap.params.analysis_dir)
        custom_fields = ['datestamp','run_id','run_number','source','null','project','user','PI','library_type','single_cell_platform','organism','platform','sequencer_model','#samples','#cells','paired_end','samples','path']
        for o,e in zip(report_projects(ap,fields=custom_fields).split('\n'),
                       expected.split('\n')):
            self.assertEqual(o,e)

    def test_report_projects_with_composite_custom_fields(self):
        """report: report run in 'projects' mode with composite custom fields
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": 87,
                       "sequencer_model": "MiSeq" },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "RNA-seq",
                        "Single cell platform": "10xGenomics Chromium 3'v3",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Sequencer model": "MiSeq" },
                "CDE": { "User": "Charles David Edwards",
                         "Library type": "ChIP-seq",
                         "Organism": "Mouse",
                         "PI": "Colin Delaney Eccleston",
                         "Sequencer model": "MiSeq" }
            },
            top_dir=self.dirn)
        mockdir.create()
        # Make autoprocess instance and set required metadata
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Generate projects report with composite field
        # (i.e. using '+' notation)
        expected = """MISEQ_170901#87\tAlison Bell\tHuman RNA-seq
MISEQ_170901#87\tCharles David Edwards\tMouse ChIP-seq
"""
        custom_fields = ['run_id','user','organism+library_type']
        for o,e in zip(report_projects(ap,fields=custom_fields).split('\n'),
                       expected.split('\n')):
            self.assertEqual(o,e)
        # Generate projects report with composite field with
        # null value in subfield
        expected = """MISEQ_170901#87\tAlison Bell\t10xGenomics Chromium 3'v3 RNA-seq
MISEQ_170901#87\tCharles David Edwards\tChIP-seq
"""
        custom_fields = ['run_id','user','single_cell_platform+library_type']
        for o,e in zip(report_projects(ap,fields=custom_fields).split('\n'),
                       expected.split('\n')):
            self.assertEqual(o,e)
        # Generate projects report with composite field with
        # alternative delimiter (i.e. using '[...]:' notation)
        expected = """MISEQ_170901#87\tAlison Bell\ttesting_87
MISEQ_170901#87\tCharles David Edwards\ttesting_87
"""
        custom_fields = ['run_id','user','[_]:source+run_number']
        for o,e in zip(report_projects(ap,fields=custom_fields).split('\n'),
                       expected.split('\n')):
            self.assertEqual(o,e)

    def test_report_projects_novaseq_flow_cell_mode(self):
        """report: report run in 'projects' mode with NovaSeq flow cell mode
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '221216_A00879_0047_000000000-AGEW9',
            'novaseq',
            metadata={ "source": "testing",
                       "run_number": 47,
                       "sequencer_model": "NovaSeq 6000",
                       "flow_cell_mode": "SP" },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "RNA-seq",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Sequencer model": "NovaSeq 6000" },
                "CDE": { "User": "Charles David Edwards",
                         "Library type": "ChIP-seq",
                         "Organism": "Mouse",
                         "PI": "Colin Delaney Eccleston",
                         "Sequencer model": "NovaSeq 6000" }
            },
            top_dir=self.dirn)
        mockdir.create()
        # Make autoprocess instance and set required metadata
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Generate projects report
        expected = """221216\tNOVASEQ_221216#47\t47\tSP\ttesting\t\tAB\tAlison Bell\tAudrey Bower\tRNA-seq\t\tHuman\tNOVASEQ\tNovaSeq 6000\t2\t\tyes\tAB1-2\t%s
221216\tNOVASEQ_221216#47\t47\tSP\ttesting\t\tCDE\tCharles David Edwards\tColin Delaney Eccleston\tChIP-seq\t\tMouse\tNOVASEQ\tNovaSeq 6000\t2\t\tyes\tCDE3-4\t%s
""" % (ap.params.analysis_dir,ap.params.analysis_dir)
        custom_fields = ['datestamp','run_id','run_number','flow_cell_mode','source','null','project','user','PI','library_type','single_cell_platform','organism','platform','sequencer_model','#samples','#cells','paired_end','samples','path']
        for o,e in zip(report_projects(ap,fields=custom_fields).split('\n'),
                       expected.split('\n')):
            self.assertEqual(o,e)

    def test_report_projects_single_cell(self):
        """report: report single-cell run in 'projects' mode
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": 87,
                       "sequencer_model": "MiSeq" },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "scRNA-seq",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Single cell platform": "10xGenomics Chromium 3'",
                        "Number of cells": 1311,
                        "Sequencer model": "MiSeq" },
                "CDE": { "User": "Charles David Edwards",
                         "Library type": "ChIP-seq",
                         "Organism": "Mouse",
                         "PI": "Colin Delaney Eccleston",
                         "Sequencer model": "MiSeq" }
            },
            top_dir=self.dirn)
        mockdir.create()
        # Make autoprocess instance and set required metadata
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Generate projects report
        expected = """MISEQ_170901#87\t87\ttesting\t\tAlison Bell\tAudrey Bower\tscRNA-seq\t10xGenomics Chromium 3'\tHuman\tMISEQ\t2\t1311\tyes\tAB1-2
MISEQ_170901#87\t87\ttesting\t\tCharles David Edwards\tColin Delaney Eccleston\tChIP-seq\t\tMouse\tMISEQ\t2\t\tyes\tCDE3-4
"""
        for o,e in zip(report_projects(ap).split('\n'),
                       expected.split('\n')):
            self.assertEqual(o,e)

    def test_report_projects_10x_cellplex(self):
        """report: report 10xGenomics CellPlex run in 'projects' mode
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": 87,
                       "sequencer_model": "MiSeq" },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "CellPlex",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Single cell platform": "10xGenomics Chromium 3'v3",
                        "Number of cells": 1311,
                        "Sequencer model": "MiSeq" },
                "CDE": { "User": "Charles David Edwards",
                         "Library type": "ChIP-seq",
                         "Organism": "Mouse",
                         "PI": "Colin Delaney Eccleston",
                         "Sequencer model": "MiSeq" }
            },
            top_dir=self.dirn)
        mockdir.create()
        # Add a cellranger multi config.csv file
        with open(os.path.join(mockdir.dirn,
                               "AB",
                               "10x_multi_config.csv"),'wt') as fp:
            fastq_dir = os.path.join(mockdir.dirn,
                                     "AB",
                                     "fastqs")
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
AB1,%s,any,AB1,gene expression,
AB2,%s,any,AB2,Multiplexing Capture,

[samples]
sample_id,cmo_ids,description
ABM1,CMO301,ABM1
ABM2,CMO302,ABM2
ABM3,CMO303,ABM3
ABM4,CMO304,ABM4
""" % (fastq_dir,fastq_dir))
        # Make autoprocess instance and set required metadata
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Generate projects report
        expected = """MISEQ_170901#87\t87\ttesting\t\tAlison Bell\tAudrey Bower\tCellPlex\t10xGenomics Chromium 3'v3\tHuman\tMISEQ\t4\t1311\tyes\tABM1-4
MISEQ_170901#87\t87\ttesting\t\tCharles David Edwards\tColin Delaney Eccleston\tChIP-seq\t\tMouse\tMISEQ\t2\t\tyes\tCDE3-4
"""
        for o,e in zip(report_projects(ap).split('\n'),
                       expected.split('\n')):
            self.assertEqual(o,e)

    def test_report_projects_10x_cellplex_scrnaseq(self):
        """report: report 10xGenomics CellPlex scRNA-seq run in 'projects' mode
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": 87,
                       "sequencer_model": "MiSeq" },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "CellPlex scRNA-seq",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Single cell platform": "10xGenomics Chromium 3'v3",
                        "Number of cells": 1311,
                        "Sequencer model": "MiSeq" },
                "CDE": { "User": "Charles David Edwards",
                         "Library type": "ChIP-seq",
                         "Organism": "Mouse",
                         "PI": "Colin Delaney Eccleston",
                         "Sequencer model": "MiSeq" }
            },
            top_dir=self.dirn)
        mockdir.create()
        # Add a cellranger multi config.csv file
        with open(os.path.join(mockdir.dirn,
                               "AB",
                               "10x_multi_config.csv"),'wt') as fp:
            fastq_dir = os.path.join(mockdir.dirn,
                                     "AB",
                                     "fastqs")
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
AB1,%s,any,AB1,gene expression,
AB2,%s,any,AB2,Multiplexing Capture,

[samples]
sample_id,cmo_ids,description
ABM1,CMO301,ABM1
ABM2,CMO302,ABM2
ABM3,CMO303,ABM3
ABM4,CMO304,ABM4
""" % (fastq_dir,fastq_dir))
        # Make autoprocess instance and set required metadata
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Generate projects report
        expected = """MISEQ_170901#87\t87\ttesting\t\tAlison Bell\tAudrey Bower\tCellPlex scRNA-seq\t10xGenomics Chromium 3'v3\tHuman\tMISEQ\t4\t1311\tyes\tABM1-4
MISEQ_170901#87\t87\ttesting\t\tCharles David Edwards\tColin Delaney Eccleston\tChIP-seq\t\tMouse\tMISEQ\t2\t\tyes\tCDE3-4
"""
        for o,e in zip(report_projects(ap).split('\n'),
                       expected.split('\n')):
            self.assertEqual(o,e)

    def test_report_projects_10x_flex(self):
        """report: report 10xGenomics Flex run in 'projects' mode
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '240213_A00879_0087_000000000-AGEW9',
            'novaseq',
            metadata={ "source": "testing",
                       "run_number": 87,
                       "sequencer_model": "MiSeq" },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "Flex",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Single cell platform": "10xGenomics Chromium 3'v3",
                        "Number of cells": 1311,
                        "Sequencer model": "NovaSeq" },
                "CDE": { "User": "Charles David Edwards",
                         "Library type": "ChIP-seq",
                         "Organism": "Mouse",
                         "PI": "Colin Delaney Eccleston",
                         "Sequencer model": "NovaSeq" }
            },
            top_dir=self.dirn)
        mockdir.create()
        # Add a cellranger multi config.csv file
        with open(os.path.join(mockdir.dirn,
                               "AB",
                               "10x_multi_config.csv"),'wt') as fp:
            fastq_dir = os.path.join(mockdir.dirn,
                                     "AB",
                                     "fastqs")
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A
probe-set,/data/Probe_Set_v1.0_GRCh38-2020-A.csv
no-bam,true

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB1_Flex,{fastq_dir},any,PJB1,Gene Expression,

[samples]
sample_id,probe_barcode_ids,description
ABF1,BC001,ABF1
ABF2,BC002,ABF2
ABF3,BC003,ABF3
ABF4,BC004,ABF4
""".format(fastq_dir=fastq_dir))
        # Make autoprocess instance and set required metadata
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Generate projects report
        expected = """NOVASEQ_240213#87\t87\ttesting\t\tAlison Bell\tAudrey Bower\tFlex\t10xGenomics Chromium 3'v3\tHuman\tNOVASEQ\t4\t1311\tyes\tABF1-4
NOVASEQ_240213#87\t87\ttesting\t\tCharles David Edwards\tColin Delaney Eccleston\tChIP-seq\t\tMouse\tNOVASEQ\t2\t\tyes\tCDE3-4
"""
        for o,e in zip(report_projects(ap).split('\n'),
                       expected.split('\n')):
            self.assertEqual(o,e)

    def test_report_projects_multiplexed_samples(self):
        """report: report run with multiplexed samples in 'projects' mode
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": 87,
                       "sequencer_model": "MiSeq" },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "RNA-seq",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Sequencer model": "MiSeq",
                        "Multiplexed samples": "M1,M2,M3,M4" },
                "CDE": { "User": "Charles David Edwards",
                         "Library type": "ChIP-seq",
                         "Organism": "Mouse",
                         "PI": "Colin Delaney Eccleston",
                         "Sequencer model": "MiSeq" }
            },
            top_dir=self.dirn)
        mockdir.create()
        # Make autoprocess instance and set required metadata
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Generate projects report
        expected = """MISEQ_170901#87\t87\ttesting\t\tAlison Bell\tAudrey Bower\tRNA-seq\t\tHuman\tMISEQ\t4\t\tyes\tM1-4
MISEQ_170901#87\t87\ttesting\t\tCharles David Edwards\tColin Delaney Eccleston\tChIP-seq\t\tMouse\tMISEQ\t2\t\tyes\tCDE3-4
"""
        for o,e in zip(report_projects(ap).split('\n'),
                       expected.split('\n')):
            self.assertEqual(o,e)

    def test_report_projects_unknown_multiplexed_samples(self):
        """report: report run with "unknown" multiplexed samples in 'projects' mode
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": 87,
                       "sequencer_model": "MiSeq" },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "RNA-seq",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Sequencer model": "MiSeq",
                        "Multiplexed samples": "?" },
                "CDE": { "User": "Charles David Edwards",
                         "Library type": "ChIP-seq",
                         "Organism": "Mouse",
                         "PI": "Colin Delaney Eccleston",
                         "Sequencer model": "MiSeq" }
            },
            top_dir=self.dirn)
        mockdir.create()
        # Make autoprocess instance and set required metadata
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Generate projects report
        expected = """MISEQ_170901#87\t87\ttesting\t\tAlison Bell\tAudrey Bower\tRNA-seq\t\tHuman\tMISEQ\t?\t\tyes\t?
MISEQ_170901#87\t87\ttesting\t\tCharles David Edwards\tColin Delaney Eccleston\tChIP-seq\t\tMouse\tMISEQ\t2\t\tyes\tCDE3-4
"""
        for o,e in zip(report_projects(ap).split('\n'),
                       expected.split('\n')):
            self.assertEqual(o,e)

    def test_report_no_projects(self):
        """report: report run with no projects in 'projects' mode
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": 87,
                       "sequencer_model": "MiSeq" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Make autoprocess instance and set required metadata
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Generate projects report
        expected = ""
        for o,e in zip(report_projects(ap).split('\n'),
                       expected.split('\n')):
            self.assertEqual(o,e)

    def test_report_projects_stored_path_differs_from_actual_path(self):
        """report: report run in 'projects' mode (stored path is wrong)
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": 87,
                       "sequencer_model": "MiSeq" },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "RNA-seq",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Sequencer model": "MiSeq" },
                "CDE": { "User": "Charles David Edwards",
                         "Library type": "ChIP-seq",
                         "Organism": "Mouse",
                         "PI": "Colin Delaney Eccleston",
                         "Sequencer model": "MiSeq" }
            },
            top_dir=self.dirn)
        mockdir.create()
        # Relocate (so stored paths will be wrong)
        new_dir = os.path.join(self.dirn,"elsewhere")
        os.mkdir(new_dir)
        analysis_dir_path = os.path.join(new_dir,
                                         os.path.basename(mockdir.dirn))
        os.rename(mockdir.dirn,analysis_dir_path)
        # Make autoprocess instance and set required metadata
        ap = AutoProcess(analysis_dir=analysis_dir_path)
        # Custom fields (including path)
        custom_fields = ("run_id","project","path")
        # Generate projects report
        expected = """MISEQ_170901#87\tAB\t{path}
MISEQ_170901#87\tCDE\t{path}
""".format(path=analysis_dir_path)
        for o,e in zip(report_projects(ap,fields=custom_fields).split('\n'),
                       expected.split('\n')):
            self.assertEqual(o,e)

    def test_report_projects_with_legacy_metadata(self):
        """report: report run in 'projects' mode with legacy metadata items
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": 87,
                       "sequencer_model": "MiSeq" },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "RNA-seq",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Sequencer model": "MiSeq",
                        "ICELL8 well list": None },
                "CDE": { "User": "Charles David Edwards",
                         "Library type": "ChIP-seq",
                         "Organism": "Mouse",
                         "PI": "Colin Delaney Eccleston",
                         "Sequencer model": "MiSeq",
                         "ICELL8 well list": None }
            },
            top_dir=self.dirn)
        mockdir.create()
        # Make autoprocess instance and set required metadata
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Generate projects report
        expected = """MISEQ_170901#87\t87\ttesting\t\tAlison Bell\tAudrey Bower\tRNA-seq\t\tHuman\tMISEQ\t2\t\tyes\tAB1-2
MISEQ_170901#87\t87\ttesting\t\tCharles David Edwards\tColin Delaney Eccleston\tChIP-seq\t\tMouse\tMISEQ\t2\t\tyes\tCDE3-4
"""
        for o,e in zip(report_projects(ap).split('\n'),
                       expected.split('\n')):
            self.assertEqual(o,e)

class TestReport(unittest.TestCase):
    """
    Tests for the 'report' command invoked directly
    """
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestProjectsSummary')
        # Store original location so we can get back at the end
        self.pwd = os.getcwd()
        # Move to working dir
        os.chdir(self.dirn)

    def tearDown(self):
        # Return to original dir
        os.chdir(self.pwd)
        # Remove the temporary test directory
        def del_rw(action,name,excinfo):
            # Explicitly remove read only files/
            # dirs
            os.chmod(os.path.dirname(name),0o755)
            os.chmod(name,0o655)
            os.remove(name)
        shutil.rmtree(self.dirn,onerror=del_rw)

    def test_report_to_file(self):
        """report: report run to file
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": 87,
                       "sequencer_model": "MiSeq"},
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "RNA-seq",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Sequencer model": "MiSeq" },
                "CDE": { "User": "Charles David Edwards",
                         "Library type": "ChIP-seq",
                         "Organism": "Mouse",
                         "PI": "Colin Delaney Eccleston",
                         "Sequencer model": "MiSeq" }
            },
            top_dir=self.dirn)
        mockdir.create()
        # Make autoprocess instance and set required metadata
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Generate projects report
        out_file = os.path.join(self.dirn,"projects.tsv")
        report(ap,mode=ReportingMode.PROJECTS,out_file=out_file)
        # Check the outputs
        expected = """MISEQ_170901#87\t87\ttesting\t\tAlison Bell\tAudrey Bower\tRNA-seq\t\tHuman\tMISEQ\t2\t\tyes\tAB1-2
MISEQ_170901#87\t87\ttesting\t\tCharles David Edwards\tColin Delaney Eccleston\tChIP-seq\t\tMouse\tMISEQ\t2\t\tyes\tCDE3-4
"""
        self.assertTrue(os.path.exists(out_file))
        with open(out_file,'r') as fp:
            report_contents = fp.read()
            for o,e in zip(report_contents.split('\n'),
                           expected.split('\n')):
                self.assertEqual(o,e)

class TestFetchValueFunction(unittest.TestCase):
    """
    Tests for the 'fetch_value' function
    """
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestProjectsSummary')
        # Store original location so we can get back at the end
        self.pwd = os.getcwd()
        # Move to working dir
        os.chdir(self.dirn)

    def tearDown(self):
        # Return to original dir
        os.chdir(self.pwd)
        # Remove the temporary test directory
        def del_rw(action,name,excinfo):
            # Explicitly remove read only files/
            # dirs
            os.chmod(os.path.dirname(name),0o755)
            os.chmod(name,0o655)
            os.remove(name)
        shutil.rmtree(self.dirn,onerror=del_rw)

    def test_fetch_value(self):
        """
        report: test the fetch_value function
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": 87,
                       "sequencer_model": "MiSeq",
                       "analysis_number": 2 },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "scRNA-seq",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Single cell platform": "10xGenomics Chromium 3'",
                        "Number of cells": 1311,
                        "Sequencer model": "MiSeq" },
                "CDE": { "User": "Charles David Edwards",
                         "Library type": "ChIP-seq",
                         "Organism": "Mouse",
                         "PI": "Colin Delaney Eccleston",
                         "Sequencer model": "MiSeq" }
            },
            top_dir=self.dirn)
        mockdir.create()
        # Make autoprocess instance and set required metadata
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Acquire project
        project = AnalysisProject('AB',
                                  os.path.join(mockdir.dirn,'AB'))
        # Check the returned values
        self.assertEqual(fetch_value(ap,project,'datestamp'),'170901')
        self.assertEqual(fetch_value(ap,project,'run_id'),'MISEQ_170901#87.2')
        self.assertEqual(fetch_value(ap,project,'run_number'),'87')
        self.assertEqual(fetch_value(ap,project,'analysis_number'),'2')
        self.assertEqual(fetch_value(ap,project,'source'),'testing')
        self.assertEqual(fetch_value(ap,project,'data_source'),'testing')
        self.assertEqual(fetch_value(ap,project,'analysis_dir'),mockdir.dirn)
        self.assertEqual(fetch_value(ap,project,'path'),mockdir.dirn)
        self.assertEqual(fetch_value(ap,project,'project'),'AB')
        self.assertEqual(fetch_value(ap,project,'project_name'),'AB')
        self.assertEqual(fetch_value(ap,project,'user'),'Alison Bell')
        self.assertEqual(fetch_value(ap,project,'PI'),'Audrey Bower')
        self.assertEqual(fetch_value(ap,project,'pi'),'Audrey Bower')
        self.assertEqual(fetch_value(ap,project,'application'),'scRNA-seq')
        self.assertEqual(fetch_value(ap,project,'library_type'),'scRNA-seq')
        self.assertEqual(fetch_value(ap,project,'organism'),'Human')
        self.assertEqual(fetch_value(ap,project,'sequencer_model'),'MiSeq')
        self.assertEqual(fetch_value(ap,project,'sequencer_platform'),'MISEQ')
        self.assertEqual(fetch_value(ap,project,'platform'),'MISEQ')
        self.assertEqual(fetch_value(ap,project,'no_of_samples'),'2')
        self.assertEqual(fetch_value(ap,project,'#samples'),'2')
        self.assertEqual(fetch_value(ap,project,'single_cell_platform'),
                         "10xGenomics Chromium 3'")
        self.assertEqual(fetch_value(ap,project,'no_of_cells'),'1311')
        self.assertEqual(fetch_value(ap,project,'#cells'),'1311')
        self.assertEqual(fetch_value(ap,project,'paired_end'),'yes')
        self.assertEqual(fetch_value(ap,project,'samples'),'AB1-2')
        self.assertEqual(fetch_value(ap,project,'sample_names'),'AB1-2')
        self.assertEqual(fetch_value(ap,project,'null'),'')

class TestDefaultValueFunction(unittest.TestCase):
    """
    Tests for the 'default_value' helper function
    """
    def test_default_value(self):
        """report: test the default_value function
        """
        self.assertEqual("Hello",default_value("Hello"))
        self.assertEqual("Hello",default_value("Hello",default="Goodbye"))
        self.assertEqual("",default_value(None))
        self.assertEqual("Goodbye",default_value(None,default="Goodbye"))
        self.assertEqual(0,default_value(0))

class TestGetMultiplexedSamplesFunction(unittest.TestCase):
    """
    Tests for the 'get_multiplexed_samples' function
    """
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestGetMultiplexedSamples')
        # Store original location so we can get back at the end
        self.pwd = os.getcwd()
        # Move to working dir
        os.chdir(self.dirn)

    def tearDown(self):
        # Return to original dir
        os.chdir(self.pwd)
        # Remove the temporary test directory
        def del_rw(action,name,excinfo):
            # Explicitly remove read only files/
            # dirs
            os.chmod(os.path.dirname(name),0o755)
            os.chmod(name,0o655)
            os.remove(name)
        shutil.rmtree(self.dirn,onerror=del_rw)

    def test_get_multiplexed_samples_no_samples(self):
        """
        report: test 'get_multiplexed_samples' (no multiplexed samples)
        """
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": 87,
                       "sequencer_model": "MiSeq",
                       "analysis_number": 2 },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "scRNA-seq",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Single cell platform": None,
                        "Number of cells": None,
                        "Sequencer model": "MiSeq" }
            },
            top_dir=self.dirn)
        mockdir.create()
        project = AnalysisProject(os.path.join(mockdir.dirn,'AB'))
        self.assertEqual(get_multiplexed_samples(project), None)

    def test_get_multiplexed_samples_missing_samples(self):
        """
        report: test 'get_multiplexed_samples' (missing multiplexed samples)
        """
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": 87,
                       "sequencer_model": "MiSeq",
                       "analysis_number": 2 },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "CellPlex scRNA-seq",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Single cell platform": "10xGenomics Chromium 3'v3",
                        "Number of cells": None,
                        "Sequencer model": "MiSeq" }
            },
            top_dir=self.dirn)
        mockdir.create()
        project = AnalysisProject(os.path.join(mockdir.dirn,'AB'))
        self.assertEqual(get_multiplexed_samples(project), [])

    def test_get_multiplexed_samples_10x_multi_config(self):
        """
        report: test 'get_multiplexed_samples' (use 10x multi config)
        """
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": 87,
                       "sequencer_model": "MiSeq",
                       "analysis_number": 2 },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "CellPlex scRNA-seq",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Single cell platform": "10xGenomics Chromium 3'v3",
                        "Number of cells": None,
                        "Sequencer model": "MiSeq" }
            },
            top_dir=self.dirn)
        mockdir.create()
        # Add a cellranger multi config.csv file
        with open(os.path.join(mockdir.dirn,
                               "AB",
                               "10x_multi_config.csv"),'wt') as fp:
            fastq_dir = os.path.join(mockdir.dirn,
                                     "AB",
                                     "fastqs")
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
AB1,%s,any,AB1,gene expression,
AB2,%s,any,AB2,Multiplexing Capture,

[samples]
sample_id,cmo_ids,description
ABM1,CMO301,ABM1
ABM2,CMO302,ABM2
ABM3,CMO303,ABM3
ABM4,CMO304,ABM4
""" % (fastq_dir,fastq_dir))
        project = AnalysisProject(os.path.join(mockdir.dirn,'AB'))
        self.assertEqual(get_multiplexed_samples(project),
                         ["ABM1", "ABM2", "ABM3", "ABM4"])

    def test_get_multiplexed_samples_multiple_10x_multi_configs(self):
        """
        report: test 'get_multiplexed_samples' (multiple 10x multi configs)
        """
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": 87,
                       "sequencer_model": "MiSeq",
                       "analysis_number": 2 },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "CellPlex scRNA-seq",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Single cell platform": "10xGenomics Chromium 3'v3",
                        "Number of cells": None,
                        "Sequencer model": "MiSeq" }
            },
            top_dir=self.dirn)
        mockdir.create()
        # Add cellranger multi config.csv files
        with open(os.path.join(mockdir.dirn,
                               "AB",
                               "10x_multi_config.AB1.csv"),'wt') as fp:
            fastq_dir = os.path.join(mockdir.dirn,
                                     "AB",
                                     "fastqs")
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
AB1,%s,any,AB1,gene expression,
AB1,%s,any,AB1,Multiplexing Capture,

[samples]
sample_id,cmo_ids,description
ABM1,CMO301,ABM1
ABM2,CMO302,ABM2
""" % (fastq_dir,fastq_dir))
        with open(os.path.join(mockdir.dirn,
                               "AB",
                               "10x_multi_config.AB2.csv"),'wt') as fp:
            fastq_dir = os.path.join(mockdir.dirn,
                                     "AB",
                                     "fastqs")
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
AB2,%s,any,AB2,gene expression,
AB2,%s,any,AB2,Multiplexing Capture,

[samples]
sample_id,cmo_ids,description
ABM3,CMO303,ABM3
ABM4,CMO304,ABM4
""" % (fastq_dir,fastq_dir))
        project = AnalysisProject(os.path.join(mockdir.dirn,'AB'))
        self.assertEqual(get_multiplexed_samples(project),
                         ["ABM1", "ABM2", "ABM3", "ABM4"])
