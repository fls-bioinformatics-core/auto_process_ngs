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
        # Generate concise report
        expected = """Run reference: MISEQ_170901#87
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
                        "Single cell platform": "ICELL8",
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
        # Generate concise report
        expected = """Run reference: MISEQ_170901#87
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
  SC Plat.: ICELL8
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
        # Generate concise report
        expected = """Run reference: MISEQ_170901#87
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
        # Generate concise report
        expected = """Run reference: MISEQ_170901#87
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
        # Generate concise report
        expected = """Run reference: MISEQ_170901#87
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
                        "Single cell platform": "ICELL8",
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
                         "Paired end: 'AB': Alison Bell, Human ICELL8 scRNA-seq (PI: Audrey Bower) (2 samples/1311 cells); 'CDE': Charles David Edwards, Mouse ChIP-seq (PI: Colin Delaney Eccleston) (2 samples)")

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
                         "Paired end: 'AB': Alison Bell, Human 10xGenomics Chromium 3'v3 CellPlex (PI: Audrey Bower) (4 samples/1311 cells); 'CDE': Charles David Edwards, Mouse ChIP-seq (PI: Colin Delaney Eccleston) (2 samples)")

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
                        "Single cell platform": "ICELL8",
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
- 'AB':  Alison Bell           Human scRNA-seq (ICELL8) 2 samples/1311 cells (PI Audrey Bower)           
- 'CDE': Charles David Edwards Mouse ChIP-seq           2 samples            (PI Colin Delaney Eccleston)

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
- 'AB':  Alison Bell           Human CellPlex (10xGenomics Chromium 3'v3) 4 samples/1311 cells (PI Audrey Bower)           
- 'CDE': Charles David Edwards Mouse ChIP-seq                             2 samples            (PI Colin Delaney Eccleston)

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
                        "Single cell platform": "ICELL8",
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
        expected = """MISEQ_170901#87\t87\ttesting\t\tAlison Bell\tAudrey Bower\tscRNA-seq\tICELL8\tHuman\tMISEQ\t2\t1311\tyes\tAB1-2
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
                       "sequencer_model": "MiSeq" },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "scRNA-seq",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Single cell platform": "ICELL8",
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
        self.assertEqual(fetch_value(ap,project,'run_id'),'MISEQ_170901#87')
        self.assertEqual(fetch_value(ap,project,'run_number'),'87')
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
        self.assertEqual(fetch_value(ap,project,'single_cell_platform'),'ICELL8')
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
