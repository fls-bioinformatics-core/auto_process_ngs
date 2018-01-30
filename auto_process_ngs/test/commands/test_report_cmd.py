#######################################################################
# Tests for report_cmd.py module
#######################################################################

import unittest
import tempfile
import shutil
import os
from auto_process_ngs.auto_processor import AutoProcess
from auto_process_ngs.mock import MockAnalysisDirFactory
from auto_process_ngs.commands.report_cmd import report_info
from auto_process_ngs.commands.report_cmd import report_concise
from auto_process_ngs.commands.report_cmd import report_summary
from auto_process_ngs.commands.report_cmd import report_projects

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
            os.chmod(os.path.dirname(name),0755)
            os.chmod(name,0655)
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
                       "run_number": 87,
                       "assay": "Nextera" },
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
                       "run_number": 87,
                       "assay": "Nextera" },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "scRNA-seq",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Single cell platform": "ICELL8"
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
            os.chmod(os.path.dirname(name),0755)
            os.chmod(name,0655)
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
                       "assay": "Nextera" },
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
                       "run_number": 87,
                       "assay": "Nextera" },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "scRNA-seq",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Single cell platform": "ICELL8" },
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
                         "Paired end: 'AB': Alison Bell, Human ICELL8 scRNA-seq (PI: Audrey Bower) (2 samples); 'CDE': Charles David Edwards, Mouse ChIP-seq (PI: Colin Delaney Eccleston) (2 samples)")

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
            os.chmod(os.path.dirname(name),0755)
            os.chmod(name,0655)
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
                       "assay": "Nextera" },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "RNA-seq",
                        "Organism": "Human",
                        "PI": "Audrey Bower" },
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
Run name : 170901_M00879_0087_000000000-AGEW9
Reference: MISEQ_170901#87
Platform : MISEQ
Directory: %s
Endedness: Paired end
Bcl2fastq: Unknown
Assay    : Nextera

2 projects:
- 'AB':  Alison Bell           Human RNA-seq   2 samples (PI Audrey Bower)           
- 'CDE': Charles David Edwards Mouse ChIP-seq  2 samples (PI Colin Delaney Eccleston)

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
                       "assay": "Nextera" },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "scRNA-seq",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Single cell platform": "ICELL8" },
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
Run name : 170901_M00879_0087_000000000-AGEW9
Reference: MISEQ_170901#87
Platform : MISEQ
Directory: %s
Endedness: Paired end
Bcl2fastq: Unknown
Assay    : Nextera

2 projects:
- 'AB':  Alison Bell           Human scRNA-seq (ICELL8) 2 samples (PI Audrey Bower)           
- 'CDE': Charles David Edwards Mouse ChIP-seq           2 samples (PI Colin Delaney Eccleston)

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
            os.chmod(os.path.dirname(name),0755)
            os.chmod(name,0655)
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
                       "assay": "Nextera" },
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
        # Make autoprocess instance and set required metadata
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Generate projects report
        expected = """2 projects found
MISEQ_170901#87\t87\ttesting\t\tAlison Bell\tAudrey Bower\tRNA-seq\t\tHuman\tMISEQ\t2\tyes\tAB1-2
MISEQ_170901#87\t87\ttesting\t\tCharles David Edwards\tColin Delaney Eccleston\tChIP-seq\t\tMouse\tMISEQ\t2\tyes\tCDE3-4
"""
        for o,e in zip(report_projects(ap).split('\n'),
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
                       "assay": "Nextera" },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "scRNA-seq",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Single cell platform": "ICELL8" },
                "CDE": { "User": "Charles David Edwards",
                         "Library type": "ChIP-seq",
                         "Organism": "Mouse",
                         "PI": "Colin Delaney Eccleston" }
            },
            top_dir=self.dirn)
        mockdir.create()
        # Make autoprocess instance and set required metadata
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Generate projects report
        expected = """2 projects found
MISEQ_170901#87\t87\ttesting\t\tAlison Bell\tAudrey Bower\tscRNA-seq\tICELL8\tHuman\tMISEQ\t2\tyes\tAB1-2
MISEQ_170901#87\t87\ttesting\t\tCharles David Edwards\tColin Delaney Eccleston\tChIP-seq\t\tMouse\tMISEQ\t2\tyes\tCDE3-4
"""
        for o,e in zip(report_projects(ap).split('\n'),
                       expected.split('\n')):
            self.assertEqual(o,e)
