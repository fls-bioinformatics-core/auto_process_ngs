#######################################################################
# Tests for autoprocessor.py module
#######################################################################

import unittest
import tempfile
import shutil
import os
from auto_process_ngs.auto_processor import AutoProcess
from auto_process_ngs.mock import MockAnalysisDirFactory

IEMSampleSheetContents = """[Header]
IEMFileVersion,4
Date,21/9/2015
Workflow,GenerateFASTQ
Application,FASTQ Only
Assay,Nextera XT
Description,
Chemistry,Amplicon

[Reads]
75
75

[Settings]
ReverseComplement,0
Adapter,CTGTCTCTTATACACATCT

[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
S1_A,S1_A,,,N701,TAAGGCGA,S501,TAGATCGC,,
S2_A,S2_A,,,N702,CGTACTAG,S501,TAGATCGC,,
TM3,TM3,,,N703,AGGCAGAA,S502,CTCTCTAT,,
"""

RunInfoXmlContents = """<?xml version="1.0"?>
<RunInfo xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" Version="2">
  <Run Id="150921_M00123_0045_000000000-ABC6D" Number="45">
    <Flowcell>000000000-ABC6D</Flowcell>
    <Instrument>M00123</Instrument>
    <Date>150921</Date>
    <Reads>
      <Read NumCycles="75" Number="1" IsIndexedRead="N" />
      <Read NumCycles="8" Number="2" IsIndexedRead="Y" />
      <Read NumCycles="8" Number="3" IsIndexedRead="Y" />
      <Read NumCycles="75" Number="4" IsIndexedRead="N" />
    </Reads>
    <FlowcellLayout LaneCount="1" SurfaceCount="2" SwathCount="1" TileCount="19" />
  </Run>
</RunInfo>
"""

# Helper functions for testing

def _make_mock_primary_data_dir(top_dir):
    # Creates a mock primary data directory
    data_dir = os.path.join(top_dir,
                            '150921_M00123_0045_000000000-ABC6D')
    # Make the Data/Intensities/BaseCalls dir
    os.makedirs(os.path.join(data_dir,
                             'Data','Intensities','BaseCalls'))
    # Make an IEM sample sheet
    sample_sheet = os.path.join(data_dir,
                                'Data','Intensities','BaseCalls',
                                'SampleSheet.csv')
    with open(sample_sheet,'w') as fp:
        fp.write(IEMSampleSheetContents)
    # Make a RunInfo.xml file
    run_info_xml =  os.path.join(data_dir,'RunInfo.xml')
    with open(run_info_xml,'w') as fp:
        fp.write(RunInfoXmlContents)
    return data_dir

# Unit tests

class TestAutoProcess(unittest.TestCase):

    def setUp(self):
        # Create a temporary working dir for tests
        self.this_dir = os.getcwd()
        self.wd = tempfile.mkdtemp(suffix='TestAutoProcess')
        self.data_dir = _make_mock_primary_data_dir(self.wd)
        os.chdir(self.wd)

    def tearDown(self):
        # Remove the temporary test directory
        os.chdir(self.this_dir)
        shutil.rmtree(self.wd)

    def test_autoprocess_setup(self):
        # Test the 'setup' command of AutoProcess
        ap = AutoProcess()
        ap.setup(self.data_dir)
        self.assertTrue(os.path.isdir(
            '150921_M00123_0045_000000000-ABC6D_analysis'))

class TestAutoProcessImportProject(unittest.TestCase):
    """Tests for AutoProcess.import_project

    """
    def setUp(self):
        self.dirn = tempfile.mkdtemp(suffix='TestAutoProcessImportProject')

    def tearDown(self):
        # Remove the temporary test directory
        shutil.rmtree(self.dirn)

    def test_import_project(self):
        """Check AutoProcess.import_project
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_M00879_0087_000000000-AGEW9',
            'miseq',
            top_dir=self.dirn)
        mockdir.create()
        open(os.path.join(mockdir.dirn,'auto_process.info'),'w').write(
            """analysis_dir\t%s
bases_mask\ty76,I8,I8,y76
data_dir\t/mnt/data/%s
per_lane_stats_file\tper_lane_statistics.info
primary_data_dir\t%s/primary_data/%s
project_metadata\tprojects.info
sample_sheet\t%s/custom_SampleSheet.csv
stats_file\tstatistics.info
unaligned_dir\tbcl2fastq
"""
            % (mockdir.dirn,
               os.path.basename(mockdir.dirn[:-9]),
               mockdir.dirn,os.path.basename(mockdir.dirn[:-9]),
               mockdir.dirn))
        # Make a mock project
        project_dir = os.path.join(self.dirn,'NewProj')
        os.mkdir(project_dir)
        os.mkdir(os.path.join(project_dir,'fastqs'))
        for fq in ('NP01_S1_R1_001.fastq.gz','NP01_S1_R1_001.fastq.gz'):
            open(os.path.join(project_dir,'fastqs',fq),'w').write('')
        open(os.path.join(project_dir,'README.info'),'w').write(
            """Run\t160622_NB5001234_0011_ABCDE5AFXX
Platform\tnextseq
User\tPeter Briggs
PI\tAnne Cleaves
Organism\tHuman
Library type\tRNA-seq
Paired_end\tY
Samples\t1 sample (NP01)
Comments\t1% PhiX spike in
""")
        # Check that the project is not currently present
        ap = AutoProcess(mockdir.dirn)
        self.assertFalse('NewProj' in [p.name
                                       for p in ap.get_analysis_projects()])
        self.assertFalse('NewProj' in [p.name
                                       for p in ap.get_analysis_projects_from_dirs()])
        self.assertFalse(os.path.exists(os.path.join(ap.analysis_dir,'NewProj')))
        # Import the project
        ap.import_project(project_dir)
        self.assertTrue('NewProj' in [p.name
                                      for p in ap.get_analysis_projects()])
        self.assertTrue('NewProj' in [p.name
                                      for p in ap.get_analysis_projects_from_dirs()])
        self.assertTrue(os.path.exists(os.path.join(ap.analysis_dir,'NewProj')))
