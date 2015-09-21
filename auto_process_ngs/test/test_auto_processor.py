#######################################################################
# Tests for autoprocessor.py module
#######################################################################

import unittest
import tempfile
import shutil
import os
from auto_process_ngs.auto_processor import AutoProcess

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
