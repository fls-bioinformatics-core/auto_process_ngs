#######################################################################
# Tests for barcodes/analysis.py module
#######################################################################

import os
import unittest
import tempfile
import shutil
from auto_process_ngs.barcodes.analysis import BarcodeCounter
from auto_process_ngs.barcodes.analysis import BarcodeGroup
from auto_process_ngs.barcodes.analysis import SampleSheetBarcodes
from auto_process_ngs.barcodes.analysis import Reporter
from auto_process_ngs.barcodes.analysis import report_barcodes

# BarcodeCounter
class TestBarcodeCounter(unittest.TestCase):
    def setUp(self):
        # Temporary working dir (if needed)
        self.wd = None

    def _make_working_dir(self):
        if self.wd is None:
            self.wd = tempfile.mkdtemp(suffix='.test_BarcodeCounter')

    def _make_file(self,name,contents):
        # Create a file under the working directory
        # called "name" and populated with "contents"
        # Working directory will be created if not already
        # set up
        # Returns the path to the file
        self._make_working_dir()
        filen = os.path.join(self.wd,name)
        with open(filen,'w') as fp:
            fp.write(contents)
        return filen

    def tearDown(self):
        # Remove temporary working dir
        if self.wd is not None and os.path.isdir(self.wd):
            shutil.rmtree(self.wd)

    def test_empty_counter(self):
        """BarcodeCounter: check empty counter
        """
        # Initialise counter object
        bc = BarcodeCounter()
        self.assertEqual(bc.barcodes(),[])
        self.assertEqual(bc.lanes,[])
        self.assertEqual(bc.filter_barcodes(),[])
        self.assertEqual(bc.counts("AGGCAGAATCTTACGC"),0)
        self.assertEqual(bc.counts("AGGCAGAATCTTACGC",lane=1),0)
        self.assertEqual(bc.counts_all("AGGCAGAATCTTACGC"),0)
        self.assertEqual(bc.nreads(),0)
        self.assertEqual(bc.nreads(1),0)

    def test_count_fastq_sequences(self):
        """BarcodeCounter: count barcode sequences
        """
        # Initialise counter object
        bc = BarcodeCounter()
        # Populate with sequences
        for r,incr in (((1,"AGGCAGAATCTTACGC"),102),
                       ((1,"TCCTGAGCTCTTACGC"),10),
                       ((1,"ACAGTGATTCTTTCCC"),3),
                       ((1,"ATGCTCGTCTCGCATC"),1),
                       ((2,"CGTACTAGTCTTACGC"),95),
                       ((2,"ATGTCAGATCTTTCCC"),29),
                       ((2,"AGGCAGAATCTTACGC"),12),
                       ((2,"CAGATCATTCTTTCCC"),6),
                       ((3,"GGACTCCTTCTTACGC"),75),
                       ((3,"ACCGATTCGCGCGTAG"),74),
                       ((3,"CCAGCAATATCGCGAG"),2),
                       ((3,"CCGCGTAAGCAATAGA"),1)):
            lane,seq = r
            for i in xrange(incr):
                bc.count_barcode(seq,lane=lane)
        # Check contents
        self.assertEqual(bc.barcodes(),["AGGCAGAATCTTACGC",
                                        "CGTACTAGTCTTACGC",
                                        "GGACTCCTTCTTACGC",
                                        "ACCGATTCGCGCGTAG",
                                        "ATGTCAGATCTTTCCC",
                                        "TCCTGAGCTCTTACGC",
                                        "CAGATCATTCTTTCCC",
                                        "ACAGTGATTCTTTCCC",
                                        "CCAGCAATATCGCGAG",
                                        "ATGCTCGTCTCGCATC",
                                        "CCGCGTAAGCAATAGA"])
        # Lanes
        self.assertEqual(bc.lanes,[1,2,3])
        # Counts for individual barcodes
        self.assertEqual(bc.counts("AGGCAGAATCTTACGC"),114)
        self.assertEqual(bc.counts("AGGCAGAATCTTACGC",lane=1),102)
        self.assertEqual(bc.counts("AGGCAGAATCTTACGC",lane=2),12)
        self.assertEqual(bc.counts("AGGCAGAATCTTACGC",lane=3),0)
        self.assertEqual(bc.counts_all("AGGCAGAATCTTACGC"),114)
        self.assertEqual(bc.counts("CCGCGTAAGCAATAGA"),1)
        self.assertEqual(bc.counts("CCGCGTAAGCAATAGA",lane=1),0)
        self.assertEqual(bc.counts("CCGCGTAAGCAATAGA",lane=2),0)
        self.assertEqual(bc.counts("CCGCGTAAGCAATAGA",lane=3),1)
        self.assertEqual(bc.counts_all("CCGCGTAAGCAATAGA"),1)
        # Read counts
        self.assertEqual(bc.nreads(),410)
        self.assertEqual(bc.nreads(1),116)
        self.assertEqual(bc.nreads(2),142)
        self.assertEqual(bc.nreads(3),152)
        # Lengths
        self.assertEqual(bc.barcode_lengths(),[16])
        self.assertEqual(bc.barcode_lengths(1),[16])
        self.assertEqual(bc.barcode_lengths(2),[16])
        self.assertEqual(bc.barcode_lengths(3),[16])

    def test_filter_barcodes(self):
        """BarcodeCounter: check filtering by lane and cutoff
        """
        bc = BarcodeCounter()
        bc.count_barcode("TATGCGCGGTA",lane=1,incr=285302)
        bc.count_barcode("TATGCGCGGTG",lane=1,incr=532)
        bc.count_barcode("ACCTACCGGTA",lane=1,incr=315)
        bc.count_barcode("CCCTTATGCGA",lane=1,incr=22)
	bc.count_barcode("ACCTAGCGGTA",lane=2,incr=477)
        bc.count_barcode("ACCTCTATGCT",lane=2,incr=368)
        self.assertEqual(bc.barcodes(),["TATGCGCGGTA",
                                        "TATGCGCGGTG",
                                        "ACCTAGCGGTA",
                                        "ACCTCTATGCT",
                                        "ACCTACCGGTA",
                                        "CCCTTATGCGA"])
        # No filtering
        self.assertEqual(bc.filter_barcodes(),["TATGCGCGGTA",
                                               "TATGCGCGGTG",
                                               "ACCTAGCGGTA",
                                               "ACCTCTATGCT",
                                               "ACCTACCGGTA",
                                               "CCCTTATGCGA"])
        # Filter by lane
        self.assertEqual(bc.filter_barcodes(lane=1),["TATGCGCGGTA",
                                                     "TATGCGCGGTG",
                                                     "ACCTACCGGTA",
                                                     "CCCTTATGCGA"]),
        self.assertEqual(bc.filter_barcodes(lane=2),["ACCTAGCGGTA",
                                                     "ACCTCTATGCT"])
        # Filter by cutoff
        self.assertEqual(bc.filter_barcodes(cutoff=0.5),
                         ["TATGCGCGGTA",])
        self.assertEqual(bc.filter_barcodes(cutoff=0.0015,lane=1),
                         ["TATGCGCGGTA","TATGCGCGGTG"])
        self.assertEqual(bc.filter_barcodes(cutoff=0.5,lane=2),
                         ["ACCTAGCGGTA",])

    def test_group(self):
        """BarcodeCounter: check grouping of barcode sequences
        """
        bc = BarcodeCounter()
        bc.count_barcode("TATGCGCGGTA",lane=1,incr=285302)
        bc.count_barcode("CATGCGCGGTA",lane=1,incr=8532)
        bc.count_barcode("GATGCGCGGTA",lane=1,incr=5321)
        bc.count_barcode("GCTGCGCGGTA",lane=1,incr=7853)
        bc.count_barcode("GCTGCGCGGTC",lane=1,incr=325394)
        bc.count_barcode("GTCACGCGGTA",lane=2,incr=296201)
        bc.count_barcode("GTCACGCGGTT",lane=2,incr=2853)
        bc.count_barcode("GTCACGCTGTT",lane=2,incr=278539)
        ## 2 mismatches across all lanes
        groups = bc.group(None,mismatches=2)
        ##"GCTGCGCGGTC","GCTGCGCGGTA","GATGCGCGGTA" = 338568
        ##"TATGCGCGGTA","CATGCGCGGTA" = 293834
        ##"GTCACGCGGTA","GTCACGCTGTT","GTCACGCGGTT" = 577593
        self.assertEqual(len(groups),3)
        self.assertEqual(groups[0].reference,"GTCACGCGGTA")
        self.assertEqual(groups[0].sequences,["GTCACGCGGTA",
                                              "GTCACGCTGTT",
                                              "GTCACGCGGTT"])
        self.assertEqual(groups[0].counts,577593)
        self.assertEqual(groups[1].reference,"GCTGCGCGGTC")
        self.assertEqual(groups[1].sequences,["GCTGCGCGGTC",
                                              "GCTGCGCGGTA",
                                              "GATGCGCGGTA"])
        self.assertEqual(groups[1].counts,338568)
        self.assertEqual(groups[2].reference,"TATGCGCGGTA")
        self.assertEqual(groups[2].sequences,["TATGCGCGGTA",
                                              "CATGCGCGGTA"])
        self.assertEqual(groups[2].counts,293834)
        ## 1 mismatch across all lanes
        groups = bc.group(None,mismatches=1)
        ##"TATGCGCGGTA","CATGCGCGGTA","GATGCGCGGTA" = 299155
        ##"GCTGCGCGGTC","GCTGCGCGGTA" = 333247
        ##"GTCACGCGGTA","GTCACGCGGTT" = 299054
        ##"GTCACGCTGTT" = 278539
        self.assertEqual(len(groups),4)
        self.assertEqual(groups[0].reference,"GCTGCGCGGTC")
        self.assertEqual(groups[0].sequences,["GCTGCGCGGTC",
                                              "GCTGCGCGGTA"])
        self.assertEqual(groups[0].counts,333247)
        self.assertEqual(groups[1].reference,"TATGCGCGGTA")
        self.assertEqual(groups[1].sequences,["TATGCGCGGTA",
                                              "CATGCGCGGTA",
                                              "GATGCGCGGTA"])
        self.assertEqual(groups[1].counts,299155)
        self.assertEqual(groups[2].reference,"GTCACGCGGTA")
        self.assertEqual(groups[2].sequences,["GTCACGCGGTA",
                                              "GTCACGCGGTT"])
        self.assertEqual(groups[2].counts,299054)
        self.assertEqual(groups[3].reference,"GTCACGCTGTT")
        self.assertEqual(groups[3].sequences,["GTCACGCTGTT",])
        self.assertEqual(groups[3].counts,278539)
        ## 1 mismatch in lane 1
        groups = bc.group(1,mismatches=1)
        ##"TATGCGCGGTA","CATGCGCGGTA","GATGCGCGGTA" = 299155
        ##"GCTGCGCGGTC","GCTGCGCGGTA" = 333247
        self.assertEqual(len(groups),2)
        self.assertEqual(groups[0].reference,"GCTGCGCGGTC")
        self.assertEqual(groups[0].sequences,["GCTGCGCGGTC",
                                              "GCTGCGCGGTA"])
        self.assertEqual(groups[0].counts,333247)
        self.assertEqual(groups[1].reference,"TATGCGCGGTA")
        self.assertEqual(groups[1].sequences,["TATGCGCGGTA",
                                              "CATGCGCGGTA",
                                              "GATGCGCGGTA"])
        self.assertEqual(groups[1].counts,299155)
        ## 2 mismatches across all lanes
        groups = bc.group(None,mismatches=2)
        ##"GCTGCGCGGTC","GCTGCGCGGTA","GATGCGCGGTA" = 338568
        ##"TATGCGCGGTA","CATGCGCGGTA" = 293834
        ##"GTCACGCGGTA","GTCACGCTGTT","GTCACGCGGTT" = 577593
        self.assertEqual(len(groups),3)
        self.assertEqual(groups[0].reference,"GTCACGCGGTA")
        self.assertEqual(groups[0].sequences,["GTCACGCGGTA",
                                              "GTCACGCTGTT",
                                              "GTCACGCGGTT"])
        self.assertEqual(groups[0].counts,577593)
        self.assertEqual(groups[1].reference,"GCTGCGCGGTC")
        self.assertEqual(groups[1].sequences,["GCTGCGCGGTC",
                                              "GCTGCGCGGTA",
                                              "GATGCGCGGTA"])
        self.assertEqual(groups[1].counts,338568)
        self.assertEqual(groups[2].reference,"TATGCGCGGTA")
        self.assertEqual(groups[2].sequences,["TATGCGCGGTA",
                                              "CATGCGCGGTA"])
        self.assertEqual(groups[2].counts,293834)

    def test_analyse(self):
        """BarcodeCounter: perform analysis with defaults
        """
        bc = BarcodeCounter()
        bc.count_barcode("TATGCGCGGTA",lane=1,incr=285302)
        bc.count_barcode("CATGCGCGGTA",lane=1,incr=8532)
        bc.count_barcode("GATGCGCGGTA",lane=1,incr=5321)
        bc.count_barcode("GCTGCGCGGTA",lane=1,incr=7853)
        bc.count_barcode("GCTGCGCGGTC",lane=1,incr=325394)
        analysis = bc.analyse(lane=1)
        self.assertEqual(analysis.cutoff,None)
        self.assertEqual(analysis.mismatches,0)
        self.assertEqual(analysis.total_reads,632402)
        self.assertEqual(analysis.coverage,632402)
        self.assertEqual(analysis.barcodes,["GCTGCGCGGTC",
                                            "TATGCGCGGTA",
                                            "CATGCGCGGTA",
                                            "GCTGCGCGGTA",
                                            "GATGCGCGGTA"])
        self.assertEqual(analysis.counts["GCTGCGCGGTC"].reads,325394)
        self.assertEqual(analysis.counts["TATGCGCGGTA"].reads,285302)
        self.assertEqual(analysis.counts["CATGCGCGGTA"].reads,8532)
        self.assertEqual(analysis.counts["GCTGCGCGGTA"].reads,7853)
        self.assertEqual(analysis.counts["GATGCGCGGTA"].reads,5321)
        self.assertEqual(analysis.counts["GCTGCGCGGTC"].sample,None)
        self.assertEqual(analysis.counts["TATGCGCGGTA"].sample,None)
        self.assertEqual(analysis.counts["CATGCGCGGTA"].sample,None)
        self.assertEqual(analysis.counts["GCTGCGCGGTA"].sample,None)
        self.assertEqual(analysis.counts["GATGCGCGGTA"].sample,None)
        self.assertEqual(analysis.counts["GCTGCGCGGTC"].sequences,1)
        self.assertEqual(analysis.counts["TATGCGCGGTA"].sequences,1)
        self.assertEqual(analysis.counts["CATGCGCGGTA"].sequences,1)
        self.assertEqual(analysis.counts["GCTGCGCGGTA"].sequences,1)
        self.assertEqual(analysis.counts["GATGCGCGGTA"].sequences,1)

    def test_analyse_with_cutoff(self):
        """BarcodeCounter: perform analysis with cutoff
        """
        bc = BarcodeCounter()
        bc.count_barcode("TATGCGCGGTA",lane=1,incr=285302)
        bc.count_barcode("CATGCGCGGTA",lane=1,incr=8532)
        bc.count_barcode("GATGCGCGGTA",lane=1,incr=5321)
        bc.count_barcode("GCTGCGCGGTA",lane=1,incr=7853)
        bc.count_barcode("GCTGCGCGGTC",lane=1,incr=325394)
        analysis = bc.analyse(lane=1,cutoff=0.013)
        self.assertEqual(analysis.cutoff,0.013)
        self.assertEqual(analysis.mismatches,0)
        self.assertEqual(analysis.total_reads,632402)
        self.assertEqual(analysis.coverage,619228)
        self.assertEqual(analysis.barcodes,["GCTGCGCGGTC",
                                            "TATGCGCGGTA",
                                            "CATGCGCGGTA"])
        self.assertEqual(analysis.counts["GCTGCGCGGTC"].reads,325394)
        self.assertEqual(analysis.counts["TATGCGCGGTA"].reads,285302)
        self.assertEqual(analysis.counts["CATGCGCGGTA"].reads,8532)
        self.assertEqual(analysis.counts["GCTGCGCGGTC"].sample,None)
        self.assertEqual(analysis.counts["TATGCGCGGTA"].sample,None)
        self.assertEqual(analysis.counts["CATGCGCGGTA"].sample,None)
        self.assertEqual(analysis.counts["GCTGCGCGGTC"].sequences,1)
        self.assertEqual(analysis.counts["TATGCGCGGTA"].sequences,1)
        self.assertEqual(analysis.counts["CATGCGCGGTA"].sequences,1)

    def test_analyse_with_sample_sheet(self):
        """BarcodeCounter: perform analysis with samplesheet
        """
        # Create sample sheet
        sample_sheet_file = self._make_file("SampleSheet.csv",
                                            """[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description
1,SMPL1,,,,A006,CATGCGCGGTA,,
1,SMPL2,,,,A012,GCTGCGCGGTC,,
2,SMPL3,,,,A005,ACAGTGCGGTA,,
2,SMPL4,,,,A019,GTGAAACGGTC,,
""")
        # Set up barcode counts
        bc = BarcodeCounter()
        bc.count_barcode("TATGCGCGGTA",lane=1,incr=285302)
        bc.count_barcode("CATGCGCGGTA",lane=1,incr=8532)
        bc.count_barcode("GATGCGCGGTA",lane=1,incr=5321)
        bc.count_barcode("GCTGCGCGGTA",lane=1,incr=7853)
        bc.count_barcode("GCTGCGCGGTC",lane=1,incr=325394)
        analysis = bc.analyse(lane=1,sample_sheet=sample_sheet_file)
        self.assertEqual(analysis.cutoff,None)
        self.assertEqual(analysis.mismatches,0)
        self.assertEqual(analysis.total_reads,632402)
        self.assertEqual(analysis.coverage,632402)
        self.assertEqual(analysis.barcodes,["GCTGCGCGGTC",
                                            "TATGCGCGGTA",
                                            "CATGCGCGGTA",
                                            "GCTGCGCGGTA",
                                            "GATGCGCGGTA"])
        self.assertEqual(analysis.counts["GCTGCGCGGTC"].reads,325394)
        self.assertEqual(analysis.counts["TATGCGCGGTA"].reads,285302)
        self.assertEqual(analysis.counts["CATGCGCGGTA"].reads,8532)
        self.assertEqual(analysis.counts["GCTGCGCGGTA"].reads,7853)
        self.assertEqual(analysis.counts["GATGCGCGGTA"].reads,5321)
        self.assertEqual(analysis.counts["GCTGCGCGGTC"].sample,"SMPL2")
        self.assertEqual(analysis.counts["TATGCGCGGTA"].sample,None)
        self.assertEqual(analysis.counts["CATGCGCGGTA"].sample,"SMPL1")
        self.assertEqual(analysis.counts["GCTGCGCGGTA"].sample,None)
        self.assertEqual(analysis.counts["GATGCGCGGTA"].sample,None)
        self.assertEqual(analysis.counts["GCTGCGCGGTC"].sequences,1)
        self.assertEqual(analysis.counts["TATGCGCGGTA"].sequences,1)
        self.assertEqual(analysis.counts["CATGCGCGGTA"].sequences,1)
        self.assertEqual(analysis.counts["GCTGCGCGGTA"].sequences,1)
        self.assertEqual(analysis.counts["GATGCGCGGTA"].sequences,1)

    def test_analyse_groups(self):
        """BarcodeCounter: perform analysis with grouping
        """
        bc = BarcodeCounter()
        bc.count_barcode("TATGCGCGGTA",lane=1,incr=285302)
        bc.count_barcode("CATGCGCGGTA",lane=1,incr=8532)
        bc.count_barcode("GATGCGCGGTA",lane=1,incr=5321)
        bc.count_barcode("GCTGCGCGGTA",lane=1,incr=7853)
        bc.count_barcode("GCTGCGCGGTC",lane=1,incr=325394)
        analysis = bc.analyse(lane=1,mismatches=1)
        ##"TATGCGCGGTA","CATGCGCGGTA","GATGCGCGGTA" = 299155
        ##"GCTGCGCGGTC","GCTGCGCGGTA" = 333247
        self.assertEqual(analysis.cutoff,None)
        self.assertEqual(analysis.mismatches,1)
        self.assertEqual(analysis.total_reads,632402)
        self.assertEqual(analysis.coverage,632402)
        self.assertEqual(analysis.barcodes,["GCTGCGCGGTC",
                                            "TATGCGCGGTA"])
        self.assertEqual(analysis.counts["GCTGCGCGGTC"].reads,333247)
        self.assertEqual(analysis.counts["TATGCGCGGTA"].reads,299155)
        self.assertEqual(analysis.counts["GCTGCGCGGTC"].sample,None)
        self.assertEqual(analysis.counts["TATGCGCGGTA"].sample,None)
        self.assertEqual(analysis.counts["GCTGCGCGGTC"].sequences,2)
        self.assertEqual(analysis.counts["TATGCGCGGTA"].sequences,3)

    def test_analyse_groups_with_sample_sheet(self):
        """BarcodeCounter: perform analysis with grouping and samplesheet
        """
        # Create sample sheet
        sample_sheet_file = self._make_file("SampleSheet.csv",
                                            """[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description
1,SMPL1,,,,A006,CATGCGCGGTA,,
1,SMPL2,,,,A012,GCTGCGCGGTC,,
2,SMPL3,,,,A005,ACAGTGCGGTA,,
2,SMPL4,,,,A019,GTGAAACGGTC,,
""")
        # Set up barcode counts
        bc = BarcodeCounter()
        bc.count_barcode("TATGCGCGGTA",lane=1,incr=285302)
        bc.count_barcode("CATGCGCGGTA",lane=1,incr=8532)
        bc.count_barcode("GATGCGCGGTA",lane=1,incr=5321)
        bc.count_barcode("GCTGCGCGGTA",lane=1,incr=7853)
        bc.count_barcode("GCTGCGCGGTC",lane=1,incr=325394)
        analysis = bc.analyse(lane=1,
                              mismatches=2,
                              sample_sheet=sample_sheet_file)
        ##"CATGCGCGGTA","TATGCGCGGTA","GATGCGCGGTA","GCTGCGCGGTA" = 307008
        ##"GCTGCGCGGTC" = 325394
        self.assertEqual(analysis.cutoff,None)
        self.assertEqual(analysis.mismatches,2)
        self.assertEqual(analysis.total_reads,632402)
        self.assertEqual(analysis.coverage,632402)
        self.assertEqual(analysis.barcodes,["GCTGCGCGGTC",
                                            "CATGCGCGGTA"])
        self.assertEqual(analysis.counts["GCTGCGCGGTC"].reads,325394)
        self.assertEqual(analysis.counts["CATGCGCGGTA"].reads,307008)
        self.assertEqual(analysis.counts["GCTGCGCGGTC"].sample,"SMPL2")
        self.assertEqual(analysis.counts["CATGCGCGGTA"].sample,"SMPL1")
        self.assertEqual(analysis.counts["GCTGCGCGGTC"].sequences,1)
        self.assertEqual(analysis.counts["CATGCGCGGTA"].sequences,4)

    def test_analyse_with_no_counts(self):
        """BarcodeCounter: perform analysis for zero counts
        """
        bc = BarcodeCounter()
        analysis = bc.analyse()
        self.assertEqual(analysis.cutoff,None)
        self.assertEqual(analysis.mismatches,0)
        self.assertEqual(analysis.total_reads,0)
        self.assertEqual(analysis.coverage,0)
        self.assertEqual(analysis.barcodes,[])

    def test_read_counts_file(self):
        """BarcodeCounter: read in data from '.counts' file
        """
        # Read a counts file
        counts_file = self._make_file("test.counts","""#Lane	Rank	Sequence	Count
1	1	TATGCGCGGTA	285302
1	2	TATGCGCGGTG	532
1	3	ACCTACCGGTA	315
1	4	CCCTTATGCGA	22
2	5	ACCTAGCGGTA	477
2	6	ACCTCTATGCT	368
3	7	ACCCTNCGGTA	312
3	8	ACCTTATGCGC	248""")
        # Read the file
        bc = BarcodeCounter(counts_file)
        # Check the contents
        self.assertEqual(bc.barcodes(),["TATGCGCGGTA",
                                        "TATGCGCGGTG",
                                        "ACCTAGCGGTA",
                                        "ACCTCTATGCT",
                                        "ACCTACCGGTA",
                                        "ACCCTNCGGTA",
                                        "ACCTTATGCGC",
                                        "CCCTTATGCGA"])
        # Lanes
        self.assertEqual(bc.lanes,[1,2,3])
        # Counts for individual barcodes
        self.assertEqual(bc.counts("TATGCGCGGTA"),285302)
        self.assertEqual(bc.counts("TATGCGCGGTG"),532)
        self.assertEqual(bc.counts("ACCTAGCGGTA"),477)
        self.assertEqual(bc.counts("ACCTCTATGCT"),368)
        self.assertEqual(bc.counts("ACCTACCGGTA"),315)
        self.assertEqual(bc.counts("ACCCTNCGGTA"),312)
        self.assertEqual(bc.counts("ACCTTATGCGC"),248)
        self.assertEqual(bc.counts("CCCTTATGCGA"),22)
        self.assertEqual(bc.counts("TATGCGCGGTA",lane=1),285302)
        self.assertEqual(bc.counts("TATGCGCGGTA",lane=2),0)
        self.assertEqual(bc.counts("TATGCGCGGTA",lane=3),0)
        self.assertEqual(bc.counts_all("TATGCGCGGTA"),285302)
        self.assertEqual(bc.counts("ACCTTATGCGC",lane=1),0)
        self.assertEqual(bc.counts("ACCTTATGCGC",lane=2),0)
        self.assertEqual(bc.counts("ACCTTATGCGC",lane=3),248)
        self.assertEqual(bc.counts_all("ACCTTATGCGC"),248)
        # Read counts
        self.assertEqual(bc.nreads(),287576)
        self.assertEqual(bc.nreads(1),286171)
        self.assertEqual(bc.nreads(2),845)
        self.assertEqual(bc.nreads(3),560)

    def test_read_multiple_counts_file(self):
        """BarcodeCounter: read in data from multiple '.counts' files
        """
        # Read multiple counts files
        counts_lane1 = self._make_file("lane1.counts",
                                       """#Lane	Rank	Sequence	Count
1	1	TATGCGCGGTA	285302
1	2	TATGCGCGGTG	532
1	3	ACCTACCGGTA	315
1	4	CCCTTATGCGA	22""")
        counts_lane2 = self._make_file("lane2.counts",
                                       """#Lane	Rank	Sequence	Count
2	1	ACCTAGCGGTA	477
2	2	ACCTCTATGCT	368""")
        counts_lane3 = self._make_file("lane3.counts",
                                       """#Lane	Rank	Sequence	Count
3	1	ACCCTNCGGTA	312
3	2	ACCTTATGCGC	248""")
        # Read the file
        bc = BarcodeCounter(counts_lane1,counts_lane2,counts_lane3)
        # Check the contents
        self.assertEqual(bc.barcodes(),["TATGCGCGGTA",
                                        "TATGCGCGGTG",
                                        "ACCTAGCGGTA",
                                        "ACCTCTATGCT",
                                        "ACCTACCGGTA",
                                        "ACCCTNCGGTA",
                                        "ACCTTATGCGC",
                                        "CCCTTATGCGA"])
        # Lanes
        self.assertEqual(bc.lanes,[1,2,3])
        # Counts for individual barcodes
        self.assertEqual(bc.counts("TATGCGCGGTA"),285302)
        self.assertEqual(bc.counts("TATGCGCGGTG"),532)
        self.assertEqual(bc.counts("ACCTAGCGGTA"),477)
        self.assertEqual(bc.counts("ACCTCTATGCT"),368)
        self.assertEqual(bc.counts("ACCTACCGGTA"),315)
        self.assertEqual(bc.counts("ACCCTNCGGTA"),312)
        self.assertEqual(bc.counts("ACCTTATGCGC"),248)
        self.assertEqual(bc.counts("CCCTTATGCGA"),22)
        self.assertEqual(bc.counts("TATGCGCGGTA",lane=1),285302)
        self.assertEqual(bc.counts("TATGCGCGGTA",lane=2),0)
        self.assertEqual(bc.counts("TATGCGCGGTA",lane=3),0)
        self.assertEqual(bc.counts_all("TATGCGCGGTA"),285302)
        self.assertEqual(bc.counts("ACCTTATGCGC",lane=1),0)
        self.assertEqual(bc.counts("ACCTTATGCGC",lane=2),0)
        self.assertEqual(bc.counts("ACCTTATGCGC",lane=3),248)
        self.assertEqual(bc.counts_all("ACCTTATGCGC"),248)
        # Read counts
        self.assertEqual(bc.nreads(),287576)
        self.assertEqual(bc.nreads(1),286171)
        self.assertEqual(bc.nreads(2),845)
        self.assertEqual(bc.nreads(3),560)

    def test_read_old_style_counts_file(self):
        """BarcodeCounter: read in data from old-style 3 column '.counts' file
        """
        # Read old-style 3 column counts files
        self._make_working_dir()
        old_style_counts_file = self._make_file("old_style.counts",
                                                """#Rank	Sequence	Count
1	TATGCGCGGTA	285302
2	TATGCGCGGTG	532
3	ACCTACCGGTA	315
4	CCCTTATGCGA	22""")
        # Read the file
        bc = BarcodeCounter(old_style_counts_file)
        # Check the contents
        self.assertEqual(bc.barcodes(),["TATGCGCGGTA",
                                        "TATGCGCGGTG",
                                        "ACCTACCGGTA",
                                        "CCCTTATGCGA"])
        # Lanes
        self.assertEqual(bc.lanes,[])
        # Counts for individual barcodes
        self.assertEqual(bc.counts("TATGCGCGGTA"),285302)
        self.assertEqual(bc.counts("TATGCGCGGTG"),532)
        self.assertEqual(bc.counts("ACCTACCGGTA"),315)
        self.assertEqual(bc.counts("CCCTTATGCGA"),22)
        self.assertEqual(bc.counts("TATGCGCGGTA",lane=1),0)
        self.assertEqual(bc.counts_all("TATGCGCGGTA"),285302)
        # Read counts
        self.assertEqual(bc.nreads(),286171)

    def test_write_counts_file(self):
        """BarcodeCounter: write counts to a file
        """
        # Write a file
        self._make_working_dir()
        bc = BarcodeCounter()
        bc.count_barcode("TATGCGCGGTA",lane=1,incr=285302)
        bc.count_barcode("TATGCGCGGTG",lane=1,incr=532)
        bc.count_barcode("ACCTACCGGTA",lane=1,incr=315)
        bc.count_barcode("CCCTTATGCGA",lane=1,incr=22)
        bc.count_barcode("ACCTAGCGGTA",lane=2,incr=477)
        bc.count_barcode("ACCTCTATGCT",lane=2,incr=368)
        bc.count_barcode("ACCCTNCGGTA",lane=3,incr=312)
        bc.count_barcode("ACCTTATGCGC",lane=3,incr=248)
        counts_file = os.path.join(self.wd,"out.counts")
        bc.write(counts_file)
        expected_contents = """#Lane	Rank	Sequence	Count
1	1	TATGCGCGGTA	285302
1	2	TATGCGCGGTG	532
1	3	ACCTACCGGTA	315
1	4	CCCTTATGCGA	22
2	1	ACCTAGCGGTA	477
2	2	ACCTCTATGCT	368
3	1	ACCCTNCGGTA	312
3	2	ACCTTATGCGC	248
"""
        self.assertTrue(os.path.exists(counts_file))
        self.assertEqual(open(counts_file,'r').read(),
                         expected_contents)

# BarcodeGroup
class TestBarcodeGroup(unittest.TestCase):
    def test_barcodegroup(self):
        """BarcodeGroup: check making a new instance
        """
        # Create a new BarcodeGroup
        grp = BarcodeGroup("CTAAGCCT",2894178)
        self.assertEqual(grp.reference,"CTAAGCCT")
        self.assertEqual(grp.sequences,["CTAAGCCT"])
        self.assertEqual(grp.counts,2894178)
        self.assertEqual(len(grp),1)

    def test_barcodegroup_add(self):
        """BarcodeGroup: check adding a sequence to the group
        """
        # Add sequences to a BarcodeGroup
        grp = BarcodeGroup("CTAAGCCT",2894178)
        grp.add("CTAAGCCA",92417)
        self.assertEqual(grp.reference,"CTAAGCCT")
        self.assertEqual(grp.sequences,["CTAAGCCT","CTAAGCCA"])
        self.assertEqual(grp.counts,2986595)
        self.assertEqual(len(grp),2)

    def test_barcodegroup_match(self):
        """BarcodeGroup: check matching sequences against the reference
        """
        # Match sequence against the group
        grp = BarcodeGroup("CTAAGCCT",2894178)
        # 1 mismatch allowed
        self.assertTrue(grp.match("CTAAGCCA",mismatches=1))
        self.assertFalse(grp.match("CGAAGCCA",mismatches=1))
        self.assertFalse(grp.match("CGATGCCA",mismatches=1))
        # Default (2 mismatches)
        self.assertTrue(grp.match("CTAAGCCA"))
        self.assertTrue(grp.match("CGAAGCCA"))
        self.assertFalse(grp.match("CGATGCCA"))
        # Differing lengths of barcode
        # -- Too short
        self.assertFalse(grp.match("CGATGCC"))
        # -- Too long
        self.assertFalse(grp.match("CGATGCCGG"))

# SampleSheetBarcodes
class TestSampleSheetBarcodes(unittest.TestCase):
    def setUp(self):
        # Test data
        sample_sheet_header = "[Header]\nIEMFileVersion,4\n\n[Reads]\n150\n150\n\n[Settings]\nAdapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA\nAdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT\n\n"
        sample_sheet_single_index_with_lanes = """
"""
        # Temporary working dir
        self.wd = tempfile.mkdtemp(suffix='.test_SampleSheetBarcodes')
        # Create files
        sample_sheet_single_index_no_lanes = """[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description
ES1,,,,A006,GCCAAT,,
EP1,,,,A012,CTTGTA,,
ES2,,,,A005,ACAGTG,,
EP2,,,,A019,GTGAAA,,"""
        self.single_index_no_lanes = \
            os.path.join(self.wd,"single_index_no_lanes.csv")
        with open(self.single_index_no_lanes,"w") as fp:
            fp.write(sample_sheet_header +
                     sample_sheet_single_index_no_lanes)
        #
        sample_sheet_dual_index_no_lanes = """[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
SW1,SW1,,,N701,TAAGGCGA,S517,TCTTACGC,,
SW2,SW2,,,N702,CGTACTAG,S517,TCTTACGC,,
SW3,SW3,,,N703,AGGCAGAA,S517,TCTTACGC,,
SW4,SW4,,,N704,TCCTGAGC,S517,TCTTACGC,,
SW5,SW5,,,N705,GGACTCCT,S517,TCTTACGC,,
SW6,SW6,,,N706,TAGGCATG,S517,TCTTACGC,,"""
        self.dual_index_no_lanes = \
            os.path.join(self.wd,"dual_index_no_lanes.csv")
        with open(self.dual_index_no_lanes,"w") as fp:
            fp.write(sample_sheet_header +
                     sample_sheet_dual_index_no_lanes)
        #
        sample_sheet_single_index_with_lanes = """[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description
1,ES1,,,,A006,GCCAAT,,
1,EP1,,,,A012,CTTGTA,,
2,ES2,,,,A005,ACAGTG,,
2,EP2,,,,A019,GTGAAA,,"""
        self.single_index_with_lanes = \
            os.path.join(self.wd,"single_index_with_lanes.csv")
        with open(self.single_index_with_lanes,"w") as fp:
            fp.write(sample_sheet_header +
                     sample_sheet_single_index_with_lanes)
        #
        sample_sheet_dual_index_with_lanes = """[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,DI1,DI1,,,D701,CGTGTAGG,D501,GACCTGTA,HO,
1,DI2,DI2,,,D702,CGTGTAGG,D501,ATGTAACT,HO,
1,DI3,DI3,,,D703,CGTGTAGG,D501,GTTTCAGA,HO,
1,DI4,DI4,,,D704,CGTGTAGG,D501,CACAGGAT,HO,
2,E11,E11,,,D703,GCCAATAT,D503,TCTTTCCC,FL,
2,E12,E12,,,D704,CAGATCAT,D504,TCTTTCCC,FL,
3,AT1,AT1,,,D701,GGATTCGC,D501,TAGTAGCC,JF,
3,AT3,AT3,,,D702,TACCAGCG,D502,CGCTGCTG,JF,
3,AT4,AT4,,,D703,ACCGATTC,D503,GCGCGTAG,JF,
3,AT5,AT5,,,D704,CCGCGTAA,D504,GCAATAGA,JF,
3,AEx,AEx,,,D705,ATGCTCGT,D505,CTCGCATC,JF,
4,AEx,AEx,,,D705,ATGCTCGT,D505,CTCGCATC,JF,
4,AD5,AD5,,,D707,CCAGCAAT,D507,ATCGCGAG,JF,
4,D3K1,D3K1,,,D701,ACAGTGAT,D501,TCTTTCCC,FL,
4,D3K2,D3K2,,,D702,ATGTCAGA,D502,TCTTTCCC,FL,"""
        self.dual_index_with_lanes = \
            os.path.join(self.wd,"dual_index_with_lanes.csv")
        with open(self.dual_index_with_lanes,"w") as fp:
            fp.write(sample_sheet_header +
                     sample_sheet_dual_index_with_lanes)

    def tearDown(self):
        # Remove temporary working dir
        if os.path.isdir(self.wd):
            shutil.rmtree(self.wd)

    def test_single_index_no_lanes(self):
        """SampleSheetBarcodes: single index sample sheet with no lanes defined
        """
        s = SampleSheetBarcodes(self.single_index_no_lanes)
        # Check all barcodes
        self.assertEqual(s.barcodes(),["ACAGTG","CTTGTA",
                                       "GCCAAT","GTGAAA"])
        # Check all samples
        self.assertEqual(s.samples(),["EP1","EP2","ES1","ES2"])
        # Look up sample names for barcodes in each lane
        self.assertEqual(s.lookup_sample("GCCAAT"),"ES1")
        self.assertEqual(s.lookup_sample("CTTGTA"),"EP1")
        self.assertEqual(s.lookup_sample("ACAGTG"),"ES2")
        self.assertEqual(s.lookup_sample("GTGAAA"),"EP2")
        # Look up barcode matching sample names in each lane
        self.assertEqual(s.lookup_barcode("ES1"),"GCCAAT")
        self.assertEqual(s.lookup_barcode("EP1"),"CTTGTA")
        self.assertEqual(s.lookup_barcode("ES2"),"ACAGTG")
        self.assertEqual(s.lookup_barcode("EP2"),"GTGAAA")

    def test_dual_index_no_lanes(self):
        """SampleSheetBarcodes: dual index sample sheet with no lanes defined
        """
        s = SampleSheetBarcodes(self.dual_index_no_lanes)
        # Check all barcodes
        self.assertEqual(s.barcodes(),["AGGCAGAATCTTACGC",
                                       "CGTACTAGTCTTACGC",
                                       "GGACTCCTTCTTACGC",
                                       "TAAGGCGATCTTACGC",
                                       "TAGGCATGTCTTACGC",
                                       "TCCTGAGCTCTTACGC"])
        # Check all samples
        self.assertEqual(s.samples(),["SW1","SW2","SW3","SW4","SW5","SW6"])
        # Look up sample names for barcodes in each lane
        self.assertEqual(s.lookup_sample("TAAGGCGATCTTACGC"),"SW1")
        self.assertEqual(s.lookup_sample("CGTACTAGTCTTACGC"),"SW2")
        self.assertEqual(s.lookup_sample("AGGCAGAATCTTACGC"),"SW3")
        self.assertEqual(s.lookup_sample("TCCTGAGCTCTTACGC"),"SW4")
        self.assertEqual(s.lookup_sample("GGACTCCTTCTTACGC"),"SW5")
        self.assertEqual(s.lookup_sample("TAGGCATGTCTTACGC"),"SW6")
        # Look up barcode matching sample names in each lane
        self.assertEqual(s.lookup_barcode("SW1"),"TAAGGCGATCTTACGC")
        self.assertEqual(s.lookup_barcode("SW2"),"CGTACTAGTCTTACGC")
        self.assertEqual(s.lookup_barcode("SW3"),"AGGCAGAATCTTACGC")
        self.assertEqual(s.lookup_barcode("SW4"),"TCCTGAGCTCTTACGC")
        self.assertEqual(s.lookup_barcode("SW5"),"GGACTCCTTCTTACGC")
        self.assertEqual(s.lookup_barcode("SW6"),"TAGGCATGTCTTACGC")

    def test_single_index_with_lanes(self):
        """SampleSheetBarcodes: single index sample sheet with lanes defined
        """
        s = SampleSheetBarcodes(self.single_index_with_lanes)
        # Check all barcodes
        self.assertEqual(s.barcodes(),["ACAGTG","CTTGTA",
                                       "GCCAAT","GTGAAA"])
        # Check barcodes in each lane
        self.assertEqual(s.barcodes(1),["CTTGTA","GCCAAT"])
        self.assertEqual(s.barcodes(2),["ACAGTG","GTGAAA"])
        # Check all samples
        self.assertEqual(s.samples(),["EP1","EP2","ES1","ES2"])
        # Check samples in each lane
        self.assertEqual(s.samples(1),["EP1","ES1"])
        self.assertEqual(s.samples(2),["EP2","ES2"])
        # Look up sample names for barcodes in each lane
        self.assertEqual(s.lookup_sample("GCCAAT",1),"ES1")
        self.assertEqual(s.lookup_sample("CTTGTA",1),"EP1")
        self.assertEqual(s.lookup_sample("ACAGTG",2),"ES2")
        self.assertEqual(s.lookup_sample("GTGAAA",2),"EP2")
        # Look up barcode matching sample names in each lane
        self.assertEqual(s.lookup_barcode("ES1",1),"GCCAAT")
        self.assertEqual(s.lookup_barcode("EP1",1),"CTTGTA")
        self.assertEqual(s.lookup_barcode("ES2",2),"ACAGTG")
        self.assertEqual(s.lookup_barcode("EP2",2),"GTGAAA")

    def test_dual_index_with_lanes(self):
        """SampleSheetBarcodes: dual index sample sheet with lanes defined
        """
        s = SampleSheetBarcodes(self.dual_index_with_lanes)
        # Check all barcodes
        self.assertEqual(s.barcodes(),["ACAGTGATTCTTTCCC",
                                       "ACCGATTCGCGCGTAG",
                                       "ATGCTCGTCTCGCATC",
                                       "ATGCTCGTCTCGCATC",
                                       "ATGTCAGATCTTTCCC",
                                       "CAGATCATTCTTTCCC",
                                       "CCAGCAATATCGCGAG",
                                       "CCGCGTAAGCAATAGA",
                                       "CGTGTAGGATGTAACT",
                                       "CGTGTAGGCACAGGAT",
                                       "CGTGTAGGGACCTGTA",
                                       "CGTGTAGGGTTTCAGA",
                                       "GCCAATATTCTTTCCC",
                                       "GGATTCGCTAGTAGCC",
                                       "TACCAGCGCGCTGCTG"])
        # Check barcodes in each lane
        self.assertEqual(s.barcodes(1),["CGTGTAGGATGTAACT",
                                        "CGTGTAGGCACAGGAT",
                                        "CGTGTAGGGACCTGTA",
                                        "CGTGTAGGGTTTCAGA"])
        self.assertEqual(s.barcodes(2),["CAGATCATTCTTTCCC",
                                        "GCCAATATTCTTTCCC"])
        self.assertEqual(s.barcodes(3),["ACCGATTCGCGCGTAG",
                                        "ATGCTCGTCTCGCATC",
                                        "CCGCGTAAGCAATAGA",
                                        "GGATTCGCTAGTAGCC",
                                        "TACCAGCGCGCTGCTG"])
        self.assertEqual(s.barcodes(4),["ACAGTGATTCTTTCCC",
                                        "ATGCTCGTCTCGCATC",
                                        "ATGTCAGATCTTTCCC",
                                        "CCAGCAATATCGCGAG"])
        # Check all samples
        self.assertEqual(s.samples(),["AD5","AEx","AEx",
                                      "AT1","AT3","AT4","AT5",
                                      "D3K1","D3K2",
                                      "DI1","DI2","DI3","DI4",
                                      "E11","E12"])
        # Check samples in each lane
        self.assertEqual(s.samples(1),["DI1","DI2","DI3","DI4"])
        self.assertEqual(s.samples(2),["E11","E12"])
        self.assertEqual(s.samples(3),["AEx","AT1","AT3","AT4","AT5"])
        self.assertEqual(s.samples(4),["AD5","AEx","D3K1","D3K2"])
        # Look up sample names for barcodes in each lane
        self.assertEqual(s.lookup_sample("CGTGTAGGGACCTGTA",1),"DI1")
        self.assertEqual(s.lookup_sample("CGTGTAGGATGTAACT",1),"DI2")
        self.assertEqual(s.lookup_sample("CGTGTAGGGTTTCAGA",1),"DI3")
        self.assertEqual(s.lookup_sample("CGTGTAGGCACAGGAT",1),"DI4")
        self.assertEqual(s.lookup_sample("GCCAATATTCTTTCCC",2),"E11")
        self.assertEqual(s.lookup_sample("CAGATCATTCTTTCCC",2),"E12")
        self.assertEqual(s.lookup_sample("GGATTCGCTAGTAGCC",3),"AT1")
        self.assertEqual(s.lookup_sample("TACCAGCGCGCTGCTG",3),"AT3")
        self.assertEqual(s.lookup_sample("ACCGATTCGCGCGTAG",3),"AT4")
        self.assertEqual(s.lookup_sample("CCGCGTAAGCAATAGA",3),"AT5")
        self.assertEqual(s.lookup_sample("ATGCTCGTCTCGCATC",3),"AEx")
        self.assertEqual(s.lookup_sample("ATGCTCGTCTCGCATC",4),"AEx")
        self.assertEqual(s.lookup_sample("CCAGCAATATCGCGAG",4),"AD5")
        self.assertEqual(s.lookup_sample("ACAGTGATTCTTTCCC",4),"D3K1")
        self.assertEqual(s.lookup_sample("ATGTCAGATCTTTCCC",4),"D3K2")
        # Look up barcode matching sample names in each lane
        self.assertEqual(s.lookup_barcode("DI1",1),"CGTGTAGGGACCTGTA")
        self.assertEqual(s.lookup_barcode("DI2",1),"CGTGTAGGATGTAACT")
        self.assertEqual(s.lookup_barcode("DI3",1),"CGTGTAGGGTTTCAGA")
        self.assertEqual(s.lookup_barcode("DI4",1),"CGTGTAGGCACAGGAT")
        self.assertEqual(s.lookup_barcode("E11",2),"GCCAATATTCTTTCCC")
        self.assertEqual(s.lookup_barcode("E12",2),"CAGATCATTCTTTCCC")
        self.assertEqual(s.lookup_barcode("AT1",3),"GGATTCGCTAGTAGCC")
        self.assertEqual(s.lookup_barcode("AT3",3),"TACCAGCGCGCTGCTG")
        self.assertEqual(s.lookup_barcode("AT4",3),"ACCGATTCGCGCGTAG")
        self.assertEqual(s.lookup_barcode("AT5",3),"CCGCGTAAGCAATAGA")
        self.assertEqual(s.lookup_barcode("AEx",3),"ATGCTCGTCTCGCATC")
        self.assertEqual(s.lookup_barcode("AEx",4),"ATGCTCGTCTCGCATC")
        self.assertEqual(s.lookup_barcode("AD5",4),"CCAGCAATATCGCGAG")
        self.assertEqual(s.lookup_barcode("D3K1",4),"ACAGTGATTCTTTCCC")
        self.assertEqual(s.lookup_barcode("D3K2",4),"ATGTCAGATCTTTCCC")

    def test_request_non_existent_lane(self):
        """SampleSheetBarcodes: handle request for non-existent lane
        """
        # Bad lane for sample sheet with lanes
        s = SampleSheetBarcodes(self.dual_index_with_lanes)
        self.assertRaises(KeyError,s.barcodes,5)
        self.assertRaises(KeyError,s.samples,5)
        # Any lane for sample sheet with no lanes
        s = SampleSheetBarcodes(self.dual_index_no_lanes)
        self.assertEqual(s.barcodes(1),
                         ["AGGCAGAATCTTACGC",
                          "CGTACTAGTCTTACGC",
                          "GGACTCCTTCTTACGC",
                          "TAAGGCGATCTTACGC",
                          "TAGGCATGTCTTACGC",
                          "TCCTGAGCTCTTACGC"])
        self.assertEqual(s.samples(1),
                         ["SW1","SW2","SW3","SW4","SW5","SW6"])

# Reporter
class TestReporter(unittest.TestCase):
    def setUp(self):
        # Temporary working dir (if needed)
        self.wd = None

    def _make_working_dir(self):
        if self.wd is None:
            self.wd = tempfile.mkdtemp(suffix='.test_Reporter')

    def tearDown(self):
        # Remove temporary working dir
        if self.wd is not None and os.path.isdir(self.wd):
            shutil.rmtree(self.wd)

    def test_add(self):
        """Reporter: can add content
        """
        reporter = Reporter()
        self.assertEqual(len(reporter),0)
        self.assertFalse(reporter)
        self.assertEqual(str(reporter),"")
        reporter.add("Title text",title=True)
        reporter.add("Some words")
        self.assertEqual(len(reporter),2)
        self.assertTrue(reporter)
        self.assertEqual(str(reporter),
                         """Title text
==========
Some words""")

    def test_write(self):
        """Reporter: can write to a text file
        """
        self._make_working_dir()
        reporter = Reporter()
        reporter.add("Test Document",title=True)
        reporter.add("Lorem ipsum")
        report_file = os.path.join(self.wd,"report.txt")
        reporter.write(filen=report_file,title="My Report")
        self.assertTrue(os.path.isfile(report_file))
        expected_contents = """My Report
*********

Test Document
=============
Lorem ipsum
"""
        self.assertTrue(os.path.exists(report_file))
        self.assertEqual(open(report_file,'r').read(),
                         expected_contents)

    def test_write_xls(self):
        """Reporter: can write to an XLS file
        """
        self._make_working_dir()
        reporter = Reporter()
        reporter.add("Test Document",title=True)
        reporter.add("Lorem ipsum")
        report_xls = os.path.join(self.wd,"report.xls")
        reporter.write_xls(report_xls)
        self.assertTrue(os.path.isfile(report_xls))

    def test_write_html(self):
        """Reporter: can write to a HTML file
        """
        self._make_working_dir()
        reporter = Reporter()
        reporter.add("This is a Test",title=True)
        reporter.add("Lorem ipsum")
        reporter.add("Column1\tColumn2",heading=True)
        reporter.add("1\t2")
        reporter.add("3\t4")
        reporter.add("Lorem more ipsum")
        report_html = os.path.join(self.wd,"report.html")
        reporter.write_html(report_html,
                            title="Test Document",
                            no_styles=True)
        self.assertTrue(os.path.isfile(report_html))
        expected_contents = """<html>
<head>
<title>Test Document</title>
</head>
<body>
<h1>Test Document</h1>
<div id='toc'>
<h2>Contents</h2>
<ul><li><a href='#This_is_a_Test'>This is a Test</a></li></ul>
</div>
<div id='This_is_a_Test'>
<h2>This is a Test</h2>
<p>Lorem ipsum</p>
</div>
<div>
<table>
<tr><th>Column1</th><th>Column2</th></tr>
<tr><td>1</td><td>2</td></tr>
<tr><td>3</td><td>4</td></tr>
</table>
</div>
<div>
<p>Lorem more ipsum</p>
</div></body>
</html>
"""
        for expected,actual in zip(expected_contents.split('\n'),
                                   open(report_html,'r').read().split('\n')):
            self.assertEqual(expected,actual)

# report_barcodes
class TestReportBarcodesFunction(unittest.TestCase):
    def setUp(self):
        # Temporary working dir (if needed)
        self.wd = None

    def _make_working_dir(self):
        if self.wd is None:
            self.wd = tempfile.mkdtemp(suffix='.test_report_barcodes')

    def tearDown(self):
        # Remove temporary working dir
        if self.wd is not None and os.path.isdir(self.wd):
            shutil.rmtree(self.wd)

    def _make_file(self,name,contents):
        # Create a file under the working directory
        # called "name" and populated with "contents"
        # Working directory will be created if not already
        # set up
        # Returns the path to the file
        self._make_working_dir()
        filen = os.path.join(self.wd,name)
        with open(filen,'w') as fp:
            fp.write(contents)
        return filen

    def test_report_barcodes(self):
        """report_barcodes: check output for mismatches and sample sheet
        """
        # Create sample sheet
        sample_sheet_file = self._make_file("SampleSheet.csv",
                                            """[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description
1,SMPL1,,,,A006,CATGCGCGGTA,,
1,SMPL2,,,,A012,GCTGCGCGGTC,,
""")
        # Set up barcode counts
        bc = BarcodeCounter()
        bc.count_barcode("TATGCGCGGTA",lane=1,incr=285302)
        bc.count_barcode("CATGCGCGGTA",lane=1,incr=8532)
        bc.count_barcode("GATGCGCGGTA",lane=1,incr=5321)
        bc.count_barcode("GCTGCGCGGTA",lane=1,incr=7853)
        bc.count_barcode("GCTGCGCGGTC",lane=1,incr=325394)
        analysis = bc.analyse(lane=1,mismatches=2,
                              sample_sheet=sample_sheet_file)
        ##"CATGCGCGGTA","TATGCGCGGTA","GATGCGCGGTA","GCTGCGCGGTA" = 307008
        ##"GCTGCGCGGTC" = 325394
        self.assertEqual(analysis.cutoff,None)
        self.assertEqual(analysis.mismatches,2)
        self.assertEqual(analysis.total_reads,632402)
        self.assertEqual(analysis.coverage,632402)
        self.assertEqual(analysis.barcodes,["GCTGCGCGGTC",
                                            "CATGCGCGGTA"])
        self.assertEqual(analysis.counts["GCTGCGCGGTC"].reads,325394)
        self.assertEqual(analysis.counts["CATGCGCGGTA"].reads,307008)
        self.assertEqual(analysis.counts["GCTGCGCGGTC"].sample,"SMPL2")
        self.assertEqual(analysis.counts["CATGCGCGGTA"].sample,"SMPL1")
        self.assertEqual(analysis.counts["GCTGCGCGGTC"].sequences,1)
        self.assertEqual(analysis.counts["CATGCGCGGTA"].sequences,4)
        # Create report
        reporter = report_barcodes(bc,
                                   lane=1,
                                   mismatches=2,
                                   sample_sheet=sample_sheet_file)
        # Check content
        self.assertEqual(str(reporter),
                         """Barcode analysis for lane #1
============================
 * Barcodes have been grouped by allowing 2 mismatches

#Rank	Index	Sample	N_seqs	N_reads	%reads	(%Total_reads)
    1	GCTGCGCGGTC	SMPL2	1	325394	51.5%	(51.5%)
    2	CATGCGCGGTA	SMPL1	4	307008	48.5%	(100.0%)""")

    def test_report_barcodes_for_no_counts(self):
        """report_barcodes: check output when there are no counts
        """
        bc = BarcodeCounter()
        analysis = bc.analyse()
        reporter = report_barcodes(bc)
        # Check content
        self.assertEqual(str(reporter),
                         """Barcode analysis for all lanes
==============================
 * No mismatches were allowed (exact matches only)
No barcodes counted""")
