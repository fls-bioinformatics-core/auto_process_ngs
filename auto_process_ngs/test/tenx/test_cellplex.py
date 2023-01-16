#######################################################################
# Tests for tenx_genomics/multiome.py module
#######################################################################

import unittest
import os
import shutil
import tempfile
from auto_process_ngs.tenx.cellplex import *

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

class TestCellrangerMultiConfigCsv(unittest.TestCase):
    """
    Tests for the 'CellrangerMultiConfigCsv' class
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestCellrangerMultiConfigCsv')
    def tearDown(self):
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    def test_cellranger_multi_config_csv(self):
        """
        CellrangerMultiConfigCsv: check data is extracted from config.csv
        """
        with open(os.path.join(self.wd,"10x_multi_config.csv"),'wt') as fp:
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB1_GEX,/data/runs/fastqs_gex,any,PJB1,gene expression,
PJB2_MC,/data/runs/fastqs_mc,any,PJB2,Multiplexing Capture,

[samples]
sample_id,cmo_ids,description
PBA,CMO301,PBA
PBB,CMO302,PBB
""")
        config_csv = CellrangerMultiConfigCsv(
            os.path.join(self.wd,
                         "10x_multi_config.csv"))
        self.assertEqual(config_csv.sample_names,['PBA','PBB'])
        self.assertEqual(config_csv.reference_data_path,
                         "/data/refdata-cellranger-gex-GRCh38-2020-A")
        self.assertEqual(config_csv.gex_libraries,
                         ['PJB1_GEX',])
        self.assertEqual(config_csv.gex_library('PJB1_GEX'),
                         { 'fastqs': '/data/runs/fastqs_gex',
                           'lanes': 'any',
                           'library_id': 'PJB1',
                           'feature_type': 'gene expression',
                           'subsample_rate': ''
                         })
        self.assertEqual(config_csv.fastq_dirs,
                         { 'PJB1_GEX': '/data/runs/fastqs_gex',
                           'PJB2_MC': '/data/runs/fastqs_mc'
                         })
        self.assertEqual(config_csv.pretty_print_samples(),
                         "PBA, PBB")

    def test_cellranger_multi_config_csv_reduced_fields(self):
        """
        CellrangerMultiConfigCsv: check with 'reduced' config.csv
        """
        with open(os.path.join(self.wd,"10x_multi_config.csv"),'wt') as fp:
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A

[libraries]
fastq_id,fastqs,feature_types
PJB1_GEX,/data/runs/fastqs_gex,gene expression
PJB2_MC,/data/runs/fastqs_mc,Multiplexing Capture

[samples]
sample_id,cmo_ids
PBA,CMO301
PBB,CMO302
""")
        config_csv = CellrangerMultiConfigCsv(
            os.path.join(self.wd,
                         "10x_multi_config.csv"))
        self.assertEqual(config_csv.sample_names,['PBA','PBB'])
        self.assertEqual(config_csv.reference_data_path,
                         "/data/refdata-cellranger-gex-GRCh38-2020-A")
        self.assertEqual(config_csv.gex_libraries,
                         ['PJB1_GEX',])
        self.assertEqual(config_csv.gex_library('PJB1_GEX'),
                         { 'fastqs': '/data/runs/fastqs_gex',
                           'lanes': 'any',
                           'library_id': None,
                           'feature_type': 'gene expression',
                           'subsample_rate': ''
                         })
        self.assertEqual(config_csv.fastq_dirs,
                         { 'PJB1_GEX': '/data/runs/fastqs_gex',
                           'PJB2_MC': '/data/runs/fastqs_mc'
                         })
        self.assertEqual(config_csv.pretty_print_samples(),
                         "PBA, PBB")
