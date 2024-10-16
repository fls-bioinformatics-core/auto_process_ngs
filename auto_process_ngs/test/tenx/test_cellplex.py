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
PJB1_GEX,/data/runs/fastqs_gex,any,PJB1,Gene Expression,
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
        self.assertEqual(config_csv.sections,['gene-expression',
                                              'libraries',
                                              'samples'])
        self.assertEqual(config_csv.reference_data_path,
                         "/data/refdata-cellranger-gex-GRCh38-2020-A")
        self.assertEqual(config_csv.probe_set_path,None)
        self.assertEqual(config_csv.feature_reference_path,None)
        self.assertEqual(config_csv.vdj_reference_path,None)
        self.assertEqual(config_csv.feature_types,
                         ['gene expression','multiplexing capture'])
        self.assertEqual(config_csv.gex_libraries,
                         ['PJB1_GEX',])
        self.assertEqual(config_csv.gex_library('PJB1_GEX'),
                         { 'fastqs': '/data/runs/fastqs_gex',
                           'lanes': 'any',
                           'library_id': 'PJB1',
                           'feature_type': 'Gene Expression',
                           'subsample_rate': ''
                         })
        self.assertEqual(config_csv.fastq_dirs,
                         { 'PJB1_GEX': '/data/runs/fastqs_gex',
                           'PJB2_MC': '/data/runs/fastqs_mc'
                         })
        self.assertEqual(config_csv.pretty_print_samples(),
                         "PBA, PBB")
        self.assertTrue(config_csv.is_valid)

    def test_cellranger_multi_config_csv_frp(self):
        """
        CellrangerMultiConfigCsv: check FRP data is extracted from config.csv
        """
        with open(os.path.join(self.wd,"10x_multi_config.csv"),'wt') as fp:
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A
probe-set,/data/Probe_Set_v1.0_GRCh38-2020-A.csv
no-bam,true

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB1_Flex,/data/runs/fastqs,any,PJB1,gene expression,

[samples]
sample_id,probe_barcode_ids,description
PB1,BC001,PB1
PB2,BC002,PB2
""")
        config_csv = CellrangerMultiConfigCsv(
            os.path.join(self.wd,
                         "10x_multi_config.csv"))
        self.assertEqual(config_csv.sample_names,['PB1','PB2'])
        self.assertEqual(config_csv.sections,['gene-expression',
                                              'libraries',
                                              'samples'])
        self.assertEqual(config_csv.reference_data_path,
                         "/data/refdata-cellranger-gex-GRCh38-2020-A")
        self.assertEqual(config_csv.probe_set_path,
                         "/data/Probe_Set_v1.0_GRCh38-2020-A.csv")
        self.assertEqual(config_csv.feature_reference_path,None)
        self.assertEqual(config_csv.vdj_reference_path,None)
        self.assertEqual(config_csv.feature_types,
                         ['gene expression'])
        self.assertEqual(config_csv.gex_libraries,
                         ['PJB1_Flex',])
        self.assertEqual(config_csv.gex_library('PJB1_Flex'),
                         { 'fastqs': '/data/runs/fastqs',
                           'lanes': 'any',
                           'library_id': 'PJB1',
                           'feature_type': 'gene expression',
                           'subsample_rate': ''
                         })
        self.assertEqual(config_csv.fastq_dirs,
                         {
                             'PJB1_Flex': '/data/runs/fastqs',
                         })
        self.assertEqual(config_csv.pretty_print_samples(),
                         "PB1-2")
        self.assertTrue(config_csv.is_valid)

    def test_cellranger_multi_config_csv_feature_ref(self):
        """
        CellrangerMultiConfigCsv: check feature reference is extracted from config.csv
        """
        with open(os.path.join(self.wd,"10x_multi_config.csv"),'wt') as fp:
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A

[feature]
reference,/data/feature_ref.csv

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB1_GEX,/data/runs/fastqs_gex,any,PJB1,Gene Expression,
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
        self.assertEqual(config_csv.sections,['feature',
                                              'gene-expression',
                                              'libraries',
                                              'samples'])
        self.assertEqual(config_csv.reference_data_path,
                         "/data/refdata-cellranger-gex-GRCh38-2020-A")
        self.assertEqual(config_csv.probe_set_path,None)
        self.assertEqual(config_csv.feature_reference_path,
                         "/data/feature_ref.csv")
        self.assertEqual(config_csv.vdj_reference_path,None)
        self.assertEqual(config_csv.feature_types,
                         ['gene expression','multiplexing capture'])
        self.assertEqual(config_csv.gex_libraries,
                         ['PJB1_GEX',])
        self.assertEqual(config_csv.gex_library('PJB1_GEX'),
                         { 'fastqs': '/data/runs/fastqs_gex',
                           'lanes': 'any',
                           'library_id': 'PJB1',
                           'feature_type': 'Gene Expression',
                           'subsample_rate': ''
                         })
        self.assertEqual(config_csv.fastq_dirs,
                         { 'PJB1_GEX': '/data/runs/fastqs_gex',
                           'PJB2_MC': '/data/runs/fastqs_mc'
                         })
        self.assertEqual(config_csv.pretty_print_samples(),
                         "PBA, PBB")
        self.assertTrue(config_csv.is_valid)

    def test_cellranger_multi_config_csv_vdj_t_ref(self):
        """
        CellrangerMultiConfigCsv: check V(D)J TCR data are extracted from config.csv
        """
        with open(os.path.join(self.wd,"10x_multi_config.csv"),'wt') as fp:
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A

[vdj]
reference,/data/vdj_ref.csv

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB1_GEX,/data/runs/fastqs_gex,any,PJB1,Gene Expression,
PJB2_TCR,/data/runs/fastqs_vdj_t,any,PJB2,VDJ-T,
""")
        config_csv = CellrangerMultiConfigCsv(
            os.path.join(self.wd,
                         "10x_multi_config.csv"))
        self.assertEqual(config_csv.sample_names,[])
        self.assertEqual(config_csv.sections,['gene-expression',
                                              'libraries',
                                              'vdj'])
        self.assertEqual(config_csv.reference_data_path,
                         "/data/refdata-cellranger-gex-GRCh38-2020-A")
        self.assertEqual(config_csv.probe_set_path,None)
        self.assertEqual(config_csv.feature_reference_path,None)
        self.assertEqual(config_csv.vdj_reference_path,
                         "/data/vdj_ref.csv")
        self.assertEqual(config_csv.feature_types,
                         ['gene expression','vdj-t'])
        self.assertEqual(config_csv.gex_libraries,
                         ['PJB1_GEX',])
        self.assertEqual(config_csv.gex_library('PJB1_GEX'),
                         { 'fastqs': '/data/runs/fastqs_gex',
                           'lanes': 'any',
                           'library_id': 'PJB1',
                           'feature_type': 'Gene Expression',
                           'subsample_rate': ''
                         })
        self.assertEqual(config_csv.libraries('VDJ-T'),
                         ['PJB2_TCR',])
        self.assertEqual(config_csv.library('VDJ-T','PJB2_TCR'),
                         { 'fastqs': '/data/runs/fastqs_vdj_t',
                           'lanes': 'any',
                           'library_id': 'PJB2',
                           'feature_type': 'VDJ-T',
                           'subsample_rate': ''
                         })
        self.assertEqual(config_csv.fastq_dirs,
                         { 'PJB1_GEX': '/data/runs/fastqs_gex',
                           'PJB2_TCR': '/data/runs/fastqs_vdj_t'
                         })
        self.assertEqual(config_csv.pretty_print_samples(),"")
        self.assertTrue(config_csv.is_valid)

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
        self.assertEqual(config_csv.sections,['gene-expression',
                                              'libraries',
                                              'samples'])
        self.assertEqual(config_csv.reference_data_path,
                         "/data/refdata-cellranger-gex-GRCh38-2020-A")
        self.assertEqual(config_csv.probe_set_path,None)
        self.assertEqual(config_csv.feature_reference_path,None)
        self.assertEqual(config_csv.vdj_reference_path,None)
        self.assertEqual(config_csv.feature_types,
                         ['gene expression','multiplexing capture'])
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
        self.assertTrue(config_csv.is_valid)

    def test_cellranger_multi_config_csv_duplicated_sample(self):
        """
        CellrangerMultiConfigCsv: handle for duplicated sample
        """
        with open(os.path.join(self.wd,"10x_multi_config.csv"),'wt') as fp:
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A

[vdj]
reference,/data/vdj_ref.csv

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB1_GEX,/data/runs/fastqs_gex,any,PJB1,Gene Expression,
PJB2_CML,/data/runs/fastqs_cml,any,PJB2,Multiplexing Capture,

[samples]
PB1,CMO1,PB1
PB2,CMO2,PB2
PB1,CMO3,PB1
""")
        # Raise exception by default
        self.assertRaises(Exception,
                          CellrangerMultiConfigCsv,
                          os.path.join(self.wd,"10x_multi_config.csv"))
        # Invalid if loaded with 'strict' turned off
        config_csv = CellrangerMultiConfigCsv(
            os.path.join(self.wd,
                         "10x_multi_config.csv"), strict=False)
        self.assertFalse(config_csv.is_valid)
        self.assertEqual(len(config_csv.get_errors()), 1)

    def test_cellranger_multi_config_csv_duplicated_cmo(self):
        """
        CellrangerMultiConfigCsv: handle duplicated CMO
        """
        with open(os.path.join(self.wd,"10x_multi_config.csv"),'wt') as fp:
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A

[vdj]
reference,/data/vdj_ref.csv

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB1_GEX,/data/runs/fastqs_gex,any,PJB1,Gene Expression,
PJB2_CML,/data/runs/fastqs_cml,any,PJB2,Multiplexing Capture,

[samples]
PB1,CMO1|CMO2,PB1
PB2,CMO2|CMO3,PB2
""")
        # Raise exception by default
        self.assertRaises(Exception,
                          CellrangerMultiConfigCsv,
                          os.path.join(self.wd,"10x_multi_config.csv"))
        # Invalid if loaded with 'strict' turned off
        config_csv = CellrangerMultiConfigCsv(
            os.path.join(self.wd,
                         "10x_multi_config.csv"), strict=False)
        self.assertFalse(config_csv.is_valid)
        self.assertEqual(len(config_csv.get_errors()), 1)

    def test_cellranger_multi_config_csv_bad_cmo(self):
        """
        CellrangerMultiConfigCsv: handle "bad" CMO
        """
        with open(os.path.join(self.wd,"10x_multi_config.csv"),'wt') as fp:
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A

[vdj]
reference,/data/vdj_ref.csv

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB1_GEX,/data/runs/fastqs_gex,any,PJB1,Gene Expression,
PJB2_CML,/data/runs/fastqs_cml,any,PJB2,Multiplexing Capture,

[samples]
PB1,CMO1|CMO2...,PB1
""")
        # Raise exception by default
        self.assertRaises(Exception,
                          CellrangerMultiConfigCsv,
                          os.path.join(self.wd,"10x_multi_config.csv"))
        # Invalid if loaded with 'strict' turned off
        config_csv = CellrangerMultiConfigCsv(
            os.path.join(self.wd,
                         "10x_multi_config.csv"), strict=False)
        self.assertFalse(config_csv.is_valid)
        self.assertEqual(len(config_csv.get_errors()), 1)

    def test_cellranger_multi_config_csv_unknown_feature_type(self):
        """
        CellrangerMultiConfigCsv: handle unknown feature type
        """
        with open(os.path.join(self.wd,"10x_multi_config.csv"),'wt') as fp:
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A

[vdj]
reference,/data/vdj_ref.csv

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB1_GEX,/data/runs/fastqs_gex,any,PJB1,Gene Expression,
PJB2_UNKNOWN,/data/runs/fastqs_unknown,any,PJB2,Unknown,
""")
        # Raise exception by default
        self.assertRaises(Exception,
                          CellrangerMultiConfigCsv,
                          os.path.join(self.wd,"10x_multi_config.csv"))
        # Invalid if loaded with 'strict' turned off
        config_csv = CellrangerMultiConfigCsv(
            os.path.join(self.wd,
                         "10x_multi_config.csv"), strict=False)
        self.assertFalse(config_csv.is_valid)
        self.assertEqual(len(config_csv.get_errors()), 1)

    def test_cellranger_multi_config_csv_duplicated_feature_type(self):
        """
        CellrangerMultiConfigCsv: handle duplicated feature type
        """
        with open(os.path.join(self.wd,"10x_multi_config.csv"),'wt') as fp:
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A

[vdj]
reference,/data/vdj_ref.csv

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB1_GEX,/data/runs/fastqs,any,PJB2,Gene Expression,
PJB1_CML,/data/runs/fastqs,any,PJB2,Multiplexing capture,
PJB2_GEX,/data/runs/fastqs,any,PJB2,Gene Expression,
PJB2_CML,/data/runs/fastqs,any,PJB2,Multiplexing capture,
""")
        # Raise exception by default
        self.assertRaises(Exception,
                          CellrangerMultiConfigCsv,
                          os.path.join(self.wd,"10x_multi_config.csv"))
        # Invalid if loaded with 'strict' turned off
        config_csv = CellrangerMultiConfigCsv(
            os.path.join(self.wd,
                         "10x_multi_config.csv"), strict=False)
        self.assertFalse(config_csv.is_valid)
        self.assertEqual(len(config_csv.get_errors()), 2)
