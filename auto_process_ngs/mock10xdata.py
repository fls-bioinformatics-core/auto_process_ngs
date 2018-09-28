#     mock.py: module providing mock 10xGenomics data for testing
#     Copyright (C) University of Manchester 2018 Peter Briggs
#
########################################################################

QC_SUMMARY_JSON = """{
    "10x_software_version": "cellranger 2.2.0", 
    "PhiX_aligned": null, 
    "PhiX_error_worst_tile": null, 
    "bcl2fastq_args": "bcl2fastq --minimum-trimmed-read-length 8 --mask-short-adapter-reads 8 --create-fastq-for-index-reads --ignore-missing-positions --ignore-missing-filter --ignore-missing-bcls --use-bases-mask=y26,I8,y98 -R /net/lustre/mqbsspbd/180926_MN00218_0021_A000H2GYCG_analysis/primary_data/180926_MN00218_0021_A000H2GYCG --output-dir=/scratch/mqbsspbd/180926_MN00218_0021_A000H2GYCG_analysis/bcl2fastq --interop-dir=/scratch/mqbsspbd/180926_MN00218_0021_A000H2GYCG_analysis/000H2GYCG/MAKE_FASTQS_CS/MAKE_FASTQS/BCL2FASTQ_WITH_SAMPLESHEET/fork0/chnk0-u4f20ace4f5/files/interop_path --sample-sheet=/scratch/mqbsspbd/180926_MN00218_0021_A000H2GYCG_analysis/000H2GYCG/MAKE_FASTQS_CS/MAKE_FASTQS/PREPARE_SAMPLESHEET/fork0/chnk0-u4f20ace4cb/files/samplesheet.csv -p 6 -d 6 -r 6 -w 6", 
    "bcl2fastq_version": "2.17.1.14", 
    "experiment_name": "", 
    "fwhm_A": null, 
    "fwhm_C": null, 
    "fwhm_G": null, 
    "fwhm_T": null, 
    "index1_PhiX_error_by_cycle": null, 
    "index1_mean_phasing": null, 
    "index1_mean_prephasing": null, 
    "index1_q20_fraction": null, 
    "index1_q20_fraction_by_cycle": null, 
    "index1_q30_fraction": null, 
    "index1_q30_fraction_by_cycle": null, 
    "index2_PhiX_error_by_cycle": null, 
    "index2_mean_phasing": null, 
    "index2_mean_prephasing": null, 
    "index2_q20_fraction": null, 
    "index2_q20_fraction_by_cycle": null, 
    "index2_q30_fraction": null, 
    "index2_q30_fraction_by_cycle": null, 
    "intensity_A": null, 
    "intensity_C": null, 
    "intensity_G": null, 
    "intensity_T": null, 
    "lanecount": 1, 
    "mean_cluster_density": null, 
    "mean_cluster_density_pf": null, 
    "num_clusters": null, 
    "percent_pf_clusters": null, 
    "read1_PhiX_error_by_cycle": null, 
    "read1_mean_phasing": null, 
    "read1_mean_prephasing": null, 
    "read1_q20_fraction": null, 
    "read1_q20_fraction_by_cycle": null, 
    "read1_q30_fraction": null, 
    "read1_q30_fraction_by_cycle": null, 
    "read2_PhiX_error_by_cycle": null, 
    "read2_mean_phasing": null, 
    "read2_mean_prephasing": null, 
    "read2_q20_fraction": null, 
    "read2_q20_fraction_by_cycle": null, 
    "read2_q30_fraction": null, 
    "read2_q30_fraction_by_cycle": null, 
    "rta_version": "unknown", 
    "run_id": "180926_MN00218_0021_A000H2GYCG", 
    "sample_qc": {
        "smpl1": {
            "1": {
                "barcode_exact_match_ratio": 0.940013016010283, 
                "barcode_q30_base_ratio": 0.9468684602040338, 
                "bc_on_whitelist": 0.9634447859958355, 
                "gem_count_estimate": 68214, 
                "mean_barcode_qscore": 35.47116735755179, 
                "number_reads": 3388135, 
                "read1_q30_base_ratio": 0.9460734878931605, 
                "read2_q30_base_ratio": 0.7687690301946499
            }, 
            "all": {
                "barcode_exact_match_ratio": 0.940013016010283, 
                "barcode_q30_base_ratio": 0.9468684602040338, 
                "bc_on_whitelist": 0.9634447859958355, 
                "gem_count_estimate": 68214, 
                "mean_barcode_qscore": 35.47116735755179, 
                "number_reads": 3388135, 
                "read1_q30_base_ratio": 0.9460734878931605, 
                "read2_q30_base_ratio": 0.7687690301946499
            }
        }, 
        "smpl2": {
            "1": {
                "barcode_exact_match_ratio": 0.9412661883398474, 
                "barcode_q30_base_ratio": 0.9483778169720642, 
                "bc_on_whitelist": 0.9653697937148712, 
                "gem_count_estimate": 21399, 
                "mean_barcode_qscore": 35.50137749763865, 
                "number_reads": 2990424, 
                "read1_q30_base_ratio": 0.9474780353080827, 
                "read2_q30_base_ratio": 0.7802990769663007
            }, 
            "all": {
                "barcode_exact_match_ratio": 0.9412661883398474, 
                "barcode_q30_base_ratio": 0.9483778169720642, 
                "bc_on_whitelist": 0.9653697937148712, 
                "gem_count_estimate": 21399, 
                "mean_barcode_qscore": 35.50137749763865, 
                "number_reads": 2990424, 
                "read1_q30_base_ratio": 0.9474780353080827, 
                "read2_q30_base_ratio": 0.7802990769663007
            }
        }
    }, 
    "signoise_ratio": null, 
    "start_datetime": "None", 
    "surfacecount": 2, 
    "swathcount": 3, 
    "tilecount": 6, 
    "total_cluster_density": null, 
    "total_cluster_density_pf": null, 
    "yield": null, 
    "yield_pf": null, 
    "yield_pf_q30": null
}
"""

CELLRANGER_QC_SUMMARY = """<html>
<head>
<title>mkfastq QC report</title>
<style type="text/css">

h1    { background-color: #42AEC2;
        color: white;
        padding: 5px 10px; }
h2    { background-color: #8CC63F;
        color: white;
        display: inline-block;
        padding: 5px 15px;
        margin: 0;
        border-top-left-radius: 20px;
        border-bottom-right-radius: 20px; }
h3, h4 { background-color: grey;
         color: white;
         display: block;
         padding: 5px 15px;
         margin: 0;
         border-top-left-radius: 20px;
         border-bottom-right-radius: 20px; }
table { margin: 10 10;
        border: solid 1px grey;
        background-color: white; }
th    { background-color: grey;
        color: white;
        padding: 2px 5px; }
td    { text-align: left;
        vertical-align: top;
        padding: 2px 5px;
        border-bottom: solid 1px lightgray; }
td.param { background-color: grey;
           color: white;
           padding: 2px 5px;
           font-weight: bold; }
p.footer { font-style: italic;
           font-size: 70%; }

table { font-size: 80%;
        font-family: sans-serif; }
td { text-align: right; }</style>
</head>
<body>
<h1>mkfastq QC report</h1>
<div id='toc'>
<h2>Contents</h2>
<ul><li><a href='#General_information'>General information</a></li><li><a href='#Sample_smpl1'>Sample: smpl1</a></li><li><a href='#Sample_smpl2'>Sample: smpl2</a></li></ul>
</div>
<div id='General_information'>
<h2>General information</h2>
<table>
<tr><th>Parameter</th><th>Value</th></tr>
<tr><td>run_id</td><td>180926_MN00218_0021_A000H2GYCG</td></tr>
<tr><td>experiment_name</td><td></td></tr>
<tr><td>10x_software_version</td><td>cellranger 2.2.0</td></tr>
<tr><td>bcl2fastq_version</td><td>2.17.1.14</td></tr>
<tr><td>bcl2fastq_args</td><td>bcl2fastq --minimum-trimmed-read-length 8 --mask-short-adapter-reads 8 --create-fastq-for-index-reads --ignore-missing-positions --ignore-missing-filter --ignore-missing-bcls --use-bases-mask=y26,I8,y98 -R /net/lustre/mqbsspbd/180926_MN00218_0021_A000H2GYCG_analysis/primary_data/180926_MN00218_0021_A000H2GYCG --output-dir=/scratch/mqbsspbd/180926_MN00218_0021_A000H2GYCG_analysis/bcl2fastq --interop-dir=/scratch/mqbsspbd/180926_MN00218_0021_A000H2GYCG_analysis/000H2GYCG/MAKE_FASTQS_CS/MAKE_FASTQS/BCL2FASTQ_WITH_SAMPLESHEET/fork0/chnk0-u4f20ace4f5/files/interop_path --sample-sheet=/scratch/mqbsspbd/180926_MN00218_0021_A000H2GYCG_analysis/000H2GYCG/MAKE_FASTQS_CS/MAKE_FASTQS/PREPARE_SAMPLESHEET/fork0/chnk0-u4f20ace4cb/files/samplesheet.csv -p 6 -d 6 -r 6 -w 6</td></tr>
<tr><td>rta_version</td><td>unknown</td></tr>
</table>
</div>
<div id='Sample_smpl1'>
<h2>Sample: smpl1</h2>
<table>
<tr><th></th><th>1</th><th>all</th></tr>
<tr><td>bc_on_whitelist</td><td>0.966331905673</td><td>0.966331905673</td></tr>
<tr><td>barcode_exact_match_ratio</td><td>0.942857435466</td><td>0.942857435466</td></tr>
<tr><td>read2_q30_base_ratio</td><td>0.772020178873</td><td>0.772020178873</td></tr>
<tr><td>gem_count_estimate</td><td>39745</td><td>39745</td></tr>
<tr><td>number_reads</td><td>2831672</td><td>2831672</td></tr>
<tr><td>read1_q30_base_ratio</td><td>0.946857262157</td><td>0.946857262157</td></tr>
<tr><td>mean_barcode_qscore</td><td>35.4841823647</td><td>35.4841823647</td></tr>
<tr><td>barcode_q30_base_ratio</td><td>0.947610600507</td><td>0.947610600507</td></tr>
</table>
</div>
<div id='Sample_smpl2'>
<h2>Sample: smpl2</h2>
<table>
<tr><th></th><th>1</th><th>all</th></tr>
<tr><td>bc_on_whitelist</td><td>0.967838002215</td><td>0.967838002215</td></tr>
<tr><td>barcode_exact_match_ratio</td><td>0.943587582063</td><td>0.943587582063</td></tr>
<tr><td>read2_q30_base_ratio</td><td>0.772472408134</td><td>0.772472408134</td></tr>
<tr><td>gem_count_estimate</td><td>13063</td><td>13063</td></tr>
<tr><td>number_reads</td><td>2972732</td><td>2972732</td></tr>
<tr><td>read1_q30_base_ratio</td><td>0.946434513574</td><td>0.946434513574</td></tr>
<tr><td>mean_barcode_qscore</td><td>35.4759437333</td><td>35.4759437333</td></tr>
<tr><td>barcode_q30_base_ratio</td><td>0.947180621102</td><td>0.947180621102</td></tr>
</table>
</div></body>
</html>
"""
