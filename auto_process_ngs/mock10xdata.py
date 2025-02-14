#     mock.py: module providing mock 10xGenomics data for testing
#     Copyright (C) University of Manchester 2018-2025 Peter Briggs
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

# CellRanger
METRICS_SUMMARY = """Estimated Number of Cells,Mean Reads per Cell,Median Genes per Cell,Number of Reads,Valid Barcodes,Reads Mapped Confidently to Transcriptome,Reads Mapped Confidently to Exonic Regions,Reads Mapped Confidently to Intronic Regions,Reads Mapped Confidently to Intergenic Regions,Reads Mapped Antisense to Gene,Sequencing Saturation,Q30 Bases in Barcode,Q30 Bases in RNA Read,Q30 Bases in Sample Index,Q30 Bases in UMI,Fraction Reads in Cells,Total Genes Detected,Median UMI Counts per Cell
"2,272","107,875","1,282","245,093,084",98.3%,69.6%,71.9%,6.1%,3.2%,4.4%,51.3%,98.5%,79.2%,93.6%,98.5%,12.0%,"16,437","2,934"
"""

# CellRanger 7.1.0
METRICS_SUMMARY_7_1_0 = """Estimated Number of Cells,Mean Reads per Cell,Median Genes per Cell,Number of Reads,Valid Barcodes,Sequencing Saturation,Q30 Bases in Barcode,Q30 Bases in RNA Read,Q30 Bases in UMI,Reads Mapped to Genome,Reads Mapped Confidently to Genome,Reads Mapped Confidently to Intergenic Regions,Reads Mapped Confidently to Intronic Regions,Reads Mapped Confidently to Exonic Regions,Reads Mapped Confidently to Transcriptome,Reads Mapped Antisense to Gene,Fraction Reads in Cells,Total Genes Detected,Median UMI Counts per Cell
"5,529",551,267,"3,045,784",97.6%,1.7%,97.4%,89.4%,97.1%,94.3%,89.0%,3.9%,30.6%,54.5%,76.5%,8.0%,88.5%,"17,781",331
"""

# Cellranger ATAC < 2.0.0
ATAC_SUMMARY = """annotated_cells,bc_q30_bases_fract,cellranger-atac_version,cells_detected,frac_cut_fragments_in_peaks,frac_fragments_nfr,frac_fragments_nfr_or_nuc,frac_fragments_nuc,frac_fragments_overlapping_peaks,frac_fragments_overlapping_targets,frac_mapped_confidently,frac_waste_chimeric,frac_waste_duplicate,frac_waste_lowmapq,frac_waste_mitochondrial,frac_waste_no_barcode,frac_waste_non_cell_barcode,frac_waste_overall_nondup,frac_waste_total,frac_waste_unmapped,median_fragments_per_cell,median_per_cell_unique_fragments_at_30000_RRPC,median_per_cell_unique_fragments_at_50000_RRPC,num_fragments,r1_q30_bases_fract,r2_q30_bases_fract,si_q30_bases_fract,total_usable_fragments,tss_enrichment_score
5682,0.925226023701,1.0.1,6748,0.512279447992,0.392368676637,0.851506103882,0.459137427245,0.556428090013,0.575082094792,0.534945791083,0.00123066129161,0.160515305655,0.0892973647982,0.00899493352094,0.352907229061,0.0135851297269,0.471714266123,0.632229571777,0.00569894772443,16119.5,5769.94794925,8809.29425158,366582587,0.947387774999,0.941378123188,0.962708567847,134818235,6.91438390781
"""

# Cellranger ATAC 2.0.0
ATAC_SUMMARY_2_0_0 = """Sample ID,Genome,Pipeline version,Estimated number of cells,Confidently mapped read pairs,Estimated bulk library complexity,Fraction of genome in peaks,Fraction of high-quality fragments in cells,Fraction of high-quality fragments overlapping TSS,Fraction of high-quality fragments overlapping peaks,Fraction of transposition events in peaks in cells,Fragments flanking a single nucleosome,Fragments in nucleosome-free regions,Mean raw read pairs per cell,Median high-quality fragments per cell,Non-nuclear read pairs,Number of peaks,Percent duplicates,Q30 bases in barcode,Q30 bases in read 1,Q30 bases in read 2,Q30 bases in sample index i1,Sequenced read pairs,Sequencing saturation,TSS enrichment score,Unmapped read pairs,Valid barcodes
ATAC_1,GRCh38,cellranger-atac-2.0.0,3582,0.9223,467083069.2000,0.0631,0.8358,0.2910,0.4856,0.4571,0.3372,0.5777,111103.6343,51354.5000,0.0034,228838,0.3622,0.9328,0.9640,0.9348,0.9632,397973218,0.4853,6.8333,0.0102,0.9747
"""

MULTIOME_SUMMARY = """Sample ID,Pipeline version,Genome,Estimated number of cells,ATAC Confidently mapped read pairs,ATAC Fraction of genome in peaks,ATAC Fraction of high-quality fragments in cells,ATAC Fraction of high-quality fragments overlapping TSS,ATAC Fraction of high-quality fragments overlapping peaks,ATAC Fraction of transposition events in peaks in cells,ATAC Mean raw read pairs per cell,ATAC Median high-quality fragments per cell,ATAC Non-nuclear read pairs,ATAC Number of peaks,ATAC Percent duplicates,ATAC Q30 bases in barcode,ATAC Q30 bases in read 1,ATAC Q30 bases in read 2,ATAC Q30 bases in sample index i1,ATAC Sequenced read pairs,ATAC TSS enrichment score,ATAC Unmapped read pairs,ATAC Valid barcodes,Feature linkages detected,GEX Fraction of transcriptomic reads in cells,GEX Mean raw reads per cell,GEX Median UMI counts per cell,GEX Median genes per cell,GEX Percent duplicates,GEX Q30 bases in UMI,GEX Q30 bases in barcode,GEX Q30 bases in read 2,GEX Q30 bases in sample index i1,GEX Q30 bases in sample index i2,GEX Reads mapped antisense to gene,GEX Reads mapped confidently to exonic regions,GEX Reads mapped confidently to genome,GEX Reads mapped confidently to intergenic regions,GEX Reads mapped confidently to intronic regions,GEX Reads mapped confidently to transcriptome,GEX Reads mapped to genome,GEX Reads with TSO,GEX Sequenced read pairs,GEX Total genes detected,GEX Valid UMIs,GEX Valid barcodes,Linked genes,Linked peaks
PB_ATAC_multiome,cellranger-arc-1.0.0,GRCh38,744,0.892,0.02691,0.54188,0.4003,0.59379,0.57403,404636.24194,8079,0.02154,80637,0.942,0.96564,0.96911,0.96414,0.97688,301049364,9.00247,0.01356,0.98281,20374,0.36971,569488.26747,5441,1490,0.9288,0.96798,0.97522,0.887,0.98441,0.97715,0.0497,0.48039,0.76483,0.03809,0.24635,0.67216,0.89641,0.2635,423699271,24699,0.99958,0.96795,3930,14251
"""

MULTIOME_SUMMARY_2_0_0 = """Sample ID,Genome,Pipeline version,Estimated number of cells,Feature linkages detected,Linked genes,Linked peaks,ATAC Confidently mapped read pairs,ATAC Fraction of genome in peaks,ATAC Fraction of high-quality fragments in cells,ATAC Fraction of high-quality fragments overlapping TSS,ATAC Fraction of high-quality fragments overlapping peaks,ATAC Fraction of transposition events in peaks in cells,ATAC Mean raw read pairs per cell,ATAC Median high-quality fragments per cell,ATAC Non-nuclear read pairs,ATAC Number of peaks,ATAC Percent duplicates,ATAC Q30 bases in barcode,ATAC Q30 bases in read 1,ATAC Q30 bases in read 2,ATAC Q30 bases in sample index i1,ATAC Sequenced read pairs,ATAC TSS enrichment score,ATAC Unmapped read pairs,ATAC Valid barcodes,GEX Fraction of transcriptomic reads in cells,GEX Mean raw reads per cell,GEX Median UMI counts per cell,GEX Median genes per cell,GEX Percent duplicates,GEX Q30 bases in UMI,GEX Q30 bases in barcode,GEX Q30 bases in read 2,GEX Reads mapped antisense to gene,GEX Reads mapped confidently to exonic regions,GEX Reads mapped confidently to genome,GEX Reads mapped confidently to intergenic regions,GEX Reads mapped confidently to intronic regions,GEX Reads mapped confidently to transcriptome,GEX Reads mapped to genome,GEX Reads with TSO,GEX Sequenced read pairs,GEX Total genes detected,GEX Valid UMIs,GEX Valid barcodes
PB_ATAC_multiome,GRCh38,cellranger-arc-2.0.0,785,2,2,2,0.7659,0.0125,0.7934,0.5388,0.5809,0.5591,127.3885,9.0000,0.6139,167,0.0624,0.8433,0.9228,0.9009,0.8869,100000,40.0000,0.0881,0.9838,0.8720,509.5541,161.0000,15.0000,0.3135,0.9698,0.9707,0.9550,0.0241,0.6598,0.8938,0.2190,0.0150,0.6503,0.9612,0.1222,400000,135,1.0000,0.9956
"""

CELLPLEX_METRICS_SUMMARY = """Library or Sample,Library Type,Grouped By,Group Name,Metric Name,Metric Value
Library,Gene Expression,Fastq ID,JR_10K_1_gex,Number of reads,"384,834,949"
Library,Gene Expression,Fastq ID,JR_10K_1_gex,Number of short reads skipped,0
Library,Gene Expression,Fastq ID,JR_10K_1_gex,Q30 RNA read,93.5%
Library,Gene Expression,Fastq ID,JR_10K_1_gex,Q30 UMI,94.3%
Library,Gene Expression,Fastq ID,JR_10K_1_gex,Q30 barcodes,95.8%
Library,Gene Expression,Physical library ID,JR_10K_gex,Cell-associated barcodes identified as multiplets,903 (6.65%)
Library,Gene Expression,Physical library ID,JR_10K_gex,Cell-associated barcodes not assigned any CMOs,1635 (12.03%)
Library,Gene Expression,Physical library ID,JR_10K_gex,Cells assigned to a sample,11050 (81.32%)
Library,Gene Expression,Physical library ID,JR_10K_gex,Confidently mapped antisense,1.29%
Library,Gene Expression,Physical library ID,JR_10K_gex,Confidently mapped to exonic regions,66.79%
Library,Gene Expression,Physical library ID,JR_10K_gex,Confidently mapped to genome,94.53%
Library,Gene Expression,Physical library ID,JR_10K_gex,Confidently mapped to intergenic regions,4.18%
Library,Gene Expression,Physical library ID,JR_10K_gex,Confidently mapped to intronic regions,23.57%
Library,Gene Expression,Physical library ID,JR_10K_gex,Confidently mapped to transcriptome,64.13%
Library,Gene Expression,Physical library ID,JR_10K_gex,Estimated number of cells,"13,588"
Library,Gene Expression,Physical library ID,JR_10K_gex,Fraction reads in cells,92.18%
Library,Gene Expression,Physical library ID,JR_10K_gex,Mapped to genome,97.28%
Library,Gene Expression,Physical library ID,JR_10K_gex,Mean reads per cell,"28,322"
Library,Gene Expression,Physical library ID,JR_10K_gex,Number of reads,"384,834,949"
Library,Gene Expression,Physical library ID,JR_10K_gex,Number of reads in the library,"384,834,949"
Library,Gene Expression,Physical library ID,JR_10K_gex,Sequencing saturation,20.38%
Library,Gene Expression,Physical library ID,JR_10K_gex,Valid UMIs,99.96%
Library,Gene Expression,Physical library ID,JR_10K_gex,Valid barcodes,98.53%
Library,Multiplexing Capture,,,Cell-associated barcodes identified as multiplets,903 (6.65%)
Library,Multiplexing Capture,,,Cells assigned to a sample,11050 (81.32%)
Library,Multiplexing Capture,,,Estimated number of cell-associated barcodes,"13,588"
Library,Multiplexing Capture,,,Median CMO UMIs per cell,"7,963"
Library,Multiplexing Capture,,,Number of samples assigned at least one cell,2
Library,Multiplexing Capture,,,Singlet capture ratio,0.87
Library,Multiplexing Capture,CMO Name,CMO301,CMO signal-to-noise ratio,5.09
Library,Multiplexing Capture,CMO Name,CMO301,Cells assigned to CMO,53.17%
Library,Multiplexing Capture,CMO Name,CMO301,Fraction reads in cell-associated barcodes,70.78%
Library,Multiplexing Capture,CMO Name,CMO302,CMO signal-to-noise ratio,3.57
Library,Multiplexing Capture,CMO Name,CMO302,Cells assigned to CMO,46.83%
Library,Multiplexing Capture,CMO Name,CMO302,Fraction reads in cell-associated barcodes,55.19%
Library,Multiplexing Capture,Fastq ID,JR_10K_1_multiplex,Number of reads,"215,531,826"
Library,Multiplexing Capture,Fastq ID,JR_10K_1_multiplex,Number of short reads skipped,0
Library,Multiplexing Capture,Fastq ID,JR_10K_1_multiplex,Q30 RNA read,96.8%
Library,Multiplexing Capture,Fastq ID,JR_10K_1_multiplex,Q30 UMI,95.5%
Library,Multiplexing Capture,Fastq ID,JR_10K_1_multiplex,Q30 barcodes,95.8%
Library,Multiplexing Capture,Fastq ID,JR_10K_2_multiplex,Number of reads,"538,688"
Library,Multiplexing Capture,Fastq ID,JR_10K_2_multiplex,Number of short reads skipped,0
Library,Multiplexing Capture,Fastq ID,JR_10K_2_multiplex,Q30 RNA read,97.1%
Library,Multiplexing Capture,Fastq ID,JR_10K_2_multiplex,Q30 UMI,95.9%
Library,Multiplexing Capture,Fastq ID,JR_10K_2_multiplex,Q30 barcodes,95.9%
Library,Multiplexing Capture,Physical library ID,JR_10K_multiplex,Cell-associated barcodes identified as multiplets,903 (6.65%)
Library,Multiplexing Capture,Physical library ID,JR_10K_multiplex,Cell-associated barcodes not assigned any CMOs,1635 (12.03%)
Library,Multiplexing Capture,Physical library ID,JR_10K_multiplex,Cells assigned to a sample,11050 (81.32%)
Library,Multiplexing Capture,Physical library ID,JR_10K_multiplex,Estimated number of cell-associated barcodes,"13,588"
Library,Multiplexing Capture,Physical library ID,JR_10K_multiplex,Fraction CMO reads,98.30%
Library,Multiplexing Capture,Physical library ID,JR_10K_multiplex,Fraction CMO reads usable,63.49%
Library,Multiplexing Capture,Physical library ID,JR_10K_multiplex,Fraction reads from multiplets,7.07%
Library,Multiplexing Capture,Physical library ID,JR_10K_multiplex,Fraction reads in cell-associated barcodes,64.85%
Library,Multiplexing Capture,Physical library ID,JR_10K_multiplex,Fraction unrecognized CMO,1.70%
Library,Multiplexing Capture,Physical library ID,JR_10K_multiplex,Mean reads per cell-associated barcode,"15,568"
Library,Multiplexing Capture,Physical library ID,JR_10K_multiplex,Median CMO UMIs per cell-associated barcode,"7,735"
Library,Multiplexing Capture,Physical library ID,JR_10K_multiplex,Number of reads,"216,070,514"
Library,Multiplexing Capture,Physical library ID,JR_10K_multiplex,Samples assigned at least one cell,2
Library,Multiplexing Capture,Physical library ID,JR_10K_multiplex,Valid UMIs,99.99%
Library,Multiplexing Capture,Physical library ID,JR_10K_multiplex,Valid barcodes,99.59%
Sample,Gene Expression,,,Cells,"5,175"
Sample,Gene Expression,,,Confidently mapped antisense,1.38%
Sample,Gene Expression,,,Confidently mapped to exonic regions,68.40%
Sample,Gene Expression,,,Confidently mapped to genome,95.05%
Sample,Gene Expression,,,Confidently mapped to intergenic regions,3.77%
Sample,Gene Expression,,,Confidently mapped to intronic regions,22.88%
Sample,Gene Expression,,,Confidently mapped to transcriptome,65.67%
Sample,Gene Expression,,,Mapped to genome,97.66%
Sample,Gene Expression,,,Median UMI counts per cell,"10,515"
Sample,Gene Expression,,,Median genes per cell,"3,086"
Sample,Gene Expression,,,Median reads per cell,"20,052"
Sample,Gene Expression,,,Number of reads assigned to the sample,"114,898,392"
Sample,Gene Expression,,,Total genes detected,"21,260"
Sample,Gene Expression,Physical library ID,JR_10K_gex,Cell-associated barcodes identified as multiplets,903 (6.65%)
Sample,Gene Expression,Physical library ID,JR_10K_gex,Cell-associated barcodes not assigned any CMOs,1635 (12.03%)
Sample,Gene Expression,Physical library ID,JR_10K_gex,Cells assigned to other samples,5875 (43.24%)
Sample,Gene Expression,Physical library ID,JR_10K_gex,Cells assigned to this sample,5175 (38.09%)
"""

CELLPLEX_METRICS_SUMMARY_7_1_0 = """Category,Library Type,Grouped By,Group Name,Metric Name,Metric Value
Cells,Gene Expression,,,Cells,"1,569"
Cells,Gene Expression,,,Confidently mapped antisense,3.66%
Cells,Gene Expression,,,Confidently mapped to exonic regions,65.79%
Cells,Gene Expression,,,Confidently mapped to genome,92.83%
Cells,Gene Expression,,,Confidently mapped to intergenic regions,4.66%
Cells,Gene Expression,,,Confidently mapped to intronic regions,22.39%
Cells,Gene Expression,,,Confidently mapped to transcriptome,84.16%
Cells,Gene Expression,,,Mapped to genome,96.27%
Cells,Gene Expression,,,Median UMI counts per cell,"6,685"
Cells,Gene Expression,,,Median genes per cell,"2,468"
Cells,Gene Expression,,,Median reads per cell,"26,198"
Cells,Gene Expression,,,Number of reads from cells called from this sample,"61,845,367"
Cells,Gene Expression,,,Total genes detected,"20,942"
Cells,Gene Expression,Physical library ID,SC_GE,Cell-associated barcodes identified as multiplets,"1,222 (8.31%)"
Cells,Gene Expression,Physical library ID,SC_GE,Cell-associated barcodes not assigned any CMOs,"4,359 (29.65%)"
Cells,Gene Expression,Physical library ID,SC_GE,Cells assigned to other samples,"7,553 (51.37%)"
Cells,Gene Expression,Physical library ID,SC_GE,Cells assigned to this sample,"1,569 (10.67%)"
Library,Gene Expression,Fastq ID,SC_GE,Number of reads,"771,660,546"
Library,Gene Expression,Fastq ID,SC_GE,Number of short reads skipped,0
Library,Gene Expression,Fastq ID,SC_GE,Q30 RNA read,92.8%
Library,Gene Expression,Fastq ID,SC_GE,Q30 UMI,95.7%
Library,Gene Expression,Fastq ID,SC_GE,Q30 barcodes,96.2%
Library,Gene Expression,Physical library ID,SC_GE,Cell-associated barcodes identified as multiplets,"1,222 (8.31%)"
Library,Gene Expression,Physical library ID,SC_GE,Cell-associated barcodes not assigned any CMOs,"4,359 (29.65%)"
Library,Gene Expression,Physical library ID,SC_GE,Cells assigned to a sample,"9,122 (62.04%)"
Library,Gene Expression,Physical library ID,SC_GE,Confidently mapped antisense,3.41%
Library,Gene Expression,Physical library ID,SC_GE,Confidently mapped reads in cells,96.92%
Library,Gene Expression,Physical library ID,SC_GE,Confidently mapped to exonic regions,68.96%
Library,Gene Expression,Physical library ID,SC_GE,Confidently mapped to genome,91.86%
Library,Gene Expression,Physical library ID,SC_GE,Confidently mapped to intergenic regions,4.04%
Library,Gene Expression,Physical library ID,SC_GE,Confidently mapped to intronic regions,18.87%
Library,Gene Expression,Physical library ID,SC_GE,Confidently mapped to transcriptome,84.09%
Library,Gene Expression,Physical library ID,SC_GE,Estimated number of cells,"14,703"
Library,Gene Expression,Physical library ID,SC_GE,Mapped to genome,95.29%
Library,Gene Expression,Physical library ID,SC_GE,Mean reads per cell,"52,483"
Library,Gene Expression,Physical library ID,SC_GE,Number of reads,"771,660,546"
Library,Gene Expression,Physical library ID,SC_GE,Number of reads in the library,"771,660,546"
Library,Gene Expression,Physical library ID,SC_GE,Sequencing saturation,69.03%
Library,Gene Expression,Physical library ID,SC_GE,Valid UMIs,99.93%
Library,Gene Expression,Physical library ID,SC_GE,Valid barcodes,97.62%
Library,Multiplexing Capture,,,Cell-associated barcodes identified as multiplets,"1,222 (8.31%)"
Library,Multiplexing Capture,,,Cells assigned to a sample,"9,122 (62.04%)"
Library,Multiplexing Capture,,,Estimated number of cell-associated barcodes,"14,703"
Library,Multiplexing Capture,,,Median CMO UMIs per cell,"2,987"
Library,Multiplexing Capture,,,Number of samples assigned at least one cell,5
Library,Multiplexing Capture,,,Singlet capture ratio,0.69
Library,Multiplexing Capture,CMO Name,CMO301,CMO signal-to-noise ratio,4.25
Library,Multiplexing Capture,CMO Name,CMO301,Cells assigned to CMO,17.20%Library,Multiplexing Capture,CMO Name,CMO301,Fraction reads in cell-associated barcodes,72.63%
Library,Multiplexing Capture,CMO Name,CMO302,CMO signal-to-noise ratio,2.70
Library,Multiplexing Capture,CMO Name,CMO302,Cells assigned to CMO,13.24%
Library,Multiplexing Capture,CMO Name,CMO302,Fraction reads in cell-associated barcodes,49.82%
Library,Multiplexing Capture,CMO Name,CMO303,CMO signal-to-noise ratio,3.50
Library,Multiplexing Capture,CMO Name,CMO303,Cells assigned to CMO,27.26%
Library,Multiplexing Capture,CMO Name,CMO303,Fraction reads in cell-associated barcodes,71.93%
Library,Multiplexing Capture,CMO Name,CMO304,CMO signal-to-noise ratio,4.60
Library,Multiplexing Capture,CMO Name,CMO304,Cells assigned to CMO,14.42%
Library,Multiplexing Capture,CMO Name,CMO304,Fraction reads in cell-associated barcodes,66.88%
Library,Multiplexing Capture,CMO Name,CMO305,CMO signal-to-noise ratio,3.34
Library,Multiplexing Capture,CMO Name,CMO305,Cells assigned to CMO,27.88%
Library,Multiplexing Capture,CMO Name,CMO305,Fraction reads in cell-associated barcodes,73.27%
Library,Multiplexing Capture,Fastq ID,SC_CML,Number of reads,"204,261,166"
Library,Multiplexing Capture,Fastq ID,SC_CML,Number of short reads skipped,0
Library,Multiplexing Capture,Fastq ID,SC_CML,Q30 RNA read,96.9%
Library,Multiplexing Capture,Fastq ID,SC_CML,Q30 UMI,96.4%
Library,Multiplexing Capture,Fastq ID,SC_CML,Q30 barcodes,96.3%
Library,Multiplexing Capture,Physical library ID,SC_CML,Cell-associated barcodes identified as multiplets,"1,222 (8.31%)"
Library,Multiplexing Capture,Physical library ID,SC_CML,Cell-associated barcodes not assigned any CMOs,"4,359 (29.65%)"
Library,Multiplexing Capture,Physical library ID,SC_CML,Cells assigned to a sample,"9,122 (62.04%)"
Library,Multiplexing Capture,Physical library ID,SC_CML,Estimated number of cell-associated barcodes,"14,703"
Library,Multiplexing Capture,Physical library ID,SC_CML,Fraction CMO reads,97.66%
Library,Multiplexing Capture,Physical library ID,SC_CML,Fraction CMO reads usable,65.16%
Library,Multiplexing Capture,Physical library ID,SC_CML,Fraction reads from multiplets,0.00%
Library,Multiplexing Capture,Physical library ID,SC_CML,Fraction reads in cell-associated barcodes,67.12%
Library,Multiplexing Capture,Physical library ID,SC_CML,Fraction unrecognized CMO,2.34%
Library,Multiplexing Capture,Physical library ID,SC_CML,Mean reads per cell-associated barcode,"13,892"
Library,Multiplexing Capture,Physical library ID,SC_CML,Median CMO UMIs per cell-associated barcode,"2,814"
Library,Multiplexing Capture,Physical library ID,SC_CML,Number of reads,"204,261,166"
Library,Multiplexing Capture,Physical library ID,SC_CML,Samples assigned at least one cell,5
Library,Multiplexing Capture,Physical library ID,SC_CML,Sequencing saturation,29.48%
Library,Multiplexing Capture,Physical library ID,SC_CML,Valid UMIs,99.99%
Library,Multiplexing Capture,Physical library ID,SC_CML,Valid barcodes,99.17%
"""
CELLPLEX_METRICS_SUMMARY_CSP_7_1_0 = """Category,Library Type,Grouped By,Group Name,Metric Name,Metric Value
Cells,Antibody Capture,,,Cells,"1,582"
Cells,Antibody Capture,,,Mean antibody reads usable per cell,"4,527"
Cells,Antibody Capture,,,Median UMI counts per cell,"2,025"
Cells,Antibody Capture,Physical library ID,SS_CSP,Cell-associated barcodes identified as multiplets,"4,167 (20.38%)"
Cells,Antibody Capture,Physical library ID,SS_CSP,Cell-associated barcodes not assigned any CMOs,"1,728 (8.45%)"
Cells,Antibody Capture,Physical library ID,SS_CSP,Cells assigned to other samples,"12,966 (63.43%)"
Cells,Antibody Capture,Physical library ID,SS_CSP,Cells assigned to this sample,"1,582 (7.74%)"
Cells,Gene Expression,,,Cells,"1,582"
Cells,Gene Expression,,,Confidently mapped antisense,8.82%
Cells,Gene Expression,,,Confidently mapped to exonic regions,66.29%
Cells,Gene Expression,,,Confidently mapped to genome,93.98%
Cells,Gene Expression,,,Confidently mapped to intergenic regions,3.99%
Cells,Gene Expression,,,Confidently mapped to intronic regions,23.69%
Cells,Gene Expression,,,Confidently mapped to transcriptome,80.52%
Cells,Gene Expression,,,Mapped to genome,97.34%
Cells,Gene Expression,,,Median UMI counts per cell,"5,789"
Cells,Gene Expression,,,Median genes per cell,"2,042"
Cells,Gene Expression,,,Median reads per cell,"18,376"
Cells,Gene Expression,,,Number of reads from cells called from this sample,"32,611,101"
Cells,Gene Expression,,,Total genes detected,"23,199"
Cells,Gene Expression,Physical library ID,SS_GEX,Cell-associated barcodes identified as multiplets,"4,167 (20.38%)"
Cells,Gene Expression,Physical library ID,SS_GEX,Cell-associated barcodes not assigned any CMOs,"1,728 (8.45%)"
Cells,Gene Expression,Physical library ID,SS_GEX,Cells assigned to other samples,"12,966 (63.43%)"
Cells,Gene Expression,Physical library ID,SS_GEX,Cells assigned to this sample,"1,582 (7.74%)"
Library,Antibody Capture,Fastq ID,SS_CSP,Number of reads,"286,620,444"
Library,Antibody Capture,Fastq ID,SS_CSP,Number of short reads skipped,0
Library,Antibody Capture,Fastq ID,SS_CSP,Q30 RNA read,96.0%
Library,Antibody Capture,Fastq ID,SS_CSP,Q30 UMI,97.3%
Library,Antibody Capture,Fastq ID,SS_CSP,Q30 barcodes,97.2%
Library,Antibody Capture,Physical library ID,SS_CSP,Cell-associated barcodes identified as multiplets,"4,167 (20.38%)"
Library,Antibody Capture,Physical library ID,SS_CSP,Cell-associated barcodes not assigned any CMOs,"1,728 (8.45%)"
Library,Antibody Capture,Physical library ID,SS_CSP,Cells assigned to a sample,"14,548 (71.16%)"
Library,Antibody Capture,Physical library ID,SS_CSP,Estimated number of cells,"20,443"
Library,Antibody Capture,Physical library ID,SS_CSP,Fraction antibody reads,44.74%
Library,Antibody Capture,Physical library ID,SS_CSP,Fraction antibody reads in aggregate barcodes,2.01%
Library,Antibody Capture,Physical library ID,SS_CSP,Fraction antibody reads usable,29.02%
Library,Antibody Capture,Physical library ID,SS_CSP,Fraction reads in cells,65.34%
Library,Antibody Capture,Physical library ID,SS_CSP,Fraction unrecognized antibody,55.26%
Library,Antibody Capture,Physical library ID,SS_CSP,Mean reads per cell,"14,020"
Library,Antibody Capture,Physical library ID,SS_CSP,Number of reads,"286,620,444"
Library,Antibody Capture,Physical library ID,SS_CSP,Sequencing saturation,37.16%
Library,Antibody Capture,Physical library ID,SS_CSP,Valid UMIs,100.00%
Library,Antibody Capture,Physical library ID,SS_CSP,Valid barcodes,99.12%
Library,Gene Expression,Fastq ID,SS_GEX,Number of reads,"510,571,324"
Library,Gene Expression,Fastq ID,SS_GEX,Number of short reads skipped,0
Library,Gene Expression,Fastq ID,SS_GEX,Q30 RNA read,93.3%
Library,Gene Expression,Fastq ID,SS_GEX,Q30 UMI,95.8%
Library,Gene Expression,Fastq ID,SS_GEX,Q30 barcodes,96.5%
Library,Gene Expression,Physical library ID,SS_GEX,Cell-associated barcodes identified as multiplets,"4,167 (20.38%)"
Library,Gene Expression,Physical library ID,SS_GEX,Cell-associated barcodes not assigned any CMOs,"1,728 (8.45%)"
Library,Gene Expression,Physical library ID,SS_GEX,Cells assigned to a sample,"14,548 (71.16%)"
Library,Gene Expression,Physical library ID,SS_GEX,Confidently mapped antisense,9.40%
Library,Gene Expression,Physical library ID,SS_GEX,Confidently mapped reads in cells,94.20%
Library,Gene Expression,Physical library ID,SS_GEX,Confidently mapped to exonic regions,65.20%
Library,Gene Expression,Physical library ID,SS_GEX,Confidently mapped to genome,92.78%
Library,Gene Expression,Physical library ID,SS_GEX,Confidently mapped to intergenic regions,3.60%
Library,Gene Expression,Physical library ID,SS_GEX,Confidently mapped to intronic regions,23.99%
Library,Gene Expression,Physical library ID,SS_GEX,Confidently mapped to transcriptome,79.14%
Library,Gene Expression,Physical library ID,SS_GEX,Estimated number of cells,"20,443"
Library,Gene Expression,Physical library ID,SS_GEX,Mapped to genome,96.67%
Library,Gene Expression,Physical library ID,SS_GEX,Mean reads per cell,"24,975"
Library,Gene Expression,Physical library ID,SS_GEX,Number of reads,"510,571,324"
Library,Gene Expression,Physical library ID,SS_GEX,Number of reads in the library,"510,571,324"
Library,Gene Expression,Physical library ID,SS_GEX,Sequencing saturation,59.89%
Library,Gene Expression,Physical library ID,SS_GEX,Valid UMIs,99.94%
Library,Gene Expression,Physical library ID,SS_GEX,Valid barcodes,97.81%
Library,Multiplexing Capture,,,Cell-associated barcodes identified as multiplets,"4,167 (20.38%)"
Library,Multiplexing Capture,,,Cells assigned to a sample,"14,548 (71.16%)"
Library,Multiplexing Capture,,,Estimated number of cell-associated barcodes,"20,443"
Library,Multiplexing Capture,,,Median CMO UMIs per cell,"1,865"
Library,Multiplexing Capture,,,Number of samples assigned at least one cell,9
Library,Multiplexing Capture,,,Singlet capture ratio,0.86
Library,Multiplexing Capture,CMO Name,CMO301,CMO signal-to-noise ratio,5.97
Library,Multiplexing Capture,CMO Name,CMO301,Cells assigned to CMO,10.87%
Library,Multiplexing Capture,CMO Name,CMO301,Fraction reads in cell-associated barcodes,66.78%
Library,Multiplexing Capture,CMO Name,CMO302,CMO signal-to-noise ratio,4.87
Library,Multiplexing Capture,CMO Name,CMO302,Cells assigned to CMO,10.21%
Library,Multiplexing Capture,CMO Name,CMO302,Fraction reads in cell-associated barcodes,63.16%
Library,Multiplexing Capture,CMO Name,CMO303,CMO signal-to-noise ratio,4.98
Library,Multiplexing Capture,CMO Name,CMO303,Cells assigned to CMO,12.87%
Library,Multiplexing Capture,CMO Name,CMO303,Fraction reads in cell-associated barcodes,73.73%
Library,Multiplexing Capture,CMO Name,CMO304,CMO signal-to-noise ratio,6.13
Library,Multiplexing Capture,CMO Name,CMO304,Cells assigned to CMO,11.60%
Library,Multiplexing Capture,CMO Name,CMO304,Fraction reads in cell-associated barcodes,63.48%
Library,Multiplexing Capture,CMO Name,CMO305,CMO signal-to-noise ratio,4.55
Library,Multiplexing Capture,CMO Name,CMO305,Cells assigned to CMO,9.99%
Library,Multiplexing Capture,CMO Name,CMO305,Fraction reads in cell-associated barcodes,65.07%
Library,Multiplexing Capture,CMO Name,CMO306,CMO signal-to-noise ratio,6.22
Library,Multiplexing Capture,CMO Name,CMO306,Cells assigned to CMO,12.41%
Library,Multiplexing Capture,CMO Name,CMO306,Fraction reads in cell-associated barcodes,64.30%
Library,Multiplexing Capture,CMO Name,CMO307,CMO signal-to-noise ratio,5.84
Library,Multiplexing Capture,CMO Name,CMO307,Cells assigned to CMO,9.61%
Library,Multiplexing Capture,CMO Name,CMO307,Fraction reads in cell-associated barcodes,59.05%
Library,Multiplexing Capture,CMO Name,CMO308,CMO signal-to-noise ratio,4.34
Library,Multiplexing Capture,CMO Name,CMO308,Cells assigned to CMO,10.17%
Library,Multiplexing Capture,CMO Name,CMO308,Fraction reads in cell-associated barcodes,56.30%
Library,Multiplexing Capture,CMO Name,CMO309,CMO signal-to-noise ratio,5.54
Library,Multiplexing Capture,CMO Name,CMO309,Cells assigned to CMO,12.27%
Library,Multiplexing Capture,CMO Name,CMO309,Fraction reads in cell-associated barcodes,60.63%
Library,Multiplexing Capture,Fastq ID,SS_CML,Number of reads,"123,420,555"
Library,Multiplexing Capture,Fastq ID,SS_CML,Number of short reads skipped,0
Library,Multiplexing Capture,Fastq ID,SS_CML,Q30 RNA read,97.0%
Library,Multiplexing Capture,Fastq ID,SS_CML,Q30 UMI,97.3%
Library,Multiplexing Capture,Fastq ID,SS_CML,Q30 barcodes,97.1%
Library,Multiplexing Capture,Physical library ID,SS_CML,Cell-associated barcodes identified as multiplets,"4,167 (20.38%)"
Library,Multiplexing Capture,Physical library ID,SS_CML,Cell-associated barcodes not assigned any CMOs,"1,728 (8.45%)"
Library,Multiplexing Capture,Physical library ID,SS_CML,Cells assigned to a sample,"14,548 (71.16%)"
Library,Multiplexing Capture,Physical library ID,SS_CML,Estimated number of cell-associated barcodes,"20,443"
Library,Multiplexing Capture,Physical library ID,SS_CML,Fraction CMO reads,98.10%
Library,Multiplexing Capture,Physical library ID,SS_CML,Fraction CMO reads usable,62.28%
Library,Multiplexing Capture,Physical library ID,SS_CML,Fraction reads from multiplets,0.00%
Library,Multiplexing Capture,Physical library ID,SS_CML,Fraction reads in cell-associated barcodes,63.75%
Library,Multiplexing Capture,Physical library ID,SS_CML,Fraction unrecognized CMO,1.90%
Library,Multiplexing Capture,Physical library ID,SS_CML,Mean reads per cell-associated barcode,"6,037"
Library,Multiplexing Capture,Physical library ID,SS_CML,Median CMO UMIs per cell-associated barcode,"2,087"
Library,Multiplexing Capture,Physical library ID,SS_CML,Number of reads,"123,420,555"
Library,Multiplexing Capture,Physical library ID,SS_CML,Samples assigned at least one cell,9
Library,Multiplexing Capture,Physical library ID,SS_CML,Sequencing saturation,14.40%
Library,Multiplexing Capture,Physical library ID,SS_CML,Valid UMIs,100.00%
Library,Multiplexing Capture,Physical library ID,SS_CML,Valid barcodes,99.51%
"""

CELLPLEX_METRICS_SUMMARY_8_0_0 = """Category,Library Type,Grouped By,Group Name,Metric Name,Metric Value
Cells,Antibody Capture,,,Cells,"1,582"
Cells,Antibody Capture,,,Fraction antibody reads,44.13%
Cells,Antibody Capture,,,Fraction antibody reads in aggregate barcodes,---
Cells,Antibody Capture,,,Mean antibody reads usable per cell,"4,527"
Cells,Antibody Capture,,,Median UMI counts per cell,"2,025"
Cells,Antibody Capture,,,Number of reads from cells associated with this sample,"16,235,484"
Cells,Antibody Capture,Physical library ID,SS_CSP,Cell-associated barcodes identified as multiplets,"4,167 (20.38%)"
Cells,Antibody Capture,Physical library ID,SS_CSP,Cell-associated barcodes not assigned any CMOs,"1,728 (8.45%)"
Cells,Antibody Capture,Physical library ID,SS_CSP,Cells assigned to other samples,"12,966 (63.43%)"
Cells,Antibody Capture,Physical library ID,SS_CSP,Cells assigned to this sample,"1,582 (7.74%)"
Cells,Gene Expression,,,Cells,"1,582"
Cells,Gene Expression,,,Confidently mapped antisense,8.82%
Cells,Gene Expression,,,Confidently mapped to exonic regions,66.29%
Cells,Gene Expression,,,Confidently mapped to genome,93.98%
Cells,Gene Expression,,,Confidently mapped to intergenic regions,3.99%
Cells,Gene Expression,,,Confidently mapped to intronic regions,23.69%
Cells,Gene Expression,,,Confidently mapped to transcriptome,80.52%
Cells,Gene Expression,,,Mapped to genome,97.34%
Cells,Gene Expression,,,Mean reads per cell,"20,614"
Cells,Gene Expression,,,Median UMI counts per cell,"5,789"
Cells,Gene Expression,,,Median genes per cell,"2,042"
Cells,Gene Expression,,,Number of reads from cells called from this sample,"32,611,101"
Cells,Gene Expression,,,Total genes detected,"23,199"
Cells,Gene Expression,Physical library ID,SS_GEX,Cell-associated barcodes identified as multiplets,"4,167 (20.38%)"
Cells,Gene Expression,Physical library ID,SS_GEX,Cell-associated barcodes not assigned any CMOs,"1,728 (8.45%)"
Cells,Gene Expression,Physical library ID,SS_GEX,Cells assigned to other samples,"12,966 (63.43%)"
Cells,Gene Expression,Physical library ID,SS_GEX,Cells assigned to this sample,"1,582 (7.74%)"
Library,Antibody Capture,Fastq ID,SS_CSP,Number of reads,"286,620,444"
Library,Antibody Capture,Fastq ID,SS_CSP,Number of short reads skipped,0
Library,Antibody Capture,Fastq ID,SS_CSP,Q30 RNA read,96.0%
Library,Antibody Capture,Fastq ID,SS_CSP,Q30 UMI,97.3%
Library,Antibody Capture,Fastq ID,SS_CSP,Q30 barcodes,97.2%
Library,Antibody Capture,Physical library ID,SS_CSP,Antibody reads in cells,65.34%
Library,Antibody Capture,Physical library ID,SS_CSP,Cell-associated barcodes identified as multiplets,"4,167 (20.38%)"
Library,Antibody Capture,Physical library ID,SS_CSP,Cell-associated barcodes not assigned any CMOs,"1,728 (8.45%)"
Library,Antibody Capture,Physical library ID,SS_CSP,Cells assigned to a sample,"14,548 (71.16%)"
Library,Antibody Capture,Physical library ID,SS_CSP,Estimated number of cells,"20,443"
Library,Antibody Capture,Physical library ID,SS_CSP,Fraction antibody reads,44.74%
Library,Antibody Capture,Physical library ID,SS_CSP,Fraction antibody reads in aggregate barcodes,2.01%
Library,Antibody Capture,Physical library ID,SS_CSP,Fraction antibody reads usable,29.02%
Library,Antibody Capture,Physical library ID,SS_CSP,Mean reads per cell,"14,020"
Library,Antibody Capture,Physical library ID,SS_CSP,Number of reads,"286,620,444"
Library,Antibody Capture,Physical library ID,SS_CSP,Number of reads in the library,"286,620,444"
Library,Antibody Capture,Physical library ID,SS_CSP,Sequencing saturation,37.16%
Library,Antibody Capture,Physical library ID,SS_CSP,Valid UMIs,100.00%
Library,Antibody Capture,Physical library ID,SS_CSP,Valid barcodes,99.12%
Library,Gene Expression,Fastq ID,SS_GEX,Number of reads,"510,571,324"
Library,Gene Expression,Fastq ID,SS_GEX,Number of short reads skipped,0
Library,Gene Expression,Fastq ID,SS_GEX,Q30 RNA read,93.3%
Library,Gene Expression,Fastq ID,SS_GEX,Q30 UMI,95.8%
Library,Gene Expression,Fastq ID,SS_GEX,Q30 barcodes,96.5%
Library,Gene Expression,Physical library ID,SS_GEX,Cell-associated barcodes identified as multiplets,"4,167 (20.38%)"
Library,Gene Expression,Physical library ID,SS_GEX,Cell-associated barcodes not assigned any CMOs,"1,728 (8.45%)"
Library,Gene Expression,Physical library ID,SS_GEX,Cells assigned to a sample,"14,548 (71.16%)"
Library,Gene Expression,Physical library ID,SS_GEX,Confidently mapped antisense,9.40%
Library,Gene Expression,Physical library ID,SS_GEX,Confidently mapped reads in cells,94.20%
Library,Gene Expression,Physical library ID,SS_GEX,Confidently mapped to exonic regions,65.20%
Library,Gene Expression,Physical library ID,SS_GEX,Confidently mapped to genome,92.78%
Library,Gene Expression,Physical library ID,SS_GEX,Confidently mapped to intergenic regions,3.60%
Library,Gene Expression,Physical library ID,SS_GEX,Confidently mapped to intronic regions,23.99%
Library,Gene Expression,Physical library ID,SS_GEX,Confidently mapped to transcriptome,79.14%
Library,Gene Expression,Physical library ID,SS_GEX,Estimated number of cells,"20,443"
Library,Gene Expression,Physical library ID,SS_GEX,Mapped to genome,96.67%
Library,Gene Expression,Physical library ID,SS_GEX,Mean reads per cell,"24,975"
Library,Gene Expression,Physical library ID,SS_GEX,Number of reads,"510,571,324"
Library,Gene Expression,Physical library ID,SS_GEX,Number of reads in the library,"510,571,324"
Library,Gene Expression,Physical library ID,SS_GEX,Sequencing saturation,59.89%
Library,Gene Expression,Physical library ID,SS_GEX,Valid UMIs,99.94%
Library,Gene Expression,Physical library ID,SS_GEX,Valid barcodes,97.81%
Library,Multiplexing Capture,,,Cell-associated barcodes identified as multiplets,"4,167 (20.38%)"
Library,Multiplexing Capture,,,Cells assigned to a sample,"14,548 (71.16%)"
Library,Multiplexing Capture,,,Estimated number of cell-associated barcodes,"20,443"
Library,Multiplexing Capture,,,Median CMO UMIs per cell,"1,865"
Library,Multiplexing Capture,,,Number of samples assigned at least one cell,9
Library,Multiplexing Capture,,,Singlet capture ratio,0.86
Library,Multiplexing Capture,CMO Name,CMO301,CMO signal-to-noise ratio,5.97
Library,Multiplexing Capture,CMO Name,CMO301,Cells assigned to CMO,10.87%
Library,Multiplexing Capture,CMO Name,CMO301,Fraction reads in cell-associated barcodes,66.78%
Library,Multiplexing Capture,CMO Name,CMO302,CMO signal-to-noise ratio,4.87
Library,Multiplexing Capture,CMO Name,CMO302,Cells assigned to CMO,10.21%
Library,Multiplexing Capture,CMO Name,CMO302,Fraction reads in cell-associated barcodes,63.16%
Library,Multiplexing Capture,CMO Name,CMO303,CMO signal-to-noise ratio,4.98
Library,Multiplexing Capture,CMO Name,CMO303,Cells assigned to CMO,12.87%
Library,Multiplexing Capture,CMO Name,CMO303,Fraction reads in cell-associated barcodes,73.73%
Library,Multiplexing Capture,CMO Name,CMO304,CMO signal-to-noise ratio,6.13
Library,Multiplexing Capture,CMO Name,CMO304,Cells assigned to CMO,11.60%
Library,Multiplexing Capture,CMO Name,CMO304,Fraction reads in cell-associated barcodes,63.48%
Library,Multiplexing Capture,CMO Name,CMO305,CMO signal-to-noise ratio,4.55
Library,Multiplexing Capture,CMO Name,CMO305,Cells assigned to CMO,9.99%
Library,Multiplexing Capture,CMO Name,CMO305,Fraction reads in cell-associated barcodes,65.07%
Library,Multiplexing Capture,CMO Name,CMO306,CMO signal-to-noise ratio,6.22
Library,Multiplexing Capture,CMO Name,CMO306,Cells assigned to CMO,12.41%
Library,Multiplexing Capture,CMO Name,CMO306,Fraction reads in cell-associated barcodes,64.30%
Library,Multiplexing Capture,CMO Name,CMO307,CMO signal-to-noise ratio,5.84
Library,Multiplexing Capture,CMO Name,CMO307,Cells assigned to CMO,9.61%
Library,Multiplexing Capture,CMO Name,CMO307,Fraction reads in cell-associated barcodes,59.05%
Library,Multiplexing Capture,CMO Name,CMO308,CMO signal-to-noise ratio,4.34
Library,Multiplexing Capture,CMO Name,CMO308,Cells assigned to CMO,10.17%
Library,Multiplexing Capture,CMO Name,CMO308,Fraction reads in cell-associated barcodes,56.30%
Library,Multiplexing Capture,CMO Name,CMO309,CMO signal-to-noise ratio,5.54
Library,Multiplexing Capture,CMO Name,CMO309,Cells assigned to CMO,12.27%
Library,Multiplexing Capture,CMO Name,CMO309,Fraction reads in cell-associated barcodes,60.63%
Library,Multiplexing Capture,Fastq ID,SS_CML,Number of reads,"123,420,555"
Library,Multiplexing Capture,Fastq ID,SS_CML,Number of short reads skipped,0
Library,Multiplexing Capture,Fastq ID,SS_CML,Q30 RNA read,97.0%
Library,Multiplexing Capture,Fastq ID,SS_CML,Q30 UMI,97.3%
Library,Multiplexing Capture,Fastq ID,SS_CML,Q30 barcodes,97.1%
Library,Multiplexing Capture,Physical library ID,SS_CML,Cell-associated barcodes identified as multiplets,"4,167 (20.38%)"
Library,Multiplexing Capture,Physical library ID,SS_CML,Cell-associated barcodes not assigned any CMOs,"1,728 (8.45%)"
Library,Multiplexing Capture,Physical library ID,SS_CML,Cells assigned to a sample,"14,548 (71.16%)"
Library,Multiplexing Capture,Physical library ID,SS_CML,Estimated number of cell-associated barcodes,"20,443"
Library,Multiplexing Capture,Physical library ID,SS_CML,Fraction CMO reads,98.10%
Library,Multiplexing Capture,Physical library ID,SS_CML,Fraction CMO reads usable,62.28%
Library,Multiplexing Capture,Physical library ID,SS_CML,Fraction reads from multiplets,21.90%
Library,Multiplexing Capture,Physical library ID,SS_CML,Fraction reads in cell-associated barcodes,63.75%
Library,Multiplexing Capture,Physical library ID,SS_CML,Fraction unrecognized CMO,1.90%
Library,Multiplexing Capture,Physical library ID,SS_CML,Mean reads per cell-associated barcode,"6,037"
Library,Multiplexing Capture,Physical library ID,SS_CML,Median CMO UMIs per cell-associated barcode,"2,087"
Library,Multiplexing Capture,Physical library ID,SS_CML,Number of reads,"123,420,555"
Library,Multiplexing Capture,Physical library ID,SS_CML,Samples assigned at least one cell,9
Library,Multiplexing Capture,Physical library ID,SS_CML,Sequencing saturation,14.40%
Library,Multiplexing Capture,Physical library ID,SS_CML,Valid UMIs,100.00%
Library,Multiplexing Capture,Physical library ID,SS_CML,Valid barcodes,99.51%
"""

FLEX_METRICS_SUMMARY_8_0_0 = """Category,Library Type,Grouped By,Group Name,Metric Name,Metric Value
Cells,Gene Expression,,,Cells,"2,302"
Cells,Gene Expression,,,Confidently mapped reads in cells,93.06%
Cells,Gene Expression,,,Estimated UMIs from genomic DNA,0.14%
Cells,Gene Expression,,,Estimated UMIs from genomic DNA per unspliced probe,1
Cells,Gene Expression,,,Mean reads per cell,"27,485"
Cells,Gene Expression,,,Median UMI counts per cell,"4,202"
Cells,Gene Expression,,,Median genes per cell,"2,465"
Cells,Gene Expression,,,Number of reads from cells called from this sample,"63,271,140"
Cells,Gene Expression,,,Reads confidently mapped to filtered probe set,96.97%
Cells,Gene Expression,,,Reads confidently mapped to probe set,97.30%
Cells,Gene Expression,,,Reads half-mapped to probe set,0.18%
Cells,Gene Expression,,,Reads mapped to probe set,99.53%
Cells,Gene Expression,,,Reads split-mapped to probe set,2.04%
Cells,Gene Expression,,,Total genes detected,"15,937"
Cells,Gene Expression,Physical library ID,JMF1,Cells detected in other samples,"41,675 (94.77%)"
Cells,Gene Expression,Physical library ID,JMF1,Cells detected in this sample,"2,302 (5.23%)"
Library,Gene Expression,,,Estimated UMIs from genomic DNA,0.01%
Library,Gene Expression,,,Estimated UMIs from genomic DNA per unspliced probe,2
Library,Gene Expression,Fastq ID,JMF1,Number of reads,"1,947,449,723"
Library,Gene Expression,Fastq ID,JMF1,Number of short reads skipped,0
Library,Gene Expression,Fastq ID,JMF1,Q30 GEM barcodes,97.7%
Library,Gene Expression,Fastq ID,JMF1,Q30 RNA read,97.3%
Library,Gene Expression,Fastq ID,JMF1,Q30 UMI,97.4%
Library,Gene Expression,Fastq ID,JMF1,Q30 barcodes,97.2%
Library,Gene Expression,Fastq ID,JMF1,Q30 probe barcodes,96.3%
Library,Gene Expression,Physical library ID,JMF1,Confidently mapped reads in cells,92.87%
Library,Gene Expression,Physical library ID,JMF1,Estimated number of cells,"43,977"
Library,Gene Expression,Physical library ID,JMF1,Fraction of initial cell barcodes passing high occupancy GEM filtering,99.46%
Library,Gene Expression,Physical library ID,JMF1,Mean reads per cell,"44,283"
Library,Gene Expression,Physical library ID,JMF1,Number of reads,"1,947,449,723"
Library,Gene Expression,Physical library ID,JMF1,Number of reads in the library,"1,947,449,723"
Library,Gene Expression,Physical library ID,JMF1,Reads confidently mapped to filtered probe set,94.91%
Library,Gene Expression,Physical library ID,JMF1,Reads confidently mapped to probe set,95.23%
Library,Gene Expression,Physical library ID,JMF1,Reads half-mapped to probe set,1.29%
Library,Gene Expression,Physical library ID,JMF1,Reads mapped to probe set,98.66%
Library,Gene Expression,Physical library ID,JMF1,Reads split-mapped to probe set,2.14%
Library,Gene Expression,Physical library ID,JMF1,Sequencing saturation,64.85%
Library,Gene Expression,Physical library ID,JMF1,Valid GEM barcodes,98.81%
Library,Gene Expression,Physical library ID,JMF1,Valid UMIs,100.00%
Library,Gene Expression,Physical library ID,JMF1,Valid barcodes,96.04%
Library,Gene Expression,Physical library ID,JMF1,Valid probe barcodes,96.94%
Library,Gene Expression,Probe barcode ID,BC001,Cells per probe barcode,"2,302 (5.23%)"
Library,Gene Expression,Probe barcode ID,BC001,Sample ID,JMF1
Library,Gene Expression,Probe barcode ID,BC001,UMIs per probe barcode,"23,357,091 (3.67%)"
Library,Gene Expression,Probe barcode ID,BC002,Cells per probe barcode,"2,308 (5.25%)"
Library,Gene Expression,Probe barcode ID,BC002,Sample ID,JMF2
Library,Gene Expression,Probe barcode ID,BC002,UMIs per probe barcode,"22,251,203 (3.50%)"
Library,Gene Expression,Probe barcode ID,BC003,Cells per probe barcode,"2,863 (6.51%)"
Library,Gene Expression,Probe barcode ID,BC003,Sample ID,JMF3
Library,Gene Expression,Probe barcode ID,BC003,UMIs per probe barcode,"29,304,274 (4.61%)"
Library,Gene Expression,Probe barcode ID,BC004,Cells per probe barcode,"2,250 (5.12%)"
Library,Gene Expression,Probe barcode ID,BC004,Sample ID,JMF4
Library,Gene Expression,Probe barcode ID,BC004,UMIs per probe barcode,"31,493,332 (4.95%)"
Library,Gene Expression,Probe barcode ID,BC005,Cells per probe barcode,"2,415 (5.49%)"
Library,Gene Expression,Probe barcode ID,BC005,Sample ID,JMF5
Library,Gene Expression,Probe barcode ID,BC005,UMIs per probe barcode,"35,011,126 (5.50%)"
Library,Gene Expression,Probe barcode ID,BC006,Cells per probe barcode,"3,193 (7.26%)"
Library,Gene Expression,Probe barcode ID,BC006,Sample ID,JMF6
Library,Gene Expression,Probe barcode ID,BC006,UMIs per probe barcode,"46,457,109 (7.30%)"
Library,Gene Expression,Probe barcode ID,BC007,Cells per probe barcode,"3,446 (7.84%)"
Library,Gene Expression,Probe barcode ID,BC007,Sample ID,JMF7
Library,Gene Expression,Probe barcode ID,BC007,UMIs per probe barcode,"43,765,251 (6.88%)"
Library,Gene Expression,Probe barcode ID,BC008,Cells per probe barcode,"3,244 (7.38%)"
Library,Gene Expression,Probe barcode ID,BC008,Sample ID,JMF8
Library,Gene Expression,Probe barcode ID,BC008,UMIs per probe barcode,"42,711,979 (6.72%)"
Library,Gene Expression,Probe barcode ID,BC009,Cells per probe barcode,"3,089 (7.02%)"
Library,Gene Expression,Probe barcode ID,BC009,Sample ID,JMF9
Library,Gene Expression,Probe barcode ID,BC009,UMIs per probe barcode,"60,515,493 (9.51%)"
Library,Gene Expression,Probe barcode ID,BC010,Cells per probe barcode,"2,276 (5.18%)"
Library,Gene Expression,Probe barcode ID,BC010,Sample ID,JMF10
Library,Gene Expression,Probe barcode ID,BC010,UMIs per probe barcode,"44,437,435 (6.99%)"
Library,Gene Expression,Probe barcode ID,BC011,Cells per probe barcode,"2,492 (5.67%)"
Library,Gene Expression,Probe barcode ID,BC011,Sample ID,JMF11
Library,Gene Expression,Probe barcode ID,BC011,UMIs per probe barcode,"47,811,412 (7.52%)"
Library,Gene Expression,Probe barcode ID,BC012,Cells per probe barcode,"2,254 (5.13%)"
Library,Gene Expression,Probe barcode ID,BC012,Sample ID,JMF12
Library,Gene Expression,Probe barcode ID,BC012,UMIs per probe barcode,"32,291,662 (5.08%)"
Library,Gene Expression,Probe barcode ID,BC013,Cells per probe barcode,"2,405 (5.47%)"
Library,Gene Expression,Probe barcode ID,BC013,Sample ID,JMF13
Library,Gene Expression,Probe barcode ID,BC013,UMIs per probe barcode,"37,113,794 (5.84%)"
Library,Gene Expression,Probe barcode ID,BC014,Cells per probe barcode,"2,219 (5.05%)"
Library,Gene Expression,Probe barcode ID,BC014,Sample ID,JMF14
Library,Gene Expression,Probe barcode ID,BC014,UMIs per probe barcode,"35,921,883 (5.65%)"
Library,Gene Expression,Probe barcode ID,BC015,Cells per probe barcode,"3,560 (8.10%)"
Library,Gene Expression,Probe barcode ID,BC015,Sample ID,JMF15
Library,Gene Expression,Probe barcode ID,BC015,UMIs per probe barcode,"51,216,165 (8.05%)"
Library,Gene Expression,Probe barcode ID,BC016,Cells per probe barcode,"3,661 (8.32%)"
Library,Gene Expression,Probe barcode ID,BC016,Sample ID,JMF16
Library,Gene Expression,Probe barcode ID,BC016,UMIs per probe barcode,"52,394,225 (8.24%)"
"""

CELLPLEX_METRICS_SUMMARY_9_0_0 = """Category,Library Type,Grouped By,Group Name,Metric Name,Metric Value
Cells,Gene Expression,,,Cells,"1,571"
Cells,Gene Expression,,,Confidently mapped antisense,3.7%
Cells,Gene Expression,,,Confidently mapped to exonic regions,65.8%
Cells,Gene Expression,,,Confidently mapped to genome,92.8%
Cells,Gene Expression,,,Confidently mapped to intergenic regions,4.7%
Cells,Gene Expression,,,Confidently mapped to intronic regions,22.4%
Cells,Gene Expression,,,Confidently mapped to transcriptome,84.2%
Cells,Gene Expression,,,Mapped to genome,96.3%
Cells,Gene Expression,,,Mean reads per cell,"39,445"
Cells,Gene Expression,,,Median UMI counts per cell,"6,692"
Cells,Gene Expression,,,Median genes per cell,"2,472"
Cells,Gene Expression,,,Number of reads from cells called from this sample,"61,968,623"
Cells,Gene Expression,,,Total genes detected,"20,944"
Cells,Gene Expression,Physical library ID,SC_GE,Cell-associated barcodes identified as multiplets,"1,218 (8.3%)"
Cells,Gene Expression,Physical library ID,SC_GE,Cell-associated barcodes not assigned any CMOs,"4,345 (29.6%)"
Cells,Gene Expression,Physical library ID,SC_GE,Cells assigned to other samples,"7,545 (51.4%)"
Cells,Gene Expression,Physical library ID,SC_GE,Cells assigned to this sample,"1,571 (10.7%)"
Library,Gene Expression,Fastq ID,SC_GE,Number of reads,"771,660,546"
Library,Gene Expression,Fastq ID,SC_GE,Number of short reads skipped,0
Library,Gene Expression,Fastq ID,SC_GE,Q30 RNA read,92.8%
Library,Gene Expression,Fastq ID,SC_GE,Q30 UMI,95.7%
Library,Gene Expression,Fastq ID,SC_GE,Q30 barcodes,96.2%
Library,Gene Expression,Physical library ID,SC_GE,Cell-associated barcodes identified as multiplets,"1,218 (8.3%)"
Library,Gene Expression,Physical library ID,SC_GE,Cell-associated barcodes not assigned any CMOs,"4,345 (29.6%)"
Library,Gene Expression,Physical library ID,SC_GE,Cells,"14,679"
Library,Gene Expression,Physical library ID,SC_GE,Cells assigned to a sample,"9,116 (62.1%)"
Library,Gene Expression,Physical library ID,SC_GE,Confidently mapped antisense,3.4%
Library,Gene Expression,Physical library ID,SC_GE,Confidently mapped reads in cells,96.9%
Library,Gene Expression,Physical library ID,SC_GE,Confidently mapped to exonic regions,69.0%
Library,Gene Expression,Physical library ID,SC_GE,Confidently mapped to genome,91.9%
Library,Gene Expression,Physical library ID,SC_GE,Confidently mapped to intergenic regions,4.0%
Library,Gene Expression,Physical library ID,SC_GE,Confidently mapped to intronic regions,18.9%
Library,Gene Expression,Physical library ID,SC_GE,Confidently mapped to transcriptome,84.1%
Library,Gene Expression,Physical library ID,SC_GE,Mapped to genome,95.3%
Library,Gene Expression,Physical library ID,SC_GE,Mean reads per cell,"52,569"
Library,Gene Expression,Physical library ID,SC_GE,Number of reads,"771,660,546"
Library,Gene Expression,Physical library ID,SC_GE,Number of reads in the library,"771,660,546"
Library,Gene Expression,Physical library ID,SC_GE,Sequencing saturation,69.0%
Library,Gene Expression,Physical library ID,SC_GE,Valid UMIs,99.9%
Library,Gene Expression,Physical library ID,SC_GE,Valid barcodes,97.6%
Library,Multiplexing Capture,,,Cell-associated barcodes identified as multiplets,"1,218 (8.3%)"
Library,Multiplexing Capture,,,Cell-associated barcodes not assigned any CMOs,"4,345 (29.6%)"
Library,Multiplexing Capture,,,Cells assigned to a sample,"9,116 (62.1%)"
Library,Multiplexing Capture,,,Estimated number of cell-associated barcodes,"14,679"
Library,Multiplexing Capture,,,Median CMO UMI Counts per cell assigned to a sample,"2,987"
Library,Multiplexing Capture,,,Number of samples assigned at least one cell,5
Library,Multiplexing Capture,,,Singlet capture ratio,0.69
Library,Multiplexing Capture,CMO Name,CMO301,CMO signal-to-noise ratio,4.25
Library,Multiplexing Capture,CMO Name,CMO301,Cells per CMO,"1,571 (17.2%)"
Library,Multiplexing Capture,CMO Name,CMO301,Fraction reads in cell-associated barcodes,72.4%
Library,Multiplexing Capture,CMO Name,CMO301,Sample ID,SC1
Library,Multiplexing Capture,CMO Name,CMO302,CMO signal-to-noise ratio,2.70
Library,Multiplexing Capture,CMO Name,CMO302,Cells per CMO,"1,204 (13.2%)"
Library,Multiplexing Capture,CMO Name,CMO302,Fraction reads in cell-associated barcodes,49.7%
Library,Multiplexing Capture,CMO Name,CMO302,Sample ID,SC2
Library,Multiplexing Capture,CMO Name,CMO303,CMO signal-to-noise ratio,3.50
Library,Multiplexing Capture,CMO Name,CMO303,Cells per CMO,"2,487 (27.3%)"
Library,Multiplexing Capture,CMO Name,CMO303,Fraction reads in cell-associated barcodes,71.9%
Library,Multiplexing Capture,CMO Name,CMO303,Sample ID,SC3
Library,Multiplexing Capture,CMO Name,CMO304,CMO signal-to-noise ratio,4.60
Library,Multiplexing Capture,CMO Name,CMO304,Cells per CMO,"1,313 (14.4%)"
Library,Multiplexing Capture,CMO Name,CMO304,Fraction reads in cell-associated barcodes,66.8%
Library,Multiplexing Capture,CMO Name,CMO304,Sample ID,SC4
Library,Multiplexing Capture,CMO Name,CMO305,CMO signal-to-noise ratio,3.34
Library,Multiplexing Capture,CMO Name,CMO305,Cells per CMO,"2,541 (27.9%)"
Library,Multiplexing Capture,CMO Name,CMO305,Fraction reads in cell-associated barcodes,73.0%
Library,Multiplexing Capture,CMO Name,CMO305,Sample ID,SC5
Library,Multiplexing Capture,Fastq ID,SC_CML,Number of reads,"204,261,166"
Library,Multiplexing Capture,Fastq ID,SC_CML,Number of short reads skipped,0
Library,Multiplexing Capture,Fastq ID,SC_CML,Q30 RNA read,96.9%
Library,Multiplexing Capture,Fastq ID,SC_CML,Q30 UMI,96.4%
Library,Multiplexing Capture,Fastq ID,SC_CML,Q30 barcodes,96.3%
Library,Multiplexing Capture,Physical library ID,SC_CML,Fraction CMO reads,97.7%
Library,Multiplexing Capture,Physical library ID,SC_CML,Fraction CMO reads usable,65.0%
Library,Multiplexing Capture,Physical library ID,SC_CML,Fraction reads from multiplets,16.0%
Library,Multiplexing Capture,Physical library ID,SC_CML,Fraction reads in cell-associated barcodes,67.0%
Library,Multiplexing Capture,Physical library ID,SC_CML,Fraction unrecognized CMO,2.3%
Library,Multiplexing Capture,Physical library ID,SC_CML,Mean reads per cell-associated barcode,"13,915"
Library,Multiplexing Capture,Physical library ID,SC_CML,Number of reads,"204,261,166"
Library,Multiplexing Capture,Physical library ID,SC_CML,Sequencing saturation,29.5%
Library,Multiplexing Capture,Physical library ID,SC_CML,Valid UMIs,100.0%
Library,Multiplexing Capture,Physical library ID,SC_CML,Valid barcodes,99.2%
"""

MULTIOME_LIBRARIES = """#Local sample\tLinked sample
PB2_ATAC\tNEXTSEQ_210111/12:PB_GEX/PB2_GEX
PB1_ATAC\tNEXTSEQ_210111/12:PB_GEX/PB1_GEX
"""

MULTIOME_LIBRARIES_NO_RUN = """#Local sample\tLinked sample
PB2_ATAC\tPB_GEX/PB2_GEX
PB1_ATAC\tPB_GEX/PB1_GEX
"""
