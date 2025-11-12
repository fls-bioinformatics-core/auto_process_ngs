#!/usr/bin/env python3
#
#     bcl2fastq.protocols.py: Fastq generation protocol definitions
#     Copyright (C) University of Manchester 2025 Peter Briggs
#

"""
Defines the ``PROTOCOLS`` dictionary, which specifies the available Fastq
generation protocols which can be used in the Fastq generation pipeline.

Each protocol is defined as a name (the dictionary key) and a dictionary
of parameters (the dictionary value).

The following parameters compulsory for each protocol:

* description: free-text string describing the protocol
* pipeline_variant: a string specifying which of the implemented sub-pipelines
  the protocol will use; must be one of: "standard", "10x_cellranger",
  "10x_cellranger-atc", "10x_cellranger-arc", "10x_spaceranger"
* supported_indexes: a sequence (list or tuple) which specifies the types
  of sample sheet index that the protocol supports; must be one or more of:
  "ILLUMINA", "10X", "NONE"

The remaining parameters can be none or more of any of the possible lane
subset attributes. Most commonly:

* r1_length
* r2_length
* r3_length
* i1_length
* i2_length
* minimum_trimmed_read_lengths
* mask_short_adapter_reads
* no_lane_splitting
* create_fastq_for_index_read
* trim_adapters
"""

PROTOCOLS = {
    "standard": {
        "description": "Standard Illumina sequencing data (default)",
        "pipeline_variant": "standard",
        "supported_indexes": ("ILLUMINA", "NONE"),
    },
    "mirna": {
        # miRNA-seq protocol
        # Set minimum trimmed read length and turn off masking
        "description": "miRNA-seq data",
        "pipeline_variant": "standard",
        "supported_indexes": ("ILLUMINA",),
        "minimum_trimmed_read_length": 10,
        "mask_short_adapter_reads": 0,
    },
    "10x_chromium_sc": {
        # 10xGenomics Chromium SC (GEX/Flex)
        # -- truncate R1 to 28 bases
        # -- truncate R2 to 90 bases
        # -- truncate I1 and I2 to 10 bases
        # -- minimum trimmed read length 8bp
        # -- minimum masked read length 8bp
        # -- split by lane
        # -- create Fastqs for index read
        # -- disable adapter trimming
        "description": "10x Genomics Chromium 3' and 5' single cell "
                        "gene expression data",
        "pipeline_variant": "10x_cellranger",
        "supported_indexes": ("ILLUMINA", "10X"),
        "r1_length": 28,
        "r2_length": 90,
        "i1_length": 10,
        "i2_length": 10,
        "minimum_trimmed_read_length": 8,
        "mask_short_adapter_reads": 8,
        "no_lane_splitting": False,
        "create_fastq_for_index_read": True,
        "trim_adapters": False
    },
    "10x_atac": {
        # 10xGenomics ATAC-seq
        # -- convert I2 to R2
        # -- truncate R1 to 50 bases
        # -- truncate R2 to 16 bases
        # -- truncate R3 to 50 bases
        # -- truncate I1 to 8 bases
        # -- enable filter single index
        # -- split by lane
        # -- create Fastqs for index read
        # -- disable adapter trimming
        "description": "10x Genomics Chromium single cell ATAC-seq data",
        "pipeline_variant": "10x_cellranger-atac",
        "supported_indexes": ("10X",),
        "r1_length": 50,
        "r2_length": 16,
        "r3_length": 50,
        "i1_length": 8,
        "override_template": "RIRR",
        "tenx_filter_single_index": True,
        "no_lane_splitting": False,
        "create_fastq_for_index_read": True,
        "trim_adapters": False
    },
    "10x_multiome": {
        # 10xGenomics multiome
        # -- set bases mask to "auto"
        # -- split by lane
        # -- create Fastqs for index read
        # -- disable adapter trimming
        "description": "10x Genomics single cell multiome data "
        "(unpooled data i.e. ATAC or GEX data only in single run)",
        "pipeline_variant": "10x_cellranger-arc",
        "supported_indexes": ("10X",),
        "bases_mask": "auto",
        "no_lane_splitting": False,
        "create_fastq_for_index_read": True,
        "trim_adapters": False
    },
    "10x_multiome_atac": {
        # 10xGenomics multiome (ATAC)
        # -- convert I2 to R2
        # -- truncate R1 to 50 bases
        # -- truncate R2 to 24 bases
        # -- truncate R3 to 49 bases
        # -- truncate I1 to 8 bases
        # -- enable filter single index
        # -- split by lane
        # -- create Fastqs for index read
        # -- disable adapter trimming
        "description": "10x Genomics single cell multiome ATAC-seq data "
        "(run with pooled GEX and ATAC data)",
        "pipeline_variant": "10x_cellranger-arc",
        "supported_indexes": ("10X",),
        "r1_length": 50,
        "r2_length": 24,
        "r3_length": 49,
        "i1_length": 8,
        "override_template": "RIRR",
        "tenx_filter_single_index": True,
        "no_lane_splitting": False,
        "create_fastq_for_index_read": True,
        "trim_adapters": False
    },
    "10x_multiome_gex": {
        # 10xGenomics multiome (GEX)
        # -- truncate I1 and I2 to 10 bases
        # -- truncate R1 to 28 bases
        # -- truncate R2 to 90 bases
        # -- enable filter dual index
        # -- split by lane
        # -- create Fastqs for index read
        # -- disable adapter trimming
        "description": "10x Genomics single cell multiome GEX data "
        "(run with pooled GEX and ATAC data)",
        "pipeline_variant": "10x_cellranger-arc",
        "supported_indexes": ("10X",),
        "r1_length": 28,
        "r2_length": 90,
        "i1_length": 10,
        "i2_length": 10,
        "tenx_filter_dual_index": True,
        "no_lane_splitting": False,
        "create_fastq_for_index_read": True,
        "trim_adapters": False
    },
    "10x_visium": {
        # 10xGenomics Visium
        # -- truncate R1 to 28 bases
        # -- truncate R2 to 50 bases
        # -- truncate I1 and I2 to 10 bases
        # -- minimum trimmed read length 8bp
        # -- minimum masked read length 8bp
        # -- split by lane
        # -- create Fastqs for index read
        # -- disable adapter trimming
        "description": "10x Genomics Visium CytAssist FFPE, Fresh Frozen, Fixed "
        "Frozen spatial GEX or FFPE PEX data",
        "pipeline_variant": "10x_spaceranger",
        "supported_indexes": ("ILLUMINA", "10X"),
        "r1_length": 28,
        "r2_length": 50,
        "i1_length": 10,
        "i2_length": 10,
        "minimum_trimmed_read_length": 8,
        "mask_short_adapter_reads": 8,
        "no_lane_splitting": False,
        "create_fastq_for_index_read": True,
        "trim_adapters": False
    },
    "10x_visium_v1": {
        # 10xGenomics Visium v1
        # -- truncate R1 to 28 bases
        # -- truncate R2 to 90 bases
        # -- truncate I1 and I2 to 10 bases
        # -- minimum trimmed read length 8bp
        # -- minimum masked read length 8bp
        # -- split by lane
        # -- create Fastqs for index read
        # -- disable adapter trimming
        "description": "10x Genomics Visium Fresh Frozen Spatial GEX (v1) "
        "data (no CytAssist)",
        "pipeline_variant": "10x_spaceranger",
        "supported_indexes": ("ILLUMINA", "10X"),
        "r1_length": 28,
        "r2_length": 90,
        "i1_length": 10,
        "i2_length": 10,
        "minimum_trimmed_read_length": 8,
        "mask_short_adapter_reads": 8,
        "no_lane_splitting": False,
        "create_fastq_for_index_read": True,
        "trim_adapters": False
    },
    "10x_visium_hd": {
        # 10xGenomics Visium (HD)
        # -- truncate R1 to 43 bases
        # -- truncate R2 to 50 bases
        # -- truncate I1 and I2 to 10 bases
        # -- minimum trimmed read length 8bp
        # -- minimum masked read length 8bp
        # -- no lane splitting
        # -- create Fastqs for index read
        # -- disable adapter trimming
        "description": "10x Genomics Visium CytAssist FFPE HD spatial GEX "
        "data",
        "pipeline_variant": "10x_spaceranger",
        "supported_indexes": ("ILLUMINA", "10X"),
        "r1_length": 43,
        "r2_length": 50,
        "i1_length": 10,
        "i2_length": 10,
        "minimum_trimmed_read_length": 8,
        "mask_short_adapter_reads": 8,
        "no_lane_splitting": True,
        "create_fastq_for_index_read": True,
        "trim_adapters": False
    },
    "10x_visium_hd_3prime": {
        # 10xGenomics Visium (HD 3')
        # -- truncate R1 to 43 bases
        # -- truncate R2 to 75 bases
        # -- truncate I1 and I2 to 10 bases
        # -- minimum trimmed read length 8bp
        # -- minimum masked read length 8bp
        # -- split by lane
        # -- create Fastqs for index read
        # -- disable adapter trimming
        "description": "10x Visium CytAssist FFPE HD 3' spatial GEX data",
        "pipeline_variant": "10x_spaceranger",
        "supported_indexes": ("ILLUMINA", "10X"),
        "r1_length": 43,
        "r2_length": 75,
        "i1_length": 10,
        "i2_length": 10,
        "minimum_trimmed_read_length": 8,
        "mask_short_adapter_reads": 8,
        "no_lane_splitting": False,
        "create_fastq_for_index_read": True,
        "trim_adapters": False
    },
    "parse_evercode": {
        # Parse Evercode
        # Disable adapter trimming
        "description": "Parse Evercode single cell data",
        "pipeline_variant": "standard",
        "supported_indexes": ("ILLUMINA",),
        "trim_adapters": False,
    },
    "biorad_ddseq": {
        # Bio-Rad ddSEQ
        # Disable adapter trimming
        "description": "Bio-Rad ddSEQ single cell data",
        "pipeline_variant": "standard",
        "supported_indexes": ("ILLUMINA",),
        "trim_adapters": False,
    }
}