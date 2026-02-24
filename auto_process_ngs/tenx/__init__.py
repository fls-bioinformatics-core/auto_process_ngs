#######################################################################
# Data
#######################################################################

# Permissible values for cellranger count --chemistry option
# See https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count
CELLRANGER_ASSAY_CONFIGS = {
    'auto': 'autodetection',
    'threeprime': 'Single Cell 3\'',
    'fiveprime': 'Single Cell 5\'',
    'SC3Pv1': 'Single Cell 3\' v1',
    'SC3Pv2': 'Single Cell 3\' v2',
    'SC3Pv3': 'Single Cell 3\' v3',
    'SC3Pv4': 'Single Cell 3\' v4',
    'SC3Pv3LT': 'Single Cell 3\' v3 LT',
    'SC3Pv3HT': 'Single Cell 3\' v3 HT',
    'SC5P-PE': 'Single Cell 5\' paired-end (both R1 and R2 are used for alignment)',
    'SC5P-PE-v3': 'Single Cell 5\' paired-end (both R1 and R2 are used for alignment) v3',
    'SC5P-R2': 'Single Cell 5\' R2-only (where only R2 is used for alignment)',
    'SC5P-R2-v3': 'Single Cell 5\' R2-only (where only R2 is used for alignment)',
    'ARC-v1': 'Single Cell Multiome (ATAC+GEX) v1 (cannot be autodetected)'
}

# Default Cellranger version
DEFAULT_CELLRANGER_VERSION = "10.0.0"

# Feature types for different library type extensions for Cellranger multi
CELLRANGER_MULTI_EXTENSIONS_TO_FEATURE_TYPES = {
    # Canonical (used in 10x product descriptions)
    "CSP": "Antibody Capture",
    "VDJ-T": "VDJ-T",
    "VDJ-B": "VDJ-B",
    # Non-canonical (for local usage)
    "Antibody Capture": "Antibody Capture",
    "Feature Barcode": "Antibody Capture",
    "TCR": "VDJ-T",
    "BCR": "VDJ-B",
}

# Sample name endings to feature types for Cellranger multi config templates
# Expect names to be of the form e.g. "<SAMPLE>_CSP"
# (These are local conventions not canonical names used by 10x)
CELLRANGER_MULTI_FEATURE_TYPES = {
    "GEX": "Gene Expression",
    "GE": "Gene Expression",
    "FLEX": "Gene Expression",
    "CML": "Multiplexing Capture",
    "CSP": "Antibody Capture",
    "TCR": "VDJ-T",
    "BCR": "VDJ-B",
}