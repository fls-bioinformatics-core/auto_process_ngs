#######################################################################
# Data
#######################################################################

# Permissible values for 10xGenomics platforms
PLATFORMS = (
    "10xGenomics Chromium",
    "10xGenomics Chromium 3'",
    "10xGenomics Chromium 3'v2",
    "10xGenomics Chromium 3'v3",
    "10xGenomics Chromium 3'v3.1",
    "10xGenomics Chromium 5'",
    "10xGenomics Single Cell ATAC",
    "10xGenomics Visium",
    "10xGenomics CytAssist Visium",
    "10xGenomics Single Cell Multiome",
)

# Permissible values for cellranger count --chemistry option
# See https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count
CELLRANGER_ASSAY_CONFIGS = {
    'auto': 'autodetection',
    'threeprime': 'Single Cell 3\'',
    'fiveprime': 'Single Cell 5\'',
    'SC3Pv1': 'Single Cell 3\' v1',
    'SC3Pv2': 'Single Cell 3\' v2',
    'SC3Pv3': 'Single Cell 3\' v3',
    'SC5P-PE': 'Single Cell 5\' paired-end (both R1 and R2 are used for alignment)',
    'SC5P-R2': 'Single Cell 5\' R2-only (where only R2 is used for alignment)',
    'ARC-v1': 'Single Cell Multiome (ATAC+GEX) v1', # Not documented?
}

