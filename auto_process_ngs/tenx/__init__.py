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
    "10xGenomics Chromium GEM-X 3'v4",
    "10xGenomics Chromium Next GEM",
    "10xGenomics Chromium Next GEM 3'v3.1",
    "10xGenomics Single Cell ATAC",
    "10xGenomics Visium",
    "10xGenomics Visium (CytAssist)",
    "10xGenomics CytAssist Visium",
    "10xGenomics Single Cell Multiome",
)

# List of known library types for 10xGenomics
LIBRARIES = {
    "10xGenomics Chromium GEM-X 3'*": (
        "scRNA-seq",
    ),
    "10xGenomics Chromium Next GEM*": (
        "scRNA-seq",
        "CellPlex scRNA-seq",
        "Flex",
    ),
    "10xGenomics Chromium*": (
        "scRNA-seq",
        "snRNA-seq",
        "CellPlex",
        "CellPlex scRNA-seq",
        "Flex",
        "Single Cell Immune Profiling",
    ),
    "10xGenomics Single Cell ATAC": (
        "scATAC-seq",
        "snATAC-seq",
    ),
    "10xGenomics Single Cell Multiome": (
        "ATAC",
        "GEX",
    ),
    "10xGenomics Visium": (
        "Fresh Frozen Spatial GEX",
        "FFPE Spatial GEX",
        "Fresh Frozen Spatial Gene Expression",
        "FFPE Spatial Gene Expression",
    ),
    "10xGenomics Visium (CytAssist)": (
        "FFPE HD Spatial GEX",
        "FFPE Spatial GEX",
        "Fixed Frozen Spatial GEX",
        "Fresh Frozen Spatial GEX",
        "FFPE Spatial PEX",
        "FFPE HD Spatial Gene Expression",
        "FFPE Spatial Gene Expression",
        "Fixed Frozen Spatial Gene Expression",
        "Fresh Frozen Spatial Gene Expression",
        "FFPE Spatial Protein Expression",
    ),
    # Legacy Visium platform/application combinations kept for
    # backwards compatibility
    "10xGenomics * Visium": (
        "FFPE Spatial RNA-seq",
        "Fresh Frozen RNA-seq",
        "GEX",
        "PEX",
        "Spatial RNA-seq",
        "HD Spatial GEX",
    ),
}

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
DEFAULT_CELLRANGER_VERSION = "8.0.0"
