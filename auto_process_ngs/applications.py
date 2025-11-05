#!/usr/bin/env python3
#
#     applications.py: information about platforms and libraries
#     Copyright (C) University of Manchester 2025 Peter Briggs
#

"""
Collects information about sequencing applications (combinations of platform
and library type).

Defines the ``APPLICATIONS`` list, which contains dictionaries defining known
applications, and the ``identify_application()`` function, which can be used
to identify the appropriate application for a given platform and library type.

Each application is defined as a dictionary with the following keys:

* platforms: list of platform names (with optional wildcards) which
  correspond to the application; use ["*"] to indicate all platforms, or
  an empty list to indicate that the platform must not be set
* libraries: list of library types (with optional wildcards) which
  correspond to the application; use ["*"] to indicate all library types, or
  an empty list to indicate that the library type must not be set
* fastq_generation: name of the Fastq generation protocol to use for
  this application
* qc_protocol: name of the QC protocol to use for this application
* setup: dictionary with information about actions that should be performed
  as part of setting up analysis project directories for this
  application; contains the following keys:
  - templates: list of template names to use for this application; can be
    one or more of: "10x_multi_config", "10x_multiome_libraries".
  - directories: list of subdirectory names which will be created in the
    analysis project directory (for example "Visium_images")
* tags: optional list of tags to associate with this application; tags can
  be one or more of: "10x", "bio_rad", "parse", "single_cell", "spatial",
  "legacy". Tags are used for automated documentation generation.

The minimum required keys for each application are ``platforms``, ``libraries``,
``fastq_generation``, and ``qc_protocol``.

The module also defines the following user-facing function:

* ``identify_application``: returns the dictionary defining the application
    which matches a given platform and library type

The following functions are also defined for internal use:

* ``match_application``: determines whether a given platform and library type
  match a given application definition
* ``score_match``: returns a score for a platform/library combination
"""
from fnmatch import fnmatch

# Module specific logger
import logging
logger = logging.getLogger(__name__)

APPLICATIONS = [
    # Minimal (no defined platform or defined library)
    {
        "platforms": [],
        "libraries": [],
        "fastq_generation": "standard",
        "qc_protocol": "minimal",
        "setup": {
            "templates": [],
            "directories": [],
        }
    },
    # Minimal (any defined platform & defined library)
    {
        "platforms": ["*"],
        "libraries": ["*"],
        "fastq_generation": "standard",
        "qc_protocol": "minimal",
        "setup": {
            "templates": [],
            "directories": [],
        }
    },
    # Standard (no platform, any defined library)
    {
        "platforms": [],
        "libraries": ["*"],
        "fastq_generation": "standard",
        "qc_protocol": "standard",
        "setup": {
            "templates": [],
            "directories": [],
        }
    },
    # miRNA
    {
        "platforms": [],
        "libraries": ["miRNA-seq"],
        "fastq_generation": "mirna",
        "qc_protocol": "standard",
        "setup": {
            "templates": [],
            "directories": [],
        }
    },
    # WGS/DNA-seq/CRISPR
    {
        "platforms": [],
        "libraries": ["WGS", "DNA-seq", "Amplicon DNA-seq", "CRISPR", "CRISPR-Cas9"],
        "fastq_generation": "standard",
        "qc_protocol": "minimal",
        "setup": {
            "templates": [],
            "directories": [],
        }
},
    # 10x Chromium 3' single-cell
    {
        "platforms": ["10x Chromium 3' (v4 GEM-X)",
                      "10x Chromium 3' (v4 GEM-X) OCM",],
        "libraries": ["GEX", "scRNA-seq"],
        "extensions": ["CSP", "CRISPR"],
        "fastq_generation": "10x_chromium_sc",
        "qc_protocol": "10x_scRNAseq",
        "setup": {
            "templates": ["10x_multi_config"],
            "directories": [],
        },
        "tags": ["10x", "single_cell"]
    },
    # 10x Chromium 3' single-nuclei
    {
        "platforms": ["10x Chromium 3' (v4 GEM-X)",
                      "10x Chromium 3' (v4 GEM-X) OCM",],
        "libraries": ["snRNA-seq"],
        "fastq_generation": "10x_chromium_sc",
        "qc_protocol": "10x_snRNAseq",
        "setup": {
            "templates": ["10x_multi_config"],
            "directories": [],
        },
        "tags": ["10x", "single_cell"]
    },
    # 10x Chromium 5'
    {
        "platforms": ["10x Chromium 5' (v3 GEM-X)",
                      "10x Chromium 5' (v3 GEM-X) OCM",],
        "libraries": ["GEX", "scRNA-seq", "snRNA-seq"],
        "extensions": ["VDJ", "VDJ-T", "VDJ-B", "CSP", "CRISPR", "BEAM"],
        "fastq_generation": "10x_chromium_sc",
        "qc_protocol": "10x_ImmuneProfiling",
        "setup": {
            "templates": ["10x_multi_config"],
            "directories": [],
        },
        "tags": ["10x", "single_cell"]
    },
    # 10x Chromium Flex
    {
        "platforms": ["10x Chromium Flex (v1 GEM-X)",],
        "libraries": ["GEX", "scRNA-seq", "snRNA-seq"],
        "fastq_generation": "10x_chromium_sc",
        "qc_protocol": "10x_Flex",
        "setup": {
            "templates": ["10x_multi_config(include_probeset=true)"],
            "directories": [],
        },
        "tags": ["10x", "single_cell"]
    },
    # 10x Chromium CellPlex
    {
        "platforms": ["10x Chromium 3' (v3.1 Next GEM ST) CellPlex",
                      "10x Chromium 3' (v3.1 Next GEM HT) CellPlex",],
        "libraries": ["GEX", "scRNA-seq"],
        "fastq_generation": "10x_chromium_sc",
        "qc_protocol": "10x_CellPlex",
        "setup": {
            "templates": ["10x_multi_config"],
            "directories": [],
        },
        "tags": ["10x", "single_cell"]
    },
    # 10x Chromium ATAC
    {
        "platforms": ["10x Chromium Epi ATAC (v2)",],
        "libraries": ["ATAC"],
        "fastq_generation": "10x_atac",
        "qc_protocol": "10x_scATAC", # Should be called 10x_snATAC?
        "setup": {
            "templates": [],
            "directories": [],
        },
        "tags": ["10x", "single_cell"]
    },
    # 10x Chromium Multiome
    {
        "platforms": ["10x Chromium Epi Multiome ATAC (v1)"],
        "libraries": ["ATAC"],
        "fastq_generation": "10x_multiome_atac",
        "qc_protocol": "10x_Multiome_ATAC",
        "setup": {
            "templates": ["10x_multiome_libraries"],
            "directories": [],
        },
        "tags": ["10x", "single_cell"]
    },
    {
        "platforms": ["10x Chromium Epi Multiome ATAC (v1)"],
        "libraries": ["GEX"],
        "fastq_generation": "10x_multiome_gex",
        "qc_protocol": "10x_Multiome_GEX",
        "setup": {
            "templates": ["10x_multiome_libraries"],
            "directories": [],
        },
        "tags": ["10x", "single_cell"]
    },
    # 10x Visium
    {
        "platforms": ["10x Visium"],
        "libraries": ["Fresh Frozen Spatial GEX (v1)"],
        "fastq_generation": "10x_visium_v1",
        "qc_protocol": "10x_Visium_GEX_90bp_insert",
        "setup": {
            "templates": [],
            "directories": ["Visium_images"],
        },
        "tags": ["10x", "spatial"]
    },
    {
        "platforms": ["10x Visium (CytAssist)"],
        "libraries": ["FFPE Spatial GEX",
                      "FFPE Spatial GEX (v2)",
                      "Fresh Frozen Spatial GEX (v2)",
                      "Fixed Frozen Spatial GEX (v2)",],
        "fastq_generation": "10x_visium",
        "qc_protocol": "10x_Visium_GEX",
        "setup": {
            "templates": [],
            "directories": ["Visium_images"],
        },
        "tags": ["10x", "spatial"]
    },
    {
        "platforms": ["10x Visium (CytAssist)"],
        "libraries": ["FFPE Spatial PEX",],
        "fastq_generation": "10x_visium",
        "qc_protocol": "10x_Visium_PEX",
        "setup": {
            "templates": [],
            "directories": ["Visium_images"],
        },
        "tags": ["10x", "spatial"]
    },
    {
        "platforms": ["10x Visium (CytAssist)"],
        "libraries": ["FFPE HD Spatial GEX",],
        "fastq_generation": "10x_visium_hd",
        "qc_protocol": "10x_Visium_GEX",
        "setup": {
            "templates": [],
            "directories": ["Visium_images"],
        },
        "tags": ["10x", "spatial"]
    },
    {
        "platforms": ["10x Visium (CytAssist)"],
        "libraries": ["FFPE HD 3' Spatial GEX",],
        "fastq_generation": "10x_visium_hd_3prime",
        "qc_protocol": "10x_Visium_GEX_75bp_insert",
        "setup": {
            "templates": [],
            "directories": ["Visium_images"],
        },
        "tags": ["10x", "spatial"]
    },
    # Parse
    {
        "platforms": ["Parse Evercode"],
        "libraries": ["scRNA-seq",
                      "snRNA-seq",
                      "TCR",
                      "TCR scRNA-seq",
                      "WT",
                      "WT scRNA-seq"],
        "fastq_generation": "parse_evercode",
        "qc_protocol": "ParseEvercode",
        "setup": {
            "templates": [],
            "directories": [],
        },
        "tags": ["parse", "single_cell"]
    },
    # Bio-Rad
    {
        "platforms": ["Bio-Rad ddSEQ Single Cell 3' RNA-Seq"],
        "libraries": ["scRNA-seq",
                      "snRNA-seq"],
        "fastq_generation": "biorad_ddseq",
        "qc_protocol": "minimal", # No protocol defined?
        "setup": {
            "templates": [],
            "directories": [],
        },
        "tags": ["bio_rad", "single_cell"]
    },
    {
        "platforms": ["Bio-Rad ddSEQ Single Cell ATAC"],
        "libraries": ["scATAC-seq",
                      "snATAC-seq"],
        "fastq_generation": "biorad_ddseq",
        "qc_protocol": "BioRad_ddSEQ_ATAC",
        "setup": {
            "templates": [],
            "directories": [],
        },
        "tags": ["bio_rad", "single_cell"]
    },

    # Legacy application definitions below
    # These are retained for backward compatibility

    # 10x Chromium 3' single-cell (legacy)
    {
        "platforms": ["10xGenomics Chromium 3'*",
                      "10xGenomics Chromium GEM-X*",
                      "10xGenomics Chromium Next GEM*"],
        "libraries": ["scRNA-seq"],
        "fastq_generation": "10x_chromium_sc",
        "qc_protocol": "10x_scRNAseq",
        "setup": {
            "templates": [],
            "directories": [],
        },
        "tags": ["10x", "single_cell", "legacy"]
    },
    # 10x Chromium 3' single-nuclei (legacy)
    {
        "platforms": ["10xGenomics Chromium 3'*",
                      "10xGenomics Chromium GEM-X*",
                      "10xGenomics Chromium Next GEM*"],
        "libraries": ["snRNA-seq"],
        "fastq_generation": "10x_chromium_sc",
        "qc_protocol": "10x_snRNAseq",
        "setup": {
            "templates": [],
            "directories": [],
        }
    },
    # 10x Chromium CellPlex (legacy)
    {
        "platforms": ["10xGenomics Chromium 3'*",
                      "10xGenomics Chromium GEM-X*",
                      "10xGenomics Chromium Next GEM*"],
        "libraries": ["CellPlex",
                      "CellPlex scRNA-seq",
                      "CellPlex snRNA-seq"],
        "fastq_generation": "10x_chromium_sc",
        "qc_protocol": "10x_CellPlex",
        "setup": {
            "templates": ["10x_multi_config"],
            "directories": [],
        },
        "tags": ["10x", "single_cell", "legacy"]
    },
    # 10x Chromium 5' (legacy)
    {
        "platforms": ["10xGenomics Chromium 5'*"],
        "libraries": ["Single Cell Immune Profiling"],
        "fastq_generation": "10x_chromium_sc",
        "qc_protocol": "10x_ImmuneProfiling",
        "setup": {
            "templates": ["10x_multi_config"],
            "directories": [],
        }
    },
    # 10x Chromium Flex (legacy)
    {
        "platforms": ["10xGenomics Chromium 3'*",
                      "10xGenomics Chromium GEM-X*",
                      "10xGenomics Chromium Next GEM*"],
        "libraries": ["Flex",],
        "fastq_generation": "10x_chromium_sc",
        "qc_protocol": "10x_Flex",
        "setup": {
            "templates": ["10x_multi_config"],
            "directories": [],
        },
        "tags": ["10x", "single_cell", "legacy"]
    },
    # 10x ATAC (legacy)
    {
        "platforms": ["10xGenomics Single Cell ATAC",],
        "libraries": ["scATAC-seq",
                      "snATAC-seq"],
        "fastq_generation": "10x_atac",
        "qc_protocol": "10x_scATAC", # Should be called 10x_snATAC?
        "setup": {
            "templates": [],
            "directories": [],
        },
        "tags": ["10x", "single_cell", "legacy"]
    },
    # 10x Chromium Multiome (legacy)
    {
        "platforms": ["10xGenomics Single Cell Multiome"],
        "libraries": ["ATAC"],
        "fastq_generation": "10x_multiome_atac",
        "qc_protocol": "10x_Multiome_ATAC",
        "setup": {
            "templates": ["10x_multiome_libraries"],
            "directories": [],
        },
        "tags": ["10x", "single_cell", "legacy"]
    },
    {
        "platforms": ["10xGenomics Single Cell Multiome"],
        "libraries": ["GEX"],
        "fastq_generation": "10x_multiome_gex",
        "qc_protocol": "10x_Multiome_GEX",
        "setup": {
            "templates": ["10x_multiome_libraries"],
            "directories": [],
        },
        "tags": ["10x", "single_cell", "legacy"]
    },
    # 10x Visium (legacy)
    {
        "platforms": ["10xGenomics Visium",
                      "10xGenomics Visium (CytAssist)",
                      "10xGenomics CytAssist Visium"],
        "libraries": ["FFPE Spatial PEX",
                      "FFPE Spatial Protein Expression"],
        "fastq_generation": "10x_visium",
        "qc_protocol": "10x_Visium_PEX",
        "setup": {
            "templates": [],
            "directories": ["Visium_images"],
        },
        "tags": ["10x", "spatial", "legacy"]
    },
    {
        "platforms": ["10xGenomics Visium"],
        "libraries": ["Fresh Frozen Spatial GEX",
                      "Fresh Frozen Spatial Gene Expression"],
        "fastq_generation": "10x_visium_v1",
        "qc_protocol": "10x_Visium_GEX_90bp_insert",
        "setup": {
            "templates": [],
            "directories": ["Visium_images"],
        },
        "tags": ["10x", "spatial", "legacy"]
    },
    {
        "platforms": ["10xGenomics Visium",
                      "10xGenomics CytAssist Visium"],
        "libraries": ["Spatial RNA-seq",
                      "spatial RNA-seq"],
        "fastq_generation": "10x_visium",
        "qc_protocol": "10x_Visium_legacy",
        "setup": {
            "templates": [],
            "directories": ["Visium_images"],
        },
        "tags": ["10x", "spatial", "legacy"]
    },
    # Catch-all for legacy 10x Visium not covered by more
    # specific applications above
    {
        "platforms": ["10xGenomics Visium",
                      "10xGenomics Visium (CytAssist)",
                      "10xGenomics CytAssist Visium"],
        "libraries": ["*"],
        "fastq_generation": "10x_visium",
        "qc_protocol": "10x_Visium_GEX",
        "setup": {
            "templates": [],
            "directories": ["Visium_images"],
        },
        "tags": ["10x", "spatial", "legacy"]
    },
]

def identify_application(platform_name, library_type):
    """
    Returns information about an application

    Applications are combinations of platforms and libraries.

    Arguments:
        platform_name (str): name of the platform
        library_type (str): name of the library

    Returns:
        dict: application-specific information
    """
    # Collect all possible matches
    matching_applications = []
    for platform in APPLICATIONS:
        match = match_application(platform, platform_name, library_type)
        if match:
            matching_applications.append(match)
    if not matching_applications:
        # Failed to find any matches
        raise Exception(f"No matching platforms found for '{platform_name}'/'{library_type}'")
    # Sort matches by score
    matching_applications = sorted(matching_applications, key=lambda m: score_match(m[1], m[2]))
    if len(matching_applications) == 1:
        # If there is a single match then return it
        return matching_applications[0][0]
    # Extract unique scores from matches
    scores = sorted(set([score_match(m[1], m[2]) for m in matching_applications]))
    # Loop over scores in ascending order looking for unique matches
    for score in scores:
        filter_matches = [m for m in matching_applications if score_match(m[1], m[2]) == score]
        if len(filter_matches) == 1:
            # Only one match at this score so return it
            return filter_matches[0][0]
    # No unique match found for any score
    raise Exception(f"No unique matching platform found for '{platform_name}'/'{library_type}'")

def split_library_type(library_type):
    """
    Splits a library type into its components

    Library types are expected to consist of a "base"
    library type followed by none or more optional
    "extensions", which are identified by a preceeding
    '+' character.

    For example: "GEX" has a base library type
    with no extensions; "GEX+CSP" has the base type
    "GEX" with extensions ["CSP"]; and "GEX+CSP+VDJ"
    has the base type "GEX" with extensions
    ["CSP", "VDJ"].

    This function returns a tuple of the form:

    (BASE, EXTENSIONS)

    Arguments:
        library_type (str): name of the library

    Returns:
        tuple: tuple with two elements, first is
        the base library type, the second is
        a list of the extensions.
    """
    if library_type is None:
        return (None, None)
    library_components = library_type.split("+")
    library = library_components[0]
    try:
        return (library, [c.strip() for c in library_components[1:]])
    except IndexError:
        return (library, None)

def match_application(application_info, platform_name, library_type):
    """
    Determine if platform and library type matches the supplied application

    Given information about an application (supplied as a dictionary with
    elements ``platforms`` and ``libraries``), determines whether the
    supplied platform and library match that information.

    FIXME doesn't currently include matching against the library extensions

    Arguments:
         application_info (dict): information about the application
         platform_name (str): name of the platform
         library_type (str): name of the library

    Returns:
        tuple: if the platform and library match the application then
        returns a tuple of the form (application, list of platform matches,
        list of library matches); otherwise returns None
    """
    # See if platform name is a match
    platform_matches = []
    if platform_name is None:
        # No platform specified so only applications which specify all (i.e. '*')
        # or no platforms can be a match
        if not application_info["platforms"] or application_info["platforms"] == ["*"]:
            logger.debug(f"==> matched empty input platform {platform_name}")
            try:
                platform_matches.append(application_info["platforms"][0])
            except IndexError:
                platform_matches.append(None)
        else:
            # Application specifies a platform so it can't be a match
            logger.debug("=> Platforms defined for application, but platform must be empty")
            return None
    else:
        # Platform was specified so only applications which specify at least
        # one platform can be a match
        logger.debug(f"=> looking for match to platform '{platform_name}'")
        if not application_info["platforms"]:
            # Application doesn't specify any platforms so cannot be a match
            logger.debug("=> No platforms defined for application, but looking for a platform")
            return None
        # Iterate over platform names looking for matches
        for platform in application_info["platforms"]:
            logger.debug(f"==> seeing if '{platform}' matches input '{platform_name}'...")
            if fnmatch(platform_name, platform):
                logger.debug("--- yes")
                platform_matches.append(platform)
    # Check if any platforms are a match
    if not platform_matches:
        # No matches in the application for the platform, so reject
        logger.debug(f"=> No matching platforms found for input '{platform_name}'")
        return None
    # See if base library type is a match
    library_matches = []
    library_type, extensions = split_library_type(library_type)
    if library_type is None:
        # No library specified so only applications which specify all (i.e. "*") or
        # no library types can be a match
        if not application_info["libraries"] or application_info["libraries"] == ["*"]:
            logger.debug(f"==> matched empty input library type {library_type}")
            try:
                library_matches.append(application_info["libraries"][0])
            except IndexError:
                library_matches.append(None)
        else:
            # Application defines libraries so it can't be a match
            logger.debug("=> Libraries defined for application, but library type must be empty")
            return None
    else:
        # Library was specified so only applications which specify at least one library
        # type can be a match
        if not application_info["libraries"]:
            return None
        # Iterate over libraries looking for matches
        for library in application_info["libraries"]:
            logger.debug(f"==> seeing if '{library}' matches input '{library_type}'...")
            if fnmatch(library_type, library):
                logger.debug("--- yes")
                library_matches.append(library)
    # Check if any libraries are a match
    if not library_matches:
        # No matches in the application for the library type, so reject
        return None
    # FIXME Ignore matching extensions for now
    # Return matches
    return (application_info,
            [m for m in platform_matches if m is not None],
            [m for m in library_matches if m is not None])


def score_match(platforms, libraries):
    """
    Return a score for a platform/library combination

    The score is calculated as the sum of the minimum number of wildcard
    characters ('*') in the platform and library lists. A lower score
    indicates a more specific match.

    Arguments:
        platforms (list): list of platforms
        libraries (list): list of libraries
    """
    score = 0
    for lst in (platforms, libraries):
        if lst:
            score += min([str(x).count("*") for x in lst if x is not None])
    return score