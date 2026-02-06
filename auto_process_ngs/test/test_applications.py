import unittest
from unittest import TestCase
from auto_process_ngs.applications import identify_application
from auto_process_ngs.applications import fetch_application_data
from auto_process_ngs.applications import split_library_type
from auto_process_ngs.applications import match_application
from auto_process_ngs.applications import score_match


class TestIdentifyApplication(TestCase):
    def test_identify_application_undefined_platform_and_library(self):
        """
        identify_application: identify undefined platform and library
        """
        application = identify_application(None, None)
        self.assertEqual(application["fastq_generation"],"standard")
        self.assertEqual(application["qc_protocol"], "minimal")

    def test_identify_application_rnaseq(self):
        """
        identify_application: identify RNA-seq
        """
        application = identify_application(None, "RNA-seq")
        self.assertEqual(application["fastq_generation"], "standard")
        self.assertEqual(application["qc_protocol"], "standard")

    def test_identify_application_wgs(self):
        """
        identify_application: identify WGS/DNA-seq/CRISPR
        """
        for library_type in ["WGS", "DNA-seq", "Amplicon DNA-seq", "CRISPR", "CRISPR-Cas9"]:
            application = identify_application(None, library_type)
            self.assertEqual(application["fastq_generation"], "standard")
            self.assertEqual(application["qc_protocol"], "minimal")

    def test_identify_application_mirnaseq(self):
        """
        identify_application: identify miRNA-seq
        """
        application = identify_application(None, "miRNA-seq")
        self.assertEqual(application["fastq_generation"], "mirna")
        self.assertEqual(application["qc_protocol"], "standard")

    def test_identify_application_10x_chromium_3prime_scrnaseq(self):
        """
        identify_application: identify 10x Chromium 3' scRNA-seq
        """
        platform = "10x Chromium 3'"
        for library_type in ["scRNA-seq",
                             "scRNA-seq+CSP",
                             "scRNA-seq+CRISPR"]:
            application = identify_application(platform, library_type)
            self.assertEqual(application["fastq_generation"], "10x_chromium_sc")
            self.assertEqual(application["qc_protocol"], "10x_scRNAseq")

    def test_identify_application_10x_chromium_3prime_ocm_scrnaseq(self):
        """
        identify_application: identify 10x Chromium 3' OCM scRNA-seq
        """
        platform = "10x Chromium 3' OCM"
        for library_type in ["scRNA-seq",
                             "scRNA-seq+CSP",
                             "scRNA-seq+CRISPR"]:
            application = identify_application(platform, library_type)
            self.assertEqual(application["fastq_generation"], "10x_chromium_sc")
            self.assertEqual(application["qc_protocol"], "10x_scRNAseq")

    def test_identify_application_10x_chromium_3prime_snrnaseq(self):
        """
        identify_application: identify 10x Chromium 3' snRNA-seq
        """
        platform = "10x Chromium 3'"
        for library_type in ["snRNA-seq",
                             "snRNA-seq+CSP",
                             "snRNA-seq+CRISPR"]:
            application = identify_application(platform, library_type)
            self.assertEqual(application["fastq_generation"], "10x_chromium_sc")
            self.assertEqual(application["qc_protocol"], "10x_snRNAseq")

    def test_identify_application_10x_chromium_3prime_ocm_snrnaseq(self):
        """
        identify_application: identify 10x Chromium 3' OCM snRNA-seq
        """
        platform = "10x Chromium 3' OCM"
        for library_type in ["snRNA-seq",
                             "snRNA-seq+CSP",
                             "snRNA-seq+CRISPR"]:
            application = identify_application(platform, library_type)
            self.assertEqual(application["fastq_generation"], "10x_chromium_sc")
            self.assertEqual(application["qc_protocol"], "10x_snRNAseq")

    def test_identify_application_10x_chromium_5prime_immune_profiling(self):
        """
        identify_application: identify 10x Chromium 5' Immune Profiling
        """
        platform = "10x Chromium 5'"
        for library_type in ["Immune Profiling",
                             "Immune Profiling+CSP",
                             "Immune Profiling+VDJ",
                             "Immune Profiling+VDJ+CSP"]:
            application = identify_application(platform, library_type)
            self.assertEqual(application["fastq_generation"], "10x_chromium_sc")
            self.assertEqual(application["qc_protocol"], "10x_ImmuneProfiling")

    def test_identify_application_10x_chromium_5prime_ocm_immune_profiling(self):
        """
        identify_application: identify 10x Chromium 5' OCM Immune Profiling
        """
        platform = "10x Chromium 5' OCM"
        for library_type in ["Immune Profiling",
                             "Immune Profiling+CSP",
                             "Immune Profiling+VDJ",
                             "Immune Profiling+VDJ+CSP"]:
            application = identify_application(platform, library_type)
            self.assertEqual(application["fastq_generation"], "10x_chromium_sc")
            self.assertEqual(application["qc_protocol"], "10x_ImmuneProfiling")

    def test_identify_application_10x_chromium_flex(self):
        """
        identify_application: identify 10x Chromium Flex sc/snRNA-seq
        """
        platform = "10x Chromium Flex"
        for library_type in ["GEX", "GEX+PEX"]:
            application = identify_application(platform, library_type)
            self.assertEqual(application["fastq_generation"], "10x_chromium_sc")
            self.assertEqual(application["qc_protocol"], "10x_Flex")

    def test_identify_application_10x_chromium_3prime_cellplex(self):
        """
        identify_application: identify 10x Chromium 3' CellPlex
        """
        platform = "10x Chromium 3' CellPlex"
        library_type = "scRNA-seq"
        application = identify_application(platform, library_type)
        self.assertEqual(application["fastq_generation"], "10x_chromium_sc")
        self.assertEqual(application["qc_protocol"], "10x_CellPlex")

    def test_identify_application_10x_single_cell_atac(self):
        """
        identify_application: identify 10x Single Cell ATAC
        """
        platform = "10x Single Cell ATAC"
        library_type = "snATAC-seq"
        application = identify_application(platform, library_type)
        self.assertEqual(application["fastq_generation"], "10x_atac")
        self.assertEqual(application["qc_protocol"], "10x_scATAC")

    def test_identify_application_10x_single_cell_multiome_atac(self):
        """
        identify_application: identify 10x Single Cell Multiome ATAC
        """
        platform = "10x Single Cell Multiome"
        library_type = "ATAC"
        application = identify_application(platform, library_type)
        self.assertEqual(application["fastq_generation"], "10x_multiome_atac")
        self.assertEqual(application["qc_protocol"], "10x_Multiome_ATAC")

    def test_identify_application_10x_single_cell_multiome_gex(self):
        """
        identify_application: identify 10x Single Cell Multiome GEX
        """
        platform = "10x Single Cell Multiome"
        library_type = "GEX"
        application = identify_application(platform, library_type)
        self.assertEqual(application["fastq_generation"], "10x_multiome_gex")
        self.assertEqual(application["qc_protocol"], "10x_Multiome_GEX")

    def test_identify_application_10x_visium(self):
        """
        identify_application: identify 10x Visium applications
        """
        application = identify_application("10x Visium", "Fresh Frozen Spatial RNA-seq")
        self.assertEqual(application["fastq_generation"], "10x_visium_v1")
        self.assertEqual(application["qc_protocol"], "10x_Visium_GEX_90bp_insert")

    def test_identify_application_10x_visium_cytassist_gex(self):
        """
        identify_application: identify 10x Visium CytAssist GEX applications
        """
        for library_type in ["FFPE Spatial RNA-seq",
                             "Fresh Frozen Spatial RNA-seq",
                             "Fixed Frozen Spatial RNA-seq"]:
            application = identify_application("10x Visium (CytAssist)", library_type)
            self.assertEqual(application["fastq_generation"], "10x_visium")
            self.assertEqual(application["qc_protocol"], "10x_Visium_GEX")

    def test_identify_application_10x_visium_cytassist_pex(self):
        """
        identify_application: identify 10x Visium PEX applications
        """
        application = identify_application("10x Visium (CytAssist)", "FFPE Spatial PEX")
        self.assertEqual(application["fastq_generation"], "10x_visium")
        self.assertEqual(application["qc_protocol"], "10x_Visium_PEX")

    def test_identify_application_10x_visium_cytassist_hd_gex(self):
        """
        identify_application: identify 10x Visium CytAssist HD GEX applications
        """
        application = identify_application("10x Visium (CytAssist)", "FFPE HD Spatial RNA-seq")
        self.assertEqual(application["fastq_generation"], "10x_visium_hd")
        self.assertEqual(application["qc_protocol"], "10x_Visium_GEX")

    def test_identify_application_10x_visium_cytassist_hd_3prime_gex(self):
        """
        identify_application: identify 10x Visium CytAssist HD 3' GEX applications
        """
        application = identify_application("10x Visium (CytAssist)", "Fresh Frozen HD 3' Spatial RNA-seq")
        self.assertEqual(application["fastq_generation"], "10x_visium_hd_3prime")
        self.assertEqual(application["qc_protocol"], "10x_Visium_GEX")

    def test_identify_application_parse_evercode(self):
        """
        identify_application: identify Parse Evercode applications
        """
        for library_type in ["scRNA-seq", "snRNA-seq", "TCR", "TCR scRNA-seq", "WT", "WT scRNA-seq"]:
            application = identify_application("Parse Evercode", library_type)
            self.assertEqual(application["fastq_generation"], "parse_evercode")
            self.assertEqual(application["qc_protocol"], "ParseEvercode")

    def test_identify_application_biorad_ddseq_scrna_seq(self):
        """
        identify_application: identify BioRad ddSEQ scRNA-seq applications
        """
        for library_type in ["scRNA-seq", "snRNA-seq"]:
            application = identify_application("Bio-Rad ddSEQ Single Cell 3' RNA-Seq", library_type)
            self.assertEqual(application["fastq_generation"], "biorad_ddseq")
            self.assertEqual(application["qc_protocol"], "minimal")

    def test_identify_application_biorad_ddseq_atac(self):
        """
        identify_application: identify BioRad ddSEQ ATAC applications
        """
        for library_type in ["scATAC-seq", "snATAC-seq"]:
            application = identify_application("Bio-Rad ddSEQ Single Cell ATAC", library_type)
            self.assertEqual(application["fastq_generation"], "biorad_ddseq")
            self.assertEqual(application["qc_protocol"], "BioRad_ddSEQ_ATAC")

    def test_identify_application_legacy_platforms(self):
        """
        identify_application: identify legacy applications
        """
        # 10x Chromium 3' single cell
        for platform in ["10xGenomics Chromium 3'", "10xGenomics Chromium 3'v3", "10xGenomics Chromium 3'v3.1",
                         "10xGenomics Chromium GEM-X", "10xGenomics Chromium GEM-X 3'",
                         "10xGenomics Chromium Next GEM"]:
            application = identify_application(platform, "scRNA-seq")
            self.assertEqual(application["fastq_generation"], "10x_chromium_sc")
            self.assertEqual(application["qc_protocol"], "10x_scRNAseq")

        # 10x Chromium 3' single nuclei
        for platform in ["10xGenomics Chromium 3'", "10xGenomics Chromium 3'v3", "10xGenomics Chromium 3'v3.1",
                         "10xGenomics Chromium GEM-X", "10xGenomics Chromium GEM-X 3'",
                         "10xGenomics Chromium Next GEM"]:
            application = identify_application(platform, "snRNA-seq")
            self.assertEqual(application["fastq_generation"], "10x_chromium_sc")
            self.assertEqual(application["qc_protocol"], "10x_snRNAseq")

        # 10x Chromium CellPlex
        for platform in ["10xGenomics Chromium 3'", "10xGenomics Chromium 3'v3", "10xGenomics Chromium 3'v3.1",
                         "10xGenomics Chromium GEM-X", "10xGenomics Chromium GEM-X 3'",
                         "10xGenomics Chromium Next GEM"]:
            for library_type in ["CellPlex", "CellPlex scRNA-seq", "CellPlex snRNA-seq"]:
                application = identify_application(platform, library_type)
                self.assertEqual(application["fastq_generation"], "10x_chromium_sc")
                self.assertEqual(application["qc_protocol"], "10x_CellPlex")

        # 10x Chromium 5'
        application = identify_application("10xGenomics Chromium 5'", "Single Cell Immune Profiling")
        self.assertEqual(application["fastq_generation"], "10x_chromium_sc")
        self.assertEqual(application["qc_protocol"], "10x_ImmuneProfiling")

        # 10x Chromium Flex
        for platform in ["10xGenomics Chromium 3'", "10xGenomics Chromium 3'v3", "10xGenomics Chromium 3'v3.1",
                         "10xGenomics Chromium GEM-X", "10xGenomics Chromium GEM-X 3'",
                         "10xGenomics Chromium Next GEM"]:
            application = identify_application(platform, "Flex")
            self.assertEqual(application["fastq_generation"], "10x_chromium_sc")
            self.assertEqual(application["qc_protocol"], "10x_Flex")

        # 10x ATAC
        for library_type in ["scATAC-seq", "snATAC-seq"]:
            application = identify_application("10xGenomics Single Cell ATAC", library_type)
            self.assertEqual(application["fastq_generation"], "10x_atac")
            self.assertEqual(application["qc_protocol"], "10x_scATAC")

        # 10x Multiome
        application = identify_application("10xGenomics Single Cell Multiome", "ATAC")
        self.assertEqual(application["fastq_generation"], "10x_multiome_atac")
        self.assertEqual(application["qc_protocol"], "10x_Multiome_ATAC")

        application = identify_application("10xGenomics Single Cell Multiome", "GEX")
        self.assertEqual(application["fastq_generation"], "10x_multiome_gex")
        self.assertEqual(application["qc_protocol"], "10x_Multiome_GEX")

        # 10x Visium
        for platform in ["10xGenomics Visium", "10xGenomics Visium (CytAssist)", "10xGenomics CytAssist Visium"]:
            for library_type in ["FFPE Spatial PEX", "FFPE Spatial Protein Expression"]:
                application = identify_application(platform, library_type)
                self.assertEqual(application["fastq_generation"], "10x_visium")
                self.assertEqual(application["qc_protocol"], "10x_Visium_PEX")

        for library_type in ["Fresh Frozen Spatial GEX", "Fresh Frozen Spatial Gene Expression"]:
            application = identify_application("10xGenomics Visium", library_type)
            self.assertEqual(application["fastq_generation"], "10x_visium_v1")
            self.assertEqual(application["qc_protocol"], "10x_Visium_GEX_90bp_insert")

        for platform in ["10xGenomics Visium", "10xGenomics CytAssist Visium"]:
            for library_type in ["Spatial RNA-seq", "spatial RNA-seq"]:
                application = identify_application(platform, library_type)
                self.assertEqual(application["fastq_generation"], "10x_visium")
                self.assertEqual(application["qc_protocol"], "10x_Visium_legacy")

        for platform in ["10xGenomics Visium", "10xGenomics Visium (CytAssist)", "10xGenomics CytAssist Visium"]:
            for library_type in ["Spatial GEX", "FFPE Spatial GEX"]:
                application = identify_application(platform, library_type)
                self.assertEqual(application["fastq_generation"], "10x_visium")
                self.assertEqual(application["qc_protocol"], "10x_Visium_GEX")

    def test_identify_application_unrecognised_platform(self):
        """
        identify_application: handle "unrecognised" platform
        """
        application = identify_application("Unknown", "Unknown")
        self.assertEqual(application["fastq_generation"], "standard")
        self.assertEqual(application["qc_protocol"], "minimal")


class TestFetchApplicationData(TestCase):
    def test_fetch_application_data_no_matching_tag(self):
        """
        fetch_application_data: no matches for nonexistent tag
        """
        no_applications = fetch_application_data(["nonexistent"])
        self.assertEqual(len(no_applications), 0)

    def test_fetch_application_data_single_tag(self):
        """
        fetch_application_data: match single tag
        """
        legacy_applications = fetch_application_data(["legacy"])
        for application in legacy_applications:
            self.assertTrue("legacy" in application["tags"])

    def test_fetch_application_data_exclude_single_tag(self):
        """
        fetch_application_data: exclude single tag
        """
        legacy_applications = fetch_application_data(["!legacy"])
        for application in legacy_applications:
            self.assertFalse("legacy" in application["tags"] if "tags" in application else False)

    def test_fetch_application_data_multiple_tags(self):
        """
        fetch_application_data: match multiple tags
        """
        multiple_tags = fetch_application_data(["10x", "single_cell"])
        for application in multiple_tags:
            self.assertTrue("10x" in application["tags"])
            self.assertTrue("single_cell" in application["tags"])

    def test_fetch_application_data_mix_include_and_exclude_tags(self):
        multiple_tags = fetch_application_data(["!10x", "single_cell"])
        for application in multiple_tags:
            self.assertFalse("10x" in application["tags"])
            self.assertTrue("single_cell" in application["tags"])


class TestSplitLibraryType(TestCase):
    def test_split_library_type(self):
        """
        split_library_type: no extensions
        """
        self.assertEqual(split_library_type("RNA-seq"), ("RNA-seq", []))

    def test_split_library_type_single_extension(self):
        """
        split_library_type: single extension
        """
        self.assertEqual(split_library_type("GEX+CSP"), ("GEX", ["CSP"]))

    def test_split_library_type_multiple_extensions(self):
        """
        split_library_type: multiple extensions
        """
        self.assertEqual(split_library_type("GEX+CSP+CRISPR"), ("GEX", ["CSP", "CRISPR"]))

    def test_split_library_type_none_type(self):
        """
        split_library_type: library is 'None'
        """
        self.assertEqual(split_library_type(None), (None, None))

    def test_split_library_type_empty(self):
        """
        split_library_type: library is empty string'None'
        """
        self.assertEqual(split_library_type(""), ("", []))

class TestMatchApplication(TestCase):
    def test_match_application_exact_library_no_platform_no_library(self):
        """
        match_application: match exact library (no platform or library)
        """
        platform_info = {
            "platforms": [],
            "libraries": []
        }
        # Matches if no platform or library specified
        self.assertEqual(match_application(platform_info, None, None),
                         (platform_info, [], []))
        # Doesn't match if platform specified
        self.assertEqual(match_application(platform_info, "10x Chromium 3'", None),
                         None)
        # Doesn't match if library specified
        self.assertEqual(match_application(platform_info, None, "RNA-seq"),
                         None)

    def test_match_application_exact_library_no_platform(self):
        """
        match_application: match exact library (no platform)
        """
        platform_info = {
            "platforms": [],
            "libraries": ["RNA-seq"]
        }
        # Matches if no platform specified
        self.assertEqual(match_application(platform_info, None, "RNA-seq"),
                         (platform_info, [], ["RNA-seq"]))
        # Doesn't match for different library
        self.assertEqual(match_application(platform_info, None, "ATAC-seq"),
                         None)
        # Doesn't match if platform specified
        self.assertEqual(match_application(platform_info, "10x Chromium 3'", "RNA-seq"),
                         None)

    def test_match_application_exact_library_exact_platform(self):
        """
        match_application: match exact library (exact platform)
        """
        platform_info = {
            "platforms": ["10x Chromium 3'",],
            "libraries": ["scRNA-seq"]
        }
        # Matches if exact library and platform specified
        self.assertEqual(match_application(platform_info, "10x Chromium 3'", "scRNA-seq"),
                         (platform_info, ["10x Chromium 3'"], ["scRNA-seq"]))
        # Doesn't match for different library
        self.assertEqual(match_application(platform_info, "10x Chromium 3'", "snRNA-seq"),
                         None)
        # Doesn't match if no platform specified
        self.assertEqual(match_application(platform_info, None, "scRNA-seq"),
                         None)

    def test_match_application_exact_library_exact_platform_multiple_options(self):
        """
        match_application: match exact library (exact platform) (multiple options)
        """
        platform_info = {
            "platforms": ["10x Chromium 3' (GEM-X)",
                          "10x Chromium 3' (Next GEM)"],
            "libraries": ["scRNA-seq", "GEX"]
        }
        # Matches if exact libraries and platforms specified
        self.assertEqual(match_application(platform_info, "10x Chromium 3' (GEM-X)", "scRNA-seq"),
                         (platform_info, ["10x Chromium 3' (GEM-X)"], ["scRNA-seq"]))
        self.assertEqual(match_application(platform_info, "10x Chromium 3' (GEM-X)", "GEX"),
                         (platform_info, ["10x Chromium 3' (GEM-X)"], ["GEX"]))
        self.assertEqual(match_application(platform_info, "10x Chromium 3' (Next GEM)", "scRNA-seq"),
                         (platform_info, ["10x Chromium 3' (Next GEM)"], ["scRNA-seq"]))
        self.assertEqual(match_application(platform_info, "10x Chromium 3' (Next GEM)", "GEX"),
                         (platform_info, ["10x Chromium 3' (Next GEM)"], ["GEX"]))
        # Doesn't match for different library
        self.assertEqual(match_application(platform_info, "10x Chromium 3'", "snRNA-seq"),
                         None)
        # Doesn't match if no platform specified
        self.assertEqual(match_application(platform_info, None, "scRNA-seq"),
                         None)

    def test_match_application_any_library_no_platform(self):
        """
        match_application: match any library (no platform)
        """
        platform_info = {
            "platforms": [],
            "libraries": ["*"]
        }
        # Matches if no platform specified
        self.assertEqual(match_application(platform_info, None, "RNA-seq"),
                         (platform_info, [], ["*"]))
        # Doesn't match if platform specified
        self.assertEqual(match_application(platform_info, "10x Chromium 3'", "scRNA-seq"),
                         None)

    def test_match_application_any_library_any_platform(self):
        """
        match_application: match any library (any platform)
        """
        platform_info = {
            "platforms": ["*"],
            "libraries": ["*"]
        }
        # Matches for no library and no platform specified
        self.assertEqual(match_application(platform_info, None, None),
                         (platform_info, ["*"], ["*"]))
        # Matches for library but no platform specified
        self.assertEqual(match_application(platform_info, None, "RNA-seq"),
                         (platform_info, ["*"], ["*"]))

    def test_match_application_ignore_library_extensions(self):
        """
        match_application: ignore library extensions
        """
        # No platform
        application_info = {
            "platforms": ["10x Chromium 3'",],
            "libraries": ["GEX"]
        }
        # Matches if no extensions specified
        self.assertEqual(match_application(application_info, "10x Chromium 3'", "GEX"),
                         (application_info, ["10x Chromium 3'"], ["GEX"]))
        # Matches with extensions
        self.assertEqual(match_application(application_info, "10x Chromium 3'", "GEX+CSP"),
                         (application_info, ["10x Chromium 3'"], ["GEX"]))

class TestScoreMatch(unittest.TestCase):
    def test_score_match(self):
        """
        score_match: score single matches
        """
        self.assertEqual(score_match(["10x Chromium 3'"], ["scRNA-seq"]), 0)
        self.assertEqual(score_match([], ["RNA-seq"]), 0)
        self.assertEqual(score_match(["10x Chromium 3'"], ["s*RNA-seq"]), 1)
        self.assertEqual(score_match(["10x Chromium 3'*"], ["s*RNA-seq"]), 2)
        self.assertEqual(score_match([], ["*"]), 1)
        self.assertEqual(score_match(["*"], ["*"]), 2)

    def test_score_match_multiple_matches(self):
        """
        score_match: score multiple matches
        """
        self.assertEqual(score_match(["10x Chromium 3' (GEM-X)", "10x Chromium 3'*"],
                                     ["scRNA-seq"]), 0)
        self.assertEqual(score_match(["10x Chromium 3' (GEM-X)", "10x Chromium 3'*"],
                                     ["scRNA-seq", "s*RNA-seq"]), 0)
        self.assertEqual(score_match(["10x Chromium 3'*"],
                                     ["scRNA-seq", "s*RNA-seq"]), 1)
        self.assertEqual(score_match(["10x Chromium 3' (GEM-X)", "10x Chromium 3'*"],
                                     ["s*RNA-seq"]), 1)