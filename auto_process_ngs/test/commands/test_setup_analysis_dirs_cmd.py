#######################################################################
# Tests for setup_analysis_dirs_cmd.py module
#######################################################################

import unittest
import os
import tempfile
import shutil
from auto_process_ngs.auto_processor import AutoProcess
from auto_process_ngs.mock import MockAnalysisDirFactory
from auto_process_ngs.commands.setup_analysis_dirs_cmd import setup_analysis_dirs

# Unit tests

class TestSetupAnalysisDirs(unittest.TestCase):
    """
    Tests for the 'setup_analysis_dirs' command
    """
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestSetupAnalysisDirsCommand')
        # Store original location so we can get back at the end
        self.pwd = os.getcwd()
        # Move to working dir
        os.chdir(self.dirn)

    def tearDown(self):
        # Return to original dir
        os.chdir(self.pwd)
        # Remove the temporary test directory
        def del_rw(action,name,excinfo):
            # Explicitly remove read only files/
            # dirs
            os.chmod(os.path.dirname(name),0o755)
            os.chmod(name,0o655)
            os.remove(name)
        shutil.rmtree(self.dirn,onerror=del_rw)

    def test_setup_analysis_dirs(self):
        """
        setup_analysis_dirs: test create new analysis dirs
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Check project dirs don't exist
        for project in ("AB","CDE"):
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertFalse(os.path.exists(project_dir_path))
        # Add required metadata to 'projects.info'
        projects_info = os.path.join(mockdir.dirn,"projects.info")
        with open(projects_info,"w") as fp:
            fp.write(
"""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1,AB2\tAlan Brown\tRNA-seq\t.\tHuman\tAudrey Benson\t1% PhiX
CDE\tCDE3,CDE4\tClive David Edwards\tChIP-seq\t.\tMouse\tClaudia Divine Eccleston\t1% PhiX
""")
        # Expected data
        projects = {
            "AB": ["AB1_S1_R1_001.fastq.gz",
                   "AB1_S1_R2_001.fastq.gz",
                   "AB2_S2_R1_001.fastq.gz",
                   "AB2_S2_R2_001.fastq.gz"],
            "CDE": ["CDE3_S3_R1_001.fastq.gz",
                    "CDE3_S3_R2_001.fastq.gz",
                    "CDE4_S4_R1_001.fastq.gz",
                    "CDE4_S4_R2_001.fastq.gz"],
            "undetermined": ["Undetermined_S0_R1_001.fastq.gz",
                             "Undetermined_S0_R1_001.fastq.gz"]
        }
        # Setup the project dirs
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        setup_analysis_dirs(ap)
        # Check project dirs and contents
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertTrue(os.path.exists(project_dir_path))
            # Check README.info file
            readme_file = os.path.join(project_dir_path,
                                       "README.info")
            self.assertTrue(os.path.exists(readme_file))
            # Check Fastqs
            fastqs_dir = os.path.join(project_dir_path,
                                      "fastqs")
            self.assertTrue(os.path.exists(fastqs_dir))
            for fq in projects[project]:
                fastq = os.path.join(fastqs_dir,fq)
                self.assertTrue(os.path.exists(fastq))

    def test_setup_analysis_dirs_with_identifier(self):
        """
        setup_analysis_dirs: test create new analysis dirs with an identifier
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Check project dirs don't exist
        for project in ("AB","CDE"):
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertFalse(os.path.exists(project_dir_path))
        # Add required metadata to 'projects.info'
        projects_info = os.path.join(mockdir.dirn,"projects.info")
        with open(projects_info,"w") as fp:
            fp.write(
"""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1,AB2\tAlan Brown\tRNA-seq\t.\tHuman\tAudrey Benson\t1% PhiX
CDE\tCDE3,CDE4\tClive David Edwards\tChIP-seq\t.\tMouse\tClaudia Divine Eccleston\t1% PhiX
""")
        # Expected data
        projects = {
            "AB": ["AB1_S1_R1_001.fastq.gz",
                   "AB1_S1_R2_001.fastq.gz",
                   "AB2_S2_R1_001.fastq.gz",
                   "AB2_S2_R2_001.fastq.gz"],
            "CDE": ["CDE3_S3_R1_001.fastq.gz",
                    "CDE3_S3_R2_001.fastq.gz",
                    "CDE4_S4_R1_001.fastq.gz",
                    "CDE4_S4_R2_001.fastq.gz"],
            "undetermined": ["Undetermined_S0_R1_001.fastq.gz",
                             "Undetermined_S0_R1_001.fastq.gz"]
        }
        # Setup the project dirs
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        setup_analysis_dirs(ap,name="test_id")
        # Check project dirs and contents
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,
                                            "%s_%s" % (project,
                                                       "test_id"))
            self.assertTrue(os.path.exists(project_dir_path))
            # Check README.info file
            readme_file = os.path.join(project_dir_path,
                                       "README.info")
            self.assertTrue(os.path.exists(readme_file))
            # Check Fastqs
            fastqs_dir = os.path.join(project_dir_path,
                                      "fastqs")
            self.assertTrue(os.path.exists(fastqs_dir))
            for fq in projects[project]:
                fastq = os.path.join(fastqs_dir,fq)
                self.assertTrue(os.path.exists(fastq))

    def test_setup_analysis_dirs_ignore_commented_projects(self):
        """
        setup_analysis_dirs: ignore commented line in 'projects.info'
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Check project dirs don't exist
        for project in ("AB","CDE"):
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertFalse(os.path.exists(project_dir_path))
        # Add required metadata to 'projects.info'
        projects_info = os.path.join(mockdir.dirn,"projects.info")
        with open(projects_info,"w") as fp:
            fp.write(
"""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
#AB\tAB1,AB2\tAlan Brown\tRNA-seq\t.\tHuman\tAudrey Benson\t1% PhiX
CDE\tCDE3,CDE4\tClive David Edwards\tChIP-seq\t.\tMouse\tClaudia Divine Eccleston\t1% PhiX
""")
        # Expected data
        projects = {
            "CDE": ["CDE3_S3_R1_001.fastq.gz",
                    "CDE3_S3_R2_001.fastq.gz",
                    "CDE4_S4_R1_001.fastq.gz",
                    "CDE4_S4_R2_001.fastq.gz"],
            "undetermined": ["Undetermined_S0_R1_001.fastq.gz",
                             "Undetermined_S0_R1_001.fastq.gz"]
        }
        # Setup the project dirs
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        setup_analysis_dirs(ap)
        # Check project dirs and contents
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertTrue(os.path.exists(project_dir_path))
            # Check README.info file
            readme_file = os.path.join(project_dir_path,
                                       "README.info")
            self.assertTrue(os.path.exists(readme_file))
            # Check Fastqs
            fastqs_dir = os.path.join(project_dir_path,
                                      "fastqs")
            self.assertTrue(os.path.exists(fastqs_dir))
            for fq in projects[project]:
                fastq = os.path.join(fastqs_dir,fq)
                self.assertTrue(os.path.exists(fastq))
        # Check that 'AB' and '#AB' projects weren't created
        for project in ('AB','#AB'):
            self.assertFalse(os.path.exists(
                os.path.join(mockdir.dirn,project)))

    def test_setup_analysis_dirs_icell8_atac(self):
        """
        setup_analysis_dirs: test create new analysis dir for ICELL8 ATAC
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            paired_end=True,
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Add required metadata to 'projects.info'
        projects_info = os.path.join(mockdir.dirn,"projects.info")
        with open(projects_info,"w") as fp:
            fp.write(
"""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1,AB2\tAlan Brown\tscATAC-seq\tICELL8 ATAC\tHuman\tAudrey Benson\t1% PhiX
""")
        # Add ICELL8 ATAC outputs
        xlsx_file = os.path.join(mockdir.dirn,
                                 "bcl2fastq",
                                 "Reports",
                                 "icell8_atac_stats.xlsx")
        with open(xlsx_file,'w') as fp:
            fp.write("")
        # Expected data
        projects = {
            "AB": ["AB1_S1_R1_001.fastq.gz",
                   "AB1_S1_R2_001.fastq.gz",
                   "AB2_S2_R1_001.fastq.gz",
                   "AB2_S2_R2_001.fastq.gz"],
            "undetermined": ["Undetermined_S0_R1_001.fastq.gz",
                             "Undetermined_S0_R2_001.fastq.gz"]
        }
        # Check project dirs don't exist
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertFalse(os.path.exists(project_dir_path))
        # Setup the project dirs
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        setup_analysis_dirs(ap)
        # Check project dirs and contents
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertTrue(os.path.exists(project_dir_path))
            # Check README.info file
            readme_file = os.path.join(project_dir_path,
                                       "README.info")
            self.assertTrue(os.path.exists(readme_file))
            # Check Fastqs
            fastqs_dir = os.path.join(project_dir_path,
                                      "fastqs")
            self.assertTrue(os.path.exists(fastqs_dir))
            for fq in projects[project]:
                fastq = os.path.join(fastqs_dir,fq)
                self.assertTrue(os.path.exists(fastq))
        # Check extra data for ICELL8 ATAC
        icell8_atac_xlsx = os.path.join(mockdir.dirn,
                                        "AB",
                                        "icell8_atac_stats.xlsx")
        self.assertTrue(os.path.exists(icell8_atac_xlsx))

    def test_setup_analysis_dirs_10x_visium(self):
        """
        setup_analysis_dirs: test create new analysis dir for 10x Visium
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            reads=('R1','R2','R3','I1'),
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq")))
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq","AB")))
        # Add required metadata to 'projects.info'
        projects_info = os.path.join(mockdir.dirn,"projects.info")
        with open(projects_info,"w") as fp:
            fp.write(
"""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1,AB2\tAlan Brown\tFFPE Spatial RNA-seq\t10xGenomics Visium\tHuman\tAudrey Benson\t1% PhiX
""")
        # Expected data
        projects = {
            "AB": ["AB1_S1_R1_001.fastq.gz",
                   "AB1_S1_R2_001.fastq.gz",
                   "AB1_S1_R3_001.fastq.gz",
                   "AB1_S1_I1_001.fastq.gz",
                   "AB2_S2_R1_001.fastq.gz",
                   "AB2_S2_R2_001.fastq.gz",
                   "AB2_S2_R3_001.fastq.gz",
                   "AB2_S2_I1_001.fastq.gz"],
            "undetermined": ["Undetermined_S0_R1_001.fastq.gz",]
        }
        # Check project dirs don't exist
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertFalse(os.path.exists(project_dir_path))
        # Setup the project dirs
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        setup_analysis_dirs(ap)
        # Check project dirs and contents
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertTrue(os.path.exists(project_dir_path))
            # Check README.info file
            readme_file = os.path.join(project_dir_path,
                                       "README.info")
            self.assertTrue(os.path.exists(readme_file))
            # Check Fastqs
            fastqs_dir = os.path.join(project_dir_path,
                                      "fastqs")
            self.assertTrue(os.path.exists(fastqs_dir))
            for fq in projects[project]:
                fastq = os.path.join(fastqs_dir,fq)
                self.assertTrue(os.path.exists(fastq))
        # Check that Visium project includes directory
        # for images
        self.assertTrue(os.path.isdir(os.path.join(mockdir.dirn,
                                                   "AB",
                                                   "Visium_images")))

    def test_setup_analysis_dirs_10x_cytassist_visium(self):
        """
        setup_analysis_dirs: test create new analysis dir for 10x CytAssist Visium
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            reads=('R1','R2','R3','I1'),
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq")))
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq","AB")))
        # Add required metadata to 'projects.info'
        projects_info = os.path.join(mockdir.dirn,"projects.info")
        with open(projects_info,"w") as fp:
            fp.write(
"""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1,AB2\tAlan Brown\tFFPE Spatial RNA-seq\t10xGenomics CytAssist Visium\tHuman\tAudrey Benson\t1% PhiX
""")
        # Expected data
        projects = {
            "AB": ["AB1_S1_R1_001.fastq.gz",
                   "AB1_S1_R2_001.fastq.gz",
                   "AB1_S1_R3_001.fastq.gz",
                   "AB1_S1_I1_001.fastq.gz",
                   "AB2_S2_R1_001.fastq.gz",
                   "AB2_S2_R2_001.fastq.gz",
                   "AB2_S2_R3_001.fastq.gz",
                   "AB2_S2_I1_001.fastq.gz"],
            "undetermined": ["Undetermined_S0_R1_001.fastq.gz",]
        }
        # Check project dirs don't exist
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertFalse(os.path.exists(project_dir_path))
        # Setup the project dirs
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        setup_analysis_dirs(ap)
        # Check project dirs and contents
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertTrue(os.path.exists(project_dir_path))
            # Check README.info file
            readme_file = os.path.join(project_dir_path,
                                       "README.info")
            self.assertTrue(os.path.exists(readme_file))
            # Check Fastqs
            fastqs_dir = os.path.join(project_dir_path,
                                      "fastqs")
            self.assertTrue(os.path.exists(fastqs_dir))
            for fq in projects[project]:
                fastq = os.path.join(fastqs_dir,fq)
                self.assertTrue(os.path.exists(fastq))
        # Check that Visium project includes directory
        # for images
        self.assertTrue(os.path.isdir(os.path.join(mockdir.dirn,
                                                   "AB",
                                                   "Visium_images")))

    def test_setup_analysis_dirs_10x_multiome(self):
        """
        setup_analysis_dirs: test create new analysis dir for 10x Multiome
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            reads=('R1','R2','I1'),
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq")))
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq","AB")))
        # Add required metadata to 'projects.info'
        projects_info = os.path.join(mockdir.dirn,"projects.info")
        with open(projects_info,"w") as fp:
            fp.write(
"""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1,AB2\tAlan Brown\tGEX\t10xGenomics Single Cell Multiome\tHuman\tAudrey Benson\t1% PhiX
""")
        # Expected data
        projects = {
            "AB": ["AB1_S1_R1_001.fastq.gz",
                   "AB1_S1_R2_001.fastq.gz",
                   "AB1_S1_I1_001.fastq.gz",
                   "AB2_S2_R1_001.fastq.gz",
                   "AB2_S2_R2_001.fastq.gz",
                   "AB2_S2_I1_001.fastq.gz"],
            "undetermined": ["Undetermined_S0_R1_001.fastq.gz",]
        }
        # Check project dirs don't exist
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertFalse(os.path.exists(project_dir_path))
        # Setup the project dirs
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        setup_analysis_dirs(ap)
        # Check project dirs and contents
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertTrue(os.path.exists(project_dir_path))
            # Check README.info file
            readme_file = os.path.join(project_dir_path,
                                       "README.info")
            self.assertTrue(os.path.exists(readme_file))
            # Check Fastqs
            fastqs_dir = os.path.join(project_dir_path,
                                      "fastqs")
            self.assertTrue(os.path.exists(fastqs_dir))
            for fq in projects[project]:
                fastq = os.path.join(fastqs_dir,fq)
                self.assertTrue(os.path.exists(fastq))
            # Check template 10x_multiome_libraries.info
            template_libraries_file = os.path.join(
                project_dir_path,
                "10x_multiome_libraries.info.template")
            if project == "AB":
                # Template file should exist
                self.assertTrue(os.path.exists(template_libraries_file),
                                "Missing %s" % template_libraries_file)
            else:
                # No template file
                self.assertFalse(os.path.exists(template_libraries_file),
                                 "Found %s" % template_libraries_file)

    def test_setup_analysis_dirs_10x_cellplex_710(self):
        """
        setup_analysis_dirs: create new analysis dir for 10x CellPlex (7.1.0)
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901",
                "processing_software": {
                    "cellranger" : (
                        "/usr/local/cellranger/7.1.0/cellranger",
                        "cellranger",
                        "7.1.0"
                    )
                }
            },
            reads=('R1','R2','I1','I2'),
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq")))
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq","AB")))
        # Add required metadata to 'projects.info'
        projects_info = os.path.join(mockdir.dirn,"projects.info")
        with open(projects_info,"w") as fp:
            fp.write(
"""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1,AB2\tAlan Brown\tCellPlex scRNA-seq\t10xGenomics Chromium 3'v3\tHuman\tAudrey Benson\t1% PhiX
""")
        # Expected data
        projects = {
            "AB": ["AB1_S1_R1_001.fastq.gz",
                   "AB1_S1_R2_001.fastq.gz",
                   "AB1_S1_I1_001.fastq.gz",
                   "AB1_S1_I2_001.fastq.gz",
                   "AB2_S2_R1_001.fastq.gz",
                   "AB2_S2_R2_001.fastq.gz",
                   "AB2_S2_I1_001.fastq.gz",
                   "AB2_S2_I2_001.fastq.gz"],
            "undetermined": ["Undetermined_S0_R1_001.fastq.gz",]
        }
        # Check project dirs don't exist
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertFalse(os.path.exists(project_dir_path))
        # Setup the project dirs
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        setup_analysis_dirs(ap)
        # Check project dirs and contents
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertTrue(os.path.exists(project_dir_path))
            # Check README.info file
            readme_file = os.path.join(project_dir_path,
                                       "README.info")
            self.assertTrue(os.path.exists(readme_file))
            # Check Fastqs
            fastqs_dir = os.path.join(project_dir_path,
                                      "fastqs")
            self.assertTrue(os.path.exists(fastqs_dir))
            for fq in projects[project]:
                fastq = os.path.join(fastqs_dir,fq)
                self.assertTrue(os.path.exists(fastq))
            # Check template 10x_multi_config.csv
            template_multi_config_file = os.path.join(
                project_dir_path,
                "10x_multi_config.csv.template")
            if project == "AB":
                # Template file should exist
                self.assertTrue(os.path.exists(template_multi_config_file),
                                "Missing %s" % template_multi_config_file)
                # Template should contain '#no-bam'
                with open(template_multi_config_file,'rt') as fp:
                    template = fp.read()
                    self.assertTrue(template.find("\n#no-bam,") != -1)
            else:
                # No template file
                self.assertFalse(os.path.exists(template_multi_config_file),
                                 "Found %s" % template_multi_config_file)

    def test_setup_analysis_dirs_10x_cellplex_800(self):
        """
        setup_analysis_dirs: create new analysis dir for 10x CellPlex (8.0.0)
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901",
                "processing_software": {
                    "cellranger" : (
                        "/usr/local/cellranger/8.0.0/cellranger",
                        "cellranger",
                        "8.0.0"
                    )
                }
            },
            reads=('R1','R2','I1','I2'),
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq")))
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq","AB")))
        # Add required metadata to 'projects.info'
        projects_info = os.path.join(mockdir.dirn,"projects.info")
        with open(projects_info,"w") as fp:
            fp.write(
"""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1,AB2\tAlan Brown\tCellPlex scRNA-seq\t10xGenomics Chromium 3'v3\tHuman\tAudrey Benson\t1% PhiX
""")
        # Expected data
        projects = {
            "AB": ["AB1_S1_R1_001.fastq.gz",
                   "AB1_S1_R2_001.fastq.gz",
                   "AB1_S1_I1_001.fastq.gz",
                   "AB1_S1_I2_001.fastq.gz",
                   "AB2_S2_R1_001.fastq.gz",
                   "AB2_S2_R2_001.fastq.gz",
                   "AB2_S2_I1_001.fastq.gz",
                   "AB2_S2_I2_001.fastq.gz"],
            "undetermined": ["Undetermined_S0_R1_001.fastq.gz",]
        }
        # Check project dirs don't exist
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertFalse(os.path.exists(project_dir_path))
        # Setup the project dirs
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        setup_analysis_dirs(ap)
        # Check project dirs and contents
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertTrue(os.path.exists(project_dir_path))
            # Check README.info file
            readme_file = os.path.join(project_dir_path,
                                       "README.info")
            self.assertTrue(os.path.exists(readme_file))
            # Check Fastqs
            fastqs_dir = os.path.join(project_dir_path,
                                      "fastqs")
            self.assertTrue(os.path.exists(fastqs_dir))
            for fq in projects[project]:
                fastq = os.path.join(fastqs_dir,fq)
                self.assertTrue(os.path.exists(fastq))
            # Check template 10x_multi_config.csv
            template_multi_config_file = os.path.join(
                project_dir_path,
                "10x_multi_config.csv.template")
            if project == "AB":
                # Template file should exist
                self.assertTrue(os.path.exists(template_multi_config_file),
                                "Missing %s" % template_multi_config_file)
                # Template should contain 'create-bam'
                with open(template_multi_config_file,'rt') as fp:
                    template = fp.read()
                    self.assertTrue(template.find("\ncreate-bam,") != -1)
            else:
                # No template file
                self.assertFalse(os.path.exists(template_multi_config_file),
                                 "Found %s" % template_multi_config_file)

    def test_setup_analysis_dirs_10x_flex_710(self):
        """
        setup_analysis_dirs: create new analysis dir for 10x Flex (7.1.0)
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={
                "instrument_datestamp": "170901",
                "processing_software": {
                    "cellranger" : (
                        "/usr/local/cellranger/7.1.0/cellranger",
                        "cellranger",
                        "7.1.0"
                    )
                }
            },
            reads=('R1','R2','I1','I2'),
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq")))
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq","AB")))
        # Add required metadata to 'projects.info'
        projects_info = os.path.join(mockdir.dirn,"projects.info")
        with open(projects_info,"w") as fp:
            fp.write(
"""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1\tAlan Brown\tFlex\t10xGenomics Chromium 3'v3\tHuman\tAudrey Benson\t1% PhiX
""")
        # Expected data
        projects = {
            "AB": ["AB1_S1_R1_001.fastq.gz",
                   "AB1_S1_R2_001.fastq.gz",
                   "AB1_S1_I1_001.fastq.gz",
                   "AB1_S1_I2_001.fastq.gz",],
            "undetermined": ["Undetermined_S0_R1_001.fastq.gz",]
        }
        # Check project dirs don't exist
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertFalse(os.path.exists(project_dir_path))
        # Setup the project dirs
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        setup_analysis_dirs(ap)
        # Check project dirs and contents
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertTrue(os.path.exists(project_dir_path))
            # Check README.info file
            readme_file = os.path.join(project_dir_path,
                                       "README.info")
            self.assertTrue(os.path.exists(readme_file))
            # Check Fastqs
            fastqs_dir = os.path.join(project_dir_path,
                                      "fastqs")
            self.assertTrue(os.path.exists(fastqs_dir))
            for fq in projects[project]:
                fastq = os.path.join(fastqs_dir,fq)
                self.assertTrue(os.path.exists(fastq))
            # Check template 10x_multi_config.csv
            template_multi_config_file = os.path.join(
                project_dir_path,
                "10x_multi_config.csv.template")
            if project == "AB":
                # Template file should exist
                self.assertTrue(os.path.exists(template_multi_config_file),
                                "Missing %s" % template_multi_config_file)
                # Template should contain '#no-bam'
                with open(template_multi_config_file,'rt') as fp:
                    template = fp.read()
                    self.assertTrue(template.find("\n#no-bam,") != -1)
            else:
                # No template file
                self.assertFalse(os.path.exists(template_multi_config_file),
                                 "Found %s" % template_multi_config_file)

    def test_setup_analysis_dirs_10x_flex_800(self):
        """
        setup_analysis_dirs: create new analysis dir for 10x Flex (8.0.0)
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={
                "instrument_datestamp": "170901",
                "processing_software": {
                    "cellranger" : (
                        "/usr/local/cellranger/8.0.0/cellranger",
                        "cellranger",
                        "8.0.0"
                    )
                }
            },
            reads=('R1','R2','I1','I2'),
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq")))
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq","AB")))
        # Add required metadata to 'projects.info'
        projects_info = os.path.join(mockdir.dirn,"projects.info")
        with open(projects_info,"w") as fp:
            fp.write(
"""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1\tAlan Brown\tFlex\t10xGenomics Chromium 3'v3\tHuman\tAudrey Benson\t1% PhiX
""")
        # Expected data
        projects = {
            "AB": ["AB1_S1_R1_001.fastq.gz",
                   "AB1_S1_R2_001.fastq.gz",
                   "AB1_S1_I1_001.fastq.gz",
                   "AB1_S1_I2_001.fastq.gz",],
            "undetermined": ["Undetermined_S0_R1_001.fastq.gz",]
        }
        # Check project dirs don't exist
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertFalse(os.path.exists(project_dir_path))
        # Setup the project dirs
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        setup_analysis_dirs(ap)
        # Check project dirs and contents
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertTrue(os.path.exists(project_dir_path))
            # Check README.info file
            readme_file = os.path.join(project_dir_path,
                                       "README.info")
            self.assertTrue(os.path.exists(readme_file))
            # Check Fastqs
            fastqs_dir = os.path.join(project_dir_path,
                                      "fastqs")
            self.assertTrue(os.path.exists(fastqs_dir))
            for fq in projects[project]:
                fastq = os.path.join(fastqs_dir,fq)
                self.assertTrue(os.path.exists(fastq))
            # Check template 10x_multi_config.csv
            template_multi_config_file = os.path.join(
                project_dir_path,
                "10x_multi_config.csv.template")
            if project == "AB":
                # Template file should exist
                self.assertTrue(os.path.exists(template_multi_config_file),
                                "Missing %s" % template_multi_config_file)
                # Template should contain 'create-bam'
                with open(template_multi_config_file,'rt') as fp:
                    template = fp.read()
                    self.assertTrue(template.find("\ncreate-bam,") != -1)
            else:
                # No template file
                self.assertFalse(os.path.exists(template_multi_config_file),
                                 "Found %s" % template_multi_config_file)

    def test_setup_analysis_dirs_10x_immune_profiling_710(self):
        """
        setup_analysis_dirs: create new analysis dir for 10x Single Cell Immune Profiling (7.1.0)
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901",
                "processing_software": {
                    "cellranger" : (
                        "/usr/local/cellranger/7.1.0/cellranger",
                        "cellranger",
                        "7.1.0"
                    )
                }
            },
            reads=('R1','R2','I1','I2'),
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq")))
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq","AB")))
        # Add required metadata to 'projects.info'
        projects_info = os.path.join(mockdir.dirn,"projects.info")
        with open(projects_info,"w") as fp:
            fp.write(
"""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1\tAlan Brown\tSingle Cell Immune Profiling\t10xGenomics Chromium 5'\tMouse\tAudrey Benson\t1% PhiX
""")
        # Expected data
        projects = {
            "AB": ["AB1_S1_R1_001.fastq.gz",
                   "AB1_S1_R2_001.fastq.gz",
                   "AB1_S1_I1_001.fastq.gz",
                   "AB1_S1_I2_001.fastq.gz",],
            "undetermined": ["Undetermined_S0_R1_001.fastq.gz",]
        }
        # Check project dirs don't exist
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertFalse(os.path.exists(project_dir_path))
        # Setup the project dirs
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        setup_analysis_dirs(ap)
        # Check project dirs and contents
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertTrue(os.path.exists(project_dir_path))
            # Check README.info file
            readme_file = os.path.join(project_dir_path,
                                       "README.info")
            self.assertTrue(os.path.exists(readme_file))
            # Check Fastqs
            fastqs_dir = os.path.join(project_dir_path,
                                      "fastqs")
            self.assertTrue(os.path.exists(fastqs_dir))
            for fq in projects[project]:
                fastq = os.path.join(fastqs_dir,fq)
                self.assertTrue(os.path.exists(fastq))
            # Check template 10x_multi_config.csv
            template_multi_config_file = os.path.join(
                project_dir_path,
                "10x_multi_config.csv.template")
            if project == "AB":
                # Template file should exist
                self.assertTrue(os.path.exists(template_multi_config_file),
                                "Missing %s" % template_multi_config_file)
                # Template should contain '#no-bam'
                with open(template_multi_config_file,'rt') as fp:
                    template = fp.read()
                    self.assertTrue(template.find("\n#no-bam,") != -1)
            else:
                # No template file
                self.assertFalse(os.path.exists(template_multi_config_file),
                                 "Found %s" % template_multi_config_file)

    def test_setup_analysis_dirs_10x_immune_profiling_800(self):
        """
        setup_analysis_dirs: create new analysis dir for 10x Single Cell Immune Profiling (8.0.0)
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901",
                "processing_software": {
                    "cellranger" : (
                        "/usr/local/cellranger/8.0.0/cellranger",
                        "cellranger",
                        "8.0.0"
                    )
                }
            },
            reads=('R1','R2','I1','I2'),
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq")))
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq","AB")))
        # Add required metadata to 'projects.info'
        projects_info = os.path.join(mockdir.dirn,"projects.info")
        with open(projects_info,"w") as fp:
            fp.write(
"""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1\tAlan Brown\tSingle Cell Immune Profiling\t10xGenomics Chromium 5'\tMouse\tAudrey Benson\t1% PhiX
""")
        # Expected data
        projects = {
            "AB": ["AB1_S1_R1_001.fastq.gz",
                   "AB1_S1_R2_001.fastq.gz",
                   "AB1_S1_I1_001.fastq.gz",
                   "AB1_S1_I2_001.fastq.gz",],
            "undetermined": ["Undetermined_S0_R1_001.fastq.gz",]
        }
        # Check project dirs don't exist
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertFalse(os.path.exists(project_dir_path))
        # Setup the project dirs
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        setup_analysis_dirs(ap)
        # Check project dirs and contents
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertTrue(os.path.exists(project_dir_path))
            # Check README.info file
            readme_file = os.path.join(project_dir_path,
                                       "README.info")
            self.assertTrue(os.path.exists(readme_file))
            # Check Fastqs
            fastqs_dir = os.path.join(project_dir_path,
                                      "fastqs")
            self.assertTrue(os.path.exists(fastqs_dir))
            for fq in projects[project]:
                fastq = os.path.join(fastqs_dir,fq)
                self.assertTrue(os.path.exists(fastq))
            # Check template 10x_multi_config.csv
            template_multi_config_file = os.path.join(
                project_dir_path,
                "10x_multi_config.csv.template")
            if project == "AB":
                # Template file should exist
                self.assertTrue(os.path.exists(template_multi_config_file),
                                "Missing %s" % template_multi_config_file)
                # Template should contain 'create-bam'
                with open(template_multi_config_file,'rt') as fp:
                    template = fp.read()
                    self.assertTrue(template.find("\ncreate-bam,") != -1)
            else:
                # No template file
                self.assertFalse(os.path.exists(template_multi_config_file),
                                 "Found %s" % template_multi_config_file)

    def test_setup_analysis_dirs_missing_metadata(self):
        """
        setup_analysis_dirs: raise exception if metadata not set
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Attempt to set up the project dirs
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        self.assertRaises(Exception,
                          setup_analysis_dirs,ap)

    def test_setup_analysis_dirs_ignore_missing_metadata(self):
        """
        setup_analysis_dirs: test create new analysis dirs (ignore missing metadata)
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Expected data
        projects = {
            "AB": ["AB1_S1_R1_001.fastq.gz",
                   "AB1_S1_R2_001.fastq.gz",
                   "AB2_S2_R1_001.fastq.gz",
                   "AB2_S2_R2_001.fastq.gz"],
            "CDE": ["CDE3_S3_R1_001.fastq.gz",
                    "CDE3_S3_R2_001.fastq.gz",
                    "CDE4_S4_R1_001.fastq.gz",
                    "CDE4_S4_R2_001.fastq.gz"],
            "undetermined": ["Undetermined_S0_R1_001.fastq.gz",
                             "Undetermined_S0_R1_001.fastq.gz"]
        }
        # Check project dirs don't exist
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertFalse(os.path.exists(project_dir_path))
        # Setup the project dirs
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        setup_analysis_dirs(ap,ignore_missing_metadata=True)
        # Check project dirs and contents
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertTrue(os.path.exists(project_dir_path))
            # Check README.info file
            readme_file = os.path.join(project_dir_path,
                                       "README.info")
            self.assertTrue(os.path.exists(readme_file))
            # Check Fastqs
            fastqs_dir = os.path.join(project_dir_path,
                                      "fastqs")
            self.assertTrue(os.path.exists(fastqs_dir))
            for fq in projects[project]:
                fastq = os.path.join(fastqs_dir,fq)
                self.assertTrue(os.path.exists(fastq))

    def test_setup_analysis_dirs_bad_single_cell_platform(self):
        """
        setup_analysis_dirs: raise exception for bad single cell platform
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Add required metadata to 'projects.info'
        projects_info = os.path.join(mockdir.dirn,"projects.info")
        with open(projects_info,"w") as fp:
            fp.write(
"""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1,AB2\tAlan Brown\tRNA-seq\t.\tHuman\tAudrey Benson\t1% PhiX
CDE\tCDE3,CDE4\tClive David Edwards\tscRNA-seq\t11xGenomics Chromium 5'v0\tMouse\tClaudia Divine Eccleston\t1% PhiX
""")
        # Attempt to set up the project dirs
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        self.assertRaises(Exception,
                          setup_analysis_dirs,ap)
