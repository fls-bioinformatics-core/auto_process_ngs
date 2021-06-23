#!/usr/bin/env python
#
#     tenx_genomics_utils.py: utility functions for handling 10xGenomics data
#     Copyright (C) University of Manchester 2017-2021 Peter Briggs
#

"""
tenx_genomics_utils.py

Utility classes and functions for processing the outputs from 10xGenomics
platforms:

- MetricSummary
- GexSummary
- AtacSummary
- MultiomeSummary
- MultiomeLibraries
- flow_cell_id
- has_10x_indices
- has_chromium_sc_indices
- cellranger_info
- spaceranger_info
- make_qc_summary_html
- set_cell_count_for_project
- add_cellranger_args
"""

#######################################################################
# Imports
#######################################################################

import os
import re
import json
import shutil
from builtins import range
from bcftbx.IlluminaData import SampleSheet
from bcftbx.IlluminaData import IlluminaData
from bcftbx.IlluminaData import IlluminaDataError
from bcftbx.IlluminaData import split_run_name_full
from bcftbx.JobRunner import SimpleJobRunner
from bcftbx.TabFile import TabFile
from bcftbx.utils import mkdirs
from bcftbx.utils import find_program
from .analysis import locate_project
from .analysis import split_sample_reference
from .command import Command
from .docwriter import Document
from .docwriter import List
from .docwriter import Link
from .docwriter import Table
from .analysis import AnalysisProject
from .bcl2fastq.utils import get_bases_mask
from .metadata import AnalysisProjectQCDirInfo
from .utils import get_numbered_subdir
from .utils import ZipArchive
from . import css_rules

# Initialise logging
import logging
logger = logging.getLogger(__name__)

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
    'SC5P-PE': 'Single Cell 5\' paired-end (both R1 and R2 are used for alignment)',
    'SC5P-R2': 'Single Cell 5\' R2-only (where only R2 is used for alignment)',
}

#######################################################################
# Classes
#######################################################################

class MetricsSummary(TabFile):
    """
    Base class for extracting data from cellranger* count
    *summary.csv files

    The files consists of two lines: the first is a
    header line, the second consists of corresponding
    data values.

    In addition: in some variants (e.g.
    'metrics_summary.csv'), integer data values are formatted
    to use commas to separate thousands (e.g. 2,272) and
    values which contain commas are enclosed in double
    quotes.

    For example:

    Estimated Number of Cells,Mean Reads per Cell,...
    "2,272","107,875","1,282","245,093,084",98.3%,...

    This class extracts the data values and where
    possible converts them to integers.
    """
    def __init__(self,f):
        """
        Create a new MetricsSummary instance

        Arguments:
          f (str): path to the 'metrics_summary.csv' file
        """
        # Read in data from the file
        with open(f,'rt') as fp:
            s = fp.read()
        self._data = dict()
        s = s.strip().split('\n')
        if len(s) != 2:
            raise Exception("%s: MetricsSummary expects 2 lines"
                            % f)
        # Set up the tabfile instance
        TabFile.__init__(self,
                         column_names=self._tokenise(s[0]),
                         delimiter=',')
        # Add the data
        self.append(data=self._tokenise(s[1]))
    def _tokenise(self,line):
        """
        Internal: process line from *summary.csv

        Arguments:
          line (str): line to tokenise

        Returns:
          List: list of tokens extracted from the
            supplied line, with enclosing quotes
            removed and converted to integers
            where possible.
        """
        # Split the line into tokens
        # (taking account of quotes)
        tokens = list()
        this_token = ''
        quoted = False
        for c in line:
            if c == '"':
                # Start or end of a quote
                quoted = (not quoted)
            elif c == ',':
                if not quoted:
                    # End of a token
                    tokens.append(this_token)
                    this_token = ''
                    continue
            # Append to current token
            this_token += c
        # Deal with the last token
        tokens.append(this_token)
        # Strip quotes
        tokens = [t.strip('"') for t in tokens]
        # Convert to integer where possible
        # (i.e. remove commas from e.g. "2,272")
        for i in range(len(tokens)):
            try:
                tokens[i] = int(tokens[i].replace(',',''))
            except ValueError:
                pass
        return tokens
    def fetch(self,field):
        """
        Fetch data associated with an arbitrary field
        """
        return self[0][field]

class GexSummary(MetricsSummary):
    """
    Extract data from metrics_summary.csv file for scRNA-seq

    Utility class for extracting data from a
    'metrics_summary.csv' file output from running
    'cellranger count'.

    The file consists of two lines: the first is a
    header line, the second consists of corresponding
    data values.

    The following properties are available:

    - estimated_number_of_cells
    - mean_reads_per_cell
    - median_genes_per_cell
    - frac_reads_in_cells
    """
    def __init__(self,f):
        """
        Create a new GexSummary instance

        Arguments:
          f (str): path to the 'metrics_summary.csv' file
        """
        MetricsSummary.__init__(self,f)
    @property
    def estimated_number_of_cells(self):
        """
        Return the estimated number of cells
        """
        return self.fetch('Estimated Number of Cells')
    @property
    def mean_reads_per_cell(self):
        """
        Return the mean reads per cell
        """
        return self.fetch('Mean Reads per Cell')
    @property
    def median_genes_per_cell(self):
        """
        Return the median genes per cell
        """
        return self.fetch('Median Genes per Cell')
    @property
    def frac_reads_in_cells(self):
        """
        Return the fraction of reads in cells
        """
        return self.fetch('Fraction Reads in Cells')

class AtacSummary(MetricsSummary):
    """
    Extract data from summary.csv file for scATAC-seq

    Utility class for extracting data from a 'summary.csv'
    file output from running 'cellranger-atac count'.

    The file consists of two lines: the first is a
    header line, the second consists of corresponding
    data values.

    The following properties are available:

    - cells_detected
    - annotated_cells
    - median_fragments_per_cell
    - frac_fragments_overlapping_targets
    """
    def __init__(self,f):
        """
        Create a new AtacSummary instance

        Arguments:
          f (str): path to the 'summary.csv' file
        """
        MetricsSummary.__init__(self,f)
    @property
    def cells_detected(self):
        """
        Return the number of cells detected

        Only supported for Cellranger ATAC < 2.0.0; raises
        AttributeError otherwise.
        """
        if self.version != "2.0.0":
            # Only supported for pre-2.0.0
            return self.fetch('cells_detected')
        else:
            # Not supported
            raise AttributeError
    @property
    def annotated_cells(self):
        """
        Return the number of annotated cells

        Only supported for Cellranger ATAC < 2.0.0; raises
        AttributeError otherwise.
        """
        if self.version != "2.0.0":
            # Only supported for pre-2.0.0
            return self.fetch('annotated_cells')
        else:
            # Not supported
            raise AttributeError
    @property
    def estimated_number_of_cells(self):
        """
        Return the estimated number of cells

        Only supported for Cellranger ATAC < 2.0.0; raises
        AttributeError otherwise.
        """
        if self.version == "2.0.0":
            # Only supported for 2.0.0
            return self.fetch('Estimated number of cells')
        else:
            # Not supported
            raise AttributeError
    @property
    def median_fragments_per_cell(self):
        """
        Return the median fragments per cell
        """
        try:
            return self.fetch('median_fragments_per_cell')
        except KeyError:
            return self.fetch('Median high-quality fragments per cell')
    @property
    def frac_fragments_overlapping_targets(self):
        """
        Return the fraction of fragments overlapping targets
        """
        try:
            return self.fetch('frac_fragments_overlapping_targets')
        except KeyError:
            return self.fetch('Fraction of high-quality fragments overlapping TSS')
    @property
    def frac_fragments_overlapping_peaks(self):
        """
        Return the fraction of fragments overlapping targets
        """
        try:
            return self.fetch('frac_fragments_overlapping_peaks')
        except KeyError:
            return self.fetch('Fraction of high-quality fragments overlapping peaks')
    @property
    def version(self):
        """
        Return the pipeline version
        """
        try:
            version = self.fetch('cellranger-atac_version')
        except KeyError:
            version = self.fetch('Pipeline version')
        if version.startswith('cellranger-atac-'):
            version = version[len('cellranger-atac-'):]
        return version

class MultiomeSummary(MetricsSummary):
    """
    Extract data from summary.csv file for multiome GEX-ATAC

    Utility class for extracting data from a 'summary.csv'
    file output from running 'cellranger-arc count'.

    The file consists of two lines: the first is a
    header line, the second consists of corresponding
    data values.

    The following properties are available:

    - estimated_number_of_cells
    - atac_median_high_quality_fragments_per_cell
    - gex_median_cells_per_gene
    """
    def __init__(self,f):
        """
        Create a new MultiomeSummary instance

        Arguments:
          f (str): path to the 'summary.csv' file
        """
        MetricsSummary.__init__(self,f)
    @property
    def estimated_number_of_cells(self):
        """
        Return the estimated number of cells
        """
        return self.fetch('Estimated number of cells')
    @property
    def atac_median_high_quality_fragments_per_cell(self):
        """
        Return the median high-quality fragments per cell
        for ATAC data
        """
        return self.fetch(
            'ATAC Median high-quality fragments per cell')
    @property
    def gex_median_cells_per_gene(self):
        """
        Return the median genes per cell for GEX data
        """
        return self.fetch('GEX Median genes per cell')

class MultiomeLibraries(object):
    """
    Class to handle '10x_multiome_libraries.info' files

    These files link sample names in an analysis project
    with those in another project. They consist of
    tab-delimited lines of the form:

    <LOCAL_SAMPLE>   <REMOTE_SAMPLE>

    where LOCAL_SAMPLE should be the name of a sample in
    the local project directory, and REMOTE_SAMPLE is a
    compound ID for a sample in a different project, of
    the form:

    [RUN:]PROJECT/SAMPLE_NAME

    In turn, RUN can be a run reference ID, a run name, or
    path to an analysis directory.
    """
    def __init__(self,filen):
        """
        Create new MultiomeLibraries instance

        Arguments:
          filen (str): path to a 10x Multiome libraries
            info file
        """
        self._filen = os.path.abspath(filen)
        self._samples = dict()
        # Read data in from libraries file
        with open(self._filen,'rt') as fp:
            for line in fp:
                # Ignore comment lines
                if line.startswith('#'):
                    continue
                # Lines should be 'local_sample<TAB>remote_sample'
                local_sample,remote_sample = line.split()
                try:
                    self._samples[local_sample].append(remote_sample)
                except KeyError:
                    self._samples[local_sample] = [remote_sample,]

    def _library_description(self,library_type):
        """
        Internal: get cellranger-arc description of library type
        """
        if library_type == 'ATAC':
            return "Chromatin Accessibility"
        elif library_type == 'GEX':
            return "Gene Expression"
        else:
            raise Exception("Unsupported library: '%s'"
                            % library_type)

    @property
    def local_samples(self):
        """
        List the local sample names
        """
        return sorted(list(self._samples.keys()))

    def linked_samples(self,sample):
        """
        List the remote samples associated with a local sample
        """
        return self._samples[sample]

    def linked_projects(self):
        """
        List the projects linked to this one
        """
        projects = list()
        for sample in self.local_samples:
            for linked_sample in self.linked_samples(sample):
                # Extract and locate the location of the remote sample
                run,project,remote_sample = split_sample_reference(
                    linked_sample)
                linked_project = locate_project(
                    "%s:%s" % (run,project),
                    start_dir=os.path.dirname(self._filen),
                    ascend=True)
                if not linked_project:
                    raise Exception("Failed to locate project for "
                                    "linked sample: %s"
                                    % linked_sample)
                else:
                    # Check linked project not already in the list
                    if linked_project.dirn not in \
                       [p.dirn for p in projects]:
                        projects.append(linked_project)
        # Return list of projects
        return sorted(projects,key=lambda p: p.name)

    def write_libraries_csv(self,sample,fastq_dir,library_type,
                            filen=None):
        """
        Create a cellranger-arc libraries.csv file

        The format is a header line of the form:

        fastqs,sample,library_type

        See https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/using/using/count#libraries

        Arguments:
          sample (str): local sample name
          fastq_dir (str): path to directory with the
            Fastqs for the sample
          library_type (str): name of the library type
            for the sample
          filen (str): optional, path to the output CSV
            file (defaults to 'libraries.SAMPLE.csv')
        """
        # Determine output file
        if filen is None:
            filen = "libraries.%s.csv" % sample
        filen = os.path.abspath(filen)
        # Get linked samples
        linked_samples = self._samples[sample]
        # Write libraries.csv file for cellranger-arc
        with open(filen,'wt') as fp:
            fp.write("fastqs,sample,library_type\n")
            # Data for local sample
            fp.write("%s,%s,%s\n" %
                     (fastq_dir,
                      sample,
                      self._library_description(library_type)))
            # Data for linked sample(s)
            for sample_id in linked_samples:
                print("Locating linked sample: '%s'" % sample_id)
                run,project,linked_sample = \
                    split_sample_reference(sample_id)
                project = locate_project(
                    "%s:%s" % (run,project),
                    start_dir=os.path.dirname(self._filen),
                    ascend=True)
                if not project:
                    raise Exception("Failed to locate project for "
                                    "linked sample: %s"
                                    % sample_id)
                fp.write("%s,%s,%s\n" %
                         (project.fastq_dir,
                          linked_sample,
                          self._library_description(
                              project.info.library_type)))
            print("Generated %s" % filen)

#######################################################################
# Functions
#######################################################################

def flow_cell_id(run_name):
    """
    Extract the flow cell ID from a run name

    For example for run name "170426_K00311_0033_AHJCY7BBXX"
    the extracted flow cell ID will be "HJCY7BBXX".

    Arguments:
      run_name (str): path to the run name to extract
        flow cell ID from

    Returns:
      String: the extracted flow cell ID.
    """
    ds,inst,run,prefix,flow_cell_id = split_run_name_full(run_name)
    return flow_cell_id

def has_10x_indices(sample_sheet):
    """
    Check if a sample sheet contains 10xGenomics-format indices

    The Chromium SC 3'v2 indices are of the form:

    SI-GA-[A-H][1-12]

    e.g. 'SI-GA-B11' (see
    https://support.10xgenomics.com/permalink/27rGqWvNYYuqkgeS66sksm)

    For scATAC-seq the indices are assumed to be of the form:

    SI-NA-[A-H][1-12]

    e.g. 'SI-NA-G9'

    For Visium data the indices are assumed to be of the form:

    SI-TT-[A-H][1-12]

    e.g. 'SI-TT-B1'

    Arguments:
      sample_sheet (str): path to the sample sheet CSV
        file to check

    Returns:
      Boolean: True if the sample sheet contains at least
        one 10xGenomics-style index, False if not.
    """
    index_pattern = re.compile(r"SI-(GA|NA|TT)-[A-H](1[0-2]|[1-9])$")
    s = SampleSheet(sample_sheet)
    for line in s:
        try:
            if index_pattern.match(line['index']):
                return True
        except KeyError:
            pass
    return False

def has_chromium_sc_indices(sample_sheet):
    """
    Wrapper for 'has_10x_indices'.

    Maintained for backwards compatibility
    """
    return has_10x_indices(sample_sheet)

def get_bases_mask_10x_atac(runinfo_xml):
    """
    Acquire a bases mask for 10xGenomics scATAC-seq

    Generates an initial bases mask based on the run
    contents, and then updates this so that:

    1. Only the first 8 bases of the first index read
       are actually used, and
    2. The second index read is converted to a data
       read.

    For example: if the initial bases mask is
    'Y50,I16,I16,Y50' then the scATAC-seq bases mask
    will be 'Y50,I8nnnnnnnn,Y16,Y50'.

    Arguments:
      runinfo_xml (str): path to the RunInfo.xml for
        the sequencing run

    Returns:
      String: 10xGenomics scATAC-seq bases mask string
    """
    bases_mask = get_bases_mask(runinfo_xml).lower().split(',')
    # Check there are four reads defined
    if len(bases_mask) != 4:
        raise Exception("Bases mask '%s' should have 4 reads "
                        "defined (has %d)" % (bases_mask,
                                              len(bases_mask)))
    # First read
    r1_mask = bases_mask[0]
    # Update first index to restrict to 8 bases
    num_cycles = int(bases_mask[1][1:])
    if num_cycles < 8:
        raise Exception("Index read < 8 bases")
    i1_mask = "I8%s" % ('n'*(num_cycles-8),)
    # Update second index to second read
    r2_mask = bases_mask[2].replace('i','y')
    # Keep last read as is
    r3_mask = bases_mask[3]
    # Reassemble and return
    return ','.join((r1_mask,i1_mask,r2_mask,r3_mask,))

def make_qc_summary_html(json_file,html_file):
    """
    Make HTML report for cellranger mkfastqs processing stats

    Arguments:
      json_file (str): path to JSON file output from
        cellranger mkfastq command
      html_file (str): path to output HTML file
    """
    # Load the JSON data
    with open(json_file,'r') as fp:
        data = json.load(fp)
    # Initialise the HTML report
    qc_summary = Document("mkfastq QC report")
    qc_summary.add_css_rule(css_rules.QC_REPORT_CSS_RULES)
    qc_summary.add_css_rule("table { font-size: 80%;\n"
                            "        font-family: sans-serif; }")
    qc_summary.add_css_rule("td { text-align: right; }")
    # Add table of contents
    toc = qc_summary.add_section("Contents",name="toc")
    toc_list = List()
    toc.add(toc_list)
    # General information
    general_info = qc_summary.add_section("General information")
    toc_list.add_item(Link("General information",general_info))
    data_items = ['run_id',
                  'experiment_name',
                  '10x_software_version',
                  'bcl2fastq_version',
                  'bcl2fastq_args',
                  'rta_version',]
    tbl = Table(columns=['Parameter','Value'])
    for item in data_items:
        tbl.add_row(Parameter=item,Value=data[item])
    general_info.add(tbl)
    # Get the sample names
    sample_names = list(data['sample_qc'].keys())
    # Get names of the associated data items
    sample0 = sample_names[0]
    item_names = data['sample_qc'][sample0]['all'].keys()
    # Report QC for each sample in tables
    for sample in sample_names:
        sample_qc = qc_summary.add_section("Sample: %s" % sample)
        toc_list.add_item(Link("Sample: %s" % sample,sample_qc))
        # Set up the table
        tbl = Table(['items',],items="")
        for item in item_names:
            tbl.add_row(items=item)
        # Lanes
        lanes = data['sample_qc'][sample].keys()
        for lane in lanes:
            column = "%s" % lane 
            tbl.append_columns(column)
            # Add the data
            for i,item in enumerate(item_names):
                tbl.set_value(i,column,data['sample_qc'][sample][lane][item])
        # Add to the document
        sample_qc.add(tbl)
    # Write the report
    qc_summary.write(html_file)

def set_cell_count_for_project(project_dir,qc_dir=None):
    """
    Set the total number of cells for a project

    Sums the number of cells for each sample in a project
    (as determined from 'cellranger count' and extracted
    from the 'metrics_summary.csv' file for scRNA-seq, or
    from 'cellranger-atac count' and extracted from the
    'summary.csv' file for scATAC-seq).

    The final count is written to the 'number_of_cells'
    metadata item for the project.

    Arguments:
      project_dir (str): path to the project directory
      qc_dir (str): path to QC directory (if not the default
        QC directory for the project)

    Returns:
      Integer: exit code, non-zero values indicate problems
        were encountered.
    """
    # Set up basic info
    project = AnalysisProject(os.path.basename(project_dir),
                              project_dir)
    if qc_dir is None:
        qc_dir = project.qc_dir
    qc_dir = os.path.abspath(qc_dir)
    # Determine which 10x pipeline was used
    pipeline = None
    single_cell_platform = project.info.single_cell_platform
    if single_cell_platform:
        if single_cell_platform.startswith("10xGenomics Chromium 3'"):
            pipeline = "cellranger"
        elif single_cell_platform == "10xGenomics Single Cell ATAC":
            pipeline = "cellranger-atac"
        elif single_cell_platform == "10xGenomics Single Cell Multiome":
            pipeline = "cellranger-arc"
    if not pipeline:
        raise NotImplementedError("Not implemented for platform '%s'"
                                  % single_cell_platform)
    # Determine possible locations for outputs
    count_dirs = []
    # New-style with version and reference data levels
    qc_info_file = os.path.join(qc_dir,"qc.info")
    if os.path.exists(qc_info_file):
        qc_info = AnalysisProjectQCDirInfo(filen=qc_info_file)
        try:
            cellranger_refdata = qc_info['cellranger_refdata']
        except KeyError:
            cellranger_refdata = None
        try:
            cellranger_version = qc_info['cellranger_version']
        except KeyError:
            cellranger_version = None
        if cellranger_version and cellranger_refdata:
            count_dirs.append(os.path.join(qc_dir,
                                           "cellranger_count",
                                           cellranger_version,
                                           os.path.basename(
                                               cellranger_refdata)))
    # Old-style without additional subdirectories
    count_dirs.append(os.path.join(qc_dir,
                                   "cellranger_count"))
    # Check each putative output location in turn
    for count_dir in count_dirs:
        # Check that the directory exists
        if not os.path.exists(count_dir):
            logger.warning("%s: not found" % count_dir)
            continue
        print("Looking for outputs under %s" % count_dir)
        # Loop over samples and collect cell numbers for each
        number_of_cells = 0
        try:
            for sample in project.samples:
                print("- %s" % sample)
                outs_dir = os.path.join(count_dir,
                                        sample.name,
                                        "outs")
                if pipeline == "cellranger":
                    # Single cell/nuclei RNA-seq output
                    metrics_summary_csv = os.path.join(
                        outs_dir,
                        "metrics_summary.csv")
                    if not os.path.exists(metrics_summary_csv):
                        raise Exception("Failed to add cell count "
                                        "for sample '%s': missing "
                                        "file '%s'" %
                                        (sample.name,
                                         metrics_summary_csv))
                    # Extract cell numbers
                    metrics = GexSummary(metrics_summary_csv)
                    number_of_cells += \
                        metrics.estimated_number_of_cells
                elif pipeline == "cellranger-atac":
                    # Single cell ATAC-seq output
                    summary_csv = os.path.join(
                        outs_dir,
                        "summary.csv")
                    if not os.path.exists(summary_csv):
                        raise Exception("Failed to add cell count "
                                        "for sample '%s': missing "
                                        "file '%s'" %
                                        (sample.name,
                                         summary_csv))
                    # Extract cell numbers
                    try:
                        # Cellranger ATAC pre-2.0.0 uses 'annotated_cells'
                        number_of_cells += AtacSummary(summary_csv).\
                                           annotated_cells
                    except AttributeError:
                        # Cellranger ATAC 2.0.0 uses 'Estimated number of
                        # cells'
                        number_of_cells += AtacSummary(summary_csv).\
                                           estimated_number_of_cells
                elif pipeline == "cellranger-arc":
                    # Single cell multiome output
                    summary_csv = os.path.join(
                        outs_dir,
                        "summary.csv")
                    if not os.path.exists(summary_csv):
                        raise Exception("Failed to add cell count "
                                        "for sample '%s': missing "
                                        "file '%s'" %
                                        (sample.name,
                                         summary_csv))
                    # Extract cell numbers
                    number_of_cells += MultiomeSummary(summary_csv).\
                                       estimated_number_of_cells
            # Store in the project metadata
            project.info['number_of_cells'] = number_of_cells
            project.info.save()
            return 0
        except Exception as ex:
            logger.warning("Unable to set cell count from data in "
                           "%s: %s" % (count_dir,ex))
    # Reached the end without setting count
    return 1

def cellranger_info(path=None,name=None):
    """
    Retrieve information on the cellranger software

    If called without any arguments this will locate the first
    cellranger executable that is available on the user's PATH,
    and attempts to extract the version.

    Alternatively if the path to an executable is supplied then
    the version will be determined from that instead.

    If no version is identified then the script path is still
    returned, but without any version info.

    If a 'path' is supplied then the package name will be taken
    from the basename; otherwise the package name can be supplied
    via the 'name' argument. If neither are supplied then the
    package name defaults to 'cellranger'.

    Returns:
      Tuple: tuple consisting of (PATH,PACKAGE,VERSION) where PATH
        is the full path for the cellranger program, PACKAGE is
        'cellranger', and VERSION is the package version. If any
        value can't be determined then it will be returned as an
        empty string.
    """
    # Initialise
    cellranger_path = ''
    if name is None:
        if path:
            name = os.path.basename(path)
        else:
            name = 'cellranger'
    package_name = name
    package_version = ''
    # Locate the core script
    if not path:
        cellranger_path = find_program(package_name)
    else:
        cellranger_path = os.path.abspath(path)
    # Identify the version
    if os.path.basename(cellranger_path) == package_name:
        # Run the program to get the version
        version_cmd = Command(cellranger_path,'--version')
        output = version_cmd.subprocess_check_output()[1]
        for line in output.split('\n'):
            if package_name in ('cellranger',
                                'cellranger-atac',):
                try:
                    if line.startswith("%s %s-" % (package_name,
                                                   package_name)):
                        # Extract version from line of the form
                        # cellranger cellranger-5.0.1
                        package_version = line.split('-')[-1]
                    elif line.startswith("%s " % package_name) and \
                         line.endswith(")"):
                        # Extract version from line of the form
                        # cellranger ... (2.0.1)
                        package_version = line.split('(')[-1].strip(')')
                    else:
                        # Raise an exception
                        raise("unrecognised version format")
                except Exception as ex:
                    logger.warning("Unable to get version from '%s': "
                                   "%s" % (line,ex))
            elif package_name == 'cellranger-arc':
                # Extract version from line of the form
                # cellranger-arc cellranger-arc-1.0.0
                try:
                    package_version = line.split('-')[-1]
                except Exception as ex:
                    logger.warning("Unable to get version from '%s': "
                                   "%s" % (line,ex))
    else:
        # No package supplied or located
        logger.warning("Unable to identify %s package from '%s'" %
                       (name,cellranger_path))
    # Return what we found
    return (cellranger_path,package_name,package_version)

def spaceranger_info(path=None,name=None):
    """
    Retrieve information on the spaceranger software

    If called without any arguments this will locate the first
    spaceranger executable that is available on the user's PATH,
    and attempts to extract the version.

    Alternatively if the path to an executable is supplied then
    the version will be determined from that instead.

    If no version is identified then the script path is still
    returned, but without any version info.

    If a 'path' is supplied then the package name will be taken
    from the basename; otherwise the package name can be supplied
    via the 'name' argument. If neither are supplied then the
    package name defaults to 'cellranger'.

    Returns:
      Tuple: tuple consisting of (PATH,PACKAGE,VERSION) where PATH
        is the full path for the spaceranger program, PACKAGE is
        'spaceranger', and VERSION is the package version. If any
        value can't be determined then it will be returned as an
        empty string.
    """
    # Initialise
    spaceranger_path = ''
    if name is None:
        if path:
            name = os.path.basename(path)
        else:
            name = 'spaceranger'
    package_name = name
    package_version = ''
    # Locate the core script
    if not path:
        spaceranger_path = find_program(package_name)
    else:
        spaceranger_path = os.path.abspath(path)
    # Identify the version
    if os.path.basename(spaceranger_path) == package_name:
        # Run the program to get the version
        version_cmd = Command(spaceranger_path,'--version')
        output = version_cmd.subprocess_check_output()[1]
        # Extract version from line of the from
        # spaceranger 1.1.0
        try:
            package_version = output.split()[-1]
        except Exception as ex:
            logger.warning("Unable to get version from '%s': %s" %
                           (output,ex))
    else:
        # No package supplied or located
        logger.warning("Unable to identify spaceranger package "
                       "from '%s'" % spaceranger_path)
    # Return what we found
    return (spaceranger_path,package_name,package_version)

def add_cellranger_args(cellranger_cmd,
                        jobmode=None,
                        maxjobs=None,
                        mempercore=None,
                        jobinterval=None,
                        localcores=None,
                        localmem=None,
                        disable_ui=False):
    """
    Configure options for cellranger

    Given a Command instance for running cellranger,
    add the appropriate options (e.g. --jobmode)
    according to the supplied arguments.

    Arguments:
      cellranger_cmd (Command): Command instance for
        running cellranger
      jobmode (str): if specified, will be passed to the
        --jobmode option
      maxjobs (int): if specified, will be passed to the
        --mempercore option
      mempercore (int): if specified, will be passed to
        the --maxjobs option (only if jobmode is not
        "local")
      jobinterval (int):  if specified, will be passed to
        the --jobinterval option
      localcores (int): if specified, will be passed to
        the --localcores option (only if jobmode is
        "local")
      localmem (int): if specified, will be passed to the
        the --localmem option (only if jobmode is
        "local")
      disable_ui (bool): if True, add the --disable-ui
        option (default is not to add it)

    Returns:
      Command: the original command updated with the
        appropriate options.
    """
    if jobmode is not None:
        cellranger_cmd.add_args("--jobmode=%s" % jobmode)
    if jobmode == "local":
        if localcores is not None:
            cellranger_cmd.add_args("--localcores=%s" %
                                    localcores)
        if localmem is not None:
            cellranger_cmd.add_args("--localmem=%s" %
                                    localmem)
    else:
        if mempercore is not None:
            cellranger_cmd.add_args("--mempercore=%s" %
                                    mempercore)
    if maxjobs is not None:
        cellranger_cmd.add_args("--maxjobs=%s" % maxjobs)
    if jobinterval is not None:
        cellranger_cmd.add_args("--jobinterval=%s" % jobinterval)
    if disable_ui:
        cellranger_cmd.add_args("--disable-ui")
    return cellranger_cmd
