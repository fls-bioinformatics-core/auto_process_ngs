#!/usr/bin/env python
#
#     tenx/metrics.py: classes for handling 10xGenomics metric summaries
#     Copyright (C) University of Manchester 2023 Peter Briggs
#

"""
Utility classes for handling the metric summary files produced by
various 10x Genomics pipelines:

- MetricSummary
- GexSummary
- AtacSummary
- MultiomeSummary
- MultiplexSummary
"""

#######################################################################
# Imports
#######################################################################

from bcftbx.TabFile import TabFile
from ..utils import parse_version

#######################################################################
# Classes
#######################################################################

class MetricsSummary(TabFile):
    """
    Base class for extracting data from cellranger* count
    *summary.csv files

    By default the files consists of two lines: the first
    is a header line, the second consists of corresponding
    data values. There is also a multi-line variation (e.g.
    from cellranger multi) with a header line followed by
    multiple lines of data (use the 'multiline' argument
    to indicate this is the expected format).

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
    def __init__(self,f,multiline=False):
        """
        Create a new MetricsSummary instance

        Arguments:
          f (str): path to the 'metrics_summary.csv' file
          multiline (bool): if True then expect multiple lines
            of data (default is to expect a single line of
            data)
        """
        # Store file
        self._filename = f
        # Multi-line summary?
        self._multiline = multiline
        # Read in data from the file
        with open(f,'rt') as fp:
            s = fp.read()
        self._data = dict()
        s = s.strip().split('\n')
        if not self._multiline and len(s) != 2:
            raise Exception("%s: MetricsSummary expects 2 lines (or specify "
                            "multi-line mode)" % f)
        # Set up the tabfile instance
        TabFile.__init__(self,
                         column_names=self._tokenise(s[0]),
                         delimiter=',')
        # Add the data
        for line in s[1:]:
            self.append(data=self._tokenise(line))
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
        if self._multiline:
            raise NotImplementedError("Superclass should implement "
                                      "'fetch' method for multi-line "
                                      "metrics files")
        return self[0][field]
    @property
    def metrics_file(self):
        return self._filename

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
        if parse_version(self.version) < parse_version("2.0.0"):
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
        if parse_version(self.version) < parse_version("2.0.0"):
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
        if parse_version(self.version) >= parse_version("2.0.0"):
            # Only supported for 2.0.0 and later
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
            # Not supported
            raise AttributeError
    @property
    def tss_enrichment_score(self):
        """
        Return the TSS enrichment score
        """
        try:
            return self.fetch('TSS enrichment score')
        except KeyError:
            return self.fetch('tss_enrichment_score')
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
    - gex_median_genes_per_cell
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
    def gex_median_genes_per_cell(self):
        """
        Return the median genes per cell for GEX data
        """
        return self.fetch('GEX Median genes per cell')

class MultiplexSummary(MetricsSummary):
    """
    Extract data from metrics_summary.csv file for CellPlex

    Utility class for extracting data from a
    'metrics_summary.csv' file output from running
    'cellranger multi'.

    The file consists of a header line followed by multiple
    data lines, with one set of values per line.

    The following properties are available:

    - cells
    - median_reads_per_cell
    - median_genes_per_cell
    - total_genes_detected
    - median_umi_counts_per_cell

    NB the returned values for these properties are all
    from the gene expression data.

    Values for other library types can be fetched
    directly by using the 'fetch' method and specifying
    the library type, for example:

    >>> MultiplexSummary(f).fetch('Cells','Antibody Capture')
    """
    def __init__(self,f):
        """
        Create a new MultiomeSummary instance

        Arguments:
          f (str): path to the 'summary.csv' file
        """
        MetricsSummary.__init__(self,f,multiline=True)
    def fetch(self,name,library_type='Gene Expression'):
        """
        Fetch data associated with an arbitrary field

        By default data associated with 'Gene Expression'
        are returned; data associated with other library
        types can be fetched by specifying a different
        value for the 'library_type' argument (either
        'Multiplexing Capture' or 'Antibody Capture').

        Arguments:
          name (str): name of the metric
          library_type (str): library type to fetch
            the metric for (default: 'Gene Expression')
        """
        metric = self.lookup('Metric Name',name)
        if not metric:
            raise Exception("Failed to lookup metric '%s'" % name)
        # Pick out the specific library
        for m in metric:
            if m['Library Type'] == library_type:
                return m['Metric Value']
        # No matching library
        raise Exception("No value for metric '%s' associated with "
                        "library type '%s'" % (name,library_type))
    @property
    def cells(self):
        """
        Returns the number of cells
        """
        return self.fetch('Cells')
    @property
    def median_reads_per_cell(self):
        """
        Returns the median reads per cell
        """
        return self.fetch('Median reads per cell')
    @property
    def median_genes_per_cell(self):
        """
        Returns the median genes per cell
        """
        return self.fetch('Median genes per cell')
    @property
    def total_genes_detected(self):
        """
        Returns the total genes detected
        """
        return self.fetch('Total genes detected')
    @property
    def median_umi_counts_per_cell(self):
        """
        Returns the median UMI counts per cell
        """
        return self.fetch('Median UMI counts per cell')
