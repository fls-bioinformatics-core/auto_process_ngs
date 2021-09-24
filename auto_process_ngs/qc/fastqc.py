#!/usr/bin/env python
#
# fastqc library
import os
from bcftbx.TabFile import TabFile
from bcftbx.htmlpagewriter import PNGBase64Encoder
from ..docwriter import Table
from ..docwriter import Link

"""
Example Fastqc summary text file (FASTQ_fastqc/summary.txt):

PASS	Basic Statistics	ES1_GTCCGC_L008_R1_001.fastq.gz
PASS	Per base sequence quality	ES1_GTCCGC_L008_R1_001.fastq.gz
PASS	Per tile sequence quality	ES1_GTCCGC_L008_R1_001.fastq.gz
PASS	Per sequence quality scores	ES1_GTCCGC_L008_R1_001.fastq.gz
FAIL	Per base sequence content	ES1_GTCCGC_L008_R1_001.fastq.gz
WARN	Per sequence GC content	ES1_GTCCGC_L008_R1_001.fastq.gz
PASS	Per base N content	ES1_GTCCGC_L008_R1_001.fastq.gz
PASS	Sequence Length Distribution	ES1_GTCCGC_L008_R1_001.fastq.gz
FAIL	Sequence Duplication Levels	ES1_GTCCGC_L008_R1_001.fastq.gz
PASS	Overrepresented sequences	ES1_GTCCGC_L008_R1_001.fastq.gz
PASS	Adapter Content	ES1_GTCCGC_L008_R1_001.fastq.gz
FAIL	Kmer Content	ES1_GTCCGC_L008_R1_001.fastq.gz

Head of the FastQC data file (FASTQ_fastqc/fastqc_data.txt), which
contains raw numbers for the plots etc):

##FastQC	0.11.3
>>Basic Statistics	pass
#Measure	Value
Filename	ES1_GTCCGC_L008_R1_001.fastq.gz
File type	Conventional base calls
Encoding	Sanger / Illumina 1.9
Total Sequences	12317096
Sequences flagged as poor quality	0
Sequence length	101
%GC	50
>>END_MODULE
>>Per base sequence quality	pass
#Base	Mean	Median	Lower Quartile	Upper Quartile	10th Percentile	90th Pe
rcentile
1	32.80553403172306	33.0	33.0	33.0	33.0	33.0
...

"""

class Fastqc(object):
    """
    Wrapper class for handling outputs from FastQC

    The ``Fastqc`` object gives access to various
    aspects of the outputs of the FastQC program.

    """
    # Base names for plots in the 'Images' subdir
    plot_names = ('adapter_content',
                  'duplication_levels',
                  'kmer_profiles',
                  'per_base_n_content',
                  'per_base_quality',
                  'per_base_sequence_content',
                  'per_sequence_gc_content',
                  'per_sequence_quality',
                  'per_tile_quality',
                  'sequence_length_distribution',
                  )
    def __init__(self,fastqc_dir):
        """
        Create a new Fastqc instance

        Arguments:
          fastqc_dir (str): path to the top-level
            output directory from a FastQC run.
        """
        self._fastqc_dir = os.path.abspath(fastqc_dir)
        self._fastqc_summary = FastqcSummary(
            summary_file=os.path.join(self._fastqc_dir,
                                      'summary.txt'))
        self._fastqc_data = FastqcData(
            os.path.join(self._fastqc_dir,
                         'fastqc_data.txt'))
        self._html_report = self._fastqc_dir + '.html'
        self._zip = self._fastqc_dir + '.zip'

    @property
    def version(self):
        """
        Version of FastQC that was used
        """
        return self.data.version

    @property
    def dir(self):
        """
        Path to the directory with the FastQC outputs
        """
        return self._fastqc_dir

    @property
    def html_report(self):
        """
        Path to the associated HTML report file
        """
        return self._html_report

    @property
    def zip(self):
        """
        Path to the associated ZIP archive
        """
        return self._zip

    @property
    def summary(self):
        """
        Return a FastqcSummary instance
        """
        return self._fastqc_summary

    @property
    def data(self):
        """
        Return a FastqcData instance
        """
        return self._fastqc_data

    def plot(self,module,inline=False):
        """
        """
        # Normalise name
        name = module.lower().replace(' ','_')
        plot_png = os.path.join(self._fastqc_dir,
                                'Images',
                                '%s.png' % name)
        # Check png exists
        if not os.path.exists(plot_png):
            return None
        # Return requested format
        if inline:
            return "data:image/png;base64," + \
                PNGBase64Encoder().encodePNG(plot_png)
        else:
            return plot_png

    def quality_boxplot(self,inline=False):
        """
        """
        return self.plot('per_base_quality',
                         inline=inline)

    def adapter_content_plot(self,inline=False):
        """
        """
        return self.plot('Adapter Content',
                         inline=inline)

class FastqcSummary(TabFile):
    """
    Class representing data from a Fastqc summary file
    """
    # Names for Fastqc modules
    module_names = (
        'Basic Statistics',
        'Per base sequence quality',
        'Per tile sequence quality',
        'Per sequence quality scores',
        'Per base sequence content',
        'Per sequence GC content',
        'Per base N content',
        'Sequence Length Distribution',
        'Sequence Duplication Levels',
        'Overrepresented sequences',
        'Adapter Content',
        'Kmer Content',
    )
    def __init__(self,summary_file=None):
        """
        Create a new FastqcSummary instance

        Arguments:
          summary_file (str): path to the FastQC
            ``summary.txt`` file
        """
        TabFile.__init__(self,
                         column_names=('Status',
                                       'Module',
                                       'File',))
        if summary_file:
            summary_file = os.path.abspath(summary_file)
            with open(summary_file,'r') as fp:
                for line in fp:
                    line = line.strip()
                    self.append(tabdata=line)
        self._summary_file = summary_file

    @property
    def path(self):
        # Path to the summary file
        return self._summary_file

    @property
    def modules(self):
        # Return list of modules
        return [r['Module'] for r in self]

    def status(self,name):
        # Return status for module 'name'
        return self.lookup('Module',name)[0]['Status']

    @property
    def passes(self):
        """
        Return modules with passes

        Returns a list with the names of the modules that
        have status 'PASS'.
        """
        return [r['Module'] for r in filter(lambda x: x['Status'] == 'PASS',
                                            self)]

    @property
    def warnings(self):
        """
        Return modules with warnings

        Returns a list with the names of the modules that
        have status 'WARN'.
        """
        return [r['Module'] for r in filter(lambda x: x['Status'] == 'WARN',
                                            self)]

    @property
    def failures(self):
        """
        Return modules with failures

        Returns a list with the names of the modules that
        have status 'FAIL'.
        """
        return [r['Module'] for r in filter(lambda x: x['Status'] == 'FAIL',
                                            self)]

    def link_to_module(self,name,full_path=True,relpath=None):
        """
        Return link to the result of a specified FastQC module

        Arguments:
          name (str): name of the module (e.g. 'Basic Statistics')
          full_path (boolean): optional, if True then return the
            full path; otherwise return just the anchor (e.g. '#M1')
          relpath (str): optional, if supplied then specifies the
            path that full paths will be made relative to (implies
            full_path is True)

        """
        i = self.module_names.index(name)
        link = "#M%d" % i
        if full_path:
            link = self.html_report() + link
            if relpath:
                link = os.path.relpath(link,relpath)
        return link

    def html_report(self):
        """
        Return the path of the HTML report from FastQC
        """
        return os.path.dirname(self.path)+'.html'

    def html_table(self,relpath=None):
        """
        Generate HTML table for FastQC summary

        Arguments:
          relpath (str): optional, if supplied then links in the
            table will be relative to this path

        """
        tbl = Table(('module','status'),
                    module='FastQC test',status='Outcome')
        tbl.add_css_classes('fastqc_summary','summary')
        for name in self.modules:
            tbl.add_row(module=Link(name,self.link_to_module(name,
                                                             relpath=relpath)),
                        status="<span class='%s'>%s</span>" % (
                            self.status(name),
                            self.status(name)))
        return tbl.html()

class FastqcData(object):
    """
    Class representing data from a Fastqc data file

    Reads in the data from a ``fastqc_data.txt`` file
    and makes it available programmatically.

    To create a new FastqcData instance:

    >>> fqc = FastqcData('fastqc_data.txt')

    To access a field in the 'Basic Statistics' module:

    >>> nreads = fqc.basic_statistics('Total Sequences')
    """
    def __init__(self,data_file):
        """
        Create a new FastqcData instance

        Arguments:
          data_file (str): path to a ``fastqc_data.txt``
            file which will be read in and processed
        """
        self._data_file = os.path.abspath(data_file)
        self._fastqc_version = None
        self._modules = {}
        if data_file:
            fastqc_module = None
            with open(data_file,'r') as fp:
                for line in fp:
                    line = line.strip()
                    if fastqc_module is None:
                        if line.startswith('##FastQC'):
                            self._fastqc_version = line.split()[-1]
                        elif line.startswith('>>'):
                            fastqc_module = line.split('\t')[0][2:]
                            self._modules[fastqc_module] = []
                    elif line.startswith('>>END_MODULE'):
                        fastqc_module = None
                    else:
                        self._modules[fastqc_module].append(line)

    @property
    def version(self):
        """
        FastQC version number
        """
        return self._fastqc_version

    @property
    def path(self):
        """
        Path to the fastqc_data.txt file
        """
        return self._data_file

    @property
    def modules(self):
        """
        List of the modules in the raw data
        """
        return self._modules

    def data(self,module):
        """
        Return the raw data for a module

        Returns the data for the specified module as a list
        of lines.

        The first list item/line is the header line; data items
        within each line are tab-delimited.

        For example:

        >>> Fastqc('myfastq_fastq').data.data('Sequence Length Distribution')
        ['#Length\tCount',
         '35\t8826.0',
         '36\t2848.0',
         '37\t4666.0',
         '38\t4524.0']
        """
        if module in self._modules:
            return self._modules[module]
        return None

    def basic_statistics(self,measure):
        """
        Access a data item in the ``Basic Statistics`` section

        Possible values include:

        - Filename
        - File type
        - Encoding
        - Total Sequences
        - Sequences flagged as poor quality
        - Sequence length
        - %GC

        Arguments:
          measure (str): key corresponding to a 'measure'
            in the ``Basic Statistics`` section.

        Returns:
          String: value of the requested 'measure'

        Raises:
          KeyError: if 'measure' is not found.
        """
        for line in self.data('Basic Statistics'):
            key,value = line.split('\t')
            if key == measure:
                return value
        raise KeyError("No measure '%s'" % measure)
