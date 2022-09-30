#!/usr/bin/env python
#
# picard library
import os
from bcftbx.utils import AttributeDictionary

"""
Example CollectInsertSizeMetrics output (SAMPLE.insert_size_metrics.txt):

## htsjdk.samtools.metrics.StringHeader
# CollectInsertSizeMetrics --Histogram_FILE <SAMPLE>.insert_size_histogram.pdf --INPUT <BAM> --OUTPUT <BASENAME>.insert_size_metrics.txt --DEVIATIONS 10.0 --MINIMUM_PCT 0.05 --METRIC_ACCUMULATION_LEVEL ALL_READS --INCLUDE_DUPLICATES false --ASSUME_SORTED true --STOP_AFTER 0 --VERBOSITY INFO --QUIET false --VALIDATION_STRINGENCY STRICT --COMPRESSION_LEVEL 5 --MAX_RECORDS_IN_RAM 500000 --CREATE_INDEX false --CREATE_MD5_FILE false --GA4GH_CLIENT_SECRETS client_secrets.json --help false --version false --showHidden false --USE_JDK_DEFLATER false --USE_JDK_INFLATER false
## htsjdk.samtools.metrics.StringHeader
# Started on: Mon Oct 04 12:53:24 BST 2021

## METRICS CLASS	picard.analysis.InsertSizeMetrics
MEDIAN_INSERT_SIZE	MODE_INSERT_SIZE	MEDIAN_ABSOLUTE_DEVIATION	MIN_INSERT_SIZE	MAX_INSERT_SIZE	MEAN_INSERT_SIZE	STANDARD_DEVIATION	READ_PAIRS	PAIR_ORIENTATION	WIDTH_OF_10_PERCENT	WIDTH_OF_20_PERCENT	WIDTH_OF_30_PERCENT	WIDTH_OF_40_PERCENT	WIDTH_OF_50_PERCENT	WIDTH_OF_60_PERCENT	WIDTH_OF_70_PERCENT	WIDTH_OF_80_PERCENT	WIDTH_OF_90_PERCENT	WIDTH_OF_95_PERCENT	WIDTH_OF_99_PERCENT	SAMPLE	LIBRARY	READ_GROUP
139	103	37	28	753074	153.754829	69.675347	175847	FR	15	31	47	61	75	93	115	157	323	975	31983

## HISTOGRAM	java.lang.Integer
insert_size	All_Reads.fr_count
28	1
29	1
30	1
31	1
32	1
33	2
34	3
...

"""

class CollectInsertSizeMetrics(object):
    """
    Wrapper class for outputs from CollectInsertSizeMetrics

    The ``CollectInsertSizeMetrics`` object gives access to
    various aspects of the outputs of the Picard
    ``CollectInsertSizeMetrics`` utility.

    The following properties are available:

    - metrics (dict): dictionary mapping metrics to values
    - histogram (dict): dictionary holding histogram data

    """
    def __init__(self,metrics_file):
        """
        Create a new CollectInsertSizeMetrics instance

        Arguments:
          metrics_file (str): path to the .insert_size_metrics.txt
            output file from collectInsertSizeMetrics
        """
        # Initialise
        self._metrics_file = os.path.abspath(metrics_file)
        self._metrics = dict()
        self._histogram = dict()
        # Read in raw data for each section from metrics file
        raw_data = dict()
        with open(self._metrics_file,'rt') as fp:
            section = None
            for line in fp:
                if line.startswith("## METRICS CLASS"):
                    section = 'metrics'
                elif line.startswith("## HISTOGRAM"):
                    section = 'histogram'
                elif line.startswith('#') or not line.strip():
                    section = None
                elif section:
                    if section not in raw_data:
                        raw_data[section] = list()
                    raw_data[section].append(line.rstrip('\n'))
                else:
                    pass
        # Convert metrics to AttributeDict
        self._metrics = AttributeDictionary()
        for key,value in zip(raw_data['metrics'][0].split('\t'),
                             raw_data['metrics'][1].split('\t')):
            try:
                value = int(value)
            except ValueError:
                try:
                    value = float(value)
                except ValueError:
                    pass
            self._metrics[key] = value
        # Convert histogram to dictionary
        self._histogram = dict()
        for line in raw_data['histogram'][1:]:
            position,count = line.split('\t')
            self._histogram[int(position)] = int(count)

    @property
    def metrics_file(self):
        """
        Return path to source metrics file
        """
        return self._metrics_file

    @property
    def metrics(self):
        """
        Dictionary mapping metrics to values

        To get the value associated with a metric,
        use e.g.:

        >>> mean = insertsizes.metrics['MEAN_INSERT_SIZE']

        To see what metrics are available, use e.g.:

        >>> insertsizes.metrics.keys()
        """
        return self._metrics

    @property
    def histogram(self):
        """
        Histogram data for insert sizes

        Dictionary mapping insert sizes (keys) to
        associated number of alignments.
        """
        return self._histogram
