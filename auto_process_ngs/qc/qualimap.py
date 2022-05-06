#!/usr/bin/env python
#
# qualimap library
import os

"""
Example 'rnaseq_qc_results.txt' file:

RNA-Seq QC report
-----------------------------------

>>>>>>> Input

    bam file = /path/to/SMP1_S1_R1_001.bam
    gff file = /path/to/annotation.gtf
    counting algorithm = uniquely-mapped-reads
    protocol = strand-specific-reverse


>>>>>>> Reads alignment

    reads aligned (left/right) = 175,847 / 175,847
    read pairs aligned  = 175,847
    total alignments = 728,678
    secondary alignments = 376,984
    non-unique alignments = 357,582
    aligned to genes  = 63,784
    ambiguous alignments = 160,764
    no feature assigned = 19,130
    not aligned = 0


>>>>>>> Reads genomic origin

    exonic =  63,784 (76.93%)
    intronic = 4,778 (5.76%)
    intergenic = 14,352 (17.31%)
    overlapping exon = 6,523 (7.87%)
    rRNA = 0 (0%)


>>>>>>> Transcript coverage profile

    5' bias = 0.73
    3' bias = 0.41
    5'-3' bias = 1.52


>>>>>>> Junction analysis

    reads at junctions = 28,238

    AGAT : 3.22%
    AGGT : 3.18%
    ATCT : 2.79%
    ACCT : 2.79%
    ...
"""

class QualimapRnaseq(object):
    """
    Wrapper class for handling outputs from Qualimap 'rnaseq'

    The ``QualimapRnaseq`` object gives access to various
    aspects of the outputs of the Qualimap ``rnaseq`` module.

    The following properties are available:

    - html_report (str): path to the HTML report
    - results_text (str): path to the results .txt file
    - input (dict): access values in the 'Input' section
    - reads_genomic_origin (dict): access values in the
      'Reads Genomic Origin' section

    """
    def __init__(self,qualimap_dir):
        """
        """
        # Initialise
        self._qualimap_dir = os.path.abspath(qualimap_dir)
        self._raw_results = dict()
        self._input = None
        self._reads_genomic_origin = None
        self._raw_coverage_profile_along_genes_total = None
        # Read in data from 'rnaseq_qc_results.txt'
        self._results_txt = os.path.join(self._qualimap_dir,
                                         "rnaseq_qc_results.txt")
        with open(self._results_txt,'rt') as fp:
            section = None
            for line in fp:
                line = line.strip()
                if line.startswith(">>>>>>> "):
                    # New section
                    section = ' '.join(line.split(' ')[1:])
                    if section not in self._raw_results:
                        self._raw_results[section] = list()
                elif section and line:
                    self._raw_results[section].append(line)
                else:
                    pass

    @property
    def html_report(self):
        """
        Return path to the 'qualimapReport.html' HTML report
        """
        return os.path.join(self._qualimap_dir,
                            "qualimapReport.html")

    @property
    def results_txt(self):
        """
        Return path to the 'rnaseq_qc_results.txt' file
        """
        return self._results_txt

    @property
    def input(self):
        """
        Provides access to values in the 'Input' section
        """
        if not self._input:
            data = dict()
            for line in self._raw_results['Input']:
                # Line of the form
                # 'bam file = /path/to/SMP1_S1_R1_001.bam'
                item,value = line.split('=')
                item = item.strip()
                data[item] = value.strip()
            self._input = data
        return self._input

    @property
    def reads_genomic_origin(self):
        """
        Provides access to values in the 'Reads genomic origin' section
        """
        if self._reads_genomic_origin is None:
            data = dict()
            for line in self._raw_results['Reads genomic origin']:
                # Line of the form 'exonic =  63,784 (76.93%)'
                item,values = line.split('=')
                item = item.strip()
                count,perc = values[:-2].split('(')
                count = int(count.strip().replace(',',''))
                perc = float(perc)
                data[item] = (count,perc)
            self._reads_genomic_origin = data
        return self._reads_genomic_origin

    @property
    def raw_coverage_profile_along_genes_total(self):
        """
        Provides access to values in the 'raw_data_qualimapReport/coverage_profile_along_genes_\(total\).txt' file
        """
        if self._raw_coverage_profile_along_genes_total is None:
            f = os.path.join(self._qualimap_dir,
                             "raw_data_qualimapReport",
                             "coverage_profile_along_genes_(total).txt")
            # FIXME use TabFile instead?
            data = dict()
            with open(f,'rt') as fp:
                for line in fp:
                    if line.startswith('#'):
                        # Ignore header line
                        continue
                    # Line of the form '
                    x,y = line.rstrip().split('\t')
                    data[int(float(x))] = float(y)
            self._raw_coverage_profile_along_genes_total = data
        return self._raw_coverage_profile_along_genes_total
