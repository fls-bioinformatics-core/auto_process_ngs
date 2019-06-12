#!/usr/bin/env python
#
#     barcodes/analysis.py: classes and functions for analysing barcodes
#     Copyright (C) University of Manchester 2016-2019 Peter Briggs
#
########################################################################
#
# barcode_analysis.py
#
#########################################################################

"""
Classes and functions for analysing barcodes (i.e. index sequences) from
FASTQ read headers:

- BarcodeCounter: utility class for counting barcode sequences
- BarcodeGroup: class for sorting groups of related barcodes
- SampleSheetBarcode: class for handling barcode information from a
  sample sheet
- Reporter: class for generating reports of barcode statistics
- report_barcodes: populate Reporter with analysis of barcode stats
- detect_barcodes_warnings: check whether reports include warnings
- make_title: turn a string into a Markdown/rst title
- percent: return values as percentage
"""

#######################################################################
# Imports
#######################################################################

import sys
from itertools import izip
from bcftbx.IlluminaData import SampleSheet
from bcftbx.IlluminaData import samplesheet_index_sequence
from bcftbx.IlluminaData import normalise_barcode
from bcftbx.utils import AttributeDictionary
from bcftbx.simple_xls import XLSWorkBook
from bcftbx.simple_xls import XLSStyle
from ..docwriter import Document
from ..docwriter import Table
from ..docwriter import List
from ..docwriter import Link
from ..docwriter import Img
from ..docwriter import Para
from ..docwriter import WarningIcon

#######################################################################
# Data
#######################################################################

PROBLEMS_DETECTED_TEXT = "One or more lanes have problems"

#######################################################################
# Classes
#######################################################################

class BarcodeCounter(object):
    """
    Utility class to mange counts of barcode sequences

    The BarcodeCounter class manages the counting of
    barcode index sequences and facilitates subsequent
    analysis of their frequency.

    Usage
    -----

    To set up a new (empty) BarcodeCounter:

    >>> bc = BarcodeCounter()

    To count the sequences from a fastq file (keeping the
    lane association):

    >>> for read in FastqIterator("example.fq"):
    ...   bc.count_barcode(read.seqid.index_sequence,
    ...                    lane=read.seqid.flowcell_lane)

    To get a list of the barcodes sorted by order of
    frequency (highest to lowest) across all lanes:

    >>> bc.barcodes()

    To get the count for a barcode in lane 1:

    >>> bc.counts('TAGGCATGTAGCCTCT',1)

    Reading and writing counts
    --------------------------

    Counts can be output to a file:

    >>> bc.write("counts.out")

    and subsequently loaded into a new BarcodeCounter:

    >>> bc2 = BarcodeCounter("counts.out")

    Multiple counts files can be combined:

    >>> bc_all = BarcodeCounter("counts1.out","counts2.out")

    Grouping barcodes
    -----------------

    To put barcodes into groups use the 'group' method,
    e.g. to group all barcodes from lane 2:

    >>> groups = bc.group(2)

    which produces a list of BarcodeGroup instances.

    """
    def __init__(self,*counts_files):
        """
        Create a new BarcodeCounter instance

        Arguments:
          counts_files: (optionally) one or more count
           files; if there are multiple files then the
           data will be combined for all input files

        """
        self._seqs = {}
        self._seqs_all = {}
        for counts_file in counts_files:
            self.read(counts_file)

    def count_barcode(self,barcode,lane=None,incr=1):
        """
        Increment count of a barcode sequence

        Arguments:
          barcode (str): barcode sequence to count
          lane (int): lane that the barcode appears
            in (None if unknown)
          incr (int): increment the count for the
            barcode in the lane by this amount
            (defaults to 1)

        """
        # Normalise barcode
        barcode = normalise_barcode(barcode)
        # Store by lane
        try:
            self._seqs[lane][barcode] += incr
        except KeyError:
            try:
                self._seqs[lane][barcode] = incr
            except KeyError:
                self._seqs[lane] = { barcode: incr }
        # Store overall
        try:
            self._seqs_all[barcode] += incr
        except KeyError:
            self._seqs_all[barcode] = incr

    def barcode_lengths(self,lane=None):
        """
        Return lengths of barcode sequences

        Returns a list of the barcode sequence lengths.

        Arguments:
          lane (int): if specified then restricts the
            list to barcodes that appear in the named
            lane (default is to get lengths from all
            barcodes in all lanes)
        """
        lengths = set()
        for barcode in self.barcodes(lane=lane):
            lengths.add(len(normalise_barcode(barcode)))
        return sorted(list(lengths))

    @property
    def lanes(self):
        """
        List of lane numbers that barcodes are counted in

        Returns:
          List: list of integer lane numbers in
            ascending order.

        """
        lanes = self._seqs.keys()
        if lanes == [None]:
            return []
        else:
            return sorted(lanes)

    def barcodes(self,lane=None):
        """
        List barcodes sorted into order of frequency
        
        Optionally restrict the list to frequency in
        a specified lane

        Arguments:
          lane (int): if not None then only list
            barcodes that appear in the specified lane
            (and ordering only uses the frequencies in
            that lane)

        Returns:
          List: list of barcodes.

        """
        if lane is None:
            return sorted([s for s in self._seqs_all],
                          cmp=lambda x,y: cmp(self.counts_all(y),
                                              self.counts_all(x)))
        else:
            return sorted([s for s in self._seqs[lane]],
                          cmp=lambda x,y: cmp(self.counts(y,lane),
                                              self.counts(x,lane)))

    def filter_barcodes(self,cutoff=None,lane=None):
        """
        Return subset of index sequences filtered by specified criteria

        Arguments:
          cutoff (float): barcodes must account for at least this
            fraction of all reads to be included. Total reads are
            all reads if lane is 'None', or else total for the
            specified lane
          lane (integer): barcodes must appear in this lane

        Returns:
          List: list of barcodes.

        """
        # Initialise
        bc = self.barcodes(lane=lane)
        nreads = self.nreads(lane=lane)
        # Apply cutoff if specified
        # i.e. exclude reads that are less than specified
        # fraction of total reads
        if cutoff is not None:
            cutoff_reads = int(float(nreads)*cutoff)
            bc = [s for s in bc if self.counts(s,lane) >= cutoff_reads]
        return bc

    def counts(self,barcode,lane=None):
        """
        Return the number of counts for the barcode

        If 'lane' is None then return counts across all
        lanes.

        """
        try:
            if lane is None:
                return self._seqs_all[barcode]
            else:
                return self._seqs[lane][barcode]
        except KeyError:
            return 0

    def counts_all(self,barcode):
        """
        Number of counts for the barcode across all lanes

        """
        return self.counts(barcode)

    def nreads(self,lane=None):
        """
        Number of reads counted

        If lane is None then is the number of reads
        across all lanes (i.e. all reads from all files).

        Arguments:
          lane (int): only report number of reads for
            the specified lane

        Returns:
          Integer: number of reads.

        """
        if lane is None:
            return sum([self._seqs_all[barcode]
                        for barcode in self._seqs_all])
        else:
            try:
                return sum([self._seqs[lane][barcode]
                            for barcode in self._seqs[lane]])
            except KeyError:
                return 0

    def read(self,filen):
        """
        Read count data from a file

        The format of the 'counts' file is four column
        tab-delimited file:

        - Column 1: lane
        - Column 2: rank (ignored when reading in)
        - Column 3: barcode sequence
        - Column 4: counts

        Older 'counts' files are three column tab-delimited
        with the following columns:

        - Column 1: rank (ignored when reading in)
        - Column 2: barcode sequence
        - Column 3: counts

        In either case if the lane is missing or cannot be
        interpreted as an integer then it's set to be 'None'.

        """
        print "Reading count data from file '%s'" % filen
        counts = {}
        with open(filen,'r') as fp:
            for line in fp:
                # Skip comment line
                if line.startswith('#'):
                    continue
                # Split into elements
                items = line.split('\t')
                if len(items) == 3:
                    # Old style counts file (no lane column)
                    lane = None
                    barcode = items[1]
                    counts = int(items[2])
                elif len(items) == 4:
                    # New style (includes lane)
                    try:
                        lane = int(items[0])
                    except ValueError:
                        lane = None
                    barcode = items[2]
                    counts = int(items[3])
                else:
                    raise ValueError("Wrong number of items in line")
                # Store the data
                self.count_barcode(barcode,lane,counts)

    def write(self,filen):
        """
        Write barcode data to a file

        """
        print "Writing all counts to file '%s'" % filen
        with open(filen,'w') as fp:
            fp.write("#Lane\tRank\tSequence\tCount\n")
            for lane in self.lanes:
                for i,seq in enumerate(self.barcodes(lane)):
                    fp.write("%d\t%d\t%s\t%d\n" % (lane,
                                                   i+1,
                                                   seq,
                                                   self.counts(seq,lane)))

    def group(self,lane,mismatches=2,n=None,cutoff=None,
              seed_barcodes=None,exclude_reads=0.000001):
        """
        Put barcodes into groups of similar sequences

        Returns a list of BarcodeGroup instances

        Arguments:
          lane: lane number to restrict the pool of barcodes to
            (set to None to use barcodes from all lanes)
          mismatches: number of mismatches to allow when creating
            groups
          cutoff: minimum number of reads as a fraction of all
            reads that a group must contain to be included
            (set to None to disable cut-off)
          exclude_reads: speed-up parameter, excludes barcodes with
            less than this fraction of associated reads. Speeds up
            the grouping calculation at the cost of some precision
          seed_barcodes (list): optional, set of barcode sequences
            (typically, expected index sequences from a sample sheet)
            which will be used to build groups around even if they
            have low associated counts

        """
        # Initialise
        barcodes = self.filter_barcodes(lane=lane,cutoff=exclude_reads)
        nreads = self.nreads(lane=lane)
        groups = []
        # Update barcode list if 'seed' barcodes were provided
        if seed_barcodes:
            promoted_barcodes = []
            for barcode in seed_barcodes:
                try:
                    barcodes.remove(barcode)
                    promoted_barcodes.append(barcode)
                except ValueError:
                    # Barcode doesn't appear in the list
                    pass
            barcodes = promoted_barcodes + barcodes
        # Cutoff
        if cutoff is not None:
            cutoff_reads = int(float(nreads)*cutoff)
        else:
            cutoff_reads = 0
        # Iteratively assign barcodes to groups
        # until we run out
        cumulative = 0
        while barcodes:
            # Fetch next reference sequence
            group = BarcodeGroup(barcodes[0],
                                 self.counts(barcodes[0],lane))
            # Save non-matching sequences
            rejected = []
            # Iterate through the remaining sequences
            # looking for matches
            for seq in barcodes[1:]:
                if group.match(seq,mismatches):
                    group.add(seq,self.counts(seq,lane))
                else:
                    rejected.append(seq)
            # Finished checking sequences for this group
            # Reset sequences to check in next iteration
            barcodes = rejected
            # Check cutoff
            if group.counts >= cutoff_reads:
                groups.append(group)
            else:
                # Discard group
                pass
        # Sort groups into order of total counts
        groups = sorted(groups,cmp=lambda x,y: cmp(y.counts,x.counts))
        return groups

    def analyse(self,lane=None,sample_sheet=None,cutoff=None,
                mismatches=0):
        """
        Analyse barcode frequencies

        Returns a dictionary with the following keys:

        - barcodes: list of barcodes (or reference barcodes,
          if mismatches > 0)
        - cutoff: the specified cutoff fraction
        - mismatches: the specified number of mismatches to
          allow
        - total_reads: the total number of reads for the
          specified lane (or all reads, if no lane was
          specified)
        - coverage: the number of reads after cutoffs have
          been applied
        - counts: dictionary with barcodes from the 'barcodes'
          list as keys; each key points to a dictionary with
          keys:
          * reads: number of reads associated with this barcode
            (or group, if mismatches > 0)
          * sample: name of the associated sample (if a sample
            sheet was supplied, otherwise 'None')
          * sequences: number of sequences in the group (always
            1 if mismatches == 0)

        Arguments:
          lane (integer): lane to restrict analysis to (None
            analyses all lanes)
          sample_sheet (str): sample sheet file to compare
            barcodes against (None skips comparison)
          cutoff (float): if mismatches == 0 then barcodes must
            have at least this fraction of reads to be included;
            (if mismatches > 0 then this condition is applied to
            groups instead)

        """
        sample_lookup = {}
        if sample_sheet is not None:
            sample_sheet = SampleSheetBarcodes(sample_sheet)
            sample_sheet_barcodes = sample_sheet.barcodes(lane)
        else:
            sample_sheet_barcodes = None
        if not mismatches:
            groups = None
            barcodes = self.filter_barcodes(cutoff=cutoff,lane=lane)
        else:
            groups = self.group(lane,mismatches=mismatches,
                                seed_barcodes=sample_sheet_barcodes,
                                cutoff=cutoff)
            barcodes = [grp.reference for grp in groups]
        analysis = AttributeDictionary(
            barcodes=barcodes,
            cutoff=cutoff,
            counts=dict(),
            total_reads=self.nreads(lane=lane),
            mismatches=mismatches
        )
        cum_reads = 0
        if groups:
            for group in groups:
                barcode = group.reference
                barcode_reads = group.counts
                cum_reads += barcode_reads
                try:
                    # Exact match
                    sample = sample_sheet.lookup_sample(barcode,lane)
                except KeyError:
                    # Closest match(es)
                    sample = []
                    for seq in sample_sheet.barcodes(lane):
                        if group.match(seq,mismatches):
                            sample.append(sample_sheet.lookup_sample(seq,lane))
                    if sample:
                        sample = ','.join(sample)
                    else:
                        sample = None
                except AttributeError:
                    # No sample sheet
                    sample = None
                analysis.counts[barcode] = AttributeDictionary(
                    reads=barcode_reads,
                    sample=sample,
                    sequences=len(group)
                )
        else:
            for barcode in barcodes:
                barcode_reads = self.counts(barcode,lane)
                cum_reads += barcode_reads
                try:
                    sample = sample_sheet.lookup_sample(barcode,lane)
                except (KeyError,AttributeError):
                    sample = None
                analysis.counts[barcode] = AttributeDictionary(
                    reads=barcode_reads,
                    sample=sample,
                    sequences=1
                )
        analysis['coverage'] = cum_reads
        return analysis

class BarcodeGroup(object):
    """
    Class for storing groups of related barcodes

    A group stores a representative 'reference' sequence
    and none or more related sequences, along with the
    total counts for all the sequences combined.

    Create a group:

    >>> grp = BarcodeGroup('GCTACGCTCTAAGCCT',2894178)

    Add a barcode to the group if it's related:

    >>> if grp.match('GCTCCGCTCTAAGCCT'):
    ...   grp.add('GCTCCGCTCTAAGCCT',94178)

    List the barcodes in the group:

    >>> grp.sequences

    Get the total number of counts:

    >>> grp.counts

    Retrieve the reference barcode:

    >>> grp.reference
    
    """

    def __init__(self,barcode,counts=0):
        """
        Create a new BarcodeGroup instance

        Arguments:
          barcode (str): representative reference
            barcode for the group
          counts (int): number of counts for the
            reference (defaults to zero)

        """
        self._barcode = barcode
        self._sequences = [barcode]
        self._counts = counts

    def add(self,seq,counts):
        """
        Add a sequence to the group

        Arguments:
          seq (str): sequence to add to the group
          counts (int): count for that sequence (will be
            added to the total for the group)

        """
        self._sequences.append(seq)
        self._counts += counts

    def match(self,seq,mismatches=2):
        """
        Check if a sequence is related to the reference

        Note that if sequences differ in length then they
        automatically fail to match.

        Arguments:
          seq (str): sequence to check against the reference
          mismatches (int): maximum number of mismatches that
            are allowed for the sequences to be considered as
            related (default is 2). Note that 'N's in either
            sequence automatically count as a mismatch.

        Returns:
          Boolean: True if sequences match within the
            specified tolerance, False otherwise (or if
            sequence lengths differ)

        """
        if len(self._barcode) != len(seq):
            return False
        m = 0
        for c1,c2 in izip(self._barcode,seq):
            if c1 == 'N' or c2 == 'N' or c1 != c2:
                m += 1
                if m > mismatches:
                    return False
        return True

    @property
    def reference(self):
        """
        Return the reference sequence for the group

        """
        return self._barcode

    @property
    def sequences(self):
        """
        List the barcode sequences in the group

        """
        return [b for b in self._sequences]

    @property
    def counts(self):
        """
        Return the total number of counts for all sequences

        """
        return self._counts

    def __len__(self):
        return len(self._sequences)

class SampleSheetBarcodes(object):
    """
    Class for handling index sequence information from a sample sheet

    Given a SampleSheet.csv file this class can extract
    the index sequences (aka barcodes) corresponding to
    sample names, and provides methods to look up one from
    the other.

    Note that the sequences are 'normalised' i.e. dual
    indexes are concatenated with no intermediate character
    (so 'AGCCCTT' and 'GTTACAT' becomes 'AGCCCTTGTTACAT',
    not 'AGCCCTT-GTTACAT' or 'AGCCCTT+GTTACAT')

    Create an initial lookup object:

    >>> s = SampleSheetBarcodes('SampleSheet.csv')

    Get a list of barcodes in lane 2:

    >>> s.barcodes(1)

    Look up the sample name in lane 2 matching a barcode:

    >>> s.lookup_sample('ATTGTG',2)

    Look up the sample name matching a barcode in all lanes
    (e.g. if the sample sheet doesn't include explicit lane
    information):

    >>> s.lookup_sample('ATTGTG')

    """
    def __init__(self,sample_sheet_file):
        """
        Create a new SampleSheetBarcodes instance

        Arguments:
          sample_sheet_file (str): path of a SampleSheet.csv
            file

        """
        self._sample_sheet = SampleSheet(sample_sheet_file)
        self._sample_lookup = {}
        self._barcode_lookup = {}
        self._lanes = []
        sample_id = self._sample_sheet.sample_id_column
        for line in self._sample_sheet.data:
            if self._sample_sheet.has_lanes:
                lane = line['Lane']
            else:
                lane = None
            if lane not in self._lanes:
                self._lanes.append(lane)
                self._sample_lookup[lane] = {}
                self._barcode_lookup[lane] = {}
            sample = line[sample_id]
            index_seq = normalise_barcode(samplesheet_index_sequence(line))
            self._sample_lookup[lane][index_seq] = sample
            self._barcode_lookup[lane][sample] = index_seq

    def barcodes(self,lane=None):
        """
        Return a list of index sequences

        If a lane is specified then a list of normalised
        barcode sequences for that lane is returned; if no
        lane is specified (or lane is 'None') then all
        barcode sequences are returned.

        Arguments:
          lane (int): lane to restrict barcodes to, or None
            to get all barcode sequences

        """
        if not self._sample_sheet.has_lanes:
            # Special case: sample sheet doesn't define any lanes
            lane = None
        if lane in self._lanes:
            barcodes = self._sample_lookup[lane].keys()
        elif lane is None:
            barcodes = []
            for l in self._lanes:
                barcodes.extend(self.barcodes(l))
        else:
            raise KeyError("Lane %s not in sample sheet" % lane)
        return sorted(barcodes)

    def samples(self,lane=None):
        """
        Return a list of sample names

        If a lane is specified then a list of sample names
        that lane is returned; if no lane is specified (or
        lane is 'None') then all barcode sequences are
        returned.

        Arguments:
          lane (int): lane to restrict sample names to, or
            None to get all sample names

        """
        if not self._sample_sheet.has_lanes:
            # Special case:sample sheet doesn't define any lanes
            lane = None
        if lane in self._lanes:
            samples = self._barcode_lookup[lane].keys()
        elif lane is None:
            samples = []
            for l in self._lanes:
                samples.extend(self.samples(l))
        else:
            raise KeyError("Lane %s not in sample sheet" % lane)
        return sorted(samples)

    def lookup_sample(self,barcode,lane=None):
        """
        Return sample name matching barcode sequence

        Arguments:
          barcode (str): normalised barcode sequence
            to get sample name for
          lane (int): optional, lane to look for
            matching sample in

        """
        if not self._sample_sheet.has_lanes:
            # Special case:sample sheet doesn't define any lanes
            lane = None
        if lane in self._lanes:
            return self._sample_lookup[lane][barcode]
        elif lane is None:
            samples = []
            for l in self._lanes:
                try:
                    samples.append(self.lookup_sample(barcode,l))
                except KeyError:
                    pass
            return ','.join(samples)
        else:
            raise KeyError("Lane %s not in sample sheet" % lane)

    def lookup_barcode(self,sample,lane=None):
        """
        Return normalised barcode sequence matching
        sample name

        Arguments:
          sample (str): sample name to get normalised
            barcode sequence for
          lane (int): optional, lane to look for
            matching barcode in

        """
        if not self._sample_sheet.has_lanes:
            # Special case:sample sheet doesn't define any lanes
            lane = None
        if lane in self._lanes:
            return self._barcode_lookup[lane][sample]
        elif lane is None:
            barcodes = []
            for l in self._lanes:
                try:
                    barcodes.append(self.lookup_barcode(sample,l))
                except KeyError:
                    pass
            return ','.join(barcodes)
        else:
            raise KeyError("Lane %s not in sample sheet" % lane)

class Reporter(object):
    """
    Class for generating reports of barcode statistics

    Add arbitrary blocks of text with optional keyword
    'tags', which can then be written to a text file,
    stream, or as an XLS file.

    Usage:

    Make a new Reporter:

    >>> r = Reporter()

    Add a title:

    >>> r.add("This is the title",title=True)

    Add a heading:

    >>> r.add("A heading",heading=True)

    Add some text:

    >>> r.add("Lorem ipsum")

    Write to file:

    >>> r.write("report.txt")

    Write as XLS:

    >>> r.write_xls("report.xls")

    Write as HTML:

    >>> r.write_html("report.html",title="Barcodes")

    """
    def __init__(self):
        """
        Create new Reporter instance
        """
        self._content = []

    def __len__(self):
        return len(self._content)

    def __nonzero__(self):
        return bool(self._content)

    def __str__(self):
        text = []
        for item in self._content:
            content = item[0]
            attrs = item[1]
            if attrs.get('title',False):
                content = make_title(content,'=')
            text.append("%s" % content)
        return '\n'.join(text)

    @property
    def has_warnings(self):
        """
        Report whether warnings were found in analysis
        """
        for item in self._content:
            attrs = item[1]
            if attrs.get('warning',False):
                return True
        return False

    def add(self,content,**kws):
        """
        Add content to the report

        Supplied content is appended to the existing
        content.

        Also arbitrary keyword-value parts can be
        associated with the content.

        """
        for line in content.split('\n'):
            self._content.append((line,dict(**kws)))

    def write(self,fp=None,filen=None,title=None):
        """
        Write the report to a file or stream

        """
        if title is None:
            title = "Barcodes Report"
        if fp is None:
            if filen is None:
                fp = sys.stdout
            else:
                fp = open(filen,'w')
        fp.write("%s\n\n" % make_title(title,'*'))
        if self.has_warnings:
            fp.write("*** %s ***\n\n" % PROBLEMS_DETECTED_TEXT)
        for item in self._content:
            content = item[0]
            attrs = item[1]
            if attrs.get('title',False):
                content = make_title(content,'=')
            fp.write("%s\n" % content)
        if filen is not None:
            fp.close()

    def write_xls(self,xls_file,title=None):
        """
        Write the report to an XLS file

        """
        if title is None:
            title = "Barcodes Report"
        wb = XLSWorkBook(title)
        ws = wb.add_work_sheet("barcodes")
        for item in self._content:
            content = item[0]
            attrs = item[1]
            style = None
            if attrs.get('title',False):
                style = XLSStyle(bold=True,
                                 color='white',
                                 bgcolor='gray50')
            elif attrs.get('heading',False):
                style = XLSStyle(bold=True,
                                 bgcolor='gray25')
            elif attrs.get('strong',False):
                style = XLSStyle(bold=True)
            elif attrs.get('warning',False):
                style = XLSStyle(bold=True,
                                 color='red')
            ws.append_row(data=content.split('\t'),style=style)
        wb.save_as_xls(xls_file)

    def write_html(self,html_file,title=None,no_styles=False):
        """
        Write the report to a HTML file
        """
        if title is None:
            title = "Barcodes Report"
        html = Document(title)
        if self.has_warnings:
            warnings = html.add_section(css_classes=('warnings',))
            warnings.add(Para(WarningIcon(size=50),PROBLEMS_DETECTED_TEXT))
        toc = html.add_section(title="Contents",name="toc")
        toc_list = List()
        toc.add(toc_list)
        section = None
        table = None
        lst = None
        for item in self._content:
            content = item[0]
            attrs = item[1]
            style = None
            # Check if it's tabular data
            is_tabular = (len(content.split('\t')) > 1)
            if is_tabular:
                # Deal with table data
                items = content.lstrip('\t').split('\t')
                if table is None:
                    # New table
                    header = items
                    table = Table(columns=header)
                    section = html.add_section()
                    section.add(table)
                else:
                    # Append to existing table
                    if attrs.get('warning',False):
                        items[0] = Para(WarningIcon(size=20),items[0])
                    table.add_row(**dict(zip(header,items)))
            elif content.startswith(" * "):
                # List
                item = content[3:]
                if lst is None:
                    # New list
                    lst = List()
                    section.add(lst)
                lst.add_item(item)
            else:
                # Not a table or a list
                if attrs.get('title',False):
                    # New section with title
                    section = html.add_section(title=content)
                    if attrs.get('warning',False):
                        toc_list.add_item(WarningIcon(size=20),
                                          Link(section.title,section))
                    else:
                        toc_list.add_item(Link(section.title,section))
                    continue
                if attrs.get('warning',False):
                    section = html.add_section(css_classes=('warning',))
                    section.add(Para(WarningIcon(),content))
                    continue
                if table is not None:
                    # New section after table (no title)
                    section = html.add_section()
                    table = None
                if lst is not None:
                    # Clear list
                    lst = None
                if content:
                    section.add(content)
        # Add styles
        if not no_styles:
            html.add_css_rule("h1 { background-color: #42AEC2;\n"
                              "     color: white;\n"
                              "     padding: 5px 10px; }")
            html.add_css_rule("h2 { background-color: #8CC63F;\n"
                              "     color: white;\n"
                              "     display: inline-block;\n"
                              "     padding: 5px 15px;\n"
                              "     margin: 0;\n"
                              "     border-top-left-radius: 20px;\n"
                              "     border-bottom-right-radius: 20px; }")
            html.add_css_rule("table { border: solid 1px grey;\n"
                              "        font-size: 80%;\n"
                              "        font-family: sans-serif;\n"
                              "        margin: 5px 5px 20px 20px; }")
            html.add_css_rule("table th { background-color: grey;\n"
                              "           color: white;\n"
                              "           padding: 2px 5px; }")
            html.add_css_rule("table td { text-align: right;\n"
                              "           padding: 2px 5px;\n"
                              "           border-bottom: solid 1px lightgray; }")
            html.add_css_rule("div.warning  { padding: 5px;\n"
                              "               border: solid 1px red;\n"
                              "               background-color: F5BCA9;\n"
                              "               color: red;\n"
                              "               font-weight: bold;\n"
                              "               border-radius: 10px;\n"
                              "               display: inline-block; }")
            html.add_css_rule("div.warnings { padding: 2px;\n"
                              "               border: solid 3px red;\n"
                              "               background-color: F5BCA9;\n"
                              "               color: red;\n"
                              "               font-weight: bold;\n"
                              "               margin: 10px;\n"
                              "               border-radius: 10px;\n"
                              "               display: inline-block; }")
            html.add_css_rule("img { vertical-align: middle; }")
        # Write to file
        html.write(html_file)

#######################################################################
# Functions
#######################################################################

def report_barcodes(counts,lane=None,sample_sheet=None,cutoff=None,
                    mismatches=0,reporter=None):
    """
    Report barcode statistics

    Arguments:
      counts (BarcodeCounter): BarcodeCounter instance
        with barcode data to report on
      lane (integer): optional, restrict report to the
        specified lane ('None' reports across all lanes)
      sample_sheet (str): optional, sample sheet file
        to match actual barcodes/groups against
      mismatches (integer): optional, maximum number of
        mismatches to allow for grouping; zero means there
        will be no grouping (default)
      cutoff (float): optional, decimal fraction
        representing minimum percentage of total reads
        that must be associated with a barcode in order
        to be included in analyses (e.g. 0.001 = 0.1%).
        Default is to include all barcodes
      reporter (Reporter): Reporter instance to write
        results to for reporting (optional, default is to
        write to stdout)

    Returns:
      Reporter: Reporter instance
    """
    # Initialise report content
    if reporter is None:
        reporter = Reporter()
    # Get analysis
    analysis = counts.analyse(lane=lane,cutoff=cutoff,
                              sample_sheet=sample_sheet,
                              mismatches=mismatches)
    # Check for overrepresented sequences, and missing and
    # underrepresented samples
    underrepresented = []
    overrepresented = []
    missing = []
    if sample_sheet is not None:
        sample_sheet = SampleSheetBarcodes(sample_sheet)
        found_samples = filter(lambda s: s is not None,
                               [analysis.counts[bc].sample
                                for bc in analysis.barcodes])
        # Get the underrepresented and missing sample names
        for sample in sample_sheet.samples(lane):
            if sample not in found_samples:
                barcode = sample_sheet.lookup_barcode(sample,lane)
                try:
                    count = counts.counts(barcode,lane)
                except KeyError:
                    count = 0
                if count > 0:
                    underrepresented.append(
                        {
                            'name': sample,
                            'barcode': barcode,
                            'counts': count,
                        })
                else:
                    missing.append(
                        {
                            'name': sample,
                            'barcode': barcode,
                        })
        # Get the overrepresented barcodes not associated
        # with a sample name
        lowest_ranked_sample = None
        unassigned_indexes = {}
        for i,barcode in enumerate(analysis.barcodes):
            sample_name = analysis.counts[barcode].sample
            if sample_name is not None:
                lowest_ranked_sample = i
            else:
                unassigned_indexes[barcode] = i
        if lowest_ranked_sample is not None:
            for barcode in unassigned_indexes:
                if unassigned_indexes[barcode] < lowest_ranked_sample:
                    overrepresented.append(
                        {
                            'barcode': barcode,
                            'rank': unassigned_indexes[barcode],
                            'counts': counts.counts(barcode,lane),
                        })
    # Add separator line if the reporter already has content
    if reporter:
        reporter.add('')
    # Check lanes
    warning = (missing or underrepresented or overrepresented)
    if lane is not None:
        reporter.add("Barcode analysis for lane #%d" % lane,
                     title=True,
                     warning=warning)
    else:
        reporter.add("Barcode analysis for all lanes",
                     title=True,
                     warning=warning)
    # Report settings
    if cutoff is not None:
        reporter.add(" * Barcodes which cover less than %.1f%% of reads "
                     "have been excluded" % (cutoff*100.0))
        reporter.add(" * Reported barcodes cover %.1f%% of the data "
                     "(%d/%d)" % (percent(analysis.coverage,
                                          analysis.total_reads),
                                  analysis.coverage,
                                  analysis.total_reads))
    if mismatches:
        reporter.add(" * Barcodes have been grouped by allowing %d "
                     "mismatch%s" % (mismatches,
                                     ('' if mismatches == 1 else 'es')))
    else:
        reporter.add(" * No mismatches were allowed (exact matches only)")
    # Check there are results
    if analysis.total_reads == 0:
        reporter.add("No barcodes counted",warning=True)
        return reporter
    # Warning about specific problems
    if underrepresented or missing or overrepresented:
        reporter.add("")
        reporter.add("Problems detected:",warning=True)
        if underrepresented:
            reporter.add(" * Underrepresented samples")
        if missing:
            reporter.add(" * Missing samples")
        if overrepresented:
            reporter.add(" * Overrepresented barcodes not "
                         "assigned to any samples")
    # Report information on the top barcodes
    cumulative_reads = 0
    reporter.add("")
    reporter.add("%s" % '\t'.join(("#Rank",
                                   "Index",
                                   "Sample",
                                   "N_seqs",
                                   "N_reads",
                                   "%reads",
                                   "(%Total_reads)")),heading=True)
    for i,barcode in enumerate(analysis.barcodes):
        cumulative_reads += analysis.counts[barcode].reads
        sample_name = analysis.counts[barcode].sample
        if sample_name is None:
            sample_name = ''
        reporter.add(
            "%s" % '\t'.join(
                [
                    str(x) for x in
                    ('% 5d' % (i+1),
                     barcode,
                     sample_name,
                     analysis.counts[barcode].sequences,
                     analysis.counts[barcode].reads,
                     '%.1f%%' % (percent(analysis.counts[barcode].reads,
                                         analysis['total_reads'])),
                     '(%.1f%%)' % (percent(cumulative_reads,
                                           analysis['total_reads'])))
                ]),
            warning=(barcode in [b['barcode'] for b in overrepresented]))
    # Report underrepresented samples
    if underrepresented:
        # Sort into order of highest to lowest counts
        underrepresented = sorted(underrepresented,
                                  key=lambda x: x['counts'],
                                  reverse=True)
        # Report
        reporter.add("")
        reporter.add("The following samples are underrepresented:")
        reporter.add("")
        reporter.add("\t#Sample\tIndex\tN_reads\t%reads",heading=True)
        for sample in underrepresented:
            reporter.add("\t%s\t%s\t%d\t%.2f%%" %
                         (sample['name'],
                          sample['barcode'],
                          sample['counts'],
                          percent(sample['counts'],
                                  analysis['total_reads'])))
    # Report missing samples
    if missing:
        # Sort into alphabetical order
        missing = sorted(missing,
                         key=lambda x: x['name'])
        # Report
        reporter.add("")
        reporter.add("The following samples had no counts:")
        reporter.add("")
        reporter.add("\t#Sample\tIndex",heading=True)
        for sample in missing:
            reporter.add("\t%s\t%s" % (sample['name'],
                                       sample['barcode']))
    # Report overrepresented sequences
    if overrepresented:
        reporter.add("")
        reporter.add("The following unassigned barcodes are "
                     "overrepresented compared to the assigned "
                     "barcodes:")
        reporter.add("#Index\tN_reads\t%reads",heading=True)
        for barcode in overrepresented:
            reporter.add("%s\t%d\t%.2f%%" %
                         (barcode['barcode'],
                          barcode['counts'],
                          percent(barcode['counts'],
                                  analysis['total_reads'])))
    return reporter

def detect_barcodes_warnings(report_file):
    """
    Look for warning text in barcode.report file

    Arguments:
      report_file (str): path to barcode report file

    Returns:
      Boolean: True if warnings were found, False if not.
    """
    with open(report_file) as fp:
        for line in fp:
            if PROBLEMS_DETECTED_TEXT in line:
                return True
                break
    return False

def make_title(text,underline="="):
    """
    Turn a string into a Markdown/rst title

    Arguments:
      text (str): text to make into title
      underline (str): underline character (defaults to
        '=')

    Returns:
      String: title text.
    """
    return "%s\n%s" % (text,underline*len(text))

def percent(num,denom):
    """
    Return values as percentage

    Arguments:
      num (float): number to express as percentage
      denom (float): denominator

    Returns:
      Float: value expressed as a percentage.
    """
    return float(num)/float(denom)*100.0
