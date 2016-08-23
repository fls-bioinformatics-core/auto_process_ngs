#!/usr/bin/env python
#
#     analyse_barcodes.py: analyse index sequences from Illumina FASTQs
#     Copyright (C) University of Manchester 2016 Peter Briggs
#
"""
analyse_barcodes.py

Counts the index sequences (i.e. barcodes) for all reads in one or more
Fastq files, and reports the most numerous.

"""

import optparse
import sys
import os
from itertools import izip
from bcftbx.IlluminaData import IlluminaData
from bcftbx.IlluminaData import SampleSheet
from bcftbx.IlluminaData import IlluminaDataError
from bcftbx.IlluminaData import samplesheet_index_sequence
from bcftbx.IlluminaData import normalise_barcode
from bcftbx.FASTQFile import FastqIterator
from bcftbx.utils import AttributeDictionary
from bcftbx.simple_xls import XLSWorkBook
from bcftbx.simple_xls import XLSStyle

__version__ = "0.0.1"

class BarcodeCounter(object):
    """
    Utility class to track counts of barcode sequences

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

    def filter_barcodes(self,cutoff=None,coverage=None,lane=None):
        """
        Return subset of index sequences filtered by specified criteria

        Arguments:
          cutoff (float): barcodes must account for at least this
            fraction of all reads to be included. Total reads are
            all reads if lane is 'None', or else total for the
            specified lane
          coverage (float): barcodes must be in top fraction of
            reads up to this limit
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
        # Limit coverage if specified
        # i.e. only include reads up to the specified total
        # fraction of reads
        if coverage is not None:
            limit = int(float(nreads)*coverage)
            cumulative = 0
            filtered = []
            for s in bc:
                cumulative += self.counts(s,lane)
                filtered.append(s)
                if cumulative >= limit:
                    break
            bc = filtered
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

    def group(self,lane,mismatches=2,n=None,cutoff=None,coverage=None,
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
          coverage: include groups to cover up to this fraction of
            all reads (set to None to cover all reads)
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
        # Cutoff and coverage
        if cutoff is not None:
            cutoff_reads = int(float(nreads)*cutoff)
        else:
            cutoff_reads = 0
        if coverage is not None:
            limit = int(float(nreads)*coverage)
        else:
            limit = 0
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
                # Check if coverage limit has been exceeded
                cumulative += group.counts
                if limit and cumulative >= limit:
                    break
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

        Arguments:
          seq (str): sequence to check against the reference
          mismatches (int): maximum number of mismatches that
            are allowed for the sequences to be considered as
            related (default is 2). Note that 'N's in either
            sequence automatically count as a mismatch.

        """
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
    Class for index sequence information from a sample sheet

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

    def write(self,fp=None,filen=None):
        """
        Write the report to a file or stream

        """
        if fp is None:
            if filen is None:
                fp = sys.stdout
            else:
                fp = open(filen,'w')
        for item in self._content:
            content = item[0]
            attrs = item[1]
            if attrs.get('title',False):
                content = make_title(content,'=')
            fp.write("%s\n" % content)
        if filen is not None:
            fp.close()

    def write_xls(self,xls_file):
        """
        Write the report to an XLS file

        """
        wb = XLSWorkBook("Barcodes Report")
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
            ws.append_row(data=content.split('\t'),style=style)
        wb.save_as_xls(xls_file)

def count_barcodes_bcl2fastq(dirn):
    """
    Count the barcodes from bcl2fastq output

    """
    try:
        unaligned = os.path.basename(dirn.rstrip(os.sep))
        dirn = os.path.dirname(os.path.abspath(dirn.rstrip(os.sep)))
        illumina_data = IlluminaData(dirn,unaligned_dir=unaligned)
    except IlluminaDataError,ex:
        print "%s: not an Illumina output directory?" % dirn
        raise ex
    fqs = []
    for p in illumina_data.projects:
        for s in p.samples:
            for fq in s.fastq_subset(read_number=1,full_path=True):
                fqs.append(fq)
    if illumina_data.undetermined:
        for s in illumina_data.undetermined.samples:
            for fq in s.fastq_subset(read_number=1,full_path=True):
                fqs.append(fq)
    return count_barcodes(fqs)

def count_barcodes(fastqs):
    """
    Count the barcodes from multiple fastqs

    """
    print "Reading in %s fastq%s" % (len(fastqs),
                                     ('' if len(fastqs) == 1
                                      else 's'))
    counts = BarcodeCounter()
    for fq in fastqs:
        print "%s" % os.path.basename(fq)
        for r in FastqIterator(fq):
            seq = r.seqid.index_sequence
            lane = int(r.seqid.flowcell_lane)
            counts.count_barcode(seq,lane)
    return counts

def make_title(text,underline="="):
    return "%s\n%s" % (text,underline*len(text))

def percent(num,denom):
    return float(num)/float(denom)*100.0

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
      coverage (float): optional, decimal fraction
        representing percentage of reads that will be
        included in analyses (e.g. 0.9 = 90%). Default
        is to include all barcodes
      reporter (Reporter): Reporter instance to write
        results to for reporting (optional, default is to
        write to stdout)

    """
    # Initialise report content
    if reporter is None:
        reporter = Reporter()
    # Add separator line if the reporter already has content
    if reporter:
        reporter.add('')
    # Check lanes
    if lane is not None:
        reporter.add("Barcode analysis for lane #%d" % lane,title=True)
    else:
        reporter.add("Barcode analysis for all lanes",title=True)
    # Get analysis
    analysis = counts.analyse(lane=lane,cutoff=cutoff,
                              sample_sheet=sample_sheet,
                              mismatches=mismatches)
    # Report settings
    if cutoff is not None:
        reporter.add("Barcodes which cover less than %.1f%% of reads "
                     "have been excluded" % (cutoff*100.0))
        reporter.add("Reported barcodes cover %.1f%% of the data "
                     "(%d/%d)" % (percent(analysis['coverage'],
                                          analysis['total_reads']),
                                  analysis['coverage'],
                                  analysis['total_reads']))
    if mismatches:
        reporter.add("Barcodes have been grouped by allowing %d mismatch%s" %
                     (mismatches,
                      ('' if mismatches == 1 else 'es')))
    cumulative_reads = 0
    # Report information on the top barcodes
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
        reporter.add("%s" % '\t'.join(
            [str(x) for x in
             ('% 5d' % (i+1),
              barcode,
              sample_name,
              analysis.counts[barcode].sequences,
              analysis.counts[barcode].reads,
              '%.1f%%' % (percent(analysis.counts[barcode].reads,analysis['total_reads'])),
              '(%.1f%%)' % (percent(cumulative_reads,
                                    analysis['total_reads'])))]))
    # Report "missing" samples
    if sample_sheet is not None:
        sample_sheet = SampleSheetBarcodes(sample_sheet)
        found_samples = filter(lambda s: s is not None,
                               [analysis.counts[bc].sample
                                for bc in analysis.barcodes])
        found_samples = ','.join(found_samples).split(',')
        missing = []
        missing_no_counts = []
        for sample in sample_sheet.samples(lane):
            if sample not in found_samples:
                barcode = sample_sheet.lookup_barcode(sample,lane)
                try:
                    missing.append({'name': sample,
                                    'barcode': barcode,
                                    'counts': counts.counts(barcode,lane)
                                })
                except KeyError:
                    missing_no_counts.append({'name': sample,
                                              'barcode': barcode,
                                          })
        if missing:
            # Sort into order of highest to lowest counts
            missing = sorted(missing,key=lambda x: x.counts,
                             reverse=True)
            # Report
            reporter.add("")
            reporter.add("The following samples had too few counts to "
                         "appear in the results:",strong=True)
            reporter.add("")
            reporter.add("\t#Sample\tIndex\tN_reads\t%reads",heading=True)
            for sample in missing:
                reporter.add("\t%s\t%s\t%d\t%.2f%%" %
                             (sample['name'],
                              sample['barcode'],
                              sample.counts,
                              percent(sample.counts,
                                      analysis['total_reads'])))
        if missing_no_counts:
            # Sort into alphabetical order
            missing_no_counts = sorted(missing_no_counts,
                                       key=lambda x: x['name'])
            # Report
            reporter.add("")
            reporter.add("The following samples had no counts:",
                         strong=True)
            reporter.add("")
            reporter.add("\t#Sample\tIndex",heading=True)
            for sample in missing_no_counts:
                reporter.add("\t%s\t%s" % (sample['name'],
                                           sample['barcode']))
    return reporter

def parse_lanes_expression(lanes):
    """
    Parse a string and return list of lane numbers

    Process a string specifying one or more lane numbers
    expressed as a single digit and returns a sorted
    list of unique lane numbers.

    For example:

    >>> parse_lanes_expression('1')
    [1]

    or as a list of digits:

    >>> parse_lanes_expression('1,3')
    [1, 3]

    or as a range of digits:

    >>> parse_lanes_expression('5-7')
    [5, 6, 7]

    or as a mixture of both:

    >>> parse_lanes_expression('1,3,5-7')
    [1, 3, 5, 6, 7]
    
    """
    if lanes is None:
        return None
    lane_numbers = []
    for l in lanes.split(','):
        l1 = l.split('-')
        if len(l1) == 1:
            lane_numbers.append(int(l))
        else:
            lane_numbers.extend(xrange(int(l1[0]),int(l1[1])+1))
    lane_numbers = sorted(set(lane_numbers))
    return lane_numbers

# Unit tests

import unittest
import tempfile
import shutil

# BarcodeCounter
class TestBarcodeCounter(unittest.TestCase):
    def setUp(self):
        # Temporary working dir (if needed)
        self.wd = None

    def _make_working_dir(self):
        if self.wd is None:
            self.wd = tempfile.mkdtemp(suffix='.test_BarcodeCounter')

    def _make_file(self,name,contents):
        # Create a file under the working directory
        # called "name" and populated with "contents"
        # Working directory will be created if not already
        # set up
        # Returns the path to the file
        self._make_working_dir()
        filen = os.path.join(self.wd,name)
        with open(filen,'w') as fp:
            fp.write(contents)
        return filen

    def tearDown(self):
        # Remove temporary working dir
        if self.wd is not None and os.path.isdir(self.wd):
            shutil.rmtree(self.wd)

    def test_empty_counter(self):
        """BarcodeCounter: check empty counter
        """
        # Initialise counter object
        bc = BarcodeCounter()
        self.assertEqual(bc.barcodes(),[])
        self.assertEqual(bc.lanes,[])
        self.assertEqual(bc.filter_barcodes(),[])
        self.assertEqual(bc.counts("AGGCAGAATCTTACGC"),0)
        self.assertEqual(bc.counts("AGGCAGAATCTTACGC",lane=1),0)
        self.assertEqual(bc.counts_all("AGGCAGAATCTTACGC"),0)
        self.assertEqual(bc.nreads(),0)
        self.assertEqual(bc.nreads(1),0)

    def test_count_fastq_sequences(self):
        """BarcodeCounter: count barcode sequences
        """
        # Initialise counter object
        bc = BarcodeCounter()
        # Populate with sequences
        for r,incr in (((1,"AGGCAGAATCTTACGC"),102),
                       ((1,"TCCTGAGCTCTTACGC"),10),
                       ((1,"ACAGTGATTCTTTCCC"),3),
                       ((1,"ATGCTCGTCTCGCATC"),1),
                       ((2,"CGTACTAGTCTTACGC"),95),
                       ((2,"ATGTCAGATCTTTCCC"),29),
                       ((2,"AGGCAGAATCTTACGC"),12),
                       ((2,"CAGATCATTCTTTCCC"),6),
                       ((3,"GGACTCCTTCTTACGC"),75),
                       ((3,"ACCGATTCGCGCGTAG"),74),
                       ((3,"CCAGCAATATCGCGAG"),2),
                       ((3,"CCGCGTAAGCAATAGA"),1)):
            lane,seq = r
            for i in xrange(incr):
                bc.count_barcode(seq,lane=lane)
        # Check contents
        self.assertEqual(bc.barcodes(),["AGGCAGAATCTTACGC",
                                        "CGTACTAGTCTTACGC",
                                        "GGACTCCTTCTTACGC",
                                        "ACCGATTCGCGCGTAG",
                                        "ATGTCAGATCTTTCCC",
                                        "TCCTGAGCTCTTACGC",
                                        "CAGATCATTCTTTCCC",
                                        "ACAGTGATTCTTTCCC",
                                        "CCAGCAATATCGCGAG",
                                        "ATGCTCGTCTCGCATC",
                                        "CCGCGTAAGCAATAGA"])
        # Lanes
        self.assertEqual(bc.lanes,[1,2,3])
        # Counts for individual barcodes
        self.assertEqual(bc.counts("AGGCAGAATCTTACGC"),114)
        self.assertEqual(bc.counts("AGGCAGAATCTTACGC",lane=1),102)
        self.assertEqual(bc.counts("AGGCAGAATCTTACGC",lane=2),12)
        self.assertEqual(bc.counts("AGGCAGAATCTTACGC",lane=3),0)
        self.assertEqual(bc.counts_all("AGGCAGAATCTTACGC"),114)
        self.assertEqual(bc.counts("CCGCGTAAGCAATAGA"),1)
        self.assertEqual(bc.counts("CCGCGTAAGCAATAGA",lane=1),0)
        self.assertEqual(bc.counts("CCGCGTAAGCAATAGA",lane=2),0)
        self.assertEqual(bc.counts("CCGCGTAAGCAATAGA",lane=3),1)
        self.assertEqual(bc.counts_all("CCGCGTAAGCAATAGA"),1)
        # Read counts
        self.assertEqual(bc.nreads(),410)
        self.assertEqual(bc.nreads(1),116)
        self.assertEqual(bc.nreads(2),142)
        self.assertEqual(bc.nreads(3),152)

    def test_filter_barcodes(self):
        """BarcodeCounter: check filtering by lane and cutoff
        """
        bc = BarcodeCounter()
        bc.count_barcode("TATGCGCGGTA",lane=1,incr=285302)
        bc.count_barcode("TATGCGCGGTG",lane=1,incr=532)
        bc.count_barcode("ACCTACCGGTA",lane=1,incr=315)
        bc.count_barcode("CCCTTATGCGA",lane=1,incr=22)
	bc.count_barcode("ACCTAGCGGTA",lane=2,incr=477)
        bc.count_barcode("ACCTCTATGCT",lane=2,incr=368)
        self.assertEqual(bc.barcodes(),["TATGCGCGGTA",
                                        "TATGCGCGGTG",
                                        "ACCTAGCGGTA",
                                        "ACCTCTATGCT",
                                        "ACCTACCGGTA",
                                        "CCCTTATGCGA"])
        # No filtering
        self.assertEqual(bc.filter_barcodes(),["TATGCGCGGTA",
                                               "TATGCGCGGTG",
                                               "ACCTAGCGGTA",
                                               "ACCTCTATGCT",
                                               "ACCTACCGGTA",
                                               "CCCTTATGCGA"])
        # Filter by lane
        self.assertEqual(bc.filter_barcodes(lane=1),["TATGCGCGGTA",
                                                     "TATGCGCGGTG",
                                                     "ACCTACCGGTA",
                                                     "CCCTTATGCGA"]),
        self.assertEqual(bc.filter_barcodes(lane=2),["ACCTAGCGGTA",
                                                     "ACCTCTATGCT"])
        # Filter by cutoff
        self.assertEqual(bc.filter_barcodes(cutoff=0.5),
                         ["TATGCGCGGTA",])
        self.assertEqual(bc.filter_barcodes(cutoff=0.0015,lane=1),
                         ["TATGCGCGGTA","TATGCGCGGTG"])
        self.assertEqual(bc.filter_barcodes(cutoff=0.5,lane=2),
                         ["ACCTAGCGGTA",])

    def test_group(self):
        """BarcodeCounter: check grouping of barcode sequences
        """
        bc = BarcodeCounter()
        bc.count_barcode("TATGCGCGGTA",lane=1,incr=285302)
        bc.count_barcode("CATGCGCGGTA",lane=1,incr=8532)
        bc.count_barcode("GATGCGCGGTA",lane=1,incr=5321)
        bc.count_barcode("GCTGCGCGGTA",lane=1,incr=7853)
        bc.count_barcode("GCTGCGCGGTC",lane=1,incr=325394)
        bc.count_barcode("GTCACGCGGTA",lane=2,incr=296201)
        bc.count_barcode("GTCACGCGGTT",lane=2,incr=2853)
        bc.count_barcode("GTCACGCTGTT",lane=2,incr=278539)
        ## 2 mismatches across all lanes
        groups = bc.group(None,mismatches=2)
        ##"GCTGCGCGGTC","GCTGCGCGGTA","GATGCGCGGTA" = 338568
        ##"TATGCGCGGTA","CATGCGCGGTA" = 293834
        ##"GTCACGCGGTA","GTCACGCTGTT","GTCACGCGGTT" = 577593
        self.assertEqual(len(groups),3)
        self.assertEqual(groups[0].reference,"GTCACGCGGTA")
        self.assertEqual(groups[0].sequences,["GTCACGCGGTA",
                                              "GTCACGCTGTT",
                                              "GTCACGCGGTT"])
        self.assertEqual(groups[0].counts,577593)
        self.assertEqual(groups[1].reference,"GCTGCGCGGTC")
        self.assertEqual(groups[1].sequences,["GCTGCGCGGTC",
                                              "GCTGCGCGGTA",
                                              "GATGCGCGGTA"])
        self.assertEqual(groups[1].counts,338568)
        self.assertEqual(groups[2].reference,"TATGCGCGGTA")
        self.assertEqual(groups[2].sequences,["TATGCGCGGTA",
                                              "CATGCGCGGTA"])
        self.assertEqual(groups[2].counts,293834)
        ## 1 mismatch across all lanes
        groups = bc.group(None,mismatches=1)
        ##"TATGCGCGGTA","CATGCGCGGTA","GATGCGCGGTA" = 299155
        ##"GCTGCGCGGTC","GCTGCGCGGTA" = 333247
        ##"GTCACGCGGTA","GTCACGCGGTT" = 299054
        ##"GTCACGCTGTT" = 278539
        self.assertEqual(len(groups),4)
        self.assertEqual(groups[0].reference,"GCTGCGCGGTC")
        self.assertEqual(groups[0].sequences,["GCTGCGCGGTC",
                                              "GCTGCGCGGTA"])
        self.assertEqual(groups[0].counts,333247)
        self.assertEqual(groups[1].reference,"TATGCGCGGTA")
        self.assertEqual(groups[1].sequences,["TATGCGCGGTA",
                                              "CATGCGCGGTA",
                                              "GATGCGCGGTA"])
        self.assertEqual(groups[1].counts,299155)
        self.assertEqual(groups[2].reference,"GTCACGCGGTA")
        self.assertEqual(groups[2].sequences,["GTCACGCGGTA",
                                              "GTCACGCGGTT"])
        self.assertEqual(groups[2].counts,299054)
        self.assertEqual(groups[3].reference,"GTCACGCTGTT")
        self.assertEqual(groups[3].sequences,["GTCACGCTGTT",])
        self.assertEqual(groups[3].counts,278539)
        ## 1 mismatch in lane 1
        groups = bc.group(1,mismatches=1)
        ##"TATGCGCGGTA","CATGCGCGGTA","GATGCGCGGTA" = 299155
        ##"GCTGCGCGGTC","GCTGCGCGGTA" = 333247
        self.assertEqual(len(groups),2)
        self.assertEqual(groups[0].reference,"GCTGCGCGGTC")
        self.assertEqual(groups[0].sequences,["GCTGCGCGGTC",
                                              "GCTGCGCGGTA"])
        self.assertEqual(groups[0].counts,333247)
        self.assertEqual(groups[1].reference,"TATGCGCGGTA")
        self.assertEqual(groups[1].sequences,["TATGCGCGGTA",
                                              "CATGCGCGGTA",
                                              "GATGCGCGGTA"])
        self.assertEqual(groups[1].counts,299155)
        ## 2 mismatches across all lanes
        groups = bc.group(None,mismatches=2)
        ##"GCTGCGCGGTC","GCTGCGCGGTA","GATGCGCGGTA" = 338568
        ##"TATGCGCGGTA","CATGCGCGGTA" = 293834
        ##"GTCACGCGGTA","GTCACGCTGTT","GTCACGCGGTT" = 577593
        self.assertEqual(len(groups),3)
        self.assertEqual(groups[0].reference,"GTCACGCGGTA")
        self.assertEqual(groups[0].sequences,["GTCACGCGGTA",
                                              "GTCACGCTGTT",
                                              "GTCACGCGGTT"])
        self.assertEqual(groups[0].counts,577593)
        self.assertEqual(groups[1].reference,"GCTGCGCGGTC")
        self.assertEqual(groups[1].sequences,["GCTGCGCGGTC",
                                              "GCTGCGCGGTA",
                                              "GATGCGCGGTA"])
        self.assertEqual(groups[1].counts,338568)
        self.assertEqual(groups[2].reference,"TATGCGCGGTA")
        self.assertEqual(groups[2].sequences,["TATGCGCGGTA",
                                              "CATGCGCGGTA"])
        self.assertEqual(groups[2].counts,293834)

    def test_analyse(self):
        """BarcodeCounter: perform analysis with defaults
        """
        bc = BarcodeCounter()
        bc.count_barcode("TATGCGCGGTA",lane=1,incr=285302)
        bc.count_barcode("CATGCGCGGTA",lane=1,incr=8532)
        bc.count_barcode("GATGCGCGGTA",lane=1,incr=5321)
        bc.count_barcode("GCTGCGCGGTA",lane=1,incr=7853)
        bc.count_barcode("GCTGCGCGGTC",lane=1,incr=325394)
        analysis = bc.analyse(lane=1)
        self.assertEqual(analysis.cutoff,None)
        self.assertEqual(analysis.mismatches,0)
        self.assertEqual(analysis.total_reads,632402)
        self.assertEqual(analysis.coverage,632402)
        self.assertEqual(analysis.barcodes,["GCTGCGCGGTC",
                                            "TATGCGCGGTA",
                                            "CATGCGCGGTA",
                                            "GCTGCGCGGTA",
                                            "GATGCGCGGTA"])
        self.assertEqual(analysis.counts["GCTGCGCGGTC"].reads,325394)
        self.assertEqual(analysis.counts["TATGCGCGGTA"].reads,285302)
        self.assertEqual(analysis.counts["CATGCGCGGTA"].reads,8532)
        self.assertEqual(analysis.counts["GCTGCGCGGTA"].reads,7853)
        self.assertEqual(analysis.counts["GATGCGCGGTA"].reads,5321)
        self.assertEqual(analysis.counts["GCTGCGCGGTC"].sample,None)
        self.assertEqual(analysis.counts["TATGCGCGGTA"].sample,None)
        self.assertEqual(analysis.counts["CATGCGCGGTA"].sample,None)
        self.assertEqual(analysis.counts["GCTGCGCGGTA"].sample,None)
        self.assertEqual(analysis.counts["GATGCGCGGTA"].sample,None)
        self.assertEqual(analysis.counts["GCTGCGCGGTC"].sequences,1)
        self.assertEqual(analysis.counts["TATGCGCGGTA"].sequences,1)
        self.assertEqual(analysis.counts["CATGCGCGGTA"].sequences,1)
        self.assertEqual(analysis.counts["GCTGCGCGGTA"].sequences,1)
        self.assertEqual(analysis.counts["GATGCGCGGTA"].sequences,1)

    def test_analyse_with_cutoff(self):
        """BarcodeCounter: perform analysis with cutoff
        """
        bc = BarcodeCounter()
        bc.count_barcode("TATGCGCGGTA",lane=1,incr=285302)
        bc.count_barcode("CATGCGCGGTA",lane=1,incr=8532)
        bc.count_barcode("GATGCGCGGTA",lane=1,incr=5321)
        bc.count_barcode("GCTGCGCGGTA",lane=1,incr=7853)
        bc.count_barcode("GCTGCGCGGTC",lane=1,incr=325394)
        analysis = bc.analyse(lane=1,cutoff=0.013)
        self.assertEqual(analysis.cutoff,0.013)
        self.assertEqual(analysis.mismatches,0)
        self.assertEqual(analysis.total_reads,632402)
        self.assertEqual(analysis.coverage,619228)
        self.assertEqual(analysis.barcodes,["GCTGCGCGGTC",
                                            "TATGCGCGGTA",
                                            "CATGCGCGGTA"])
        self.assertEqual(analysis.counts["GCTGCGCGGTC"].reads,325394)
        self.assertEqual(analysis.counts["TATGCGCGGTA"].reads,285302)
        self.assertEqual(analysis.counts["CATGCGCGGTA"].reads,8532)
        self.assertEqual(analysis.counts["GCTGCGCGGTC"].sample,None)
        self.assertEqual(analysis.counts["TATGCGCGGTA"].sample,None)
        self.assertEqual(analysis.counts["CATGCGCGGTA"].sample,None)
        self.assertEqual(analysis.counts["GCTGCGCGGTC"].sequences,1)
        self.assertEqual(analysis.counts["TATGCGCGGTA"].sequences,1)
        self.assertEqual(analysis.counts["CATGCGCGGTA"].sequences,1)

    def test_analyse_with_sample_sheet(self):
        """BarcodeCounter: perform analysis with samplesheet
        """
        # Create sample sheet
        sample_sheet_file = self._make_file("SampleSheet.csv",
                                            """[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description
1,SMPL1,,,,A006,CATGCGCGGTA,,
1,SMPL2,,,,A012,GCTGCGCGGTC,,
2,SMPL3,,,,A005,ACAGTGCGGTA,,
2,SMPL4,,,,A019,GTGAAACGGTC,,
""")
        # Set up barcode counts
        bc = BarcodeCounter()
        bc.count_barcode("TATGCGCGGTA",lane=1,incr=285302)
        bc.count_barcode("CATGCGCGGTA",lane=1,incr=8532)
        bc.count_barcode("GATGCGCGGTA",lane=1,incr=5321)
        bc.count_barcode("GCTGCGCGGTA",lane=1,incr=7853)
        bc.count_barcode("GCTGCGCGGTC",lane=1,incr=325394)
        analysis = bc.analyse(lane=1,sample_sheet=sample_sheet_file)
        self.assertEqual(analysis.cutoff,None)
        self.assertEqual(analysis.mismatches,0)
        self.assertEqual(analysis.total_reads,632402)
        self.assertEqual(analysis.coverage,632402)
        self.assertEqual(analysis.barcodes,["GCTGCGCGGTC",
                                            "TATGCGCGGTA",
                                            "CATGCGCGGTA",
                                            "GCTGCGCGGTA",
                                            "GATGCGCGGTA"])
        self.assertEqual(analysis.counts["GCTGCGCGGTC"].reads,325394)
        self.assertEqual(analysis.counts["TATGCGCGGTA"].reads,285302)
        self.assertEqual(analysis.counts["CATGCGCGGTA"].reads,8532)
        self.assertEqual(analysis.counts["GCTGCGCGGTA"].reads,7853)
        self.assertEqual(analysis.counts["GATGCGCGGTA"].reads,5321)
        self.assertEqual(analysis.counts["GCTGCGCGGTC"].sample,"SMPL2")
        self.assertEqual(analysis.counts["TATGCGCGGTA"].sample,None)
        self.assertEqual(analysis.counts["CATGCGCGGTA"].sample,"SMPL1")
        self.assertEqual(analysis.counts["GCTGCGCGGTA"].sample,None)
        self.assertEqual(analysis.counts["GATGCGCGGTA"].sample,None)
        self.assertEqual(analysis.counts["GCTGCGCGGTC"].sequences,1)
        self.assertEqual(analysis.counts["TATGCGCGGTA"].sequences,1)
        self.assertEqual(analysis.counts["CATGCGCGGTA"].sequences,1)
        self.assertEqual(analysis.counts["GCTGCGCGGTA"].sequences,1)
        self.assertEqual(analysis.counts["GATGCGCGGTA"].sequences,1)

    def test_analyse_groups(self):
        """BarcodeCounter: perform analysis with grouping
        """
        bc = BarcodeCounter()
        bc.count_barcode("TATGCGCGGTA",lane=1,incr=285302)
        bc.count_barcode("CATGCGCGGTA",lane=1,incr=8532)
        bc.count_barcode("GATGCGCGGTA",lane=1,incr=5321)
        bc.count_barcode("GCTGCGCGGTA",lane=1,incr=7853)
        bc.count_barcode("GCTGCGCGGTC",lane=1,incr=325394)
        analysis = bc.analyse(lane=1,mismatches=1)
        ##"TATGCGCGGTA","CATGCGCGGTA","GATGCGCGGTA" = 299155
        ##"GCTGCGCGGTC","GCTGCGCGGTA" = 333247
        self.assertEqual(analysis.cutoff,None)
        self.assertEqual(analysis.mismatches,1)
        self.assertEqual(analysis.total_reads,632402)
        self.assertEqual(analysis.coverage,632402)
        self.assertEqual(analysis.barcodes,["GCTGCGCGGTC",
                                            "TATGCGCGGTA"])
        self.assertEqual(analysis.counts["GCTGCGCGGTC"].reads,333247)
        self.assertEqual(analysis.counts["TATGCGCGGTA"].reads,299155)
        self.assertEqual(analysis.counts["GCTGCGCGGTC"].sample,None)
        self.assertEqual(analysis.counts["TATGCGCGGTA"].sample,None)
        self.assertEqual(analysis.counts["GCTGCGCGGTC"].sequences,2)
        self.assertEqual(analysis.counts["TATGCGCGGTA"].sequences,3)

    def test_analyse_groups_with_sample_sheet(self):
        """BarcodeCounter: perform analysis with grouping and samplesheet
        """
        # Create sample sheet
        sample_sheet_file = self._make_file("SampleSheet.csv",
                                            """[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description
1,SMPL1,,,,A006,CATGCGCGGTA,,
1,SMPL2,,,,A012,GCTGCGCGGTC,,
2,SMPL3,,,,A005,ACAGTGCGGTA,,
2,SMPL4,,,,A019,GTGAAACGGTC,,
""")
        # Set up barcode counts
        bc = BarcodeCounter()
        bc.count_barcode("TATGCGCGGTA",lane=1,incr=285302)
        bc.count_barcode("CATGCGCGGTA",lane=1,incr=8532)
        bc.count_barcode("GATGCGCGGTA",lane=1,incr=5321)
        bc.count_barcode("GCTGCGCGGTA",lane=1,incr=7853)
        bc.count_barcode("GCTGCGCGGTC",lane=1,incr=325394)
        analysis = bc.analyse(lane=1,
                              mismatches=2,
                              sample_sheet=sample_sheet_file)
        ##"CATGCGCGGTA","TATGCGCGGTA","GATGCGCGGTA","GCTGCGCGGTA" = 307008
        ##"GCTGCGCGGTC" = 325394
        self.assertEqual(analysis.cutoff,None)
        self.assertEqual(analysis.mismatches,2)
        self.assertEqual(analysis.total_reads,632402)
        self.assertEqual(analysis.coverage,632402)
        self.assertEqual(analysis.barcodes,["GCTGCGCGGTC",
                                            "CATGCGCGGTA"])
        self.assertEqual(analysis.counts["GCTGCGCGGTC"].reads,325394)
        self.assertEqual(analysis.counts["CATGCGCGGTA"].reads,307008)
        self.assertEqual(analysis.counts["GCTGCGCGGTC"].sample,"SMPL2")
        self.assertEqual(analysis.counts["CATGCGCGGTA"].sample,"SMPL1")
        self.assertEqual(analysis.counts["GCTGCGCGGTC"].sequences,1)
        self.assertEqual(analysis.counts["CATGCGCGGTA"].sequences,4)

    def test_read_counts_file(self):
        """BarcodeCounter: read in data from '.counts' file
        """
        # Read a counts file
        counts_file = self._make_file("test.counts","""#Lane	Rank	Sequence	Count
1	1	TATGCGCGGTA	285302
1	2	TATGCGCGGTG	532
1	3	ACCTACCGGTA	315
1	4	CCCTTATGCGA	22
2	5	ACCTAGCGGTA	477
2	6	ACCTCTATGCT	368
3	7	ACCCTNCGGTA	312
3	8	ACCTTATGCGC	248""")
        # Read the file
        bc = BarcodeCounter(counts_file)
        # Check the contents
        self.assertEqual(bc.barcodes(),["TATGCGCGGTA",
                                        "TATGCGCGGTG",
                                        "ACCTAGCGGTA",
                                        "ACCTCTATGCT",
                                        "ACCTACCGGTA",
                                        "ACCCTNCGGTA",
                                        "ACCTTATGCGC",
                                        "CCCTTATGCGA"])
        # Lanes
        self.assertEqual(bc.lanes,[1,2,3])
        # Counts for individual barcodes
        self.assertEqual(bc.counts("TATGCGCGGTA"),285302)
        self.assertEqual(bc.counts("TATGCGCGGTG"),532)
        self.assertEqual(bc.counts("ACCTAGCGGTA"),477)
        self.assertEqual(bc.counts("ACCTCTATGCT"),368)
        self.assertEqual(bc.counts("ACCTACCGGTA"),315)
        self.assertEqual(bc.counts("ACCCTNCGGTA"),312)
        self.assertEqual(bc.counts("ACCTTATGCGC"),248)
        self.assertEqual(bc.counts("CCCTTATGCGA"),22)
        self.assertEqual(bc.counts("TATGCGCGGTA",lane=1),285302)
        self.assertEqual(bc.counts("TATGCGCGGTA",lane=2),0)
        self.assertEqual(bc.counts("TATGCGCGGTA",lane=3),0)
        self.assertEqual(bc.counts_all("TATGCGCGGTA"),285302)
        self.assertEqual(bc.counts("ACCTTATGCGC",lane=1),0)
        self.assertEqual(bc.counts("ACCTTATGCGC",lane=2),0)
        self.assertEqual(bc.counts("ACCTTATGCGC",lane=3),248)
        self.assertEqual(bc.counts_all("ACCTTATGCGC"),248)
        # Read counts
        self.assertEqual(bc.nreads(),287576)
        self.assertEqual(bc.nreads(1),286171)
        self.assertEqual(bc.nreads(2),845)
        self.assertEqual(bc.nreads(3),560)

    def test_read_multiple_counts_file(self):
        """BarcodeCounter: read in data from multiple '.counts' files
        """
        # Read multiple counts files
        counts_lane1 = self._make_file("lane1.counts",
                                       """#Lane	Rank	Sequence	Count
1	1	TATGCGCGGTA	285302
1	2	TATGCGCGGTG	532
1	3	ACCTACCGGTA	315
1	4	CCCTTATGCGA	22""")
        counts_lane2 = self._make_file("lane2.counts",
                                       """#Lane	Rank	Sequence	Count
2	1	ACCTAGCGGTA	477
2	2	ACCTCTATGCT	368""")
        counts_lane3 = self._make_file("lane3.counts",
                                       """#Lane	Rank	Sequence	Count
3	1	ACCCTNCGGTA	312
3	2	ACCTTATGCGC	248""")
        # Read the file
        bc = BarcodeCounter(counts_lane1,counts_lane2,counts_lane3)
        # Check the contents
        self.assertEqual(bc.barcodes(),["TATGCGCGGTA",
                                        "TATGCGCGGTG",
                                        "ACCTAGCGGTA",
                                        "ACCTCTATGCT",
                                        "ACCTACCGGTA",
                                        "ACCCTNCGGTA",
                                        "ACCTTATGCGC",
                                        "CCCTTATGCGA"])
        # Lanes
        self.assertEqual(bc.lanes,[1,2,3])
        # Counts for individual barcodes
        self.assertEqual(bc.counts("TATGCGCGGTA"),285302)
        self.assertEqual(bc.counts("TATGCGCGGTG"),532)
        self.assertEqual(bc.counts("ACCTAGCGGTA"),477)
        self.assertEqual(bc.counts("ACCTCTATGCT"),368)
        self.assertEqual(bc.counts("ACCTACCGGTA"),315)
        self.assertEqual(bc.counts("ACCCTNCGGTA"),312)
        self.assertEqual(bc.counts("ACCTTATGCGC"),248)
        self.assertEqual(bc.counts("CCCTTATGCGA"),22)
        self.assertEqual(bc.counts("TATGCGCGGTA",lane=1),285302)
        self.assertEqual(bc.counts("TATGCGCGGTA",lane=2),0)
        self.assertEqual(bc.counts("TATGCGCGGTA",lane=3),0)
        self.assertEqual(bc.counts_all("TATGCGCGGTA"),285302)
        self.assertEqual(bc.counts("ACCTTATGCGC",lane=1),0)
        self.assertEqual(bc.counts("ACCTTATGCGC",lane=2),0)
        self.assertEqual(bc.counts("ACCTTATGCGC",lane=3),248)
        self.assertEqual(bc.counts_all("ACCTTATGCGC"),248)
        # Read counts
        self.assertEqual(bc.nreads(),287576)
        self.assertEqual(bc.nreads(1),286171)
        self.assertEqual(bc.nreads(2),845)
        self.assertEqual(bc.nreads(3),560)

    def test_read_old_style_counts_file(self):
        """BarcodeCounter: read in data from old-style 3 column '.counts' file
        """
        # Read old-style 3 column counts files
        self._make_working_dir()
        old_style_counts_file = self._make_file("old_style.counts",
                                                """#Rank	Sequence	Count
1	TATGCGCGGTA	285302
2	TATGCGCGGTG	532
3	ACCTACCGGTA	315
4	CCCTTATGCGA	22""")
        # Read the file
        bc = BarcodeCounter(old_style_counts_file)
        # Check the contents
        self.assertEqual(bc.barcodes(),["TATGCGCGGTA",
                                        "TATGCGCGGTG",
                                        "ACCTACCGGTA",
                                        "CCCTTATGCGA"])
        # Lanes
        self.assertEqual(bc.lanes,[])
        # Counts for individual barcodes
        self.assertEqual(bc.counts("TATGCGCGGTA"),285302)
        self.assertEqual(bc.counts("TATGCGCGGTG"),532)
        self.assertEqual(bc.counts("ACCTACCGGTA"),315)
        self.assertEqual(bc.counts("CCCTTATGCGA"),22)
        self.assertEqual(bc.counts("TATGCGCGGTA",lane=1),0)
        self.assertEqual(bc.counts_all("TATGCGCGGTA"),285302)
        # Read counts
        self.assertEqual(bc.nreads(),286171)

    def test_write_counts_file(self):
        """BarcodeCounter: write counts to a file
        """
        # Write a file
        self._make_working_dir()
        bc = BarcodeCounter()
        bc.count_barcode("TATGCGCGGTA",lane=1,incr=285302)
        bc.count_barcode("TATGCGCGGTG",lane=1,incr=532)
        bc.count_barcode("ACCTACCGGTA",lane=1,incr=315)
        bc.count_barcode("CCCTTATGCGA",lane=1,incr=22)
        bc.count_barcode("ACCTAGCGGTA",lane=2,incr=477)
        bc.count_barcode("ACCTCTATGCT",lane=2,incr=368)
        bc.count_barcode("ACCCTNCGGTA",lane=3,incr=312)
        bc.count_barcode("ACCTTATGCGC",lane=3,incr=248)
        counts_file = os.path.join(self.wd,"out.counts")
        bc.write(counts_file)
        expected_contents = """#Lane	Rank	Sequence	Count
1	1	TATGCGCGGTA	285302
1	2	TATGCGCGGTG	532
1	3	ACCTACCGGTA	315
1	4	CCCTTATGCGA	22
2	1	ACCTAGCGGTA	477
2	2	ACCTCTATGCT	368
3	1	ACCCTNCGGTA	312
3	2	ACCTTATGCGC	248
"""
        self.assertTrue(os.path.exists(counts_file))
        self.assertEqual(open(counts_file,'r').read(),
                         expected_contents)

# BarcodeGroup
class TestBarcodeGroup(unittest.TestCase):
    def test_barcodegroup(self):
        """BarcodeGroup: check making a new instance
        """
        # Create a new BarcodeGroup
        grp = BarcodeGroup("CTAAGCCT",2894178)
        self.assertEqual(grp.reference,"CTAAGCCT")
        self.assertEqual(grp.sequences,["CTAAGCCT"])
        self.assertEqual(grp.counts,2894178)
        self.assertEqual(len(grp),1)

    def test_barcodegroup_add(self):
        """BarcodeGroup: check adding a sequence to the group
        """
        # Add sequences to a BarcodeGroup
        grp = BarcodeGroup("CTAAGCCT",2894178)
        grp.add("CTAAGCCA",92417)
        self.assertEqual(grp.reference,"CTAAGCCT")
        self.assertEqual(grp.sequences,["CTAAGCCT","CTAAGCCA"])
        self.assertEqual(grp.counts,2986595)
        self.assertEqual(len(grp),2)

    def test_barcodegroup_match(self):
        """BarcodeGroup: check matching sequences against the reference
        """
        # Match sequence against the group
        grp = BarcodeGroup("CTAAGCCT",2894178)
        # 1 mismatch allowed
        self.assertTrue(grp.match("CTAAGCCA",mismatches=1))
        self.assertFalse(grp.match("CGAAGCCA",mismatches=1))
        self.assertFalse(grp.match("CGATGCCA",mismatches=1))
        # Default (2 mismatches)
        self.assertTrue(grp.match("CTAAGCCA"))
        self.assertTrue(grp.match("CGAAGCCA"))
        self.assertFalse(grp.match("CGATGCCA"))

# SampleSheetBarcodes
class TestSampleSheetBarcodes(unittest.TestCase):
    def setUp(self):
        # Test data
        sample_sheet_header = "[Header]\nIEMFileVersion,4\n\n[Reads]\n150\n150\n\n[Settings]\nAdapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA\nAdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT\n\n"
        sample_sheet_single_index_with_lanes = """
"""
        # Temporary working dir
        self.wd = tempfile.mkdtemp(suffix='.test_SampleSheetBarcodes')
        # Create files
        sample_sheet_single_index_no_lanes = """[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description
ES1,,,,A006,GCCAAT,,
EP1,,,,A012,CTTGTA,,
ES2,,,,A005,ACAGTG,,
EP2,,,,A019,GTGAAA,,"""
        self.single_index_no_lanes = \
            os.path.join(self.wd,"single_index_no_lanes.csv")
        with open(self.single_index_no_lanes,"w") as fp:
            fp.write(sample_sheet_header +
                     sample_sheet_single_index_no_lanes)
        #
        sample_sheet_dual_index_no_lanes = """[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
SW1,SW1,,,N701,TAAGGCGA,S517,TCTTACGC,,
SW2,SW2,,,N702,CGTACTAG,S517,TCTTACGC,,
SW3,SW3,,,N703,AGGCAGAA,S517,TCTTACGC,,
SW4,SW4,,,N704,TCCTGAGC,S517,TCTTACGC,,
SW5,SW5,,,N705,GGACTCCT,S517,TCTTACGC,,
SW6,SW6,,,N706,TAGGCATG,S517,TCTTACGC,,"""
        self.dual_index_no_lanes = \
            os.path.join(self.wd,"dual_index_no_lanes.csv")
        with open(self.dual_index_no_lanes,"w") as fp:
            fp.write(sample_sheet_header +
                     sample_sheet_dual_index_no_lanes)
        #
        sample_sheet_single_index_with_lanes = """[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description
1,ES1,,,,A006,GCCAAT,,
1,EP1,,,,A012,CTTGTA,,
2,ES2,,,,A005,ACAGTG,,
2,EP2,,,,A019,GTGAAA,,"""
        self.single_index_with_lanes = \
            os.path.join(self.wd,"single_index_with_lanes.csv")
        with open(self.single_index_with_lanes,"w") as fp:
            fp.write(sample_sheet_header +
                     sample_sheet_single_index_with_lanes)
        #
        sample_sheet_dual_index_with_lanes = """[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,DI1,DI1,,,D701,CGTGTAGG,D501,GACCTGTA,HO,
1,DI2,DI2,,,D702,CGTGTAGG,D501,ATGTAACT,HO,
1,DI3,DI3,,,D703,CGTGTAGG,D501,GTTTCAGA,HO,
1,DI4,DI4,,,D704,CGTGTAGG,D501,CACAGGAT,HO,
2,E11,E11,,,D703,GCCAATAT,D503,TCTTTCCC,FL,
2,E12,E12,,,D704,CAGATCAT,D504,TCTTTCCC,FL,
3,AT1,AT1,,,D701,GGATTCGC,D501,TAGTAGCC,JF,
3,AT3,AT3,,,D702,TACCAGCG,D502,CGCTGCTG,JF,
3,AT4,AT4,,,D703,ACCGATTC,D503,GCGCGTAG,JF,
3,AT5,AT5,,,D704,CCGCGTAA,D504,GCAATAGA,JF,
3,AEx,AEx,,,D705,ATGCTCGT,D505,CTCGCATC,JF,
4,AEx,AEx,,,D705,ATGCTCGT,D505,CTCGCATC,JF,
4,AD5,AD5,,,D707,CCAGCAAT,D507,ATCGCGAG,JF,
4,D3K1,D3K1,,,D701,ACAGTGAT,D501,TCTTTCCC,FL,
4,D3K2,D3K2,,,D702,ATGTCAGA,D502,TCTTTCCC,FL,"""
        self.dual_index_with_lanes = \
            os.path.join(self.wd,"dual_index_with_lanes.csv")
        with open(self.dual_index_with_lanes,"w") as fp:
            fp.write(sample_sheet_header +
                     sample_sheet_dual_index_with_lanes)

    def tearDown(self):
        # Remove temporary working dir
        if os.path.isdir(self.wd):
            shutil.rmtree(self.wd)

    def test_single_index_no_lanes(self):
        """SampleSheetBarcodes: single index sample sheet with no lanes defined
        """
        s = SampleSheetBarcodes(self.single_index_no_lanes)
        # Check all barcodes
        self.assertEqual(s.barcodes(),["ACAGTG","CTTGTA",
                                       "GCCAAT","GTGAAA"])
        # Check all samples
        self.assertEqual(s.samples(),["EP1","EP2","ES1","ES2"])
        # Look up sample names for barcodes in each lane
        self.assertEqual(s.lookup_sample("GCCAAT"),"ES1")
        self.assertEqual(s.lookup_sample("CTTGTA"),"EP1")
        self.assertEqual(s.lookup_sample("ACAGTG"),"ES2")
        self.assertEqual(s.lookup_sample("GTGAAA"),"EP2")
        # Look up barcode matching sample names in each lane
        self.assertEqual(s.lookup_barcode("ES1"),"GCCAAT")
        self.assertEqual(s.lookup_barcode("EP1"),"CTTGTA")
        self.assertEqual(s.lookup_barcode("ES2"),"ACAGTG")
        self.assertEqual(s.lookup_barcode("EP2"),"GTGAAA")

    def test_dual_index_no_lanes(self):
        """SampleSheetBarcodes: dual index sample sheet with no lanes defined
        """
        s = SampleSheetBarcodes(self.dual_index_no_lanes)
        # Check all barcodes
        self.assertEqual(s.barcodes(),["AGGCAGAATCTTACGC",
                                       "CGTACTAGTCTTACGC",
                                       "GGACTCCTTCTTACGC",
                                       "TAAGGCGATCTTACGC",
                                       "TAGGCATGTCTTACGC",
                                       "TCCTGAGCTCTTACGC"])
        # Check all samples
        self.assertEqual(s.samples(),["SW1","SW2","SW3","SW4","SW5","SW6"])
        # Look up sample names for barcodes in each lane
        self.assertEqual(s.lookup_sample("TAAGGCGATCTTACGC"),"SW1")
        self.assertEqual(s.lookup_sample("CGTACTAGTCTTACGC"),"SW2")
        self.assertEqual(s.lookup_sample("AGGCAGAATCTTACGC"),"SW3")
        self.assertEqual(s.lookup_sample("TCCTGAGCTCTTACGC"),"SW4")
        self.assertEqual(s.lookup_sample("GGACTCCTTCTTACGC"),"SW5")
        self.assertEqual(s.lookup_sample("TAGGCATGTCTTACGC"),"SW6")
        # Look up barcode matching sample names in each lane
        self.assertEqual(s.lookup_barcode("SW1"),"TAAGGCGATCTTACGC")
        self.assertEqual(s.lookup_barcode("SW2"),"CGTACTAGTCTTACGC")
        self.assertEqual(s.lookup_barcode("SW3"),"AGGCAGAATCTTACGC")
        self.assertEqual(s.lookup_barcode("SW4"),"TCCTGAGCTCTTACGC")
        self.assertEqual(s.lookup_barcode("SW5"),"GGACTCCTTCTTACGC")
        self.assertEqual(s.lookup_barcode("SW6"),"TAGGCATGTCTTACGC")

    def test_single_index_with_lanes(self):
        """SampleSheetBarcodes: single index sample sheet with lanes defined
        """
        s = SampleSheetBarcodes(self.single_index_with_lanes)
        # Check all barcodes
        self.assertEqual(s.barcodes(),["ACAGTG","CTTGTA",
                                       "GCCAAT","GTGAAA"])
        # Check barcodes in each lane
        self.assertEqual(s.barcodes(1),["CTTGTA","GCCAAT"])
        self.assertEqual(s.barcodes(2),["ACAGTG","GTGAAA"])
        # Check all samples
        self.assertEqual(s.samples(),["EP1","EP2","ES1","ES2"])
        # Check samples in each lane
        self.assertEqual(s.samples(1),["EP1","ES1"])
        self.assertEqual(s.samples(2),["EP2","ES2"])
        # Look up sample names for barcodes in each lane
        self.assertEqual(s.lookup_sample("GCCAAT",1),"ES1")
        self.assertEqual(s.lookup_sample("CTTGTA",1),"EP1")
        self.assertEqual(s.lookup_sample("ACAGTG",2),"ES2")
        self.assertEqual(s.lookup_sample("GTGAAA",2),"EP2")
        # Look up barcode matching sample names in each lane
        self.assertEqual(s.lookup_barcode("ES1",1),"GCCAAT")
        self.assertEqual(s.lookup_barcode("EP1",1),"CTTGTA")
        self.assertEqual(s.lookup_barcode("ES2",2),"ACAGTG")
        self.assertEqual(s.lookup_barcode("EP2",2),"GTGAAA")

    def test_dual_index_with_lanes(self):
        """SampleSheetBarcodes: dual index sample sheet with lanes defined
        """
        s = SampleSheetBarcodes(self.dual_index_with_lanes)
        # Check all barcodes
        self.assertEqual(s.barcodes(),["ACAGTGATTCTTTCCC",
                                       "ACCGATTCGCGCGTAG",
                                       "ATGCTCGTCTCGCATC",
                                       "ATGCTCGTCTCGCATC",
                                       "ATGTCAGATCTTTCCC",
                                       "CAGATCATTCTTTCCC",
                                       "CCAGCAATATCGCGAG",
                                       "CCGCGTAAGCAATAGA",
                                       "CGTGTAGGATGTAACT",
                                       "CGTGTAGGCACAGGAT",
                                       "CGTGTAGGGACCTGTA",
                                       "CGTGTAGGGTTTCAGA",
                                       "GCCAATATTCTTTCCC",
                                       "GGATTCGCTAGTAGCC",
                                       "TACCAGCGCGCTGCTG"])
        # Check barcodes in each lane
        self.assertEqual(s.barcodes(1),["CGTGTAGGATGTAACT",
                                        "CGTGTAGGCACAGGAT",
                                        "CGTGTAGGGACCTGTA",
                                        "CGTGTAGGGTTTCAGA"])
        self.assertEqual(s.barcodes(2),["CAGATCATTCTTTCCC",
                                        "GCCAATATTCTTTCCC"])
        self.assertEqual(s.barcodes(3),["ACCGATTCGCGCGTAG",
                                        "ATGCTCGTCTCGCATC",
                                        "CCGCGTAAGCAATAGA",
                                        "GGATTCGCTAGTAGCC",
                                        "TACCAGCGCGCTGCTG"])
        self.assertEqual(s.barcodes(4),["ACAGTGATTCTTTCCC",
                                        "ATGCTCGTCTCGCATC",
                                        "ATGTCAGATCTTTCCC",
                                        "CCAGCAATATCGCGAG"])
        # Check all samples
        self.assertEqual(s.samples(),["AD5","AEx","AEx",
                                      "AT1","AT3","AT4","AT5",
                                      "D3K1","D3K2",
                                      "DI1","DI2","DI3","DI4",
                                      "E11","E12"])
        # Check samples in each lane
        self.assertEqual(s.samples(1),["DI1","DI2","DI3","DI4"])
        self.assertEqual(s.samples(2),["E11","E12"])
        self.assertEqual(s.samples(3),["AEx","AT1","AT3","AT4","AT5"])
        self.assertEqual(s.samples(4),["AD5","AEx","D3K1","D3K2"])
        # Look up sample names for barcodes in each lane
        self.assertEqual(s.lookup_sample("CGTGTAGGGACCTGTA",1),"DI1")
        self.assertEqual(s.lookup_sample("CGTGTAGGATGTAACT",1),"DI2")
        self.assertEqual(s.lookup_sample("CGTGTAGGGTTTCAGA",1),"DI3")
        self.assertEqual(s.lookup_sample("CGTGTAGGCACAGGAT",1),"DI4")
        self.assertEqual(s.lookup_sample("GCCAATATTCTTTCCC",2),"E11")
        self.assertEqual(s.lookup_sample("CAGATCATTCTTTCCC",2),"E12")
        self.assertEqual(s.lookup_sample("GGATTCGCTAGTAGCC",3),"AT1")
        self.assertEqual(s.lookup_sample("TACCAGCGCGCTGCTG",3),"AT3")
        self.assertEqual(s.lookup_sample("ACCGATTCGCGCGTAG",3),"AT4")
        self.assertEqual(s.lookup_sample("CCGCGTAAGCAATAGA",3),"AT5")
        self.assertEqual(s.lookup_sample("ATGCTCGTCTCGCATC",3),"AEx")
        self.assertEqual(s.lookup_sample("ATGCTCGTCTCGCATC",4),"AEx")
        self.assertEqual(s.lookup_sample("CCAGCAATATCGCGAG",4),"AD5")
        self.assertEqual(s.lookup_sample("ACAGTGATTCTTTCCC",4),"D3K1")
        self.assertEqual(s.lookup_sample("ATGTCAGATCTTTCCC",4),"D3K2")
        # Look up barcode matching sample names in each lane
        self.assertEqual(s.lookup_barcode("DI1",1),"CGTGTAGGGACCTGTA")
        self.assertEqual(s.lookup_barcode("DI2",1),"CGTGTAGGATGTAACT")
        self.assertEqual(s.lookup_barcode("DI3",1),"CGTGTAGGGTTTCAGA")
        self.assertEqual(s.lookup_barcode("DI4",1),"CGTGTAGGCACAGGAT")
        self.assertEqual(s.lookup_barcode("E11",2),"GCCAATATTCTTTCCC")
        self.assertEqual(s.lookup_barcode("E12",2),"CAGATCATTCTTTCCC")
        self.assertEqual(s.lookup_barcode("AT1",3),"GGATTCGCTAGTAGCC")
        self.assertEqual(s.lookup_barcode("AT3",3),"TACCAGCGCGCTGCTG")
        self.assertEqual(s.lookup_barcode("AT4",3),"ACCGATTCGCGCGTAG")
        self.assertEqual(s.lookup_barcode("AT5",3),"CCGCGTAAGCAATAGA")
        self.assertEqual(s.lookup_barcode("AEx",3),"ATGCTCGTCTCGCATC")
        self.assertEqual(s.lookup_barcode("AEx",4),"ATGCTCGTCTCGCATC")
        self.assertEqual(s.lookup_barcode("AD5",4),"CCAGCAATATCGCGAG")
        self.assertEqual(s.lookup_barcode("D3K1",4),"ACAGTGATTCTTTCCC")
        self.assertEqual(s.lookup_barcode("D3K2",4),"ATGTCAGATCTTTCCC")

    def test_non_existent_lane_raises_exception(self):
        """SampleSheetBarcodes: request non-existent lane raises KeyError
        """
        # Bad lane for sample sheet with lanes
        s = SampleSheetBarcodes(self.dual_index_with_lanes)
        self.assertRaises(KeyError,s.barcodes,5)
        self.assertRaises(KeyError,s.samples,5)
        # Any lane for sample sheet with no lanes
        s = SampleSheetBarcodes(self.dual_index_no_lanes)
        self.assertRaises(KeyError,s.barcodes,1)
        self.assertRaises(KeyError,s.samples,1)

# Reporter
class TestReporter(unittest.TestCase):
    def setUp(self):
        # Temporary working dir (if needed)
        self.wd = None

    def _make_working_dir(self):
        if self.wd is None:
            self.wd = tempfile.mkdtemp(suffix='.test_Reporter')

    def tearDown(self):
        # Remove temporary working dir
        if self.wd is not None and os.path.isdir(self.wd):
            shutil.rmtree(self.wd)

    def test_add(self):
        """Reporter: can add content
        """
        reporter = Reporter()
        self.assertEqual(len(reporter),0)
        self.assertFalse(reporter)
        self.assertEqual(str(reporter),"")
        reporter.add("Title text",title=True)
        reporter.add("Some words")
        self.assertEqual(len(reporter),2)
        self.assertTrue(reporter)
        self.assertEqual(str(reporter),
                         """Title text
==========
Some words""")

    def test_write(self):
        """Reporter: can write to a text file
        """
        self._make_working_dir()
        reporter = Reporter()
        reporter.add("Test Document",title=True)
        reporter.add("Lorem ipsum")
        report_file = os.path.join(self.wd,"report.txt")
        reporter.write(filen=report_file)
        self.assertTrue(os.path.isfile(report_file))
        expected_contents = """Test Document
=============
Lorem ipsum
"""
        self.assertTrue(os.path.exists(report_file))
        self.assertEqual(open(report_file,'r').read(),
                         expected_contents)

    def test_write_xls(self):
        """Reporter: can write to an XLS file
        """
        self._make_working_dir()
        reporter = Reporter()
        reporter.add("Test Document",title=True)
        reporter.add("Lorem ipsum")
        report_xls = os.path.join(self.wd,"report.xls")
        reporter.write_xls(report_xls)
        self.assertTrue(os.path.isfile(report_xls))

# report_barcodes
class TestReportBarcodesFunction(unittest.TestCase):
    def setUp(self):
        # Temporary working dir (if needed)
        self.wd = None

    def _make_working_dir(self):
        if self.wd is None:
            self.wd = tempfile.mkdtemp(suffix='.test_report_barcodes')

    def tearDown(self):
        # Remove temporary working dir
        if self.wd is not None and os.path.isdir(self.wd):
            shutil.rmtree(self.wd)

    def _make_file(self,name,contents):
        # Create a file under the working directory
        # called "name" and populated with "contents"
        # Working directory will be created if not already
        # set up
        # Returns the path to the file
        self._make_working_dir()
        filen = os.path.join(self.wd,name)
        with open(filen,'w') as fp:
            fp.write(contents)
        return filen

    def test_report_barcodes(self):
        """report_barcodes: check output for mismatches and sample sheet
        """
        # Create sample sheet
        sample_sheet_file = self._make_file("SampleSheet.csv",
                                            """[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description
1,SMPL1,,,,A006,CATGCGCGGTA,,
1,SMPL2,,,,A012,GCTGCGCGGTC,,
""")
        # Set up barcode counts
        bc = BarcodeCounter()
        bc.count_barcode("TATGCGCGGTA",lane=1,incr=285302)
        bc.count_barcode("CATGCGCGGTA",lane=1,incr=8532)
        bc.count_barcode("GATGCGCGGTA",lane=1,incr=5321)
        bc.count_barcode("GCTGCGCGGTA",lane=1,incr=7853)
        bc.count_barcode("GCTGCGCGGTC",lane=1,incr=325394)
        analysis = bc.analyse(lane=1,mismatches=2,
                              sample_sheet=sample_sheet_file)
        ##"CATGCGCGGTA","TATGCGCGGTA","GATGCGCGGTA","GCTGCGCGGTA" = 307008
        ##"GCTGCGCGGTC" = 325394
        self.assertEqual(analysis.cutoff,None)
        self.assertEqual(analysis.mismatches,2)
        self.assertEqual(analysis.total_reads,632402)
        self.assertEqual(analysis.coverage,632402)
        self.assertEqual(analysis.barcodes,["GCTGCGCGGTC",
                                            "CATGCGCGGTA"])
        self.assertEqual(analysis.counts["GCTGCGCGGTC"].reads,325394)
        self.assertEqual(analysis.counts["CATGCGCGGTA"].reads,307008)
        self.assertEqual(analysis.counts["GCTGCGCGGTC"].sample,"SMPL2")
        self.assertEqual(analysis.counts["CATGCGCGGTA"].sample,"SMPL1")
        self.assertEqual(analysis.counts["GCTGCGCGGTC"].sequences,1)
        self.assertEqual(analysis.counts["CATGCGCGGTA"].sequences,4)
        # Create report
        reporter = report_barcodes(bc,
                                   lane=1,
                                   mismatches=2,
                                   sample_sheet=sample_sheet_file)
        # Check content
        self.assertEqual(str(reporter),
                         """Barcode analysis for lane #1
============================
Barcodes have been grouped by allowing 2 mismatches

#Rank	Index	Sample	N_seqs	N_reads	%reads	(%Total_reads)
    1	GCTGCGCGGTC	SMPL2	1	325394	51.5%	(51.5%)
    2	CATGCGCGGTA	SMPL1	4	307008	48.5%	(100.0%)""")

# Main program
if __name__ == '__main__':
    p = optparse.OptionParser(usage=
                              "\n\t%prog FASTQ [FASTQ...]\n"
                              "\t%prog DIR\n"
                              "\t%prog -c COUNTS_FILE [COUNTS_FILE...]",
                              version="%%prog %s" % __version__,
                              description="Collate and report counts and "
                              "statistics for Fastq index sequences (aka "
                              "barcodes). If multiple Fastq files are "
                              "supplied then sequences will be pooled before "
                              "being analysed. If a single directory is "
                              "supplied then this will be assumed to "
                              "be an output directory from CASAVA or "
                              "bclToFastq and files will be processed on a "
                              "per-lane basis. If the -c option is "
                              "supplied then the input must be one or more "
                              "file of barcode counts generated previously "
                              "using the -o option.")
    p.add_option('-c','--counts',
                 action='store_true',dest='use_counts',default=False,
                 help="input is one or more counts files generated by "
                 "previous runs")
    p.add_option('-o','--output',
                 action='store',dest='counts_file_out',default=None,
                 help="output all counts to tab-delimited file "
                 "COUNTS_FILE. This can be used again in another "
                 "run by specifying the '-c' option.")
    p.add_option('-l','--lanes',action='store',dest='lanes',default=None,
                 help="restrict analysis to the specified lane numbers "
                 "(default is to process all lanes). Multiple lanes "
                 "can be specified using ranges (e.g. '2-3'), comma-"
                 "separated list ('5,7') or a mixture ('2-3,5,7')")
    p.add_option('-m','--mismatches',action='store',dest='mismatches',
                 default=0,type='int',
                 help="maximum number of mismatches to use when "
                 "grouping similar barcodes (default is 0)")
    p.add_option('--cutoff',action='store',dest='cutoff',
                 default=0.001,type='float',
                 help="exclude barcodes with a smaller fraction of "
                 "associated reads than CUTOFF, e.g. '0.01' excludes "
                 "barcodes with < 0.01% of reads (default is 0.1%)")
    p.add_option('--coverage',action='store',dest='coverage',
                 default=None,type='float',
                 help="include most numerous barcodes to cover only "
                 "the fraction of reads specified by COVERAGE, e.g. "
                 "'0.9' limits barcodes to those associated with 90% "
                 " of total reads (default is to include all barcodes)")
    p.add_option('-s','--sample-sheet',
                 action='store',dest='sample_sheet',default=None,
                 help="report best matches against barcodes in "
                 "SAMPLE_SHEET")
    p.add_option('-r','--report',
                 action='store',dest='report_file',default=None,
                 help="write report to REPORT_FILE (otherwise write to "
                 "stdout)")
    p.add_option('-x','--xls',
                 action='store',dest='xls_file',default=None,
                 help="write XLS version of report to XLS_FILE")
    p.add_option('-n','--no-report',
                 action='store_true',dest='no_report',default=None,
                 help="suppress reporting (overrides --report)")
    # Report name and version
    p.print_version()
    # Process command line
    opts,args = p.parse_args()
    # Determine subset of lanes to examine
    lanes = parse_lanes_expression(opts.lanes)
    # Determine mode
    if opts.use_counts:
        # Read counts from counts file(s)
        counts = BarcodeCounter(*args)
    elif len(args) == 1 and os.path.isdir(args[0]):
        # Generate counts from bcl2fastq output
        counts = count_barcodes_bcl2fastq(args[0])
    else:
        # Generate counts from fastq files
        counts = count_barcodes(args)
    # Deal with cutoff
    if opts.cutoff == 0.0:
        cutoff = None
    else:
        cutoff = opts.cutoff
    # Report the counts
    if not opts.no_report:
        reporter = Reporter()
        if lanes is None:
            report_barcodes(counts,
                            cutoff=cutoff,
                            sample_sheet=opts.sample_sheet,
                            mismatches=opts.mismatches,
                            reporter=reporter)
        else:
            for lane in lanes:
                if lane not in counts.lanes:
                    logging.error("Requested analysis for lane %d but "
                                  "only have counts for lanes %s" %
                                  (lane,
                                   ','.join([str(l) for l in counts.lanes])))
                    sys.exit(1)
                report_barcodes(counts,
                                lane=lane,
                                cutoff=cutoff,
                                sample_sheet=opts.sample_sheet,
                                mismatches=opts.mismatches,
                                reporter=reporter)
        if opts.report_file is not None:
            print "Writing report to %s" % opts.report_file
            reporter.write(filen=opts.report_file)
        else:
            reporter.write()
        if opts.xls_file is not None:
            print "Writing XLS to %s" % opts.xls_file
            reporter.write_xls(opts.xls_file)
    # Output counts if requested
    if opts.counts_file_out is not None:
        counts.write(opts.counts_file_out)
