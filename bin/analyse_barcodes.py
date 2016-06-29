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
from bcftbx.FASTQFile import FastqIterator

__version__ = "0.0.1"

class BarcodeCounter:
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
        return sorted(self._seqs.keys())

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
        if lane is None:
            return self._seqs_all[barcode]
        else:
            return self._seqs[lane][barcode]

    def counts_all(self,barcode):
        """
        Number of counts for the barcode across all lanes

        """
        return self._seqs_all[barcode]

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
            return sum([self._seqs[lane][barcode]
                        for barcode in self._seqs[lane]])

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
              exclude_reads=0.000001):
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

        """
        # Initialise
        barcodes = self.filter_barcodes(lane=lane,cutoff=exclude_reads)
        nreads = self.nreads(lane=lane)
        groups = []
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

        """
        if not mismatches:
            groups = None
            barcodes = self.filter_barcodes(cutoff=cutoff,lane=lane)
        else:
            groups = self.group(lane,mismatches=mismatches,
                                cutoff=cutoff)
            barcodes = [grp.reference for grp in groups]
        analysis = {
            'barcodes': barcodes,
            'cutoff': cutoff,
            'counts': dict(),
            'total_reads': self.nreads(lane=lane),
            'mismatches': mismatches,
        }
        sample_lookup = {}
        if sample_sheet is not None:
            sample_sheet = SampleSheet(sample_sheet)
            sample_id = sample_sheet.sample_id_column
            for line in sample_sheet.data:
                sample = line[sample_id]
                index_seq = normalise_barcode(samplesheet_index_sequence(line))
                sample_lookup[index_seq] = sample
        cum_reads = 0
        if groups:
            for group in groups:
                barcode = group.reference
                barcode_reads = group.counts
                cum_reads += barcode_reads
                try:
                    sample = sample_lookup[barcode]
                except KeyError:
                    sample = None
                analysis['counts'][barcode] = { 'reads': barcode_reads,
                                                'sample': sample,
                                                'sequences': len(group) }
        else:
            for barcode in barcodes:
                barcode_reads = self.counts(barcode,lane)
                cum_reads += barcode_reads
                try:
                    sample = sample_lookup[barcode]
                except KeyError:
                    sample = None
                analysis['counts'][barcode] = { 'reads': barcode_reads,
                                                'sample': sample,
                                                'sequences': 1 }
        analysis['coverage'] = cum_reads
        return analysis

class BarcodeGroup:
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

def normalise_barcode(seq):
    """
    Return normalised version of barcode sequence

    This standardises the sequence so that:

    - all bases are uppercase
    - dual index barcodes have '-' and '+' removed

    """
    return str(seq).upper().replace('-','').replace('+','')

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

def print_title(text,underline):
    print "\n%s\n%s" % (text,underline*len(text))

def write_title(fp,text,underline):
    fp.write("\n%s\n%s\n" % (text,underline*len(text)))

def report_barcodes(counts,lane=None,sample_sheet=None,cutoff=None,
                    mismatches=0,fp=None):
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
      fp (File): optional, file object opened for writing
        (defaults to stdout)

    """
    # Where to write the output
    if fp is None:
        fp = sys.stdout
    # Check lanes
    if lane is not None:
        write_title(fp,"Barcode analysis for lane #%d" % lane,'=')
    else:
        write_title(fp,"Barcode analysis for all lanes",'=')
    # Get analysis
    analysis = counts.analyse(lane=lane,cutoff=cutoff,
                              sample_sheet=sample_sheet,
                              mismatches=mismatches)
    # Report settings
    if cutoff is not None:
        fp.write("Barcodes which cover less than %.1f%% of reads have been "
                 "excluded\n" % (cutoff*100.0))
    fp.write("Reported barcodes cover %.1f%% of the data\n" %
             (float(analysis['coverage'])/float(analysis['total_reads'])*100.0))
    if mismatches:
        fp.write("Barcodes have been grouped by allowing %d mismatch%s\n" %
                 (mismatches,
                  ('' if mismatches == 1 else 'es')))
    cumulative_reads = 0
    fp.write("\n%s\n" % '\t'.join(("#Rank",
                                   "Index",
                                   "Sample",
                                   "N_seqs",
                                   "N_reads",
                                   "%reads",
                                   "(%Total_reads)")))
    for i,barcode in enumerate(analysis['barcodes']):
        cumulative_reads += analysis['counts'][barcode]['reads']
        sample_name = analysis['counts'][barcode]['sample']
        if sample_name is None:
            sample_name = ''
        fp.write("%s\n" % '\t'.join([str(x) for x in
                                     ('% 5d' % (i+1),
                                      barcode,
                                      sample_name,
                                      analysis['counts'][barcode]['sequences'],
                                      analysis['counts'][barcode]['reads'],
                                      '%.1f%%' % (float(analysis['counts'][barcode]['reads'])/float(analysis['total_reads'])*100.0),
                                      '(%.1f%%)' % (float(cumulative_reads)/float(analysis['total_reads'])*100.0))]))

def samplesheet_index_sequence(line):
    """
    Return the index sequence for a sample sheet line

    Arguments:
      line (TabDataLine): line from a SampleSheet instance

    """
    # Index sequence
    try:
        # Try dual-indexed IEM4 format
        return "%s-%s" % (line['index'].strip(),
                          line['index2'].strip())
    except KeyError:
        pass
    # Try single indexed IEM4 (no index2)
    try:
        return line['index'].strip()
    except KeyError:
        pass
    # Try CASAVA format
    indx = line['Index'].strip()
    if not indx:
        indx = None
    return indx

def match_against_samplesheet(groups,sample_sheet,lanes,mismatches,nreads,
                              fp=None):
    """
    Find best matches to barcodes in samplesheet

    """
    # Where to write the output
    if fp is None:
        fp = sys.stdout
    write_title(fp,"Best matches of sample sheet against groups",'-')
    # Get barcodes from sample sheet
    sample_sheet = SampleSheet(sample_sheet)
    sample_id = sample_sheet.sample_id_column
    sample = {}
    index_sequences = []
    if sample_sheet.has_lanes:
        for line in sample_sheet.data:
            if line['Lane'] in lanes:
                samples[samplesheet_index_sequence(line)] = line[sample_id]
                index_sequences.append(samplesheet_index_sequence(line))
    else:
        for line in sample_sheet.data:
            samples[samplesheet_index_sequence(line)] = line[sample_id]
            index_sequences.append(samplesheet_index_sequence(line))
    # Normalise sequences
    index_sequences = [normalise_barcode(s) for s in index_sequences]
    normalised_samples = {}
    for sample in samples:
        normalised_samples[normalise_barcode(sample)] = samples[sample]
    samples = normalised_samples
    # Match barcodes
    fp.write("%d mismatch%s allowed\n" % (mismatches,
                                          ('' if mismatches == 1
                                           else 'es')))
    fp.write("%12s\t%18s\t%18s\t%8s\t%6s\n" % ("Sample",
                                               "Sample_index",
                                               "Group_seq",
                                               "Nreads",
                                               "%reads"))
    for seq in index_sequences:
        matched = False
        rejected = []
        for i,group in enumerate(groups):
            if not matched and group.match(seq,mismatches):
                group_counts = group.counts
                fp.write("%12s\t% 18s\t% 18s\t% 8d\t% 5.1f%%\n" %
                         (samples[seq],seq,
                          group.reference,
                          group_counts,
                          float(group_counts)/float(nreads)*100.0))
                matched = True
                rejected.extend(groups[i+1:])
                break
            else:
                rejected.append(group)
        if not matched:
            fp.write("%12s\t%18s\t%18s\t%8s\t%6s\n" % (samples[seq],
                                                         seq,
                                                         "No match",
                                                         "-",
                                                         "-"))
        groups = rejected

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
                 default=None,type='float',
                 help="exclude barcodes with a smaller fraction of "
                 "associated reads than CUTOFF, e.g. '0.001' excludes "
                 "barcodes with < 0.1% of reads (default is to include "
                 "all barcodes)")
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
    # Report the counts
    if not opts.no_report:
        if opts.report_file is not None:
            print "Writing report to %s" % opts.report_file
            fp = open(opts.report_file,'w')
        else:
            fp = sys.stdout
        if lanes is None:
            report_barcodes(counts,
                            cutoff=opts.cutoff,
                            sample_sheet=opts.sample_sheet,
                            mismatches=opts.mismatches,
                            fp=fp)
        else:
            for lane in lanes:
                report_barcodes(counts,
                                lane=lane,
                                cutoff=opts.cutoff,
                                sample_sheet=opts.sample_sheet,
                                mismatches=opts.mismatches,
                                fp=fp)
    # Output counts if requested
    if opts.counts_file_out is not None:
        counts.write(opts.counts_file_out)
