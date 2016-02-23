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

    Grouping barcodes
    -----------------

    To put barcodes into groups use the 'group' method,
    e.g. to group all barcodes from lane 2:

    >>> groups = bc.group(2)

    which produces a list of BarcodeGroup instances.

    """

    def __init__(self,counts_file=None):
        """
        Create a new BarcodeCounter instance

        """
        self._seqs = {}
        if counts_file is not None:
            self.read(counts_file)

    def count_barcode(self,barcode,lane=None,incr=1):
        """
        Increment count of a barcode sequence

        Arguments:
          lane (int): lane that the barcode appears
            in (None if unknown)
          barcode (str): barcode sequence to count
          incr (int): increment the count for the
            barcode in the lane by this amount
            (defaults to 1)

        """
        try:
            self._seqs[lane][barcode] += incr
        except KeyError:
            try:
                self._seqs[lane][barcode] = incr
            except KeyError:
                self._seqs[lane] = { barcode: incr }

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
            bc = []
            for i in self._seqs.keys():
                bc.extend([s for s in self._seqs[i]])
            return sorted(bc,cmp=lambda x,y: cmp(self.counts(y),
                                                 self.counts(x)))
        else:
            bc = [s for s in self._seqs[lane]]
            return sorted(bc,cmp=lambda x,y: cmp(self.counts(y,lane),
                                                 self.counts(x,lane)))

    def counts(self,barcode,lane):
        """
        Return the number of counts for the barcode

        """
        return self._seqs[lane][barcode]

    def counts_all(self,barcode):
        """
        Number of counts for the barcode across all lanes

        """
        return sum(map(lambda i: 0 if barcode not in self._seqs[i]
                       else self._seqs[i][barcode],self.lanes))

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
            return sum([sum([self._seqs[i][b] for b in self._seqs[i]])
                        for i in self._seqs.keys()])
        else:
            return sum([self._seqs[lane][barcode]
                        for barcode in self._seqs[lane]])

    def read(self,filen):
        """
        Read count data from a file

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

    def group(self,lane,mismatches=2,n=None):
        """
        Put barcodes into groups of similar sequences

        Returns a list of BarcodeGroup instances

        """
        # Initialise
        barcodes = self.barcodes(lane)
        if n is not None:
            barcodes = barcodes[:n]
        groups = []
        # Iteratively assign barcodes to groups
        # until we run out
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
            groups.append(group)
            # Reset sequences to check in next iteration
            barcodes = rejected
        # Sort groups into order of total counts
        groups = sorted(groups,cmp=lambda x,y: cmp(y.counts,x.counts))
        return groups

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

def report(counts,lane,n,mismatches,sample_sheet=None,fp=None):
    """
    Report barcodes

    """
    # Where to write the output
    if fp is None:
        fp = sys.stdout
    write_title(fp,"Barcode analysis for lane #%d" % lane,'=')
    # Calculate total number of reads
    nreads = counts.nreads(lane)
    write_title(fp,"Top barcode sequences in lane %d" % lane,'-')
    for i,barcode in enumerate(counts.barcodes(lane)[:n]):
        fp.write("% 8d\t%s\t% 8d\t% 5.1f%%\n" %
                 (i+1,
                  barcode,
                  counts.counts(barcode,lane),
                  float(counts.counts(barcode,lane))/nreads*100.0))
    # Group matching barcodes
    mismatches = 2
    write_title(fp,"Groups",'-')
    groups = counts.group(lane,mismatches)
    fp.write("%d group%s (%d mismatch%s allowed)\n" %
             (len(groups),
              ('' if len(groups) == 1
               else 's'),
              mismatches,
              ('' if mismatches == 1
               else 'es')))
    for i,group in enumerate(groups[:n]):
        ngroup = group.counts
        fp.write("% 6d\t%s\t% 6d\t% 8d\t% 5.1f%%\n" %
                 (i+1,
                  group.reference,
                  len(group),
                  ngroup,
                  float(ngroup)/nreads*100.0))
    # Match against sample sheet
    if sample_sheet is not None:
        match_against_samplesheet(groups,sample_sheet,lanes,mismatches,
                                  fp=fp)

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

def match_against_samplesheet(groups,sample_sheet,lanes,mismatches,fp=None):
    """
    Find best matches to barcodes in samplesheet

    """
    # Where to write the output
    if fp is None:
        fp = sys.stdout
    write_title(fp,"Best matches against sample sheet",'-')
    # Get barcodes from sample sheet
    sample_sheet = SampleSheet(sample_sheet)
    sample_id = sample_sheet.sample_id_column
    if sample_sheet.has_lanes:
        samples = { samplesheet_index_sequence(line): line[sample_id]
                    for line in sample_sheet.data
                    if line['Lane'] in lanes }
        index_sequences = [samplesheet_index_sequence(line)
                           for line in sample_sheet.data
                           if line['Lane'] in lanes]
    else:
        samples = { samplesheet_index_sequence(line): line[sample_id]
                    for line in sample_sheet.data }
        index_sequences = [samplesheet_index_sequence(line)
                           for line in sample_sheet.data]
    # Normalise sequences
    index_sequences = [s.replace('-','') for s in index_sequences]
    samples = { s.replace('-',''): samples[s] for s in samples }
    # Get groups
    ##groups = counts.group(lane,mismatches)
    # Match barcodes
    fp.write("%d mismatch%s allowed\n" % (mismatches,
                                          ('' if mismatches == 1
                                           else 'es')))
    for seq in index_sequences:
        matched = False
        rejected = []
        for i,group in enumerate(groups):
            if not matched and group.match(seq):
                fp.write("%s\t%s\t%s\t%d\n" % (samples[seq],seq,
                                               group.reference,group.counts))
                matched = True
                rejected.extend(groups[i+1:])
                break
            else:
                rejected.append(group)
        if not matched:
            fp.write("%s\t%s\tNo match\t\n" % (samples[seq],seq))
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
                              "\t%prog -c COUNTS_FILE",
                              version="%%prog %s" % __version__)
    p.add_option('-c','--counts',
                 action='store',dest='counts_file_in',default=None,
                 help="counts file generated by a previous run")
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
                 default=1,type='int',
                 help="maximum number of mismatches to use when "
                 "grouping similar barcodes (default is 1)")
    p.add_option('-s','--sample-sheet',
                 action='store',dest='sample_sheet',default=None,
                 help="report best matches against barcodes in "
                 "SAMPLE_SHEET")
    p.add_option('-r','--report',
                 action='store',dest='report_file',default=None,
                 help="write report to REPORT_FILE (otherwise write to "
                 "stdout)")
    # Report name and version
    p.print_version()
    # Process command line
    opts,args = p.parse_args()
    # Determine subset of lanes to examine
    lanes = parse_lanes_expression(opts.lanes)
    # Determine mode
    if opts.counts_file_in is not None:
        # Read counts from counts file
        counts = BarcodeCounter(opts.counts_file_in)
    elif len(args) == 1 and os.path.isdir(args[0]):
        # Generate counts from bcl2fastq output
        counts = count_barcodes_bcl2fastq(args[0])
    else:
        # Generate counts from fastq files
        counts = count_barcodes(args)
    # Report the counts
    if opts.report_file is not None:
        print "Writing report to %s" % opts.report_file
        fp = open(opts.report_file,'w')
    else:
        fp = sys.stdout
    if lanes is None:
        lanes = counts.lanes
    for lane in lanes:
        report(counts,lane,20,opts.mismatches,opts.sample_sheet,fp=fp)
    # Output counts if requested
    if opts.counts_file_out is not None:
        counts.write(opts.counts_file_out)

        
        
