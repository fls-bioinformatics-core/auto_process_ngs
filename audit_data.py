#!/bin/env python
#
#     audit_data.py: audit sequencing data and disk usage
#     Copyright (C) University of Manchester 2014 Peter Briggs
#
########################################################################
#
# audit_data.py
#
#########################################################################


"""Perform auditing operations for sequencing data directories

"""

#######################################################################
# Module metadata
#######################################################################

__version__ = '0.0.10'

#######################################################################
# Import modules that this module depends on
#######################################################################

import platforms
import applications
import auto_process_utils
import bcf_utils
import logging
import optparse
import sys
import os
import time
import data_manager
import simple_xls

#######################################################################
# Classes
#######################################################################

class SeqDataSizes:
    def __init__(self,name,dirn,year=None,platform=None,is_subdir=False,
                 include_subdirs=True):
        """Create a new SeqDataDirectorySizes object

        name: name for the object
        dirn: path to the directory

        """
        self.name = name
        self.dirn = dirn
        self.year = year
        self.platform = platform
        self.du_total = None
        self.du_fastqs = None
        self.du_fastqgzs = None
        self.du_solid = None
        self.du_bam = None
        self.du_compressed_bams = None
        self.du_sams = None
        self.du_compressed_sams = None
        self.du_beds = None
        self.is_subdir = is_subdir
        self.include_subdirs = include_subdirs
        self.subdirs = []
        if not is_subdir and include_subdirs:
            for d in auto_process_utils.list_dirs(self.dirn):
                self.subdirs.append(SeqDataSizes(d,os.path.join(self.dirn,d),
                                                 year=self.year,
                                                 platform=self.platform,
                                                 is_subdir=True))

    @property
    def parent_dir(self):
        """Return the 'parent' directory

        For non-subdir instances, the parent directory is the directory
        above the run directory, with the platform and year stripped off
        (if provided).

        For subdir instances, the parent directory is the basename of the
        parent run directory.

        """
        parent = os.path.dirname(self.dirn)
        if self.is_subdir:
            parent = os.path.basename(parent)
        else:
            if self.platform is not None and parent.endswith(self.platform):
                parent = parent[:-len(self.platform)].rstrip(os.sep)
                if self.year is not None and parent.endswith(str(year)):
                    parent = parent[:-len(str(year))].rstrip(os.sep)
        return parent

    @property
    def du_other(self):
        try:
            return (self.du_total -
                    self.du_fastqs -
                    self.du_fastqgzs - 
                    self.du_solid -
                    self.du_bams -
                    self.du_compressed_bams -
                    self.du_sams -
                    self.du_compressed_sams -
                    self.du_beds)
        except TypeError:
            return None

    @property
    def users(self):
        """
        """
        return data_manager.DataDir(self.dirn).users

    def get_disk_usage(self):
        d = data_manager.DataDir(self.dirn)
        self.du_total = d.get_size()
        self.du_fastqs = d.get_size(pattern='.*\.fastq$')
        self.du_fastqgzs = d.get_size(pattern='.*\.fastq.gz$')
        self.du_solid = d.get_size(pattern='.*\.(csfasta|qual)$')
        self.du_bams = d.get_size(pattern='.*\.bam$')
        self.du_sams = d.get_size(pattern='.*\.sam$')
        self.du_compressed_bams = d.get_size(pattern='.*\.bam\.(gz|bz2)$')
        self.du_compressed_sams = d.get_size(pattern='.*\.sam\.(gz|bz2)$')
        self.du_beds = d.get_size(pattern='.*\.bed$')
        for subdir in self.subdirs:
            subdir.get_disk_usage()

    def dataline(self,formatter):
        if not self.is_subdir:
            dataline = [self.year,
                        self.platform,
                        self.parent_dir,
                        self.name]
        else:
            dataline = ['','','',self.name]
        dataline.extend([formatter.format(self.du_solid),
                         formatter.format(self.du_fastqgzs),
                         formatter.format(self.du_fastqs),
                         formatter.format(self.du_compressed_bams),
                         formatter.format(self.du_bams),
                         formatter.format(self.du_compressed_sams),
                         formatter.format(self.du_sams),
                         formatter.format(self.du_beds),
                         formatter.format(self.du_other),
                         formatter.format(self.du_total)])
        if self.is_subdir:
            dataline.append(', '.join(self.users))
        return dataline

class UsageFormatter:
    """Utility class for formatting disk usage sizes

    Create a UsageFormatter and configure via the 'units' argument.

    >>> f = UsageFormatter()

    creates a formatter that returns byte values in human-reable form,
    e.g. 4.0K, 165.0G etc.

    >>> f = UsageFormatter(units='G')

    creates a formatter that returns byte values in gigabytes, with no
    trailing unit identifier.

    To format values use e.g.:

    >>> f.format(1022)

    """
    def __init__(self,units=None):
        """Create a new UsageFormatter object

        'units' controls how the object formats sizes:

        None: returns byte values in human-reable form, with trailing
              unit identifier e.g. 4.0K, 165.0G
        'K','M','G','T': returns values in Kb, Mb, Gb or Tb's (no
              trailing unit identifier)

        """
        self.__units = units

    def format(self,size):
        """Return 'size' (in bytes) in the appropriate format

        """
        if self.__units is None:
            # Human-readable
            return bcf_utils.format_file_size(size)
        elif self.__units == 'bytes':
            # Raw bytes
            return size
        else:
            return bcf_utils.format_file_size(size,units=self.__units)[0:-1]

#######################################################################
# Functions
#######################################################################

def get_seqdir_dataline(seqdir,formatter):
    if not seqdir.is_subdir:
        dataline = [seqdir.year,
                    seqdir.platform,
                    seqdir.parent_dir,
                    seqdir.name]
    else:
        dataline = ['','','',seqdir.name]
    dataline.extend([formatter.format(seqdir.du_solid),
                     formatter.format(seqdir.du_fastqgzs),
                     formatter.format(seqdir.du_fastqs),
                     formatter.format(seqdir.du_compressed_bams),
                     formatter.format(seqdir.du_bams),
                     formatter.format(seqdir.du_compressed_sams),
                     formatter.format(seqdir.du_sams),
                     formatter.format(seqdir.du_beds),
                     formatter.format(seqdir.du_other),
                     formatter.format(seqdir.du_total)])
    if seqdir.is_subdir:
        dataline.append(', '.join(seqdir.users))
    return dataline

def plot_usage_data(usage_data,show=False,out_file="audit_data.png"):
    # Create a stacked bar plot
    import numpy as np
    # Need to specify a non-interactive backend
    import matplotlib
    matplotlib.use('Agg')
    # Import pyplot
    import matplotlib.pyplot as plt
    # Make the stacked bar plot
    ind = np.arange(len(years))    # the x locations for the groups
    width = 0.70       # the width of the bars: can also be len(x) sequence
    colours = ['blue',
               'orange',
               'yellow',
               'maroon',
               'green',
               'cyan',
               'black']
    data = {}
    for platform in platform_list:
        data[platform] = []
        for year in years:
            data[platform].append(usage_data[year][platform])
    plots = {}
    culmulative = []
    for i,platform in enumerate(platform_list):
        if not culmulative:
            plots[platform] = plt.bar(ind,data[platform],width,color=colours[i],
                                      linewidth=0)
            for j,year in enumerate(years):
                culmulative.append(data[platform][j])
        else:
            plots[platform] = plt.bar(ind,data[platform],width,color=colours[i],
                                      linewidth=0,
                                      bottom=culmulative)
            for j,year in enumerate(years):
                culmulative[j] += data[platform][j]

    plt.ylabel('Usage (Gb)')
    plt.title('Year')
    plt.xticks(ind+width/2.,years)
    plt.yticks(np.arange(0.,16000,2000))
    legend = plt.legend([plots[p][0] for p in platform_list],platform_list)
    legend.draw_frame(False)
    plt.savefig(out_file,format='png')
    if show:
        plt.show()

#######################################################################
# Main program
#######################################################################

if __name__ == "__main__":

    # Process command line
    p = optparse.OptionParser(usage="%prog DIR",
                              version="%prog "+__version__,
                              description="Audit sequence data directories under DIR.")
    p.add_option("--log-dir",action='store',dest='log_dir',default=os.getcwd(),
                 help="write log files from archiving operations to LOG_DIR (defaults "
                 "to cwd)")
    p.add_option("-o","--output-file",action='store',dest ='output_file',default=None,
                 help="specify output file to write results to (otherwise written to stdout)")
    p.add_option("--human-readable",action='store_true',dest ='human_readable',
                 default=False,help="output disk usage in human-readable format (otherwise "
                 "output in Gb)")
    p.add_option("--bytes",action="store_true",dest='bytes',default=False,
                 help="output sizes in bytes, rather than Gb")
    p.add_option("--years",action="store",dest='years',default=None,
                 help="examine data for specific year(s); can be a single year, or multiple "
                 "years expressed as either a range (e.g. '2011-2013') or a comma-separated "
                 "list (e.g.'2010,2014'). Default is current year.")
    p.add_option("--platforms",action="store",dest='platforms',default=None,
                 help="examine data for specific platform(s); can be a single platform, or "
                 "a comma-separated list. Default is all known platforms i.e. %s"
                 % ', '.join(platforms.list_platforms()))
    p.add_option("--include-subdirs",action='store_true',dest='include_subdirs',default=False,
                 help="also collect and report disk usage information for first level of "
                 "subdirectories in each data directory")
    p.add_option("--plot",action='store',dest='plot_file',default=None,
                 help="plot a stacked barchart of the usage and output as a PNG to PLOT_FILE")
    p.add_option("--debug",action='store_true',dest='debug',
                 help="turn on debugging output (nb very verbose!)")
    options,args = p.parse_args()
    if len(args) < 1:
        p.error("Need to supply at least one directory")

    # Debugging output
    if options.debug:
        logging.getLogger().setLevel(logging.DEBUG)

    # Formatter
    if options.human_readable:
        formatter = UsageFormatter()
    elif options.bytes:
        formatter = UsageFormatter('bytes')
    else:
        formatter = UsageFormatter('G')

    # Years and platforms
    if options.years is None:
        years = [int(time.strftime("%Y"))]
    else:
        if '-' in options.years:
            try:
                years = options.years.split('-')
                start_year = int(years[0])
                if years[1] != '':
                    end_year = int(years[1])
                else:
                    end_year = int(time.strftime("%Y"))
                years = range(start_year,end_year+1)
            except Exception, ex:
                p.error("Bad year range supplied to --years option")
        elif ',' in options.years:
            try:
                years = [int(x) for x in options.years.split(',')]
            except Exception, ex:
                p.error("Bad year list supplied to --years option")
        else:
            years = [int(options.years)]
    if options.platforms is None:
        platform_list = platforms.list_platforms()
    else:
        platform_list = options.platforms.split(',')
        for platform in platform_list:
            if platform not in platforms.list_platforms():
                p.error("Unknown platform '%s' supplied to --platform option" % platform)

    # Report before starting
    print "Years:\t%s" % ','.join([str(x) for x in years])
    print "Platforms:\t%s" % ','.join(platform_list)

    # Process directories
    seqdirs = []
    for d in args:
        print "Processing %s..." % d
        for year in years:
            print "- Year %s" % year
            dirn = os.path.join(d,"%s" % year)
            if not os.path.isdir(dirn):
                continue
            for platform in platform_list:
                platform_dirn = os.path.join(dirn,platform)
                if not os.path.isdir(platform_dirn):
                    continue
                print "-- Platform %s" % platform
                # Get all the directories for this year/platform combination
                for run in auto_process_utils.list_dirs(platform_dirn):
                    print "--- Run %s" % run
                    run_dir = os.path.join(platform_dirn,run)
                    if os.path.islink(run_dir):
                        logging.warning("%s: is link, ignoring" % run_dir)
                    else:
                        seqdir = SeqDataSizes(run,run_dir,year=year,platform=platform,
                                              include_subdirs=options.include_subdirs)
                        seqdir.get_disk_usage()
                        seqdirs.append(seqdir)

    # Calculate totals for each year
    usage = dict()
    for year in years:
        usage[year] = dict()
        usage[year]['total'] = 0
        for platform in platform_list:
            usage[year][platform] = 0
            for seqdir in seqdirs:
                if seqdir.year == year and seqdir.platform == platform:
                    usage[year]['total'] += seqdir.du_total
                    usage[year][platform] += seqdir.du_total

    # Report usage
    xls = simple_xls.XLSWorkBook("Usage stats")
    header_style = simple_xls.XLSStyle(bold=True)
    # Summary of all years and platforms
    summary = xls.add_work_sheet("Summary")
    header = ['#Year']
    for platform in platform_list:
        header.append(platform)
    header.append('Total')
    summary.insert_row(1,header,style=header_style)
    for year in years:
        line = [year]
        for platform in platform_list:
            line.append(formatter.format(usage[year][platform]))
        line.append(formatter.format(usage[year]['total']))
        summary.append_row(data=line)
    print summary.render_as_text()
    # Detailed breakdown
    for year in years:
        report = xls.add_work_sheet(year)
        report.insert_row(1,data=['#Year','Platform','Directory','Run',
                                  'Csfasta/qual',
                                  'Fastq.gz','Fastq',
                                  'Bam.(gz|bz2)','Bam',
                                  'Sam.(gz|bz2)','Sam',
                                  'Bed','Other',
                                  'Total Size (GB)'],
                          style=header_style)
        for seqdir in seqdirs:
            if seqdir.year != year:
                continue
            report.append_row(data=seqdir.dataline(formatter))
            for subdir in seqdir.subdirs:
                report.append_row(data=subdir.dataline(formatter))
            report.append_row(data=('.',))
        print report.render_as_text()
    # Save to file
    xls.save_as_xls("usage_stats.xls")

    # Make a stacked bar char plot
    if options.plot_file is not None:
        plot_usage_data(usage,out_file=options.plot_file)
