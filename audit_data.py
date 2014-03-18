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

__version__ = '0.0.1'

KNOWN_PLATFORMS = ['solid4',
                   'solid5500',
                   'illumina-ga2x',
                   'miseq',
                   'hiseq',
                   'other']

#######################################################################
# Import modules that this module depends on
#######################################################################

import applications
import simple_scheduler
import auto_process_utils
import JobRunner
import bcf_utils
import logging
import optparse
import sys
import os
import time

#######################################################################
# Classes
#######################################################################

class SeqDataSizes:
    """Class for acquiring & storing sequence data size (i.e disk usage) info

    """
    def __init__(self,name,dirn,sched,year=None,platform=None,is_subdir=False):
        """Create a new SeqDataDirectorySizes object

        name: name for the object
        dirn: path to the directory
        sched: SimpleScheduler instance

        """
        self.name = name
        self.dirn = dirn
        self.scheduler = sched
        self.year = year
        self.platform = platform
        self.du_total = None
        self.du_fastqs = None
        self.du_fastqgzs = None
        self.is_subdir = is_subdir
        self.subdirs = []
        if not is_subdir:
            for d in auto_process_utils.list_dirs(self.dirn):
                self.subdirs.append(SeqDataSizes(d,os.path.join(self.dirn,d),
                                                 self.scheduler,
                                                 year=self.year,
                                                 platform=self.platform,
                                                 is_subdir=True))

    @property
    def du_other(self):
        try:
            return self.du_total - self.du_fastqs - self.du_fastqgzs
        except TypeError:
            return None

    @property
    def has_data(self):
        has_data = self.du_total is not None and \
                   self.du_fastqs is not None and \
                   self.du_fastqgzs is not None
        if not self.is_subdir:
            for subdir in self.subdirs:
                has_data = has_data and subdir.has_data
        return has_data

    def get_disk_usage(self):
        if self.is_subdir:
            name = "%s.%s" % (os.path.basename(os.path.dirname(self.dirn)),self.name)
        else:
            name = self.name
        job_grp = self.scheduler.group(name,callbacks=(self.process_directory_sizes,))
        job_grp.add(applications.Command('du','-sb',self.dirn),
                    name="du.%s.total" % name)
        job_grp.add(applications.Command('find',self.dirn,
                                         '-name','"*.fastq.gz"',
                                         '-exec','du','-b','{}','\;'),
                    name="du.%s.fastqgzs" % name)
        job_grp.add(applications.Command('find',self.dirn,
                                         '-name','"*.fastq"',
                                         '-exec','du','-b','{}','\;'),
                    name="du.%s.fastqs" % name)
        job_grp.submit()
        for subdir in self.subdirs:
            subdir.get_disk_usage()

    def process_directory_sizes(self,name,job_groups,sched):
        # Callback function
        logging.debug("### Handling %s ###" % name)
        for group in job_groups:
            for job in group.jobs:
                size = get_total_size(job.log)
                if job.name.endswith("total"):
                    self.du_total = size
                elif job.name.endswith("fastqgzs"):
                    self.du_fastqgzs = size
                elif job.name.endswith("fastqs"):
                    self.du_fastqs = size
                else:
                    raise Exception,"Unknown job '%s'" % job.name
                logging.debug("Deleting file %s" % job.log)
                os.remove(job.log)

    def report(self,no_header=False,pretty_print=True):
        report = []
        if not no_header:
            report.append("#Year\tPlatform\tDirectory\tTotal Size\tFastq.gz size\tFastq size\tOther files")
        if pretty_print:
            report.append("%s\t%s\t%s\t%s\t%s\t%s\t%s" % (
                self.year,
                self.platform,
                self.name,
                bcf_utils.format_file_size(self.du_total),
                bcf_utils.format_file_size(self.du_fastqgzs),
                bcf_utils.format_file_size(self.du_fastqs),
                bcf_utils.format_file_size(self.du_other)))
        else:
            report.append("%s\t%s\t%s\t%d\t%d\t%d\t%d" % (
                self.year,
                self.platform,
                self.name,
                self.du_total,
                self.du_fastqgzs,
                self.du_fastqs,
                self.du_other))
        report = '\n'.join(report)
        for subdir in self.subdirs:
            report += '\n' + subdir.report(no_header=True,pretty_print=pretty_print)
        return report

#######################################################################
# Functions
#######################################################################

def get_total_size(du_log_file):
    """Read 'du' output from log_file and return the total 

    """
    total_size = 0
    for line in open(du_log_file,'r'):
        size = int(line.split()[0])
        total_size += size
    return total_size

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
    p.add_option("--use-grid-engine",action='store_true',dest='use_grid_engine',default=False,
                 help="submit 'du' jobs to Grid Engine")
    p.add_option("-o","--output-file",action='store',dest ='output_file',default=None,
                 help="specify output file to write results to (otherwise written to stdout)")
    p.add_option("--human-readable",action='store_true',dest ='human_readable',
                 default=False,help="output disk usage in human-readable format (otherwise "
                 "output in bytes)")
    p.add_option("--years",action="store",dest='years',default=None,
                 help="examine data for specific year(s); can be a single year, or multiple "
                 "years expressed as either a range (e.g. '2011-2013') or a comma-separated "
                 "list (e.g.'2010,2014'). Default is current year.")
    p.add_option("--platforms",action="store",dest='platforms',default=None,
                 help="examine data for specific platform(s); can be a single platform, or "
                 "a comma-separated list. Default is all known platforms i.e. %s"
                 % KNOWN_PLATFORMS)
    p.add_option("--debug",action='store_true',dest='debug',
                 help="turn on debugging output (nb very verbose!)")
    options,args = p.parse_args()
    if len(args) < 1:
        p.error("Need to supply at least one directory")

    # Debugging output
    if options.debug:
        logging.getLogger().setLevel(logging.DEBUG)

    # Make a job runner for computationally-intensive (CI) jobs
    if options.use_grid_engine:
        ci_runner = JobRunner.GEJobRunner(log_dir=options.log_dir,
                                          ge_extra_args=['-j','y'])
    else:
        ci_runner = JobRunner.SimpleJobRunner(log_dir=options.log_dir,
                                              join_logs=True)

    # Start scheduler
    s = simple_scheduler.SimpleScheduler(runner=ci_runner,
                                         max_concurrent=16)
    s.start()

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
            years = int(options.years)
    if options.platforms is None:
        platforms = KNOWN_PLATFORMS
    else:
        platforms = options.platforms.split(',')
        for platform in platforms:
            if platform not in KNOWN_PLATFORMS:
                p.error("Unknown platform '%s' supplied to --platform option" % platform)

    # Report before starting
    print "Years:\t%s" % ','.join([str(x) for x in years])
    print "Platforms:\t%s" % ','.join(platforms)

    # Process directories
    seqdirs = []
    for d in args:
        print "Processing %s..." % d
        for year in years:
            dirn = os.path.join(d,"%s" % year)
            if not os.path.isdir(dirn):
                continue
            for platform in platforms:
                platform_dirn = os.path.join(dirn,platform)
                if not os.path.isdir(platform_dirn):
                    continue
                print "%s" % platform_dirn
                # Get all the directories for this year/platform combination
                for run in auto_process_utils.list_dirs(platform_dirn):
                    print "%s" % run
                    run_dir = os.path.join(platform_dirn,run)
                    seqdir = SeqDataSizes(run,run_dir,s,year=year,platform=platform)
                    seqdir.get_disk_usage()
                    seqdirs.append(seqdir)
    # Wait for all data to be processed
    while not reduce(lambda x,y: x and y.has_data,seqdirs,True):
        time.sleep(5)

    # Report usage for all directories
    if options.output_file is not None:
        fp = open(options.output_file,'w')
    else:
        fp = sys.stdout
    no_header = False
    for seqdir in seqdirs:
        fp.write(seqdir.report(pretty_print=options.human_readable,no_header=no_header)+'\n')
        if not no_header:
            no_header = True

    # Report usage by year and platform
    usage = dict()
    for year in years:
        usage[year] = dict()
        usage[year]['total'] = 0
        for platform in platforms:
            usage[year][platform] = 0
            for seqdir in seqdirs:
                if seqdir.year == year and seqdir.platform == platform:
                    usage[year]['total'] += seqdir.du_total
                    usage[year][platform] += seqdir.du_total
    header = ['#Year']
    for platform in platforms:
        header.append(platform)
    header.append('Total')
    print '\t'.join(header)
    for year in years:
        line = [str(year)]
        for platform in platforms:
            if options.human_readable:
                line.append(bcf_utils.format_file_size(usage[year][platform]))
            else:
                line.append(str(usage[year][platform]))
        if options.human_readable:
            line.append(bcf_utils.format_file_size(usage[year]['total']))
        else:
            line.append(str(usage[year]['total']))
        print '\t'.join(line)

        

