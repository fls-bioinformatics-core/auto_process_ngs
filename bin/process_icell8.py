#!/usr/bin/env python
#
#     process_icell8.py: perform processing of Wafergen iCell8 data
#     Copyright (C) University of Manchester 2017 Peter Briggs
#
"""
process_icell8.py

Utility to peform initial processing of data from Wafergen iCell8
platform.
"""

######################################################################
# Imports
######################################################################

import os
import sys
import logging
import argparse
import glob
import time
import shutil
import uuid
from bcftbx.utils import mkdir
from bcftbx.utils import strip_ext
from bcftbx.IlluminaData import IlluminaData
from bcftbx.IlluminaData import IlluminaDataError
from bcftbx.JobRunner import fetch_runner
from auto_process_ngs.applications import Command
from auto_process_ngs.simple_scheduler import SimpleScheduler
from auto_process_ngs.simple_scheduler import SchedulerReporter
from auto_process_ngs.simple_scheduler import SchedulerGroup
from auto_process_ngs.fastq_utils import pair_fastqs
from auto_process_ngs.utils import AnalysisFastq
import auto_process_ngs.envmod as envmod

# Fetch configuration settings
import auto_process_ngs.settings
__settings = auto_process_ngs.settings.Settings()

######################################################################
# Magic numbers
######################################################################

DEFAULT_BATCH_SIZE = 5000000

######################################################################
# Generic pipeline base classes
######################################################################

class PipelineStage(object):
    """
    """
    def __init__(self,name):
        self._name = str(name)
        self._commands = []
    def name(self):
        return self._name.lower().replace(' ','_')
    def add(self,pipeline_job):
        self._commands.append(pipeline_job)
    def run(self,sched=None,runner=None,working_dir=None,log_dir=None,
            scripts_dir=None,wait_for=(),async=True,use_wrapper=False):
        # Generate commands to run
        cmds = []
        for command in self._commands:
            if use_wrapper:
                script_file = command.make_wrapper_script(scripts_dir=scripts_dir)
                cmd = Command('/bin/bash',script_file)
            else:
                cmd = command.cmd()
            cmds.append(cmd)
        # Run the commands
        if sched is not None:
            # Use the scheduler
            use_group = (len(cmds)!=1)
            if use_group:
                # Run as a group
                group = sched.group(self.name())
                for cmd in cmds:
                    name = "%s.%s" % (self.name(),uuid.uuid4())
                    group.add(cmd,
                              wd=working_dir,
                              name=name,
                              runner=runner,
                              log_dir=log_dir,
                              wait_for=wait_for)
                group.close()
                callback_name = group.name
                callback_function = check_status
            else:
                # Run a single job
                cmd = cmds[0]
                name = "%s.%s" % (self.name(),uuid.uuid4())
                job = sched.submit(cmd,
                                   wd=working_dir,
                                   name=name,
                                   runner=runner,
                                   log_dir=log_dir,
                                   wait_for=wait_for)
                callback_name = job.name
                callback_function = check_status
            # If asynchronous then setup callback and
            # return immediately
            if async:
                sched.callback("%s" % self._name,
                               callback_function,
                               wait_for=(callback_name,))
            else:
                # Wait for job or group to complete before returning
                sched.wait_for((callback_name,))
                if use_group:
                    exit_code = max([j.exit_code for j in group.jobs])
                else:
                    exit_code = job.exit_code
                if exit_code != 0:
                    logging.warning("%s: failed (exit code %s)" %
                                    (self._name,exit_code))
            if use_group:
                return group
            else:
                return job
        else:
            # Run each stage locally
            for cmd in cmds:
                cmd.run_subprocess(working_dir=working_dir)

class PipelineCommand(object):
    """
    """
    def __init__(self,name):
        self._name = str(name)
    def cmd(self):
        # Build the command
        # Must be implemented by the subclass and return a
        # Command instance
        raise NotImplementedError("Subclass must implement 'cmd' method")
    def name(self):
        return self._name.lower().replace(' ','_')
    def make_wrapper_script(self,scripts_dir=None,shell="/bin/bash"):
        # Wrap in a script
        if scripts_dir is None:
            scripts_dir = os.getcwd()
        script_file = os.path.join(scripts_dir,"%s.%s.sh" % (self.name(),
                                                             uuid.uuid4()))
        self.cmd().make_wrapper_script(filen=script_file,
                                       shell=shell)
        return script_file

######################################################################
# ICell8 pipeline commands
######################################################################

class ICell8Statistics(PipelineCommand):
    """
    """
    def __init__(self,name,fastqs,stats_file,well_list=None,
                 suffix=None,append=False):
        PipelineCommand.__init__(self,name)
        self._fastqs = fastqs
        self._stats_file = os.path.abspath(stats_file)
        self._well_list = well_list
        if self._well_list is not None:
            self._well_list = os.path.abspath(self._well_list)
        self._append = append
        self._suffix = suffix
    def cmd(self):
        # Build command
        cmd = Command('icell8_stats.py',
                      '-f',self._stats_file)
        if self._well_list:
            cmd.add_args('-w',self._well_list)
        if self._suffix:
            cmd.add_args('--suffix',self._suffix)
        if self._append:
            cmd.add_args('--append',)
        cmd.add_args(*self._fastqs)
        return cmd

class SplitAndFilterFastqPair(PipelineCommand):
    """
    """
    def __init__(self,name,fastq_pair,filter_dir,well_list=None,
                 basename=None,mode='none',filter=False):
        PipelineCommand.__init__(self,name)
        self._fastq_pair = fastq_pair
        self._filter_dir = os.path.abspath(filter_dir)
        self._well_list = well_list
        if self._well_list is not None:
            self._well_list = os.path.abspath(self._well_list)
        self._basename = basename
        self._mode = mode
        self._filter = filter
    def cmd(self):
        # Build command
        cmd = Command('split_icell8_fastqs.py',
                      '-o',self._filter_dir,
                      '-b',self._basename)
        if self._well_list:
            cmd.add_args('-w',self._well_list)
        if self._mode:
            cmd.add_args('-m',self._mode)
        if self._filter:
            cmd.add_args('--filter')
        cmd.add_args(*self._fastq_pair)
        return cmd

class BatchFastqs(PipelineCommand):
    """
    """
    def __init__(self,name,fastqs,batch_dir,basename,
                 batch_size=None,merge=False):
        PipelineCommand.__init__(self,name)
        self._fastqs = fastqs
        self._batch_dir = os.path.abspath(batch_dir)
        self._basename = basename
        self._batch_size = batch_size
        self._merge = merge
    def cmd(self):
        # Build command
        cmd = Command('batch_fastqs.py',
                      '-o',self._batch_dir,
                      '-b',self._basename)
        if self._batch_size:
            cmd.add_args('-s',self._batch_size)
        if self._merge:
            cmd.add_args('-m')
        cmd.add_args(*self._fastqs)
        return cmd

class TrimFastqPair(PipelineCommand):
    """
    """
    def __init__(self,name,fastq_pair,trim_dir):
        PipelineCommand.__init__(self,name)
        self._fastq_pair = fastq_pair
        self._trim_dir = os.path.abspath(trim_dir)
    def cmd(self):
        # Generate output file pair names
        fastq_pair_out = [os.path.join(self._trim_dir,
                                       strip_ext(os.path.basename(fq),'.fastq')
                                       + '.trimmed.fastq')
                          for fq in self._fastq_pair]
        # Build command
        cmd = Command(
            'cutadapt',
            '-a','AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT',
            '-a','AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT',
            '-a','AAAAAAAA',
            '-a','TTTTTTTT',
            '-m',20,
            '--trim-n',
            '--max-n',0.7,
            '-q',25)
        # NB reverse R1 and R2 for input and output
        cmd.add_args('-o',fastq_pair_out[1],
                     '-p',fastq_pair_out[0])
        cmd.add_args(self._fastq_pair[1],
                     self._fastq_pair[0])
        return cmd

class ContaminantFilterFastqPair(PipelineCommand):
    """
    """
    def __init__(self,name,fastq_pair,filter_dir,
                 preferred_conf,contaminants_conf,
                 aligner=None,threads=None):
        PipelineCommand.__init__(self,name)
        self._fastq_pair = fastq_pair
        self._filter_dir = os.path.abspath(filter_dir)
        self._preferred_conf = os.path.abspath(preferred_conf)
        self._contaminants_conf = os.path.abspath(contaminants_conf)
        self._aligner = aligner
        self._threads = threads
    def cmd(self):
        # Build the command
        cmd = Command(
            'icell8_contamination_filter.py',
            '-p',self._preferred_conf,
            '-c',self._contaminants_conf,
            '-o',self._filter_dir)
        if self._threads:
            cmd.add_args('-n',self._threads)
        if self._aligner is not None:
            cmd.add_args('-a',self._aligner)
        cmd.add_args(*self._fastq_pair)
        return cmd

######################################################################
# ICell8 pipeline stages
######################################################################

class GetICell8Stats(PipelineStage):
    """
    """
    def __init__(self,name,*args,**kws):
        PipelineStage.__init__(self,name)
        self.add(ICell8Statistics(self._name,
                                  *args,
                                  **kws))

class SplitFastqsIntoBatches(PipelineStage):
    """
    """
    def __init__(self,name,*args,**kws):
        PipelineStage.__init__(self,name)
        self.add(BatchFastqs(self._name,
                             *args,
                             **kws))

class FilterICell8Fastqs(PipelineStage):
    """
    """
    def __init__(self,name,fastqs,filter_dir,well_list=None,
                 mode='none',filter=False):
        PipelineStage.__init__(self,name)
        fastq_pairs = pair_fastqs(fastqs)[0]
        for fastq_pair in fastq_pairs:
            basename = os.path.basename(fastq_pair[0])[:-len(".r1.fastq")]
            self.add(SplitAndFilterFastqPair(self._name,
                                             fastq_pair,
                                             filter_dir,
                                             well_list=well_list,
                                             basename=basename,
                                             mode=mode,
                                             filter=filter))

class TrimReads(PipelineStage):
    """
    """
    def __init__(self,name,fastqs,trim_dir):
        PipelineStage.__init__(self,name)
        fastq_pairs = pair_fastqs(fastqs)[0]
        for fastq_pair in fastq_pairs:
            self.add(TrimFastqPair(self._name,
                                   fastq_pair,
                                   trim_dir))

class FilterContaminatedReads(PipelineStage):
    """
    """
    def __init__(self,name,fastqs,filter_dir,preferred_conf,
                 contaminants_conf,aligner=None,threads=None):
        PipelineStage.__init__(self,name)
        fastq_pairs = pair_fastqs(fastqs)[0]
        for fastq_pair in fastq_pairs:
            self.add(ContaminantFilterFastqPair(self._name,
                                                fastq_pair,
                                                filter_dir,
                                                preferred_conf,
                                                contaminants_conf,
                                                aligner=aligner,
                                                threads=threads))

class SplitByBarcodes(PipelineStage):
    """
    """
    def __init__(self,name,fastqs,barcodes_dir):
        PipelineStage.__init__(self,name)
        fastq_pairs = pair_fastqs(fastqs)[0]
        for fastq_pair in fastq_pairs:
            basename = os.path.basename(fastq_pair[0])[:-len(".r1.fastq")+1]
            self.add(SplitAndFilterFastqPair(self._name,
                                             fastq_pair,
                                             barcodes_dir,
                                             basename=basename,
                                             mode="barcodes"))

class MergeFastqs(PipelineStage):
    """
    """
    def __init__(self,name,fastqs,unassigned_fastqs,
                 failed_barcode_fastqs,failed_umi_fastqs,
                 merge_dir,basename,batch_size=25):
        PipelineStage.__init__(self,name)
        # Extract the barcodes from the fastq names
        barcodes = set()
        for fq in barcoded_fastqs:
            barcode = os.path.basename(fq).split('.')[-3]
            barcodes.add(barcode)
        barcodes = sorted(list(barcodes))
        # Group files by barcode
        fastq_groups = dict()
        for barcode in barcodes:
            fqs = filter(lambda fq: (fq.endswith("%s.r1.fastq" % barcode) or
                                     fq.endswith("%s.r2.fastq" % barcode)),
                         fastqs)
            fastq_groups[barcode] = fqs
        # Group barcodes into batches
        barcode_batches = [barcodes[i:i+batch_size]
                           for i in xrange(0,len(barcodes),batch_size)]
        # Concat fastqs
        for i,barcode_batch in enumerate(barcode_batches):
            batch_name = "barcodes%06d" % i
            print "Barcode batch: %s" % batch_name
            fastq_pairs = []
            for barcode in barcode_batch:
                print "-- %s" % barcode
                fastq_pairs.extend(fastq_groups[barcode])
            self.add(SplitAndFilterFastqPair(self._name,
                                             fastq_pairs,
                                             merge_dir,
                                             basename=basename,
                                             mode="barcodes"))
        # Handle unassigned and failed quality reads
        for fqs in (unassigned_fastqs,
                    failed_barcode_fastqs,
                    failed_umi_fastqs):
            if not fqs:
                continue
            self.add(BatchFastqs(self._name,fqs,merge_dir,basename))

######################################################################
# Functions
######################################################################

def collect_fastqs(dirn,pattern):
    """
    Return names of Fastqs in a directory which match a glob pattern

    Arguments:
      dirn (str): path to a directory containing the files
      pattern (str): a glob pattern to match

    Returns:
      List: list of matching files
    """
    return sorted(glob.glob(os.path.join(os.path.abspath(dirn),pattern)))

def check_status(name,jobs,sched):
    """
    Check the exit status of a set of groups

    This function is not called directly, instead it should be
    passed as part of a callback to check the exit status of a
    job or group of jobs from the scheduler.

    For example:

    >>> sched.callback("My job",check_status,wait_for(('my_job',))

    If any of the jobs have a non-zero exit status then
    the program is terminated.

    Arguments:
      name (str): name for the callback
      jobs (list): list of SchedulerJob instances
      sched (SimpleScheduler): scheduler instance

    """
    print "*** %s completed ***" % name
    for job in jobs:
        exit_code = job.exit_code
        if exit_code != 0:
            logging.critical("Job '%s' failed: exit code %s"
                             % (job.name,exit_code))
            sched.stop()
            sys.exit(1)
    print "%s: ok" % name

######################################################################
# Main
######################################################################

if __name__ == "__main__":
    # Handle the command line
    p = argparse.ArgumentParser(
        description="Perform initial QC on FASTQs from Wafergen "
        "ICell8: assign to barcodes, filter on barcode & UMI quality, "
        "trim reads, perform contaminant filtering and split by "
        "barcode.")
    p.add_argument("well_list",metavar="WELL_LIST",help="Well list file")
    p.add_argument("fastqs",nargs='*',metavar="FASTQ_R1 FASTQ_R2",
                   help="FASTQ file pairs")
    p.add_argument("--unaligned",
                   dest="unaligned_dir",default="bcl2fastq",
                   help="'unaligned' dir with output from "
                   "bcl2fastq")
    p.add_argument("-o","--outdir",
                   dest="outdir",default="icell8",
                   help="directory to write outputs to "
                   "(default: 'CWD/icell8')")
    p.add_argument("-p","--preferred",
                   dest="preferred_conf",
                   help="fastq_screen 'conf' file with the "
                   "'preferred' genome indices")
    p.add_argument("-c","--contaminants",
                   dest="contaminants_conf",
                   help="fastq_screen 'conf' file with the "
                   "'contaminant' genome indices")
    p.add_argument("-a","--aligner",
                   dest="aligner",default=None,
                   choices=["bowtie","bowtie2"],
                   help="aligner to use with fastq_screen (default: "
                   "don't specify the aligner)")
    p.add_argument("-n","--threads",type=int,
                   dest="threads",default=1,
                   help="number of threads to use with fastq_screen "
                   "(default: 1)")
    p.add_argument("-r","--runner",metavar="STAGE=RUNNER",
                   action="append",dest="runners",default=list(),
                   help="explicitly specify runner definitions for "
                   "running pipeline jobs at each stage. STAGE "
                   "can be one of 'default','contaminant_filter'. "
                   "RUNNER must be a valid job runner specification "
                   "e.g. 'GEJobRunner(-j y)'. Multiple --runner "
                   "arguments can be specified (default: "
                   "SimpleJobRunner)")
    p.add_argument("-s","--size",type=int,
                   dest="batch_size",default=DEFAULT_BATCH_SIZE,
                   help="number of reads per batch when splitting "
                   "FASTQ files for processing (default: %s)" %
                   DEFAULT_BATCH_SIZE)
    p.add_argument("-m","--max-jobs",type=int,
                   dest="max_jobs",
                   default= __settings.general.max_concurrent_jobs,
                   help="maxiumum number of concurrent jobs to run "
                   "(default: %d)"
                   % __settings.general.max_concurrent_jobs)
    p.add_argument('--modulefiles',action='store',
                   dest='modulefiles',default=None,
                   help="comma-separated list of environment "
                   "modules to load before executing commands "
                   "(overrides any modules specified in the global "
                   "settings)")
    p.add_argument('--force',action='store_true',
                   dest='force',default=False,
                   help="force overwrite of existing outputs")
    args = p.parse_args()

    # Deal with module files
    if args.modulefiles is not None:
        modulefiles = args.modulefiles.split(',')
        for modulefile in modulefiles:
            envmod.load(modulefile)

    # Deal with job runners
    stages = ('default','contaminant_filter')
    runners = dict()
    for runner in args.runners:
        try:
            stage,runner_spec = runner.split('=')
        except ValueError: # too few values to unpack
            stage = 'default'
            runner_spec = runner
        if stage not in stages:
            logging.fatal("Bad stage for --runner option: %s" % stage)
            sys.exit(1)
        runners[stage] = fetch_runner(runner_spec)
    try:
        default_runner = runners['default']
    except KeyError:
        default_runner = fetch_runner('SimpleJobRunner')
    for stage in stages:
        if stage not in runners:
            runners[stage] = default_runner

    # Other settings
    well_list = os.path.abspath(args.well_list)
    max_jobs = args.max_jobs

    # Report settings
    print "Unaligned dir     : %s" % args.unaligned_dir
    print "Well list file    : %s" % well_list
    print "Output dir        : %s" % args.outdir
    print "Batch size (reads): %s" % args.batch_size
    print "Preferred genomes screen: %s" % args.preferred_conf
    with open(args.preferred_conf) as fp:
        for line in fp:
            if line.startswith("DATABASE"):
                print "-- %s" % line.split('\t')[1]
    print "Contaminants screen     : %s" % args.contaminants_conf
    with open(args.contaminants_conf) as fp:
        for line in fp:
            if line.startswith("DATABASE"):
                print "-- %s" % line.split('\t')[1]
    print "Fastq_screen aligner    : %s" % args.aligner
    print "Fastq_screen threads    : %s" % args.threads
    print "Maximum concurrent jobs : %s" % max_jobs
    print "Job runners:"
    for stage in stages:
        print "-- %s: %s" % (stage,runners[stage])
    if args.modulefiles is not None:
        print "Environment modules:"
        for modulefile in modulefiles:
            print "-- %s" % modulefile

    # Get the input FASTQ file pairs
    fastqs = []
    try:
        illumina_data = IlluminaData(os.getcwd(),
                                     unaligned_dir=args.unaligned_dir)
        for project in illumina_data.projects:
            for sample in project.samples:
                for fq in sample.fastq:
                    fastqs.append(os.path.join(sample.dirn,fq))
    except IlluminaDataError:
        logging.warning("Couldn't find FASTQS in directory '%s'" %
                        args.unaligned_dir)
    for fq in args.fastqs:
        fastqs.append(os.path.abspath(fq))
    if not fastqs:
        logging.fatal("No FASTQs found")
        sys.exit(1)

    # Basename for output fastqs and job names etc
    basename = AnalysisFastq(fastqs[0]).sample_name

    # Set up a scheduler for running jobs
    sched_reporter = SchedulerReporter(
        job_start="Started  #%(job_number)d: %(job_name)s:\n-- %(command)s",
        job_end=  "Finished #%(job_number)d: %(job_name)s"
    )
    sched = SimpleScheduler(runner=runners['default'],
                            max_concurrent=max_jobs,
                            reporter=sched_reporter)
    sched.start()

    # Make top-level output dirs
    icell8_dir = os.path.abspath(args.outdir)
    if os.path.exists(icell8_dir):
        if not args.force:
            logging.fatal("Output destination '%s': already exists "
                          "(remove or use --force to overwrite)" %
                          icell8_dir)
            sys.exit(1)
        logging.warning("Removing existing output destination '%s'" %
                        icell8_dir)
        shutil.rmtree(icell8_dir)
    log_dir = os.path.join(icell8_dir,"logs")
    stats_dir = os.path.join(icell8_dir,"stats")
    scripts_dir = os.path.join(icell8_dir,"scripts")
    for dirn in (icell8_dir,log_dir,stats_dir,scripts_dir):
        mkdir(dirn)

    # Initial stats
    stats_file = os.path.join(stats_dir,"icell8_stats.tsv")
    initial_stats = GetICell8Stats("Initial statistics",
                                   fastqs,stats_file,well_list)
    initial_stats = initial_stats.run(sched,
                                      working_dir=icell8_dir,
                                      log_dir=log_dir)

    # Split fastqs into batches
    batch_dir = os.path.join(icell8_dir,"_fastqs.batched")
    mkdir(batch_dir)
    batch_fastqs = SplitFastqsIntoBatches("Batch Fastqs",fastqs,
                                          batch_dir,basename,
                                          batch_size=args.batch_size)
    batch_fastqs = batch_fastqs.run(sched,
                                    working_dir=icell8_dir,
                                    log_dir=log_dir,
                                    async=False)
    print "*** Fastq batching stage completed ***"

    # Collect the batched files for processing
    batched_fastqs = collect_fastqs(batch_dir,"*.B*.r*.fastq")

    # Setup the quality filter jobs as a group
    filter_dir = os.path.join(icell8_dir,"_fastqs.quality_filter")
    mkdir(filter_dir)
    quality_filter = FilterICell8Fastqs("Quality filter Fastqs",
                                        batched_fastqs,
                                        filter_dir,
                                        well_list=well_list,
                                        mode='none',
                                        filter=True)
    quality_filter = quality_filter.run(sched,
                                        working_dir=filter_dir,
                                        log_dir=log_dir,
                                        async=False)
    print "*** Quality filtering stage completed ***"

    # Collect the post-filtering files for stats and read trimming
    filtered_fastqs = collect_fastqs(filter_dir,"*.B*.filtered.r*.fastq")
    
    # Collect the files with unassigned and failed quality reads
    unassigned_fastqs = collect_fastqs(filter_dir,"*.unassigned.r*.fastq")
    failed_barcode_fastqs = collect_fastqs(filter_dir,"*.failed_barcode.r*.fastq")
    failed_umi_fastqs = collect_fastqs(filter_dir,"*.failed_umi.r*.fastq")
    
    # Post quality filter stats
    filter_stats = GetICell8Stats("Post-quality filter statistics",
                                  filtered_fastqs,stats_file,
                                  suffix="_quality_filtered",
                                  append=True)
    filter_stats = filter_stats.run(sched,
                                    working_dir=icell8_dir,
                                    log_dir=log_dir,
                                    scripts_dir=scripts_dir,
                                    use_wrapper=True,
                                    wait_for=(initial_stats.name,))

    # Set up the cutadapt jobs as a group
    trim_dir = os.path.join(icell8_dir,"_fastqs.trim_reads")
    mkdir(trim_dir)
    trim_reads = TrimReads("Read trimming",filtered_fastqs,trim_dir)
    trim_reads = trim_reads.run(sched,
                                working_dir=trim_dir,
                                log_dir=log_dir,
                                async=False)
    print "*** Read trimming stage completed ***"

    # Collect the files for stats and contamination filter step
    trimmed_fastqs = collect_fastqs(trim_dir,"*.trimmed.fastq")

    # Post read trimming stats
    trim_stats = GetICell8Stats("Post-trimming statistics",
                                trimmed_fastqs,stats_file,
                                suffix="_trimmed",
                                append=True)
    trim_stats = trim_stats.run(sched,
                                working_dir=icell8_dir,
                                log_dir=log_dir,
                                scripts_dir=scripts_dir,
                                use_wrapper=True,
                                wait_for=(filter_stats.name,))

    # Set up the contaminant filter jobs as a group
    contaminant_filter_dir = os.path.join(icell8_dir,
                                          "_fastqs.contaminant_filter")
    mkdir(contaminant_filter_dir)
    contaminant_filter = FilterContaminatedReads("Contaminant filtering",
                                                 trimmed_fastqs,
                                                 contaminant_filter_dir,
                                                 args.preferred_conf,
                                                 args.contaminants_conf,
                                                 aligner=args.aligner,
                                                 threads=args.threads)
    contaminant_filter = contaminant_filter.run(
        sched,
        working_dir=contaminant_filter_dir,
        scripts_dir=scripts_dir,
        log_dir=log_dir,
        runner=runners['contaminant_filter'],
        use_wrapper=True,
        async=False)
    print "*** Contaminant filtering stage completed ***"

    # Post contaminant filter stats
    filtered_fastqs = collect_fastqs(contaminant_filter_dir,"*.trimmed.filtered.fastq")
    final_stats = GetICell8Stats("Post-contaminant filter statistics",
                                 filtered_fastqs,stats_file,
                                 suffix="_contaminant_filtered",
                                 append=True)
    final_stats = final_stats.run(sched,
                                  working_dir=icell8_dir,
                                  log_dir=log_dir,
                                  scripts_dir=scripts_dir,
                                  use_wrapper=True,
                                  wait_for=(trim_stats.name,))

    # Rebatch reads by barcode
    # First: split each batch by barcode
    barcoded_fastqs_dir = os.path.join(icell8_dir,"_fastqs.barcodes")
    mkdir(barcoded_fastqs_dir)
    split_barcodes = SplitByBarcodes("Split by barcodes",
                                     filtered_fastqs,
                                     barcoded_fastqs_dir)
    split_barcodes = split_barcodes.run(sched,
                                        working_dir=barcoded_fastqs_dir,
                                        log_dir=log_dir,
                                        scripts_dir=scripts_dir,
                                        use_wrapper=True,
                                        async=False)
    print "*** Splitting by barcodes stage completed ***"
    # Second: merge across batches
    barcoded_fastqs = collect_fastqs(barcoded_fastqs_dir,"*.B*.*.fastq")
    # Merge (concat) fastqs into single pairs per barcode
    final_fastqs_dir = os.path.join(icell8_dir,"fastqs")
    mkdir(final_fastqs_dir)
    merge_fastqs = MergeFastqs("Merge Fastqs",
                               barcoded_fastqs,
                               unassigned_fastqs,
                               failed_barcode_fastqs,
                               failed_umi_fastqs,
                               final_fastqs_dir,
                               basename)
    merge_fastqs = merge_fastqs.run(sched,
                                    working_dir=icell8_dir,
                                    log_dir=log_dir,
                                    use_wrapper=True,
                                    async=False)
    print "*** Merging barcoded FASTQs stage completed ***"

    # Finish
    sched.wait()
    print "All jobs completed"
    sched.stop()
