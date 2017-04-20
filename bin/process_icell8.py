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
from collections import Iterator
from bcftbx.utils import mkdir
from bcftbx.utils import strip_ext
from bcftbx.utils import AttributeDictionary
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

class FileCollection(Iterator):
    """
    """
    def __init__(self,dirn,pattern):
        self._dirn = os.path.abspath(dirn)
        self._pattern = pattern
        self._files = None
        self._idx = None
    def next(self):
        if self._files is None:
            self._idx = 0
            self._files = collect_fastqs(self._dirn,self._pattern)
        else:
            self._idx += 1
        try:
            return self._files[self._idx]
        except IndexError:
            self._files = None
            self._idx = None
            raise StopIteration

class Pipeline(object):
    """
    Class to define and run a 'pipeline' of 'tasks'

    A 'pipeline' in this case is a set of 'tasks', some of
    which may depend on other tasks in the pipeline.

    Example usage:

    >> p = Pipeline()
    >> t1 = p.add_task(Task1())
    >> t2 = p.add_task(Task2(),dependencies=(t1,))
    >> ...
    >> p.run()

    Tasks will only run when all dependenices have
    completed (or will run immediately if they don't have
    any dependencies).
    """
    def __init__(self,default_runner=None):
        self._pending = []
        self._running = []
        self._finished = []
        self._default_runner = default_runner
    def add_task(self,task,dependencies=(),**kws):
        self._pending.append((task,dependencies,kws))
        print "Adding task '%s'" % task.name()
        if dependencies:
            for dep in dependencies:
                if dep.name() not in [t[0].name() for t in self._pending]:
                    print "-> Adding dependency '%s'" % dep.name()
                    self.add_task(dep)
        return task
    def run(self,sched=None,log_dir=None,scripts_dir=None):
        # Execute the pipeline
        print "PIPELINE: started"
        # Run while there are still pending or running tasks
        update = True
        while self._pending or self._running:
            # Report the current running and pending tasks
            if update:
                if self._running:
                    print "PIPELINE: %d running tasks:" % len(self._running)
                    for t in self._running: print "- %s" % t.name()
                if self._pending:
                    print "PIPELINE: %d pending tasks:" % len(self._pending)
                    for t in self._pending: print "- %s" % t[0].name()
                update = False
            # Check for running tasks that have completed
            running = []
            for task in self._running:
                if task.completed:
                    print "PIPELINE: finished %s" % task.name()
                    self._finished.append(task)
                    update = True
                else:
                    running.append(task)
            self._running = running
            # Check for pending tasks that can start
            pending = []
            for task,dependencies,kws in self._pending:
                run_task = False
                if not dependencies:
                    # No dependencies - start it
                    run_task = True
                else:
                    # Check dependencies
                    run_task = reduce(lambda x,y: x and y.completed,
                                      dependencies,True)
                if run_task:
                    print "PIPELINE: started %s" % task.name()
                    if 'runner' not in kws:
                        kws['runner'] = self._default_runner
                    task.run(sched=sched,
                             log_dir=log_dir,
                             scripts_dir=scripts_dir,
                             **kws)
                    self._running.append(task)
                    update = True
                else:
                    pending.append((task,dependencies,kws))
            self._pending = pending
            # Pause before checking again
            if not update:
                time.sleep(5)
        # Finished
        print "PIPELINE: completed"

class PipelineTask(object):
    """
    Base class defining a 'task' to run as part of a pipeline

    A 'task' wraps one or more external programs which can
    be run concurrently, and which produces a set of outputs.
    Individual programs should be wrapped in instances of the
    'PipelineCommand' class.

    This class should be subclassed to implement the 'setup' and
    'output' methods. The 'add_cmd' method can be used within
    'setup' to add one or 'PipelineCommand' instances.

    For example:

    >>> class LsDir(PipelineTask):
    ...    def setup(self,dirn):
    ...      self.add_cmd(LsCommand(dirn))
    ...    def outputs(self):
    ...      return None

    """
    def __init__(self,name,*args,**kws):
        self._name = str(name)
        self._args = args
        self._kws = kws
        self._commands = []
        self._use_wrapper = False
        self._task_name = "%s.%s" % (self._name.lower().replace(' ','_'),
                                     uuid.uuid4())
        self._completed = False
    @property
    def completed(self):
        return self._completed
    def name(self):
        return self._task_name
    def check_status(self,name,jobs,sched):
        self._completed = True
        print "TASK: %s completed" % name
        ##return check_status(name,jobs,sched)
    def add_cmd(self,pipeline_job):
        self._commands.append(pipeline_job)
    def use_wrapper(self,use_wrapper):
        self._use_wrapper = bool(use_wrapper)
    def run(self,sched=None,runner=None,working_dir=None,log_dir=None,
            scripts_dir=None,wait_for=(),async=True,use_wrapper=None):
        # Do setup
        self.setup(*self._args,**self._kws)
        # Sort out defaults
        if use_wrapper is None:
            use_wrapper = self._use_wrapper
        # Generate commands to run
        cmds = []
        for command in self._commands:
            print "TASK: %s" % command.cmd()
            if use_wrapper:
                script_file = command.make_wrapper_script(scripts_dir=scripts_dir)
                cmd = Command('/bin/bash',script_file)
                print "TASK: wrapper script %s" % script_file
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
                for j,cmd in enumerate(cmds):
                    name = "%s#%s" % (self.name(),j)
                    group.add(cmd,
                              wd=working_dir,
                              name=name,
                              runner=runner,
                              log_dir=log_dir,
                              wait_for=wait_for)
                group.close()
                callback_name = group.name
                callback_function = self.check_status
            else:
                # Run a single job
                cmd = cmds[0]
                name = self.name()
                job = sched.submit(cmd,
                                   wd=working_dir,
                                   name=name,
                                   runner=runner,
                                   log_dir=log_dir,
                                   wait_for=wait_for)
                callback_name = job.name
                callback_function = self.check_status
            # If asynchronous then setup callback and
            # return immediately
            if async:
                sched.callback("%s" % self._name,
                               callback_function,
                               wait_for=(callback_name,))
            else:
                # Wait for job or group to complete before returning
                sched.wait_for((callback_name,))
                self._completed = True
                if use_group:
                    exit_code = max([j.exit_code for j in group.jobs])
                else:
                    exit_code = job.exit_code
                if exit_code != 0:
                    logging.warning("%s: failed (exit code %s)" %
                                    (self._name,exit_code))
        else:
            # Run each stage locally
            for cmd in cmds:
                cmd.run_subprocess(working_dir=working_dir)
            self._completed = True
        return self
    def setup(self,*args,**kws):
        raise NotImplementedError("Subclass must implement 'setup' method")
    def output(self):
        raise NotImplementedError("Subclass must implement 'output' method")

class PipelineCommand(object):
    """
    Base class for constructing program command lines

    This class should be subclassed to implement the 'cmd'
    methods, which should return a 'Command' instance.
    For example, to wrap the 'ls' command:

    >>> class LsCommand(PipelineCommand):
    ...    def __init__(self,dirn):
    ...      PipelineCommand.__init__(self,name)
    ...      self._dirn = dirn
    ...    def cmd(self):
    ...      return Command('ls',self._dirn)

    """
    def __init__(self,name):
        self._name = str(name)
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
    def cmd(self):
        # Build the command
        # Must be implemented by the subclass and return a
        # Command instance
        raise NotImplementedError("Subclass must implement 'cmd' method")

######################################################################
# ICell8 pipeline commands
######################################################################

class ICell8Statistics(PipelineCommand):
    """
    Build command to run the 'icell8_stats.py' utility
    """
    def __init__(self,name,fastqs,stats_file,well_list=None,
                 suffix=None,append=False):
        """
        Create new ICell8Statistics instance

        Arguments:
          name (str): description of the command
          fastqs (list): list of FASTQ file names
          stats_file (str): path to output file
          well_list (str): path to 'well list' file (optional)
          suffix (str): suffix to append to columns with read and
            UMI counts (optional)
          append (bool): if True then append columns to existing
            output file (by default creates new output file)
        """
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
    Build command to run the 'split_icell8_fastqs.py' utility
    """
    def __init__(self,name,fastq_pair,filter_dir,well_list=None,
                 basename=None,mode='none',filter=False):
        """
        Create a new SplitAndFilterFastqPair instance

        Arguments:
          name (str): description of the command
          fastq_pair (list): R1/R2 FASTQ file pair
          filter_dir (str): destination directory to
            write output files to
          well_list (str): 'well list' file to use
            (optional)
          basename (str): basename to use for output
            FASTQ files (optional)
          mode (str): mode to run the utility in
          filter (bool): if True then also do barcode
            and UMI quality filtering (no filtering is
            performed by default)
        """
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
    Build command to run the 'batch_fastqs.py' utility
    """
    def __init__(self,name,fastqs,batch_dir,basename,
                 batch_size=None,merge=False):
        """
        Create a new BatchFastqs instance

        Arguments:
          name (str): description of the command
          fastqs (list): list of input FASTQ files
          batch_dir (str): destination directory to
            write output files to
          batch_size (int): number of reads per output
            FASTQ (in batch mode) (optional)
          merge (bool): if True then merge all reads
            into a single FASTQ file (default is to
            split into multiple FASTQs)
        """
        PipelineCommand.__init__(self,name)
        self._fastqs = fastqs
        self._batch_dir = os.path.abspath(batch_dir)
        self._basename = basename
        self._batch_size = batch_size
        self._merge = merge
    def cmd(self):
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
    Build command to run 'cutadapt' with ICell8 settings
    """
    def __init__(self,name,fastq_pair,trim_dir):
        """
        Create a new TrimFastqPair instance

        Arguments:
          name (str): description of the command
          fastq_pair (list): R1/R1 FASTQ file pair
          trim_dir (str): destination directory to
            write output files to
        """
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
    Build command to run 'icell8_contaminantion_filter.py' utility
    """
    def __init__(self,name,fastq_pair,filter_dir,
                 mammalian_conf,contaminants_conf,
                 aligner=None,threads=None):
        """
        Create a new TrimFastqPair instance

        Arguments:
          name (str): description of the command
          fastq_pair (list): R1/R1 FASTQ file pair
          filter_dir (str): destination directory to
            write output files to
          mammalian_conf (str): path to FastqScreen
            .conf file with mammalian genome indexes
          contaminants_conf (str): path FastqScreen
            .conf file with contaminant genome indexes
          aligner (str): explicitly specify name of
            aligner to use with FastqScreen (e.g.
            'bowtie2') (optional)
          threads (int): explicitly specify number of
            threads to run FastqScreen using
            (optional)
        """
        PipelineCommand.__init__(self,name)
        self._fastq_pair = fastq_pair
        self._filter_dir = os.path.abspath(filter_dir)
        self._mammalian_conf = os.path.abspath(mammalian_conf)
        self._contaminants_conf = os.path.abspath(contaminants_conf)
        self._aligner = aligner
        self._threads = threads
    def cmd(self):
        # Build the command
        cmd = Command(
            'icell8_contamination_filter.py',
            '-p',self._mammalian_conf,
            '-c',self._contaminants_conf,
            '-o',self._filter_dir)
        if self._threads:
            cmd.add_args('-n',self._threads)
        if self._aligner is not None:
            cmd.add_args('-a',self._aligner)
        cmd.add_args(*self._fastq_pair)
        return cmd

######################################################################
# ICell8 pipeline tasks
######################################################################

class GetICell8Stats(PipelineTask):
    """
    """
    def setup(self,*args,**kws):
        self.add_cmd(ICell8Statistics(self._name,
                                      *args,
                                      **kws))
        self.use_wrapper(True)
    def output(self):
        stats_file = self._args[1]
        return stats_file

class SplitFastqsIntoBatches(PipelineTask):
    """
    """
    def setup(self,*args,**kws):
        mkdir(self._args[1])
        self.add_cmd(BatchFastqs(self._name,
                                 *args,
                                 **kws))
    def output(self):
        out_dir = self._args[1]
        return FileCollection(out_dir,"*.B*.r*.fastq")

class FilterICell8Fastqs(PipelineTask):
    """
    """
    def setup(self,fastqs,filter_dir,well_list=None,
              mode='none',filter=False):
        mkdir(filter_dir)
        fastq_pairs = pair_fastqs(fastqs)[0]
        for fastq_pair in fastq_pairs:
            basename = os.path.basename(fastq_pair[0])[:-len(".r1.fastq")]
            self.add_cmd(SplitAndFilterFastqPair(self._name,
                                                 fastq_pair,
                                                 filter_dir,
                                                 well_list=well_list,
                                                 basename=basename,
                                                 mode=mode,
                                                 filter=filter))
        self.use_wrapper(True)
    def output(self):
        out_dir = self._args[1]
        return AttributeDictionary(
            assigned=FileCollection(out_dir,"*.B*.filtered.r*.fastq"),
            unassigned=FileCollection(out_dir,"*.unassigned.r*.fastq"),
            failed_barcodes=FileCollection(out_dir,"*.failed_barcode.r*.fastq"),
            failed_umis=FileCollection(out_dir,"*.failed_umi.r*.fastq")
        )

class TrimReads(PipelineTask):
    """
    """
    def setup(self,fastqs,trim_dir):
        mkdir(trim_dir)
        fastq_pairs = pair_fastqs(fastqs)[0]
        for fastq_pair in fastq_pairs:
            self.add_cmd(TrimFastqPair(self._name,
                                       fastq_pair,
                                       trim_dir))
    def output(self):
        out_dir = self._args[1]
        return FileCollection(out_dir,"*.trimmed.fastq")

class FilterContaminatedReads(PipelineTask):
    """
    """
    def setup(self,fastqs,filter_dir,mammalian_conf,
                 contaminants_conf,aligner=None,threads=None):
        mkdir(filter_dir)
        fastq_pairs = pair_fastqs(fastqs)[0]
        for fastq_pair in fastq_pairs:
            self.add_cmd(ContaminantFilterFastqPair(self._name,
                                                    fastq_pair,
                                                    filter_dir,
                                                    mammalian_conf,
                                                    contaminants_conf,
                                                    aligner=aligner,
                                                    threads=threads))
        self.use_wrapper(True)
    def output(self):
        out_dir = self._args[1]
        return FileCollection(out_dir,"*.trimmed.filtered.fastq")

class SplitByBarcodes(PipelineTask):
    """
    """
    def setup(self,fastqs,barcodes_dir):
        mkdir(barcodes_dir)
        fastq_pairs = pair_fastqs(fastqs)[0]
        for fastq_pair in fastq_pairs:
            basename = os.path.basename(fastq_pair[0])[:-len(".r1.fastq")+1]
            self.add_cmd(SplitAndFilterFastqPair(self._name,
                                                 fastq_pair,
                                                 barcodes_dir,
                                                 basename=basename,
                                                 mode="barcodes"))
        self.use_wrapper(True)
    def output(self):
        out_dir = self._args[1]
        return FileCollection(out_dir,"*.r*.fastq")

class MergeFastqs(PipelineTask):
    """
    """
    def setup(self,fastqs,unassigned_fastqs,
              failed_barcode_fastqs,failed_umi_fastqs,
              merge_dir,basename,batch_size=25):
        # Make the output directory
        mkdir(merge_dir)
        # Extract the barcodes from the fastq names
        barcodes = set()
        for fq in fastqs:
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
            self.add_cmd(SplitAndFilterFastqPair(self._name,
                                                 fastq_pairs,
                                                 merge_dir,
                                                 basename=basename,
                                                 mode="barcodes"))
        # Handle unassigned and failed quality reads
        for name,fqs in (('unassigned',unassigned_fastqs),
                         ('failed_barcodes',failed_barcode_fastqs),
                         ('failed_umis',failed_umi_fastqs)):
            if not fqs:
                continue
            self.add_cmd(BatchFastqs(self._name,fqs,
                                     merge_dir,
                                     "%s.%s" % (basename,name),
                                     merge=True))
        # Force jobs to use a wrapper script
        self.use_wrapper(True)
    def output(self):
        out_dir = self._args[4]
        return AttributeDictionary(
            assigned=FileCollection(out_dir,"*.[ACGT]*.r*.fastq"),
            unassigned=FileCollection(out_dir,"*.unassigned.r*.fastq"),
            failed_barcodes=FileCollection(out_dir,"*.failed_barcodes.r*.fastq"),
            failed_umis=FileCollection(out_dir,"*.failed_umis.r*.fastq"),
        )

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
    p.add_argument("-u","--unaligned",
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
        job_start="SCHEDULER: Started  #%(job_number)d: %(job_name)s:\n-- %(command)s",
        job_end=  "SCHEDULER: Finished #%(job_number)d: %(job_name)s"
    )
    sched_reporter = SchedulerReporter()
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

    # Set up a pipeline
    ppl = Pipeline()

    # Initial stats
    initial_stats = GetICell8Stats("Initial statistics",
                                   fastqs,
                                   os.path.join(stats_dir,"icell8_stats.tsv"),
                                   well_list)
    ##ppl.add_task(initial_stats)

    # Split fastqs into batches
    batch_dir = os.path.join(icell8_dir,"_fastqs.batched")
    batch_fastqs = SplitFastqsIntoBatches("Batch Fastqs",fastqs,
                                          batch_dir,basename,
                                          batch_size=args.batch_size)
    ##ppl.add_task(batch_fastqs)

    # Setup the quality filter jobs as a group
    filter_dir = os.path.join(icell8_dir,"_fastqs.quality_filter")
    quality_filter = FilterICell8Fastqs("Quality filter Fastqs",
                                        batch_fastqs.output(),
                                        filter_dir,
                                        well_list=well_list,
                                        mode='none',
                                        filter=True)
    ppl.add_task(quality_filter,dependencies=(batch_fastqs,))
    
    # Post quality filter stats
    filter_stats = GetICell8Stats("Post-quality filter statistics",
                                  quality_filter.output().assigned,
                                  initial_stats.output(),
                                  suffix="_quality_filtered",
                                  append=True)
    ppl.add_task(filter_stats,dependencies=(initial_stats,quality_filter))

    # Set up the cutadapt jobs as a group
    trim_dir = os.path.join(icell8_dir,"_fastqs.trim_reads")
    trim_reads = TrimReads("Read trimming",
                           quality_filter.output().assigned,
                           trim_dir)
    ppl.add_task(trim_reads,dependencies=(quality_filter,))

    # Post read trimming stats
    trim_stats = GetICell8Stats("Post-trimming statistics",
                                trim_reads.output(),
                                initial_stats.output(),
                                suffix="_trimmed",
                                append=True)
    ppl.add_task(trim_stats,dependencies=(trim_reads,))

    # Set up the contaminant filter jobs as a group
    contaminant_filter_dir = os.path.join(icell8_dir,
                                          "_fastqs.contaminant_filter")
    contaminant_filter = FilterContaminatedReads("Contaminant filtering",
                                                 trim_reads.output(),
                                                 contaminant_filter_dir,
                                                 args.preferred_conf,
                                                 args.contaminants_conf,
                                                 aligner=args.aligner,
                                                 threads=args.threads)
    ppl.add_task(contaminant_filter,dependencies=(trim_reads,),
                 runner=runners['contaminant_filter'])

    # Post contaminant filter stats
    final_stats = GetICell8Stats("Post-contaminant filter statistics",
                                 contaminant_filter.output(),
                                 initial_stats.output(),
                                 suffix="_contaminant_filtered",
                                 append=True)
    ppl.add_task(final_stats,dependencies=(contaminant_filter,trim_stats))

    # Rebatch reads by barcode
    # First: split each batch by barcode
    barcoded_fastqs_dir = os.path.join(icell8_dir,"_fastqs.barcodes")
    split_barcodes = SplitByBarcodes("Split by barcodes",
                                     contaminant_filter.output(),
                                     barcoded_fastqs_dir)
    ppl.add_task(split_barcodes,dependencies=(contaminant_filter,))
    # Merge (concat) fastqs into single pairs per barcode
    final_fastqs_dir = os.path.join(icell8_dir,"fastqs")
    merge_fastqs = MergeFastqs("Merge Fastqs",
                               split_barcodes.output(),
                               quality_filter.output().unassigned,
                               quality_filter.output().failed_barcodes,
                               quality_filter.output().failed_umis,
                               final_fastqs_dir,
                               basename)
    ppl.add_task(merge_fastqs,dependencies=(split_barcodes,))

    # Final stats for verification
    final_barcode_stats = GetICell8Stats(
        "Post-barcode splitting and merging statistics",
        merge_fastqs.output().assigned,
        initial_stats.output(),
        suffix="_final",
        append=True)
    ppl.add_task(final_barcode_stats,dependencies=(merge_fastqs,))

    # Execute the pipeline
    ppl.run(sched=sched,log_dir=log_dir,scripts_dir=scripts_dir)

    # Finish
    print "All jobs completed"
    sched.stop()
