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
import inspect
import traceback
from collections import Iterator
from cStringIO import StringIO
from bcftbx.utils import mkdir
from bcftbx.utils import strip_ext
from bcftbx.utils import AttributeDictionary
from bcftbx.IlluminaData import IlluminaData
from bcftbx.IlluminaData import IlluminaDataError
from bcftbx.FASTQFile import FastqIterator
from bcftbx.JobRunner import fetch_runner
from bcftbx.TabFile import TabFile
from auto_process_ngs.applications import Command
from auto_process_ngs.simple_scheduler import SimpleScheduler
from auto_process_ngs.simple_scheduler import SchedulerReporter
from auto_process_ngs.simple_scheduler import SchedulerGroup
from auto_process_ngs.fastq_utils import pair_fastqs
from auto_process_ngs.fastq_utils import get_read_number
from auto_process_ngs.utils import BaseFastqAttrs
from auto_process_ngs.utils import AnalysisFastq
from auto_process_ngs.utils import AnalysisProject
from auto_process_ngs.qc.illumina_qc import check_qc_outputs
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
    Class to return set of files based on glob pattern
    """
    def __init__(self,dirn,pattern):
        self._dirn = os.path.abspath(dirn)
        self._pattern = pattern
        self._files = None
        self._idx = None
    def __len__(self):
        if self._files is None:
            self._files = collect_fastqs(self._dirn,self._pattern)
            self._idx = -1
        return len(self._files)
    def next(self):
        if self._files is None:
            self._files = collect_fastqs(self._dirn,self._pattern)
        if self._idx is None:
            self._idx = 0
        else:
            self._idx += 1
        try:
            return self._files[self._idx]
        except IndexError:
            self._files = None
            self._idx = None
            raise StopIteration

# Capture stdout from a function call - see
# http://stackoverflow.com/a/16571630/579925
class Capturing(list):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = StringIO()
        return self
    def __exit__(self,*args):
        self.extend(self._stringio.getvalue().splitlines())
        del self._stringio    # free up some memory
        sys.stdout = self._stdout

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
    def __init__(self,name="PIPELINE",default_runner=None,
                 working_dir=None):
        self._name = str(name)
        self._pending = []
        self._running = []
        self._finished = []
        self._default_runner = default_runner
        if working_dir is None:
            working_dir = os.getcwd()
        self._working_dir = os.path.abspath(working_dir)
    def add_task(self,task,dependencies=(),**kws):
        self._pending.append((task,dependencies,kws))
        self.report("Adding task '%s'" % task.name())
        if dependencies:
            for dep in dependencies:
                if dep.name() not in [t[0].name() for t in self._pending]:
                    self.report("-> Adding dependency '%s'" % dep.name())
                    self.add_task(dep)
        return task
    def report(self,s):
        print "%s [%s] %s" % (time.strftime("%Y-%m-%d %H:%M:%S"),
                              self._name,s)
    def run(self,sched=None,log_dir=None,scripts_dir=None):
        # Execute the pipeline
        self.report("Started")
        # Run while there are still pending or running tasks
        update = True
        while self._pending or self._running:
            # Report the current running and pending tasks
            if update:
                if self._running:
                    self.report("%d running tasks:"
                                % len(self._running))
                    for t in self._running:
                        self.report("- %s" % t.name())
                if self._pending:
                    self.report("%d pending tasks:"
                                % len(self._pending))
                    for t in self._pending:
                        self.report("- %s" % t[0].name())
                update = False
            # Check for running tasks that have completed
            running = []
            failed = []
            for task in self._running:
                if task.completed:
                    self.report("finished %s"
                                % task.name())
                    self._finished.append(task)
                    update = True
                    # Check if task failed
                    if task.exit_code != 0:
                        failed.append(task)
                else:
                    running.append(task)
            self._running = running
            # Check for finished tasks that have failed
            if failed:
                self.report("Following tasks failed:")
                for task in failed:
                    self.report("- %s" % task.name())
                self.report("Terminating prematurely")
                return 1
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
                    self.report("started %s" % task.name())
                    if 'runner' not in kws:
                        kws['runner'] = self._default_runner
                    if 'working_dir' not in kws:
                        kws['working_dir'] = self._working_dir
                    try:
                        task.run(sched=sched,
                                 log_dir=log_dir,
                                 scripts_dir=scripts_dir,
                                 **kws)
                    except Exception as ex:
                        self.report("Failed to start task '%s': %s" %
                                    (task.name(),ex))
                        logging.critical("Failed to start task '%s': %s" %
                                         (task.name(),ex))
                        return 1
                    self._running.append(task)
                    update = True
                else:
                    pending.append((task,dependencies,kws))
            self._pending = pending
            # Pause before checking again
            if not update:
                time.sleep(5)
        # Finished
        self.report("Completed")
        return 0

class PipelineTask(object):
    """
    Base class defining a 'task' to run as part of a pipeline

    A 'task' wraps one or more external programs which can
    be run concurrently, and which produces a set of outputs.
    Individual programs should be wrapped in instances of the
    'PipelineCommand' class.

    This class should be subclassed to implement the 'setup',
    'finish' and 'output' methods.

    The 'add_cmd' method can be used within 'setup' to add one
    or 'PipelineCommand' instances.

    For example:

    >>> class LsDir(PipelineTask):
    ...    def init(self,dirn):
    ...      pass
    ...    def setup(self):
    ...      self.add_cmd(LsCommand(self.args.dirn))
    ...    def outputs(self):
    ...      return None

    """
    def __init__(self,name,*args,**kws):
        self._name = str(name)
        self._args = args
        self._kws = kws
        self._commands = []
        self._task_name = "%s.%s" % (self._name.lower().replace(' ','_'),
                                     uuid.uuid4())
        self._completed = False
        self._exit_code = 0
        # Deal with subclass arguments
        self._callargs = inspect.getcallargs(self.init,*args,**kws)
        try:
            del(self._callargs['self'])
        except KeyError:
            pass
    @property
    def args(self):
        return AttributeDictionary(**self._callargs)
    @property
    def completed(self):
        return self._completed
    @property
    def exit_code(self):
        if not self.completed:
            return None
        else:
            return self._exit_code
    def name(self):
        return self._task_name
    def report(self,s):
        print "%s [Task: %s] %s" % (time.strftime("%Y-%m-%d %H:%M:%S"),
                                    self._name,s)
    def invoke(self,f,args=None,kws=None):
        try:
            with Capturing() as output:
                if args is None:
                    f()
                else:
                    f(*args,**kws)
            self.report("done '%s'" % f.__name__)
            for line in output:
                self.report("%s STDOUT: %s" % (f.__name__,line))
        except NotImplementedError:
            pass
        except Exception as ex:
            self.report("exception invoking '%s': %s" %
                        (f.__name__,ex))
            traceback.print_exc(ex)
            self._exit_code += 1
    def task_completed(self,name,jobs,sched):
        # Callback method invoked when scheduled
        # jobs in the task finish
        # Arguments:
        # name (str): name for the callback
        # jobs (list): list of SchedulerJob instances
        # sched (SimpleScheduler): scheduler instance
        for job in jobs:
            try:
                if job.exit_code != 0:
                    self._exit_code += 1
            except AttributeError:
                # Assume it's a group
                for j in job.jobs:
                    if j.exit_code != 0:
                        self._exit_code += 1
        if self._exit_code != 0:
            logging.critical("%s failed: exit code %s"
                             % (self._name,self._exit_code))
        else:
            # Execute 'finish', if implemented
            self.invoke(self.finish)
        # Flag job as completed
        self._completed = True
        self.report("%s completed" % name)
    def add_cmd(self,pipeline_job):
        self._commands.append(pipeline_job)
    def run(self,sched=None,runner=None,working_dir=None,log_dir=None,
            scripts_dir=None,wait_for=(),async=True):
        # Do init and setup
        self.invoke(self.init,self._args,self._kws)
        self.invoke(self.setup)
        # Generate commands to run
        cmds = []
        for command in self._commands:
            self.report("%s" % command.cmd())
            script_file = command.make_wrapper_script(scripts_dir=scripts_dir)
            cmd = Command('/bin/bash',script_file)
            self.report("wrapper script %s" % script_file)
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
                callback_function = self.task_completed
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
                callback_function = self.task_completed
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
        else:
            # Run each stage locally
            for cmd in cmds:
                cmd.run_subprocess(working_dir=working_dir)
            self.invoke(self.finish)
            self._completed = True
        return self
    def init(self,*args,**kws):
        pass
    def setup(self):
        raise NotImplementedError("Subclass must implement 'setup' method")
    def finish(self):
        raise NotImplementedError("Subclass must implement 'finish' method")
    def output(self):
        raise NotImplementedError("Subclass must implement 'output' method")

class PipelineCommand(object):
    """
    Base class for constructing program command lines

    This class should be subclassed to implement the 'init'
    and 'cmd' methods.

    The 'init' method should do any preprocessing and
    caching of arguments to be used in the 'cmd' method;
    the 'cmd' method should use these to construct and
    return a 'Command' instance.

    For example, to wrap the 'ls' command:

    >>> class LsCommand(PipelineCommand):
    ...    def init(self,dirn):
    ...      self._dirn = dirn
    ...    def cmd(self):
    ...      return Command('ls',self._dirn)

    """
    def __init__(self,*args,**kws):
        # Set internal name
        self._name = self.__class__.__name__
        # Invoke the 'init' method
        self.init(*args,**kws)
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
    def init(self):
        # Initialise and store parameters
        # Must be implemented by the subclass
        raise NotImplementedError("Subclass must implement 'init' method")
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
    def init(self,fastqs,stats_file,well_list=None,
             suffix=None,append=False,nprocs=1):
        """
        Create new ICell8Statistics instance

        Arguments:
          fastqs (list): list of FASTQ file names
          stats_file (str): path to output file
          well_list (str): path to 'well list' file (optional)
          suffix (str): suffix to append to columns with read and
            UMI counts (optional)
          append (bool): if True then append columns to existing
            output file (by default creates new output file)
          nprocs (int): number of cores available for stats
            (default: 1)
        """
        self._fastqs = fastqs
        self._stats_file = os.path.abspath(stats_file)
        self._well_list = well_list
        if self._well_list is not None:
            self._well_list = os.path.abspath(self._well_list)
        self._append = append
        self._suffix = suffix
        self._nprocs = nprocs
    def cmd(self):
        # Build command
        cmd = Command('icell8_stats.py',
                      '-f',self._stats_file,
                      '-n',self._nprocs)
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
    def init(self,fastq_pair,out_dir,well_list=None,
             basename=None,mode='none',
             discard_unknown_barcodes=False,
             quality_filter=False):
        """
        Create a new SplitAndFilterFastqPair instance

        Arguments:
          fastq_pair (list): R1/R2 FASTQ file pair
          out_dir (str): destination directory to
            write output files to
          well_list (str): 'well list' file to use
            (optional)
          basename (str): basename to use for output
            FASTQ files (optional)
          mode (str): mode to run the utility in
          discard_unknown_barcodes (bool): if True
            then discard read pairs where the barcode
            doesn't match one of those in the well
            list file (nb well list file must also
            be supplied in this case) (all reads are
            kept by default)
          quality_filter (bool): if True then also
            do filtering based on barcode- and
            UMI-quality (no filtering is performed
            by default)
        """
        self._fastq_pair = fastq_pair
        self._out_dir = os.path.abspath(out_dir)
        self._well_list = well_list
        self._basename = basename
        self._mode = mode
        self._discard_unknown_barcodes = discard_unknown_barcodes
        self._quality_filter = quality_filter
        if self._well_list is not None:
            self._well_list = os.path.abspath(self._well_list)
    def cmd(self):
        cmd = Command('split_icell8_fastqs.py',
                      '-o',self._out_dir,
                      '-b',self._basename)
        if self._well_list:
            cmd.add_args('-w',self._well_list)
        if self._mode:
            cmd.add_args('-m',self._mode)
        if self._discard_unknown_barcodes:
            cmd.add_args('--discard-unknown-barcodes')
        if self._quality_filter:
            cmd.add_args('--quality-filter')
        cmd.add_args(*self._fastq_pair)
        return cmd

class BatchFastqs(PipelineCommand):
    """
    Split reads from Fastqs into batches using (z)cat/split

    Given a list of Fastq files, combines them and then
    splits into batches of a specified number of reads by
    running a combination of '(z)cat' and 'split' commands.

    Fastqs can be gzipped, but must have the same read number
    (i.e. R1 or R2).
    """
    def init(self,fastqs,batch_dir,basename,
             batch_size=DEFAULT_BATCH_SIZE):
        """
        Create a new BatchFastqs instance

        Arguments:
          fastqs (list): list of input Fastq files
          batch_dir (str): destination directory to
            write output files to
          basename (str): basename for output Fastqs
          batch_size (int): number of reads per output
            FASTQ (in batch mode) (optional)
        """
        # Store inputs
        self._fastqs = fastqs
        self._batch_dir = os.path.abspath(batch_dir)
        self._basename = basename
        self._batch_size = batch_size
        # Determine if fastqs are gzipped
        first_fastq = self._fastqs[0]
        self._gzipped = first_fastq.endswith('.gz')
        # Determine read number
        self._read_number = get_read_number(first_fastq)
    def cmd(self):
        # Constructs command line of the form:
        # zcat FASTQ | split -l BATCH_SIZE*4 -d -a 3 \
        #   --additional-suffix=.r1.fastq - BASENAME.B
        if self._gzipped:
            cmd = Command('zcat')
        else:
            cmd = Command('cat')
        cmd.add_args(*self._fastqs)
        cmd.add_args('|',
                     'split',
                     '-l',self._batch_size*4,
                     '-d',
                     '-a',3,
                     '--additional-suffix=.r%d.fastq' %
                     self._read_number,
                     '-',
                     os.path.join(self._batch_dir,
                                  '%s.B' % self._basename))
        return cmd

class ConcatFastqs(PipelineCommand):
    """
    Concatenate reads from multiple Fastqs into a single file

    Given a list of Fastq files, combines them into a single
    Fastq using the 'cat' utility.

    FASTQs cannot be gzipped, and must all be same read number
    (i.e. R1 or R2).
    """
    def init(self,fastqs,concat_dir,fastq_out):
        """
        Create a new ConcatFastqs instance

        Arguments:
          fastqs (list): list of input FASTQ files
          concat_dir (str): destination directory to
            write output file to
          fastq_out (str): name of output Fastq file
        """
        # Store inputs
        self._fastqs = fastqs
        self._concat_dir = os.path.abspath(concat_dir)
        self._fastq_out = fastq_out
    def cmd(self):
        cmd = Command('cat')
        cmd.add_args(*self._fastqs)
        cmd.add_args('>',
                     os.path.join(self._concat_dir,
                                  self._fastq_out))
        return cmd

class TrimFastqPair(PipelineCommand):
    """
    Build command to run 'cutadapt' with ICell8 settings
    """
    def init(self,fastq_pair,trim_dir):
        """
        Create a new TrimFastqPair instance

        Arguments:
          fastq_pair (list): R1/R1 FASTQ file pair
          trim_dir (str): destination directory to
            write output files to
        """
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

class FilterPolyGReads(PipelineCommand):
    """
    'cutadapt' to fetch reads with poly-G regions
    """
    def init(self,fastq_pair,out_dir):
        """
        Create a new GetPolyGReads instance

        Arguments:
          fastq_pair (list): R1/R1 FASTQ file pair
          out_dir (str): destination directory to
            write output files to
        """
        self._fastq_pair = fastq_pair
        self._out_dir = os.path.abspath(out_dir)
    def cmd(self):
        # Generate output file pair names
        fastq_pair_out = [os.path.join(self._out_dir,
                                       strip_ext(os.path.basename(fq),'.fastq')
                                       + '.poly_g.fastq')
                          for fq in self._fastq_pair]
        # Build command
        cmd = Command(
            'cutadapt',
            '-a','GGGGGGG',
            '--discard-untrimmed')
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
    def init(self,fastq_pair,filter_dir,
             mammalian_conf,contaminants_conf,
             aligner=None,threads=None):
        """
        Create a new TrimFastqPair instance

        Arguments:
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
            '-m',self._mammalian_conf,
            '-c',self._contaminants_conf,
            '-o',self._filter_dir)
        if self._threads:
            cmd.add_args('-n',self._threads)
        if self._aligner is not None:
            cmd.add_args('-a',self._aligner)
        cmd.add_args(*self._fastq_pair)
        return cmd

class IlluminaQC(PipelineCommand):
    """
    Run the 'illumina_qc.sh' script on one or more Fastqs
    """
    def init(self,fastqs,nthreads=1,working_dir=None):
        """
        Set up parameters
        """
        self.fastqs = fastqs
        self.nthreads = nthreads
        if working_dir is None:
            working_dir = os.getcwd()
        self.working_dir = os.path.abspath(working_dir)
    def cmd(self):
        """
        Build the command
        """
        cmd = Command('cd',self.working_dir)
        for fastq in self.fastqs:
            cmd.add_args('&&','illumina_qc.sh',fastq)
            if self.nthreads > 1:
                cmd.add_args('--threads',self.nthreads)
        return cmd

class RemoveDirectory(PipelineCommand):
    """
    Command to remove a directory and its contents
    """
    def init(self,dirn):
        """
        Set up parameters
        """
        self.dirn = os.path.abspath(dirn)
    def cmd(self):
        """
        Build the command
        """
        # Does "rm -f DIRN/* && rmdir DIRN"
        return Command(
            "rm","-f","%s" % os.path.join(self.dirn,'*'),
            "&&",
            "rmdir","%s" % self.dirn)

######################################################################
# ICell8 pipeline tasks
######################################################################

class GetICell8Stats(PipelineTask):
    """
    """
    def init(self,fastqs,stats_file,well_list=None,
             suffix=None,append=False,nprocs=1):
        pass
    def setup(self):
        self.add_cmd(ICell8Statistics(self.args.fastqs,
                                      self.args.stats_file,
                                      well_list=self.args.well_list,
                                      suffix=self.args.suffix,
                                      append=self.args.append,
                                      nprocs=self.args.nprocs))
    def output(self):
        return self.args.stats_file

class GetICell8PolyGStats(GetICell8Stats):
    """
    """
    def finish(self):
        # Method invoked once commands have run
        # Adds another column to the stats file
        # with the percentage of 'unfiltered' reads
        # with poly-G regions
        print "Add number of reads with poly-G regions as percentage"
        stats_file = self.args.stats_file
        stats = TabFile(stats_file,first_line_is_header=True)
        # Add and populate the new column
        stats.appendColumn("%reads_poly_g")
        for line in stats:
            try:
                perc_poly_g = (float(line['Nreads_poly_g'])/
                               float(line['Nreads_filtered'])*100.0)
            except ValueError:
                perc_poly_g = 0.0
            line["%reads_poly_g"] = ("%.2f" % perc_poly_g)
        # Write out the updated stats file
        stats.write(stats_file,include_header=True)

class SplitFastqsIntoBatches(PipelineTask):
    """
    """
    def init(self,fastqs,batch_dir,basename,
             batch_size=DEFAULT_BATCH_SIZE):
        pass
    def setup(self):
        # Make the output directory
        mkdir(self.args.batch_dir)
        # Set up the commands
        fastq_pairs = pair_fastqs(self.args.fastqs)[0]
        fastqs_r1 = [p[0] for p in fastq_pairs]
        self.add_cmd(BatchFastqs(fastqs_r1,
                                 self.args.batch_dir,
                                 self.args.basename,
                                 batch_size=self.args.batch_size))
        fastqs_r2 = [p[1] for p in fastq_pairs]
        self.add_cmd(BatchFastqs(fastqs_r2,
                                 self.args.batch_dir,
                                 self.args.basename,
                                 batch_size=self.args.batch_size))
    def output(self):
        out_dir = self.args.batch_dir
        return FileCollection(out_dir,"*.B*.r*.fastq")

class FilterICell8Fastqs(PipelineTask):
    """
    """
    def init(self,fastqs,filter_dir,well_list=None,
              mode='none',discard_unknown_barcodes=False,
              quality_filter=False):
        pass
    def setup(self):
        mkdir(self.args.filter_dir)
        fastq_pairs = pair_fastqs(self.args.fastqs)[0]
        for fastq_pair in fastq_pairs:
            basename = os.path.basename(fastq_pair[0])[:-len(".r1.fastq")]
            self.add_cmd(SplitAndFilterFastqPair(
                fastq_pair,
                self.args.filter_dir,
                well_list=self.args.well_list,
                basename=basename,
                mode=self.args.mode,
                discard_unknown_barcodes=self.args.discard_unknown_barcodes,
                quality_filter=self.args.quality_filter))
    def output(self):
        out_dir = self.args.filter_dir
        return AttributeDictionary(
            assigned=FileCollection(out_dir,"*.B*.filtered.r*.fastq"),
            unassigned=FileCollection(out_dir,"*.unassigned.r*.fastq"),
            failed_barcodes=FileCollection(out_dir,"*.failed_barcode.r*.fastq"),
            failed_umis=FileCollection(out_dir,"*.failed_umi.r*.fastq")
        )

class TrimReads(PipelineTask):
    """
    """
    def init(self,fastqs,trim_dir):
        pass
    def setup(self):
        mkdir(self.args.trim_dir)
        fastq_pairs = pair_fastqs(self.args.fastqs)[0]
        for fastq_pair in fastq_pairs:
            self.add_cmd(TrimFastqPair(fastq_pair,
                                       self.args.trim_dir))
    def output(self):
        out_dir = self.args.trim_dir
        return FileCollection(out_dir,"*.trimmed.fastq")

class GetReadsWithPolyGRegions(PipelineTask):
    """
    """
    def init(self,fastqs,poly_g_regions_dir):
        pass
    def setup(self):
        mkdir(self.args.poly_g_regions_dir)
        fastq_pairs = pair_fastqs(self.args.fastqs)[0]
        for fastq_pair in fastq_pairs:
            self.add_cmd(FilterPolyGReads(fastq_pair,
                                          self.args.poly_g_regions_dir))
    def output(self):
        out_dir = self.args.poly_g_regions_dir
        return FileCollection(out_dir,"*.poly_g.fastq")

class FilterContaminatedReads(PipelineTask):
    """
    """
    def init(self,fastqs,filter_dir,mammalian_conf,
             contaminants_conf,aligner=None,threads=None):
        pass
    def setup(self):
        mkdir(self.args.filter_dir)
        fastq_pairs = pair_fastqs(self.args.fastqs)[0]
        for fastq_pair in fastq_pairs:
            self.add_cmd(ContaminantFilterFastqPair(
                fastq_pair,
                self.args.filter_dir,
                self.args.mammalian_conf,
                self.args.contaminants_conf,
                aligner=self.args.aligner,
                threads=self.args.threads))
    def output(self):
        out_dir = self.args.filter_dir
        return FileCollection(out_dir,"*.trimmed.filtered.fastq")

class SplitByBarcodes(PipelineTask):
    """
    """
    def init(self,fastqs,barcodes_dir):
        pass
    def setup(self):
        mkdir(self.args.barcodes_dir)
        fastq_pairs = pair_fastqs(self.args.fastqs)[0]
        for fastq_pair in fastq_pairs:
            basename = os.path.basename(fastq_pair[0])[:-len(".r1.fastq")+1]
            self.add_cmd(SplitAndFilterFastqPair(
                fastq_pair,
                self.args.barcodes_dir,
                basename=basename,
                mode="barcodes"))
    def output(self):
        out_dir = self.args.barcodes_dir
        return FileCollection(out_dir,"*.r*.fastq")

class MergeFastqs(PipelineTask):
    """
    """
    def init(self,fastqs,unassigned_fastqs,
             failed_barcode_fastqs,failed_umi_fastqs,
             merge_dir,basename,batch_size=25):
        pass
    def setup(self):
        # Move existing output directory
        if os.path.exists(self.args.merge_dir):
            print "Moving existing dir '%s'" % self.args.merge_dir
            os.rename(self.args.merge_dir,
                      "%s.save" % self.args.merge_dir)
        # Make the output directory
        mkdir(self.args.merge_dir)
        # Extract the barcodes from the fastq names
        barcodes = set()
        for fq in self.args.fastqs:
            barcode = os.path.basename(fq).split('.')[-3]
            barcodes.add(barcode)
        barcodes = sorted(list(barcodes))
        # Group files by barcode
        fastq_groups = dict()
        for barcode in barcodes:
            fqs = filter(lambda fq: (fq.endswith("%s.r1.fastq" % barcode) or
                                     fq.endswith("%s.r2.fastq" % barcode)),
                         self.args.fastqs)
            fastq_groups[barcode] = fqs
        # Group barcodes into batches
        barcode_batches = [barcodes[i:i+self.args.batch_size]
                           for i in xrange(0,len(barcodes),
                                           self.args.batch_size)]
        # Concat fastqs
        for i,barcode_batch in enumerate(barcode_batches):
            batch_name = "barcodes%06d" % i
            print "Barcode batch: %s" % batch_name
            fastq_pairs = []
            for barcode in barcode_batch:
                print "-- %s" % barcode
                fastq_pairs.extend(fastq_groups[barcode])
            self.add_cmd(SplitAndFilterFastqPair(fastq_pairs,
                                                 self.args.merge_dir,
                                                 basename=self.args.basename,
                                                 mode="barcodes"))
        # Handle unassigned and failed quality reads
        for name,fqs in (('unassigned',self.args.unassigned_fastqs),
                         ('failed_barcodes',self.args.failed_barcode_fastqs),
                         ('failed_umis',self.args.failed_umi_fastqs)):
            if not fqs:
                continue
            fastq_pairs = pair_fastqs(fqs)[0]
            fqs_r1 = [p[0] for p in fastq_pairs]
            self.add_cmd(ConcatFastqs(fqs_r1,
                                      self.args.merge_dir,
                                      "%s.%s.r1.fastq" %
                                      (self.args.basename,name)))
            fqs_r2 = [p[1] for p in fastq_pairs]
            self.add_cmd(ConcatFastqs(fqs_r2,
                                      self.args.merge_dir,
                                      "%s.%s.r2.fastq" %
                                      (self.args.basename,name)))
    def output(self):
        out_dir = self.args.merge_dir
        return AttributeDictionary(
            assigned=FileCollection(out_dir,"*.[ACGT]*.r*.fastq"),
            unassigned=FileCollection(out_dir,"*.unassigned.r*.fastq"),
            failed_barcodes=FileCollection(out_dir,"*.failed_barcodes.r*.fastq"),
            failed_umis=FileCollection(out_dir,"*.failed_umis.r*.fastq"),
        )

class RunQC(PipelineTask):
    """
    """
    def init(self,project_dir,nthreads=1,batch_size=25):
        self.qc_dir = os.path.join(self.args.project_dir,'qc')
        self.qc_report = None
    def setup(self):
        # Gather Fastqs
        project = AnalysisProject(
            os.path.basename(self.args.project_dir),
            self.args.project_dir,
            fastq_attrs=ICell8FastqAttrs)
        batch_size = self.args.batch_size
        fastqs = []
        for sample in project.samples:
            fastqs.extend(sample.fastq)
        # Make the output qc directory
        mkdir(self.qc_dir)
        # Set up QC run for batches of fastqs
        while fastqs:
            self.add_cmd(IlluminaQC(fastqs[:batch_size],
                                    nthreads=self.args.nthreads,
                                    working_dir=self.args.project_dir))
            fastqs = fastqs[batch_size:]
    def finish(self):
        # Verify the QC outputs for each fastq
        project = AnalysisProject(
            os.path.basename(self.args.project_dir),
            self.args.project_dir,
            fastq_attrs=ICell8FastqAttrs)
        if not project.verify_qc():
            print "Failed to verify QC"
            self._exit_code = 1
        else:
            # Generate QC report
            print "Generating QC report"
            self.qc_report = project.qc_report(force=True)
            if self.qc_report is None:
                print "Failed to generate QC report"
                self._exit_code = 1
            else:
                print "QC report: %s" % self.qc_report
    def output(self):
        return AttributeDictionary(
            qc_dir=self.qc_dir,
            report_zip=self.qc_report,
        )

class CleanupDirectory(PipelineTask):
    """
    """
    def init(self,dirn):
        pass
    def setup(self):
        if not os.path.isdir(self.args.dirn):
            self.report("No directory '%s'" % self.args.dirn)
        else:
            self.add_cmd(RemoveDirectory(self.args.dirn))

######################################################################
# Classes
######################################################################

class ICell8FastqAttrs(BaseFastqAttrs):
    """
    Extract attributes from Icell8 pipeline Fastq names

    Used as a custom Fastq name attribute parser for output
    Fastqs from the ICell8 pipeline.

    Fastq names for the final results look like e.g.:

    icell8.CCAGTTCAGGA.r1.fastq
    """
    def __init__(self,fastq):
        BaseFastqAttrs.__init__(self,fastq)
        name = self.basename.split('.')
        self.sample_name = '.'.join(name[0:2])
        self.barcode_sequence = name[1]
        self.read_number = int(name[2][1:])
    def __repr__(self):
        return "%s.r%s" % (self.sample_name,
                           self.read_number)

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
    p.add_argument("-m","--mammalian",
                   dest="mammalian_conf",
                   help="fastq_screen 'conf' file with the "
                   "'mammalian' genome indices")
    p.add_argument("-c","--contaminants",
                   dest="contaminants_conf",
                   help="fastq_screen 'conf' file with the "
                   "'contaminant' genome indices")
    p.add_argument("-a","--aligner",
                   dest="aligner",default=None,
                   choices=["bowtie","bowtie2"],
                   help="aligner to use with fastq_screen (default: "
                   "don't specify the aligner)")
    p.add_argument("--no-quality-filter",action='store_true',
                   dest="no_quality_filter",
                   help="turn off the barcode/UMI quality checks "
                   "(recommended for NextSeq data)")
    p.add_argument("--no-cleanup",action='store_true',
                   dest="no_cleanup",
                   help="don't remove intermediate Fastq files "
                   "(default is to delete intermediate Fastqs once "
                   "no longer needed)")
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
                   "arguments can be specified (default: '%s')" %
                   __settings.general.default_runner)
    p.add_argument("-s","--size",type=int,
                   dest="batch_size",default=DEFAULT_BATCH_SIZE,
                   help="number of reads per batch when splitting "
                   "FASTQ files for processing (default: %s)" %
                   DEFAULT_BATCH_SIZE)
    p.add_argument("-j","--max-jobs",type=int,
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
        default_runner = __settings.general.default_runner
    for stage in stages:
        if stage not in runners:
            runners[stage] = default_runner

    # Other settings
    well_list = os.path.abspath(args.well_list)
    max_jobs = args.max_jobs
    do_quality_filter = (not args.no_quality_filter)
    do_clean_up = (not args.no_cleanup)

    # Report settings
    print "Unaligned dir     : %s" % args.unaligned_dir
    print "Well list file    : %s" % well_list
    print "Output dir        : %s" % args.outdir
    print "Batch size (reads): %s" % args.batch_size
    print "Quality filter barcodes/UMIs: %s" % \
        ('yes' if do_quality_filter else 'no')
    print "Mammalian genome panel  : %s" % args.mammalian_conf
    with open(args.mammalian_conf) as fp:
        for line in fp:
            if line.startswith("DATABASE"):
                print "-- %s" % line.split('\t')[1]
    print "Contaminant genome panel: %s" % args.contaminants_conf
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
    print "Clean-up intermediate Fastqs: %s" % \
        ('yes' if do_clean_up else 'no')

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
    ppl = Pipeline(name="Process ICell8")

    # Initial stats
    initial_stats = GetICell8Stats("Initial statistics",
                                   fastqs,
                                   os.path.join(stats_dir,"icell8_stats.tsv"),
                                   well_list,
                                   nprocs=args.threads)
    ppl.add_task(initial_stats,
                 runner=runners['contaminant_filter'])

    # Split fastqs into batches
    batch_dir = os.path.join(icell8_dir,"_fastqs.batched")
    batch_fastqs = SplitFastqsIntoBatches("Batch Fastqs",fastqs,
                                          batch_dir,basename,
                                          batch_size=args.batch_size)
    ppl.add_task(batch_fastqs)

    # Setup the filtering jobs as a group
    filter_dir = os.path.join(icell8_dir,"_fastqs.quality_filter")
    filter_fastqs = FilterICell8Fastqs("Filter Fastqs",
                                       batch_fastqs.output(),
                                       filter_dir,
                                       well_list=well_list,
                                       mode='none',
                                       discard_unknown_barcodes=True,
                                       quality_filter=do_quality_filter)
    ppl.add_task(filter_fastqs,dependencies=(batch_fastqs,))
    
    # Post filtering stats
    filter_stats = GetICell8Stats("Post-filtering statistics",
                                  filter_fastqs.output().assigned,
                                  initial_stats.output(),
                                  suffix="_filtered",
                                  append=True,
                                  nprocs=args.threads)
    ppl.add_task(filter_stats,dependencies=(initial_stats,filter_fastqs),
                 runner=runners['contaminant_filter'])

    # Use cutadapt to find reads with poly-G regions
    poly_g_dir = os.path.join(icell8_dir,"_fastqs.poly_g")
    get_poly_g_reads = GetReadsWithPolyGRegions(
        "Find reads with poly-G regions",
        filter_fastqs.output().assigned,
        poly_g_dir)
    ppl.add_task(get_poly_g_reads,dependencies=(filter_fastqs,))
    poly_g_stats = GetICell8PolyGStats("Poly-G region statistics",
                                       get_poly_g_reads.output(),
                                       initial_stats.output(),
                                       suffix="_poly_g",
                                       append=True,
                                       nprocs=args.threads)
    ppl.add_task(poly_g_stats,dependencies=(get_poly_g_reads,filter_stats),
                 runner=runners['contaminant_filter'])

    # Set up the cutadapt jobs as a group
    trim_dir = os.path.join(icell8_dir,"_fastqs.trim_reads")
    trim_reads = TrimReads("Read trimming",
                           filter_fastqs.output().assigned,
                           trim_dir)
    ppl.add_task(trim_reads,dependencies=(filter_fastqs,))

    # Post read trimming stats
    trim_stats = GetICell8Stats("Post-trimming statistics",
                                trim_reads.output(),
                                initial_stats.output(),
                                suffix="_trimmed",
                                append=True,
                                nprocs=args.threads)
    ppl.add_task(trim_stats,dependencies=(trim_reads,poly_g_stats),
                 runner=runners['contaminant_filter'])

    # Set up the contaminant filter jobs as a group
    contaminant_filter_dir = os.path.join(icell8_dir,
                                          "_fastqs.contaminant_filter")
    contaminant_filter = FilterContaminatedReads("Contaminant filtering",
                                                 trim_reads.output(),
                                                 contaminant_filter_dir,
                                                 args.mammalian_conf,
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
                                 append=True,
                                 nprocs=args.threads)
    ppl.add_task(final_stats,dependencies=(contaminant_filter,trim_stats),
                 runner=runners['contaminant_filter'])

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
                               filter_fastqs.output().unassigned,
                               filter_fastqs.output().failed_barcodes,
                               filter_fastqs.output().failed_umis,
                               final_fastqs_dir,
                               basename)
    ppl.add_task(merge_fastqs,dependencies=(split_barcodes,))

    # Final stats for verification
    final_barcode_stats = GetICell8Stats(
        "Post-barcode splitting and merging statistics",
        merge_fastqs.output().assigned,
        initial_stats.output(),
        suffix="_final",
        append=True,
        nprocs=args.threads)
    ppl.add_task(final_barcode_stats,dependencies=(merge_fastqs,final_stats),
                 runner=runners['contaminant_filter'])

    # Run the QC
    run_qc = RunQC("Run QC",
                   args.outdir,
                   nthreads=args.threads)
    ppl.add_task(run_qc,dependencies=(merge_fastqs,),
                 runner=runners['contaminant_filter'])

    # Cleanup outputs
    cleanup_batch_fastqs = CleanupDirectory("Remove batched Fastqs",
                                            batch_dir)
    cleanup_quality_filter = CleanupDirectory("Remove filtered Fastqs",
                                              filter_dir)
    cleanup_poly_g = CleanupDirectory("Remove poly-G region stats data",
                                      poly_g_dir)
    cleanup_trim_reads = CleanupDirectory("Remove trimmed Fastqs",
                                          trim_dir)
    cleanup_contaminant_filtered = CleanupDirectory("Remove contaminant "
                                                    "filtered Fastqs",
                                                    contaminant_filter_dir)
    cleanup_split_barcodes = CleanupDirectory("remove barcode split Fastqs",
                                              barcoded_fastqs_dir)
    if do_clean_up:
        ppl.add_task(cleanup_batch_fastqs,dependencies=(filter_fastqs,))
        ppl.add_task(cleanup_quality_filter,dependencies=(trim_reads,
                                                          get_poly_g_reads,
                                                          merge_fastqs,
                                                          filter_stats))
        ppl.add_task(cleanup_poly_g,dependencies=(poly_g_stats,))
        ppl.add_task(cleanup_trim_reads,dependencies=(contaminant_filter,
                                                      trim_stats))
        ppl.add_task(cleanup_contaminant_filtered,dependencies=(final_stats,
                                                                split_barcodes))
        ppl.add_task(cleanup_split_barcodes,dependencies=(merge_fastqs,))

    # Execute the pipeline
    exit_status = ppl.run(sched=sched,log_dir=log_dir,scripts_dir=scripts_dir)

    # Finish
    sched.stop()
    if exit_status != 0:
        logging.critical("Pipeline failed: exit status %s" % exit_status)
    else:
        print "Pipeline completed ok"
    sys.exit(exit_status)
