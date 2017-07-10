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
import string
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
from bcftbx.simple_xls import XLSWorkBook
from auto_process_ngs.applications import Command
from auto_process_ngs.simple_scheduler import SimpleScheduler
from auto_process_ngs.simple_scheduler import SchedulerReporter
from auto_process_ngs.simple_scheduler import SchedulerGroup
from auto_process_ngs.fastq_utils import pair_fastqs
from auto_process_ngs.fastq_utils import get_read_number
from auto_process_ngs.utils import BaseFastqAttrs
from auto_process_ngs.utils import AnalysisFastq
from auto_process_ngs.utils import AnalysisProject
from auto_process_ngs.icell8_utils import ICell8WellList
from auto_process_ngs.qc.illumina_qc import check_qc_outputs
import auto_process_ngs.envmod as envmod

# Fetch configuration settings
import auto_process_ngs.settings
__settings = auto_process_ngs.settings.Settings()

# Module specific logger
logger = logging.getLogger(__name__)

######################################################################
# Pipeline infrastructure constants
######################################################################

ALLOWED_CHARS = string.lowercase + string.digits + "._-"

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
    >> t2 = p.add_task(Task2(),requires=(t1,))
    >> ...
    >> p.run()

    Tasks will only run when all requirements have
    completed (or will run immediately if they don't have
    any requirements).
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
    def add_task(self,task,requires=(),**kws):
        self._pending.append((task,requires,kws))
        self.report("Adding task '%s'" % task.name())
        if requires:
            for req in requires:
                if req.name() not in [t[0].name() for t in self._pending]:
                    self.report("-> Adding requirement '%s'" % req.name())
                    self.add_task(req)
        return task
    def report(self,s):
        print "%s [%s] %s" % (time.strftime("%Y-%m-%d %H:%M:%S"),
                              self._name,s)
    def terminate(self):
        # Terminates the pipeline
        self.report("Terminating pipeline")
        # Dump all pending tasks
        self._pending = []
        # Stop all running tasks
        if self._running:
            self.report("Stopping running tasks:")
            for task in self._running:
                self.report("- %s" % task.name())
                task.terminate()
        return 1
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
            # Check for pending tasks that can start
            pending = []
            for task,requirements,kws in self._pending:
                run_task = False
                if not requirements:
                    # No requirements - start it
                    run_task = True
                else:
                    # Check requirements
                    run_task = reduce(lambda x,y: x and y.completed,
                                      requirements,True)
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
                        logger.critical("Failed to start task '%s': %s" %
                                        (task.name(),ex))
                        failed.append(task)
                    self._running.append(task)
                    update = True
                else:
                    pending.append((task,requirements,kws))
            self._pending = pending
            # Check for tasks that have failed
            if failed:
                self.report("Following tasks failed:")
                for task in failed:
                    self.report("- %s" % task.name())
                return self.terminate()
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
    def __init__(self,_name,*args,**kws):
        self._name = str(_name)
        self._args = args
        self._kws = kws
        self._commands = []
        self._task_name = "%s.%s" % (sanitize_name(self._name),
                                     uuid.uuid4())
        self._completed = False
        self._stdout_files = []
        self._exit_code = 0
        # Running jobs
        self._jobs = []
        self._groups = []
        # Deal with subclass arguments
        try:
            self._callargs = inspect.getcallargs(self.init,*args,**kws)
        except Exception as ex:
            logger.error("Exception setting up args for task '%s' (%s): %s"
                         % (self._name,self.__class__,ex))
            raise ex
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
    @property
    def stdout(self):
        stdout = []
        for f in self._stdout_files:
            with open(f,'r') as fp:
                stdout.append(fp.read())
        return '\n'.join(stdout)
    def name(self):
        return self._task_name
    def fail(self,exit_code=1,message=None):
        # Register the task as failing
        if message:
            self.report("failed: %s" % message)
        self.report("failed: exit code set to %s" % exit_code)
        self._completed = True
        self._exit_code = exit_code
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
                self._stdout_files.append(job.log)
            except AttributeError:
                # Assume it's a group
                for j in job.jobs:
                    if j.exit_code != 0:
                        self._exit_code += 1
                    self._stdout_files.append(j.log)
        if self._exit_code != 0:
            logger.critical("%s failed: exit code %s"
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
                self._groups.append(group)
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
                self._jobs.append(job)
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
                    logger.warning("%s: failed (exit code %s)" %
                                   (self._name,exit_code))
        else:
            # Run each stage locally
            for cmd in cmds:
                cmd.run_subprocess(working_dir=working_dir)
            self.invoke(self.finish)
            self._completed = True
        return self
    def terminate(self):
        # Terminate the job
        if self.completed:
            return
        for group in self._groups:
            for job in group.jobs:
                job.terminate()
        for job in self._jobs:
            job.terminate()
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

    To make an instance using this wrapper class:

    >>> ls_command = LsCommand(path)
    """
    def __init__(self,*args,**kws):
        # Set internal name
        self._name = self.__class__.__name__
        # Invoke the 'init' method
        self.init(*args,**kws)
    def name(self):
        return sanitize_name(self._name)
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

class PipelineCommandWrapper(PipelineCommand):
    """
    Class for constructing program command lines

    This class is based on the PipelineCommand class but
    can be used directly (rather than needing to be
    subclassed).

    For example, to wrap the 'ls' command directly:

    >>> ls_command = PipelineCommandWrapper('ls',dirn)

    It is also possible to extend the command line
    using the 'add_args' method, for example:

    >>> ls_command = PipelineCommandWrapper('ls')
    >>> ls.command.add_args(dirn)
    """
    def __init__(self,*args):
        PipelineCommand.__init__(self,*args)
        self._name = self.__class__.__name__
        self._cmd = None
        if args:
            self._cmd = Command(*args)
    def add_args(self,*args):
        # Add additional arguments to extend
        # the command being built
        if self._cmd is None:
            self._cmd = Command(*args)
        else:
            self._cmd.add_args(*args)
    def init(self,*args):
        # Dummy init which does nothing
        pass
    def cmd(self):
        # Implement the 'cmd' method
        return self._cmd

######################################################################
# Generic pipeline functions
######################################################################

def sanitize_name(s):
    """
    Convert string to lowercase and replace special characters
    """
    name = []
    # Convert to lower case and replace special characters
    for c in str(s).lower():
        if c not in ALLOWED_CHARS:
            if len(name) < 2 or name[-2:] != '__':
                name.append('_')
        else:
            name.append(c)
    return ''.join(name)

######################################################################
# Magic numbers
######################################################################

DEFAULT_BATCH_SIZE = 5000000

######################################################################
# ICell8 pipeline commands
######################################################################

class ICell8Statistics(PipelineCommand):
    """
    Build command to run the 'icell8_stats.py' utility
    """
    def init(self,fastqs,stats_file,well_list=None,
             suffix=None,unassigned=False,append=False,
             nprocs=1):
        """
        Create new ICell8Statistics instance

        Arguments:
          fastqs (list): list of FASTQ file names
          stats_file (str): path to output file
          well_list (str): path to 'well list' file (optional)
          suffix (str): suffix to append to columns with read and
            UMI counts (optional)
          unassigned (bool): if True then also collect stats for
            read pairs that don't match any of the expected
            barcodes from the well list or existing stats file
            (by default unassigned stats are not collected)
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
        self._unassigned = unassigned
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
            cmd.add_args('--append')
        if self._unassigned:
            cmd.add_args('--unassigned')
        cmd.add_args(*self._fastqs)
        return cmd

class SplitAndFilterFastqPair(PipelineCommand):
    """
    Build command to run the 'split_icell8_fastqs.py' utility
    """
    def init(self,fastq_pair,out_dir,well_list=None,
             basename=None,mode='none',
             discard_unknown_barcodes=False,
             quality_filter=False,
             compress=False):
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
          compress (bool): if True then gzip the
            output files (FASTQs are uncompressed
            by default)
        """
        self._fastq_pair = fastq_pair
        self._out_dir = os.path.abspath(out_dir)
        self._well_list = well_list
        self._basename = basename
        self._mode = mode
        self._discard_unknown_barcodes = discard_unknown_barcodes
        self._quality_filter = quality_filter
        self._compress = compress
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
        if self._compress:
            cmd.add_args('--compress')
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

    If the output FASTQ names end with .gz then they will be
    automatically compressed with gzip after concatenation.

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
        compress = self._fastq_out.endswith('.gz')
        if compress:
            fastq_out = '.'.join(self._fastq_out.split('.')[:-1])
        else:
            fastq_out = self._fastq_out
        fastq_out = os.path.join(self._concat_dir,fastq_out)
        cmd = Command('cat')
        cmd.add_args(*self._fastqs)
        cmd.add_args('>',fastq_out)
        if compress:
            cmd.add_args('&&','gzip',fastq_out)
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
    def init(self,fastqs,nthreads=1,working_dir=None,qc_dir=None):
        """
        Set up parameters
        """
        self.fastqs = fastqs
        self.nthreads = nthreads
        if working_dir is None:
            working_dir = os.getcwd()
        self.working_dir = os.path.abspath(working_dir)
        self.qc_dir = qc_dir
    def cmd(self):
        """
        Build the command
        """
        cmd = Command('cd',self.working_dir)
        for fastq in self.fastqs:
            cmd.add_args('&&','illumina_qc.sh',fastq)
            if self.nthreads > 1:
                cmd.add_args('--threads',self.nthreads)
            if self.qc_dir is not None:
                cmd.add_args('--qc_dir',self.qc_dir)
        return cmd

class MultiQC(PipelineCommand):
    """
    Run the MultiQC program on a set of QC outputs
    """
    def init(self,qc_dir,out_file,title):
        """
        Set up parameters
        """
        self.qc_dir = qc_dir
        self.out_file = out_file
        self.title = title
    def cmd(self):
        """
        Build the command
        """
        return Command('multiqc',
                       '--title',self.title,
                       '--filename',self.out_file,
                       '--force',
                       self.qc_dir)

######################################################################
# ICell8 pipeline tasks
######################################################################

class GetICell8Stats(PipelineTask):
    """
    """
    def init(self,fastqs,stats_file,well_list=None,
             suffix=None,unassigned=False,append=False,
             nprocs=1):
        pass
    def setup(self):
        self.add_cmd(ICell8Statistics(self.args.fastqs,
                                      self.args.stats_file,
                                      well_list=self.args.well_list,
                                      suffix=self.args.suffix,
                                      append=self.args.append,
                                      unassigned=self.args.unassigned,
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
            except ZeroDivisionError:
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
    def finish(self):
        print self.stdout
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

class MergeBarcodeFastqs(PipelineTask):
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
                      "%s.unprocessed" % self.args.merge_dir)
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
                                                 mode="barcodes",
                                                 compress=True))
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
                                      "%s.%s.r1.fastq.gz" %
                                      (self.args.basename,name)))
            fqs_r2 = [p[1] for p in fastq_pairs]
            self.add_cmd(ConcatFastqs(fqs_r2,
                                      self.args.merge_dir,
                                      "%s.%s.r2.fastq.gz" %
                                      (self.args.basename,name)))
    def output(self):
        out_dir = self.args.merge_dir
        return AttributeDictionary(
            assigned=FileCollection(out_dir,"*.[ACGT]*.r*.fastq.gz"),
            unassigned=FileCollection(out_dir,"*.unassigned.r*.fastq.gz"),
            failed_barcodes=FileCollection(out_dir,"*.failed_barcodes.r*.fastq.gz"),
            failed_umis=FileCollection(out_dir,"*.failed_umis.r*.fastq.gz"),
        )

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

class MergeBarcodeFastqs(PipelineTask):
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
                      "%s.unprocessed" % self.args.merge_dir)
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
                                                 mode="barcodes",
                                                 compress=True))
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
                                      "%s.%s.r1.fastq.gz" %
                                      (self.args.basename,name)))
            fqs_r2 = [p[1] for p in fastq_pairs]
            self.add_cmd(ConcatFastqs(fqs_r2,
                                      self.args.merge_dir,
                                      "%s.%s.r2.fastq.gz" %
                                      (self.args.basename,name)))
    def output(self):
        out_dir = self.args.merge_dir
        return AttributeDictionary(
            assigned=FileCollection(out_dir,"*.[ACGT]*.r*.fastq.gz"),
            unassigned=FileCollection(out_dir,"*.unassigned.r*.fastq.gz"),
            failed_barcodes=FileCollection(out_dir,"*.failed_barcodes.r*.fastq.gz"),
            failed_umis=FileCollection(out_dir,"*.failed_umis.r*.fastq.gz"),
        )

class MergeSampleFastqs(PipelineTask):
    """
    """
    def init(self,fastqs,well_list,merge_dir):
        pass
    def setup(self):
        # Make the output directory
        mkdir(self.args.merge_dir)
        # Group fastqs by sample
        well_list = ICell8WellList(self.args.well_list)
        fastq_groups = dict()
        for fq in self.args.fastqs:
            barcode = os.path.basename(fq).split('.')[-3]
            sample = well_list.sample(barcode)
            try:
                fastq_groups[sample].append(fq)
            except KeyError:
                fastq_groups[sample] = [fq,]
        # Set up merge for fastq pairs in each sample
        for sample in fastq_groups:
            fastq_pairs = pair_fastqs(fastq_groups[sample])[0]
            fqs_r1 = [p[0] for p in fastq_pairs]
            self.add_cmd(ConcatFastqs(fqs_r1,
                                      self.args.merge_dir,
                                      "%s.r1.fastq.gz" % sample))
            fqs_r2 = [p[1] for p in fastq_pairs]
            self.add_cmd(ConcatFastqs(fqs_r2,
                                      self.args.merge_dir,
                                      "%s.r2.fastq.gz" % sample))
    def output(self):
        out_dir = self.args.merge_dir
        return AttributeDictionary(
            fastqs=FileCollection(out_dir,"*.r*.fastq.gz"),
        )

class RunQC(PipelineTask):
    """
    """
    def init(self,project_dir,nthreads=1,batch_size=25,
             fastq_dir='fastqs',qc_dir='qc'):
        self.qc_dir = None
        self.qc_report = None
        self.fastq_attrs = ICell8FastqAttrs
    def setup(self):
        # Gather Fastqs
        project = AnalysisProject(
            os.path.basename(self.args.project_dir),
            self.args.project_dir,
            fastq_dir=self.args.fastq_dir,
            fastq_attrs=self.fastq_attrs)
        batch_size = self.args.batch_size
        # Make the output qc directory
        self.qc_dir = project.setup_qc_dir(self.args.qc_dir)
        print "QC dir: %s" % self.qc_dir
        # Gather fastqs to run QC on
        fastqs = []
        for sample in project.samples:
            for fq in sample.fastq:
                if not sample.verify_qc(self.qc_dir,fq):
                    fastqs.append(fq)
        # Set up QC run for batches of fastqs
        while fastqs:
            self.add_cmd(IlluminaQC(fastqs[:batch_size],
                                    nthreads=self.args.nthreads,
                                    working_dir=self.args.project_dir,
                                    qc_dir=self.qc_dir))
            fastqs = fastqs[batch_size:]
    def finish(self):
        # Verify the QC outputs for each fastq
        project = AnalysisProject(
            os.path.basename(self.args.project_dir),
            self.args.project_dir,
            fastq_dir=self.args.fastq_dir,
            fastq_attrs=self.fastq_attrs)
        if not project.verify_qc(qc_dir=self.qc_dir):
            print "Failed to verify QC"
            self.fail(exit_code=1)
        else:
            # Generate QC report
            print "Generating QC report"
            qc_base = os.path.basename(self.qc_dir)
            out_file = os.path.join(project.dirn,
                                    '%s_report.html' % qc_base)
            if project.info.run is not None:
                title = "%s/%s" % (project.info.run,
                                   project.name)
            else:
                title = "%s" % project.name
            if self.args.fastq_dir is not None:
                title = "%s (%s)" % (title,self.args.fastq_dir)
            title = "%s: QC report" % title
            print "Title : '%s'" % title
            print "Report: %s" % out_file
            self.qc_report = project.qc_report(qc_dir=self.qc_dir,
                                               title=title,
                                               report_html=out_file,
                                               force=True)
            if self.qc_report is None:
                print "Failed to generate QC report"
                self.fail(exit_code=1)
            else:
                print "QC report: %s" % self.qc_report
    def output(self):
        return AttributeDictionary(
            qc_dir=self.qc_dir,
            report_zip=self.qc_report,
        )

class RunMultiQC(PipelineTask):
    """
    """
    def init(self,project_dir,fastq_dir='fastqs',qc_dir='qc'):
        self.multiqc_out = None
        self.fastq_attrs = ICell8FastqAttrs
    def setup(self):
        project = AnalysisProject(
            os.path.basename(self.args.project_dir),
            self.args.project_dir,
            fastq_dir=self.args.fastq_dir,
            fastq_attrs=self.fastq_attrs)
        project.use_qc_dir(self.args.qc_dir)
        multiqc_out = os.path.join(project.dirn,
                                   "multi%s_report.html" % \
                                   os.path.basename(project.qc_dir))
        if project.info.run is not None:
            title = "%s/%s" % (project.info.run,
                               project.name)
        else:
            title = "%s" % project.name
        self.add_cmd(MultiQC(project.qc_dir,
                             multiqc_out,
                             title))
        self.multiqc_out = multiqc_out
    def output(self):
        return AttributeDictionary(
            multiqc_out=self.multiqc_out,
        )

class CheckICell8Barcodes(PipelineTask):
    """
    Check the barcodes are consistent

    This is a sanity check: ensure that the inline barcodes
    for all reads in the R1 Fastq for the barcode Fastq pairs
    matches the assigned barcode.
    """
    def init(self,fastqs):
        self.bad_barcodes = None
    def setup(self):
        fastq_attrs =ICell8FastqAttrs
        batch_size = 25
        # Reduce fastq list to just R1 files
        fastqs = filter(lambda fq:
                        fastq_attrs(fq).read_number == 1,
                        self.args.fastqs)
        # Set up verification on batches of fastqs
        while fastqs:
            cmd = PipelineCommandWrapper()
            for fq in fastqs[:batch_size]:
                if fastq_attrs(fq).extension.endswith('.gz'):
                    cat = 'zcat'
                else:
                    cat = 'cat'
                barcode = fastq_attrs(fq).barcode_sequence
                if cmd.cmd() is not None:
                    cmd.add_args("&&")
                cmd.add_args(
                    "echo","-n","%s:" % barcode,"&&",
                    "%s" % cat,fq,"|",
                    "sed","-n","'2~4p'","|",
                    "grep","-v","^%s" % barcode,"|",
                    "wc","-l"
                )
            self.add_cmd(cmd)
            fastqs = fastqs[batch_size:]
    def finish(self):
        # Read the stdout looking for barcodes
        # with non-zero counts
        print self.stdout
        self.bad_barcodes = []
        for line in self.stdout.split('\n'):
            try:
                barcode,count = line.split(':')
                count = int(count)
                if count > 0:
                    self.bad_barcodes.append((barcode,count))
            except ValueError:
                pass
        # If there are bad barcodes then fail
        if self.bad_barcodes:
            for barcode in self.bad_barcodes:
                print "ERROR %s: %d non-matching reads" % (barcode[0],
                                                           barcode[1])
            self.fail()
    def output(self):
        return AttributeDictionary(
            bad_barcodes=self.bad_barcodes
        )

class ConvertStatsToXLSX(PipelineTask):
    """
    Convert the stats file to XLSX format
    """
    def init(self,stats_file,xlsx_file):
        pass
    def setup(self):
        convert_to_xlsx(self.args.stats_file,
                        self.args.xlsx_file,
                        title="ICell8 stats",
                        freeze_header=True)
    def output(self):
        return AttributeDictionary(
            xlsx_file=self.args.xlsx_file
        )

class ReportProcessing(PipelineTask):
    """
    Generate an HTML report on the processing
    """
    def init(self,stats_file,out_file=None,name=None):
        self.out_file = None
    def setup(self):
        if self.args.out_file is None:
            self.out_file = "icell8_processing.html"
        else:
            self.out_file = self.args.out_file
        self.out_file = os.path.abspath(self.out_file)
        cmd = PipelineCommandWrapper(
            'icell8_report.py')
        if self.args.name is not None:
            cmd.add_args('--name',self.args.name)
        cmd.add_args(self.args.stats_file,
                     self.out_file)
        self.add_cmd(cmd)
    def output(self):
        return AttributeDictionary(
            report_html=self.out_file
        )

class CleanupDirectory(PipelineTask):
    """
    """
    def init(self,dirn):
        pass
    def setup(self):
        dirn = os.path.abspath(self.args.dirn)
        if not os.path.isdir(dirn):
            self.report("No directory '%s'" % self.args.dirn)
        else:
            self.add_cmd(
                PipelineCommandWrapper(
                    "rm","-f","%s" % os.path.join(dirn,'*'),
                    "&&",
                    "rmdir","%s" % dirn))

######################################################################
# Classes
######################################################################

class ICell8FastqAttrs(BaseFastqAttrs):
    """
    Extract attributes from Icell8 pipeline Fastq names

    Used as a custom Fastq name attribute parser for output
    Fastqs from the ICell8 pipeline.

    Fastq names for the final barcoded files look like e.g.:

    icell8.CCAGTTCAGGA.r1.fastq

    For the samples e.g.:

    d1.2.r1.fastq.gz
    """
    def __init__(self,fastq):
        BaseFastqAttrs.__init__(self,fastq)
        name = self.basename.split('.')
        # Handle the file extension
        extension = []
        if name[-1] == "gz":
            extension.append("gz")
            name = name[:-1]
        if name[-1] in ("fastq","fq"):
            extension.insert(0,name[-1])
            name = name[:-1]
        if extension:
            self.extension = '.'.join(extension)
        # Get the read number
        if name[-1].startswith('r'):
            try:
                self.read_number = int(name[-1][1:])
                name = name[:-1]
            except ValueError:
                # Can't get a read number
                self.sample_name = '.'.join(name)
                return
        # Get the barcode sequence
        if (len(name) > 1) and \
           (len(name[-1]) > 0) and \
           (not filter(lambda c: c not in "AGCT",name[-1])):
            self.barcode_sequence = name[-1]
        # Assume sample name is whatever's left over
        self.sample_name = '.'.join(name)
    def __repr__(self):
        name = [self.sample_name]
        if self.read_number:
            name.append("r%s" % self.read_number)
        return '.'.join(name)

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

def convert_to_xlsx(tsv_file,xlsx_file,title=None,freeze_header=False):
    """
    Convert a tab-delimited file to an XLSX file

    Arguments:
      tsv_file (str): path to the input TSV file
      xlsx_file (str): path to the output XLSX file
      title (str): optional, name to give the worksheet in
        the output XLSX file (defaults to the input file name)
      freeze_header (bool): optional, if True then 'freezes'
        the first line of the XLSX file (default is not to
        freeze the first line)
    """
    if title is None:
        title = os.path.basename(tsv_file)
    wb = XLSWorkBook(title)
    ws = wb.add_work_sheet(title)
    with open(tsv_file,'r') as stats:
        for line in stats:
            ws.append_row(data=line.rstrip('\n').split('\t'))
    # Freeze the top row
    if freeze_header:
        ws.freeze_panes = 'A2'
    wb.save_as_xlsx(xlsx_file)

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
                   dest="unaligned_dir",default=None,
                   help="process FASTQs from 'unaligned' dir with output "
                   "from bcl2fastq (NB cannot be used with -p option)")
    p.add_argument("-p","--project",metavar="NAME",
                   dest="project",default=None,
                   help="process FASTQS from project directory NAME (NB "
                   "if -o not specified then this will also be used as "
                   "the output directory; cannot be used with -u option)")
    p.add_argument("-o","--outdir",
                   dest="outdir",default=None,
                   help="directory to write outputs to "
                   "(default: 'CWD/icell8', or project dir if -p "
                   "is specified)")
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
            logger.fatal("Bad stage for --runner option: %s" % stage)
            sys.exit(1)
        runners[stage] = fetch_runner(runner_spec)
    try:
        default_runner = runners['default']
    except KeyError:
        default_runner = __settings.general.default_runner
    for stage in stages:
        if stage not in runners:
            runners[stage] = default_runner

    # Check for clashing -u/-p
    if args.project and args.unaligned_dir:
        logger.fatal("Cannot specify -u and -p together")
        sys.exit(1)

    # Output dir
    if args.outdir is None:
        if args.project:
            outdir = args.project
        else:
            outdir = "icell8"
    else:
        outdir = args.outdir

    # Other settings
    well_list = os.path.abspath(args.well_list)
    max_jobs = args.max_jobs
    do_quality_filter = (not args.no_quality_filter)
    do_clean_up = (not args.no_cleanup)

    # Report settings
    print "Unaligned dir     : %s" % args.unaligned_dir
    print "Project           : %s" % args.project
    print "Well list file    : %s" % well_list
    print "Output dir        : %s" % outdir
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

    # Check well list file
    try:
        ICell8WellList(well_list).barcodes()
    except Exception as ex:
        logger.fatal("Couldn't load data from well list file '%s'"
                     % well_list)
        sys.exit(1)

    # Get the input FASTQ file pairs
    fastqs = []
    # Collect files from command line
    for fq in args.fastqs:
        fastqs.append(os.path.abspath(fq))
    # Collect files from unaligned dir
    if fastqs and args.unaligned_dir is not None:
        logger.warning("Ignoring unaligned dir '%s'" %
                       args.unaligned_dir)
    elif args.unaligned_dir:
        try:
            illumina_data = IlluminaData(
                os.getcwd(),
                unaligned_dir=args.unaligned_dir)
            for project in illumina_data.projects:
                for sample in project.samples:
                    for fq in sample.fastq:
                        fastqs.append(os.path.join(sample.dirn,fq))
        except IlluminaDataError:
            logger.fatal("Couldn't find FASTQS in directory '%s'" %
                          args.unaligned_dir)
    # Collect files from project
    analysis_project = None
    if fastqs and args.project is not None:
        logger.warning("Ignoring project '%s'" % args.project)
    elif args.project:
        analysis_project = AnalysisProject(args.project,
                                           args.project)
        for sample in analysis_project.samples:
            for fq in sample.fastq:
                fastqs.append(os.path.join(
                    analysis_project.fastq_dir,
                    fq))
    if not fastqs:
        logger.fatal("No FASTQs found")
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
    icell8_dir = os.path.abspath(outdir)
    if os.path.exists(icell8_dir) and args.project is None:
        if not args.force:
            logger.fatal("Output destination '%s': already exists "
                         "(remove or use --force to overwrite)" %
                         icell8_dir)
            sys.exit(1)
        logger.warning("Removing existing output destination '%s'" %
                       icell8_dir)
        shutil.rmtree(icell8_dir)
    log_dir = os.path.join(icell8_dir,"logs")
    stats_dir = os.path.join(icell8_dir,"stats")
    scripts_dir = os.path.join(icell8_dir,"scripts")
    for dirn in (icell8_dir,log_dir,stats_dir,scripts_dir):
        mkdir(dirn)

    # Copy well list file into output directory
    shutil.copy(well_list,outdir)
    well_list = os.path.join(outdir,os.path.basename(well_list))

    # Set up a pipeline
    ppl = Pipeline(name="Process ICell8")

    # Initial stats
    initial_stats = GetICell8Stats("Initial statistics",
                                   fastqs,
                                   os.path.join(stats_dir,"icell8_stats.tsv"),
                                   well_list,
                                   unassigned=True,
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
    ppl.add_task(filter_fastqs,requires=(batch_fastqs,))
    
    # Post filtering stats
    filter_stats = GetICell8Stats("Post-filtering statistics",
                                  filter_fastqs.output().assigned,
                                  initial_stats.output(),
                                  suffix="_filtered",
                                  append=True,
                                  nprocs=args.threads)
    ppl.add_task(filter_stats,requires=(initial_stats,filter_fastqs),
                 runner=runners['contaminant_filter'])

    # Use cutadapt to find reads with poly-G regions
    poly_g_dir = os.path.join(icell8_dir,"_fastqs.poly_g")
    get_poly_g_reads = GetReadsWithPolyGRegions(
        "Find reads with poly-G regions",
        filter_fastqs.output().assigned,
        poly_g_dir)
    ppl.add_task(get_poly_g_reads,requires=(filter_fastqs,))
    poly_g_stats = GetICell8PolyGStats("Poly-G region statistics",
                                       get_poly_g_reads.output(),
                                       initial_stats.output(),
                                       suffix="_poly_g",
                                       append=True,
                                       nprocs=args.threads)
    ppl.add_task(poly_g_stats,requires=(get_poly_g_reads,filter_stats),
                 runner=runners['contaminant_filter'])

    # Set up the cutadapt jobs as a group
    trim_dir = os.path.join(icell8_dir,"_fastqs.trim_reads")
    trim_reads = TrimReads("Read trimming",
                           filter_fastqs.output().assigned,
                           trim_dir)
    ppl.add_task(trim_reads,requires=(filter_fastqs,))

    # Post read trimming stats
    trim_stats = GetICell8Stats("Post-trimming statistics",
                                trim_reads.output(),
                                initial_stats.output(),
                                suffix="_trimmed",
                                append=True,
                                nprocs=args.threads)
    ppl.add_task(trim_stats,requires=(trim_reads,poly_g_stats),
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
    ppl.add_task(contaminant_filter,requires=(trim_reads,),
                 runner=runners['contaminant_filter'])

    # Post contaminant filter stats
    final_stats = GetICell8Stats("Post-contaminant filter statistics",
                                 contaminant_filter.output(),
                                 initial_stats.output(),
                                 suffix="_contaminant_filtered",
                                 append=True,
                                 nprocs=args.threads)
    ppl.add_task(final_stats,requires=(contaminant_filter,trim_stats),
                 runner=runners['contaminant_filter'])

    # Rebatch reads by barcode and sample
    # First: split each batch by barcode
    barcoded_fastqs_dir = os.path.join(icell8_dir,"_fastqs.split_barcodes")
    split_barcodes = SplitByBarcodes("Split batches by barcode",
                                     contaminant_filter.output(),
                                     barcoded_fastqs_dir)
    ppl.add_task(split_barcodes,requires=(contaminant_filter,))
    # Merge (concat) fastqs into single pairs per barcode
    final_fastqs_dir = os.path.join(icell8_dir,"fastqs.barcodes")
    merge_fastqs = MergeBarcodeFastqs("Assemble reads by barcode",
                                      split_barcodes.output(),
                                      filter_fastqs.output().unassigned,
                                      filter_fastqs.output().failed_barcodes,
                                      filter_fastqs.output().failed_umis,
                                      final_fastqs_dir,
                                      basename)
    ppl.add_task(merge_fastqs,requires=(split_barcodes,))
    # Merge (concat) fastqs into single pairs per barcode
    sample_fastqs_dir = os.path.join(icell8_dir,"fastqs.samples")
    sample_fastqs = MergeSampleFastqs("Assemble reads by sample",
                                      split_barcodes.output(),
                                      well_list,
                                      sample_fastqs_dir)
    ppl.add_task(sample_fastqs,requires=(split_barcodes,))

    # Final stats for verification
    final_barcode_stats = GetICell8Stats(
        "Post-barcode splitting and merging statistics",
        merge_fastqs.output().assigned,
        initial_stats.output(),
        suffix="_final",
        append=True,
        nprocs=args.threads)
    ppl.add_task(final_barcode_stats,requires=(merge_fastqs,final_stats),
                 runner=runners['contaminant_filter'])

    # Verify that barcodes are okay
    check_barcodes = CheckICell8Barcodes(
        "Verify barcodes are consistent",
        merge_fastqs.output().assigned)
    ppl.add_task(check_barcodes,requires=(merge_fastqs,))

    # Generate XLSX version of stats
    xlsx_stats = ConvertStatsToXLSX(
        "Convert statistics to XLSX",
        final_barcode_stats.output(),
        os.path.join(stats_dir,"icell8_stats.xlsx"))
    ppl.add_task(xlsx_stats,requires=(final_barcode_stats,))

    # Run the QC
    run_qc_barcodes = RunQC("Run QC for barcodes",
                            outdir,
                            nthreads=args.threads,
                            fastq_dir="fastqs.barcodes",
                            qc_dir="qc.barcodes")
    multiqc_barcodes = RunMultiQC("Run MultiQC for barcodes",
                                  outdir,
                                  fastq_dir="fastqs.barcodes",
                                  qc_dir="qc.barcodes")
    ppl.add_task(run_qc_barcodes,requires=(merge_fastqs,),
                 runner=runners['contaminant_filter'])
    ppl.add_task(multiqc_barcodes,requires=(run_qc_barcodes,))
    run_qc_samples = RunQC("Run QC for samples",
                           outdir,
                           nthreads=args.threads,
                           fastq_dir="fastqs.samples",
                           qc_dir="qc.samples")
    multiqc_samples = RunMultiQC("Run MultiQC for samples",
                                 outdir,
                                 fastq_dir="fastqs.samples",
                                 qc_dir="qc.samples")
    ppl.add_task(run_qc_samples,requires=(sample_fastqs,),
                 runner=runners['contaminant_filter'])
    ppl.add_task(multiqc_samples,requires=(run_qc_samples,))

    # Final report
    if analysis_project is not None:
        report_suffix = ".%s" % analysis_project.name
        if analysis_project.info.run is not None:
            report_suffix += ".%s" % analysis_project.info.run
    else:
        report_suffix = None
    final_report = ReportProcessing("Generate processing report",
                                    final_barcode_stats.output(),
                                    os.path.join(outdir,
                                                 "icell8_processing.html"),
                                    name=report_suffix)
    ppl.add_task(final_report,requires=(final_barcode_stats,))

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
        ppl.add_task(cleanup_batch_fastqs,requires=(filter_fastqs,))
        ppl.add_task(cleanup_quality_filter,requires=(trim_reads,
                                                      get_poly_g_reads,
                                                      merge_fastqs,
                                                      filter_stats))
        ppl.add_task(cleanup_poly_g,requires=(poly_g_stats,))
        ppl.add_task(cleanup_trim_reads,requires=(contaminant_filter,
                                                  trim_stats))
        ppl.add_task(cleanup_contaminant_filtered,requires=(final_stats,
                                                            split_barcodes))
        ppl.add_task(cleanup_split_barcodes,requires=(merge_fastqs,
                                                      sample_fastqs))

    # Execute the pipeline
    exit_status = ppl.run(sched=sched,log_dir=log_dir,scripts_dir=scripts_dir)

    # Finish
    sched.stop()
    if exit_status != 0:
        logger.critical("Pipeline failed: exit status %s" % exit_status)
    else:
        print "Pipeline completed ok"
    sys.exit(exit_status)
