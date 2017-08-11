#!/usr/bin/env python
#
#     pipeliner.py: utilities for building simple pipelines of tasks
#     Copyright (C) University of Manchester 2017 Peter Briggs
#
"""
pipeliner.py
============

Module providing utility classes and functions for building simple
'pipelines' of tasks.

The core classes are:

- Pipeline: class for building and executing pipelines
- PipelineTask: class for defining pipeline tasks
- PipelineCommand: class for defining commands that can be used in tasks
  (nb the ``PipelineCommandWrapper`` class is recommended over subclassing
  ``PipelineCommand``)

Additional supporting classes:

- PipelineCommandWrapper: shortcut alternative to PipelineCommand
- FileCollector: returning collections of files based on glob patterns

There are some underlying classes and functions that are intended for
internal use:

- Capturing: capture stdout from a Python function
- sanitize_name: clean up task and command names for use in pipeline
- collect_files: collect files based on glob patterns

Overview
--------

For the purposes of the ``pipeliner`` module, a 'pipeline' consists
of a set of 'tasks': a task can be independent of any other task, or
it may depend on one or more tasks being completed before it can
start.

Defining tasks: examples
------------------------

Tasks must be defined by subclassing the ``PipelineTask`` class and
implementing the following methods:

- init: used to declare any input parameters for the task
- setup: perform any set up (e.g. creating output directories) or
  other arbitary actions; typically it also defines any commands that
  will be sent to the scheduler, however this is optional.
- finish: perform final actions after setup and any defined commands
  have completed (nb ``finish`` is optional)
- output: return any outputs from the task once it has completed.

For example: let's define a task which runs the ``fastqc`` program
on a collection of Fastq files one at a time::

    class RunFastqc(PipelineTask):
        def init(self,fastqs,out_dir):
            self.out_files = list()
        def setup(self):
            if not os.path.exists(self.args.out_dir):
                os.mkdir(self.args.out_dir)
            for fq in self.args.fastqs:
                self.add_cmd(
                    PipelineCommandWrapper("Run FastQC",
                                           "fastqc",
                                           "-o",self.args.out_dir,
                                           fq))
        def finish(self):
            for fq in self.args.fastqs:
                if fq.endswith(".gz"):
                    fq = os.path.splitext(fq)[0]
                out_file = os.path.join(
                    self.args.out_dir,
                    os.path.splitext(
                        os.path.basename(fq))[0]+"_fastqc.html")
                if not os.path.exists(out_file):
                    self.fail(message="Missing output file: %s" % out_file)
                else:
                    self.out_files.append(out_file)
        def output(self):
            return self.out_files

The key features are:

1. The argument list of the ``init`` method defines an arbitrary
   set of parameters which are made available to the other methods
   via the ``self.args`` object.

   The ``init`` method also creates an attribute that is used to
   store the outputs of the task.

2. The ``setup`` method creates the output directory if it doesn't
   already exist, and then calls the ``add_cmd`` method to add a
   command for each Fastq file supplied via the ``fastqs`` argument
   (accessed as ``self.args.fastqs``), which will run ``fastqc`` on
   that Fastq.

3. The command to run ``fastqc`` is created by creating a
   ``PipelineCommandWrapper`` instance, which names the command
   and specifies the command and arguments to execute.

4. The commands defined via the ``add_cmd`` method are not guaranteed
   to run sequentially. If there are additional commands that rely on
   the first set of commands finishing, then either append these to the
   commands using the shell ``&&`` notation, or put them into a
   separate task which should be executed afterwards.

5. The ``finish`` method constructs the names of the expected output
   Fastqc HTML files, and adds them to the ``out_files`` list which
   was originally initialised within the ``init`` method.

6. The task can explicitly indicate a failure by calling the ``fail``
   method, as in the ``finish`` method above. Raising an exception
   will also implicitly indicate a failure.

7. The ``output`` method returns the ``out_files`` list object.

Another example is a task which does a simple-minded read count on
a set of Fastq files::

    class CountReads(PipelineTask):
        def init(self,fastqs):
            self.counts = dict()
        def setup(self):
            for fq in self.args.fastqs:
                if os.path.splitext(fq)[1] == ".gz":
                   cat = "zcat"
                else:
                   cat = "cat"
                self.add_cmd(
                    PipelineCommandWrapper("Count reads",
                                           "echo","-n",fq,"' '","&&",
                                            cat,fq,"|",
                                            "wc","-l"))
        def finish(self):
            for line in self.stdout.split('\n'):
                if not line:
                    continue
                fq = line.split()[0]
                read_count = int(line.split()[1])/4
                self.counts[fq] = read_count
        def output(self):
            return self.counts

The key features are:

1. The ``init`` method initialises an internal task-specific variable
   which is used in other methods.

2. The standard output from the task is available via the ``stdout``
   property of the instance.

3. The ``finish`` method is implemented to extract the line count data
   from standard output and convert this to a read count which is
   stored against the Fastq name.

A final example is a task which filters out the Fastq files which have
non-zero read counts::

    class FilterEmptyFastqs(PipelineTask):
        def init(self,read_counts):
            self.filtered_fastqs = list()
        def setup(self):
            for fq in self.args.read_counts:
                if self.args.read_counts[fq] > 0:
                     self.filtered_fastqs.append(fq)
        def output(self):
            return self.filtered_fastqs

In this case all the processing is performed by the ``setup`` method;
no commands are defined.

Task outputs
------------

Each task subclass needs to implement its own ``output`` method, which
is used to pass the outputs of the task to other tasks in the pipeline.

In principle the ``output`` method can return anything.

It is important to note that when building a pipeline, the objects
returned by ``output`` are likely to be passed to other tasks before
the task has completed. There are a number of implications:

1. Tasks that receive output from a preceeding task cannot assume that
   those outputs are ready or complete at the point of initialisation.
   It's therefore recommended that tasks don't attempt to use those
   outputs in their 'init' method - all processing should be deferred
   until the 'setup' method (when preceeding tasks will have completed).

2. Outputs from tasks should be passed by object references (e.g.
   via a list or dictionary, which can be updated by the task after
   being passed, or via specialised classes such as ``FileCollector``).

As an example, consider the following task which reverses the order of
a list of items::

    class ReverseList(PipelineTask):
        def init(self,items):
            # Create a list to use as output
            self.reversed = list()
        def setup(self):
            # Generate the reversed list
            for item in self.args.items[::-1]:
                self.reversed.append(item)
        def output(self):
            # Return the list
            return self.reversed

This might be used as follows::

    reverse_list = ReverseList("Reverse order of list",[1,2,3])

Subsequently ``reverse_list.output()`` will return the reference to
the output list, which can then be passed to another task. The list
will be populated when the reverse task runs, at which point the
output will be available to the following tasks.

The ``FileCollector`` class is a specialised class which enables
the collection of files matching a glob-type pattern, and which can
be used as an alternative where appropriate. For example::

    class MakeFiles(PipelineTask):
        def init(self,d,filenames):
            pass
        def setup(self):
            # Create a directory and "touch" the files
            os.mkdir(self.args.d)
            for f in self.args.filenames:
                with open(os.path.join(self.args.d,f) as fp:
                    fp.write()
        def output():
            return FileCollector(self.args.d,"*")

Building and running a pipeline
-------------------------------

An empty pipeline is created by making a ``Pipeline`` instance::

    ppl = Pipeline()

Specific task instances are then created from PipelineTask subclasses,
for example using the tasks defined previously::

    read_counts = CountReads("Count the reads",fastqs)
    filter_empty_fastqs = FilterEmptyFastqs("Filter empty Fastqs",
                                            read_counts.output())
    run_fastqc = RunFastqc("Run Fastqc",
                           filter_empty_fastqs.output())

Note that when instantiating a task, it must be given a name; this
can be any arbitrary text and is intended to help the end user
distinguish between different task instances in the pipeline.

The tasks are then added to the pipeline using the ``add_task``
method::

    ppl.add_task(read_counts)
    ppl.add_task(filter_empty_fastqs,
                 requires=(read_counts,))
    ppl.add_task(run_fastqc,
                 requires=(filter_empty_fastqs,))

The pipeline is then run using the ``run`` method::

    ppl.run(...)

Notes:

1. Tasks will only be executed once any tasks they depend on have
   completed successfully; these are specified via the ``requires``
   argument (which must be a list or tuple of task instances).

2. Tasks that fail (i.e. complete with non-zero exit status) will
   cause the pipeline to halt at that point.

3. The ``run`` method blocks and returns the exit status of the
   pipeline execution.

PipelineCommand versus PipelineCommandWrapper
---------------------------------------------

In the first version of the pipeliner code, commands within the
``setup`` method of ``PipelineTasks`` had to be generated using
subclasses of the ``PipelineCommand`` class.

A simple example to run Fastqc::

    class Fastqc(PipelineCommand):
        def init(self,fastq,out_dir):
            self._fastq = fastq
            self._out_dir = out_dir
        def cmd(self):
            return Command("fastqc",
                           "-o",self._out_dir,
                           self._fastq)

which can then be used within a task, for example::

    class RunFastqc(PipelineTask):
        ...
        def setup(self):
            ...
            for fq in self.args.fastqs:
                self.add_cmd(Fastqc(fq,self.args.out_dir))

as an alternative to the example ``RunFastqc`` example task given
previously, which used the ``PipelineCommandWrapper`` class to
explicitly generate the Fastqc command.

Both approaches are valid. However: using ``PipelineCommand`` is
probably better suited to situations where the same command was used
in more than one distinct tasks. In cases where the command is only
used in one task, using ``PipelineCommandWrapper`` is recommended.

"""

######################################################################
# Imports
######################################################################

import os
import sys
import logging
import shutil
import time
import glob
import uuid
import inspect
import traceback
import string
from collections import Iterator
from cStringIO import StringIO
from bcftbx.utils import mkdir
from bcftbx.utils import AttributeDictionary
from auto_process_ngs.applications import Command

# Module specific logger
logger = logging.getLogger(__name__)

######################################################################
# Pipeline infrastructure constants
######################################################################

ALLOWED_CHARS = string.lowercase + string.digits + "._-"

######################################################################
# Generic pipeline base classes
######################################################################

class FileCollector(Iterator):
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
            self._files = collect_files(self._dirn,self._pattern)
            self._idx = -1
        return len(self._files)
    def next(self):
        if self._files is None:
            self._files = collect_files(self._dirn,self._pattern)
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

    A pipeline consists of a set of tasks (defined by
    instantiating subclasses of PipelineTask) with
    simple dependency relationships (i.e. a task will
    depend on none, one or more other tasks in the
    pipeline to complete before it can be executed).

    Example usage:

    >> p = Pipeline()
    >> t1 = p.add_task(Task1())
    >> t2 = p.add_task(Task2(),requires=(t1,))
    >> ...
    >> p.run()

    Tasks will only run when all requirements have
    completed (or will run immediately if they don't
    have any requirements).
    """
    def __init__(self,name="PIPELINE"):
        """
        Create a new Pipeline instance

        Arguments:
          name (str): optional name for the pipeline
        """
        self._name = str(name)
        self._pending = []
        self._running = []
        self._finished = []

    def report(self,s):
        """
        Internal: report messages from the pipeline
        """
        print "%s [%s] %s" % (time.strftime("%Y-%m-%d %H:%M:%S"),
                              self._name,s)

    def terminate(self):
        """
        Internal: terminate a running pipeline
        """
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

    def add_task(self,task,requires=(),**kws):
        """
        Add a task to the pipeline

        Arguments:
          task (PipelineTask): task instance to add to
            the pipeline
          requires (List): list or tuple of task instances
            which need to complete before this task will
            start
          kws (Dictionary): a dictionary of keyword-value
            pairs which will be passed to the task at
            run time (see the ``run`` method of
            PipelineTask for valid options)
        """
        self._pending.append((task,requires,kws))
        self.report("Adding task '%s'" % task.name())
        if requires:
            for req in requires:
                if req.name() not in [t[0].name() for t in self._pending]:
                    self.report("-> Adding requirement '%s'" % req.name())
                    self.add_task(req)
        return task

    def run(self,working_dir=None,log_dir=None,scripts_dir=None,
            sched=None,default_runner=None):
        """
        Run the tasks in the pipeline

        Arguments:
          working_dir (str): optional path to a working
            directory (defaults to the current directory)
          log_dir (str): path of directory where log files
            will be written to
          scripts_dir (str): path of directory where script
            files will be written to
          sched (SimpleScheduler): a scheduler to use for
            running commands generated by each task
          default_runner (JobRunner): optional default
            job runner to use
        """
        # Execute the pipeline
        self.report("Started")
        # Deal with working directory
        if working_dir is None:
            working_dir = os.getcwd()
        working_dir = os.path.abspath(working_dir)
        self.report("Working directory: %s" % working_dir)
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
                    run_task = reduce(lambda x,y:
                                      x and y.completed
                                      and y.exit_code == 0,
                                      requirements,True)
                if run_task:
                    self.report("started %s" % task.name())
                    if 'runner' not in kws:
                        kws['runner'] = default_runner
                    if 'working_dir' not in kws:
                        kws['working_dir'] = working_dir
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

    This class should be subclassed to implement the 'init',
    'setup', 'finish' (optionally) and 'output' methods.

    The 'add_cmd' method can be used within 'setup' to add one
    or 'PipelineCommand' instances.

    """
    def __init__(self,_name,*args,**kws):
        """
        Create a new PipelineTask instance

        Arguments:
          name (str): an arbitrary user-friendly name for the
            task instance
          args (List): list of arguments to be supplied to
            the subclass (must match those defined in the
            'init' method)
          kws (Dictionary): dictionary of keyword-value pairs
            to be supplied to the subclass (must match those
            defined in the 'init' method)
        """
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
        # Execute the init method
        self.invoke(self.init,self._args,self._kws)

    @property
    def args(self):
        """
        Fetch parameters supplied to the instance
        """
        return AttributeDictionary(**self._callargs)

    @property
    def completed(self):
        """
        Check if the task has completed
        """
        return self._completed

    @property
    def exit_code(self):
        """
        Get the exit code for completed task

        Returns:
          Integer: exit code, or 'None' if task hasn't completed
        """
        if not self.completed:
            return None
        else:
            return self._exit_code

    @property
    def stdout(self):
        """
        Get the standard output from the task

        Returns:
          String: standard output from the task.
        """
        stdout = []
        for f in self._stdout_files:
            with open(f,'r') as fp:
                stdout.append(fp.read())
        return '\n'.join(stdout)

    def name(self):
        """
        Get the name of the task within the pipeline

        Returns:
          String: a name consisting of a 'sanitized' version
            of the supplied name appended with a unique id
            code
        """
        return self._task_name

    def fail(self,exit_code=1,message=None):
        """
        Register the task as failing

        Intended to be invoked from the subclassed 'setup'
        or 'finish' methods, to terminate the task and
        indicate that it has failed.

        Arguments:
          exit_code (int): optional, specifies the exit code
            to return (defaults to 1)
          message (str): optional, error message to report to
            the pipeline user
        """
        if message:
            self.report("failed: %s" % message)
        self.report("failed: exit code set to %s" % exit_code)
        self._completed = True
        self._exit_code = exit_code

    def report(self,s):
        """
        Internal: report messages from the task
        """
        print "%s [Task: %s] %s" % (time.strftime("%Y-%m-%d %H:%M:%S"),
                                    self._name,s)
    def invoke(self,f,args=None,kws=None):
        """
        Internal: invoke arbitrary method on the task

        Arguments:
          f (function): method to invoke (e.g. 'self.init')
          args (list): arguments to invoke function with
          kws (dictionary): keyworded parameters to invoke
            function with
        """
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
        """
        Internal: callback method

        This is a callback method which is invoked when
        scheduled jobs in the task finish

        Arguments:
          name (str): name for the callback
          jobs (list): list of SchedulerJob instances
          sched (SimpleScheduler): scheduler instance
        """
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
        """
        Add a PipelineCommand to the task

        Arguments:
           pipeline_job (PipelineCommand): a PipelineCommand
             instance to be executed by the task when it
             runs
        """
        self._commands.append(pipeline_job)

    def run(self,sched=None,runner=None,working_dir=None,log_dir=None,
            scripts_dir=None,wait_for=(),async=True):
        """
        Run the task

        This method is not normally invoked directly; instead
        it's called by the pipeline that the task has been
        added to.

        Arguments:
          sched (SimpleScheduler): scheduler to submit jobs to
          runner (JobRunner): job runner to use when running
            jobs via the scheduler
          working_dir (str): path to the working directory to use
          log_dir (str): path to the directory to write logs to
          scripts_dir (str): path to the directory to write
            scripts to
          wait_for (list): deprecated: list of scheduler jobs to
            wait for before running jobs from this task
          async (bool): deprecated: if False then block until the
            task has completed
        """
        # Do setup
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
        """
        Internal: terminate the task
        """
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
    """
    def __init__(self,*args,**kws):
        """
        Create a new PipelineCommand instance

        Arguments:
          args (List): list of arguments to be supplied to
            the subclass (must match those defined in the
            'init' method)
          kws (Dictionary): dictionary of keyword-value pairs
            to be supplied to the subclass (must match those
            defined in the 'init' method)
        """
        # Set internal name
        self._name = self.__class__.__name__
        # Invoke the 'init' method
        self.init(*args,**kws)

    def name(self):
        """
        Return a "sanitized" version of the class name
        """
        return sanitize_name(self._name)

    def make_wrapper_script(self,scripts_dir=None,shell="/bin/bash"):
        """
        Generate a uniquely-named wrapper script to run the command

        Arguments:
          scripts_dir (str): path of directory to write
            the wrapper scripts to
          shell (str): shell to use (defaults to '/bin/bash')

        Returns:
          String: name of the wrapper script.
        """
        # Wrap in a script
        if scripts_dir is None:
            scripts_dir = os.getcwd()
        script_file = os.path.join(scripts_dir,"%s.%s.sh" % (self.name(),
                                                             uuid.uuid4()))
        self.cmd().make_wrapper_script(filen=script_file,
                                       shell=shell)
        return script_file

    def init(self):
        """
        Initialise and store parameters

        Must be implemented by the subclass
        """
        raise NotImplementedError("Subclass must implement 'init' method")

    def cmd(self):
        """
        Build the command

        Must be implemented by the subclass and return a
        Command instance
        """
        raise NotImplementedError("Subclass must implement 'cmd' method")

class PipelineCommandWrapper(PipelineCommand):
    """
    Class for constructing program command lines

    This class is based on the PipelineCommand class but
    can be used directly (rather than needing to be
    subclassed).

    For example, to wrap the 'ls' command directly:

    >>> ls_command = PipelineCommandWrapper("List directory",'ls',dirn)

    It is also possible to extend the command line
    using the 'add_args' method, for example:

    >>> ls_command = PipelineCommandWrapper("List directory",'ls')
    >>> ls.command.add_args(dirn)
    """
    def __init__(self,name,*args):
        """
        Create a new PipelineCommandWrapper instance

        Arguments:
          name (str): arbitrary name for the command
          args  (List): initial list of arguments making
            up the command
        """
        PipelineCommand.__init__(self,*args)
        self._name = str(name)
        self._cmd = None
        if args:
            self._cmd = Command(*args)

    def add_args(self,*args):
        """
        Add additional arguments to extend the command being built

        Arguments:
          args  (List): one or more arguments to append to
            the command
        """
        if self._cmd is None:
            self._cmd = Command(*args)
        else:
            self._cmd.add_args(*args)

    def init(self,*args):
        """
        Internal: dummy init which does nothing
        """
        pass

    def cmd(self):
        """
        Internal: implement the 'cmd' method
        """
        return self._cmd

######################################################################
# Generic pipeline functions
######################################################################

def sanitize_name(s):
    """
    Convert string to lowercase and replace special characters

    Arguments:
      s (str): string to sanitize

    Returns:
      String: sanitized string.
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

def collect_files(dirn,pattern):
    """
    Return names of files in a directory which match a glob pattern

    Arguments:
      dirn (str): path to a directory containing the files
      pattern (str): a glob pattern to match

    Returns:
      List: list of matching files
    """
    return sorted(glob.glob(os.path.join(os.path.abspath(dirn),pattern)))
