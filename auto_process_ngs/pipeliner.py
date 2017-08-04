#!/usr/bin/env python
#
#     pipeliner.py: utilities for building simple pipelines of tasks
#     Copyright (C) University of Manchester 2017 Peter Briggs
#
"""
pipeliner.py

Module providing utility classes and functions for building simple
'pipelines' of tasks.

The core classes are:

- Pipeline: class for building and executing pipelines
- PipelineTask: class for defining pipeline tasks
- PipelineCommand: class for defining commands that can be used in tasks

Additional supporting classes:

- PipelineCommandWrapper: shortcut alternative to PipelineCommand
- FileCollection: returning collections of files based on glob patterns

There are some underlying classes and functions that are intended for
internal use:

- Capturing: capture stdout from a Python function
- sanitize_name: clean up task and command names for use in pipeline
- collect_files: collect files based on glob patterns
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
    def __init__(self,name,*args):
        PipelineCommand.__init__(self,*args)
        self._name = str(name)
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
