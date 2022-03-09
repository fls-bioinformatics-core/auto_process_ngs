#######################################################################
# Tests for pipeliner.py module
#######################################################################

import unittest
import tempfile
import shutil
import time
import os
import io
import getpass
import platform
import cloudpickle
from builtins import range
import auto_process_ngs.envmod as envmod
from auto_process_ngs.mock import MockConda
from auto_process_ngs.simple_scheduler import SimpleScheduler
from auto_process_ngs.command import Command
from auto_process_ngs.pipeliner import Pipeline
from auto_process_ngs.pipeliner import PipelineTask
from auto_process_ngs.pipeliner import PipelineFunctionTask
from auto_process_ngs.pipeliner import PipelineCommand
from auto_process_ngs.pipeliner import PipelineCommandWrapper
from auto_process_ngs.pipeliner import PipelineScriptWrapper
from auto_process_ngs.pipeliner import PipelineParam
from auto_process_ngs.pipeliner import PipelineFailure
from auto_process_ngs.pipeliner import FileCollector
from auto_process_ngs.pipeliner import Dispatcher
from auto_process_ngs.pipeliner import BaseParam
from auto_process_ngs.pipeliner import ListParam
from auto_process_ngs.pipeliner import PathJoinParam
from auto_process_ngs.pipeliner import PathExistsParam
from auto_process_ngs.pipeliner import FunctionParam
from auto_process_ngs.pipeliner import PipelineError
from auto_process_ngs.pipeliner import resolve_parameter
from bcftbx.JobRunner import SimpleJobRunner
from bcftbx.Pipeline import Job

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

# Check if Job instance has 'err' attribute
JOB_HAS_ERR = hasattr(
    Job(SimpleJobRunner(),'dummy','.','script',[]),
    'err')

# Helpers

class _Mock:
    """
    Helper class for mocking executables for testing
    """

    @staticmethod
    def fastqc(bin_dir):
        """
        Make a mock FastQC executable for testing

        Arguments:
          bin_dir (str): path to directory to put mock
           'fastqc' executable into (must already exist)

        Returns:
           String: path to mock 'fastqc' executable.
        """
        fastqc_ = os.path.join(bin_dir,"fastqc")
        with io.open(fastqc_,'wt') as fp:
            fp.write(u"""#!/bin/bash
echo $1
exit 0
""")
        os.chmod(os.path.join(bin_dir,"fastqc"),0o775)
        return fastqc_

# Unit tests

class TestPipeline(unittest.TestCase):

    def setUp(self):
        # Placeholder for scheduler instance
        self.sched = None
        # Make a temporary working dir
        self.working_dir = tempfile.mkdtemp(
            suffix='TestPipeline')

    def _get_scheduler(self):
        # Set up a scheduler
        self.sched = SimpleScheduler(poll_interval=0.01)
        self.sched.start()
        return self.sched

    def tearDown(self):
        # Stop the scheduler
        if self.sched is not None:
            self.sched.stop()
        # Remove temp dir
        if os.path.exists(self.working_dir):
            shutil.rmtree(self.working_dir)

    def test_simple_pipeline(self):
        """
        Pipeline: define and run a simple pipeline
        """
        # Define a reusable task
        # Appends item to a list
        class Append(PipelineTask):
            def init(self,l,s):
                self.add_output('list',list())
            def setup(self):
                for item in self.args.l:
                    self.output.list.append(item)
                self.output.list.append(self.args.s)
        # Build the pipeline
        ppl = Pipeline()
        task1 = Append("Append 1",(),"item1")
        task2 = Append("Append 2",task1.output.list,"item2")
        ppl.add_task(task2,requires=(task1,))
        # Run the pipeline
        exit_status = ppl.run(working_dir=self.working_dir,
                              poll_interval=0.1)
        # Check the outputs
        self.assertEqual(exit_status,0)
        self.assertEqual(task1.output.list,["item1"])
        self.assertEqual(task2.output.list,["item1","item2"])

    def test_pipeline_with_commands(self):
        """
        Pipeline: define and run pipeline with commands
        """
        # Define a task
        # Echoes/appends text to a file
        class Echo(PipelineTask):
            def init(self,f,s):
                self.add_output('file',f)
            def setup(self):
                self.add_cmd(
                    PipelineCommandWrapper(
                        "Echo text to file",
                        "echo",self.args.s,
                        ">>",self.args.f))
        # Build the pipeline
        ppl = Pipeline()
        out_file = os.path.join(self.working_dir,"out.txt")
        task1 = Echo("Write item1",out_file,"item1")
        task2 = Echo("Write item2",task1.output.file,"item2")
        ppl.add_task(task2,requires=(task1,))
        # Run the pipeline
        exit_status = ppl.run(working_dir=self.working_dir,
                              poll_interval=0.1)
        # Check the outputs
        self.assertEqual(exit_status,0)
        self.assertTrue(os.path.exists(out_file))
        with open(out_file,'rt') as fp:
            self.assertEqual(fp.read(),"item1\nitem2\n")

    def test_pipeline_with_multiple_commands(self):
        """
        Pipeline: define and run pipeline with multiple commands
        """
        # Define a task
        # Echoes/appends text to a file
        class EchoMany(PipelineTask):
            def init(self,*s):
                self.add_output('files',list())
            def setup(self):
                for f,s in self.args.s:
                    self.add_cmd(
                        PipelineCommandWrapper(
                            "Echo text to file",
                            "echo",s,
                            ">>",f))
            def finish(self):
                for f,s in self.args.s:
                    self.output.files.append(f)
        # Build the pipeline
        out_files = [os.path.join(self.working_dir,f)
                     for f in ("out1.txt","out2.txt")]
        inputs = [(f,"item") for f in out_files]
        ppl = Pipeline()
        task = EchoMany("Write items",*inputs)
        ppl.add_task(task)
        # Run the pipeline
        exit_status = ppl.run(working_dir=self.working_dir,
                              poll_interval=0.1)
        # Check the outputs
        self.assertEqual(exit_status,0)
        for out_file in out_files:
            self.assertTrue(os.path.exists(out_file))
            with open(out_file,'rt') as fp:
                self.assertEqual(fp.read(),"item\n")

    def test_pipeline_with_batched_commands(self):
        """
        Pipeline: define and run pipeline with batched commands
        """
        # Define a task
        # Echoes/appends text to a file
        class EchoMany(PipelineTask):
            def init(self,*s):
                self.add_output('files',list())
            def setup(self):
                for f,s in self.args.s:
                    self.add_cmd(
                        PipelineCommandWrapper(
                            "Echo text to file",
                            "echo",s,
                            ">>",f))
            def finish(self):
                for f,s in self.args.s:
                    self.output.files.append(f)
        # Build the pipeline
        ppl = Pipeline()
        out_files = [os.path.join(self.working_dir,f)
                     for f in ("out1.txt",
                               "out2.txt",
                               "out3.txt",
                               "out4.txt",
                               "out5.txt")]
        inputs = [(f,"item") for f in out_files]
        task = EchoMany("Write items",*inputs)
        ppl.add_task(task)
        # Run the pipeline
        exit_status = ppl.run(working_dir=self.working_dir,
                              poll_interval=0.1,
                              batch_size=2)
        # Check the outputs
        self.assertEqual(exit_status,0)
        for out_file in out_files:
            self.assertTrue(os.path.exists(out_file),
                            "Missing %s" % out_file)
            with open(out_file,'rt') as fp:
                self.assertEqual(fp.read(),"item\n")

    def test_pipeline_with_batch_limit(self):
        """
        Pipeline: define and run pipeline with batch limit set
        """
        # Define a task
        # Echoes/appends text to a file
        class EchoMany(PipelineTask):
            def init(self,*s):
                self.add_output('files',list())
            def setup(self):
                for f,s in self.args.s:
                    self.add_cmd(
                        PipelineCommandWrapper(
                            "Echo text to file",
                            "echo",s,
                            ">>",f))
            def finish(self):
                for f,s in self.args.s:
                    self.output.files.append(f)
        # Build the pipeline
        ppl = Pipeline()
        out_files = [os.path.join(self.working_dir,f)
                     for f in ("out1.txt",
                               "out2.txt",
                               "out3.txt",
                               "out4.txt",
                               "out5.txt")]
        inputs = [(f,"item") for f in out_files]
        task = EchoMany("Write items",*inputs)
        ppl.add_task(task)
        # Run the pipeline
        exit_status = ppl.run(working_dir=self.working_dir,
                              poll_interval=0.1,
                              batch_limit=3)
        # Check the outputs
        self.assertEqual(exit_status,0)
        for out_file in out_files:
            self.assertTrue(os.path.exists(out_file),
                            "Missing %s" % out_file)
            with open(out_file,'rt') as fp:
                self.assertEqual(fp.read(),"item\n")

    def test_pipeline_sets_pipelineparam_as_task_input(self):
        """
        Pipeline: handle PipelineParam passed as task input
        """
        # Define a task
        # Echoes/appends text to a file
        class Echo(PipelineTask):
            def init(self,f,s):
                self.add_output('file',f)
            def setup(self):
                self.add_cmd(
                    PipelineCommandWrapper(
                        "Echo text to file",
                        "echo",self.args.s,
                        ">>",self.args.f))
        # Create a pipeline param instance
        s = PipelineParam("hello")
        self.assertEqual(s.value,"hello")
        # Build the pipeline
        out_file = os.path.join(self.working_dir,"out.txt")
        ppl = Pipeline()
        task = Echo("Write item1",out_file,s)
        ppl.add_task(task)
        # Update the pipeline param
        s.set("goodbye")
        # Run the pipeline
        exit_status = ppl.run(working_dir=self.working_dir,
                              poll_interval=0.1)
        # Check the outputs
        self.assertEqual(exit_status,0)
        self.assertTrue(os.path.exists(out_file))
        with open(out_file,'rt') as fp:
            self.assertEqual(fp.read(),"goodbye\n")

    def test_pipeline_with_external_scheduler(self):
        """
        Pipeline: run pipeline using user-defined scheduler
        """
        # Define a task
        # Echoes/appends text to a file
        class Echo(PipelineTask):
            def init(self,f,s):
                self.add_output('file',f)
            def setup(self):
                self.add_cmd(
                    PipelineCommandWrapper(
                        "Echo text to file",
                        "echo",self.args.s,
                        ">>",self.args.f))
        # Build the pipeline
        ppl = Pipeline()
        out_file = os.path.join(self.working_dir,"out.txt")
        task1 = Echo("Write item1",out_file,"item1")
        task2 = Echo("Write item2",task1.output.file,"item2")
        ppl.add_task(task2,requires=(task1,))
        # Get a scheduler
        self._get_scheduler()
        # Run the pipeline
        exit_status = ppl.run(sched=self.sched,
                              working_dir=self.working_dir,
                              poll_interval=0.1)
        # Check the outputs
        self.assertEqual(exit_status,0)
        self.assertTrue(os.path.exists(out_file))
        with open(out_file,'rt') as fp:
            self.assertEqual(fp.read(),"item1\nitem2\n")

    def test_pipeline_task_working_dirs_isolation(self):
        """
        Pipeline: check pipeline with isolation of task working directories
        """
        # Define tasks
        # Echoes/appends working dir to a file via shell command
        class EchoWorkDir(PipelineTask):
            def init(self,f):
                self.add_output('file',f)
            def setup(self):
                self.add_cmd(
                    PipelineCommandWrapper(
                        "Echo working dir to file",
                        "echo","$(pwd)",
                        ">>",self.args.f))
        # Build the pipeline
        out_file = os.path.join(self.working_dir,"dirs.txt")
        ppl = Pipeline()
        task1 = EchoWorkDir("First task work dir",out_file)
        task2 = EchoWorkDir("Second task work dir",task1.output.file)
        ppl.add_task(task2,requires=(task1,))
        # Run the pipeline
        exit_status = ppl.run(working_dir=self.working_dir,
                              isolate_tasks=True,
                              poll_interval=0.1)
        # Check the outputs
        self.assertEqual(exit_status,0)
        self.assertTrue(os.path.exists(out_file))
        with open(out_file,'rt') as fp:
            work_dirs = fp.read().strip().split()
            # Directories for each task should be different
            self.assertNotEqual(work_dirs[0],work_dirs[1])

    def test_pipeline_task_working_dirs_no_isolation(self):
        """
        Pipeline: check pipeline with no isolation of task working directories
        """
        # Define tasks
        # Echoes/appends working dir to a file via shell command
        class EchoWorkDir(PipelineTask):
            def init(self,f):
                self.add_output('file',f)
            def setup(self):
                self.add_cmd(
                    PipelineCommandWrapper(
                        "Echo working dir to file",
                        "echo","$(pwd)",
                        ">>",self.args.f))
        # Build the pipeline
        out_file = os.path.join(self.working_dir,"dirs.txt")
        ppl = Pipeline()
        task1 = EchoWorkDir("First task work dir",out_file)
        task2 = EchoWorkDir("Second task work dir",task1.output.file)
        ppl.add_task(task2,requires=(task1,))
        # Run the pipeline
        exit_status = ppl.run(working_dir=self.working_dir,
                              isolate_tasks=False,
                              poll_interval=0.1)
        # Check the outputs
        self.assertEqual(exit_status,0)
        self.assertTrue(os.path.exists(out_file))
        with open(out_file,'rt') as fp:
            for work_dir in fp.read().strip().split():
                # Directories should be top-level working dir
                self.assertEqual(work_dir,self.working_dir)

    def test_pipeline_task_working_dirs_isolation_set_work_dir(self):
        """
        Pipeline: check pipeline sets top-level tasks work dir (with isolation)
        """
        # Define tasks
        # Echoes/appends working dir to a file via shell command
        class EchoWorkDir(PipelineTask):
            def init(self,f):
                self.add_output('file',f)
            def setup(self):
                self.add_cmd(
                    PipelineCommandWrapper(
                        "Echo working dir to file",
                        "echo","$(pwd)",
                        ">>",self.args.f))
        # Build the pipeline
        out_file = os.path.join(self.working_dir,"dirs.txt")
        ppl = Pipeline()
        task1 = EchoWorkDir("First task work dir",out_file)
        task2 = EchoWorkDir("Second task work dir",task1.output.file)
        ppl.add_task(task2,requires=(task1,))
        # Run the pipeline
        tasks_work_dir = os.path.join(self.working_dir,
                                      "tasks_work_dirs")
        exit_status = ppl.run(working_dir=self.working_dir,
                              tasks_work_dir="tasks_work_dirs",
                              isolate_tasks=True,
                              poll_interval=0.1)
        # Check the outputs
        self.assertEqual(exit_status,0)
        self.assertTrue(os.path.exists(out_file))
        with open(out_file,'rt') as fp:
            work_dirs = fp.read().strip().split()
            # Directories for each task should be different
            self.assertNotEqual(work_dirs[0],work_dirs[1])
            # Task dirs should be subdirs of tasks_work_dir
            for work_dir in work_dirs:
                self.assertEqual(os.path.dirname(work_dir),
                                 tasks_work_dir)

    def test_pipeline_task_working_dirs_no_isolation_set_work_dir(self):
        """
        Pipeline: check pipeline sets top-level tasks work dir (no isolation)
        """
        # Define tasks
        # Echoes/appends working dir to a file via shell command
        class EchoWorkDir(PipelineTask):
            def init(self,f):
                self.add_output('file',f)
            def setup(self):
                self.add_cmd(
                    PipelineCommandWrapper(
                        "Echo working dir to file",
                        "echo","$(pwd)",
                        ">>",self.args.f))
        # Build the pipeline
        out_file = os.path.join(self.working_dir,"dirs.txt")
        ppl = Pipeline()
        task1 = EchoWorkDir("First task work dir",out_file)
        task2 = EchoWorkDir("Second task work dir",task1.output.file)
        ppl.add_task(task2,requires=(task1,))
        # Run the pipeline
        tasks_work_dir = os.path.join(self.working_dir,
                                      "tasks_work_dirs")
        exit_status = ppl.run(working_dir=self.working_dir,
                              tasks_work_dir="tasks_work_dirs",
                              isolate_tasks=False,
                              poll_interval=0.1)
        # Check the outputs
        self.assertEqual(exit_status,0)
        self.assertTrue(os.path.exists(out_file))
        with open(out_file,'rt') as fp:
            for work_dir in fp.read().strip().split():
                # Directories should be top-level tasks working dir
                self.assertEqual(work_dir,tasks_work_dir)

    def test_pipeline_stops_on_task_failure(self):
        """
        Pipeline: immediate exit on task failure in 'immediate' mode
        """
        # Define a reusable task
        # Appends item to a list
        class Append(PipelineTask):
            def init(self,l,s):
                self.add_output('list',list())
            def setup(self):
                for item in self.args.l:
                    self.output.list.append(item)
                self.output.list.append(self.args.s)
        # Define a version of the 'append' task that
        # always fails
        class Failure(Append):
            def setup(self):
                self.fail(message="Automatic fail")
        # Build a failing pipeline
        ppl = Pipeline()
        task1 = Append("Append 1",(),"1")
        task2 = Append("Append 2",(),"2")
        task3 = Failure("Failing append 3",(),"3")
        task2_1 = Append("Append 2_1",task2.output.list,"2_1")
        task2_2 = Append("Append 2_2",task2_1.output.list,"2_2")
        task3_1 = Append("Append 3_1",task3.output.list,"3_1")
        task3_2 = Append("Append 3_2",task3_1.output.list,"3_2")
        ppl.add_task(task2,requires=(task1,))
        ppl.add_task(task2_1,requires=(task2,))
        ppl.add_task(task2_2,requires=(task2_1,))
        ppl.add_task(task3,requires=(task1,))
        ppl.add_task(task3_1,requires=(task3,))
        ppl.add_task(task3_2,requires=(task3_1,))
        # Run the pipeline
        exit_status = ppl.run(working_dir=self.working_dir,
                              poll_interval=0.1,
                              exit_on_failure=PipelineFailure.IMMEDIATE)
        # Check the outputs
        self.assertEqual(exit_status,1)
        self.assertEqual(task1.output.list,["1"])
        self.assertEqual(task2.output.list,["2"])
        self.assertEqual(task2_1.output.list,[])
        self.assertEqual(task2_2.output.list,[])
        self.assertEqual(task3_1.output.list,[])
        self.assertEqual(task3_2.output.list,[])
        # Check the exit codes
        self.assertEqual(task1.exit_code,0)
        self.assertEqual(task2.exit_code,0)
        self.assertEqual(task3.exit_code,1)
        self.assertEqual(task2_1.exit_code,None)
        self.assertEqual(task2_2.exit_code,None)
        self.assertEqual(task3_1.exit_code,None)
        self.assertEqual(task3_2.exit_code,None)
        self.assertEqual(exit_status,1)

    def test_pipeline_defers_task_failure(self):
        """
        Pipeline: deferred exit on task failure in 'deferred' mode
        """
        # Define a reusable task
        # Appends item to a list
        class Append(PipelineTask):
            def init(self,l,s):
                self.add_output('list',list())
            def setup(self):
                for item in self.args.l:
                    self.output.list.append(item)
                self.output.list.append(self.args.s)
        # Define a version of the 'append' task that
        # always fails
        class Failure(Append):
            def setup(self):
                self.fail(message="Automatic fail")
        # Build a failing pipeline
        ppl = Pipeline()
        task1 = Append("Append 1",(),"1")
        task2 = Append("Append 2",(),"2")
        task3 = Failure("Failing append 3",(),"3")
        task2_1 = Append("Append 2_1",task2.output.list,"2_1")
        task2_2 = Append("Append 2_2",task2_1.output.list,"2_2")
        task3_1 = Append("Append 3_1",task3.output.list,"3_1")
        task3_2 = Append("Append 3_2",task3_1.output.list,"3_2")
        ppl.add_task(task2,requires=(task1,))
        ppl.add_task(task2_1,requires=(task2,))
        ppl.add_task(task2_2,requires=(task2_1,))
        ppl.add_task(task3,requires=(task1,))
        ppl.add_task(task3_1,requires=(task3,))
        ppl.add_task(task3_2,requires=(task3_1,))
        # Run the pipeline
        exit_status = ppl.run(working_dir=self.working_dir,
                              poll_interval=0.1,
                              exit_on_failure=PipelineFailure.DEFERRED)
        # Check the outputs
        self.assertEqual(exit_status,1)
        self.assertEqual(task1.output.list,["1"])
        self.assertEqual(task2.output.list,["2"])
        self.assertEqual(task2_1.output.list,["2","2_1"])
        self.assertEqual(task2_2.output.list,["2","2_1","2_2"])
        self.assertEqual(task3_1.output.list,[])
        self.assertEqual(task3_2.output.list,[])
        # Check the exit codes
        self.assertEqual(task1.exit_code,0)
        self.assertEqual(task2.exit_code,0)
        self.assertEqual(task3.exit_code,1)
        self.assertEqual(task2_1.exit_code,0)
        self.assertEqual(task2_2.exit_code,0)
        self.assertEqual(task3_1.exit_code,None)
        self.assertEqual(task3_2.exit_code,None)

    def test_pipeline_built_in_parameters(self):
        """
        Pipeline: test the built-in parameters are set
        """
        # Define a task to output values passed in
        class OutputValues(PipelineTask):
            def init(self,working_dir,batch_size,verbose):
                self.add_output('builtins',dict())
            def setup(self):
                self.output.builtins['working_dir'] = self.args.working_dir
                self.output.builtins['batch_size'] = self.args.batch_size
                self.output.builtins['verbose'] = self.args.verbose
        # Make a pipeline
        ppl = Pipeline()
        task = OutputValues("Output the built-in parameter values",
                            working_dir=ppl.params.WORKING_DIR,
                            batch_size=ppl.params.BATCH_SIZE,
                            verbose=ppl.params.VERBOSE)
        ppl.add_task(task)
        self.assertTrue("WORKING_DIR" in ppl.params)
        self.assertTrue("BATCH_SIZE" in ppl.params)
        self.assertTrue("VERBOSE" in ppl.params)
        # Run the pipeline
        exit_status = ppl.run(working_dir=self.working_dir,
                              batch_size=100,
                              verbose=True)
        self.assertEqual(task.output.builtins['working_dir'],
                         self.working_dir)
        self.assertEqual(task.output.builtins['batch_size'],100)
        self.assertEqual(task.output.builtins['verbose'],True)

    def test_pipeline_add_param(self):
        """
        Pipeline: test the 'add_param' method
        """
        # Make an empty pipeline
        ppl = Pipeline()
        # Add a parameter
        self.assertFalse('ncores' in ppl.params)
        ppl.add_param('ncores',value=1,type=int)
        self.assertTrue('ncores' in ppl.params)
        self.assertTrue(isinstance(ppl.params.ncores,
                                   PipelineParam))
        self.assertEqual(ppl.params.ncores.value,1)

    def test_pipeline_add_param_twice_raises_exception(self):
        """
        Pipeline: 'add_param' raises exception when parameter is added twice
        """
        # Make an empty pipeline
        ppl = Pipeline()
        # Add a parameter
        self.assertFalse('ncores' in ppl.params)
        ppl.add_param('ncores',value=1,type=int)
        self.assertTrue('ncores' in ppl.params)
        self.assertRaises(KeyError,
                          ppl.add_param,
                          'ncores')

    def test_pipeline_param_passed_at_runtime(self):
        """
        Pipeline: check parameter is passed at runtime
        """
        # Define task to echoes/append text to a file
        class Echo(PipelineTask):
            def init(self,f,s):
                self.add_output('file',f)
            def setup(self):
                self.add_cmd(
                    PipelineCommandWrapper(
                        "Echo text to file",
                        "echo",self.args.s,
                        ">>",self.args.f))
        # Build the pipeline with a parameter
        ppl = Pipeline()
        ppl.add_param('out_file')
        task1 = Echo("Echo item1",
                     f=ppl.params.out_file,
                     s="item1")
        task2 = Echo("Echo item2",
                     f=ppl.params.out_file,
                     s="item2")
        ppl.add_task(task2,requires=(task1,))
        # Check the parameter before setting a value
        self.assertEqual(ppl.params.out_file.value,None)
        # Run the pipeline setting the output file
        out_file = os.path.join(self.working_dir,"out.txt")
        exit_status = ppl.run(working_dir=self.working_dir,
                              poll_interval=0.1,
                              params={ 'out_file': out_file, })
        # Check the outputs
        self.assertEqual(exit_status,0)
        self.assertTrue(os.path.exists(out_file))
        with open(out_file,'rt') as fp:
            self.assertEqual(fp.read(),"item1\nitem2\n")

    def test_pipeline_add_runner(self):
        """
        Pipeline: test the 'add_runner' method
        """
        # Make an empty pipeline
        ppl = Pipeline()
        # Add a runner definition
        self.assertFalse('test_runner' in ppl.runners)
        ppl.add_runner('test_runner')
        print(ppl.runners)
        self.assertTrue('test_runner' in ppl.runners)
        self.assertTrue(isinstance(ppl.runners['test_runner'],
                                   PipelineParam))
        self.assertTrue(isinstance(ppl.runners['test_runner'].value,
                                   SimpleJobRunner))

    def test_pipeline_add_runner_twice_raises_exception(self):
        """
        Pipeline: 'add_runner' raises exception when runner is added twice
        """
        # Make an empty pipeline
        ppl = Pipeline()
        # Add a runner definition
        self.assertFalse('test_runner' in ppl.runners)
        ppl.add_runner('test_runner')
        self.assertTrue('test_runner' in ppl.runners)
        self.assertRaises(KeyError,
                          ppl.add_runner,
                          'test_runner')

    def test_pipeline_runner_passed_at_runtime(self):
        """
        Pipeline: check runner is passed at runtime
        """
        # Define task to echoes/append text to a file
        class Echo(PipelineTask):
            def init(self,f,s):
                self.add_output('file',f)
            def setup(self):
                self.add_cmd(
                    PipelineCommandWrapper(
                        "Echo text to file",
                        "echo",self.args.s,
                        ">>",self.args.f))
        # Define custom runner for testing
        class TestJobRunner(SimpleJobRunner):
            def __init__(self,*args,**kws):
                SimpleJobRunner.__init__(self,*args,**kws)
        # Build the pipeline
        out_file = os.path.join(self.working_dir,"out.txt")
        ppl = Pipeline()
        ppl.add_runner("test_runner")
        task1 = Echo("Echo item1",
                     f=out_file,
                     s="item1")
        ppl.add_task(task1,
                     runner=ppl.runners["test_runner"])
        # Run the pipeline setting the runner
        exit_status = ppl.run(working_dir=self.working_dir,
                              poll_interval=0.1,
                              runners={ 'test_runner':
                                        TestJobRunner(), })
        # Check the outputs
        self.assertEqual(exit_status,0)
        self.assertTrue(os.path.exists(out_file))
        with open(out_file,'r') as fp:
            self.assertEqual(fp.read(),"item1\n")

    def test_pipeline_runner_set_default_runner_at_runtime(self):
        """
        Pipeline: check setting default runner at runtime
        """
        # Define task to echoes/append text to a file
        class Echo(PipelineTask):
            def init(self,f,s):
                self.add_output('file',f)
            def setup(self):
                self.add_cmd(
                    PipelineCommandWrapper(
                        "Echo text to file",
                        "echo",self.args.s,
                        ">>",self.args.f))
        # Define custom runners for testing
        class DefaultJobRunner(SimpleJobRunner):
            def __init__(self,*args,**kws):
                SimpleJobRunner.__init__(self,*args,**kws)
        class TestJobRunner(SimpleJobRunner):
            def __init__(self,*args,**kws):
                SimpleJobRunner.__init__(self,*args,**kws)
        # Build the pipeline
        out_file = os.path.join(self.working_dir,"out.txt")
        ppl = Pipeline()
        ppl.add_runner("test_runner1")
        ppl.add_runner("test_runner2")
        task1 = Echo("Echo item1",
                     f=out_file,
                     s="item1")
        task2 = Echo("Echo item2",
                     f=out_file,
                     s="item2")
        ppl.add_task(task1,
                     runner=ppl.runners["test_runner1"])
        ppl.add_task(task2,
                     requires=(task1,),
                     runner=ppl.runners["test_runner2"])
        # Run the pipeline setting the default runner
        # and only one of the two runners
        exit_status = ppl.run(working_dir=self.working_dir,
                              poll_interval=0.1,
                              default_runner=DefaultJobRunner(),
                              runners={ 'test_runner1':
                                        TestJobRunner(), })
        # Check the outputs
        self.assertEqual(exit_status,0)
        self.assertTrue(os.path.exists(out_file))
        with open(out_file,'rt') as fp:
            self.assertEqual(fp.read(),"item1\nitem2\n")

    @unittest.skipIf(not envmod.__ENVMODULES__,
                     "Environment modules not available")
    def test_pipeline_with_envmodules(self):
        """
        Pipeline: define and run pipeline with environment modules
        """
        # Set up mock Fastqc
        bin_dir = os.path.join(self.working_dir,"apps","bin")
        os.mkdir(os.path.join(self.working_dir,"apps"))
        os.mkdir(bin_dir)
        _Mock.fastqc(bin_dir)
        # Set up mock environment module
        modules_dir = os.path.join(self.working_dir,"modulefiles")
        os.mkdir(modules_dir)
        modules = "apps/fastqc/1.0"
        os.mkdir(os.path.join(modules_dir,"apps"))
        os.mkdir(os.path.join(modules_dir,"apps","fastqc"))
        with io.open(os.path.join(modules_dir,"apps","fastqc","1.0"),'wt') \
             as fp:
            fp.write(u"""#%%Module1.0
prepend-path PATH %s
""" % bin_dir)
        os.environ['MODULEPATH'] = modules_dir
        # Define a task
        class RunFastqc(PipelineTask):
            def init(self,*files):
                self.add_output('files',list())
            def setup(self):
                for f in self.args.files:
                    self.add_cmd(
                        PipelineCommandWrapper(
                            "Run fastqc for %s" % f,
                            "fastqc",f))
            def finish(self):
                for f in self.args.files:
                    self.output.files.append(f)
        # Build the pipeline
        ppl = Pipeline()
        ppl.add_envmodules("fastqc")
        task = RunFastqc("Run Fastqc",
                         "sample1.fastq","sample2.fastq")
        ppl.add_task(task,
                     envmodules=ppl.envmodules["fastqc"])
        # Run the pipeline
        exit_status = ppl.run(working_dir=self.working_dir,
                              envmodules={
                                  'fastqc': modules,
                              },
                              poll_interval=0.1,
                              verbose=True)
        # Check the outputs
        self.assertEqual(exit_status,0)

    @unittest.skipIf(not envmod.__ENVMODULES__,
                     "Environment modules not available")
    def test_pipeline_fails_with_missing_modules_environment(self):
        """
        Pipeline: check pipeline fails for missing 'modules' environment
        """
        # Missing module file
        modules = "apps/fastqc/1.0"
        # Define a task
        class RunFastqc(PipelineTask):
            def init(self,*files):
                self.add_output('files',list())
            def setup(self):
                for f in self.args.files:
                    self.add_cmd(
                        PipelineCommandWrapper(
                            "Run fastqc for %s" % f,
                            "fastqc",f))
            def finish(self):
                for f in self.args.files:
                    self.output.files.append(f)
        # Build the pipeline
        ppl = Pipeline()
        ppl.add_envmodules("fastqc")
        task = RunFastqc("Run Fastqc",
                         "sample1.fastq","sample2.fastq")
        ppl.add_task(task,
                     envmodules=ppl.envmodules["fastqc"])
        # Run the pipeline
        exit_status = ppl.run(working_dir=self.working_dir,
                              envmodules={
                                  'fastqc': modules,
                              },
                              poll_interval=0.1,
                              verbose=True)
        # Check the pipeline failed (non-zero exit)
        self.assertNotEqual(exit_status,0)

    @unittest.skipIf(not envmod.__ENVMODULES__,
                     "Environment modules not available")
    def test_pipeline_add_envmodules_twice_raises_exception(self):
        """
        Pipeline: 'add_envmodules' raises exception when environment is added twice
        """
        # Make an empty pipeline
        ppl = Pipeline()
        # Add a environment definition
        self.assertFalse('test_env' in ppl.envmodules)
        ppl.add_envmodules('test_env')
        self.assertTrue('test_env' in ppl.envmodules)
        self.assertRaises(KeyError,
                          ppl.add_envmodules,
                          'test_env')

    def test_pipeline_with_conda(self):
        """
        Pipeline: define and run pipeline with conda dependency resolution
        """
        # Set up mock conda installation
        conda_dir = MockConda.create(os.path.join(self.working_dir,"conda"))
        conda_ = os.path.join(conda_dir,"bin","conda")
        # Define a task
        class RunFastqc(PipelineTask):
            def init(self,*files):
                self.conda("fastqc=0.11.3")
                self.add_output('files',list())
            def setup(self):
                for f in self.args.files:
                    self.add_cmd(
                        PipelineCommandWrapper(
                            "Run fastqc for %s" % f,
                            "fastqc",f))
            def finish(self):
                for f in self.args.files:
                    self.output.files.append(f)
        # Build the pipeline
        ppl = Pipeline()
        task = RunFastqc("Run Fastqc",
                         "sample1.fastq","sample2.fastq")
        ppl.add_task(task)
        # Run the pipeline
        exit_status = ppl.run(working_dir=self.working_dir,
                              enable_conda=True,
                              conda=conda_,
                              poll_interval=0.1,
                              verbose=True)
        # Check the outputs
        self.assertEqual(exit_status,0)

    def test_pipeline_with_conda_existing_environment(self):
        """
        Pipeline: run with existing conda environment
        """
        # Set up mock conda installation
        conda_dir = MockConda.create(os.path.join(self.working_dir,"conda"))
        conda_ = os.path.join(conda_dir,"bin","conda")
        # Define a task
        class RunFastqc(PipelineTask):
            def init(self,*files):
                self.conda("fastqc=0.11.3")
                self.add_output('files',list())
            def setup(self):
                for f in self.args.files:
                    self.add_cmd(
                        PipelineCommandWrapper(
                            "Run fastqc for %s" % f,
                            "fastqc",f))
            def finish(self):
                for f in self.args.files:
                    self.output.files.append(f)
        # Mock an existing conda environment (with a mock
        # 'fastqc' exe) for the task
        env_path = os.path.join(self.working_dir,
                                "__conda",
                                "envs",
                                "fastqc@0.11.3")
        os.makedirs(env_path)
        with open(os.path.join(env_path,"fastqc"),'wt') as fp:
            fp.write("#!/bin/bash\necho $1\exit 0\n")
        os.chmod(os.path.join(env_path,"fastqc"),0o755)
        # Build the pipeline
        ppl = Pipeline()
        task = RunFastqc("Run Fastqc",
                         "sample1.fastq","sample2.fastq")
        ppl.add_task(task)
        # Run the pipeline
        exit_status = ppl.run(working_dir=self.working_dir,
                              enable_conda=True,
                              conda=conda_,
                              poll_interval=0.1,
                              verbose=True)
        # Check the outputs
        self.assertEqual(exit_status,0)

    def test_pipeline_with_conda_custom_env_dir(self):
        """
        Pipeline: use custom env directory with conda dependency resolution
        """
        # Set up mock conda
        conda_dir = MockConda.create(os.path.join(self.working_dir,"conda"))
        conda_ = os.path.join(conda_dir,"bin","conda")
        # Custom env dir
        custom_env_dir = os.path.join(self.working_dir,"__conda_envs")
        # Check custom env dir doesn't exist
        self.assertFalse(os.path.exists(custom_env_dir))
        # Define a task
        class RunFastqc(PipelineTask):
            def init(self,*files):
                self.conda("fastqc=0.11.3")
                self.add_output('files',list())
            def setup(self):
                for f in self.args.files:
                    self.add_cmd(
                        PipelineCommandWrapper(
                            "Run fastqc for %s" % f,
                            "fastqc",f))
            def finish(self):
                for f in self.args.files:
                    self.output.files.append(f)
        # Build the pipeline
        ppl = Pipeline()
        task = RunFastqc("Run Fastqc",
                         "sample1.fastq","sample2.fastq")
        ppl.add_task(task)
        # Run the pipeline
        exit_status = ppl.run(working_dir=self.working_dir,
                              enable_conda=True,
                              conda=conda_,
                              conda_env_dir=custom_env_dir,
                              poll_interval=0.1,
                              verbose=True)
        # Check the outputs
        self.assertEqual(exit_status,0)
        # Check custom env dir exists
        self.assertTrue(os.path.isdir(custom_env_dir))

    def test_pipeline_with_conda_fail_to_create_environment(self):
        """
        Pipeline: handle failure with conda dependency resolution
        """
        # Set up mock conda with failing create command
        conda_dir = MockConda.create(os.path.join(self.working_dir,"conda"),
                                     create_fails=True)
        conda_ = os.path.join(conda_dir,"bin","conda")
        # Define a task
        class RunFastqc(PipelineTask):
            def init(self,*files):
                self.conda("fastqc=0.11.3")
                self.add_output('files',list())
            def setup(self):
                for f in self.args.files:
                    self.add_cmd(
                        PipelineCommandWrapper(
                            "Run fastqc for %s" % f,
                            "fastqc",f))
            def finish(self):
                for f in self.args.files:
                    self.output.files.append(f)
        # Build the pipeline
        ppl = Pipeline()
        task = RunFastqc("Run Fastqc",
                         "sample1.fastq","sample2.fastq")
        ppl.add_task(task)
        # Run the pipeline
        exit_status = ppl.run(working_dir=self.working_dir,
                              enable_conda=True,
                              conda=conda_,
                              poll_interval=0.1,
                              verbose=True)
        # Check the outputs
        self.assertEqual(exit_status,1)

    def test_pipeline_with_conda_fail_to_activate_new_environment(self):
        """
        Pipeline: handle failure to activate new conda environment
        """
        # Set up mock conda with failing create command
        conda_dir = MockConda.create(os.path.join(self.working_dir,"conda"),
                                     activate_fails=True)
        conda_ = os.path.join(conda_dir,"bin","conda")
        # Define a task
        class RunFastqc(PipelineTask):
            def init(self,*files):
                self.conda("fastqc=0.11.3")
                self.add_output('files',list())
            def setup(self):
                for f in self.args.files:
                    self.add_cmd(
                        PipelineCommandWrapper(
                            "Run fastqc for %s" % f,
                            "fastqc",f))
            def finish(self):
                for f in self.args.files:
                    self.output.files.append(f)
        # Build the pipeline
        ppl = Pipeline()
        task = RunFastqc("Run Fastqc",
                         "sample1.fastq","sample2.fastq")
        ppl.add_task(task)
        # Run the pipeline
        exit_status = ppl.run(working_dir=self.working_dir,
                              enable_conda=True,
                              conda=conda_,
                              poll_interval=0.1,
                              verbose=True)
        # Check the outputs
        self.assertEqual(exit_status,1)

    def test_pipeline_with_conda_fail_to_activate_existing_environment(self):
        """
        Pipeline: handle failure to activate existing conda environment
        """
        # Set up mock conda with failing create command
        conda_dir = MockConda.create(os.path.join(self.working_dir,"conda"),
                                     activate_fails=True)
        conda_ = os.path.join(conda_dir,"bin","conda")
        # Define a task
        class RunFastqc(PipelineTask):
            def init(self,*files):
                self.conda("fastqc=0.11.3")
                self.add_output('files',list())
            def setup(self):
                for f in self.args.files:
                    self.add_cmd(
                        PipelineCommandWrapper(
                            "Run fastqc for %s" % f,
                            "fastqc",f))
            def finish(self):
                for f in self.args.files:
                    self.output.files.append(f)
        # Mock an existing conda environment for the task
        os.makedirs(os.path.join(self.working_dir,
                                 "__conda",
                                 "envs",
                                 "fastqc@0.11.3"))
        # Build the pipeline
        ppl = Pipeline()
        task = RunFastqc("Run Fastqc",
                         "sample1.fastq","sample2.fastq")
        ppl.add_task(task)
        # Run the pipeline
        exit_status = ppl.run(working_dir=self.working_dir,
                              enable_conda=True,
                              conda=conda_,
                              poll_interval=0.1,
                              verbose=True)
        # Check the outputs
        self.assertEqual(exit_status,1)

    def test_pipeline_define_outputs(self):
        """
        Pipeline: test defining pipeline outputs
        """
        # Define a reusable task
        # Appends item to a string
        class AppendString(PipelineTask):
            def init(self,s1,s2):
                self.add_output('string',PipelineParam(type=str))
            def setup(self):
                self.output.string.set(str(self.args.s1) + str(self.args.s2))
        # Build the pipeline
        ppl = Pipeline()
        task1 = AppendString("Append 1","This ","is ")
        task2 = AppendString("Append 2",task1.output.string,"the full string")
        ppl.add_task(task2,requires=(task1,))
        # Define outputs
        ppl.add_output('result',task2.output.string)
        # Run the pipeline
        exit_status = ppl.run(working_dir=self.working_dir,
                              poll_interval=0.1)
        # Check the outputs
        self.assertEqual(exit_status,0)
        self.assertTrue(isinstance(task1.output.string,PipelineParam))
        self.assertEqual(task1.output.string.value,"This is ")
        self.assertTrue(isinstance(task2.output.string,PipelineParam))
        self.assertEqual(task2.output.string.value,"This is the full string")
        self.assertEqual(ppl.output.result,"This is the full string")

    def test_pipeline_dont_finalize_outputs(self):
        """
        Pipeline: test not finalizing pipeline outputs
        """
        # Define a reusable task
        # Appends item to a string
        class AppendString(PipelineTask):
            def init(self,s1,s2):
                self.add_output('string',PipelineParam(type=str))
            def setup(self):
                self.output.string.set(str(self.args.s1) + str(self.args.s2))
        # Build the pipeline
        ppl = Pipeline()
        task1 = AppendString("Append 1","This ","is ")
        task2 = AppendString("Append 2",task1.output.string,"the full string")
        ppl.add_task(task2,requires=(task1,))
        # Define outputs
        ppl.add_output('result',task2.output.string)
        # Run the pipeline with output finalization turned off
        exit_status = ppl.run(finalize_outputs=False,
                              working_dir=self.working_dir,
                              poll_interval=0.1)
        # Check the outputs
        self.assertEqual(exit_status,0)
        self.assertTrue(isinstance(task1.output.string,PipelineParam))
        self.assertEqual(task1.output.string.value,"This is ")
        self.assertTrue(isinstance(task2.output.string,PipelineParam))
        self.assertEqual(task2.output.string.value,"This is the full string")
        self.assertTrue(isinstance(ppl.output.result,PipelineParam))
        self.assertEqual(ppl.output.result.value,"This is the full string")

    def test_pipeline_method_task_list(self):
        """
        Pipeline: test the 'task_list' method
        """
        # Define a reusable task
        # Appends item to a list
        class Append(PipelineTask):
            def init(self,l,s):
                self.add_output('list',list())
            def setup(self):
                for item in self.args.l:
                    self.output.list.append(item)
                self.output.list.append(self.args.s)
        # Make an empty pipeline
        ppl = Pipeline()
        self.assertEqual(ppl.task_list(),[])
        # Add a task
        task1 = Append("Append 1",(),"item1")
        task2 = Append("Append 2",task1.output.list,"item2")
        ppl.add_task(task2,requires=(task1,))
        # Check the task list
        task_list = ppl.task_list()
        self.assertEqual(len(task_list),2)
        self.assertTrue(task1.id() in task_list)
        self.assertTrue(task2.id() in task_list)

    def test_pipeline_method_get_task(self):
        """
        Pipeline: test the 'get_task' method
        """
        # Define a reusable task
        # Appends item to a list
        class Append(PipelineTask):
            def init(self,l,s):
                self.add_output('list',list())
            def setup(self):
                for item in self.args.l:
                    self.output.list.append(item)
                self.output.list.append(self.args.s)
        # Make a pipeline
        ppl = Pipeline()
        task1 = Append("Append 1",(),"item1")
        task2 = Append("Append 2",task1.output.list,"item2")
        ppl.add_task(task2,requires=(task1,))
        # Fetch task data
        task1_data = ppl.get_task(task1.id())
        self.assertEqual(task1_data[0],task1)
        self.assertEqual(task1_data[1],())
        self.assertEqual(task1_data[2],{})
        task2_data = ppl.get_task(task2.id())
        self.assertEqual(task2_data[0],task2)
        self.assertEqual(task2_data[1],(task1,))
        self.assertEqual(task2_data[2],{})

    def test_pipeline_method_rank_tasks(self):
        """
        Pipeline: test the 'rank_tasks' method
        """
        # Define a reusable task
        # Appends item to a list
        class Append(PipelineTask):
            def init(self,l,s):
                self.add_output('list',list())
            def setup(self):
                for item in self.args.l:
                    self.output.list.append(item)
                self.output.list.append(self.args.s)
        # Make a pipeline
        ppl = Pipeline()
        task1 = Append("Append 1",(),"item1")
        task2 = Append("Append 2",task1.output.list,"item2")
        task3 = Append("Append 3",task1.output.list,"item3")
        task4 = Append("Append 4",task3.output.list,"item4")
        ppl.add_task(task2,requires=(task1,))
        ppl.add_task(task3,requires=(task1,))
        ppl.add_task(task4,requires=(task3,))
        # Rank the tasks
        ranked_tasks = ppl.rank_tasks()
        # Should be 3 ranks
        self.assertEqual(len(ranked_tasks),3)
        # Check the ranks
        self.assertEqual(ranked_tasks[0],[task1.id()])
        self.assertEqual(sorted(ranked_tasks[1]),
                         sorted([task2.id(),task3.id()]))
        self.assertEqual(ranked_tasks[2],[task4.id()])

    def test_pipeline_method_initial_tasks(self):
        """
        Pipeline: test the 'initial_tasks' method
        """
        # Define a reusable task
        # Appends item to a list
        class Append(PipelineTask):
            def init(self,l,s):
                self.add_output('list',list())
            def setup(self):
                for item in self.args.l:
                    self.output.list.append(item)
                self.output.list.append(self.args.s)
        # Make a pipeline
        ppl = Pipeline()
        task1 = Append("Append 1",(),"item1")
        task2 = Append("Append 2",task1.output.list,"item2")
        task3 = Append("Append 3",task1.output.list,"item3")
        task4 = Append("Append 4",task3.output.list,"item4")
        ppl.add_task(task2,requires=(task1,))
        ppl.add_task(task3,requires=(task1,))
        ppl.add_task(task4,requires=(task3,))
        # Check the initial tasks
        self.assertEqual(ppl.initial_tasks,[task1])

    def test_pipeline_method_final_tasks(self):
        """
        Pipeline: test the 'final_tasks' method
        """
        # Define a reusable task
        # Appends item to a list
        class Append(PipelineTask):
            def init(self,l,s):
                self.add_output('list',list())
            def setup(self):
                for item in self.args.l:
                    self.output.list.append(item)
                self.output.list.append(self.args.s)
        # Make a pipeline
        ppl = Pipeline()
        task1 = Append("Append 1",(),"item1")
        task2 = Append("Append 2",task1.output.list,"item2")
        task3 = Append("Append 3",task1.output.list,"item3")
        task4 = Append("Append 4",task3.output.list,"item4")
        ppl.add_task(task2,requires=(task1,))
        ppl.add_task(task3,requires=(task1,))
        ppl.add_task(task4,requires=(task3,))
        # Check the initial tasks
        self.assertEqual(ppl.final_tasks,sorted([task2,task4],
                                                key=lambda x: x.id()))

    def test_pipeline_method_get_dependent_tasks(self):
        """
        Pipeline: test the 'get_dependent_tasks' method
        """
        # Define a reusable task
        # Appends item to a list
        class Append(PipelineTask):
            def init(self,l,s):
                self.add_output('list',list())
            def setup(self):
                for item in self.args.l:
                    self.output.list.append(item)
                self.output.list.append(self.args.s)
        # Make a pipeline
        ppl = Pipeline()
        task1 = Append("Append 1",(),"item1")
        task2 = Append("Append 2",task1.output.list,"item2")
        task3 = Append("Append 3",task1.output.list,"item3")
        task4 = Append("Append 4",task3.output.list,"item4")
        ppl.add_task(task2,requires=(task1,))
        ppl.add_task(task3,requires=(task1,))
        ppl.add_task(task4,requires=(task3,))
        # Check the dependent tasks
        self.assertEqual(sorted(ppl.get_dependent_tasks(task1.id())),
                         sorted([task2.id(),task3.id(),task4.id()]))
        self.assertEqual(ppl.get_dependent_tasks(task2.id()),[])
        self.assertEqual(ppl.get_dependent_tasks(task3.id()),[task4.id()])
        self.assertEqual(ppl.get_dependent_tasks(task4.id()),[])

    def test_pipeline_append_pipeline(self):
        """
        Pipeline: append one pipeline to another
        """
        # Define a reusable task
        # Appends item to a list
        class Append(PipelineTask):
            def init(self,l,s):
                self.add_output('list',list())
            def setup(self):
                for item in self.args.l:
                    self.output.list.append(item)
                self.output.list.append(self.args.s)
        # Make first pipeline
        ppl1 = Pipeline()
        ppl1.add_param("param1")
        ppl1.add_runner("runner1")
        task1 = Append("Append 1",(),"item1")
        task2 = Append("Append 2",task1.output.list,"item2")
        task3 = Append("Append 3",task1.output.list,"item3")
        task4 = Append("Append 4",task3.output.list,"item4")
        ppl1.add_task(task2,requires=(task1,))
        ppl1.add_task(task3,requires=(task1,))
        ppl1.add_task(task4,requires=(task3,))
        self.assertEqual(len(ppl1.task_list()),4)
        # Make second pipeline
        ppl2 = Pipeline()
        ppl2.add_param("param2")
        ppl2.add_runner("runner2")
        task5 = Append("Append 5",task1.output.list,"item5")
        task6 = Append("Append 6",task3.output.list,"item6")
        task7 = Append("Append 7",task3.output.list,"item7")
        ppl2.add_task(task6,requires=(task5,))
        ppl2.add_task(task7,requires=(task6,))
        self.assertEqual(len(ppl2.task_list()),3)
        # Append second pipeline to the first
        ppl1.append_pipeline(ppl2)
        self.assertEqual(len(ppl1.task_list()),7)
        # Check requirements on first task of pipeline 2
        # have been updated
        self.assertEqual(
            sorted(ppl1.get_task(task5.id())[1],key=lambda t: t.id()),
            sorted([task2,task4,],key=lambda t: t.id())
        )
        # Check params from both pipelines are defined
        self.assertTrue('param1' in ppl1.params)
        self.assertTrue('param2' in ppl1.params)
        # Check runners from both pipelines are defined
        self.assertTrue('runner1' in ppl1.runners)
        self.assertTrue('runner2' in ppl1.runners)

    def test_pipeline_merge_pipeline(self):
        """
        Pipeline: merge one pipeline into another
        """
        # Define a reusable task
        # Appends item to a list
        class Append(PipelineTask):
            def init(self,l,s):
                self.add_output('list',list())
            def setup(self):
                for item in self.args.l:
                    self.output.list.append(item)
                self.output.list.append(self.args.s)
        # Make first pipeline
        ppl1 = Pipeline()
        ppl1.add_param("param1")
        ppl1.add_runner("runner1")
        task1 = Append("Append 1",(),"item1")
        task2 = Append("Append 2",task1.output.list,"item2")
        task3 = Append("Append 3",task1.output.list,"item3")
        task4 = Append("Append 4",task3.output.list,"item4")
        ppl1.add_task(task2,requires=(task1,))
        ppl1.add_task(task3,requires=(task1,))
        ppl1.add_task(task4,requires=(task3,))
        self.assertEqual(len(ppl1.task_list()),4)
        # Make second pipeline
        ppl2 = Pipeline()
        ppl2.add_param("param2")
        ppl2.add_runner("runner2")
        task5 = Append("Append 5",task1.output.list,"item5")
        task6 = Append("Append 6",task3.output.list,"item6")
        task7 = Append("Append 7",task3.output.list,"item7")
        ppl2.add_task(task6,requires=(task5,))
        ppl2.add_task(task7,requires=(task6,))
        self.assertEqual(len(ppl2.task_list()),3)
        # Merge second pipeline into the first
        ppl1.merge_pipeline(ppl2)
        self.assertEqual(len(ppl1.task_list()),7)
        # Check params from both pipelines are defined
        self.assertTrue('param1' in ppl1.params)
        self.assertTrue('param2' in ppl1.params)
        # Check runners from both pipelines are defined
        self.assertTrue('runner1' in ppl1.runners)
        self.assertTrue('runner2' in ppl1.runners)

    def test_pipeline_add_pipeline(self):
        """
        Pipeline: add one pipeline into another
        """
        # Define a reusable task
        # Appends item to a list
        class Append(PipelineTask):
            def init(self,l,s):
                self.add_output('list',list())
            def setup(self):
                for item in self.args.l:
                    self.output.list.append(item)
                self.output.list.append(self.args.s)
        # Make first pipeline
        ppl1 = Pipeline()
        ppl1.add_param("param1")
        ppl1.add_runner("runner1")
        task1 = Append("Append 1",(),"item1")
        task2 = Append("Append 2",task1.output.list,"item2")
        task3 = Append("Append 3",task1.output.list,"item3")
        task4 = Append("Append 4",task3.output.list,"item4")
        ppl1.add_task(task2,requires=(task1,))
        ppl1.add_task(task3,requires=(task1,))
        ppl1.add_task(task4,requires=(task3,))
        self.assertEqual(len(ppl1.task_list()),4)
        # Make second pipeline
        ppl2 = Pipeline()
        ppl2.add_param("param2")
        ppl2.add_runner("runner2")
        task5 = Append("Append 5",task1.output.list,"item5")
        task6 = Append("Append 6",task3.output.list,"item6")
        task7 = Append("Append 7",task3.output.list,"item7")
        ppl2.add_task(task6,requires=(task5,))
        ppl2.add_task(task7,requires=(task6,))
        self.assertEqual(len(ppl2.task_list()),3)
        # Merge second pipeline into the first
        ppl1.add_pipeline(ppl2)
        self.assertEqual(len(ppl1.task_list()),7)
        # Check params from both pipelines are defined
        self.assertTrue('param1' in ppl1.params)
        self.assertTrue('param2' in ppl1.params)
        # Check runners from both pipelines are defined
        self.assertTrue('runner1' in ppl1.runners)
        self.assertTrue('runner2' in ppl1.runners)

class TestPipelineTask(unittest.TestCase):

    def setUp(self):
        # Set up a scheduler
        self.sched = SimpleScheduler(poll_interval=0.01)
        self.sched.start()
        # Make a temporary working dir
        self.working_dir = tempfile.mkdtemp(
            suffix='TestPipeline')
        # Store PATH
        self.path = os.environ['PATH']

    def tearDown(self):
        # Stop the scheduler
        if self.sched is not None:
            self.sched.stop()
        # Remove temp dir
        if os.path.exists(self.working_dir):
            shutil.rmtree(self.working_dir)
        # Restore PATH
        os.environ['PATH'] = self.path

    def _user(self):
        # Internal function to determine user
        return getpass.getuser()

    def _hostname(self):
        # Internal function to determine hostname
        try:
            return os.environ['HOSTNAME']
        except KeyError:
            # HOSTNAME not defined in the
            # environment, try 'platform'
            # module instead
            return platform.node()

    def test_pipelinetask_invocations(self):
        """
        PipelineTask: check task methods are invoked
        """
        # Define a simplistic task which does nothing
        class CheckInvocations(PipelineTask):
            def init(self):
                self.add_output('invocations',list())
                self.output.invocations.append("init")
            def setup(self):
                self.output.invocations.append("setup")
            def finish(self):
                self.output.invocations.append("finish")
        # Make a task instance
        task = CheckInvocations("Check method invocations")
        # Check initial state
        self.assertFalse(task.completed)
        self.assertEqual(task.exit_code,None)
        self.assertEqual(task.output.invocations,
                         ["init"])
        # Run the task
        task.run(sched=self.sched,
                 working_dir=self.working_dir,
                 asynchronous=False)
        # Check final state
        self.assertTrue(task.completed)
        self.assertEqual(task.exit_code,0)
        self.assertEqual(task.output.invocations,
                         ["init","setup","finish"])

    def test_pipelinetask_init(self):
        """
        PipelineTask: check task 'init' invocations
        """
        # Define a task for testing
        class CheckInit(PipelineTask):
            def init(self,a,b,c='hello',d=13,e=None):
                self.add_output('results',list())
            def setup(self):
                result = "a=%s b=%s c=%s d=%s e=%s" \
                         % (self.args.a,
                            self.args.b,
                            self.args.c,
                            self.args.d,
                            self.args.e)
                self.output.results.append(result)
        # Make a task instance with minimal arglist
        task = CheckInit("Minimal arglist","a","b")
        self.assertEqual(task.args.a,"a")
        self.assertEqual(task.args.b,"b")
        self.assertEqual(task.args.c,"hello")
        self.assertEqual(task.args.d,13)
        self.assertEqual(task.args.e,None)
        # Make a task instance with named minimal arglist
        task = CheckInit("Named minimal arglist",a="a",b="b")
        self.assertEqual(task.args.a,"a")
        self.assertEqual(task.args.b,"b")
        self.assertEqual(task.args.c,"hello")
        self.assertEqual(task.args.d,13)
        self.assertEqual(task.args.e,None)
        # Make a task instance with named minimal arglist (reversed)
        task = CheckInit("Named minimal arglist reversed",
                         b="a",a="b")
        self.assertEqual(task.args.a,"b")
        self.assertEqual(task.args.b,"a")
        self.assertEqual(task.args.c,"hello")
        self.assertEqual(task.args.d,13)
        self.assertEqual(task.args.e,None)
        # Make a task instance with args and subset of keywords
        task = CheckInit("Args and subset of keywords",
                         "a","b",e=True,d=12)
        self.assertEqual(task.args.a,"a")
        self.assertEqual(task.args.b,"b")
        self.assertEqual(task.args.c,"hello")
        self.assertEqual(task.args.d,12)
        self.assertEqual(task.args.e,True)
        # Make a task instance with full arglist with keywords
        task = CheckInit("Full arglist with keywords",
                         "a","b",c="goodbye",d=12,e=True)
        self.assertEqual(task.args.a,"a")
        self.assertEqual(task.args.b,"b")
        self.assertEqual(task.args.c,"goodbye")
        self.assertEqual(task.args.d,12)
        self.assertEqual(task.args.e,True)
        # Make a task instance with full arglist no keywords
        task = CheckInit("Full arglist no keywords",
                         "a","b","goodbye",12,True)
        self.assertEqual(task.args.a,"a")
        self.assertEqual(task.args.b,"b")
        self.assertEqual(task.args.c,"goodbye")
        self.assertEqual(task.args.d,12)
        self.assertEqual(task.args.e,True)
        # Check task raises exception if init fails
        # Define a task for testing
        class FailInit(PipelineTask):
            def init(self,a,b,c='hello',d=13,e=None):
                raise Exception("Forced init to fail")
            def setup(self):
                result = "a=%s b=%s c=%s d=%s e=%s" \
                         % (self.args.a,
                            self.args.b,
                            self.args.c,
                            self.args.d,
                            self.args.e)
                self.output.results.append(result)
        self.assertRaises(PipelineError,
                          FailInit,
                          "This will fail on init",
                          "a",
                          "b")

    def test_pipelinetask_requirements(self):
        """
        PipelineTask: check task requirements
        """
        # Define task for testing
        class AppendTask(PipelineTask):
            def init(self,*inputs):
                self.add_output('result',list())
            def setup(self):
                for x in self.args.inputs:
                    self.output.results.append(x)
        # Instantiate tasks
        t1 = AppendTask("Task1",1,2)
        t2 = AppendTask("Task2",3,4)
        t3 = AppendTask("Task3",5,6)
        # Check requirements on all tasks
        self.assertEqual(t1.required_task_ids,[])
        self.assertEqual(t2.required_task_ids,[])
        self.assertEqual(t3.required_task_ids,[])
        # Make second task depend on first
        t2.requires(t1)
        self.assertEqual(t1.required_task_ids,[])
        self.assertEqual(t2.required_task_ids,[t1.id()])
        self.assertEqual(t3.required_task_ids,[])
        # Make third task depend on first and second
        t3.requires(t1,t2)
        self.assertEqual(t1.required_task_ids,[])
        self.assertEqual(t2.required_task_ids,[t1.id()])
        self.assertEqual(t3.required_task_ids,
                         sorted([t2.id(),t1.id()]))

    def test_pipelinetask_requirements_as_ids(self):
        """
        PipelineTask: check task requirements supplied as IDs
        """
        # Define task for testing
        class AppendTask(PipelineTask):
            def init(self,*inputs):
                self.add_output('result',list())
            def setup(self):
                for x in self.args.inputs:
                    self.output.results.append(x)
        # Instantiate tasks
        t1 = AppendTask("Task1",1,2)
        t2 = AppendTask("Task2",3,4)
        t3 = AppendTask("Task3",5,6)
        # Check requirements on all tasks
        self.assertEqual(t1.required_task_ids,[])
        self.assertEqual(t2.required_task_ids,[])
        self.assertEqual(t3.required_task_ids,[])
        # Make second task depend on first
        t2.requires_id(t1.id())
        self.assertEqual(t1.required_task_ids,[])
        self.assertEqual(t2.required_task_ids,[t1.id()])
        self.assertEqual(t3.required_task_ids,[])
        # Make third task depend on first and second
        t3.requires_id(t1.id())
        t3.requires_id(t2.id())
        self.assertEqual(t1.required_task_ids,[])
        self.assertEqual(t2.required_task_ids,[t1.id()])
        self.assertEqual(t3.required_task_ids,
                         sorted([t2.id(),t1.id()]))

    def test_pipelinetask_required_by(self):
        """
        PipelineTask: check tasks required by others
        """
        # Define task for testing
        class AppendTask(PipelineTask):
            def init(self,*inputs):
                self.add_output('result',list())
            def setup(self):
                for x in self.args.inputs:
                    self.output.results.append(x)
        # Instantiate tasks
        t1 = AppendTask("Task1",1,2)
        t2 = AppendTask("Task2",3,4)
        t3 = AppendTask("Task3",5,6)
        # Check requirements on all tasks
        self.assertEqual(t1.required_task_ids,[])
        self.assertEqual(t2.required_task_ids,[])
        self.assertEqual(t3.required_task_ids,[])
        # Make second and third task depend on first
        t1.required_by(t2,t3)
        self.assertEqual(t1.required_task_ids,[])
        self.assertEqual(t2.required_task_ids,[t1.id()])
        self.assertEqual(t3.required_task_ids,[t1.id()])
        # Make third task depend on second
        t2.required_by(t3)
        self.assertEqual(t1.required_task_ids,[])
        self.assertEqual(t2.required_task_ids,[t1.id()])
        self.assertEqual(t3.required_task_ids,
                         sorted([t2.id(),t1.id()]))

    def test_pipelinetask_implied_requirement_from_input_param(self):
        """
        PipelineTask: check implied task requirements from inputs
        """
        # Define task for testing
        class AppendTask(PipelineTask):
            def init(self,*inputs,**kws):
                self.add_output('result',PipelineParam(type=list()))
            def setup(self):
                for x in self.args.inputs:
                    self.output.results.value.append(x)
        # Instantiate tasks
        t1 = AppendTask("Task1",1,2)
        t2 = AppendTask("Task2",t1.output.result,4)
        t3 = AppendTask("Task3",t1.output.result,extras=t2.output.result)
        # Check requirements on both tasks
        self.assertEqual(t1.required_task_ids,[])
        self.assertEqual(t2.required_task_ids,[t1.id()])
        self.assertEqual(t3.required_task_ids,
                         sorted([t2.id(),t1.id()]))

    def test_pipelinetask_raise_exception_for_non_task_requirement(self):
        """
        PipelineTask: raise exception if requirement is not a task
        """
        # Define stask for testing
        class AppendTask(PipelineTask):
            def init(self,*inputs):
                self.add_output('result',list())
            def setup(self):
                for x in self.args.inputs:
                    self.output.results.append(x)
        # Instantiate task
        t1 = AppendTask(1,2)
        # Check initial requirements
        self.assertEqual(t1.required_task_ids,[])
        # Raise exception by trying to adding a non-task
        # object as a requirement
        self.assertRaises(PipelineError,
                          t1.requires,
                          "not_a_task")

    def test_pipelinetask_no_commands(self):
        """
        PipelineTask: run task with no commands
        """
        # Define a task with no commands
        class Add(PipelineTask):
            def init(self,x,y):
                self.add_output('result',list())
            def setup(self):
                self.output.result.append(self.args.x+self.args.y)
        # Make a task instance
        task = Add("Add two numbers",1,2)
        # Check initial state
        self.assertEqual(task.args.x,1)
        self.assertEqual(task.args.y,2)
        self.assertFalse(task.completed)
        self.assertEqual(task.exit_code,None)
        self.assertEqual(task.output.result,[])
        # Run the task
        task.run(sched=self.sched,
                 working_dir=self.working_dir,
                 asynchronous=False)
        # Check final state
        self.assertTrue(task.completed)
        self.assertEqual(task.exit_code,0)
        self.assertEqual(task.output.result,[3])
        self.assertEqual(task.stdout,"")

    def test_pipelinetask_with_commands(self):
        """
        PipelineTask: run task with shell command
        """
        # Define a task with a command
        # Echoes text via shell command
        class Echo(PipelineTask):
            def init(self,s):
                pass
            def setup(self):
                self.add_cmd(
                    PipelineCommandWrapper(
                        "Echo text","echo",self.args.s))
        # Make a task instance
        task = Echo("Echo string","Hello!")
        # Check initial state
        self.assertEqual(task.args.s,"Hello!")
        self.assertFalse(task.completed)
        self.assertEqual(task.exit_code,None)
        self.assertFalse(task.output)
        # Run the task
        task.run(sched=self.sched,
                 working_dir=self.working_dir,
                 asynchronous=False)
        # Check final state
        self.assertTrue(task.completed)
        self.assertEqual(task.exit_code,0)
        self.assertFalse(task.output)
        # Check stdout
        # Should look like:
        # #### COMMAND Echo text
        # #### HOSTNAME popov
        # #### USER pjb
        # #### START Thu Aug 17 08:38:14 BST 2017
        # #### CWD /tmp/dir
        # Hello!
        # #### END Thu Aug 17 08:38:14 BST 2017
        # #### EXIT_CODE 0
        stdout = task.stdout.split("\n")
        self.assertEqual(len(stdout),9) # 9 = 8 + trailing newline
        self.assertEqual(stdout[0],"#### COMMAND Echo text")
        self.assertEqual(stdout[1],"#### HOSTNAME %s" % self._hostname())
        self.assertEqual(stdout[2],"#### USER %s" % self._user())
        self.assertTrue(stdout[3].startswith("#### START "))
        self.assertEqual(stdout[4],"#### CWD %s" % self.working_dir)
        self.assertEqual(stdout[5],"Hello!")
        self.assertTrue(stdout[6].startswith("#### END "))
        self.assertEqual(stdout[7],"#### EXIT_CODE 0")

    def test_pipelinetask_with_multiple_commands(self):
        """
        PipelineTask: run task with multiple shell commands
        """
        # Define a task with a command
        # Echoes text via shell command
        class EchoMany(PipelineTask):
            def init(self,*s):
                pass
            def setup(self):
                for s in self.args.s:
                    self.add_cmd(
                        PipelineCommandWrapper(
                            "Echo text","echo",s))
        # Make a task instance
        task = EchoMany("Echo string","Hello!","Goodbye!")
        # Check initial state
        self.assertEqual(task.args.s,("Hello!","Goodbye!"))
        self.assertFalse(task.completed)
        self.assertEqual(task.exit_code,None)
        self.assertFalse(task.output)
        # Run the task
        task.run(sched=self.sched,
                 working_dir=self.working_dir,
                 asynchronous=False)
        # Check final state
        self.assertTrue(task.completed)
        self.assertEqual(task.exit_code,0)
        self.assertFalse(task.output)
        # Check stdout
        # Should look like:
        # #### COMMAND Echo text
        # #### HOSTNAME popov
        # #### USER pjb
        # #### START Thu Aug 17 08:38:14 BST 2017
        # #### CWD /tmp/dir
        # Hello!
        # #### END Thu Aug 17 08:38:14 BST 2017
        # #### EXIT_CODE 0
        # #### COMMAND Echo text
        # #### HOSTNAME popov
        # #### USER pjb
        # #### START Thu Aug 17 08:38:14 BST 2017
        # #### CWD /tmp/dir
        # Goodbye!
        # #### END Thu Aug 17 08:38:14 BST 2017
        # #### EXIT_CODE 0
        stdout = task.stdout.split("\n")
        self.assertEqual(len(stdout),17) # 17 = 16 + trailing newline
        self.assertEqual(stdout[0],"#### COMMAND Echo text")
        self.assertEqual(stdout[1],"#### HOSTNAME %s" % self._hostname())
        self.assertEqual(stdout[2],"#### USER %s" % self._user())
        self.assertTrue(stdout[3].startswith("#### START "))
        self.assertEqual(stdout[4],"#### CWD %s" % self.working_dir)
        self.assertEqual(stdout[5],"Hello!")
        self.assertTrue(stdout[6].startswith("#### END "))
        self.assertEqual(stdout[7],"#### EXIT_CODE 0")
        self.assertEqual(stdout[8],"#### COMMAND Echo text")
        self.assertEqual(stdout[9],"#### HOSTNAME %s" % self._hostname())
        self.assertEqual(stdout[10],"#### USER %s" % self._user())
        self.assertTrue(stdout[11].startswith("#### START "))
        self.assertEqual(stdout[12],"#### CWD %s" % self.working_dir)
        self.assertEqual(stdout[13],"Goodbye!")
        self.assertTrue(stdout[14].startswith("#### END "))
        self.assertEqual(stdout[15],"#### EXIT_CODE 0")

    def test_pipelinetask_with_batched_commands(self):
        """
        PipelineTask: run task with batched shell commands
        """
        # Define a task with a command
        # Echoes text via shell command
        class EchoMany(PipelineTask):
            def init(self,*s):
                pass
            def setup(self):
                for s in self.args.s:
                    self.add_cmd(
                        PipelineCommandWrapper(
                            "Echo text","echo",s))
        # Make a task instance
        task = EchoMany("Echo string",
                        "Hello!",
                        "Bonjour!",
                        "Takk!",
                        "Wilkommen!",
                        "Benvenuto!")
        # Check initial state
        self.assertEqual(task.args.s,
                         ("Hello!",
                          "Bonjour!",
                          "Takk!",
                          "Wilkommen!",
                          "Benvenuto!"))
        self.assertFalse(task.completed)
        self.assertEqual(task.exit_code,None)
        self.assertFalse(task.output)
        # Run the task with batches
        task.run(sched=self.sched,
                 working_dir=self.working_dir,
                 batch_size=2,
                 asynchronous=False)
        # Check final state
        self.assertTrue(task.completed)
        self.assertEqual(task.exit_code,0)
        self.assertFalse(task.output)
        # Check stdout
        # Should look like:
        # #### COMMAND Batch commands for Echo string
        # #### BATCH 1
        # #### HOSTNAME popov
        # #### USER pjb
        # #### START Thu Aug 17 08:38:14 BST 2017
        # #### CWD /tmp/dir
        # Hello!
        # Bonjour!
        # #### END Thu Aug 17 08:38:14 BST 2017
        # #### EXIT_CODE 0
        # #### COMMAND Batch commands for Echo string
        # #### BATCH 2
        # #### HOSTNAME popov
        # #### USER pjb
        # #### START Thu Aug 17 08:38:14 BST 2017
        # #### CWD /tmp/dir
        # Takk!
        # Wilkommen!
        # #### END Thu Aug 17 08:38:14 BST 2017
        # #### EXIT_CODE 0
        # #### COMMAND Batch commands for Echo string
        # #### BATCH 3
        # #### HOSTNAME popov
        # #### USER pjb
        # #### START Thu Aug 17 08:38:14 BST 2017
        # #### CWD /tmp/dir
        # Benvenuto!
        # #### END Thu Aug 17 08:38:14 BST 2017
        # #### EXIT_CODE 0
        stdout = task.stdout.split("\n")
        self.assertEqual(len(stdout),30) # 30 = 29 + trailing newline
        self.assertEqual(stdout[0],"#### COMMAND Batch commands for Echo "
                         "string")
        self.assertEqual(stdout[1],"#### BATCH 1")
        self.assertEqual(stdout[2],"#### HOSTNAME %s" % self._hostname())
        self.assertEqual(stdout[3],"#### USER %s" % self._user())
        self.assertTrue(stdout[4].startswith("#### START "))
        self.assertEqual(stdout[5],"#### CWD %s" % self.working_dir)
        self.assertEqual(stdout[6],"Hello!")
        self.assertEqual(stdout[7],"Bonjour!")
        self.assertTrue(stdout[8].startswith("#### END "))
        self.assertEqual(stdout[9],"#### EXIT_CODE 0")
        self.assertEqual(stdout[10],"#### COMMAND Batch commands for Echo "
                         "string")
        self.assertEqual(stdout[11],"#### BATCH 2")
        self.assertEqual(stdout[12],"#### HOSTNAME %s" % self._hostname())
        self.assertEqual(stdout[13],"#### USER %s" % self._user())
        self.assertTrue(stdout[14].startswith("#### START "))
        self.assertEqual(stdout[15],"#### CWD %s" % self.working_dir)
        self.assertEqual(stdout[16],"Takk!")
        self.assertEqual(stdout[17],"Wilkommen!")
        self.assertTrue(stdout[18].startswith("#### END "))
        self.assertEqual(stdout[19],"#### EXIT_CODE 0")
        self.assertEqual(stdout[20],"#### COMMAND Batch commands for Echo "
                         "string")
        self.assertEqual(stdout[21],"#### BATCH 3")
        self.assertEqual(stdout[22],"#### HOSTNAME %s" % self._hostname())
        self.assertEqual(stdout[23],"#### USER %s" % self._user())
        self.assertTrue(stdout[24].startswith("#### START "))
        self.assertEqual(stdout[25],"#### CWD %s" % self.working_dir)
        self.assertEqual(stdout[26],"Benvenuto!")
        self.assertTrue(stdout[27].startswith("#### END "))
        self.assertEqual(stdout[28],"#### EXIT_CODE 0")

    def test_pipelinetask_with_batched_functions(self):
        """
        PipelineTask: run task with batched functions
        """
        # Define a task with a command
        # Echoes text via Python function
        class EchoMany(PipelineFunctionTask):
            def init(self,*s):
                pass
            def setup(self):
                for s in self.args.s:
                    self.add_call("Echo text",self.echo,s)
            def echo(self,s):
                print(s)
        # Make a task instance
        task = EchoMany("Echo string",
                        "Hello!",
                        "Bonjour!",
                        "Takk!",
                        "Wilkommen!",
                        "Benvenuto!")
        # Check initial state
        self.assertEqual(task.args.s,
                         ("Hello!",
                          "Bonjour!",
                          "Takk!",
                          "Wilkommen!",
                          "Benvenuto!"))
        self.assertFalse(task.completed)
        self.assertEqual(task.exit_code,None)
        self.assertFalse(task.output)
        # Run the task with batches
        task.run(sched=self.sched,
                 working_dir=self.working_dir,
                 batch_size=2,
                 asynchronous=False)
        # Check final state
        self.assertTrue(task.completed)
        self.assertEqual(task.exit_code,0)
        self.assertFalse(task.output)
        # Check stdout
        # Should look like:
        # #### COMMAND Batch commands for Echo string
        # #### BATCH 1
        # #### HOSTNAME popov
        # #### USER pjb
        # #### START Thu Aug 17 08:38:14 BST 2017
        # #### CWD /tmp/dir
        # Hello!
        # Bonjour!
        # #### END Thu Aug 17 08:38:14 BST 2017
        # #### EXIT_CODE 0
        # #### COMMAND Batch commands for Echo string
        # #### BATCH 2
        # #### HOSTNAME popov
        # #### USER pjb
        # #### START Thu Aug 17 08:38:14 BST 2017
        # #### CWD /tmp/dir
        # Takk!
        # Wilkommen!
        # #### END Thu Aug 17 08:38:14 BST 2017
        # #### EXIT_CODE 0
        # #### COMMAND Batch commands for Echo string
        # #### BATCH 3
        # #### HOSTNAME popov
        # #### USER pjb
        # #### START Thu Aug 17 08:38:14 BST 2017
        # #### CWD /tmp/dir
        # Benvenuto!
        # #### END Thu Aug 17 08:38:14 BST 2017
        # #### EXIT_CODE 0
        stdout = task.stdout.split("\n")
        self.assertEqual(len(stdout),30) # 30 = 29 + trailing newline
        self.assertEqual(stdout[0],"#### COMMAND Batch commands for Echo "
                         "string")
        self.assertEqual(stdout[1],"#### BATCH 1")
        self.assertEqual(stdout[2],"#### HOSTNAME %s" % self._hostname())
        self.assertEqual(stdout[3],"#### USER %s" % self._user())
        self.assertTrue(stdout[4].startswith("#### START "))
        self.assertEqual(stdout[5],"#### CWD %s" % self.working_dir)
        self.assertEqual(stdout[6],"Hello!")
        self.assertEqual(stdout[7],"Bonjour!")
        self.assertTrue(stdout[8].startswith("#### END "))
        self.assertEqual(stdout[9],"#### EXIT_CODE 0")
        self.assertEqual(stdout[10],"#### COMMAND Batch commands for Echo "
                         "string")
        self.assertEqual(stdout[11],"#### BATCH 2")
        self.assertEqual(stdout[12],"#### HOSTNAME %s" % self._hostname())
        self.assertEqual(stdout[13],"#### USER %s" % self._user())
        self.assertTrue(stdout[14].startswith("#### START "))
        self.assertEqual(stdout[15],"#### CWD %s" % self.working_dir)
        self.assertEqual(stdout[16],"Takk!")
        self.assertEqual(stdout[17],"Wilkommen!")
        self.assertTrue(stdout[18].startswith("#### END "))
        self.assertEqual(stdout[19],"#### EXIT_CODE 0")
        self.assertEqual(stdout[20],"#### COMMAND Batch commands for Echo "
                         "string")
        self.assertEqual(stdout[21],"#### BATCH 3")
        self.assertEqual(stdout[22],"#### HOSTNAME %s" % self._hostname())
        self.assertEqual(stdout[23],"#### USER %s" % self._user())
        self.assertTrue(stdout[24].startswith("#### START "))
        self.assertEqual(stdout[25],"#### CWD %s" % self.working_dir)
        self.assertEqual(stdout[26],"Benvenuto!")
        self.assertTrue(stdout[27].startswith("#### END "))
        self.assertEqual(stdout[28],"#### EXIT_CODE 0")

    def test_pipelinetask_with_batched_scripts(self):
        """
        PipelineTask: run task with batched scripts
        """
        # Define a task with a command
        # Echoes text via shell command
        class EchoMany(PipelineTask):
            def init(self,*s):
                pass
            def setup(self):
                for s in self.args.s:
                    self.add_cmd(
                        PipelineScriptWrapper(
                            "Echo text",
                            """
                            echo {s}
                            """.format(s=s)))
        # Make a task instance
        task = EchoMany("Echo string",
                        "Hello!",
                        "Bonjour!",
                        "Takk!",
                        "Wilkommen!",
                        "Benvenuto!")
        # Check initial state
        self.assertEqual(task.args.s,
                         ("Hello!",
                          "Bonjour!",
                          "Takk!",
                          "Wilkommen!",
                          "Benvenuto!"))
        self.assertFalse(task.completed)
        self.assertEqual(task.exit_code,None)
        self.assertFalse(task.output)
        # Run the task with batches
        task.run(sched=self.sched,
                 working_dir=self.working_dir,
                 batch_size=2,
                 asynchronous=False)
        # Check final state
        self.assertTrue(task.completed)
        self.assertEqual(task.exit_code,0)
        self.assertFalse(task.output)
        # Check stdout
        # Should look like:
        # #### COMMAND Batch commands for Echo string
        # #### BATCH 1
        # #### HOSTNAME popov
        # #### USER pjb
        # #### START Thu Aug 17 08:38:14 BST 2017
        # #### CWD /tmp/dir
        # Hello!
        # Bonjour!
        # #### END Thu Aug 17 08:38:14 BST 2017
        # #### EXIT_CODE 0
        # #### COMMAND Batch commands for Echo string
        # #### BATCH 2
        # #### HOSTNAME popov
        # #### USER pjb
        # #### START Thu Aug 17 08:38:14 BST 2017
        # #### CWD /tmp/dir
        # Takk!
        # Wilkommen!
        # #### END Thu Aug 17 08:38:14 BST 2017
        # #### EXIT_CODE 0
        # #### COMMAND Batch commands for Echo string
        # #### BATCH 3
        # #### HOSTNAME popov
        # #### USER pjb
        # #### START Thu Aug 17 08:38:14 BST 2017
        # #### CWD /tmp/dir
        # Benvenuto!
        # #### END Thu Aug 17 08:38:14 BST 2017
        # #### EXIT_CODE 0
        stdout = task.stdout.split("\n")
        self.assertEqual(len(stdout),30) # 30 = 29 + trailing newline
        self.assertEqual(stdout[0],"#### COMMAND Batch commands for Echo "
                         "string")
        self.assertEqual(stdout[1],"#### BATCH 1")
        self.assertEqual(stdout[2],"#### HOSTNAME %s" % self._hostname())
        self.assertEqual(stdout[3],"#### USER %s" % self._user())
        self.assertTrue(stdout[4].startswith("#### START "))
        self.assertEqual(stdout[5],"#### CWD %s" % self.working_dir)
        self.assertEqual(stdout[6],"Hello!")
        self.assertEqual(stdout[7],"Bonjour!")
        self.assertTrue(stdout[8].startswith("#### END "))
        self.assertEqual(stdout[9],"#### EXIT_CODE 0")
        self.assertEqual(stdout[10],"#### COMMAND Batch commands for Echo "
                         "string")
        self.assertEqual(stdout[11],"#### BATCH 2")
        self.assertEqual(stdout[12],"#### HOSTNAME %s" % self._hostname())
        self.assertEqual(stdout[13],"#### USER %s" % self._user())
        self.assertTrue(stdout[14].startswith("#### START "))
        self.assertEqual(stdout[15],"#### CWD %s" % self.working_dir)
        self.assertEqual(stdout[16],"Takk!")
        self.assertEqual(stdout[17],"Wilkommen!")
        self.assertTrue(stdout[18].startswith("#### END "))
        self.assertEqual(stdout[19],"#### EXIT_CODE 0")
        self.assertEqual(stdout[20],"#### COMMAND Batch commands for Echo "
                         "string")
        self.assertEqual(stdout[21],"#### BATCH 3")
        self.assertEqual(stdout[22],"#### HOSTNAME %s" % self._hostname())
        self.assertEqual(stdout[23],"#### USER %s" % self._user())
        self.assertTrue(stdout[24].startswith("#### START "))
        self.assertEqual(stdout[25],"#### CWD %s" % self.working_dir)
        self.assertEqual(stdout[26],"Benvenuto!")
        self.assertTrue(stdout[27].startswith("#### END "))
        self.assertEqual(stdout[28],"#### EXIT_CODE 0")

    def test_pipelinetask_with_batch_limit(self):
        """
        PipelineTask: run task with batch limit
        """
        # Define a task with a command
        # Echoes text via shell command
        class EchoMany(PipelineTask):
            def init(self,*s):
                pass
            def setup(self):
                for s in self.args.s:
                    self.add_cmd(
                        PipelineCommandWrapper(
                            "Echo text","echo",s))
        # Make a task instance
        task = EchoMany("Echo string",
                        "Hello!",
                        "Bonjour!",
                        "Takk!",
                        "Wilkommen!",
                        "Benvenuto!")
        # Check initial state
        self.assertEqual(task.args.s,
                         ("Hello!",
                          "Bonjour!",
                          "Takk!",
                          "Wilkommen!",
                          "Benvenuto!"))
        self.assertFalse(task.completed)
        self.assertEqual(task.exit_code,None)
        self.assertFalse(task.output)
        # Run the task with batches
        task.run(sched=self.sched,
                 working_dir=self.working_dir,
                 batch_limit=3,
                 asynchronous=False)
        # Check final state
        self.assertTrue(task.completed)
        self.assertEqual(task.exit_code,0)
        self.assertFalse(task.output)
        # Check stdout
        # Should look like:
        # #### COMMAND Batch commands for Echo string
        # #### BATCH 1
        # #### HOSTNAME popov
        # #### USER pjb
        # #### START Thu Aug 17 08:38:14 BST 2017
        # #### CWD /tmp/dir
        # Hello!
        # Bonjour!
        # #### END Thu Aug 17 08:38:14 BST 2017
        # #### EXIT_CODE 0
        # #### COMMAND Batch commands for Echo string
        # #### BATCH 2
        # #### HOSTNAME popov
        # #### USER pjb
        # #### START Thu Aug 17 08:38:14 BST 2017
        # #### CWD /tmp/dir
        # Takk!
        # Wilkommen!
        # #### END Thu Aug 17 08:38:14 BST 2017
        # #### EXIT_CODE 0
        # #### COMMAND Batch commands for Echo string
        # #### BATCH 3
        # #### HOSTNAME popov
        # #### USER pjb
        # #### START Thu Aug 17 08:38:14 BST 2017
        # #### CWD /tmp/dir
        # Benvenuto!
        # #### END Thu Aug 17 08:38:14 BST 2017
        # #### EXIT_CODE 0
        stdout = task.stdout.split("\n")
        self.assertEqual(len(stdout),30) # 30 = 29 + trailing newline
        self.assertEqual(stdout[0],"#### COMMAND Batch commands for Echo "
                         "string")
        self.assertEqual(stdout[1],"#### BATCH 1")
        self.assertEqual(stdout[2],"#### HOSTNAME %s" % self._hostname())
        self.assertEqual(stdout[3],"#### USER %s" % self._user())
        self.assertTrue(stdout[4].startswith("#### START "))
        self.assertEqual(stdout[5],"#### CWD %s" % self.working_dir)
        self.assertEqual(stdout[6],"Hello!")
        self.assertEqual(stdout[7],"Bonjour!")
        self.assertTrue(stdout[8].startswith("#### END "))
        self.assertEqual(stdout[9],"#### EXIT_CODE 0")
        self.assertEqual(stdout[10],"#### COMMAND Batch commands for Echo "
                         "string")
        self.assertEqual(stdout[11],"#### BATCH 2")
        self.assertEqual(stdout[12],"#### HOSTNAME %s" % self._hostname())
        self.assertEqual(stdout[13],"#### USER %s" % self._user())
        self.assertTrue(stdout[14].startswith("#### START "))
        self.assertEqual(stdout[15],"#### CWD %s" % self.working_dir)
        self.assertEqual(stdout[16],"Takk!")
        self.assertEqual(stdout[17],"Wilkommen!")
        self.assertTrue(stdout[18].startswith("#### END "))
        self.assertEqual(stdout[19],"#### EXIT_CODE 0")
        self.assertEqual(stdout[20],"#### COMMAND Batch commands for Echo "
                         "string")
        self.assertEqual(stdout[21],"#### BATCH 3")
        self.assertEqual(stdout[22],"#### HOSTNAME %s" % self._hostname())
        self.assertEqual(stdout[23],"#### USER %s" % self._user())
        self.assertTrue(stdout[24].startswith("#### START "))
        self.assertEqual(stdout[25],"#### CWD %s" % self.working_dir)
        self.assertEqual(stdout[26],"Benvenuto!")
        self.assertTrue(stdout[27].startswith("#### END "))
        self.assertEqual(stdout[28],"#### EXIT_CODE 0")

    def test_pipelinetask_with_commands_as_command_instances(self):
        """
        PipelineTask: run task with commands specified as 'Command' instances
        """
        # Define a task with a command
        # Echoes text via shell command
        class EchoMany(PipelineTask):
            def init(self,*s):
                pass
            def setup(self):
                for s in self.args.s:
                    self.add_cmd("Echo text",
                                 Command("echo",s))
        # Make a task instance
        task = EchoMany("Echo string","Hello!","Goodbye!")
        # Check initial state
        self.assertEqual(task.args.s,("Hello!","Goodbye!"))
        self.assertFalse(task.completed)
        self.assertEqual(task.exit_code,None)
        self.assertFalse(task.output)
        # Run the task
        task.run(sched=self.sched,
                 working_dir=self.working_dir,
                 asynchronous=False)
        # Check final state
        self.assertTrue(task.completed)
        self.assertEqual(task.exit_code,0)
        self.assertFalse(task.output)
        # Check stdout
        # Should look like:
        # #### COMMAND Echo text
        # #### HOSTNAME popov
        # #### USER pjb
        # #### START Thu Aug 17 08:38:14 BST 2017
        # #### CWD /tmp/dir
        # Hello!
        # #### END Thu Aug 17 08:38:14 BST 2017
        # #### EXIT_CODE 0
        # #### COMMAND Echo text
        # #### HOSTNAME popov
        # #### USER pjb
        # #### START Thu Aug 17 08:38:14 BST 2017
        # #### CWD /tmp/dir
        # Goodbye!
        # #### END Thu Aug 17 08:38:14 BST 2017
        # #### EXIT_CODE 0
        stdout = task.stdout.split("\n")
        self.assertEqual(len(stdout),17) # 17 = 16 + trailing newline
        self.assertEqual(stdout[0],"#### COMMAND Echo text")
        self.assertEqual(stdout[1],"#### HOSTNAME %s" % self._hostname())
        self.assertEqual(stdout[2],"#### USER %s" % self._user())
        self.assertTrue(stdout[3].startswith("#### START "))
        self.assertEqual(stdout[4],"#### CWD %s" % self.working_dir)
        self.assertEqual(stdout[5],"Hello!")
        self.assertTrue(stdout[6].startswith("#### END "))
        self.assertEqual(stdout[7],"#### EXIT_CODE 0")
        self.assertEqual(stdout[8],"#### COMMAND Echo text")
        self.assertEqual(stdout[9],"#### HOSTNAME %s" % self._hostname())
        self.assertEqual(stdout[10],"#### USER %s" % self._user())
        self.assertTrue(stdout[11].startswith("#### START "))
        self.assertEqual(stdout[12],"#### CWD %s" % self.working_dir)
        self.assertEqual(stdout[13],"Goodbye!")
        self.assertTrue(stdout[14].startswith("#### END "))
        self.assertEqual(stdout[15],"#### EXIT_CODE 0")

    def test_pipelinetask_with_commands_as_scripts(self):
        """
        PipelineTask: run task with commands specified as scripts
        """
        # Define a task with a command
        # Echoes text via shell command
        class EchoMany(PipelineTask):
            def init(self,*s):
                pass
            def setup(self):
                for s in self.args.s:
                    self.add_cmd(
                        "Echo text",
                        """
                        # Script to echo supplied text
                        TEXT={s}
                        if [ 1 == 1 ] ; then
                          echo "$TEXT"
                        fi
                        """.format(s=s))
        # Make a task instance
        task = EchoMany("Echo string","Hello!","Goodbye!")
        # Check initial state
        self.assertEqual(task.args.s,("Hello!","Goodbye!"))
        self.assertFalse(task.completed)
        self.assertEqual(task.exit_code,None)
        self.assertFalse(task.output)
        # Run the task
        task.run(sched=self.sched,
                 working_dir=self.working_dir,
                 asynchronous=False)
        # Check final state
        self.assertTrue(task.completed)
        self.assertEqual(task.exit_code,0)
        self.assertFalse(task.output)
        # Check stdout
        # Should look like:
        # #### COMMAND Echo text
        # #### HOSTNAME popov
        # #### USER pjb
        # #### START Thu Aug 17 08:38:14 BST 2017
        # #### CWD /tmp/dir
        # Hello!
        # #### END Thu Aug 17 08:38:14 BST 2017
        # #### EXIT_CODE 0
        # #### COMMAND Echo text
        # #### HOSTNAME popov
        # #### USER pjb
        # #### START Thu Aug 17 08:38:14 BST 2017
        # #### CWD /tmp/dir
        # Goodbye!
        # #### END Thu Aug 17 08:38:14 BST 2017
        # #### EXIT_CODE 0
        stdout = task.stdout.split("\n")
        self.assertEqual(len(stdout),17) # 17 = 16 + trailing newline
        self.assertEqual(stdout[0],"#### COMMAND Echo text")
        self.assertEqual(stdout[1],"#### HOSTNAME %s" % self._hostname())
        self.assertEqual(stdout[2],"#### USER %s" % self._user())
        self.assertTrue(stdout[3].startswith("#### START "))
        self.assertEqual(stdout[4],"#### CWD %s" % self.working_dir)
        self.assertEqual(stdout[5],"Hello!")
        self.assertTrue(stdout[6].startswith("#### END "))
        self.assertEqual(stdout[7],"#### EXIT_CODE 0")
        self.assertEqual(stdout[8],"#### COMMAND Echo text")
        self.assertEqual(stdout[9],"#### HOSTNAME %s" % self._hostname())
        self.assertEqual(stdout[10],"#### USER %s" % self._user())
        self.assertTrue(stdout[11].startswith("#### START "))
        self.assertEqual(stdout[12],"#### CWD %s" % self.working_dir)
        self.assertEqual(stdout[13],"Goodbye!")
        self.assertTrue(stdout[14].startswith("#### END "))
        self.assertEqual(stdout[15],"#### EXIT_CODE 0")

    def test_pipelinetask_with_failing_command(self):
        """
        PipelineTask: run task with failing shell command
        """
        # Define a task with a command
        # Attempts to run a non-existant shell command
        class Nonexistant(PipelineTask):
            def init(self):
                pass
            def setup(self):
                self.add_cmd(
                    PipelineCommandWrapper(
                        "Nonexistant","./non_existant --help"))
        # Make a task instance
        task = Nonexistant("Will fail")
        # Check initial state
        self.assertFalse(task.completed)
        self.assertEqual(task.exit_code,None)
        self.assertFalse(task.output)
        # Run the task
        task.run(sched=self.sched,
                 working_dir=self.working_dir,
                 asynchronous=False)
        # Check final state
        self.assertTrue(task.completed)
        self.assertNotEqual(task.exit_code,0)
        self.assertFalse(task.output)
        # Check stdout
        # Should look like:
        # #### COMMAND Nonexistant
        # #### HOSTNAME popov
        # #### USER pjb
        # #### START Thu Aug 17 08:38:14 BST 2017
        # #### CWD /tmp/dir
        # #### END Thu Aug 17 08:38:14 BST 2017
        # #### EXIT_CODE 127
        stdout = task.stdout.split("\n")
        self.assertEqual(len(stdout),8) # 8 = 7 + trailing newline
        self.assertEqual(stdout[0],"#### COMMAND Nonexistant")
        self.assertEqual(stdout[1],"#### HOSTNAME %s" % self._hostname())
        self.assertEqual(stdout[2],"#### USER %s" % self._user())
        self.assertTrue(stdout[3].startswith("#### START "))
        self.assertEqual(stdout[4],"#### CWD %s" % self.working_dir)
        self.assertTrue(stdout[5].startswith("#### END "))
        self.assertEqual(stdout[6],"#### EXIT_CODE 127")

    def test_pipelinetask_with_command_using_runner_nslots(self):
        """
        PipelineTask: run task with shell command using nslots from runner
        """
        # Define a task with a command
        # Echoes number of slots set by runner
        class GetNSlots(PipelineTask):
            def init(self):
                pass
            def setup(self):
                self.add_cmd(
                    PipelineCommandWrapper(
                        "Get nslots",
                        "echo","NSLOTS=%s" % self.runner_nslots))
        # Make a task instance
        task = GetNSlots("Get number of slots")
        # Run the task with a custom job runner
        task.run(sched=self.sched,
                 working_dir=self.working_dir,
                 runner=SimpleJobRunner(nslots=8),
                 asynchronous=False)
        # Check final state
        self.assertTrue(task.completed)
        self.assertEqual(task.exit_code,0)
        self.assertFalse(task.output)
        # Check stdout
        # Should look like:
        # #### COMMAND Get nslots
        # #### HOSTNAME popov
        # #### USER pjb
        # #### START Thu Aug 17 08:38:14 BST 2017
        # #### CWD /tmp/dir
        # NSLOTS=8
        # #### END Thu Aug 17 08:38:14 BST 2017
        # #### EXIT_CODE 0
        stdout = task.stdout.split("\n")
        self.assertEqual(len(stdout),9) # 9 = 8 + trailing newline
        self.assertEqual(stdout[0],"#### COMMAND Get nslots")
        self.assertEqual(stdout[1],"#### HOSTNAME %s" % self._hostname())
        self.assertEqual(stdout[2],"#### USER %s" % self._user())
        self.assertTrue(stdout[3].startswith("#### START "))
        self.assertEqual(stdout[4],"#### CWD %s" % self.working_dir)
        self.assertEqual(stdout[5],"NSLOTS=8")
        self.assertTrue(stdout[6].startswith("#### END "))
        self.assertEqual(stdout[7],"#### EXIT_CODE 0")

    def test_pipelinetask_stdout(self):
        """
        PipelineTask: check stdout recovered from task
        """
        # Define a task with a command
        # Echoes text multiple times via shell command
        class MultipleEcho(PipelineTask):
            def init(self,s,n=1):
                pass
            def setup(self):
                for i in range(self.args.n):
                    self.add_cmd(
                        PipelineCommandWrapper(
                            "Echo text","echo",self.args.s))
        # Make a task instance
        task = MultipleEcho("Echo string 3 times","Hello!",3)
        # Check initial state
        self.assertFalse(task.completed)
        self.assertEqual(task.exit_code,None)
        self.assertFalse(task.output)
        # Run the task
        task.run(sched=self.sched,
                 working_dir=self.working_dir,
                 asynchronous=False)
        # Check final state
        self.assertTrue(task.completed)
        self.assertEqual(task.exit_code,0)
        self.assertFalse(task.output)
        # Check stdout
        # Should look like:
        # #### COMMAND Echo text
        # #### HOSTNAME popov
        # #### USER pjb
        # #### START Thu Aug 17 08:38:14 BST 2017
        # #### CWD /tmp/dir
        # Hello!
        # #### END Thu Aug 17 08:38:14 BST 2017
        # #### EXIT_CODE 0
        # ...x three times
        print(task.stdout)
        stdout = task.stdout.split("\n")
        self.assertEqual(len(stdout),25) # 25 = 24 + trailing newline
        for i in range(3):
            self.assertEqual(stdout[0+i*8],"#### COMMAND Echo text")
            self.assertEqual(stdout[1+i*8],"#### HOSTNAME %s" % self._hostname())
            self.assertEqual(stdout[2+i*8],"#### USER %s" % self._user())
            self.assertTrue(stdout[3+i*8].startswith("#### START "))
            self.assertEqual(stdout[4+i*8],"#### CWD %s" % self.working_dir)
            self.assertEqual(stdout[5+i*8],"Hello!")
            self.assertTrue(stdout[6+i*8].startswith("#### END "))
            self.assertEqual(stdout[7+i*8],"#### EXIT_CODE 0")

    def test_pipelinetask_stderr(self):
        """
        PipelineTask: check stderr recovered from task
        """
        # Define a task with a command
        # Echoes text multiple times via shell command
        class MultipleEchoToStderr(PipelineTask):
            def init(self,s,n=1):
                pass
            def setup(self):
                for i in range(self.args.n):
                    self.add_cmd(
                        PipelineCommandWrapper(
                            "Echo text","echo",self.args.s,">&2"))
        # Make a task instance
        task = MultipleEchoToStderr("Echo string 3 times","Hello!",3)
        # Check initial state
        self.assertFalse(task.completed)
        self.assertEqual(task.exit_code,None)
        self.assertFalse(task.output)
        # Run the task
        task.run(sched=self.sched,
                 working_dir=self.working_dir,
                 asynchronous=False)
        # Check final state
        self.assertTrue(task.completed)
        self.assertEqual(task.exit_code,0)
        self.assertFalse(task.output)
        # Check stderr
        # Should look like:
        # Hello!
        # ...x three times
        print(task.stderr)
        stderr = task.stderr.split("\n")
        if JOB_HAS_ERR:
            self.assertEqual(len(stderr),4) # 4 = 3 + trailing newline
            for i in range(3):
                self.assertEqual(stderr[i],"Hello!")
        else:
            self.assertEqual(stderr,[""])

    def test_pipelinetask_invoke_fail(self):
        """
        PipelineTask: check task invoking 'fail' method
        """
        # Define a task which invokes 'fail'
        class FailingTask(PipelineTask):
            def init(self):
                pass
            def setup(self):
                self.fail(message="Invoked fail method",
                          exit_code=123)
        # Make a task instance
        task = FailingTask("This will fail")
        # Check initial state
        self.assertFalse(task.completed)
        self.assertEqual(task.exit_code,None)
        self.assertFalse(task.output)
        # Run the task
        task.run(sched=self.sched,
                 working_dir=self.working_dir,
                 poll_interval=0.5,
                 asynchronous=False)
        # Check final state
        self.assertTrue(task.completed)
        self.assertEqual(task.exit_code,123)
        self.assertFalse(task.output)
        self.assertEqual(task.stdout,"")

    def test_pipelinetask_cloudpickle(self):
        """
        PipelineTask: check serialization using 'cloudpickle'
        """
        # Define a task with a command
        # Echoes text via shell command
        class Echo(PipelineTask):
            def init(self,s):
                pass
            def setup(self):
                self.add_cmd(
                    PipelineCommandWrapper(
                        "Echo text","echo",self.args.s))
        # Make a task instance
        task = Echo("Echo string","Hello!")
        # Pickle it
        pickled = cloudpickle.dumps(task)
        # Unpickle it
        unpickled = cloudpickle.loads(pickled)

    def test_pipelinetask_no_conda_dependencies(self):
        """
        PipelineTask: check when no conda dependencies are declared
        """
        # Define a task without conda dependencies
        class NoCondaDeps(PipelineTask):
            def init(self,fq):
                pass
            def setup(self):
                self.add_cmd(PipelineCommandWrapper(
                    "Run FastQC","fastqc",self.args.fq))
        # Make a task instance
        task = NoCondaDeps("Test","Sample1_S1_R1_001.fastq.gz")
        # Check conda dependencies
        self.assertEqual(task.conda_dependencies,[])
        # Check conda environment name
        self.assertEqual(task.conda_env_name,None)

    def test_pipelinetask_conda_dependencies(self):
        """
        PipelineTask: check declaring conda dependencies
        """
        # Define a task with conda dependencies
        class WithCondaDeps(PipelineTask):
            def init(self,fq):
                self.conda("fastqc=0.11.3")
                self.conda("fastq-screen=0.14.0",
                           "bowtie=1.2.3")
            def setup(self):
                self.add_cmd(
                    PipelineCommandWrapper(
                        "Run FastQC","fastqc",self.args.fq))
        # Make a task instance
        task = WithCondaDeps("Test","Sample1_S1_R1_001.fastq.gz")
        # Check conda dependencies
        self.assertEqual(task.conda_dependencies,
                         ["fastqc=0.11.3",
                          "fastq-screen=0.14.0",
                          "bowtie=1.2.3"])
        # Check conda environment name
        self.assertEqual(task.conda_env_name,
                         "bowtie@1.2.3+fastq-screen@0.14.0+fastqc@0.11.3")

    def test_pipelinetask_setup_conda_env(self):
        """
        PipelineTask: set up conda environment for dependencies
        """
        # Define a task with conda dependencies
        class WithCondaDeps(PipelineTask):
            def init(self,fq):
                self.conda("fastqc=0.11.3")
                self.conda("fastq-screen=0.14.0",
                           "bowtie=1.2.3")
            def setup(self):
                self.add_cmd(
                    PipelineCommandWrapper(
                        "Run FastQC","fastqc",self.args.fq))
        # Create a mock conda installation
        conda_dir = MockConda.create(os.path.join(self.working_dir,"__conda"))
        conda_ = os.path.join(conda_dir,"bin","conda")
        # Make a task instance
        task = WithCondaDeps("Test","Sample1_S1_R1_001.fastq.gz")
        # Setup conda environment
        conda_env = task.setup_conda_env(conda_)
        # Check conda environment
        self.assertEqual(
            conda_env,
            os.path.join(conda_dir,
                         "envs",
                         "bowtie@1.2.3+fastq-screen@0.14.0+fastqc@0.11.3"))
        self.assertTrue(os.path.exists(conda_env))

    def test_pipelinetask_setup_conda_env_custom_location(self):
        """
        PipelineTask: set up conda environment for dependencies in custom location
        """
        # Define a task with conda dependencies
        class WithCondaDeps(PipelineTask):
            def init(self,fq):
                self.conda("fastqc=0.11.3")
                self.conda("fastq-screen=0.14.0",
                           "bowtie=1.2.3")
            def setup(self):
                self.add_cmd(
                    PipelineCommandWrapper(
                        "Run FastQC","fastqc",self.args.fq))
        # Create a mock conda installation
        conda_dir = MockConda.create(os.path.join(self.working_dir,"__conda"))
        conda_ = os.path.join(conda_dir,"bin","conda")
        alternative_env_dir = os.path.join(self.working_dir,
                                           '__my_conda_envs')
        os.makedirs(alternative_env_dir)
        # Make a task instance
        task = WithCondaDeps("Test","Sample1_S1_R1_001.fastq.gz")
        # Setup conda environment
        conda_env = task.setup_conda_env(conda_,
                                         env_dir=alternative_env_dir)
        # Check conda environment
        self.assertEqual(
            conda_env,
            os.path.join(alternative_env_dir,
                         "bowtie@1.2.3+fastq-screen@0.14.0+fastqc@0.11.3"))
        self.assertTrue(os.path.exists(conda_env))

    def test_pipelinetask_run_with_conda_dependencies_enabled(self):
        """
        PipelineTask: run task with conda dependencies enabled
        """
        # Define a task with conda dependencies
        class WithCondaDeps(PipelineTask):
            def init(self,fq):
                self.conda("fastqc=0.11.3")
                self.conda("fastq-screen=0.14.0",
                           "bowtie=1.2.3")
            def setup(self):
                self.add_cmd(
                    PipelineCommandWrapper(
                        "Run FastQC","fastqc",self.args.fq))
        # Make a task instance
        task = WithCondaDeps("Test","Sample1_S1_R1_001.fastq.gz")
        # Create a mock conda installation
        conda_dir = MockConda.create(os.path.join(self.working_dir,"__conda"))
        conda_ = os.path.join(conda_dir,"bin","conda")
        # Run the task
        task.run(sched=self.sched,
                 enable_conda=True,
                 conda=conda_,
                 working_dir=self.working_dir,
                 poll_interval=0.5,
                 asynchronous=False)
        # Check final state
        self.assertTrue(task.completed)
        self.assertEqual(task.exit_code,0)
        # Check conda environment
        self.assertTrue(
            os.path.exists(os.path.join(
                conda_dir,
                "envs",
                "bowtie@1.2.3+fastq-screen@0.14.0+fastqc@0.11.3")))

    def test_pipelinetask_run_with_conda_dependencies_disabled(self):
        """
        PipelineTask: run task with conda dependencies disabled
        """
        # Define a task with conda dependencies
        class WithCondaDeps(PipelineTask):
            def init(self,fq):
                self.conda("fastqc=0.11.3")
                self.conda("fastq-screen=0.14.0",
                           "bowtie=1.2.3")
            def setup(self):
                self.add_cmd(
                    PipelineCommandWrapper(
                        "Run FastQC","fastqc",self.args.fq))
        # Make a task instance
        task = WithCondaDeps("Test","Sample1_S1_R1_001.fastq.gz")
        # Create a mock conda installation
        conda_dir = MockConda.create(os.path.join(self.working_dir,"__conda"))
        conda_ = os.path.join(conda_dir,"bin","conda")
        # Create a mock FastQC instance
        fastq_bin_dir = os.path.join(self.working_dir,'__apps','bin')
        os.makedirs(fastq_bin_dir)
        _Mock.fastqc(fastq_bin_dir)
        os.environ['PATH'] = os.environ['PATH'] + os.pathsep + fastq_bin_dir
        # Run the task
        task.run(sched=self.sched,
                 enable_conda=False,
                 conda=conda_,
                 working_dir=self.working_dir,
                 poll_interval=0.5,
                 asynchronous=False)
        # Check final state
        self.assertTrue(task.completed)
        self.assertEqual(task.exit_code,0)
        # Check conda environment
        self.assertFalse(
            os.path.exists(os.path.join(
                conda_dir,
                "envs",
                "bowtie@1.2.3+fastq-screen@0.14.0+fastqc@0.11.3")))

class TestPipelineFunctionTask(unittest.TestCase):

    def setUp(self):
        # Set up a scheduler
        self.sched = SimpleScheduler(poll_interval=0.5)
        self.sched.start()
        # Make a temporary working dir
        self.working_dir = tempfile.mkdtemp(
            suffix='TestPipeline')

    def tearDown(self):
        # Stop the scheduler
        if self.sched is not None:
            self.sched.stop()
        # Remove temp dir
        if os.path.exists(self.working_dir):
            shutil.rmtree(self.working_dir)

    def test_pipelinefunctiontask(self):
        """
        PipelineFunctionTask: run task wrapping function call
        """
        # Define a task with a function call
        class Hello(PipelineFunctionTask):
            def init(self,name):
                pass
            def setup(self):
                self.add_call("Emit greeting",
                              self.hello,
                              self.args.name)
            def hello(self,name):
                return "Hello %s!" % name
        # Make a task instance
        task = Hello("Hello world","World")
        # Check initial state
        self.assertEqual(task.args.name,"World")
        self.assertFalse(task.completed)
        self.assertEqual(task.exit_code,None)
        self.assertEqual(task.result(),None)
        self.assertFalse(task.output)
        # Run the task
        task.run(sched=self.sched,
                 working_dir=self.working_dir,
                 poll_interval=0.1,
                 asynchronous=False)
        # Check final state
        self.assertTrue(task.completed)
        self.assertEqual(task.exit_code,0)
        self.assertEqual(task.result(),["Hello World!"])
        self.assertFalse(task.output)

    def test_pipelinefunctiontask_with_failing_call(self):
        """
        PipelineFunctionTask: run task with failing function call
        """
        # Define a task with a function call
        # that raises an exception
        class RaisesException(PipelineFunctionTask):
            def init(self):
                pass
            def setup(self):
                self.add_call("Raises an exception",
                              self.raises_exception)
            def raises_exception(self):
                raise Exception("Exception is raised")
        # Make a task instance
        task = RaisesException("Will raise exception")
        # Check initial state
        self.assertFalse(task.completed)
        self.assertEqual(task.exit_code,None)
        self.assertEqual(task.result(),None)
        self.assertFalse(task.output)
        # Run the task
        task.run(sched=self.sched,
                 working_dir=self.working_dir,
                 poll_interval=0.1,
                 asynchronous=False)
        # Check final state
        self.assertTrue(task.completed)
        self.assertNotEqual(task.exit_code,0)
        self.assertEqual(task.result(),None)
        self.assertFalse(task.output)

    def test_pipelinefunctiontask_with_runner_nslots(self):
        """
        PipelineFunctionTask: run task wrapping function call with runner_nslots
        """
        # Define a task with a function call
        class GetNSlots(PipelineFunctionTask):
            def init(self):
                pass
            def setup(self):
                self.add_call("Get nslots",
                              self.get_nslots)
            def get_nslots(self):
                return "NSLOTS=%s" % self.runner_nslots
        # Make a task instance
        task = GetNSlots("Get number of slots")
        # Run the task with a custom job runner
        task.run(sched=self.sched,
                 working_dir=self.working_dir,
                 poll_interval=0.1,
                 runner=SimpleJobRunner(nslots=8),
                 asynchronous=False)
        # Check final state
        self.assertTrue(task.completed)
        self.assertEqual(task.exit_code,0)
        self.assertEqual(task.result(),["NSLOTS=8"])
        self.assertFalse(task.output)

    def test_pipelinefunctiontask_cloudpickle(self):
        """
        PipelineFunctionTask: check serialization using 'cloudpickle'
        """
        # Define a task with a function call
        class Hello(PipelineFunctionTask):
            def init(self,name):
                pass
            def setup(self):
                self.add_call("Emit greeting",
                              self.hello,
                              self.args.name)
            def hello(self,name):
                return "Hello %s!" % name
        # Make a task instance
        task = Hello("Hello world","World")
        # Pickle it
        pickled = cloudpickle.dumps(task)
        # Unpickle it
        unpickled = cloudpickle.loads(pickled)

class TestPipelineCommand(unittest.TestCase):

    def setUp(self):
        # Make a temporary working dir
        self.working_dir = tempfile.mkdtemp(
            suffix='TestPipelineCommand')
        # Unset the MODULEPATH env var, if found
        if 'MODULEPATH' in os.environ:
            self.modulepath = os.environ['MODULEPATH']
            os.environ['MODULEPATH'] = ''
        else:
            self.modulepath = None

    def tearDown(self):
        # Remove temp dir
        if os.path.exists(self.working_dir):
            shutil.rmtree(self.working_dir)
        # Restore the MODULEPATH env var
        if self.modulepath:
            os.environ['MODULEPATH'] = self.modulepath

    def test_pipelinecommand(self):
        """
        PipelineCommand: check command and wrapper script
        """
        # Subclass PipelineCommand
        class EchoCmd(PipelineCommand):
            def init(self,txt):
                self._txt = txt
            def cmd(self):
                return Command(
                    "echo",
                    self._txt)
        # Make an instance
        cmd = EchoCmd("hello there")
        # Check name
        self.assertEqual(cmd.name(),"echocmd")
        # Check command
        self.assertEqual(str(cmd.cmd()),"echo hello there")
        # Check wrapper script file
        script_file = cmd.make_wrapper_script(
            scripts_dir=self.working_dir)
        self.assertTrue(os.path.isfile(script_file))
        self.assertEqual(os.path.dirname(script_file),
                         self.working_dir)
        with open(script_file,'rt') as fp:
            self.assertEqual(fp.read(),
                             "#!/bin/bash\n"
                             "echo \"#### COMMAND EchoCmd\"\n"
                             "echo \"#### HOSTNAME $HOSTNAME\"\n"
                             "echo \"#### USER $USER\"\n"
                             "echo \"#### START $(date)\"\n"
                             "echo \"#### CWD $(pwd)\"\n"
                             "echo 'hello there'\n"
                             "exit_code=$?\n"
                             "echo \"#### END $(date)\"\n"
                             "echo \"#### EXIT_CODE $exit_code\"\n"
                             "exit $exit_code")

    def test_pipelinecommand_with_modules(self):
        """
        PipelineCommand: check command and wrapper script with modules
        """
        # Subclass PipelineCommand
        class EchoCmd(PipelineCommand):
            def init(self,txt):
                self._txt = txt
            def cmd(self):
                return Command(
                    "echo",
                    self._txt)
        # Make an instance
        cmd = EchoCmd("hello there")
        # Check name
        self.assertEqual(cmd.name(),"echocmd")
        # Check command
        self.assertEqual(str(cmd.cmd()),"echo hello there")
        # Check wrapper script file
        script_file = cmd.make_wrapper_script(
            scripts_dir=self.working_dir,
            envmodules=('apps/fastq-screen/0.13.0',
                        'apps/fastqc/0.11.8',))
        self.assertTrue(os.path.isfile(script_file))
        self.assertEqual(os.path.dirname(script_file),
                         self.working_dir)
        with open(script_file,'rt') as fp:
            self.assertEqual(fp.read(),
                             "#!/bin/bash --login\n"
                             "echo \"#### COMMAND EchoCmd\"\n"
                             "echo \"#### HOSTNAME $HOSTNAME\"\n"
                             "echo \"#### USER $USER\"\n"
                             "echo \"#### START $(date)\"\n"
                             "module load apps/fastq-screen/0.13.0\n"
                             "if [ $? -ne 0 ] ; then\n"
                             "  echo Failed to load environment module >&2\n"
                             "  exit 1\n"
                             "fi\n"
                             "module load apps/fastqc/0.11.8\n"
                             "if [ $? -ne 0 ] ; then\n"
                             "  echo Failed to load environment module >&2\n"
                             "  exit 1\n"
                             "fi\n"
                             "echo \"#### CWD $(pwd)\"\n"
                             "echo 'hello there'\n"
                             "exit_code=$?\n"
                             "echo \"#### END $(date)\"\n"
                             "echo \"#### EXIT_CODE $exit_code\"\n"
                             "exit $exit_code")

    def test_pipelinecommand_with_conda_env(self):
        """
        PipelineCommand: check command and wrapper script with conda envs
        """
        # Subclass PipelineCommand
        class EchoCmd(PipelineCommand):
            def init(self,txt):
                self._txt = txt
            def cmd(self):
                return Command(
                    "echo",
                    self._txt)
        # Make an instance
        cmd = EchoCmd("hello there")
        # Check name
        self.assertEqual(cmd.name(),"echocmd")
        # Check command
        self.assertEqual(str(cmd.cmd()),"echo hello there")
        # Check wrapper script file
        script_file = cmd.make_wrapper_script(
            scripts_dir=self.working_dir,
            conda="/usr/local/conda/bin/conda",
            conda_env="/work/__conda/envs")
        self.assertTrue(os.path.isfile(script_file))
        self.assertEqual(os.path.dirname(script_file),
                         self.working_dir)
        with open(script_file,'rt') as fp:
            self.assertEqual(
                fp.read(),
                "#!/bin/bash\n"
                "echo \"#### COMMAND EchoCmd\"\n"
                "echo \"#### HOSTNAME $HOSTNAME\"\n"
                "echo \"#### USER $USER\"\n"
                "echo \"#### START $(date)\"\n"
                "echo source /usr/local/conda/bin/activate /work/__conda/envs\n"
                "source /usr/local/conda/bin/activate /work/__conda/envs\n"
                "if [ $? -ne 0 ] ; then\n"
                "  echo Failed to activate conda environment >&2\n"
                "  exit 1\n"
                "fi\n"
                "echo \"#### CWD $(pwd)\"\n"
                "echo 'hello there'\n"
                "exit_code=$?\n"
                "echo \"#### END $(date)\"\n"
                "echo \"#### EXIT_CODE $exit_code\"\n"
                "exit $exit_code")

    def test_pipelinecommand_with_working_dir(self):
        """
        PipelineCommand: check command and wrapper script with working dir
        """
        # Subclass PipelineCommand
        class EchoCmd(PipelineCommand):
            def init(self,txt):
                self._txt = txt
            def cmd(self):
                return Command(
                    "echo",
                    self._txt)
        # Make an instance
        cmd = EchoCmd("hello there")
        # Check name
        self.assertEqual(cmd.name(),"echocmd")
        # Check command
        self.assertEqual(str(cmd.cmd()),"echo hello there")
        # Check wrapper script file
        script_file = cmd.make_wrapper_script(
            scripts_dir=self.working_dir,
            working_dir="/tmp/command/wd")
        self.assertTrue(os.path.isfile(script_file))
        self.assertEqual(os.path.dirname(script_file),
                         self.working_dir)
        with open(script_file,'rt') as fp:
            self.assertEqual(fp.read(),
                             "#!/bin/bash\n"
                             "echo \"#### COMMAND EchoCmd\"\n"
                             "echo \"#### HOSTNAME $HOSTNAME\"\n"
                             "echo \"#### USER $USER\"\n"
                             "echo \"#### START $(date)\"\n"
                             "cd /tmp/command/wd\n"
                             "echo \"#### CWD $(pwd)\"\n"
                             "echo 'hello there'\n"
                             "exit_code=$?\n"
                             "echo \"#### END $(date)\"\n"
                             "echo \"#### EXIT_CODE $exit_code\"\n"
                             "exit $exit_code")

class TestPipelineCommandWrapper(unittest.TestCase):

    def test_piplinecommandwrapper(self):
        """
        PipelineCommandWrapper: check generated command
        """
        # Make a pipeline command wrapper
        cmd = PipelineCommandWrapper("Echo text","echo","hello")
        # Check name and generated command
        self.assertEqual(cmd.name(),"echo_text")
        self.assertEqual(str(cmd.cmd()),"echo hello")
        # Add argument and check updated command
        cmd.add_args("there")
        self.assertEqual(str(cmd.cmd()),"echo hello there")

class TestPipelineScriptWrapper(unittest.TestCase):

    def setUp(self):
        # Make a temporary working dir
        self.working_dir = tempfile.mkdtemp(
            suffix='TestPipelineScriptWrapper')
        # Unset the MODULEPATH env var, if found
        if 'MODULEPATH' in os.environ:
            self.modulepath = os.environ['MODULEPATH']
            os.environ['MODULEPATH'] = ''
        else:
            self.modulepath = None

    def tearDown(self):
        # Remove temp dir
        if os.path.exists(self.working_dir):
            shutil.rmtree(self.working_dir)
        # Restore the MODULEPATH env var
        if self.modulepath:
            os.environ['MODULEPATH'] = self.modulepath

    def test_pipelinescriptwrapper(self):
        """
        PipelineScriptWrapper: check generated wrapper script
        """
        # Make a pipeline script wrapper
        script = PipelineScriptWrapper(
            "Echo text",
            """
            # Example script
            TEXT={txt}
            if [ 1 == 1 ] ; then
              echo "$TEXT"
            fi
            # Finished
            """.format(txt="hello there"))
        # Check name and generated wrapper script
        self.assertEqual(script.name(),"echo_text")
        script_file = script.make_wrapper_script(
            scripts_dir=self.working_dir)
        self.assertTrue(os.path.isfile(script_file))
        self.assertEqual(os.path.dirname(script_file),
                         self.working_dir)
        with open(script_file,'rt') as fp:
            self.assertEqual(fp.read(),
                             "#!/bin/bash\n"
                             "echo \"#### COMMAND Echo text\"\n"
                             "echo \"#### HOSTNAME $HOSTNAME\"\n"
                             "echo \"#### USER $USER\"\n"
                             "echo \"#### START $(date)\"\n"
                             "echo \"#### CWD $(pwd)\"\n"
                             "# Example script\n"
                             "TEXT=hello there\n"
                             "if [ 1 == 1 ] ; then\n"
                             "  echo \"$TEXT\"\n"
                             "fi\n"
                             "# Finished\n"
                             "exit_code=$?\n"
                             "echo \"#### END $(date)\"\n"
                             "echo \"#### EXIT_CODE $exit_code\"\n"
                             "exit $exit_code")

    def test_pipelinescriptwrapper_with_blocks(self):
        """
        PipelineScriptWrapper: check for multiple script blocks
        """
        # Make a pipeline script wrapper
        script = PipelineScriptWrapper(
            "Echo text",
            """
            # Block 1
            TEXT={txt}
            """.format(txt="hello there"),
            """
            # Block 2
            if [ 1 == 1 ] ; then
              echo "$TEXT"
            fi
            """)
        # Check name and generated wrapper script
        self.assertEqual(script.name(),"echo_text")
        script_file = script.make_wrapper_script(
            scripts_dir=self.working_dir)
        self.assertTrue(os.path.isfile(script_file))
        self.assertEqual(os.path.dirname(script_file),
                         self.working_dir)
        with open(script_file,'rt') as fp:
            self.assertEqual(fp.read(),
                             "#!/bin/bash\n"
                             "echo \"#### COMMAND Echo text\"\n"
                             "echo \"#### HOSTNAME $HOSTNAME\"\n"
                             "echo \"#### USER $USER\"\n"
                             "echo \"#### START $(date)\"\n"
                             "echo \"#### CWD $(pwd)\"\n"
                             "{\n"
                             "    # Block 1\n"
                             "    TEXT=hello there\n"
                             "} && {\n"
                             "    # Block 2\n"
                             "    if [ 1 == 1 ] ; then\n"
                             "      echo \"$TEXT\"\n"
                             "    fi\n"
                             "}\n"
                             "exit_code=$?\n"
                             "echo \"#### END $(date)\"\n"
                             "echo \"#### EXIT_CODE $exit_code\"\n"
                             "exit $exit_code")

class TestBaseParam(unittest.TestCase):

    def test_baseparam_uuid(self):
        """
        BaseParam: check UUID
        """
        p = BaseParam()
        self.assertNotEqual(p.uuid,None)

    def test_associated_task(self):
        """
        BaseParam: check associated task
        """
        # Define a simple task
        class SimpleTask(PipelineTask):
            def init(self,x):
                pass
            def setup(self):
                print(self.args.x)
        t = SimpleTask("Simple task",x=12)
        # Create and test associated task in a BaseParam
        p = BaseParam()
        self.assertEqual(p.associated_task_id,None)
        p.associate_task(t)
        self.assertEqual(p.associated_task_id,t.id())

class TestPipelineParam(unittest.TestCase):

    def test_pipelineparam_no_type(self):
        """
        PipelineParam: no type specified
        """
        # No initial value
        p = PipelineParam()
        self.assertEqual(p.value,None)
        p.set("abc")
        self.assertEqual(p.value,"abc")
        p.set(123)
        self.assertEqual(p.value,123)
        self.assertEqual(repr(p),"PipelineParam(value='123')")

    def test_pipelineparam_with_initial_value(self):
        """
        PipelineParam: initial value supplied
        """
        # Specify an initial value
        p = PipelineParam("abc")
        self.assertEqual(p.value,"abc")
        p = PipelineParam(value="def")
        self.assertEqual(p.value,"def")
        self.assertEqual(repr(p),"PipelineParam(value='def')")

    def test_pipelineparam_with_type(self):
        """
        PipelineParam: type function supplied
        """
        # Specify type function as 'str'
        p = PipelineParam(type=str)
        p.set("abc")
        self.assertEqual(p.value,"abc")
        p.set(123)
        self.assertEqual(p.value,"123")
        self.assertEqual(repr(p),
                         "PipelineParam(value='123',type='%s')" % str(str))
        # Specify type function as 'int'
        p = PipelineParam(type=int)
        p.set(123)
        self.assertEqual(p.value,123)
        p.set("123")
        self.assertEqual(p.value,123)
        self.assertEqual(repr(p),
                         "PipelineParam(value='123',type='%s')" % str(int))
        # Exception for bad value
        p.set("abc")
        self.assertRaises(ValueError,lambda: p.value)
        # Specify type function as 'float'
        p = PipelineParam(type=float)
        p.set(1.23)
        self.assertEqual(p.value,1.23)
        p.set("1.23")
        self.assertEqual(p.value,1.23)
        self.assertEqual(repr(p),
                         "PipelineParam(value='1.23',type='%s')" % str(float))
        # Exception for bad value
        p.set("abc")
        self.assertRaises(ValueError,lambda: p.value)

    def test_pipelineparam_with_default(self):
        """
        PipelineParam: default function supplied
        """
        # Specify type function as 'str'
        p = PipelineParam(default=lambda: "default")
        self.assertEqual(p.value,"default")
        p.set("abc")
        self.assertEqual(p.value,"abc")

    def test_pipelineparam_with_name(self):
        """
        PipelineParam: name supplied
        """
        # Specify a name
        p = PipelineParam(name="my_param")
        self.assertEqual(p.name,"my_param")
        self.assertEqual(repr(p),
                         "PipelineParam(name='my_param',value='None')")

    def test_pipelineparam_handles_None_value(self):
        """
        PipelineParam: handle 'None' value
        """
        # Untyped parameter
        self.assertEqual(PipelineParam().value,None)
        # Typed parameter
        self.assertEqual(PipelineParam(type=str).value,None)

    def test_pipelineparam_replacement(self):
        """
        PipelineParam: replace value with another parameter
        """
        # Create initial parameter
        p = PipelineParam(name="param1",value=1)
        self.assertEqual(p.value,1)
        # Create a new parameter
        pp = PipelineParam(name="param2",value=2)
        self.assertEqual(pp.value,2)
        # Replace first parameter with the second
        p.replace_with(pp)
        self.assertEqual(p.value,2)
        # Change value of second parameter
        pp.set(3)
        self.assertEqual(pp.value,3)
        # Value of first parameter should also have changed
        self.assertEqual(p.value,3)

    def test_pipelineparam_cascade_value_resolution(self):
        """
        PipelineParam: set value to another parameter
        """
        # Create initial parameter
        p = PipelineParam(name="param1",value=1)
        self.assertEqual(p.value,1)
        # Create a new parameter
        pp = PipelineParam(name="param2",value=2)
        self.assertEqual(pp.value,2)
        # Set first parameter to the second
        p.set(pp)
        self.assertEqual(p.value,2)
        # Change value of second parameter
        pp.set(3)
        self.assertEqual(pp.value,3)
        # Value of first parameter should also have changed
        self.assertEqual(p.value,3)

class TestFileCollector(unittest.TestCase):

    def setUp(self):
        # Make a temporary working dir
        self.working_dir = tempfile.mkdtemp(
            suffix='TestFileCollector')

    def tearDown(self):
        # Remove temp dir
        if os.path.exists(self.working_dir):
            shutil.rmtree(self.working_dir)

    def test_filecollector(self):
        """
        FileCollector: collects files matching pattern
        """
        # Set up collectors
        all_files = FileCollector(self.working_dir,"*")
        txt_files = FileCollector(self.working_dir,"*.txt")
        # Put some files in
        file_list = [os.path.join(self.working_dir,f)
                     for f in ["test1.txt","test.fq",
                               "test.r1.fastq","test.r2.fastq"]]
        file_list.sort()
        for f in file_list:
            with open(os.path.join(self.working_dir,f),'w') as fp:
                fp.write("")
        # Check each collection
        self.assertEqual(len(all_files),4)
        for f1,f2 in zip(all_files,file_list):
            self.assertEqual(f1,f2)
        self.assertEqual(len(txt_files),1)
        self.assertEqual(list(txt_files),
                         [os.path.join(self.working_dir,"test1.txt")])

class TestListParam(unittest.TestCase):
    """
    Tests for the 'ListParam' class
    """
    def test_listparam_init(self):
        """
        ListParam: check initialisation
        """
        # Empty list
        self.assertEqual(ListParam().value,[])
        # Non-empty iterable
        self.assertEqual(ListParam((1,2,"hello",3)).value,
                         [1,2,"hello",3])

    def test_listparam_append(self):
        """
        ListParam: check 'append' method
        """
        l = ListParam((1,2,"hello",3))
        l.append(4)
        self.assertEqual(l.value,[1,2,"hello",3,4])

    def test_listparam_extend(self):
        """
        ListParam: check 'extend' method
        """
        l = ListParam((1,2,"hello",3))
        l.extend((4,5,"goodbye"))
        self.assertEqual(l.value,[1,2,"hello",3,4,5,"goodbye"])

    def test_listparam_len(self):
        """
        ListParam: check 'len' functionality
        """
        self.assertEqual(len(ListParam()),0)
        self.assertEqual(len(ListParam((1,2,"hello",3))),4)

    def test_listparam_with_params(self):
        """
        ListParam: handle pipeline parameters as items
        """
        l = ListParam()
        l.append(PipelineParam(value=1))
        l.append(PipelineParam(value=2))
        l.append(PipelineParam(value="hello"))
        l.append(PipelineParam(value=3))
        self.assertEqual(l.value,[1,2,"hello",3])

class TestPathJoinParam(unittest.TestCase):
    """
    Tests for the 'PathJoinParam' class
    """
    def test_pathjoinparam_with_strings(self):
        """
        PathJoinParam: handle static strings
        """
        p1 = "/mnt/data"
        p2 = "test.txt"
        pth = PathJoinParam(p1,p2)
        self.assertEqual(pth.value,"/mnt/data/test.txt")

    def test_pathjoinparam_with_pipelineparams(self):
        """
        PathJoinParam: handle PipelineParams
        """
        p1 = PipelineParam(value="/mnt/data")
        p2 = PipelineParam(value="test.txt")
        pth = PathJoinParam(p1,p2)
        self.assertEqual(pth.value,"/mnt/data/test.txt")
        p2.set("updated.txt")
        self.assertEqual(pth.value,"/mnt/data/updated.txt")

    def test_pathjoinparam_with_mixed_input(self):
        """
        PathJoinParam: handle mixture of string and PipelineParam
        """
        p1 = PipelineParam(value="/mnt/data")
        p2 = "test.txt"
        pth = PathJoinParam(p1,p2)
        self.assertEqual(pth.value,"/mnt/data/test.txt")

class TestPathExistsParam(unittest.TestCase):
    """
    Tests for the 'PathExistsParam' class
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestPathExistsParam')

    def tearDown(self):
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    def test_pathexistsparam_with_strings(self):
        """
        PathExistsParam: handle static strings
        """
        dir_path = self.wd
        file_path = os.path.join(self.wd,"file.txt")
        with open(file_path,'wt') as fp:
            fp.write("test\n")
        missing_path = os.path.join(self.wd,"missing.txt")
        self.assertTrue(PathExistsParam(dir_path).value)
        self.assertTrue(PathExistsParam(file_path).value)
        self.assertFalse(PathExistsParam(missing_path).value)

    def test_pathexistsparam_with_pipelineparams(self):
        """
        PathExistsParam: handle PipelineParams
        """
        dir_path = PipelineParam(value=self.wd)
        file_path = PipelineParam(value=os.path.join(self.wd,"file.txt"))
        with open(file_path.value,'wt') as fp:
            fp.write("test\n")
        missing_path = PipelineParam(
            value=os.path.join(self.wd,"missing.txt"))
        self.assertTrue(PathExistsParam(dir_path).value)
        self.assertTrue(PathExistsParam(file_path).value)
        self.assertFalse(PathExistsParam(missing_path).value)

    def test_pathexistsparam_with_pathjoinparams(self):
        """
        PathExistsParam: handle PathJoinParams
        """
        dir_path = PathJoinParam(self.wd)
        file_path = PathJoinParam(self.wd,"file.txt")
        with open(file_path.value,'wt') as fp:
            fp.write("test\n")
        missing_path = PathJoinParam(self.wd,"missing.txt")
        self.assertTrue(PathExistsParam(dir_path).value)
        self.assertTrue(PathExistsParam(file_path).value)
        self.assertFalse(PathExistsParam(missing_path).value)

class TestFunctionParam(unittest.TestCase):
    """
    Tests for the 'FunctionParam' class
    """
    #@unittest.skip("Skipped")
    def test_functionparam_with_static_values(self):
        """
        FunctionParam: handle static values
        """
        # Lambda function
        func_param = FunctionParam(lambda x: x,"hello")
        self.assertEqual(func_param.value,"hello")
        func_param = FunctionParam(lambda x,y: x + y,1,2)
        self.assertEqual(func_param.value,3)
        # Function with keywords
        def func(x,y,z=None):
            return (x,y,z)
        func_param = FunctionParam(func,"hello","goodbye")
        self.assertEqual(func_param.value,("hello","goodbye",None))
        func_param = FunctionParam(func,"hello","goodbye",z="surprise!")
        self.assertEqual(func_param.value,("hello","goodbye","surprise!"))

    def test_functionparam_with_parameters(self):
        """
        FunctionParam: handle pipeline parameters
        """
        # Lambda function
        p = PipelineParam(value="hello")
        func_param = FunctionParam(lambda x: x,p)
        self.assertEqual(func_param.value,"hello")
        p.set("goodbye")
        self.assertEqual(func_param.value,"goodbye")
        px = PipelineParam(value=1)
        py = PipelineParam(value=2)
        func_param = FunctionParam(lambda x,y: x + y,px,py)
        self.assertEqual(func_param.value,3)
        px.set(3)
        py.set(4)
        self.assertEqual(func_param.value,7)
        # Function with keywords
        def func(x,y,z=None):
            return (x,y,z)
        px = PipelineParam(value="hello")
        py = PipelineParam(value="goodbye")
        pz = PipelineParam(value="surprise!")
        func_param = FunctionParam(func,px,py)
        self.assertEqual(func_param.value,("hello","goodbye",None))
        func_param = FunctionParam(func,px,py,z=pz)
        self.assertEqual(func_param.value,("hello","goodbye","surprise!"))
        px.set("goodbye")
        py.set("hello")
        pz.set("backwards")
        self.assertEqual(func_param.value,("goodbye","hello","backwards"))

    def test_functionparam_trap_exceptions_from_function(self):
        """
        FunctionParam: trap exceptions from function call
        """
        # Trap type error
        exception = False
        try:
            FunctionParam(lambda x: int(x),"non integer").value
        except PipelineError:
            exception = True
        except Exception:
            pass
        self.assertTrue(exception,"Should have raised exception")
        # Trap attribute error
        exception = False
        try:
            FunctionParam(lambda x: x.missing,123).value
        except PipelineError:
            exception = True
        except Exception:
            pass
        self.assertTrue(exception,"Should have raised exception")

class TestDispatcher(unittest.TestCase):

    def setUp(self):
        # Make a temporary working dir
        self.working_dir = tempfile.mkdtemp(
            suffix='TestDispatcher')

    def tearDown(self):
        # Remove temp dir
        if os.path.exists(self.working_dir):
            shutil.rmtree(self.working_dir)

    def test_dispatcher_subprocess(self):
        """
        Dispatcher: runs a function via subprocess module
        """
        dispatcher_wd = os.path.join(self.working_dir,
                                     "test_dispatcher")
        d = Dispatcher(working_dir=dispatcher_wd)
        def hello(name):
            return "Hello %s!" % name
        cmd = d.dispatch_function_cmd(hello,"World")
        exit_code = cmd.run_subprocess(log=os.path.join(self.working_dir,
                                                        "dispatcher.log"))
        self.assertEqual(exit_code,0)
        result = d.get_result()
        self.assertEqual(result,"Hello World!")

    def test_dispatcher_simplejobrunner(self):
        """
        Dispatcher: runs a function via SimpleJobRunner
        """
        dispatcher_wd = os.path.join(self.working_dir,
                                     "test_dispatcher")
        d = Dispatcher(working_dir=dispatcher_wd)
        def hello(name):
            return "Hello %s!" % name
        cmd = d.dispatch_function_cmd(hello,"World")
        runner = SimpleJobRunner(log_dir=self.working_dir)
        job_id = runner.run("Hello world",
                            self.working_dir,
                            cmd.command,
                            cmd.args)
        while runner.isRunning(job_id):
            time.sleep(0.1)
        exit_code = runner.exit_status(job_id)
        self.assertEqual(exit_code,0)
        result = d.get_result()
        self.assertEqual(result,"Hello World!")

class TestResolveParameter(unittest.TestCase):

    def test_resolve_parameter_pipelineparam(self):
        """
        resolve_parameter: returns value for PipelineParam
        """
        self.assertEqual(
            resolve_parameter(PipelineParam(value="this is the value")),
            "this is the value")

    def test_resolve_parameter_non_pipelineparam(self):
        """
        resolve_parameter: returns original object if not PipelineParam
        """
        self.assertEqual(resolve_parameter("this is the value"),
                         "this is the value")
