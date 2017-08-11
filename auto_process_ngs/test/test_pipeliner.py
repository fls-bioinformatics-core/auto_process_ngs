#######################################################################
# Tests for pipeliner.py module
#######################################################################

import unittest
import tempfile
import shutil
import time
import os
from auto_process_ngs.simple_scheduler import SimpleScheduler
from auto_process_ngs.applications import Command
from auto_process_ngs.pipeliner import Pipeline
from auto_process_ngs.pipeliner import PipelineTask
from auto_process_ngs.pipeliner import PipelineCommand
from auto_process_ngs.pipeliner import PipelineCommandWrapper
from auto_process_ngs.pipeliner import FileCollector

# Unit tests

class TestPipeline(unittest.TestCase):

    def setUp(self):
        # Set up a scheduler
        self.sched = SimpleScheduler()
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

    def test_simple_pipeline(self):
        """
        Pipeline: define and run a simple pipeline
        """
        # Define a reusable task
        # Appends item to a list
        class Append(PipelineTask):
            def init(self,l,s):
                self.l = list()
            def setup(self):
                for item in self.args.l:
                    self.l.append(item)
                self.l.append(self.args.s)
            def output(self):
                return self.l
        # Build the pipeline
        ppl = Pipeline()
        task1 = Append("Append 1",(),"item1")
        task2 = Append("Append 2",task1.output(),"item2")
        ppl.add_task(task2,requires=(task1,))
        # Run the pipeline
        exit_status = ppl.run(sched=self.sched,
                              working_dir=self.working_dir)
        # Check the outputs
        self.assertEqual(exit_status,0)
        self.assertEqual(task1.output(),["item1"])
        self.assertEqual(task2.output(),["item1","item2"])

    def test_pipeline_with_commands(self):
        """
        Pipeline: define and run pipeline with commands
        """
        # Define a task
        # Echoes/appends text to a file
        class Echo(PipelineTask):
            def init(self,f,s):
                pass
            def setup(self):
                self.add_cmd(
                    PipelineCommandWrapper(
                        "Echo text to file",
                        "echo",self.args.s,
                        ">>",self.args.f))
            def output(self):
                return self.args.f
        # Build the pipeline
        ppl = Pipeline()
        task1 = Echo("Write item1","out.txt","item1")
        task2 = Echo("Write item2",task1.output(),"item2")
        ppl.add_task(task2,requires=(task1,))
        # Run the pipeline
        exit_status = ppl.run(sched=self.sched,
                              working_dir=self.working_dir)
        # Check the outputs
        self.assertEqual(exit_status,0)
        out_file = os.path.join(self.working_dir,"out.txt")
        self.assertTrue(os.path.exists(out_file))
        self.assertEqual(open(out_file,'r').read(),
                         "item1\nitem2\n")

    def test_pipeline_working_dir_is_respected(self):
        """
        Pipeline: check pipeline respects the working directory
        """
        # Define tasks
        # Echoes/appends text to a file via shell command
        class Echo(PipelineTask):
            def init(self,f,s):
                pass
            def setup(self):
                self.add_cmd(
                    PipelineCommandWrapper(
                        "Echo text to file",
                        "echo",self.args.s,
                        ">>",self.args.f))
            def output(self):
                return self.args.f
        # Writes text to a file via Python
        class Print(PipelineTask):
            def init(self,f,s):
                pass
            def setup(self):
                with open(self.args.f,'a') as fp:
                    fp.write("%s\n" % self.args.s)
            def output(self):
                return self.args.f
        # Build the pipeline
        ppl = Pipeline()
        task1 = Echo("Echo item1","out.txt","item1")
        task2 = Print("Print item2",task1.output(),"item2")
        ppl.add_task(task2,requires=(task1,))
        # Run the pipeline
        exit_status = ppl.run(sched=self.sched,
                              working_dir=self.working_dir)
        # Check the outputs
        self.assertEqual(exit_status,0)
        out_file = os.path.join(self.working_dir,"out.txt")
        self.assertTrue(os.path.exists(out_file))
        self.assertEqual(open(out_file,'r').read(),
                         "item1\nitem2\n")

    def test_pipeline_stops_on_task_failure(self):
        """
        Pipeline: check pipeline stops on task failure
        """
        # Define a reusable task
        # Appends item to a list
        class Append(PipelineTask):
            def init(self,l,s):
                self.l = list()
            def setup(self):
                for item in self.args.l:
                    self.l.append(item)
                self.l.append(self.args.s)
            def output(self):
                return self.l
        # Define a task that always fails
        class Failure(PipelineTask):
            def init(self):
                pass
            def setup(self):
                self.fail(message="Automatic fail")
            def output(self):
                return None
        # Build the pipeline
        ppl = Pipeline()
        task1 = Append("Append 1",(),"item1")
        task2 = Failure("Failing task")
        task3 = Append("Append 3",task1.output(),"item3")
        ppl.add_task(task2,requires=(task1,))
        ppl.add_task(task3,requires=(task2,))
        # Run the pipeline
        exit_status = ppl.run(sched=self.sched,
                              working_dir=self.working_dir)
        # Check the outputs
        self.assertEqual(exit_status,1)
        self.assertEqual(task1.output(),["item1"])
        self.assertEqual(task2.exit_code,1)
        self.assertEqual(task3.output(),[])

class TestPipelineTask(unittest.TestCase):

    def setUp(self):
        # Set up a scheduler
        self.sched = SimpleScheduler()
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

    def test_pipelinetask_invocations(self):
        """
        PipelineTask: check task methods are invoked
        """
        # Define a simplistic task which does nothing
        class CheckInvocations(PipelineTask):
            def init(self):
                self.invocations = list()
                self.invocations.append("init")
            def setup(self):
                self.invocations.append("setup")
            def finish(self):
                self.invocations.append("finish")
            def output(self):
                return self.invocations
        # Make a task instance
        task = CheckInvocations("Check method invocations")
        # Check initial state
        self.assertFalse(task.completed)
        self.assertEqual(task.exit_code,None)
        self.assertEqual(task.output(),["init"])
        # Run the task
        task.run(sched=self.sched,
                 working_dir=self.working_dir,
                 async=False)
        # Check final state
        self.assertTrue(task.completed)
        self.assertEqual(task.exit_code,0)
        self.assertEqual(task.output(),["init","setup","finish"])

    def test_pipelinetask_no_commands(self):
        """
        PipelineTask: run task with no commands
        """
        # Define a task with no commands
        class Add(PipelineTask):
            def init(self,x,y):
                self.result = list()
            def setup(self):
                self.result.append(self.args.x+self.args.y)
            def output(self):
                return self.result
        # Make a task instance
        task = Add("Add two numbers",1,2)
        # Check initial state
        self.assertEqual(task.args.x,1)
        self.assertEqual(task.args.y,2)
        self.assertFalse(task.completed)
        self.assertEqual(task.exit_code,None)
        self.assertEqual(task.output(),[])
        # Run the task
        task.run(sched=self.sched,
                 working_dir=self.working_dir,
                 async=False)
        # Check final state
        self.assertTrue(task.completed)
        self.assertEqual(task.exit_code,0)
        self.assertEqual(task.output(),[3])
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
            def output(self):
                return None
        # Make a task instance
        task = Echo("Echo string","Hello!")
        # Check initial state
        self.assertEqual(task.args.s,"Hello!")
        self.assertFalse(task.completed)
        self.assertEqual(task.exit_code,None)
        self.assertEqual(task.output(),None)
        # Run the task
        task.run(sched=self.sched,
                 working_dir=self.working_dir,
                 async=False)
        # Check final state
        self.assertTrue(task.completed)
        self.assertEqual(task.exit_code,0)
        self.assertEqual(task.output(),None)
        self.assertEqual(task.stdout,"Hello!\n")

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
                for i in xrange(self.args.n):
                    self.add_cmd(
                        PipelineCommandWrapper(
                            "Echo text","echo",self.args.s))
            def output(self):
                return None
        # Make a task instance
        task = MultipleEcho("Echo string 3 times","Hello!",3)
        # Check initial state
        self.assertFalse(task.completed)
        self.assertEqual(task.exit_code,None)
        self.assertEqual(task.output(),None)
        # Run the task
        task.run(sched=self.sched,
                 working_dir=self.working_dir,
                 async=False)
        # Check final state
        self.assertTrue(task.completed)
        self.assertEqual(task.exit_code,0)
        self.assertEqual(task.output(),None)
        self.assertEqual(task.stdout,"Hello!\n\nHello!\n\nHello!\n")

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
                self.add_cmd(
                    PipelineCommandWrapper(
                        "Echo message","echo","should not execute"))
            def output(self):
                return None
        # Make a task instance
        task = FailingTask("This will fail")
        # Check initial state
        self.assertFalse(task.completed)
        self.assertEqual(task.exit_code,None)
        self.assertEqual(task.output(),None)
        # Run the task
        task.run(sched=self.sched,
                 working_dir=self.working_dir,
                 async=False)
        # Check final state
        self.assertTrue(task.completed)
        self.assertEqual(task.exit_code,123)
        self.assertEqual(task.output(),None)
        self.assertEqual(task.stdout,"")

class TestPipelineCommand(unittest.TestCase):

    def setUp(self):
        # Make a temporary working dir
        self.working_dir = tempfile.mkdtemp(
            suffix='TestPipelineCommand')

    def tearDown(self):
        # Remove temp dir
        if os.path.exists(self.working_dir):
            shutil.rmtree(self.working_dir)

    def test_pipelinecommand(self):
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
        self.assertEqual(open(script_file,'r').read(),
                         "#!/bin/bash\necho hello there")

class TestPipelineCommandWrapper(unittest.TestCase):

    def test_piplinecommandwrapper(self):
        # Make a pipeline command wrapper
        cmd = PipelineCommandWrapper("Echo text","echo","hello")
        # Check name and generated command
        self.assertEqual(cmd.name(),"echo_text")
        self.assertEqual(str(cmd.cmd()),"echo hello")
        # Add argument and check updated command
        cmd.add_args("there")
        self.assertEqual(str(cmd.cmd()),"echo hello there")

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

