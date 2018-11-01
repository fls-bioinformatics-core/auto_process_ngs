#######################################################################
# Tests for pipeliner.py module
#######################################################################

import unittest
import tempfile
import shutil
import time
import os
import getpass
import platform
from auto_process_ngs.simple_scheduler import SimpleScheduler
from auto_process_ngs.applications import Command
from auto_process_ngs.pipeliner import Pipeline
from auto_process_ngs.pipeliner import PipelineTask
from auto_process_ngs.pipeliner import PipelineFunctionTask
from auto_process_ngs.pipeliner import PipelineCommand
from auto_process_ngs.pipeliner import PipelineCommandWrapper
from auto_process_ngs.pipeliner import PipelineFailure
from auto_process_ngs.pipeliner import FileCollector
from auto_process_ngs.pipeliner import Dispatcher
from bcftbx.JobRunner import SimpleJobRunner

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
        self.sched = SimpleScheduler(poll_interval=0.5)
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
        exit_status = ppl.run(working_dir=self.working_dir)
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
        task1 = Echo("Write item1","out.txt","item1")
        task2 = Echo("Write item2",task1.output.file,"item2")
        ppl.add_task(task2,requires=(task1,))
        # Run the pipeline
        exit_status = ppl.run(working_dir=self.working_dir)
        # Check the outputs
        self.assertEqual(exit_status,0)
        out_file = os.path.join(self.working_dir,"out.txt")
        self.assertTrue(os.path.exists(out_file))
        self.assertEqual(open(out_file,'r').read(),
                         "item1\nitem2\n")

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
        task1 = Echo("Write item1","out.txt","item1")
        task2 = Echo("Write item2",task1.output.file,"item2")
        ppl.add_task(task2,requires=(task1,))
        # Get a scheduler
        self._get_scheduler()
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
                self.add_output('file',f)
            def setup(self):
                self.add_cmd(
                    PipelineCommandWrapper(
                        "Echo text to file",
                        "echo",self.args.s,
                        ">>",self.args.f))
        # Writes text to a file via Python
        class Print(PipelineTask):
            def init(self,f,s):
                self.add_output('file',f)
            def setup(self):
                with open(self.args.f,'a') as fp:
                    fp.write("%s\n" % self.args.s)
        # Build the pipeline
        ppl = Pipeline()
        task1 = Echo("Echo item1","out.txt","item1")
        task2 = Print("Print item2",task1.output.file,"item2")
        ppl.add_task(task2,requires=(task1,))
        # Run the pipeline
        exit_status = ppl.run(working_dir=self.working_dir)
        # Check the outputs
        self.assertEqual(exit_status,0)
        out_file = os.path.join(self.working_dir,"out.txt")
        self.assertTrue(os.path.exists(out_file))
        self.assertEqual(open(out_file,'r').read(),
                         "item1\nitem2\n")

    def test_pipeline_stops_on_task_failure(self):
        """
        Pipeline: stops on task failure (IMMEDIATE mode)
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
        Pipeline: defer task failure (DEFERRED mode)
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
        self.assertEqual(ppl1.get_task(task5.id())[1],[task4,])

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
        task5 = Append("Append 5",task1.output.list,"item5")
        task6 = Append("Append 6",task3.output.list,"item6")
        task7 = Append("Append 7",task3.output.list,"item7")
        ppl2.add_task(task6,requires=(task5,))
        ppl2.add_task(task7,requires=(task6,))
        self.assertEqual(len(ppl2.task_list()),3)
        # Merge second pipeline into the first
        ppl1.merge_pipeline(ppl2)
        self.assertEqual(len(ppl1.task_list()),7)

class TestPipelineTask(unittest.TestCase):

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
                 async=False)
        # Check final state
        self.assertTrue(task.completed)
        self.assertEqual(task.exit_code,0)
        self.assertEqual(task.output.invocations,
                         ["init","setup","finish"])

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
                 async=False)
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
                 async=False)
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
        # Hello!
        # #### END Thu Aug 17 08:38:14 BST 2017
        # #### EXIT_CODE 0
        stdout = task.stdout.split("\n")
        self.assertEqual(len(stdout),8) # 8 = 7 + trailing newline
        self.assertEqual(stdout[0],"#### COMMAND Echo text")
        self.assertEqual(stdout[1],"#### HOSTNAME %s" % self._hostname())
        self.assertEqual(stdout[2],"#### USER %s" % self._user())
        self.assertTrue(stdout[3].startswith("#### START "))
        self.assertEqual(stdout[4],"Hello!")
        self.assertTrue(stdout[5].startswith("#### END "))
        self.assertEqual(stdout[6],"#### EXIT_CODE 0")

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
                 async=False)
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
        # #### END Thu Aug 17 08:38:14 BST 2017
        # #### EXIT_CODE 127
        stdout = task.stdout.split("\n")
        self.assertEqual(len(stdout),7) # 7 = 6 + trailing newline
        self.assertEqual(stdout[0],"#### COMMAND Nonexistant")
        self.assertEqual(stdout[1],"#### HOSTNAME %s" % self._hostname())
        self.assertEqual(stdout[2],"#### USER %s" % self._user())
        self.assertTrue(stdout[3].startswith("#### START "))
        self.assertTrue(stdout[4].startswith("#### END "))
        self.assertEqual(stdout[5],"#### EXIT_CODE 127")

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
        # Make a task instance
        task = MultipleEcho("Echo string 3 times","Hello!",3)
        # Check initial state
        self.assertFalse(task.completed)
        self.assertEqual(task.exit_code,None)
        self.assertFalse(task.output)
        # Run the task
        task.run(sched=self.sched,
                 working_dir=self.working_dir,
                 async=False)
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
        # Hello!
        # #### END Thu Aug 17 08:38:14 BST 2017
        # #### EXIT_CODE 0
        # ...x three times
        print task.stdout
        stdout = task.stdout.split("\n")
        self.assertEqual(len(stdout),22) # 22 = 21 + trailing newline
        for i in xrange(3):
            self.assertEqual(stdout[0+i*7],"#### COMMAND Echo text")
            self.assertEqual(stdout[1+i*7],"#### HOSTNAME %s" % self._hostname())
            self.assertEqual(stdout[2+i*7],"#### USER %s" % self._user())
            self.assertTrue(stdout[3+i*7].startswith("#### START "))
            self.assertEqual(stdout[4+i*7],"Hello!")
            self.assertTrue(stdout[5+i*7].startswith("#### END "))
            self.assertEqual(stdout[6+i*7],"#### EXIT_CODE 0")

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
        # Make a task instance
        task = FailingTask("This will fail")
        # Check initial state
        self.assertFalse(task.completed)
        self.assertEqual(task.exit_code,None)
        self.assertFalse(task.output)
        # Run the task
        task.run(sched=self.sched,
                 working_dir=self.working_dir,
                 async=False)
        # Check final state
        self.assertTrue(task.completed)
        self.assertEqual(task.exit_code,123)
        self.assertFalse(task.output)
        self.assertEqual(task.stdout,"")

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
                 async=False)
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
                 async=False)
        # Check final state
        self.assertTrue(task.completed)
        self.assertNotEqual(task.exit_code,0)
        self.assertEqual(task.result(),None)
        self.assertFalse(task.output)

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
                         "#!/bin/bash\n"
                         "echo \"#### COMMAND EchoCmd\"\n"
                         "echo \"#### HOSTNAME $HOSTNAME\"\n"
                         "echo \"#### USER $USER\"\n"
                         "echo \"#### START $(date)\"\n"
                         "echo 'hello there'\n"
                         "exit_code=$?\n"
                         "echo \"#### END $(date)\"\n"
                         "echo \"#### EXIT_CODE $exit_code\"\n"
                         "exit $exit_code")

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
