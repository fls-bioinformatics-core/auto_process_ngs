#######################################################################
# Tests for pipeliner.py module
#######################################################################

import unittest
import tempfile
import shutil
import os
from auto_process_ngs.simple_scheduler import SimpleScheduler
from auto_process_ngs.pipeliner import Pipeline
from auto_process_ngs.pipeliner import PipelineTask
from auto_process_ngs.pipeliner import PipelineCommand
from auto_process_ngs.pipeliner import PipelineCommandWrapper

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

