#######################################################################
# Tests for simple_scheduler.py module
#######################################################################
import unittest
import os
import time
import logging
from JobRunner import BaseJobRunner
from simple_scheduler import *

class MockJobRunner(BaseJobRunner):
    """Mock job runner implementation of BaseJobRunner

    """
    def __init__(self):
        self.__jobcount = 0
        self.__jobs = dict()

    def run(self,name,working_dir,script,args):
        self.__jobcount += 1
        job_id = self.__jobcount
        self.__jobs[job_id] = { 'name': name,
                                'working_dir': working_dir,
                                'script': script,
                                'args': args }
        return job_id

    def logFile(self,job_id):
        return "%s.%s.log" % (self.__jobs[job_id]['name'],job_id)

    def errFile(self,job_id):
        return self.logFile(job_id)

    def terminate(self,job_id):
        if job_id in self.__jobs:
            del(self.__jobs[job_id])

    def list(self):
        return self.__jobs.keys()

class CallbackTester:
    """Utility class for testing callbacks from scheduler

    To use: create a CallbackTester instance e.g.:

    >>> cb = CallbackTester()

    Then pass the 'call_me' method as a callback function when
    scheduling a job or adding a group e.g.

    >>> s.submit(['sleep','50'],callbacks=(cb.call_me,))

    Check if the callback has been invoked by checking the
    'invoked' attribute; other attributes ('name', 'jobs'
    and 'sched') will contain the values passed back when
    the callback was invoked.

    If the CallbackTester is created with 'raise_exception'
    set True then the 'call_me' method will also raise an
    Exception when it's invoked.

    """
    def __init__(self,raise_exception=False):
        self.raise_exception = raise_exception
        self.invoked = False
        self.name = None
        self.jobs = None
        self.sched = None

    def call_me(self,name,jobs,sched):
        self.invoked = True
        self.name = name
        self.jobs = jobs
        self.sched = sched
        if self.raise_exception:
            raise Exception,"call_me raised expected exception"

class TestSimpleScheduler(unittest.TestCase):
    """Unit tests for GetFastqFiles function

    """
    def test_simple_scheduler(self):
        """Start and stop without running any jobs

        """
        sched = SimpleScheduler(runner=MockJobRunner())
        sched.start()
        self.assertEqual(sched.n_waiting,0)
        self.assertEqual(sched.n_running,0)
        self.assertEqual(sched.n_finished,0)
        self.assertTrue(sched.is_empty())
        sched.stop()

    def test_simple_scheduler_run_single_job(self):
        """Run a single job

        """
        sched = SimpleScheduler(runner=MockJobRunner(),poll_interval=0.01)
        sched.start()
        job = sched.submit(['sleep','50'])
        # Wait for scheduler to catch up
        time.sleep(0.1)
        self.assertEqual(sched.n_waiting,0)
        self.assertEqual(sched.n_running,1)
        self.assertEqual(sched.n_finished,0)
        self.assertFalse(sched.is_empty())
        # Finish job, wait for scheduler to catch up
        job.terminate()
        time.sleep(0.1)
        self.assertEqual(sched.n_waiting,0)
        self.assertEqual(sched.n_running,0)
        self.assertEqual(sched.n_finished,1)
        self.assertTrue(sched.is_empty())
        sched.stop()

    def test_simple_scheduler_run_multiple_jobs(self):
        """Run several jobs

        """
        sched = SimpleScheduler(runner=MockJobRunner(),poll_interval=0.01)
        sched.start()
        job_1 = sched.submit(['sleep','10'])
        job_2 = sched.submit(['sleep','20'])
        job_3 = sched.submit(['sleep','30'])
        # Wait for scheduler to catch up
        time.sleep(0.1)
        self.assertEqual(sched.n_waiting,0)
        self.assertEqual(sched.n_running,3)
        self.assertFalse(sched.is_empty())
        # Finish a job, wait for scheduler to catch up
        job_1.terminate()
        time.sleep(0.1)
        self.assertEqual(sched.n_waiting,0)
        self.assertEqual(sched.n_running,2)
        self.assertFalse(sched.is_empty())
        # Finish remaining jobs, wait for scheduler to catch up
        job_2.terminate()
        job_3.terminate()
        time.sleep(0.1)
        self.assertEqual(sched.n_waiting,0)
        self.assertEqual(sched.n_running,0)
        self.assertTrue(sched.is_empty())
        sched.stop()

    def test_simple_scheduler_run_multiple_jobs_with_limit(self):
        """Run several jobs with limit on maximum concurrent jobs

        """
        sched = SimpleScheduler(runner=MockJobRunner(),poll_interval=0.01,
                                max_concurrent=2)
        sched.start()
        job_1 = sched.submit(['sleep','10'])
        job_2 = sched.submit(['sleep','20'])
        job_3 = sched.submit(['sleep','30'])
        # Wait for scheduler to catch up
        time.sleep(0.1)
        self.assertEqual(sched.n_waiting,1)
        self.assertEqual(sched.n_running,2)
        self.assertEqual(sched.n_finished,0)
        self.assertFalse(sched.is_empty())
        # Finish a job, wait for scheduler to catch up
        job_1.terminate()
        time.sleep(0.1)
        self.assertEqual(sched.n_waiting,0)
        self.assertEqual(sched.n_running,2)
        self.assertEqual(sched.n_finished,1)
        self.assertFalse(sched.is_empty())
        # Finish remaining jobs, wait for scheduler to catch up
        job_2.terminate()
        job_3.terminate()
        time.sleep(0.1)
        self.assertEqual(sched.n_waiting,0)
        self.assertEqual(sched.n_running,0)
        self.assertEqual(sched.n_finished,3)
        self.assertTrue(sched.is_empty())
        sched.stop()

    def test_simple_scheduler_run_dependent_jobs(self):
        """Run several jobs with one dependent on another

        """
        sched = SimpleScheduler(runner=MockJobRunner(),poll_interval=0.01)
        sched.start()
        job_1 = sched.submit(['sleep','10'],name="sleep_10")
        job_2 = sched.submit(['sleep','20'],name="sleep_20",wait_for=('sleep_10',))
        job_3 = sched.submit(['sleep','30'],name="sleep_30")
        # Wait for scheduler to catch up
        time.sleep(0.1)
        self.assertEqual(sched.n_waiting,1)
        self.assertEqual(sched.n_running,2)
        self.assertEqual(sched.n_finished,0)
        self.assertFalse(sched.is_empty())
        # Finish a job, wait for scheduler to catch up
        job_3.terminate()
        time.sleep(0.1)
        self.assertEqual(sched.n_waiting,1)
        self.assertEqual(sched.n_running,1)
        self.assertEqual(sched.n_finished,1)
        self.assertFalse(sched.is_empty())
        # Finish a job, wait for scheduler to catch up
        job_1.terminate()
        time.sleep(0.1)
        self.assertEqual(sched.n_waiting,0)
        self.assertEqual(sched.n_running,1)
        self.assertEqual(sched.n_finished,2)
        self.assertFalse(sched.is_empty())
        # Finish remaining job, wait for scheduler to catch up
        job_2.terminate()
        time.sleep(0.1)
        self.assertEqual(sched.n_waiting,0)
        self.assertEqual(sched.n_running,0)
        self.assertEqual(sched.n_finished,3)
        self.assertTrue(sched.is_empty())
        sched.stop()

    def test_simple_scheduler_with_group(self):
        """Run group of jobs

        """
        sched = SimpleScheduler(runner=MockJobRunner(),poll_interval=0.01)
        sched.start()
        # Add a group
        group = sched.group("grp_1")
        # Wait for scheduler to catch up
        time.sleep(0.1)
        self.assertEqual(group.group_name,"grp_1")
        self.assertFalse(group.closed)
        self.assertFalse(group.is_running)
        self.assertEqual(len(group.jobs),0)
        self.assertEqual(sched.n_waiting,0)
        self.assertEqual(sched.n_running,0)
        self.assertEqual(sched.n_finished,0)
        # Add some jobs
        job_1 = group.add(['sleep','10'],name="sleep_10")
        job_2 = group.add(['sleep','20'],name="sleep_20")
        time.sleep(0.1)
        self.assertFalse(group.closed)
        self.assertFalse(group.is_running)
        self.assertEqual(len(group.jobs),2)
        self.assertEqual(sched.n_waiting,0)
        self.assertEqual(sched.n_running,2)
        self.assertEqual(sched.n_finished,0)
        # Add another job
        job_3 = group.add(['sleep','30'],name="sleep_30")
        time.sleep(0.1)
        self.assertFalse(group.closed)
        self.assertFalse(group.is_running)
        self.assertEqual(len(group.jobs),3)
        self.assertEqual(sched.n_waiting,0)
        self.assertEqual(sched.n_running,3)
        self.assertEqual(sched.n_finished,0)
        # Close the group
        group.close()
        time.sleep(0.1)
        self.assertTrue(group.closed)
        self.assertTrue(group.is_running)
        self.assertEqual(len(group.jobs),3)
        self.assertEqual(sched.n_waiting,0)
        self.assertEqual(sched.n_running,3)
        self.assertEqual(sched.n_finished,0)
        # Try to add another job - should raise an exception
        self.assertRaises(Exception,group.add,['sleep','40'],name="sleep_40")
        time.sleep(0.1)
        self.assertTrue(group.closed)
        self.assertTrue(group.is_running)
        self.assertEqual(len(group.jobs),3)
        self.assertEqual(sched.n_waiting,0)
        self.assertEqual(sched.n_running,3)
        self.assertEqual(sched.n_finished,0)
        # Finish a job, wait for scheduler to catch up
        job_3.terminate()
        time.sleep(0.1)
        self.assertTrue(group.closed)
        self.assertTrue(group.is_running)
        self.assertEqual(sched.n_waiting,0)
        self.assertEqual(sched.n_running,2)
        self.assertEqual(sched.n_finished,1)
        # Finish remaining jobs, wait for scheduler to catch up
        job_1.terminate()
        job_2.terminate()
        time.sleep(0.1)
        self.assertTrue(group.closed)
        self.assertFalse(group.is_running)
        self.assertEqual(sched.n_waiting,0)
        self.assertEqual(sched.n_running,0)
        self.assertEqual(sched.n_finished,3)
        self.assertTrue(sched.is_empty())
        sched.stop()

    def test_callback_from_job(self):
        """Add callback function to a job

        """
        sched = SimpleScheduler(runner=MockJobRunner(),poll_interval=0.01)
        sched.start()
        # Add a job with a callback
        cb = CallbackTester()
        job = sched.submit(['sleep','50'],callbacks=(cb.call_me,))
        # Wait for scheduler to catch up
        time.sleep(0.1)
        self.assertFalse(cb.invoked)
        # Finish job, wait for scheduler to catch up
        job.terminate()
        time.sleep(0.1)
        self.assertTrue(cb.invoked)
        self.assertEqual(len(cb.jobs),1)
        self.assertEqual(cb.jobs[0],job)
        self.assertEqual(cb.sched,sched)
        self.assertTrue(sched.is_empty())
        sched.stop()

    def test_callback_from_group(self):
        """Add callback function to a group

        """
        sched = SimpleScheduler(runner=MockJobRunner(),poll_interval=0.01)
        sched.start()
        # Add a group with a callback
        cb = CallbackTester()
        group = sched.group("grp_1",callbacks=(cb.call_me,))
        job1 = group.add(['sleep','10'])
        job2 = group.add(['sleep','20'])
        group.close()
        # Wait for scheduler to catch up
        time.sleep(0.1)
        self.assertFalse(cb.invoked)
        # Finish first job, wait for scheduler to catch up
        job1.terminate()
        time.sleep(0.1)
        self.assertFalse(cb.invoked)
        # Finish second job
        job2.terminate()
        time.sleep(0.1)
        self.assertTrue(cb.invoked)
        self.assertEqual(len(cb.jobs),1)
        self.assertEqual(cb.jobs[0],group)
        self.assertEqual(cb.sched,sched)
        self.assertTrue(sched.is_empty())
        sched.stop()

    def test_callback_from_job_raises_exception(self):
        """Check 'bad' callback doesn't crash the scheduler

        """
        sched = SimpleScheduler(runner=MockJobRunner(),poll_interval=0.01)
        sched.start()
        # Add a job with a 'bad' callback
        cb = CallbackTester(raise_exception=True)
        job = sched.submit(['sleep','50'],callbacks=(cb.call_me,))
        # Wait for scheduler to catch up
        time.sleep(0.1)
        self.assertFalse(cb.invoked)
        # Finish job, wait for scheduler to catch up
        job.terminate()
        time.sleep(0.1)
        self.assertTrue(cb.invoked)
        self.assertEqual(len(cb.jobs),1)
        self.assertEqual(cb.jobs[0],job)
        self.assertEqual(cb.sched,sched)
        # Check we can still submit a job i.e. scheduler
        # is still running
        self.assertEqual(sched.n_finished,1)
        job = sched.submit(['sleep','50'])
        time.sleep(0.1)
        job.terminate()
        time.sleep(0.1)
        self.assertEqual(sched.n_finished,2)
        self.assertTrue(sched.is_empty())
        sched.stop()

#######################################################################
# Main program
#######################################################################

if __name__ == "__main__":
    # Turn off most logging output for tests
    logging.getLogger().setLevel(logging.CRITICAL)
    # Run tests
    unittest.main()
