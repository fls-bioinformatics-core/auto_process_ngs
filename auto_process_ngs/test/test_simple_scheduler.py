#######################################################################
# Tests for simple_scheduler.py module
#######################################################################
import unittest
import cStringIO
import os
import sys
import time
import logging
from bcftbx.JobRunner import BaseJobRunner
from auto_process_ngs.simple_scheduler import *

class MockJobRunner(BaseJobRunner):
    """Mock job runner implementation of BaseJobRunner

    """
    def __init__(self):
        self.__jobcount = 0
        self.__jobs = dict()
        self.__log_dirs = dict()
        BaseJobRunner.__init__(self)

    def run(self,name,working_dir,script,args):
        self.__jobcount += 1
        job_id = self.__jobcount
        self.__jobs[job_id] = { 'name': name,
                                'working_dir': working_dir,
                                'script': script,
                                'args': args }
        self.__log_dirs[job_id] = self.log_dir
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
    """Unit tests for SimpleScheduler class

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

    def test_simple_scheduler_run_single_job_set_log_dir(self):
        """Run a single job and explicitly set the log directory

        """
        sched = SimpleScheduler(runner=MockJobRunner(),poll_interval=0.01)
        sched.start()
        job = sched.submit(['sleep','50'],)
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
        self.assertFalse(group.completed)
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
        self.assertFalse(group.completed)
        self.assertEqual(len(group.jobs),2)
        self.assertEqual(sched.n_waiting,0)
        self.assertEqual(sched.n_running,2)
        self.assertEqual(sched.n_finished,0)
        # Add another job
        job_3 = group.add(['sleep','30'],name="sleep_30")
        time.sleep(0.1)
        self.assertFalse(group.closed)
        self.assertFalse(group.is_running)
        self.assertFalse(group.completed)
        self.assertEqual(len(group.jobs),3)
        self.assertEqual(sched.n_waiting,0)
        self.assertEqual(sched.n_running,3)
        self.assertEqual(sched.n_finished,0)
        # Close the group
        group.close()
        time.sleep(0.1)
        self.assertTrue(group.closed)
        self.assertTrue(group.is_running)
        self.assertFalse(group.completed)
        self.assertEqual(len(group.jobs),3)
        self.assertEqual(sched.n_waiting,0)
        self.assertEqual(sched.n_running,3)
        self.assertEqual(sched.n_finished,0)
        # Try to add another job - should raise an exception
        self.assertRaises(Exception,group.add,['sleep','40'],name="sleep_40")
        time.sleep(0.1)
        self.assertTrue(group.closed)
        self.assertTrue(group.is_running)
        self.assertFalse(group.completed)
        self.assertEqual(len(group.jobs),3)
        self.assertEqual(sched.n_waiting,0)
        self.assertEqual(sched.n_running,3)
        self.assertEqual(sched.n_finished,0)
        # Finish a job, wait for scheduler to catch up
        job_3.terminate()
        time.sleep(0.1)
        self.assertTrue(group.closed)
        self.assertTrue(group.is_running)
        self.assertFalse(group.completed)
        self.assertEqual(sched.n_waiting,0)
        self.assertEqual(sched.n_running,2)
        self.assertEqual(sched.n_finished,1)
        # Finish remaining jobs, wait for scheduler to catch up
        job_1.terminate()
        job_2.terminate()
        time.sleep(0.1)
        self.assertTrue(group.closed)
        self.assertFalse(group.is_running)
        self.assertTrue(group.completed)
        self.assertEqual(sched.n_waiting,0)
        self.assertEqual(sched.n_running,0)
        self.assertEqual(sched.n_finished,3)
        self.assertTrue(sched.is_empty())
        sched.stop()

    def test_wait_for_group_completion(self):
        """Check group completion triggers start of pending job

        """
        sched = SimpleScheduler(runner=MockJobRunner(),poll_interval=0.01)
        sched.start()
        # Add a group
        group = sched.group("grp_1")
        # Add jobs
        time.sleep(0.1)
        job_1 = group.add(['sleep','10'],name="sleep_10")
        job_2 = group.add(['sleep','20'],name="sleep_20")
        # Close the group
        time.sleep(0.1)
        group.close()
        # Add a job waiting on this group
        time.sleep(0.1)
        job_3 = sched.submit(['sleep','50'],wait_for=("grp_1",))
        # Check status
        time.sleep(0.1)
        self.assertTrue(group.is_running)
        self.assertFalse(group.completed)
        self.assertFalse(job_3.is_running)
        # Finish jobs in group, wait for group to complete
        job_1.terminate()
        job_2.terminate()
        time.sleep(0.1)
        self.assertFalse(group.is_running)
        self.assertTrue(group.completed)
        # Check that remaining job is now running
        time.sleep(0.1)
        self.assertTrue(job_3.is_running)
        # Finish
        job_3.terminate()
        time.sleep(0.1)
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

    def test_set_job_working_dir(self):
        """Explicitly specify working directory for a job

        """
        sched = SimpleScheduler(runner=MockJobRunner(),poll_interval=0.01)
        sched.start()
        job = sched.submit(['sleep','50'],wd='/tmp')
        time.sleep(0.1)
        self.assertEqual(job.working_dir,'/tmp')
        job.terminate()
        time.sleep(0.1)
        self.assertTrue(sched.is_empty())
        sched.stop()

    def test_set_job_working_dir_in_group(self):
        """Explicitly specify working directory for job in a group

        """
        sched = SimpleScheduler(runner=MockJobRunner(),poll_interval=0.01)
        sched.start()
        group = sched.group('grp1')
        job = group.add(['sleep','50'],wd='/tmp')
        group.close()
        time.sleep(0.1)
        self.assertEqual(job.working_dir,'/tmp')
        job.terminate()
        time.sleep(0.1)
        self.assertTrue(sched.is_empty())
        sched.stop()

    def test_set_job_log_dir(self):
        """Explicitly specify log directory for a job

        """
        sched = SimpleScheduler(runner=MockJobRunner(),poll_interval=0.01)
        sched.start()
        job = sched.submit(['sleep','50'],log_dir='/logs')
        time.sleep(0.1)
        self.assertEqual(job.log_dir,'/logs')
        job.terminate()
        time.sleep(0.1)
        self.assertTrue(sched.is_empty())
        sched.stop()

    def test_set_group_log_dir(self):
        """Explicitly specify log directory for a group

        """
        sched = SimpleScheduler(runner=MockJobRunner(),poll_interval=0.01)
        sched.start()
        group = sched.group('grp1',log_dir='/logs')
        job = group.add(['sleep','50'])
        group.close()
        time.sleep(0.1)
        self.assertEqual(job.log_dir,'/logs')
        job.terminate()
        time.sleep(0.1)
        self.assertTrue(sched.is_empty())
        sched.stop()

    def test_simple_scheduler_lookup(self):
        """Lookup groups and jobs by name

        """
        sched = SimpleScheduler(runner=MockJobRunner(),poll_interval=0.01)
        sched.start()
        group = sched.group('grp1')
        job_1 = group.add(['sleep','10'],name='sleep_10')
        group.close()
        job_2 = sched.submit(['sleep','20'],name='sleep_20')
        job_3 = sched.submit(['sleep','30'],name='sleep_30')
        # Wait for scheduler to catch up
        time.sleep(0.1)
        self.assertEqual(sched.n_waiting,0)
        self.assertEqual(sched.n_running,3)
        self.assertFalse(sched.is_empty())
        # Try lookups
        self.assertEqual(sched.lookup('grp1'),group)
        self.assertEqual(sched.lookup('sleep_10'),job_1)
        self.assertEqual(sched.lookup('sleep_20'),job_2)
        self.assertEqual(sched.lookup('sleep_30'),job_3)
        self.assertEqual(sched.lookup('sleep_40'),None)
        # Finish jobs, wait for scheduler to catch up
        job_1.terminate()
        job_2.terminate()
        job_3.terminate()
        time.sleep(0.1)
        sched.stop()

    def test_simple_scheduler_find(self):
        """Find groups and jobs by simple name matching

        """
        sched = SimpleScheduler(runner=MockJobRunner(),poll_interval=0.01)
        sched.start()
        group = sched.group('grp1')
        job_1 = group.add(['sleep','10'],name='sleep_10')
        group.close()
        job_2 = sched.submit(['sleep','20'],name='sleep_20')
        job_3 = sched.submit(['sleep','30'],name='sleep_30')
        # Wait for scheduler to catch up
        time.sleep(0.1)
        self.assertEqual(sched.n_waiting,0)
        self.assertEqual(sched.n_running,3)
        self.assertFalse(sched.is_empty())
        # Try lookups
        self.assertEqual(sched.find('grp1'),[group,])
        self.assertEqual(sched.find('sleep_10'),[job_1,])
        self.assertEqual(sched.find('sleep_20'),[job_2,])
        self.assertEqual(sched.find('sleep_30'),[job_3,])
        self.assertEqual(sched.find('sleep_.*'),[job_1,job_2,job_3])
        self.assertEqual(sched.find('.*'),[group,job_1,job_2,job_3])
        self.assertEqual(sched.find('sleep_40'),[])
        # Finish jobs, wait for scheduler to catch up
        job_1.terminate()
        job_2.terminate()
        job_3.terminate()
        time.sleep(0.1)
        sched.stop()

class TestSchedulerJob(unittest.TestCase):
    """Unit tests for SchedulerJob class

    """
    def test_scheduler_job_basic(self):
        """Basic test of SchedulerJob

        """
        job = SchedulerJob(MockJobRunner(),['sleep','50'])
        self.assertEqual(job.job_name,None)
        self.assertEqual(job.job_number,None)
        self.assertEqual(job.log_dir,None)
        self.assertFalse(job.is_running)
        self.assertFalse(job.completed)
        job.start()
        self.assertTrue(job.is_running)
        self.assertFalse(job.completed)
        job.terminate()
        self.assertFalse(job.is_running)
        self.assertTrue(job.completed)

    def test_scheduler_job_set_log_dir(self):
        """Set explicit log_dir for SchedulerJob

        """
        runner = MockJobRunner()
        self.assertEqual(runner.log_dir,None)
        job = SchedulerJob(runner,['sleep','50'],log_dir='/logs')
        self.assertEqual(job.job_name,None)
        self.assertEqual(job.job_number,None)
        self.assertEqual(job.log_dir,'/logs')
        self.assertFalse(job.is_running)
        self.assertFalse(job.completed)
        self.assertEqual(runner.log_dir,None)
        job.start()
        self.assertEqual(job.log_dir,'/logs')
        self.assertTrue(job.is_running)
        self.assertFalse(job.completed)
        self.assertEqual(runner.log_dir,None)
        job.terminate()
        self.assertFalse(job.is_running)
        self.assertTrue(job.completed)

class TestSchedulerReporter(unittest.TestCase):
    """Unit tests for SchedulerReporter class

    """
    def test_empty_scheduler_reporter(self):
        """'Empty' SchedulerReporter returns no output for scheduler status
        """
        fp = cStringIO.StringIO()
        reporter = SchedulerReporter(fp=fp)
        sched = SimpleScheduler()
        reporter.scheduler_status(sched)
        self.assertEqual('',fp.getvalue())

    def test_scheduler_reporter_scheduler_status(self):
        """SchedulerReporter returns correct output for scheduler status
        """
        fp = cStringIO.StringIO()
        reporter = SchedulerReporter(fp=fp,scheduler_status="%(n_running)s jobs")
        sched = SimpleScheduler()
        reporter.scheduler_status(sched)
        self.assertEqual('0 jobs\n',fp.getvalue())

    def test_scheduler_reporter_job_scheduled(self):
        """SchedulerReporter returns correct output when job is scheduled
        """
        fp = cStringIO.StringIO()
        reporter = SchedulerReporter(
            fp=fp,
            job_scheduled="Job scheduled: #%(job_number)d: \"%(job_name)s\""
        )
        job = SchedulerJob(MockJobRunner(),['sleep','50'],
                           job_number=2,name='test',wait_for=[])
        reporter.job_scheduled(job)
        self.assertEqual('Job scheduled: #2: "test"\n',fp.getvalue())

    def test_scheduler_reporter_job_start(self):
        """SchedulerReporter returns correct output when job is started
        """
        fp = cStringIO.StringIO()
        reporter = SchedulerReporter(
            fp=fp,
            job_start="Job started: #%(job_number)d (%(job_id)s): \"%(job_name)s\""
        )
        job = SchedulerJob(MockJobRunner(),['sleep','50'],
                           job_number=2,name='test',wait_for=[])
        job.start()
        reporter.job_start(job)
        self.assertEqual('Job started: #2 (1): "test"\n',fp.getvalue())

    def test_scheduler_reporter_job_end(self):
        """SchedulerReporter returns correct output when job ends
        """
        fp = cStringIO.StringIO()
        reporter = SchedulerReporter(
            fp=fp,
            job_end="Job completed: #%(job_number)d (%(job_id)s): \"%(job_name)s\""
        )
        job = SchedulerJob(MockJobRunner(),['sleep','50'],
                           job_number=2,name='test',wait_for=[])
        job.start()
        job.terminate()
        reporter.job_end(job)
        self.assertEqual('Job completed: #2 (1): "test"\n',fp.getvalue())

    def test_scheduler_reporter_group_added(self):
        """SchedulerReporter returns correct output when group is added
        """
        fp = cStringIO.StringIO()
        reporter = SchedulerReporter(
            fp=fp,
            group_added="Group has been added: #%(group_id)d: \"%(group_name)s\""
        )
        group = SchedulerGroup('test',1,SimpleScheduler())
        reporter.group_added(group)
        self.assertEqual('Group has been added: #1: "test"\n',fp.getvalue())

    def test_scheduler_reporter_group_end(self):
        """SchedulerReporter returns correct output when group finishes
        """
        fp = cStringIO.StringIO()
        reporter = SchedulerReporter(
            fp=fp,
            group_end="Group completed: #%(group_id)d: \"%(group_name)s\""
        )
        group = SchedulerGroup('test',1,SimpleScheduler())
        job = group.add(['sleep','50'],name='sleep',wait_for=[])
        group.close()
        job.terminate()
        reporter.group_end(group)
        self.assertEqual('Group completed: #1: "test"\n',fp.getvalue())
