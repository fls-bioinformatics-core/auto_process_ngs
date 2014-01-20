#!/bin/env python
#
#     simple_scheduler.py: provide basic scheduler capability
#     Copyright (C) University of Manchester 2013 Peter Briggs
#
########################################################################
#
# simple_scheduler.py
#
#########################################################################


"""Python module to provide scheduler capability for running external
programs

"""

#######################################################################
# Module metadata
#######################################################################

__version__ = "0.0.3"

#######################################################################
# Import modules that this module depends on
#######################################################################

import JobRunner
from Pipeline import Job
import time
import os
import threading
import Queue
import logging

#######################################################################
# Classes
#######################################################################

class SimpleScheduler(threading.Thread):
    """Lightweight scheduler class

    Class providing simple job control functionality. Requests to
    run external programs are submitted to the scheduler via the
    'submit' method; the scheduler handles actually executing them.

    Submitted jobs can have dependencies on earlier jobs, in which
    case they will not be executed until those dependencies have
    completed.

    The scheduler runs in its own thread.

    Usage:
    
    >>> s = SimpleScheduler()
    >>> s.start()
    >>> s.submit(['fastq_screen',...],name='fastq_screen.model')
    >>> s.submit(['fastq_screen',...],name='fastq_screen.other')
    >>> s.submit(['fastq_screen',...],name='fastq_screen.rRNA')
    >>> s.submit(['fastqc',...],wait_for=('fastq_screen.model',...))
    
    """

    def __init__(self,runner=JobRunner.SimpleJobRunner(),max_concurrent=None,
                 poll_interval=5):
        """Create a new SimpleScheduler instance

        Arguments:
          runner: optional, default job runner to use
          max_concurrent: optional, maximum number of concurrent
            processes the scheduler will run (default: no limit)
          poll_interval: optional, number of seconds to wait in
            between checking for completed jobs etc in the scheduler
            loop (default: 5 seconds)

        """

        threading.Thread.__init__(self)
        self.setDaemon(1)
        # Default job runner
        self.__runner = runner
        # Maximum number of concurrent jobs
        self.__max_concurrent = max_concurrent
        # Length of time to wait between checking jobs
        self.__poll_interval = poll_interval
        # Internal job id counter
        self.__job_count = 0
        # Queue to add jobs
        self.__submitted = Queue.Queue()
        # List of scheduled (i.e.waiting) jobs
        self.__scheduled = []
        # List of running jobs
        self.__running = []
        # Handle names
        self.__names = []
        self.__finished_names = []
        # Handle groups
        self.__groups = []
        # Flag controlling whether scheduler is active
        self.__active = False

    def stop(self):
        """Stop the scheduler

        """
        self.__active = False

    @property
    def n_waiting(self):
        """Return number of jobs waiting to run

        """
        return len(self.__scheduled)

    @property
    def n_running(self):
        """Return number of jobs currently running

        """
        return len(self.__running)

    @property
    def job_number(self):
        """Internal: increment and return job count

        """
        self.__job_count += 1
        return self.__job_count

    @property
    def default_runner(self):
        """Return the default runner

        """
        return self.__runner

    def has_name(self,name):
        """Test if name has already been used

        Returns True if the name already exists, False
        otherwise.

        """
        return name in self.__names

    def is_empty(self):
        """Test if the scheduler has any jobs remaining

        Returns False if there are jobs running and/or waiting,
        True otherwise.

        """
        return not (self.n_waiting or self.n_running)

    def submit(self,args,runner=None,name=None,wait_for=[]):
        """Submit a request to run a job
        
        Arguments:
          args
          runner
          name
          wait_for

        Returns:
          SchedulerJob instance for the submitted job.

        """
        # Use a queue rather than modifying the waiting list
        # directly to try and avoid
        #
        # Check names are not duplicated
        if name is not None:
            if self.has_name(name):
                raise Exception,"Name '%s' already assigned" % name
            self.__names.append(name)
        # Check we're not waiting on a non-existent name
        for job_name in wait_for:
            if not self.has_name(job_name):
                raise Exception,"Job depends on a non-existent name '%s'" % job_name
        # Use default runner if none explicitly specified
        if runner is None:
            runner = self.default_runner
        # Schedule the job
        job = SchedulerJob(runner,args,job_number=self.job_number,
                           name=name,wait_for=wait_for)
        self.__submitted.put(job)
        logging.debug("Scheduled job #%d: \"%s\"" % (job.job_number,job))
        return job

    def group(self,name,wait_for=[]):
        """Create a group of jobs
        
        Arguments:
          name
          wait_for

        Returns:
          Empty SchedulerGroup instance.

        """
        # Check names are not duplicated
        if name is not None:
            if self.has_name(name):
                raise Exception,"Name '%s' already assigned" % name
            self.__names.append(name)
        new_group = SchedulerGroup(name,self.job_number,self,wait_for=wait_for)
        self.__groups.append(new_group)
        return new_group

    def run(self):
        """Internal: run method overriding that from base Thread class

        This method implements the scheduler loop.

        Don't call this directly - it is invoked by the Thread's 'start'
        method.

        """
        logging.debug("Starting simple scheduler")
        self.__active = True
        while self.__active:
            # Check for completed jobs
            updated_running_list = []
            for job in self.__running:
                if job.is_running:
                    logging.debug("Job #%s (id %s) still running \"%s\"" % (job.job_number,
                                                                            job.job_id,
                                                                            job))
                    updated_running_list.append(job)
                else:
                    logging.debug("Job #%s (id %s) completed \"%s\"" % (job.job_number,
                                                                        job.job_id,
                                                                        job))
                    if job.job_name is not None:
                        self.__finished_names.append(job.job_name)
            # Update the list of running jobs
            self.__running = updated_running_list
            # Check for completed groups
            updated_group_list = []
            for group in self.__groups:
                if group.submitted:
                    if group.is_running:
                        logging.debug("Group #%s (id %s) still running" % (group.group_name,
                                                                           group.group_id))
                        updated_group_list.append(group)
                    else:
                        logging.debug("Group #%s (id %s) completed" % (group.group_name,
                                                                       group.group_id))
                        self.__finished_names.append(group.group_name)
            # Update the list of groups
            self.__groups = updated_group_list
            # Add submitted jobs to the waiting list
            while not self.__submitted.empty():
                job = self.__submitted.get()
                self.__scheduled.append(job)
                logging.debug("Added job #%d (%s): \"%s\"" % (job.job_number,job.name,job))
            remaining_jobs = []
            for job in self.__scheduled:
                ok_to_run = True
                if self.__max_concurrent is not None and \
                   self.n_running == self.__max_concurrent:
                    # Scheduler capacity maxed out
                    ok_to_run = False
                else:
                    # Check if job is waiting for another job to finish
                    for name in job.waiting_for:
                        ok_to_run = (ok_to_run and name in self.__finished_names)
                if ok_to_run:
                    # Start the job running
                    job.start()
                    self.__running.append(job)
                    logging.debug("Started job #%s (id %s)" % (job.job_number,job.job_id))
                else:
                    # Hold back for now
                    remaining_jobs.append(job)
            # Update the scheduled job list
            self.__scheduled = remaining_jobs
            # Wait before going round again
            time.sleep(self.__poll_interval)

class SchedulerGroup:
    """Class providing an interface to schedule a group of jobs

    SchedulerGroup instances should normally be returned by a
    call to the 'group' method of a SimpleScheduler object.
    The group should be populated by calling its 'add' method
    (note that jobs are passed directly to the scheduler).

    Once all jobs have been added the 'submit' method should
    be invoked to indicate to the scheduler that the group is
    complete. At this point no more jobs can be added to the
    group, and the scheduler will check for when the group
    has finished running.
    
    """

    def __init__(self,name,group_id,parent_scheduler,wait_for=[]):
        """Create a new SchedulerGroup instance

        """
        self.group_name = name
        self.group_id = group_id
        self.waiting_for = list(wait_for)
        self.__scheduler = parent_scheduler
        self.__submitted = False
        self.__jobs = []

    @property
    def submitted(self):
        """Test whether group has been submitted

        """
        return self.__submitted

    @property
    def is_running(self):
        """Test whether group is running

        Returns True if group was submitted and if it has jobs
        that are still running, False if not.

        """
        if not self.submitted:
            return False
        for job in self.__jobs:
            if not job.submitted:
                return True
            elif job.is_running:
                return True
        return False

    def add(self,args,runner=None,name=None,wait_for=[]):
        """Add a request to run a job
        
        Arguments:
          args
          runner
          name
          wait_for

        Returns:
          SchedulerJob instance for the added job.

        """
        # Check we can still add jobs
        if self.submitted:
            raise Exception, \
                "Can't add job to group '%s': group already submitted" % self.group_name
        # Update list of jobs that this one needs to wait for 
        waiting_for = self.waiting_for + list(wait_for)
        # Submit the job to the scheduler and keep a reference
        logging.debug("Group '%s' #%s: adding job" % (self.group_name,self.group_id))
        job = self.__scheduler.submit(args,runner=runner,name=name,wait_for=waiting_for)
        self.__jobs.append(job)
        return job

    @property
    def jobs(self):
        """Return list of jobs

        """
        return self.__jobs

    def submit(self):
        """Indicate that all jobs have been added to the group

        """
        if self.submitted:
            raise Exception, "Group '%s' already submitted" % self.group_name
        self.__submitted = True

    def wait(self,poll_interval=5):
        """Wait for the group to complete

        Arguments:
          poll_interval: optional, number of seconds to wait in
            between checking if the group has completed (default: 5
            seconds)

        """
        if not self.submitted:
            raise Exception, "Group '%s' not submitted" % self.group_name
        logging.debug("Waiting for group '%s' (#%s)..." % (self.group_name,
                                                           self.group_id))
        while self.is_running:
            time.sleep(poll_interval)
        logging.debug("Group '%s' (#%s) finished" % (self.group_name,
                                                     self.group_id))

class SchedulerJob(Job):
    """Class providing an interface to scheduled jobs

    SchedulerJob instances should normally be returned by a
    call to the 'submit' method of a SimpleScheduler object.

    """

    def __init__(self,runner,args,job_number=None,name=None,wait_for=[]):
        """Create a new SchedulerJob instance

        """

        self.job_number = job_number
        self.job_name = name
        self.waiting_for = list(wait_for)
        if name is None:
            name = args[0]
        Job.__init__(self,runner,name,os.getcwd(),args[0],args[1:])

    @property
    def is_running(self):
        """Test if a job is running

        Returns True if the job has started and is still running,
        False otherwise.

        """
        # Check if job is running
        return self.isRunning()

    def wait(self,poll_interval=5):
        """Wait for the job to complete

        Arguments:
          poll_interval: optional, number of seconds to wait in
            between checking if the job has completed (default: 5
            seconds)

        """
        logging.debug("Waiting for job #%s..." % self.job_number)
        while self.job_id is None or self.is_running:
            time.sleep(poll_interval)
        logging.debug("Job #%s finished" % self.job_number)

    def __repr__(self):
        """Return string representation of the job command line

        """
        return ' '.join([str(x) for x in ([self.script,] + self.args)])

#######################################################################
# Main program (example)
#######################################################################

from applications import Command
if __name__ == "__main__":
    # Examples
    logging.getLogger().setLevel(logging.DEBUG)
    #
    # Set up and start the scheduler
    sched = SimpleScheduler(max_concurrent=4)
    sched.start()
    #
    # Add some jobs to run but don't wait
    # for them to complete
    sched.submit(Command('sh','-c','echo $HOSTNAME'))
    #
    # Run a job where we do want to wait for
    # completion
    sched.submit(Command('sleep','50'),name="sleeper_50")
    sched.submit(Command('du','-sh','..')).wait()
    sched.submit(Command('sh','-c','echo Sleeper 50 finished'),
                 wait_for=('sleeper_50',))
    sched.submit(Command('sleep','10')).wait()
    # More jobs
    sched.submit(Command('ps','faux'))
    sched.submit(Command('du','-sh','.'))
    sched.submit(Command('sleep','20'),name='sleeper_20')
    sched.submit(Command('sh','-c','echo Both sleepers finished'),
                 wait_for=('sleeper_20','sleeper_50',))
    sched.submit(Command('ps','faux'),wait_for=('sleeper_20',))
    sched.submit(Command('du','-sh','.'),wait_for=('sleeper_20',))
    sched.submit(Command('sh','-c','echo Sleeper 20 finished'),
                 wait_for=('sleeper_20',))
    #
    # Wait for all jobs to finish
    while not sched.is_empty():
        print "n_waiting = %s" % sched.n_waiting
        print "n_running = %s" % sched.n_running
        time.sleep(10)
    #
    # Stop the scheduler
    sched.stop()
    
