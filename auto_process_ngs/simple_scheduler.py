#!/usr/bin/env python
#
#     simple_scheduler.py: provide basic scheduler capability
#     Copyright (C) University of Manchester 2013-19 Peter Briggs
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

__version__ = "0.0.15"

#######################################################################
# Import modules that this module depends on
#######################################################################

import bcftbx.JobRunner as JobRunner
from bcftbx.Pipeline import Job
import time
import os
import sys
import re
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

    def __init__(self,
                 runner=None,
                 reporter=None,
                 max_concurrent=None,
                 poll_interval=5,
                 job_interval=0.1,
                 max_restarts=1):
        """Create a new SimpleScheduler instance

        Arguments:
          runner: optional, default job runner to use
          reporter: optional, custom SchedulerReporter instance to use
            for messaging when jobs, groups etc are started and finished
          max_concurrent: optional, maximum number of concurrent
            processes the scheduler will run (default: no limit)
          poll_interval: optional, number of seconds to wait in
            between checking for completed jobs etc in the scheduler
            loop (default: 5 seconds)
          job_interval: optional, number of seconds to wait before
            submitting a job (default: 0.1 seconds)
          max_restarts: optional, if non-zero then attempt to restart
            jobs that are in an error state up to this many times. Set to
            zero to turn off restarting jobs in error states (they will
            be terminated instead). (default: 1)

        """

        threading.Thread.__init__(self)
        self.setDaemon(1)
        # Default job runner
        if runner is None:
            runner = JobRunner.SimpleJobRunner()
        self.__runner = runner
        # Maximum number of concurrent jobs
        self.__max_concurrent = max_concurrent
        # Length of time to wait between checking jobs
        self.__poll_interval = poll_interval
        # Length of time to wait before submitting job
        self.__job_interval = job_interval
        # Number of attempts to restart errored jobs
        self.__max_restarts = max_restarts
        # Internal job id counter
        self.__job_count = 0
        # Queue to add jobs
        self.__submitted = Queue.Queue()
        # List of scheduled (i.e.waiting) jobs
        self.__scheduled = []
        # List of running jobs
        self.__running = []
        # Dictionary with all jobs
        self.__jobs = dict()
        # Handle names
        self.__names = []
        self.__finished_names = []
        # Handle groups
        self.__active_groups = []
        self.__groups = dict()
        # Handle callbacks
        self.__callbacks = []
        # Error checking for jobs
        self.__enable_error_check = True
        self.__last_error_check_time = 0.0
        self.__error_check_interval = 60.0
        # Flag controlling whether scheduler is active
        self.__active = False
        # Default reporter
        if reporter is None:
            reporter = default_scheduler_reporter()
        self.__reporter = reporter

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
    def n_finished(self):
        """Return number of jobs that have completed

        """
        return len(self.__jobs) - self.n_waiting - self.n_running

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
        return not (self.n_waiting or self.n_running or not self.__submitted.empty())

    def lookup(self,name):
        """Look up and return SchedulerJob or SchedulerGroup instance

        Returns None if no matching instance was found.

        """
        if name in self.__jobs:
            return self.__jobs[name]
        if name in self.__groups:
            return self.__groups[name]
        return None

    def find(self,pattern):
        """Lookup groups and jobs from regex pattern matching

        Returns a list of jobs with names that match the supplied
        pattern (or an empty list if no matches were found).

        """
        matches = []
        regex = re.compile(pattern)
        for name in self.__groups:
            if regex.match(name):
                matches.append(self.__groups[name])
        for name in self.__jobs:
            if regex.match(name):
                matches.append(self.__jobs[name])
        matches.sort(lambda x,y: cmp(x.name,y.name))
        return matches

    def wait(self):
        """Wait until the scheduler is empty

        """
        try:
            while not self.is_empty():
                time.sleep(self.__poll_interval)
        except KeyboardInterrupt:
            print("KeyboardInterrupt")
            self.stop()
            print("Terminating running jobs:")
            for job in self.__running:
                print("\t#%d (%s): \"%s\" (%s)" % (job.job_number,
                                                   job.job_id,
                                                   job.name,
                                                   date_and_time(job.start_time)))
                job.terminate()
            print("Finished")

    def wait_for(self,names,timeout=None):
        """Wait until the named jobs/groups have completed

        Arguments:
          names: a list or tuple of job and/or group names
            which the scheduler will wait to complete
          timeout: optional, if set then is the maximum time
            in seconds that the job will be allowed to run before
            it's terminated and a SchedulerTimeout exception
            is raised

        """
        try:
            wait_time = 0
            while True:
                completed = True
                for name in names:
                    completed = (completed and (name in self.__finished_names))
                if completed:
                    return
                if timeout is not None and wait_time > timeout:
                    raise SchedulerTimeout(
                        "Timeout exceeded waiting for %s (%ss)" %
                        (names,timeout))
                time.sleep(self.__poll_interval)
        except KeyboardInterrupt:
            print("KeyboardInterrupt")
            self.stop()
            print("Terminating running jobs:")
            for job in self.__running:
                print("\t#%d (%s): \"%s\" (%s)" % (job.job_number,
                                                   job.job_id,
                                                   job.name,
                                                   date_and_time(job.start_time)))
                job.terminate()
            print("Finished")

    def submit(self,args,runner=None,name=None,wd=None,log_dir=None,
               wait_for=None,callbacks=[]):
        """Submit a request to run a job
        
        Arguments:
          args: the command to run expressed as a list or tuple of
                arguments
          runner: (optional) a JobRunner instance that will be used to
                dispatch and control the job.
          name: (optional) a name for the job. Must be unique within
                the scheduler instance.
          wd:   (optional) the working directory to execute the job in;
                defaults to the current working directory
          log_dir: (optional) explicitly specify directory for log files
          wait_for: (optional) a list or tuple of job and/or group
                names which must finish before this job can start
          callbacks: (optional) a list or tuple of functions that will
                be executed when the job completes.

        Returns:
          SchedulerJob instance for the submitted job.

        """
        # Use a queue rather than modifying the waiting list
        # directly to try and avoid
        #
        # Fetch a unique id number
        job_number = self.job_number
        # Generate a name if necessary
        if name is None:
            name = "%s.%s" % (str(args[0]),job_number)
        # Check names are not duplicated
        if self.has_name(name):
            raise Exception("Name '%s' already assigned" % name)
        self.__names.append(name)
        # Check we're not waiting on a non-existent name
        if wait_for:
            for job_name in wait_for:
                if not self.has_name(job_name):
                    raise Exception("Job depends on a non-existent name "
                                    "'%s'" % job_name)
        # Use default runner if none explicitly specified
        if runner is None:
            runner = self.default_runner
        # Pause before submitting
        time.sleep(self.__job_interval)
        # Schedule the job
        job = SchedulerJob(runner,args,job_number=job_number,
                           name=name,working_dir=wd,log_dir=log_dir,
                           wait_for=wait_for)
        self.__submitted.put(job)
        self.__jobs[job.job_name] = job
        # Deal with callbacks
        for function in callbacks:
            self.callback("callback.%s" % job.job_name,
                          function,wait_for=(job.job_name,))
        self.__reporter.job_scheduled(job)
        logging.debug("%s" % job)
        return job

    def group(self,name,log_dir=None,wait_for=None,callbacks=[]):
        """Create a group of jobs
        
        Arguments:
          name: a name for the group. Must be unique within the
                scheduler instance.
          wait_for: (optional) a list or tuple of job and/or group
                names which must finish before this job can start
          log_dir: (optional) explicitly specify directory for log files
          callbacks: (optional) a list or tuple of functions that will
                be executed when the job completes.

        Returns:
          Empty SchedulerGroup instance.

        """
        # Fetch a unique id number
        job_number = self.job_number
        # Generate a name if necessary
        if name is None:
            name = "group.%s" % job_number
        # Check names are not duplicated
        if self.has_name(name):
            raise Exception("Name '%s' already assigned" % name)
        self.__names.append(name)
        new_group = SchedulerGroup(name,job_number,self,log_dir=log_dir,
                                   wait_for=wait_for)
        self.__groups[name] = new_group
        self.__active_groups.append(name)
        # Deal with callbacks
        for function in callbacks:
            self.callback("callback.%s" % new_group.group_name,
                          function,
                          wait_for=(new_group.group_name,))
        self.__reporter.group_added(new_group)
        return new_group

    def callback(self,name,callback,wait_for):
        """Add a callback function

        Add a function 'callback' which will be invoked when
        all the jobs listed in 'wait_for' have completed.

        The function is invoked as:

           callback(name,(job1,...),sched)

        where name is the submitted name, job is the
        SchedulerJob that has completed, and sched is the
        scheduler instance.

        Arguments:
          name: a name for the callback. Must be unique within the
                scheduler instance; if None then a unique name will
                be generated.
          callback: the function that will be executed when the job
                completes.
          wait_for: a list or tuple of job and/or group names which
                will trigger the callback when they finish

        Returns:
          SchedulerCallback instance.

        """
        # Fetch a unique id number
        job_number = self.job_number
        # Generate a name if necessary
        if name is None:
            name = "callback.%s" % job_number
        # Check names are not duplicated
        if self.has_name(name):
            raise Exception("Name '%s' already assigned" % name)
        new_callback = SchedulerCallback(name,callback,wait_for=wait_for)
        self.__callbacks.append(new_callback)
        return new_callback

    def run(self):
        """Internal: run method overriding that from base Thread class

        This method implements the scheduler loop.

        Don't call this directly - it is invoked by the Thread's 'start'
        method.

        """
        logging.debug("Starting simple scheduler")
        self.__active = True
        while self.__active:
            # Flag to indicate status should be reported
            report_status = False
            # Flag to indicate whether to error check jobs
            if self.__enable_error_check \
               and ((time.time() - self.__last_error_check_time) >
                    self.__error_check_interval):
                self.__last_error_check_time = time.time()
                do_error_check = True
            # Check for completed jobs
            updated_running_list = []
            for job in self.__running:
                if job.is_running:
                    # Check if job is error
                    if (not do_error_check) or (not job.in_error_state):
                        logging.debug("Job #%s (id %s) running \"%s\""
                                      % (job.job_number,
                                         job.job_id,
                                         job))
                        updated_running_list.append(job)
                        continue
                    # Handle job in error state
                    logging.debug("Job #%s (id %s) in error state \"%s\""
                                  % (job.job_number,
                                     job.job_id,
                                     job))
                    logging.warning("Job #%s (id %s) in error state"
                                    % (job.job_number,
                                       job.job_id))
                    if not self.__max_restarts:
                        job.terminate()
                        logging.warning("Job #%s (id %s) terminated"
                                        % (job.job_number,
                                           job.job_id))
                    elif job.restart(max_tries=self.__max_restarts):
                        logging.warning("Job #%s (id %s) restarted"
                                        % (job.job_number,
                                           job.job_id))
                        self.__reporter.job_start(job)
                        updated_running_list.append(job)
                    else:
                        logging.warning("Job #%s (id %s) failed to "
                                        "restart" % (job.job_number,
                                                     job.job_id))
                        self.__reporter.job_end(job)
                        self.__finished_names.append(job.job_name)
                    report_status = True
                else:
                    self.__reporter.job_end(job)
                    logging.debug("Job #%s (id %s) completed \"%s\"" % (job.job_number,
                                                                        job.job_id,
                                                                        job))
                    report_status = True
                    if job.job_name is not None:
                        self.__finished_names.append(job.job_name)
            # Update the list of running jobs
            self.__running = updated_running_list
            # Check for completed groups
            updated_groups = []
            for group_name in self.__active_groups:
                group = self.__groups[group_name]
                if group.closed:
                    if group.is_running:
                        logging.debug("Group #%s (id %s) still running" % (group.group_name,
                                                                           group.group_id))
                        updated_groups.append(group_name)
                    else:
                        self.__reporter.group_end(group)
                        logging.debug("Group #%s (id %s) completed" % (group.group_name,
                                                                       group.group_id))
                        self.__finished_names.append(group_name)
                        report_status = True
                else:
                    logging.debug("Group #%s (id %s) waiting for more jobs" %
                                  (group.group_name,group.group_id))
                    updated_groups.append(group_name)
            # Update the list of groups
            self.__active_groups = updated_groups
            # Handle callbacks
            updated_callback_list = []
            for callback in self.__callbacks:
                logging.debug("Checking %s" % callback.callback_name)
                logging.debug("Waiting_for: %s" % callback.waiting_for)
                invoke_callback = True
                callback_jobs = []
                for name in callback.waiting_for:
                    if name in self.__finished_names:
                        logging.debug("- '%s' has finished" % name)
                        callback_jobs.append(self.lookup(name))
                    else:
                        # Waiting for at least one job
                        logging.debug("- still waiting for '%s'" % name)
                        invoke_callback = False
                        break
                if invoke_callback:
                    logging.debug("Invoking callback '%s'" % callback.callback_name)
                    callback.invoke(tuple(callback_jobs),self)
                else:
                    updated_callback_list.append(callback)
            self.__callbacks = updated_callback_list
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
                    try:
                        job.start()
                        self.__running.append(job)
                        self.__reporter.job_start(job)
                        logging.debug("Started job #%s (id %s)" % (job.job_number,job.job_id))
                    except Exception as ex:
                        logging.error("Failed to start job #%s: %s" % (job.job_number,ex))
                        if job.job_name is not None:
                            self.__finished_names.append(job.job_name)
                    report_status = True
                else:
                    # Hold back for now
                    remaining_jobs.append(job)
            # Update the scheduled job list
            self.__scheduled = remaining_jobs
            # Report current status, if required
            if report_status:
                self.__reporter.scheduler_status(self)
            # Wait before going round again
            time.sleep(self.__poll_interval)

class SchedulerGroup(object):
    """Class providing an interface to schedule a group of jobs

    SchedulerGroup instances should normally be returned by a
    call to the 'group' method of a SimpleScheduler object.
    The group should be populated by calling its 'add' method
    (note that jobs are passed directly to the scheduler).

    Once all jobs have been added the 'close' method should
    be invoked to indicate to the scheduler that the group is
    complete. At this point no more jobs can be added to the
    group, and the scheduler will check for when the group
    has finished running.
    
    """

    def __init__(self,name,group_id,parent_scheduler,log_dir=None,
                 wait_for=None):
        """Create a new SchedulerGroup instance

        Arguments:
          name: a name for the group. Must be unique within the
                parent scheduler instance.
          group_id: a unique id number
          parent_scheduler: the SimpleScheduler instance that the
                group belongs to
          log_dir: (optional) explicitly specify directory for log files
          wait_for: (optional) a list or tuple of job and/or group
                names which must finish before this job can start

        """
        self.group_name = name
        self.group_id = group_id
        if wait_for:
            self.waiting_for = list(wait_for)
        else:
            self.waiting_for = list()
        self.log_dir = log_dir
        self.__scheduler = parent_scheduler
        self.__closed = False
        self.__jobs = []

    @property
    def name(self):
        return self.group_name

    @property
    def closed(self):
        """Test whether group has been closed

        """
        return self.__closed

    @property
    def is_running(self):
        """Test whether group is running

        Returns True if group was closed and if it has jobs
        that are still running, False if not.

        """
        if not self.closed:
            return False
        for job in self.__jobs:
            if not job.submitted:
                return True
            elif not job.completed:
                return True
        return False

    @property
    def completed(self):
        """Test whether group has completed

        Returns True if all the jobs in the group have
        completed, False otherwise.
        """
        if not self.closed:
            return False
        for job in self.__jobs:
            if not job.completed:
                return False
        return True

    @property
    def exit_code(self):
        """Return exit code for the group

        If all jobs have completed with status zero
        then returns zero, otherwise returns a count
        of jobs which have non-zero status.

        If the group hasn't completed then returns
        None.
        """
        exit_code = 0
        if self.completed:
            for job in self.__jobs:
                if not job.completed:
                    return None
                if job.exit_code != 0:
                    exit_code += 1
            return exit_code
        else:
            return None

    def add(self,args,runner=None,name=None,wd=None,log_dir=None,
            wait_for=None,callbacks=[]):
        """Add a request to run a job
        
        Arguments:
          args: the command to run expressed as a list or tuple of
                arguments
          runner: (optional) a JobRunner instance that will be used to
                dispatch and control the job.
          name: (optional) a name for the job. Must be unique within
                the scheduler instance.
          wd:   (optional) the working directory to execute the job in;
                defaults to the current working directory
          log_dir: (optional) explicitly specify directory for log files
          wait_for: (optional) a list or tuple of job and/or group
                names which must finish before this job can start
          callbacks: (optional) a list or tuple of functions that will
                be executed when the job completes.

        Returns:
          SchedulerJob instance for the added job.

        """
        # Check we can still add jobs
        if self.closed:
            raise Exception("Can't add job to group '%s': group closed "
                            "to new jobs" % self.group_name)
        # Deal with directory for log files
        if log_dir is None:
            log_dir = self.log_dir
        # Update list of jobs that this one needs to wait for
        if wait_for:
            waiting_for = self.waiting_for + list(wait_for)
        else:
            waiting_for = self.waiting_for
        # Submit the job to the scheduler and keep a reference
        logging.debug("Group '%s' #%s: adding job" % (self.group_name,self.group_id))
        job = self.__scheduler.submit(args,runner=runner,name=name,
                                      wd=wd,log_dir=log_dir,
                                      wait_for=waiting_for,
                                      callbacks=callbacks)
        self.__jobs.append(job)
        return job

    @property
    def jobs(self):
        """Return list of jobs

        """
        return self.__jobs

    def close(self):
        """Indicate that all jobs have been added to the group

        """
        if self.closed:
            raise Exception("Group '%s' already closed" % self.group_name)
        logging.debug("Group '%s' #%s closed" % (self.group_name,self.group_id))
        self.__closed = True

    def wait(self,poll_interval=5):
        """Wait for the group to complete

        Arguments:
          poll_interval: optional, number of seconds to wait in
            between checking if the group has completed (default: 5
            seconds)

        """
        if not self.closed:
            raise Exception("Group '%s' not closed" % self.group_name)
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

    def __init__(self,runner,args,job_number=None,name=None,working_dir=None,
                 log_dir=None,wait_for=None):
        """Create a new SchedulerJob instance

        """
        self.job_number = job_number
        self.job_name = name
        self.log_dir = log_dir
        if wait_for:
            self.waiting_for = list(wait_for)
        else:
            self.waiting_for = list()
        self.command = ' '.join([str(arg) for arg in args])
        if name is None:
            name = args[0]
        if working_dir is None:
            working_dir = os.getcwd()
        else:
            working_dir = os.path.abspath(working_dir)
        Job.__init__(self,runner,name,working_dir,args[0],args[1:])
        self._restarts = 0

    @property
    def is_running(self):
        """Test if a job is running

        Returns True if the job has started and is still running,
        False otherwise.

        """
        # Check if job is running
        return self.isRunning()

    @property
    def in_error_state(self):
        """Test if a job is in an error state

        Returns True if the job appears to be in an error state,
        and False if not.
        """
        return self.errorState()

    @property
    def completed(self):
        """Test if a job has completed

        Returns True if the job has finished running, False
        otherwise.
        """
        # NB need to invoke is_running to force implicit update
        # of job status
        return ((not self.is_running) and (self.end_time is not None))

    @property
    def exit_code(self):
        """Return exit code from job

        Wrapper for the 'exit_status' property provided by the
        'Job' superclass.
        """
        return self.exit_status

    def start(self):
        """Start the job running

        Overrides the 'start' method in the base 'Job' class.

        Returns:
          Id for job

        """
        if self.log_dir is not None:
            runner_log_dir = self.runner.log_dir
            self.runner.set_log_dir(self.log_dir)
        job_id = Job.start(self)
        if self.log_dir is not None:
            self.runner.set_log_dir(runner_log_dir)
        return job_id

    def wait(self,poll_interval=5,timeout=None):
        """Wait for the job to complete

        This method blocks while waiting for the job to finish
        running.

        NB if the job is not created by submission to a scheduler
        or group then it's up to the calling subprogram to ensure
        that the job is started before the 'wait' method is
        invoked.

        Arguments:
          poll_interval: optional, number of seconds to wait in
            between checking if the job has completed (default: 5
            seconds)
          timeout: optional, if set then is the maximum time
            in seconds that the job will be allowed to run before
            it's terminated and a SchedulerTimeout exception
            is raised

        """
        logging.debug("Waiting for job #%s (%s)..." % (self.job_number,
                                                       self.job_id))
        wait_time = 0
        while ((self.job_id is None) or (not self.completed) or
               (self.exit_status is None)):
            # Check for timeout
            if timeout is not None and wait_time > timeout:
                self.terminate()
                raise SchedulerTimeout(
                    "Job #%s (%s): timeout exceeded (%ss)" %
                    (self.job_number,self.job_id,timeout))
            # Wait before polling again
            time.sleep(poll_interval)
            wait_time += poll_interval
        logging.debug("Job #%s finished" % self.job_number)

    def restart(self,max_tries=3):
        """Restart running the job

        Attempts to restart a job, by terminating the current
        instance and reinvoking the 'start' method.

        The number of times that a restart should be attempted
        is limited by the 'max_tries' parameter.

        If the restart is successful then the new job id is
        returned; otherwise None is returned to indicate
        failure.

        Arguments:
          max_tries (int): maximum number of restarts
            that should be attempted on this job before
            giving up (default: 3)

        Returns:
          Id for job

        """
        # Check if job can be restarted
        if self.completed:
            logging.debug("Job #%s: already completed, cannot "
                          "restart" % self.job_number)
            return False
        # Terminate the job
        logging.debug("Job #%s: terminating before resubmission"
                      % self.job_number)
        self.terminate()
        # Check if we can restart
        self._restarts += 1
        if self._restarts > max_tries:
            logging.debug("Job #%s: maximum number of restart "
                          "attempts exceeded (%d)" %
                          (self.job_number,max_tries))
            return False
        else:
            logging.debug("Job #%s: restarted (attempt #%d)" %
                          (self.job_number,self._restarts))
            #return self.start()
            return Job.restart(self)

    def __repr__(self):
        """Return string representation of the job command line

        """
        return ' '.join([str(x) for x in ([self.script,] + self.args)])

class SchedulerCallback(object):
    """Class providing an interface to scheduled callbacks

    SchedulerJob instances should normally be returned by a
    call to the 'submit' method of a SimpleScheduler object.

    """

    def __init__(self,name,function,wait_for=[]):
        """Create a new SchedulerCallback instance

        """
        self.callback_name = name
        self.callback_function = function
        self.waiting_for = list(wait_for)

    def invoke(self,jobs,sched):
        """Invoke the callback

        """
        logging.debug("Invoking callback function: %s" % self.callback_name)
        try:
            return self.callback_function(self.callback_name,jobs,sched)
        except Exception as ex:
            logging.error("Exception invoking callback function '%s': %s (ignored)" % \
                          (self.callback_name,ex))

class SchedulerReporter(object):
    """Class to report on scheduler operations

    SimpleScheduler calls methods of the SchedulerReporter when
    jobs are scheduled, started and finished, when groups are
    added and finished, and to report the status of the scheduler.

    The methods look for associated template format strings which
    are then used to create and output strings.

    Templates are added and modified via the 'set_template' method.
    Templates should be Python format strings which should process
    a dictionary of values.

    Template name     Used when           Available dictionary keys
    ---------------------------------------------------------------
    job_scheduled     Job is submitted    job_name, job_number,job_id,
                                          waiting_for, time_stamp,
                                          command
    job_start         Job starts running  as 'job_scheduled'
    job_end           Job finishes        as 'job_scheduled'
    group_added       Group is created    group_name, group_id, time_stamp
    group_end         Group completes     as 'group_added'
    scheduler_status  Need status         n_running, n_waiting, n_finished

    An example template string for a job could be:

    'Job %(job_name)s: %(job_id)d'

    The purpose of the SchedulerReporter class to provide a way to
    customise reporting of the standard scheduler operations when
    used in an application.

    """
    def __init__(self,fp=sys.stdout,**args):
        """Create new SchedulerReporter instance

        Arguments
          fp: optional, if provided then must be a file-like
              object opened for writing (defaults to stdout)
          args: optional keyword-value pairs, where 'value'
              defines a template for the name provided by
              'keyword'

        """
        self.__template_list = ['job_scheduled',
                                'job_start',
                                'job_end',
                                'group_added',
                                'group_start',
                                'group_end',
                                'scheduler_status']
        self.__templates = {}
        self.__fp = fp
        for name in args:
            self.set_template(name,args[name])

    def _scheduler_dict(self,sched):
        """Return dictionary of keywords derived from scheduler instance

        Arguments:
          sched: SimpleScheduler instance

        Returns:
          Dictionary containing appropriate keywords

        """
        return { 'n_running'  : sched.n_running,
                 'n_waiting'  : sched.n_waiting,
                 'n_finished' : sched.n_finished,
                 'time_stamp' : date_and_time() }

    def _job_dict(self,job):
        """Return dictionary of keywords derived from job instance

        Arguments:
          job: SchedulerJob instance

        Returns:
          Dictionary containing appropriate keywords

        """
        return { 'job_number' : job.job_number,
                 'job_name'   : job.job_name,
                 'job_id'     : job.job_id,
                 'waiting_for': ', '.join(["'%s'" % s for s in job.waiting_for]) if job.waiting_for else '',
                 'time_stamp' : date_and_time(),
                 'command'    : job.command }

    def _group_dict(self,group):
        """Return dictionary of keywords derived from group instance

        Arguments:
          group: SchedulerGroup instance

        Returns:
          Dictionary containing appropriate keywords

        """
        return { 'group_name' : group.group_name,
                 'group_id'   : group.group_id,
                 'time_stamp' : date_and_time() }

    def scheduler_status(self,sched):
        """Write report string for scheduler status

        Arguments:
          sched: SimpleScheduler instance

        """
        self._report('scheduler_status',**self._scheduler_dict(sched))

    def job_scheduled(self,job):
        """Write report string when a job is scheduled

        Arguments:
          job: SchedulerJob instance

        """
        self._report('job_scheduled',**self._job_dict(job))

    def job_start(self,job):
        """Write report string when a job starts

        Arguments:
          job: SchedulerJob instance

        """
        self._report('job_start',**self._job_dict(job))

    def job_end(self,job):
        """Write report string when a job ends

        Arguments:
          job: SchedulerJob instance

        """
        self._report('job_end',**self._job_dict(job))

    def group_added(self,group):
        """Write report string when a group is added

        Arguments:
          group: SchedulerGroup instance

        """
        self._report('group_added',**self._group_dict(group))

    def group_end(self,group):
        """Write report string when a group ends

        Arguments:
          group: SchedulerGroup instance

        """
        self._report('group_end',**self._group_dict(group))

    def set_template(self,name,template):
        """Associate a template string with an operation

        Arguments:
          name: operation name e.g. 'job_start'
          template: associated template string, or None to
            suppress reporting of the operation

        """
        if name not in self.__template_list:
            raise KeyError("SchedulerReporter: name '%s' not "
                           "defined" % name)
        if not template.endswith('\n'):
            template += '\n'
        self.__templates[name] = template

    def _report(self,name,**args):
        """Generate and write report string for operation

        Internal function that should not be called directly.
        If a matching template isn't found, or if an exception
        occurs when attempting to generate the report string,
        then write nothing.

        The result is written to the stream specified when the
        SchedulerReporter instance was created.

        Arguments:
          name: operation name e.g. 'job_start'
          args: keyword-value pairs from appropriate dictionary

        """
        try:
            template = self.__templates[name]
            if template is None:
                return
            self.__fp.write(template % args)
        except KeyError as ex:
            logging.debug("SchedulerReporter: exception '%s' (ignored)" % ex)
            return

#######################################################################
# Functions
#######################################################################

def date_and_time(epoch=None):
    """Return formatted date and time information

    """
    if epoch is not None:
        return time.asctime(time.localtime(epoch))
    else:
        return time.asctime()

def default_scheduler_reporter():
    """Return a default SchedulerReporter object

    """
    return SchedulerReporter(
        scheduler_status="%(time_stamp)s: %(n_running)d running, %(n_waiting)d waiting, %(n_finished)d finished",
        job_scheduled="Job scheduled: #%(job_number)d: \"%(job_name)s\" (%(time_stamp)s)",
        job_start="Job started: #%(job_number)d (%(job_id)s): \"%(job_name)s\" (%(time_stamp)s)",
        job_end="Job completed: #%(job_number)d (%(job_id)s): \"%(job_name)s\" (%(time_stamp)s)",

        group_added="Group has been added: #%(group_id)d: \"%(group_name)s\" (%(time_stamp)s)",
        group_end="Group completed: #%(group_id)d: \"%(group_name)s\" (%(time_stamp)s)"
    )

#######################################################################
# Exception classes
#######################################################################

class SchedulerException(Exception):
    """Base class for errors with simple scheduler code"""

class SchedulerTimeout(SchedulerException):
    """Timeout limit exceeded"""

#######################################################################
# Main program (example)
#######################################################################

from applications import Command
if __name__ == "__main__":
    # Examples
    ##logging.getLogger().setLevel(logging.DEBUG)
    #
    # Set up and start the scheduler
    sched = SimpleScheduler(max_concurrent=4,
                            runner=JobRunner.SimpleJobRunner(join_logs=True))
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
    ##while not sched.is_empty():
    ##    time.sleep(10)
    sched.wait()
    #
    # Stop the scheduler
    sched.stop()
    
