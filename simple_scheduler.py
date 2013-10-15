import JobRunner
from Pipeline import Job
from applications import Command
import time
import os
import threading
import Queue

# Python module to provide scheduler capability for functions
# and programs.

class SimpleScheduler(threading.Thread):
    # SimpleScheduler
    #
    # Simple scheduler providing simple job control
    #
    # Usage:
    #
    # >>> s = SimpleScheduler()
    # >>> s.start()
    # >>> s.submit(['fastq_screen',...],name='fastq_screen.model')
    # >>> s.submit(['fastq_screen',...],name='fastq_screen.other')
    # >>> s.submit(['fastq_screen',...],name='fastq_screen.rRNA')
    # >>> s.submit(['fastqc',...],wait_for=('fastq_screen.model',...))
    #
    # 
    #
    def __init__(self,runner=JobRunner.SimpleJobRunner(),max_concurrent=None,
                 poll_interval=5):
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
        # Flag controlling whether scheduler is active
        self.__active = False

    def stop(self):
        # Stop the scheduler
        self.__active = False

    @property
    def n_waiting(self):
        # Return number of jobs left
        return len(self.__scheduled)

    @property
    def n_running(self):
        # Return number of jobs still running
        return len(self.__running)

    @property
    def job_number(self):
        # Increment and return job count
        self.__job_count += 1
        return self.__job_count

    def submit(self,args,runner=None,name=None,wait_for=[]):
        # Schedule job to be run
        #
        # Use a queue rather than modifying the waiting list
        # directly to try and avoid
        #
        # Check names are not duplicated
        if name is not None:
            if name in self.__names:
                raise Exception,"Name '%s' already assigned" % name
            self.__names.append(name)
        # Check we're not waiting on a non-existent name
        for job_name in wait_for:
            if job_name not in self.__names:
                raise Exception,"Job depends on a non-existent name '%s'" % job_name
        # Use default runner if none explicitly specified
        if runner is None:
            runner = self.__runner
        # Schedule the job
        job = SchedulerJob(runner,args,job_number=self.job_number,
                           name=name,wait_for=wait_for)
        self.__submitted.put(job)
        print "Scheduled job #%d: \"%s\"" % (job.job_number,job)
        return job

    def run(self):
        # run method overrides that from base Thread class
        # Don't call this directly
        #
        # Implements the scheduler
        #
        print "Starting simple scheduler"
        self.__active = True
        while self.__active:
            # Check for completed jobs
            updated_running_list = []
            for job in self.__running:
                if job.is_running:
                    print "Job #%s (id %s) still running \"%s\"" % (job.job_number,
                                                                    job.job_id,
                                                                    job)
                    updated_running_list.append(job)
                else:
                    print "Job #%s (id %s) completed \"%s\"" % (job.job_number,
                                                                job.job_id,
                                                                job)
                    if job.job_name is not None:
                        self.__finished_names.append(job.job_name)
            # Update the list of running jobs
            self.__running = updated_running_list
            # Add submitted jobs to the waiting list
            while not self.__submitted.empty():
                job = self.__submitted.get()
                self.__scheduled.append(job)
                print "Added job #%d (%s): \"%s\"" % (job.job_number,job.name,job)
            # Start running jobs
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
                    print "Started job #%s (id %s)" % (job.job_number,job.job_id)
                else:
                    # Hold back for now
                    remaining_jobs.append(job)
            # Update the scheduled job list
            self.__scheduled = remaining_jobs
            # Wait before going round again
            time.sleep(self.__poll_interval)

    def is_empty(self):
        # Check if the scheduler has any jobs remaining
        return not (self.n_waiting or self.n_running)

class SchedulerJob(Job):
    # Defines a job within the scheduler
    def __init__(self,runner,args,job_number=None,name=None,wait_for=[]):
        self.job_number = job_number
        self.job_name = name
        self.waiting_for = wait_for
        if name is None:
            name = args[0]
        Job.__init__(self,runner,name,os.getcwd(),args[0],args[1:])
    @property
    def is_running(self):
        # Check if job is running
        return self.isRunning()
    def wait(self,poll_interval=5):
        # Block until the job has completed
        print "Waiting for job #%s..." % self.job_number
        while self.job_id is None or self.is_running:
            time.sleep(poll_interval)
        print "Job #%s finished" % self.job_number
    def __repr__(self):
        return ' '.join([str(x) for x in ([self.script,] + self.args)])

if __name__ == "__main__":
    # Examples
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
    ##sched.stop()
    
