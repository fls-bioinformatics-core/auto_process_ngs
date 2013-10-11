import JobRunner
from Pipeline import Job
from applications import Command
import time
import os
import threading
import Queue

# Python module to provide pipelining capability for functions
# and programs.

class SimplePipeline(threading.Thread):
    # SimplePipeline
    #
    # Simple scheduler providing simple job control
    #
    # Usage:
    # >>> p = SimplePipeline()
    # >>> p.start()
    # >>> p.add_step(['fastq_screen',...],name='fastq_screen.model')
    # >>> p.add_step(['fastq_screen',...],name='fastq_screen.other')
    # >>> p.add_step(['fastq_screen',...],name='fastq_screen.rRNA')
    # >>> p.add_step(['fastqc',...],wait_for=('fastq_screen.model',...))
    #
    def __init__(self,runner=JobRunner.SimpleJobRunner(),max_concurrent=None,
                 poll_interval=5):
        threading.Thread.__init__(self)
        self.setDaemon(1)
        # Runner to use
        self.__runner = runner
        # Maximum number of concurrent jobs
        self.__max_concurrent = max_concurrent
        # Length of time to wait between checking jobs
        self.__poll_interval = poll_interval
        # List of pipeline steps
        self.__next_step_id = 0
        self.__steps = []
        # List of running steps (jobs)
        self.__running = []
        # Queue to add jobs
        self.__scheduled = Queue.Queue()
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
        # Return number of steps left
        return len(self.__steps)

    @property
    def n_running(self):
        # Return number of jobs still running
        return len(self.__running)

    @property
    def next_step_id(self):
        # Increment and return step count
        self.__next_step_id += 1
        return self.__next_step_id

    def add_step(self,args,name=None,wait_for=[]):
        # Schedule a step to be added
        # Use a queue rather than modifying the waiting list
        # directly to try and avoid
        #
        # Check names are not duplicated
        if name is not None:
            if name in self.__names:
                raise Exception,"Name '%s' already assigned" % name
            self.__names.append(name)
        # Check we're not waiting on a non-existent name
        for job in wait_for:
            if job not in self.__names:
                raise Exception,"Step depends on a non-existent name '%s'" % job
        # Schedule the step
        step = PipelineStep(args,step_id=self.next_step_id,name=name,wait_for=wait_for)
        self.__scheduled.put(step)
        print "Scheduled step #%d: \"%s\"" % (step.step_id,step)
        return step

    def run(self):
        # run method overrides that from base Thread class
        # Don't call this directly
        #
        # Implements the scheduler
        #
        # Runs in a
        print "Starting simple pipeline"
        self.__active = True
        while self.__active:
            # Check for completed jobs
            updated_running_list = []
            for step in self.__running:
                if step.is_running:
                    print "Step %s (job %s) still running \"%s\"" % (step.step_id,
                                                                     step.job_id,
                                                                     step)
                    updated_running_list.append(step)
                else:
                    print "Step %s (job %s) completed \"%s\"" % (step.step_id,
                                                                 step.job_id,
                                                                 step)
                    if step.name is not None:
                        self.__finished_names.append(step.name)
            # Update the list of running steps
            self.__running = updated_running_list
            # Add scheduled steps to the waiting list
            while not self.__scheduled.empty():
                step = self.__scheduled.get()
                self.__steps.append(step)
                print "Added step #%d (%s): \"%s\"" % (step.step_id,step.name,step)
            # Start running steps
            remaining_steps = []
            for step in self.__steps:
                ok_to_run = True
                if self.__max_concurrent is not None and \
                   self.n_running == self.__max_concurrent:
                    # Scheduler capacity maxed out
                    ok_to_run = False
                else:
                    # Check if step is waiting for another step to finish
                    for name in step.waiting_for:
                        ok_to_run = (ok_to_run and name in self.__finished_names)
                if ok_to_run:
                    # Start the step running
                    step.run(self.__runner)
                    self.__running.append(step)
                    print "Started step %s as job %s" % (step.step_id,step.job_id)
                else:
                    # Hold back for now
                    remaining_steps.append(step)
            # Update the step list
            self.__steps = remaining_steps
            # Wait before going round again
            time.sleep(self.__poll_interval)

class PipelineStep:
    # Define a step in a pipeline
    def __init__(self,args,step_id=None,name=None,wait_for=[]):
        self.step_id = step_id
        self.args = args
        self.name = name
        self.waiting_for = wait_for
        self.__job = None
        self.__runner = None
    def run(self,runner):
        # Start job and return job id
        args = self.args
        self.__job = Job(runner,args[0],os.getcwd(),args[0],args[1:])
        self.__job.start()
    @property
    def job_id(self):
        # Return job id
        if self.__job is not None:
            return self.__job.job_id
    @property
    def is_running(self):
        # Check if step is running
        try:
            return self.__job.isRunning()
        except Exception,ex:
            print "Error checking if step %s is running: %s" % (self.step_id,ex)
            return False
    def __repr__(self):
        return ' '.join([str(x) for x in self.args])

def run_scheduled_job(simple_pipeline,args,poll_interval=5):
    # Run a job via the scheduler
    # Blocks the calling subprogram until the job has run
    #
    # Note that previously scheduled jobs will continue to
    # run in background, only the calling subprogram will be
    # waiting
    step = p.add_step(args)
    print "Scheduling a job in blocking mode (#%s): waiting..." % step.step_id
    while step.job_id is None or step.is_running:
        time.sleep(poll_interval)
    print "Returning control after step #%s completed" % step.step_id
    return step

if __name__ == "__main__":
    # Examples
    #
    # Set up and start the scheduler
    p = SimplePipeline(max_concurrent=4)
    p.start()
    #
    # Add some steps to run but don't wait
    # for them to complete
    p.add_step(Command('sh','-c','echo $HOSTNAME'))
    #
    # Run a job where we do want to wait for
    # completion
    p.add_step(Command('sleep','50'),name="sleeper_50")
    run_scheduled_job(p,Command('du','-sh','..'))
    p.add_step(Command('sh','-c','echo Sleeper 50 finished'),
               wait_for=('sleeper_50',))
    run_scheduled_job(p,Command('sleep','10'))
    # More jobs
    p.add_step(Command('ps','faux'))
    p.add_step(Command('du','-sh','.'))
    p.add_step(Command('sleep','20'),name='sleeper_20')
    p.add_step(Command('sh','-c','echo Both sleepers finished'),
               wait_for=('sleeper_20','sleeper_50',))
    p.add_step(Command('ps','faux'),wait_for=('sleeper_20',))
    p.add_step(Command('du','-sh','.'),wait_for=('sleeper_20',))
    p.add_step(Command('sh','-c','echo Sleeper 20 finished'),
               wait_for=('sleeper_20',))
    #
    # Wait for all jobs to finish
    while p.n_waiting or p.n_running:
        print "n_waiting = %s" % p.n_waiting
        print "n_running = %s" % p.n_running
        time.sleep(10)
    #
    # Stop the pipeline
    ##p.stop()
    
