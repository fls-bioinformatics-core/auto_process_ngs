#!/usr/bin/env python
#
#     pipeliner.py: utilities for building simple pipelines of tasks
#     Copyright (C) University of Manchester 2017-2021 Peter Briggs
#
"""
Module providing utility classes and functions for building simple
'pipelines' of tasks.

The core classes are:

- Pipeline: class for building and executing pipelines
- PipelineTask: class for defining pipeline tasks
- PipelineCommand: class for defining commands that can be used in tasks
  (nb the ``PipelineCommandWrapper`` class is recommended over subclassing
  ``PipelineCommand``)

Additional supporting classes:

- PipelineFunctionTask: subclass of PipelineTask which enables Python
  functions to be run as external processes
- PipelineCommandWrapper: shortcut alternative to PipelineCommand
- PipelineParam: class for passing arbitrary values between tasks
- FileCollector: returning collections of files based on glob patterns
- FunctionParam: PipelineParameter-like deferred function evaluation
- PathJoinParam: PipelineParameter-like dynamic file path joiner
- PathExistsParam: PipelineParameter-like dynamic file existence checker

There are some underlying classes and functions that are intended for
internal use:

- BaseParam: base class for parameter-like classes
- Capturing: capture stdout and stderr from a Python function
- Dispatcher: run a Python function as an external process
- sanitize_name: clean up task and command names for use in pipeline
- collect_files: collect files based on glob patterns
- resolve_parameter: get the value of arbitrary parameter or object

Overview
--------

For the purposes of the ``pipeliner`` module, a 'pipeline' consists
of a set of 'tasks': a task can be independent of any other task, or
it may depend on one or more tasks being completed before it can
start.

Defining tasks: examples
------------------------

Tasks must be defined by subclassing the ``PipelineTask`` class and
implementing the following methods:

- init: used to declare any input parameters for the task
- setup: perform any set up (e.g. creating output directories) or
  other arbitary actions; typically it also defines any commands that
  will be sent to the scheduler, however this is optional.
- finish: perform final actions after setup and any defined commands
  have completed (nb ``finish`` is optional)
- output: return any outputs from the task once it has completed.

For example: let's define a task which runs the ``fastqc`` program
on a collection of Fastq files one at a time::

    class RunFastqc(PipelineTask):
        def init(self,fastqs,out_dir):
            self.add_output('out_files',list())
        def setup(self):
            if not os.path.exists(self.args.out_dir):
                os.mkdir(self.args.out_dir)
            for fq in self.args.fastqs:
                self.add_cmd(
                    PipelineCommandWrapper("Run FastQC",
                                           "fastqc",
                                           "-o",self.args.out_dir,
                                           fq))
        def finish(self):
            for fq in self.args.fastqs:
                if fq.endswith(".gz"):
                    fq = os.path.splitext(fq)[0]
                out_file = os.path.join(
                    self.args.out_dir,
                    os.path.splitext(
                        os.path.basename(fq))[0]+"_fastqc.html")
                if not os.path.exists(out_file):
                    self.fail(message="Missing output file: %s" % out_file)
                else:
                    self.output.out_files.append(out_file)

The key features are:

1. The argument list of the ``init`` method defines an arbitrary
   set of parameters which are made available to the other methods
   via the ``self.args`` object.

   The ``init`` method also adds an output called ``out_files``
   which will be used to store the outputs of the task; it can be
   accessed via the ``output`` method.

2. The ``setup`` method creates the output directory if it doesn't
   already exist, and then calls the ``add_cmd`` method to add a
   command for each Fastq file supplied via the ``fastqs`` argument
   (accessed as ``self.args.fastqs``), which will run ``fastqc`` on
   that Fastq.

3. The command to run ``fastqc`` is created by creating a
   ``PipelineCommandWrapper`` instance, which names the command
   and specifies the command and arguments to execute.

4. The commands defined via the ``add_cmd`` method are not guaranteed
   to run sequentially. If there are additional commands that rely on
   the first set of commands finishing, then either append these to the
   commands using the shell ``&&`` notation, or put them into a
   separate task which should be executed afterwards.

5. The ``finish`` method constructs the names of the expected output
   Fastqc HTML files, and adds them to the ``out_files`` list which
   was originally initialised within the ``init`` method.

6. The task can explicitly indicate a failure by calling the ``fail``
   method, as in the ``finish`` method above. Raising an exception
   will also implicitly indicate a failure.

Another example is a task which does a simple-minded read count on
a set of Fastq files::

    class CountReads(PipelineTask):
        def init(self,fastqs):
            self.add_output('counts',dict())
        def setup(self):
            for fq in self.args.fastqs:
                if os.path.splitext(fq)[1] == ".gz":
                   cat = "zcat"
                else:
                   cat = "cat"
                self.add_cmd(
                    PipelineCommandWrapper("Count reads",
                                           "echo","-n",fq,"' '","&&",
                                            cat,fq,"|",
                                            "wc","-l"))
        def finish(self):
            for line in self.stdout.split('\n'):
                if not line:
                    continue
                fq = line.split()[0]
                read_count = int(line.split()[1])/4
                self.outputs.counts[fq] = read_count

The key features are:

1. The ``init`` method initialises the ``counts`` output which will
   be populated in the ``finish`` method, and accessed via the
   ``output`` method.

2. The standard output from the task is available via the ``stdout``
   property of the instance.

3. The ``finish`` method is implemented to extract the line count data
   from standard output and convert this to a read count which is
   stored against the Fastq name.

A final example is a task which filters out the Fastq files which have
non-zero read counts::

    class FilterEmptyFastqs(PipelineTask):
        def init(self,read_counts):
            self.add_output('filtered_fastqs',list())
        def setup(self):
            for fq in self.args.read_counts:
                if self.args.read_counts[fq] > 0:
                     self.output.filtered_fastqs.append(fq)

In this case all the processing is performed by the ``setup`` method;
no commands are defined.

Defining task outputs
---------------------

Outputs should be defined within the ``init`` method, using the
``add_output`` method of the task, for example::

    self.add_output('fastqs',dict())

The outputs can be accessed via the task's `output` property,
which returns an ``AttributeDictionary`` object where the
keys/attributes are those defined by the ``add_output`` calls,
for example::

    task = FilterEmptyFastqs(counts)
    ...
    filtered_fastqs = task.output.fastqs

Typically the values of the outputs are not known before the
task has been run, so the ``setup`` and ``finish`` methods can
be used to populate the outputs when the task completes (as in
the ``FilterEmptyFastqs`` example above).

Using task output as input to another task
------------------------------------------

The key feature of a pipeline is that the outputs from an
upstream task can be used as input to one or more downstream
tasks.

In this case the objects returned by ``output`` are likely to
be passed to other tasks before the task has completed. There
are a number of implications:

1. Tasks that receive output from a preceeding task cannot assume that
   those outputs are ready or complete at the point of initialisation.
   It's therefore recommended that tasks don't attempt to use those
   outputs in their 'init' method - all processing should be deferred
   until the 'setup' method (when preceeding tasks will have completed).

2. Outputs from tasks should be passed by object references (e.g.
   via a list or dictionary, which can be updated by the task after
   being passed, via specialised classes such as ``FileCollector``
   or ``PipelineParam`` - see sections below).

As an example, consider the following task which reverses the order of
a list of items::

    class ReverseList(PipelineTask):
        def init(self,items):
            self.add_output('reversed_list',list())
        def setup(self):
            # Generate the reversed list
            for item in self.args.items[::-1]:
                self.output.reversed_list.append(item)

This might be used as follows::

    reverse = ReverseList("Reverse order of list",[1,2,3])

Subsequently ``reverse.output.reversed_list`` will return the
reference to the output list, which can then be passed to another task.
The list will be populated when the reverse task runs, at which point
the output will be available to the following tasks.

Specialised task input/output classes
-------------------------------------

There are a number of specialised classes which can be used to
pass data from the output of one task into the input of another:

``FileCollector``
*****************

The ``FileCollector`` class enables the collection of files matching
a glob-type pattern. It behaves as an iterator, with the file collection
being deferred until the iteration actually takes place; it can therefore
be used as an output for tasks where a set of files are created but the
precise number or names of the files are not known *a priori*.

For example::

    class MakeFiles(PipelineTask):
        def init(self,d,filenames):
            self.add_output('files',
                            FileCollector(self.args.d,"*"))
        def setup(self):
            # Create a directory and "touch" the files
            os.mkdir(self.args.d)
            for f in self.args.filenames:
                with open(os.path.join(self.args.d,f) as fp:
                    fp.write()

``PipelineParam``
*****************

The ``PipelineParam`` class offers a way to pass strings or numerical
types which would normally be immutable between tasks.

For example::

    import string
    from random import choice

    class RandomString(PipelineTask):
        def init(self,n):
            self.add_output('string',PipelineParam())
        def setup(self):
            # Create a random string of length 'n' characters
            allchar = string.ascii_letters + string.punctuation + string.digits
            s = "".join([choice(allchar) for x in range(self.args.n)])
            # Assign the string to the output
            self.output.string.set(s)

If passed as input to a subsequent task, the ``PipelineParam``
instance will behave as a static value when accessed via
``self.args...``.

``FunctionParam``
*****************

The ``FunctionParam`` class allows function invocations to be passed
to tasks as ``PipelineParam``-like objects, so that the function
evaluation is only performed when the task actually consumes the
object.

The function can be a ``lambda`` function, or a full function
defined elsewhere. For example:

::

    is_dir = FunctionParam(lambda d: os.path.isdir(d),out_dir)

or (more concisely):

::

    out_dir_exists = FunctionParam(os.path.isdir,out_dir)

Arguments and keywords supplied to the ``FunctionParam`` class are
passed to the function for evaluation only when the object's
``value`` method is invoked (with any parameter-like arguments
being converted to their values immediately prior to evaluation).

This can be useful when building a pipeline where functions can't
be evaluated until the pipeline is executed, and as an alternative
to creating a full ``PipelineFunctionTask`` to wrap a function.

``PathJoinParam``
*****************

The ``PathJoinParam`` class provides a parameter-like alternative
to the ``os.path.join`` function; it is useful when building
pipelines where it is not possible to construct the final paths
for files or directories until the pipeline is actually executed
(for example if some path components are supplied as parameters).

``PathExistsParam``
*******************

The ``PathExistsParam`` class provides a parameter-like alternative
to the ``os.path.exists`` function; it is useful when building
pipelines where it is not possible to check the existence of some
paths until the pipeline is actually executed (for example because
files or directories being checked are created during pipeline
execution).

Running Python functions as tasks
---------------------------------

The PipelineFunctionTask class enables a Python function to be run
as an external process, for example reimplementing the earlier
``CountReads`` example::

    from bcftbx.FASTQFile import nreads
    class CountReads(PipelineFunctionTask):
        def init(self,fastqs):
            self.add_output('counts',dict())
        def setup(self):
            self.add_call("Count reads in Fastq",
                          self.count_reads,
                          self.args.fastqs)
        def count_reads(self,fastqs):
            # Count reads in each Fastq
            counts = dict()
            for fq in fastqs:
                counts[fq] = nreads(fq)
            return counts
        def finish(self):
            for fq in self.result():
                for fq in result[fq]:
                    self.output.counts[fq] = result[fq]

The main differences between this and the PipelineTask class are:

1. The class includes a method (in this case ``make_files``) which
   implements the Python function to be executed.

2. The ``add_call`` method is used in the ``setup`` method to define
   calls to the function to run (analogous to the ``add_cmd`` method
   for normal tasks). ``add_call`` can be used multiple times, e.g.
   to break a task into several separate processes (again analogous
   to ``add_cmd``).

3. The return values from the invoked function are collected and
   made available via the ``result`` method. This returns a list
   of results (one for each ``add_call`` was used). The ``finish``
   method contains code to post-process these results.

The advantage of this approach is that intensive operations (e.g.
processing large numbers of files) previously performed by Python
functions can be farmed out to external processes without having
to be wrapped in new executable programs.

.. note::

   The ``FunctionParam`` class could also be considered as an
   alternative to ``PipelineFunctionTasks`` for light-weight
   function calls during pipeline execution.

Building and running a pipeline
-------------------------------

An empty pipeline is created by making a ``Pipeline`` instance::

    ppl = Pipeline()

Specific task instances are then created from PipelineTask subclasses,
for example using the tasks defined previously::

    read_counts = CountReads("Count the reads",fastqs)
    filter_empty_fastqs = FilterEmptyFastqs("Filter empty Fastqs",
                                            read_counts.output.filtered_fastqs)
    run_fastqc = RunFastqc("Run Fastqc",
                           filter_empty_fastqs.output.filtered_fastqs)

Note that when instantiating a task, it must be given a name; this
can be any arbitrary text and is intended to help the end user
distinguish between different task instances in the pipeline.

The tasks are then added to the pipeline using the ``add_task``
method::

    ppl.add_task(read_counts)
    ppl.add_task(filter_empty_fastqs,
                 requires=(read_counts,))
    ppl.add_task(run_fastqc,
                 requires=(filter_empty_fastqs,))

The pipeline is then run using the ``run`` method::

    ppl.run()

Notes:

1. Tasks will only be executed once any tasks they depend on have
   completed successfully; these are specified via the ``requires``
   argument of the ``add_task`` method (in which case it must be a
   list or tuple of task instances), and/or via the ``requires``
   method of a task instance.

2. Tasks that fail (i.e. complete with non-zero exit status) will
   cause the pipeline to halt at that point.

3. The ``run`` method blocks and returns the exit status of the
   pipeline execution.

Defining relationships between tasks
------------------------------------

A task in a pipeline may depend on one or more other tasks in
the pipeline to complete before it can be run; these are referred
to as "requirements" of the task.

Requirements can be specified in different ways:

1. When a task is added to a pipeline via the ``add_task``
   method then a list of required tasks can also be specified
   via the ``requires`` argument;

2. Requirements can also be added directly to a task using
   its ``requires`` method.

Note that these two approaches are not exclusive, and can be
used together on the same task to specify requirements in a
pipeline.

Additionally, requirements are automatically added implicitly for
each input that is a parameter-based output from another task.

Requirement specification can be combined with two properties
that can be useful for when adding tasks to an existing
pipeline:

* ``initial_tasks`` returns a list of all the tasks in the
  pipeline that don't have any requirements; it can be used
  when inserting tasks at the start of a pipeline, as the
  inserted task can be made into a requirement of the initial
  tasks using the idiom
  ``task.required_by(*pipeline.initial_tasks)``, and
* ``final_tasks`` returns a list of the all the tasks in the
  pipeline which no other tasks require (essentially they will
  run at the end of the pipeline); it can be used when appending
  tasks to the end of a pipeline using the idiom
  ``task.requires(*pipeline.final_tasks)``.

Scheduling and running jobs
---------------------------

By default the ``run`` method of the Pipeline instance creates
a new ``SimpleScheduler`` instance internally and uses this to
run the commands generated by the tasks in the pipeline.

The ``run`` method provides a number of arguments that can be
used to configure the scheduler, specifically:

* ``max_jobs``: this sets the maximum number of concurrent
  jobs that the scheduler will run (defaults to 1)
* ``max_slots``: this sets the maximum number of concurrent
  CPUs (aka "slots") that the scheduler will allocate work
  for (defaults to no limit).
* ``poll_interval``: the time interval that the scheduler will
  used when checking the status of running jobs (defaults to
  5 seconds)

If more control is required over the scheduler then the calling
subprogram can create and configure its own ``SimpleScheduler``
instance and pass this to the ``run`` method instead. Note that
it is the responsibility of the subprogram to start and stop
the scheduler that it provides.

Specifying job runners for tasks
--------------------------------

By default when a pipeline is executed then all jobs within
its tasks will be run using the default job runner supplied
via the ``default_runner`` argument of the pipeline's ``run``
method.

For example:

::

    ppl.run(default_runner=SpecialJobRunner())

(If no default is supplied then a ``SimpleJobRunner`` is used
as the default runner.)

Typically it is desirable to be able to exercise more granular
control over the job runners used within the pipeline. This is
possible via the ``runner`` keyword of the ``add_task`` method
of the ``Pipeline`` class, which allows a job runner to be
specified which will then used when jobs in that task are
executed.

For example:

::

    ppl.add_task(my_task,runner=SimpleJobRunner())

The ``Pipeline`` class also allows parameterisation of job
runners when building a pipeline. Placeholders for job
runners can be defined using the ``add_runner`` method,
and accessed via the ``runners`` property to be associated
with tasks.

For example:

::

    # Define runner
    ppl.add_runner('my_runner')
    ...
    # Associate runner with task
    ppl.add_task(my_task,runner=ppl.runners['my_runner'])

Job runners can then be set explicitly when the pipeline is
executed, using the ``runners`` argument of the ``run``
method to supply a mapping of runner names to job runner
instances.

For example:

::

    ppl.run(runners={ 'my_runner': SpecialJobRunner() })

Any runner names that don't have associated job runner
instances will use the default runner defined via the
``default_runner`` argument.

Dynamically setting number of CPUs/threads via job runners
----------------------------------------------------------

When job runners are created they can have a maximum number
of available CPUs (aka "slots") associated with them.

For ``SimpleJobRunner``s this has to be set explicitly via
the ``nslots`` argument, for example:

::

    runner =  SimpleJobRunner(nslots=8)

By default only a single slot is allocated. (For
``GEJobRunners`` the number of slots is set implicitly.)

The number of slots can then be accessed at runtime, so that
jobs run within a task use the appropriate number of CPUs
dynamically, by using the ``runner_nslots`` method.

For standard ``PipelineTask`` classes, this should be done
when constructing commands within the ``setup`` method. For
example: ``bowtie2`` takes a ``--threads`` option which
tells the program how many threads it should use. A minimal
task to run ``bowtie2`` with dynamically assigned number of
threads might look like:

::

   class RunBowtie2(PipelineTask):
        def init(self,fastq,index_basename,sam_out):
            pass
        def setup(self):
            self.add_cmd(
                PipelineCommandWrapper("Run bowtie",
                                       "bowtie2",
                                       "-x",self.args.index_basename,
                                       "-U",self.args.fastq,
                                       "-S",self.args.sam_out,
                                       "--threads",self.runner_nslots)

.. note::

   When using dynamic CPU assignment with ``SimpleJobRunners``,
   it may also be worth considering using the ``max_slots``
   parameter when running the pipeline.

Dealing with stdout from tasks
------------------------------

The stdout from tasks which run external commands can be
accessed via the ``stdout`` property of the task instance once
it has completed.

Where multiple jobs were run by the task, the stdout from all
jobs are concatenated and returned via this property.

The stdout for each job is topped and tailed with a standard
set of comment lines output from the wrapper scripts, of the
form::

    #### COMMAND Echo text
    #### HOSTNAME popov
    #### USER pjb
    #### START Thu Aug 17 08:38:14 BST 2017
    ...
    ...Job-specific output...
    ...
    #### END Thu Aug 17 08:38:14 BST 2017
    #### EXIT_CODE 0

When parsing the stdout it is recommended to check for these lines
using e.g. ``line.startswith("#### ")``.

Handling failed tasks in pipelines
----------------------------------

If a task in a pipeline fails (that is, completes with a
non-zero exit code) then the pipeline is considered to have
failed. In this case the pipeline can use one of a number of
strategies to handle execution of the remaining tasks:

* Pipeline execution halts immediately and all running tasks
  are terminated ('immediate' mode, the default)
* Pipeline execution continues but all tasks which depend on
  the failed tasks are removed and not executed ('deferred'
  mode)

The strategy can be set explicitly at runtime by setting the
``exit_on_failure`` argument of the pipeline ``run`` method
to one of the values defined in the ``PipelineFailure`` class.

For example::

    from pipeliner import PipelineFailure
    ...
    # Define pipeline
    ...
    # Run pipeline in 'deferred' mode
    ppl.run(exit_on_failure=PipelineFailure.DEFERRED)

Note that regardless of how the failures are handled the
pipeline will always return exit code 1 when one or more
tasks fail.

Executing pipeline commands in batches
--------------------------------------

By default when the pipeline executes the commands generated by
a task, each command is sent to the scheduler as a single job.

It is also possible to request the pipeline executes commands
in batches, by specifying a non-zero size for the ``batch_size``
option of the ``run`` method. In this case, commands are grouped
together into batches of this size, and each batch is sent to
the scheduler as a single job.

Within a batch the commands are executed sequentially, and if
one command fails then all subsequent commands in the batch
won't run.

Batch mode can also be requested on a per-task basis, by
explicitly specifying ``batch_size`` as keyword when adding
the task to the pipeline. For example::

    ppl = Pipeline()
    ...
    ppl.add_task(my_task,batch_size=5)

This will override any batch size set globally when the
pipeline is run.

Setting pipeline parameters at execution time
---------------------------------------------

When building pipelines, it is sometimes necessary or desirable
a parameter into a task where the value of the parameter isn't
known until execution time (via the ``run`` method).

For example, a task in the pipeline might need to know the
number of cores or the location of a temporary directory to be
used, which only be set this at execution time.

To handle these situations, it possible to define arbitrary
parameters within the ``Pipeline`` class at build time and then
set the values of these parameters at execution time.

Use the ``add_param`` method is used to define a parameter, for
example:

::

    ppl = Pipeline()
    ppl.add_param('ncores',value=1,type=int)
    ppl.add_param('tmpdir')

This creates a new ``PipelineParam`` instance which is associated
with the supplied name.

The parameters can be accessed via the pipeline's ``params``
property, and passed as input into tasks, for example:

::

    task = ExampleTask("This is an example",
                       ncores=ppl.params.ncores,
                       tmpdir=ppl.params.tmpdir)
    ppl.add_task(task)

The runtime values of parameters are then passed via the
``params`` argument of the pipeline's ``run`` invocation:

::

    temporary_dir = tempfile.mkdtemp()
    ppl.run(params={ 'ncores': 8,
                     'tmpdir': temporary_dir, })

Built-in parameters
-------------------

In addition to the custom parameters defined using the
``add_param`` method and outlined in the previous section, a
number of 'built-in' parameters are also available as properties
of the ``Pipeline`` instance, for use when building a pipeline.

Specifically these are:

* ``WORKING_DIR``: the working directory used by the pipeline
* ``BATCH_SIZE``: the batch size to be used when running jobs
  within pipeline tasks
* ``VERBOSE``: whether the pipeline is running in 'verbose'
  mode

These can be used in the same way as the custom parameters when
setting up tasks, for example:

::

    task = ExampleTask("This is an example",
                       ncores=ppl.params.ncores,
                       tmpdir=ppl.params.WORKING_DIR)

The values will be set when the pipeline's ``run`` method is
invoked.

Defining execution environment for a task: runners and envmodules
-----------------------------------------------------------------

It is possible to define the execution environment on a per-task
basis within a pipeline, by defining job runners and environment
modules.

Runners and environments can be declared in a parameterised
fashion when a pipline is created, using the ``add_runner`` and
``add_envmodules`` methods respectively of the ``Pipeline`` class.

For example:

::

    ppl = Pipeline()
    ppl.add_runner('4_cpus')
    ppl.add_envmodules('myenv')

This defines a ``runner`` called ``4_cpus`` and an environment
called ``myenv``.

The runners and environments are accessed via the ``runners``
and ``envmodules`` properties of the ``Pipeline`` instance, and
can be associated with tasks within the pipeline when they
are added via the ``add_task`` method, using the ``runner`` and
``envmodules`` keywords respectively).

For example:

::

    ppl.add_task(my_task,runner=ppl.runners['4_cpus'],...)

and

::

    ppl.add_task(my_task,envmodules=ppl.envmodules['myenv'],...)

Actual runners and environments can be assigned when the pipeline
is executed, via the ``runners`` and ``envmodules`` options of
the ``run`` method of the ``Pipeline`` instance - these are
mappings of the names defined previously to ``JobRunner`` instances,
and to lists of environment modules.

For example:

::

    ppl.run(runners={ '4_cpus': GEJobRunner('-pe smp.pe 4'), },
            envmodules={ 'myenv': 'apps/trimmomatic/0.38', },...)

If a runner is not explicitly set for a task then the pipeline's
default runner is used for that task; this defaults to a
``SimpleJobRunner`` instance but can be set explicitly via the
``default_runner`` argument of the ``Pipeline`` instance's ``run``
method.

Defining outputs from a pipeline
--------------------------------

It is possible to define outputs for a ``Pipeline`` instance in
the same way that outputs can be defined for individual tasks.

The ``add_output`` method of the ``Pipeline`` class allows an
arbitrary output to be defined, for example:

::

    ppl = Pipeline()
    ...
    ppl.add_output('final_result',result)
    ppl.run()

This can be accessed via the pipeline's ``output`` property:

::
    print("The result is '%s'" % ppl.output.result)

It is possible that pipeline outputs are defined as
``PipelineParam`` instances (for example, if a pipeline output is
taken from an output from one of its constituent tasks). By
default, on pipeline completion the outputs are "finalized" by
substituting the ``PipelineParam``s for their actual values. To
prevent this behaviour, set the ``finalize_outputs`` argument of
the pipeline's ``run`` method to ``False``. For example:

::

    ppl = Pipeline()
    ppl.add_output('final_result',PipelineParam())
    ...
    ppl.run(finalize_outputs=False)

It is recommended that outputs are defined as ``PipelineParam``
instances, to take advantage of the implicit task requirement
gathering mechanism.

PipelineCommand versus PipelineCommandWrapper
---------------------------------------------

In the first version of the pipeliner code, commands within the
``setup`` method of ``PipelineTasks`` had to be generated using
subclasses of the ``PipelineCommand`` class.

A simple example to run Fastqc::

    class Fastqc(PipelineCommand):
        def init(self,fastq,out_dir):
            self._fastq = fastq
            self._out_dir = out_dir
        def cmd(self):
            return Command("fastqc",
                           "-o",self._out_dir,
                           self._fastq)

which can then be used within a task, for example::

    class RunFastqc(PipelineTask):
        ...
        def setup(self):
            ...
            for fq in self.args.fastqs:
                self.add_cmd(Fastqc(fq,self.args.out_dir))

as an alternative to the example ``RunFastqc`` example task given
previously, which used the ``PipelineCommandWrapper`` class to
explicitly generate the Fastqc command.

Both approaches are valid. However: using ``PipelineCommand`` is
probably better suited to situations where the same command was used
in more than one distinct tasks. In cases where the command is only
used in one task, using ``PipelineCommandWrapper`` is recommended.

Advanced pipeline construction: combining pipelines
---------------------------------------------------

It possible to build larger pipelines out of smaller ones by using
the ``add_pipeline`` method of the ``Pipeline`` class, which pulls
tasks from one pipeline into another along with their dependency
relationships.

Optionally dependencies on tasks from the "master" pipeline can be
added to the imported pipeline tasks.

Parameters defined in the imported pipeline are also imported and
exposed with the same names; but it is also possible to override
them with with parameters defined in the master pipeline, or
parameters that are task outputs.

There are two specialised methods (``append_pipeline`` and
``merge_pipeline``) which wrap the ``add_pipeline`` method:

* ``append_pipeline`` takes all the tasks from one pipeline and
  adds them to the end of another, so that the appended tasks
  only run after the original tasks have completed;
* ``merge_pipeline`` takes all the tasks from one pipeline and
  adds them into another, without requiring that they wait until
  the original tasks have finished.

Appending can be used for building a pipeline out of distinct
'sections' of sub-pipelines; merging can be useful for running
multiple pipelines in parallel.
"""

######################################################################
# Imports
######################################################################

import os
import sys
import logging
import shutil
import time
import glob
import uuid
import inspect
import traceback
import string
import cloudpickle
import atexit
from collections import Iterator
try:
    # Python2
    from cStringIO import StringIO
except ImportError:
    # Python3
    from io import StringIO
from functools import reduce
from bcftbx.utils import mkdir
from bcftbx.utils import AttributeDictionary
from bcftbx.JobRunner import ResourceLock
from bcftbx.JobRunner import SimpleJobRunner
from .applications import Command
from .simple_scheduler import SimpleScheduler
from .simple_scheduler import SchedulerReporter

# Module specific logger
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

######################################################################
# Pipeline infrastructure constants
######################################################################

ALLOWED_CHARS = string.ascii_lowercase + string.digits + "._-"

# Pipeline failure modes
class PipelineFailure(object):
    IMMEDIATE = 0
    DEFERRED = 1

######################################################################
# Generic pipeline base classes
######################################################################

class BaseParam(object):
    """
    Provide base class for PipelineParam-type classes

    Implements core functionality that should be
    shared across all parameter-like classes, including
    assigning a UUID and enabling a task ID to be
    associated with the parameter.

    Provides the following attributes:

    - uuid
    - associated_task_id

    and the following methods:

    - associate_task
    """
    def __init__(self):
        """
        Base class for PipelineParam-type class
        """
        self._uuid = uuid.uuid4()
        self._associated_task_id = None
    def associate_task(self,task):
        """
        Associate a task with the parameter

        Arguments:
          task (PipelineTask): a task object to
            associate with the parameter
        """
        self._associated_task_id = task.id()
    @property
    def uuid(self):
        """
        Return the unique identifier (UUID) of the parameter
        """
        return self._uuid
    @property
    def associated_task_id(self):
        """
        Return the task ID of the associated task (or None)
        """
        return self._associated_task_id

class PipelineParam(BaseParam):
    """
    Class for passing arbitrary values between tasks

    The PipelineParam class offers a way to dynamically assign
    values for types which would otherwise be immutable (e.g.
    strings).

    A new PipelineParam is created using e.g.:

    >>> p = PipelineParam()

    In this form there will be no initial value, however
    one can be assigned using the ``value`` argument, e.g.:

    >>> p = PipelineParam(value="first")

    The associated value can be returned using the ``value``
    property:

    >>> p.value
    "first"

    and updated using the ``set`` method, e.g.:

    >>> p.set("last")
    >>> p.value
    "last"

    The return type can be explicitly specified on object
    instantiation using the ``type`` argument, for example
    to force that a string is always returned:

    >>> p = PipelineParam(type=str)
    >>> p.set(123)
    >>> p.value
    "123"

    If the ``default`` function is supplied then this will
    be used to generate a value if the stored value is
    ``None``:

    >>> p = PipelineParam(default=lambda: "default")
    >>> p.value
    "default"
    >>> p.set("assigned")
    >>> p.value
    "assigned"

    If a ``name`` is supplied then this will be stored and
    can be recovered via the ``name`` property:

    >>> p = PipelineParam(name="user_name")
    >>> p.name
    "user_name"

    The value of a parameter can be set to another parameter,
    in which case the value of the first parameter will be
    taken from the second:

    >>> first_param = PipelineParam(value="first")
    >>> first_param.value
    "first"
    >>> second_param = PipelineParam(value="second")
    >>> first_param.set(second_param)
    >>> first_param.value
    "second"

    A parameter can also be "replaced" with another parameter
    using its ``replace_with`` method:

    >>> p = PipelineParam(value="old")
    >>> p.value
    "old"
    >>> pp = PipelineParam(value="new")
    >>> p.replace_with(pp)
    >>> p.value
    "new"

    This feature allows parameters defined in one context
    to be matched with parameters in another (for example
    when tasks from one pipeline are imported into another).
    """
    def __init__(self,value=None,type=None,default=None,name=None):
        """
        Create a new PipelineParam instance

        Arguments:
          value (object): optional, initial value to assign
            to the instance
          type (function): optional, function used to convert
            the stored value when fetched via the `value`
            property
          default (function): optional, function which will
            return a default value if no explicit value is
            set (i.e. value is `None`)
          name (str): optional, name to associate with the
            instance
        """
        BaseParam.__init__(self)
        self._value = None
        self._type = type
        self._default = default
        if value is not None:
            self.set(value)
        self._name = str(name)
        self._replace_with = None
    def set(self,newvalue):
        """
        Update the value assigned to the instance

        Arguments:
          newvalue (object): new value to assign
        """
        self._value = newvalue
    def replace_with(self,p):
        """
        Set a parameter to replace this one with

        If a replacement parameter is set then the
        value will be taken from that parameter
        (ignoring all settings from this one)

        Arguments:
          p (PipelineParam): parameter to be used
            as a replacement on evaluation
        """
        self._replace_with = p
    @property
    def value(self):
        """
        Return the assigned value

        If a `default` function was specified on instance
        creation then this will be used to generate the
        value to return if the stored value is `None`.

        If a `type` function was also specified on instance
        creation then this will be used to convert the
        assigned value before it is returned.
        """
        if self._replace_with is not None:
            # Return the value of the parameter that
            # this should be replaced with
            try:
                return self._replace_with.value
            except AttributeError:
                raise TypeError("PipelineParameter cannot be "
                                "replaced with non-parameter "
                                "object")
        if self._value is None:
            # Try to return default value
            try:
                return self._default()
            except TypeError:
                return None
        else:
            # Return stored value
            value = resolve_parameter(self._value)
            try:
                return self._type(value)
            except TypeError:
                return value
    @property
    def name(self):
        """
        Return the name of the parameter (if supplied)
        """
        return self._name

class FileCollector(Iterator):
    """
    Class to return set of files based on glob pattern
    """
    def __init__(self,dirn,pattern):
        self._dirn = os.path.abspath(dirn)
        self._pattern = pattern
        self._files = None
        self._idx = None
    def __len__(self):
        if self._files is None:
            self._files = collect_files(self._dirn,self._pattern)
            self._idx = -1
        return len(self._files)
    def __next__(self):
        if self._files is None:
            self._files = collect_files(self._dirn,self._pattern)
        if self._idx is None:
            self._idx = 0
        else:
            self._idx += 1
        try:
            return self._files[self._idx]
        except IndexError:
            self._files = None
            self._idx = None
            raise StopIteration
    def next(self):
        """
        Implemented for Python2 compatibility
        """
        return self.__next__()

# Capture stdout and stderr from a function call
# Based on code from http://stackoverflow.com/a/16571630/579925
# Modified to handle both stdout and stderr (original version
# only handled stdout)
# Usage:
# >>> with Capturing() as output:
# ...     print("Hello!")
# >>> for line in output.stdout:
# ...     print("Line: %s" % line)
# >>> for line in output.stderr:
# ...     print("Err: %s" % line)
class Capturing(object):
    def __init__(self):
        self.stdout = list()
        self.stderr = list()
    def __enter__(self):
        self._stdout = sys.stdout
        self._stderr = sys.stderr
        sys.stdout = self._stringio = StringIO()
        sys.stderr = self._stringerr = StringIO()
        return self
    def __exit__(self,*args):
        self.stdout.extend(self._stringio.getvalue().splitlines())
        self.stderr.extend(self._stringerr.getvalue().splitlines())
        del self._stringio    # free up some memory
        del self._stringerr
        sys.stdout = self._stdout
        sys.stderr = self._stderr

class Pipeline(object):
    """
    Class to define and run a 'pipeline' of 'tasks'

    A pipeline consists of a set of tasks (defined by
    instantiating subclasses of PipelineTask) with
    simple dependency relationships (i.e. a task will
    depend on none, one or more other tasks in the
    pipeline to complete before it can be executed).

    Example usage:

    >> p = Pipeline()
    >> t1 = p.add_task(Task1())
    >> t2 = p.add_task(Task2(),requires=(t1,))
    >> ...
    >> p.run()

    Tasks will only run when all requirements have
    completed (or will run immediately if they don't
    have any requirements).
    """
    def __init__(self,name="PIPELINE"):
        """
        Create a new Pipeline instance

        Arguments:
          name (str): optional name for the pipeline
        """
        self._name = str(name)
        self._id = "%s.%s" % (sanitize_name(self._name),
                              time.strftime("%Y%m%d.%H%M%S"))
        self._tasks = {}
        self._pending = []
        self._running = []
        self._finished = []
        self._failed = []
        self._removed = []
        self._params = AttributeDictionary()
        self._runners = dict()
        self._envmodules = dict()
        self._scheduler = None
        self._log_file = None
        self._lock_manager = None
        self._exit_on_failure = PipelineFailure.IMMEDIATE
        # Initialise default runner
        self._runners['default'] = PipelineParam(value=SimpleJobRunner())
        # Initialise built-in parameters
        self.add_param("WORKING_DIR",type=str)
        self.add_param("BATCH_SIZE",type=int)
        self.add_param("VERBOSE",type=bool)
        # Initialise pipeline output
        self._output = AttributeDictionary()

    @property
    def name(self):
        """
        Return the name of the pipeline
        """
        return self._name

    @property
    def params(self):
        """
        Access the parameters defined for the pipeline
        """
        return self._params

    @property
    def runners(self):
        """
        Access the runners defined for the pipeline

        Returns the dictionary mapping runner names to
        PipelineParam instances that store the runner
        instances; so to get the runner associated with
        a name do e.g.

        >>> runner = ppl.runners['my_runner'].value
        """
        return dict(**self._runners)

    @property
    def envmodules(self):
        """
        Access the modules environments defined for the pipeline

        Returns the dictionary mapping modules environment
        names to PipelineParam instances that store the
        module lists; so to get the list of modules associated
        with an environment name do e.g.

        >>> modules = ppl.envmodules['my_env'].value
        """
        return dict(**self._envmodules)

    @property
    def output(self):
        """
        Return the output object
        """
        return self._output

    def report(self,s):
        """
        Internal: report messages from the pipeline
        """
        report = []
        for line in s.split('\n'):
            report.append("%s [%s] %s" %(time.strftime("%Y-%m-%d %H:%M:%S"),
                                         self._name,line))
        report = '\n'.join(report)
        if self._log_file:
            with open(self._log_file,'a') as log:
                log.write("%s\n" % (report,))
        print(report)

    def terminate(self):
        """
        Internal: terminate a running pipeline
        """
        self.report("Terminating pipeline: failed")
        # Dump all pending tasks
        self._pending = []
        # Stop all running tasks
        if self._running:
            self.report("Stopping running tasks:")
            for task in self._running:
                self.report("- %s" % task.id())
                task.terminate()
        # Stop the scheduler
        self.stop_scheduler()
        return 1

    def start_scheduler(self,runner=None,max_concurrent=1,
                        max_slots=None,poll_interval=5):
        """
        Internal: instantiate and start local scheduler
        """
        if self._scheduler is None:
            sched = SimpleScheduler(runner=runner,
                                    max_concurrent=max_concurrent,
                                    max_slots=max_slots,
                                    poll_interval=poll_interval,
                                    reporter=SchedulerReporter())
            sched.start()
            self._scheduler = sched
        return self._scheduler

    def stop_scheduler(self):
        """
        Internal: stop the pipeline scheduler
        """
        if self._scheduler is not None:
            self._scheduler.stop()
            self._scheduler = None

    def add_task(self,task,requires=(),**kws):
        """
        Add a task to the pipeline

        Arguments:
          task (PipelineTask): task instance to add to
            the pipeline
          requires (List): list or tuple of task instances
            which need to complete before this task will
            start
          kws (Dictionary): a dictionary of keyword-value
            pairs which will be passed to the task at
            run time (see the ``run`` method of
            PipelineTask for valid options)
        """
        self._tasks[task.id()] = (task,kws)
        for req in requires:
            task.requires(req)
            if req.id() not in self.task_list():
                self.add_task(req,())
        return task

    def add_param(self,name,value=None,type=None):
        """
        Define a new pipeline parameter

        Creates a new ``PipelineParam`` instance associated
        with the supplied name.

        Parameters can be accessed and set via the ``params``
        property of the pipeline.

        Arguments:
          name (str): name for the new parameter
          value (object): optional, initial value to assign
            to the instance
          type (function): optional, function used to convert
            the stored value when fetched via the `value`
            property
        """
        if name in self._params:
            raise KeyError("Parameter '%s' already defined"
                           % name)
        self._params[name] = PipelineParam(type=type,
                                           value=value)

    def add_runner(self,name):
        """
        Define a new runner within the pipeline

        Creates a new ``PipelineParam`` instance associated
        with the supplied runner name.

        Runner instances can be accessed and set via the
        ``runners`` property of the pipeline, for example:

        To access:

        >>> runner = ppl.runners['my_runner'].value

        To set:

        >>> ppl.runners['my_runner'].set(SimpleJobRunner())

        Arguments:
          name (str): name for the new runner
        """
        if name in self._runners:
            raise KeyError("Runner '%s' already defined" % name)
        self.report("Defining new runner '%s'" % name)
        self._runners[name] = PipelineParam(
            name=name,
            default=lambda: self.runners['default'].value)

    def add_envmodules(self,name):
        """
        Define a new environment defined by modules

        Creates a new ``PipelineParam`` instance associated
        with the supplied environment name.

        Runner instances can be accessed and set via the
        ``envmodules`` property of the pipeline, for example:

        To access:

        >>> env_modules = ppl.envmodules['my_env'].value

        To set:

        >>> ppl.envmodules['my_env'].set("")

        Arguments:
          name (str): name for the new runner
        """
        if name in self._envmodules:
            raise KeyError("Modules environment '%s' already defined"
                           % name)
        self.report("Defining new modules environment '%s'" % name)
        self._envmodules[name] = PipelineParam(
            name=name,
            value=list())

    def add_output(self,name,value):
        """
        Add an output to the pipeline

        Arguments:
          name (str): name for the output
          value (object): associated object
        """
        self._output[name] = value

    def add_pipeline(self,pipeline,params=None,requires=None):
        """
        Import tasks from another pipeline

        Adds the tasks from the supplied pipeline
        instance into this pipeline.

        Dependency relationships defined between tasks
        in the imported pipeline are preserved on
        import. The 'requires' keyword can be used to
        define new dependencies between the initial
        tasks in the imported pipeline and those in
        the destination pipeline.

        By default parameters defined in the added
        pipeline will be imported and exposed with the
        same names; the 'params' keyword can be used to
        replace them with arbitrary parameters (e.g.
        parameters defined in the master pipeline,
        outputs from tasks etc).

        Arguments:
          pipeline (Pipeline): pipeline instance with
            tasks to be added to this pipeline
          params (mapping): a dictionary or mapping
            where keys are parameter names in the
            imported pipeline and values are
            PipelineParam instances that they will
            be replaced by
          requires (list): optional list of tasks
            that the added pipeline will depend on
        """
        self.report("Importing tasks from pipeline '%s'" % pipeline.name)
        # Store direct references to imported parameters
        self.report("Importing and masking parameters")
        imported_params = pipeline.params
        # Import parameters into the pipeline
        for p in imported_params:
            if p not in self._params:
                # Imported parameter name doesn't have an equivalent
                # in this pipeline, so add a direct reference to it
                self._params[p] = imported_params[p]
                self.report("- imported parameter '%s'" % p)
            else:
                # Imported parameter name matches an existing parameter
                # Update the value
                imported_params[p].set(self.params[p].value)
            # Additionally we may need to "mask" the parameters from
            # the imported pipeline (by assigning replacement
            # parameters to them)
            imported_param = imported_params[p]
            if params and p in params:
                # Explicit replacement
                new_param = params[p]
            else:
                # Implicit replacement
                new_param = self.params[p]
            # Check that the imported parameter isn't going to
            # be replaced by itself
            try:
                if imported_param.uuid != new_param.uuid:
                    imported_param.replace_with(new_param)
                    self.report("- replaced parameter '%s'" % p)
            except Exception as ex:
                # No UUID for replacement: not a parameter?
                self.report("- imported parameter is not a "
                            "PipelineParam? (%s)" % ex)
                continue
        # Add runner definitions from the new pipeline
        for r in pipeline.runners:
            if r not in self.runners:
                self.add_runner(r)
        # Get the starting tasks from the new pipeline
        initial_tasks = pipeline.initial_tasks
        self.report("Identified initial tasks from imported pipeline:")
        for t in initial_tasks:
            self.report("- %s" % t.name())
        # Add the tasks from the new pipeline
        self.report("Importing tasks")
        imported_tasks = pipeline.task_list()
        for task_id in imported_tasks:
            task,requirements,kws = pipeline.get_task(task_id)
            # Update the runner for the task
            if 'runner' in kws:
                name = kws['runner'].name
                kws['runner'] = self.runners[name]
                self.report("- updated runner '%s' for imported task '%s'"
                            % (name,task.name()))
            # Update the environment for the task
            if 'envmod' in kws:
                name = kws['envmod'].name
                kws['envmod'] = self.envmod[name]
                self.report("- updating envmodule '%s' for imported "
                            "task '%s'" % (name,task.name()))
            # Add the task to the source pipeline
            self.add_task(task,requires=requirements,**kws)
        # Update the requirements on the initial tasks from
        # the imported pipeline
        if requires:
            self.report("Updating dependencies on imported tasks")
            for task in initial_tasks:
                task.requires(*requires)

    def append_pipeline(self,pipeline):
        """
        Append tasks from another pipeline

        Adds the tasks from the supplied pipeline
        instance into this pipeline, ensuring that the
        appended tasks depend on the original tasks
        completing.

        Dependencies which were already defined in the
        pipeline being appended are preserved when they
        are added to the first.

        Arguments:
          pipeline (Pipeline): pipeline instance with
            tasks to be appended
        """
        self.report("Appending tasks from pipeline '%s'" % pipeline.name)
        # Get the final tasks from the base pipeline
        final_tasks = self.final_tasks
        self.report("Identified final tasks from base pipeline:")
        for task in final_tasks:
            self.report("-- %s" % task.name())
        # Add the pipeline
        self.add_pipeline(pipeline,requires=final_tasks)

    def merge_pipeline(self,pipeline):
        """
        Add tasks from another pipeline

        Adds the tasks from the supplied pipeline
        instance into this pipeline, without requiring
        that the appended tasks depend on the original
        tasks.

        Dependencies which were already defined in the
        pipeline being appended are preserved when they
        are added to the first.

        Arguments:
          pipeline (Pipeline): pipeline instance with
            tasks to be added
        """
        self.add_pipeline(pipeline)

    def task_list(self):
        """
        Return a list of task ids
        """
        return list(self._tasks.keys())

    def get_task(self,task_id):
        """
        Return information on a task

        Arguments:
          task_id (str): unique identifier for
            a task

        Returns:
          Tuple: tuple of (task,requirements,kws)
            where 'task' is a PipelineTask instance,
            requirements is a list of PipelineTask
            instances that the task depends on,
            and kws is a keyword mapping
        """
        # Fetch task object and keywords
        task,kws = self._tasks[task_id]
        # Fetch dependent tasks
        requirements = tuple([self.get_task(t)[0]
                              for t in task.required_task_ids])
        return (task,requirements,kws)

    def rank_tasks(self):
        """
        Rank the tasks into order

        Returns:
          List: list of 'ranks', with each rank
            being a list of task ids.
        """
        # Create a dictionary with task ids as keys and
        # lists of the required task ids as values
        required = dict()
        for task_id in self.task_list():
            required[task_id] = [t.id()
                                 for t in self.get_task(task_id)[1]]
        # Rank the task ids
        ranks = list()
        while required:
            # Populate the current rank with tasks which don't
            # have any requirements
            current_rank = []
            for task_id in list(required.keys()):
                if not required[task_id]:
                    # Task no longer has requirements so add
                    # to the current rank
                    current_rank.append(task_id)
                    # Remove task id from the list of tasks so it
                    # isn't checked on the next iteration
                    del(required[task_id])
            # Update the requirement lists for remaining tasks
            # to remove those add to the current rank
            for task_id in list(required.keys()):
                new_reqs = []
                for req in required[task_id]:
                    if req not in current_rank:
                        new_reqs.append(req)
                required[task_id] = new_reqs
            ranks.append(current_rank)
        return ranks

    @property
    def initial_tasks(self):
        """
        Fetch a list of 'initial tasks' for the pipeline

        Returns:
          List: list of task instances from the pipeline
            which don't depend on any other tasks (i.e.
            the tasks that will run first)
        """
        ranks = self.rank_tasks()
        return [self.get_task(t)[0] for t in ranks[0]]

    @property
    def final_tasks(self):
        """
        Fetch a list of 'final tasks' for the pipeline

        Returns:
          List: list of task instances from the pipeline
            which don't have any dependent tasks (i.e.
            the tasks that will run last)
        """
        final_tasks = []
        for task_id in self.task_list():
            if not self.get_dependent_tasks(task_id):
                final_tasks.append(self.get_task(task_id)[0])
        return final_tasks

    def get_dependent_tasks(self,task_id):
        """
        Return task ids that depend on supplied task id
        """
        dependents = set()
        for id_ in self.task_list():
            task,requires,kws = self.get_task(id_)
            for req in requires:
                if req.id() == task_id:
                    dependents.update(self.get_dependent_tasks(id_))
                    dependents.add(id_)
                    break
        return list(dependents)

    def run(self,working_dir=None,log_dir=None,scripts_dir=None,
            log_file=None,sched=None,default_runner=None,max_jobs=1,
            max_slots=None,poll_interval=5,params=None,runners=None,
            envmodules=None,batch_size=None,verbose=False,
            use_locking=True,exit_on_failure=PipelineFailure.IMMEDIATE,
            finalize_outputs=True):
        """
        Run the tasks in the pipeline

        Arguments:
          working_dir (str): optional path to a working
            directory (defaults to the current directory)
          log_dir (str): path of directory where log files
            will be written to
          scripts_dir (str): path of directory where script
            files will be written to
          log_file (str): path to file to write pipeline
            log messages to (in addition to stdout)
          sched (SimpleScheduler): a scheduler to use for
            running commands generated by each task
          default_runner (JobRunner): optional default
            job runner to use
          max_jobs (int): optional maximum number of
            concurrent jobs in scheduler (defaults to 1;
            ignored if a scheduler is provided via 'sched'
            argument)
          max_slots (int): optional maximum number of 'slots'
            (i.e. concurrent threads or maximum number of
            CPUs) available to the scheduler (defaults to
            no limit)
          poll_interval (float): optional polling interval
            (seconds) to set in scheduler (if scheduler not
            provided via the 'sched' argument), and to use
            for checking if tasks have completed (defaults
            to 5s)
          params (mapping): a dictionary or mapping which
            associates parameter names with values
          runners (mapping): a dictionary or mapping which
            associates runner names with job runners
          envmodules (mapping): a dictionary or mapping which
            associates envmodule names with list of
            environment module names
          batch_size (int): if set then run commands in
            each task in batches, with each batch running
            this many commands at a time (default is to run
            one command per job)
          verbose (bool): if True then report additional
            information for diagnostics
          use_locking (book): if True then lock invocations
            of task methods
          exit_on_failure (int): either IMMEDIATE (any
            task failures cause immediate termination of
            of the pipeline; this is the default) or
            DEFERRED (the pipeline execution continues
            and only raises an error when all tasks
            have finished running)
          finalize_outputs (bool): if True then convert any
            pipeline outputs from PipelineParams to an
            actual value (this is the default)
        """
        # Deal with working directory
        if working_dir is None:
            working_dir = os.getcwd()
        working_dir = os.path.abspath(working_dir)
        # Deal with custom parameters
        if params:
            for p in params:
                self.params[p].set(params[p])
        # Deal with runners
        if runners:
            for r in runners:
                if r in self.runners:
                    self.runners[r].set(runners[r])
                else:
                    raise Exception("Undefined runner '%s'" % r)
        if default_runner:
            self.runners['default'].set(default_runner)
        # Deal with environment modules
        if envmodules:
            for name in envmodules:
                if name in self.envmodules:
                   if envmodules[name] is not None:
                       for m in envmodules[name].split(','):
                           self.envmodules[name].value.append(m)
                else:
                    self.report("WARNING ignoring undefined environment "
                                "module '%s'" % name)
        # Deal with scheduler
        if sched is None:
            # Create and start a scheduler
            sched = self.start_scheduler(runner=default_runner,
                                         max_concurrent=max_jobs,
                                         max_slots=max_slots,
                                         poll_interval=poll_interval)
        # Deal with lock manager
        if use_locking:
            self._lock_manager = ResourceLock()
        # Deal with log directory
        if log_dir is None:
            log_dir = "%s.logs" % self._id
        if not os.path.isabs(log_dir):
            log_dir = os.path.join(working_dir,log_dir)
        if not os.path.exists(log_dir):
            os.mkdir(log_dir)
        # Deal with scripts directory
        if scripts_dir is None:
            scripts_dir = "%s.scripts" % self._id
        if not os.path.isabs(scripts_dir):
            scripts_dir = os.path.join(working_dir,scripts_dir)
        if not os.path.exists(scripts_dir):
            os.mkdir(scripts_dir)
        # Deal with log file
        if log_file:
            self._log_file = os.path.abspath(log_file)
            if os.path.exists(self._log_file):
                os.remove(self._log_file)
        # How to handle task failure
        self._exit_on_failure = exit_on_failure
        # Assign built-in parameters
        self.params.WORKING_DIR.set(working_dir)
        self.params.BATCH_SIZE.set(batch_size)
        self.params.VERBOSE.set(verbose)
        # Execute the pipeline
        self.report("Started")
        self.report("-- working directory: %s" % working_dir)
        self.report("-- log directory    : %s" % log_dir)
        self.report("-- scripts directory: %s" % scripts_dir)
        self.report("-- log file         : %s" % ('<not set>'
                                                  if self._log_file is None
                                                  else self._log_file))
        self.report("-- batch size       : %s" % ('<not set>'
                                                  if batch_size is None
                                                  else batch_size))
        self.report("-- verbose output   : %s" % ('yes'
                                                  if verbose
                                                  else 'no'))
        self.report("-- concurrent jobs  : %s" % (max_jobs
                                                  if max_jobs else
                                                  'unlimited'))
        self.report("-- maximum slots    : %s" % (max_slots
                                                  if max_slots else
                                                  'unlimited'))
        self.report("-- polling interval : %ss" % ('<not set>'
                                                   if poll_interval is None
                                                   else poll_interval))
        # Report parameter settings
        if self.params:
            self.report("Pipeline parameters:")
            width = max([len(p) for p in self.params])
            for p in sorted([p for p in self.params]):
                self.report("-- %s%s: %s" % (p,
                                             ' '*(width-len(p)),
                                             ('<not set>'
                                              if self.params[p].value is None
                                              else self.params[p].value)))
        # Report runners
        self.report("Runners:")
        width = max([len(r) for r in self.runners])
        for r in sorted(self.runners):
            runner = self.runners[r].value
            self.report("-- %s%s: %s%s" % (r,
                                         ' '*(width-len(r)),
                                           runner,
                                           '' if runner.nslots == 1
                                           else ' [nslots=%s]' %
                                           runner.nslots))
        # Report modules environments
        if self.envmodules:
            self.report("Modules environments:")
            width = max([len(m) for m in self.envmodules])
            for m in sorted(self.envmodules):
                self.report("-- %s%s: %s" %
                            (m,
                             ' '*(width-len(m)),
                             ('<not set>'
                              if self.envmodules[m].value is None
                              else
                              ','.join(
                                  [str(x) for x in self.envmodules[m].value]))))
        # Sort the tasks and set up the pipeline
        self.report("Scheduling tasks...")
        for i,rank in enumerate(self.rank_tasks()):
            self.report("Task rank %d:" % i)
            for task_id in rank:
                task,requires,kws = self.get_task(task_id)
                self._pending.append((task,requires,kws))
                if verbose:
                    self.report("-- %s (%s)" % (task.name(),
                                                task.id()))
                else:
                    self.report("-- %s" % task.name())
        # Run while there are still pending or running tasks
        update = True
        while self._pending or self._running:
            pending = []
            running = []
            failed = []
            # Check for pending tasks that can start
            for task,requirements,kws in self._pending:
                run_task = False
                if not requirements:
                    # No requirements - start it
                    run_task = True
                else:
                    # Check requirements
                    run_task = reduce(lambda x,y: x and
                                      (y not in self._running) and
                                      y.completed and
                                      (y.exit_code == 0),
                                      requirements,True)
                if run_task:
                    if verbose:
                        self.report("started '%s' (%s)" % (task.name(),
                                                           task.id()))
                    else:
                        self.report("started '%s'" % task.name())
                    kws = dict(**kws)
                    if 'runner' not in kws:
                        kws['runner'] = default_runner
                    if 'working_dir' not in kws:
                        kws['working_dir'] = working_dir
                    if 'log_dir' not in kws:
                        kws['log_dir'] = log_dir
                    if 'scripts_dir' not in kws:
                        kws['scripts_dir'] = scripts_dir
                    if 'batch_size' not in kws:
                        kws['batch_size'] = batch_size
                    if 'verbose' not in kws:
                        kws['verbose'] = verbose
                    kws['log_file'] = self._log_file
                    kws['lock_manager'] = self._lock_manager
                    for k in kws:
                        # Resolve keywords which are PipelineParams etc
                        # and replace with the final value
                        kws[k] = resolve_parameter(kws[k])
                    try:
                        task.run(sched=sched,
                                 poll_interval=poll_interval,
                                 **kws)
                    except Exception as ex:
                        # Exception trying to run the task
                        task.fail(exit_code=1,
                                  message="Exception trying to run task "
                                  "'%s': %s\n%s" %
                                  (task.id(),ex,traceback.format_exc()))
                        failed.append(task)
                    self._running.append(task)
                    update = True
                else:
                    pending.append((task,requirements,kws))
            self._pending = pending
            # Check for running tasks that have completed
            for task in self._running:
                if task.completed:
                    self.report("finished '%s'"
                                % task.name())
                    self._finished.append(task)
                    update = True
                    # Check if task failed
                    if task.exit_code != 0:
                        failed.append(task)
                else:
                    # Still running
                    running.append(task)
                    # Report if status has changed
                    if task.updated:
                        njobs,ncompleted = task.njobs()
                        if njobs > 1:
                            self.report("'%s' updated (%d/%d)" %
                                        (task.name(),
                                         ncompleted,
                                         njobs))
            self._running = running
            # Check for tasks that have failed
            if failed:
                if self._exit_on_failure == PipelineFailure.IMMEDIATE:
                    # Terminate all tasks and stop immediately
                    self.report("Following tasks failed:")
                    for task in failed:
                        msg = "- '%s'" % task.name
                        if verbose:
                            msg += " (%s)" % task.id()
                        self.report(msg)
                    self._failed.extend(failed)
                    return self.terminate()
                elif self._exit_on_failure == PipelineFailure.DEFERRED:
                    # Remove any tasks waiting on the failures
                    # but defer pipeline termination
                    for task in failed:
                        self.report("Task failed: '%s' (%s)" %
                                    (task.name(),task.id()))
                        if verbose and task.stdout:
                            self.report("Task stdout:\n%s" %
                                        '\n>> '.join(
                                            task.stdout.split('\n')))
                        dependent_tasks = self.get_dependent_tasks(
                            task.id())
                        pending = []
                        for t in self._pending:
                            if t[0].id() in dependent_tasks:
                                msg = "-- removing dependent task '%s'" \
                                      % t[0].name()
                                if verbose:
                                    msg += " (%s)" % t[0].id()
                                self.report(msg)
                                self._removed.append(t[0])
                            else:
                                pending.append(t)
                        self._pending = pending
                    self._failed.extend(failed)
                    self.report("There are failed tasks but pipeline exit "
                                "will be deferred")
            # Report the current running and pending tasks
            if update:
                if self._running:
                    self.report("%d running tasks:"
                                % len(self._running))
                    for t in self._running:
                        njobs,ncompleted = t.njobs()
                        if njobs > 1:
                            self.report("- %s (%d/%d)" % (t.name(),
                                                          ncompleted,
                                                          njobs))
                        else:
                            self.report("- %s" % t.name())
                if self._pending:
                    self.report("%d pending tasks:"
                                % len(self._pending))
                    for t in self._pending:
                        self.report("- %s" % t[0].name())
                update = False
            else:
                # Pause before checking again
                time.sleep(poll_interval)
        # Finished
        if finalize_outputs:
            # Finalize the outputs
            self.report("Finalizing outputs")
            for name in self._output:
                try:
                    self._output[name] = self._output[name].value
                except AttributeError:
                    pass
                self.report("- setting '%s': %s" % (name,
                                                    self._output[name]))
        if self._failed:
            # Report failed tasks
            self.report("Pipeline completed but the following tasks failed:")
            for task in self._failed:
                msg = "- '%s'" % task.name()
                if verbose:
                    msg += " (%s)" % task.id()
                self.report(msg)
            if self._removed:
                # Report any tasks that were removed
                self.report("The following tasks were not executed:")
                for task in self._removed:
                    msg = "- '%s'" % task.name()
                    if verbose:
                        msg += " (%s)" % task.id()
                    self.report(msg)
            self.report("Completed: failed")
            return 1
        else:
            # No errors
            self.report("Completed: ok")
            return 0

class PipelineTask(object):
    """
    Base class defining a 'task' to run as part of a pipeline

    A 'task' wraps one or more external programs which can
    be run concurrently, and which produces a set of outputs.
    Individual programs should be wrapped in instances of the
    'PipelineCommand' class.

    This class should be subclassed to implement the 'init',
    'setup', 'finish' (optionally) and 'output' methods.

    The 'add_cmd' method can be used within 'setup' to add one
    or more 'PipelineCommand' instances.

    """
    def __init__(self,_name,*args,**kws):
        """
        Create a new PipelineTask instance

        Arguments:
          _name (str): an arbitrary user-friendly name for the
            task instance
          args (List): list of arguments to be supplied to
            the subclass (must match those defined in the
            'init' method)
          kws (Dictionary): dictionary of keyword-value pairs
            to be supplied to the subclass (must match those
            defined in the 'init' method)
        """
        self._name = str(_name)
        self._args = args
        self._kws = kws
        self._commands = []
        self._scripts = []
        self._required_task_ids = []
        self._task_name = "%s.%s" % (sanitize_name(self._name),
                                     uuid.uuid4())
        self._completed = False
        self._stdout_files = []
        self._exit_code = 0
        # Working directory
        self._working_dir = None
        # Running jobs
        self._runner_nslots = None
        self._jobs = []
        self._groups = []
        # Monitoring
        self._ncompleted = 0
        # Logging
        self._log_file = None
        # Output
        self._output = AttributeDictionary()
        # Locking
        self._lock_manager = None
        # Deal with subclass arguments
        try:
            self._callargs = inspect.getcallargs(self.init,*args,**kws)
        except Exception as ex:
            logger.error("Exception setting up args for task '%s' (%s): %s"
                         % (self._name,self.__class__,ex))
            raise ex
        try:
            del(self._callargs['self'])
        except KeyError:
            pass
        # Check for parameter arguments & keywords that have associated
        # tasks and add these as requirements
        for arg in self._args:
            try:
                task_id = arg.associated_task_id
                if task_id:
                    self.requires_id(task_id)
            except AttributeError as ex:
                pass
        for kw in self._kws:
            try:
                task_id = self._kws[kw].associated_task_id
                if task_id:
                    self.requires_id(task_id)
            except AttributeError:
                pass
        # Execute the init method
        self.invoke(self.init,self._args,self._kws)
        if self._exit_code != 0:
            raise Exception("Error running 'init' method for task '%s' "
                            "(%s)" % (self._name,self.__class__))

    @property
    def args(self):
        """
        Fetch parameters supplied to the instance
        """
        args = AttributeDictionary(**self._callargs)
        for a in args:
            args[a] = resolve_parameter(args[a])
        return args

    @property
    def completed(self):
        """
        Check if the task has completed
        """
        return self._completed

    @property
    def updated(self):
        """
        Check if the task has been updated since the last check
        """
        njobs,ncompleted = self.njobs()
        updated = ncompleted > self._ncompleted
        self._ncompleted = ncompleted
        return updated

    @property
    def exit_code(self):
        """
        Get the exit code for completed task

        Returns:
          Integer: exit code, or 'None' if task hasn't completed
        """
        if not self.completed:
            return None
        else:
            return self._exit_code

    @property
    def stdout(self):
        """
        Get the standard output from the task

        Returns:
          String: standard output from the task.
        """
        stdout = []
        for f in self._stdout_files:
            with open(f,'r') as fp:
                stdout.append(fp.read())
        return ''.join(stdout)

    @property
    def runner_nslots(self):
        """
        Get number of slots (i.e. available CPUs) set by runner

        This returns the number of CPUs available to the command
        as set in the job runner which will be used to run the
        job.

        Returns:
          String: environment variable to get slots from.
        """
        if self._runner_nslots:
            return self._runner_nslots
        else:
            return 1

    def name(self):
        """
        Get the name of the task

        Returns:
          String: the name supplied when the task was created
        """
        return self._name

    def id(self):
        """
        Get the name of the task within the pipeline

        Returns:
          String: a name consisting of a 'sanitized' version
            of the supplied name appended with a unique id
            code
        """
        return self._task_name

    def fail(self,exit_code=1,message=None):
        """
        Register the task as failing

        Intended to be invoked from the subclassed 'setup'
        or 'finish' methods, to terminate the task and
        indicate that it has failed.

        NB when using the 'fail' method it is recommended that
        the method that it was invoked from should return
        immediately afterwards, to avoid any unexpected side
        effects

        Arguments:
          exit_code (int): optional, specifies the exit code
            to return (defaults to 1)
          message (str): optional, error message to report to
            the pipeline user
        """
        if message:
            self.report("failed: %s" % message)
        self.report("failed: exit code set to %s" % exit_code)
        self._exit_code = exit_code
        self._completed = True

    def report(self,s):
        """
        Internal: report messages from the task
        """
        report = []
        for line in s.split('\n'):
            report.append("%s [Task: %s] %s" %(time.strftime("%Y-%m-%d %H:%M:%S"),
                                               self._name,line))
        report = '\n'.join(report)
        if self._log_file:
            with open(self._log_file,'a') as log:
                log.write("%s\n" % (report,))
        print(report)

    def invoke(self,f,args=None,kws=None):
        """
        Internal: invoke arbitrary method on the task

        Arguments:
          f (function): method to invoke (e.g. 'self.init')
          args (list): arguments to invoke function with
          kws (dictionary): keyworded parameters to invoke
            function with
        """
        # Get the lock (if lock manager exists)
        if self._lock_manager:
            lock = self._lock_manager.acquire("task.invoke")
        else:
            lock = None
        # Switch to working directory, if defined
        if self._working_dir is not None:
            current_dir = os.getcwd()
            os.chdir(self._working_dir)
        # Invoke the requested method
        caught_exception = None
        try:
            with Capturing() as output:
                if args is None:
                    f()
                else:
                    f(*args,**kws)
        except NotImplementedError:
            pass
        except Exception as ex:
            # Store exception and traceback
            caught_exception = ex
            traceback_ = traceback.format_exc()
            self._exit_code += 1
        # Report stdout and stderr
        for line in output.stdout:
            self.report("[%s] %s" % (f.__name__,line))
        for line in output.stderr:
            self.report("[%s] %s" % (f.__name__,line))
        # Report exception and diagnostics
        if caught_exception:
            self.report("exception invoking '%s': %s" %
                        (f.__name__,caught_exception))
            self.report("\n%s" % traceback_)
            self.report_diagnostics(str(caught_exception))
        # Switch back to original directory
        if self._working_dir is not None:
            os.chdir(current_dir)
        # Release the lock
        if lock:
            self._lock_manager.release(lock)

    def report_diagnostics(self,s):
        """
        Internal: report additional diagnostic information

        Reports current and working directories, current
        directory contents, scripts and script outputs;
        to be invoked on task failure.

        Arguments:
          s (str): string describing the reason for the
            diagnostics being reported
        """
        # Report diagnostic information
        self.report("**** TASK DIAGNOSTICS ****\n")
        self.report("Reason: '%s'\n" % s)
        self.report("Working directory %s" % self._working_dir)
        self.report("CWD %s" % os.getcwd())
        self.report("\nCWD contents:")
        contents = sorted(os.listdir(os.getcwd()))
        if contents:
            for item in contents:
                self.report("-- %s%s" % (item,
                                         os.sep if os.path.isdir(
                                             os.path.join(os.getcwd(),item))
                                         else ''))
        else:
            self.report("Empty")
        self.report("\nSCRIPTS:")
        if self._scripts:
            for script_file in self._scripts:
                with open(script_file,'rt') as fp:
                    self.report("%s:" % script_file)
                    for line in fp:
                        self.report("> %s" % line.rstrip('\n'))
        else:
            self.report("No scripts generated for this task")
        self.report("\nSTDOUT:")
        if self.stdout:
            for line in self.stdout.split('\n'):
                self.report("> %s" % line)
        else:
            self.report("No stdout from task scripts")
        self.report("\n**** END OF DIAGNOSTICS ****")

    def task_completed(self,name,jobs,sched):
        """
        Internal: callback method

        This is a callback method which is invoked when
        scheduled jobs in the task finish

        Arguments:
          name (str): name for the callback
          jobs (list): list of SchedulerJob instances
          sched (SimpleScheduler): scheduler instance
        """
        for job in jobs:
            try:
                if job.exit_code != 0:
                    self._exit_code += 1
                self._stdout_files.append(job.log)
            except AttributeError:
                # Assume it's a group
                for j in job.jobs:
                    if j.exit_code != 0:
                        self._exit_code += 1
                    self._stdout_files.append(j.log)
        self.finish_task()

    def finish_task(self):
        """
        Internal: perform actions to finish the task
        """
        if self._exit_code != 0:
            logger.critical("'%s' failed: exit code %s"
                            % (self.name(),self._exit_code))
        else:
            # Execute 'finish', if implemented
            self.invoke(self.finish)
        # Flag job as completed
        self._completed = True
        # Report completion
        njobs,ncompleted = self.njobs()
        status = ('ok' if self._exit_code == 0 else 'failed')
        if njobs > 1:
            self.report("completed (%d/%d): %s" % (ncompleted,
                                                   njobs,
                                                   status))
        else:
            self.report("completed: %s" % status)

    def add_cmd(self,pipeline_job):
        """
        Add a PipelineCommand to the task

        Arguments:
           pipeline_job (PipelineCommand): a PipelineCommand
             instance to be executed by the task when it
             runs
        """
        self._commands.append(pipeline_job)

    def requires(self,*tasks):
        """
        Add tasks as requirements

        Each specified task will be added to the list
        of tasks that need to complete before this
        task can run.

        Arguments:
          tasks (List): list of PipelineTask objects
            to be added to this task as requirements
        """
        for t in tasks:
            # Fetch associated ID
            try:
                task_id = t.id()
            except AttributeError:
                raise Exception("%s: not a task?" % t)
            # Add ID
            self.requires_id(task_id)

    def requires_id(self,task_id):
        """
        Add task ID to list of requirements

        Arguments:
          task_id (str): UUID of the task to be added
            as a requirement
        """
        if task_id not in self._required_task_ids:
            self._required_task_ids.append(task_id)
        self._required_task_ids = sorted(self._required_task_ids)

    @property
    def required_task_ids(self):
        """
        Return the task IDs for tasks required by this task
        """
        return self._required_task_ids

    def add_output(self,name,value):
        """
        Add an output to the task

        Arguments:
          name (str): name for the output
          value (object): associated object
        """
        self._output[name] = value
        try:
            self._output[name].associate_task(self)
        except AttributeError as ex:
            pass

    @property
    def output(self):
        """
        Return the output object
        """
        return self._output

    def run(self,sched=None,runner=None,envmodules=None,working_dir=None,
            log_dir=None,scripts_dir=None,log_file=None,wait_for=(),
            asynchronous=True,poll_interval=5,batch_size=None,
            lock_manager=None,verbose=False):
        """
        Run the task

        This method is not normally invoked directly; instead
        it's called by the pipeline that the task has been
        added to.

        Arguments:
          sched (SimpleScheduler): scheduler to submit jobs to
          runner (JobRunner): job runner to use when running
            jobs via the scheduler
          envmodules (list): list of environment modules to load when
            running jobs in the task
          working_dir (str): path to the working directory to use
            (defaults to the current working directory)
          log_dir (str): path to the directory to write logs to
            (defaults to the working directory)
          scripts_dir (str): path to the directory to write
            scripts to (defaults to the working directory)
          log_file (str): path to file to write task log messages
            to (in addition to stdout)
          wait_for (list): deprecated: list of scheduler jobs to
            wait for before running jobs from this task
          asynchronous (bool): if False then block until the task
            has completed
          poll_interval (float): interval between checks on task
            completion (in seconds) for non-asynchronous tasks
            (defaults to 5 seconds)
          batch_size (int): if set then run commands in
            each task in batches, with each batch running
            this many commands at a time (default is to run
            one command per job)
          lock_manager (ResourceLock): inter-task lock manager
          verbose (bool): if True then report additional
            information for diagnostics
        """
        # Initialise
        if working_dir is None:
            working_dir = os.getcwd()
        self._working_dir = os.path.abspath(working_dir)
        if scripts_dir is None:
            scripts_dir = self._working_dir
        if log_dir is None:
            log_dir = self._working_dir
        if log_file:
            self._log_file = os.path.abspath(log_file)
        if lock_manager:
            self._lock_manager = lock_manager
        if not runner:
            runner = sched.default_runner
        self._runner_nslots = runner.nslots
        # Do setup
        self.invoke(self.setup)
        # Generate commands to run
        cmds = []
        if not batch_size:
            # No batching: one command per job
            for command in self._commands:
                if verbose:
                    self.report("%s" % command.cmd())
                script_file = command.make_wrapper_script(
                    scripts_dir=scripts_dir,
                    envmodules=envmodules,
                    working_dir=self._working_dir)
                cmd = Command('/bin/bash')
                if envmodules:
                    cmd.add_args('-l')
                cmd.add_args(script_file)
                if verbose:
                    self.report("wrapper script %s" % script_file)
                cmds.append(cmd)
                self._scripts.append(script_file)
        else:
            # Batch multiple commands per job
            if verbose:
                self.report("Batching commands (%s commands per "
                            "batch)" % batch_size)
            remaining_cmds = self._commands
            while remaining_cmds:
                # Grab a batch of commands
                batch = remaining_cmds[:batch_size]
                # Combine batch into a single command
                batch_cmd = PipelineCommandWrapper(
                    "Batch commands for %s" % self.name(),
                    *batch[0].cmd().command_line)
                for cmd in batch[1:]:
                    batch_cmd.add_args("&&",
                                       "\\\n",
                                       *cmd.cmd().command_line)
                if verbose:
                    self.report("%s" % batch_cmd.cmd())
                script_file = batch_cmd.make_wrapper_script(
                    scripts_dir=scripts_dir,
                    envmodules=envmodules,
                    working_dir=self._working_dir)
                cmd = Command('/bin/bash')
                if envmodules:
                    cmd.add_args('-l')
                cmd.add_args(script_file)
                if verbose:
                    self.report("wrapper script %s" % script_file)
                cmds.append(cmd)
                self._scripts.append(script_file)
                # Update remaining commands to batch
                remaining_cmds = remaining_cmds[batch_size:]
        # Run the commands
        if cmds:
            # Report runner and nslots
            if verbose:
                if runner:
                    self.report("Runner: %s%s" %
                                (runner,
                                 '' if runner.nslots == 1
                                 else ' [nslots=%s]' % runner.nslots))
            use_group = (len(cmds)!=1)
            if use_group:
                # Run as a group
                group = sched.group(self.id())
                for j,cmd in enumerate(cmds):
                    name = "%s#%s" % (self.id(),j)
                    group.add(cmd,
                              wd=self._working_dir,
                              name=name,
                              runner=runner,
                              log_dir=log_dir,
                              wait_for=wait_for)
                group.close()
                callback_name = group.name
                callback_function = self.task_completed
                self._groups.append(group)
            else:
                # Run a single job
                cmd = cmds[0]
                name = self.id()
                job = sched.submit(cmd,
                                   wd=self._working_dir,
                                   name=name,
                                   runner=runner,
                                   log_dir=log_dir,
                                   wait_for=wait_for)
                callback_name = job.name
                callback_function = self.task_completed
                self._jobs.append(job)
            # Set up a callback which the scheduler will invoke
            # in background when the jobs complete
            sched.callback("%s" % self._name,
                           callback_function,
                           wait_for=(callback_name,))
            if not asynchronous:
                # Wait for job or group to complete before returning
                while not self.completed:
                    time.sleep(poll_interval)
        else:
            # No commands to execute
            self.finish_task()
        return self

    def njobs(self):
        """
        Get number of jobs associated with the task

        Returns the total number of jobs and the
        number of completed jobs associated with a
        running task, as a tuple: (njobs,ncompleted)
        """
        njobs = 0
        ncompleted = 0
        for group in self._groups:
            for job in group.jobs:
                njobs += 1
                try:
                    if job.completed:
                        ncompleted += 1
                except Exception as ex:
                    self.report("Failed to get status for job '%s': %s" %
                                (job.name,ex))
        for job in self._jobs:
            njobs += 1
            try:
                if job.completed:
                    ncompleted += 1
            except Exception as ex:
                self.report("Failed to get status for job '%s': %s" %
                            (job.name,ex))
        return (njobs,ncompleted)

    def terminate(self):
        """
        Internal: terminate the task
        """
        if self.completed:
            return
        for group in self._groups:
            for job in group.jobs:
                job.terminate()
        for job in self._jobs:
            job.terminate()

    def init(self,*args,**kws):
        pass
    def setup(self):
        raise NotImplementedError("Subclass must implement 'setup' method")
    def finish(self):
        raise NotImplementedError("Subclass must implement 'finish' method")

class PipelineFunctionTask(PipelineTask):
    """
    Class enabling a Python function to be run as a task

    A 'function task' enables one or more instances of a
    Python function or class method to be run concurrently
    as external programs, and for the return values of the
    to recovered on task completion.

    This class should be subclassed to implement the 'init',
    'setup', 'finish' (optionally) and 'output' methods.

    The 'add_call' method can be used within 'setup' to add
    one or more function calls to be invoked.
    """

    def __init__(self,_name,*args,**kws):
        """
        Create a new PipelineFunctionTask instance

        Arguments:
          _name (str): an arbitrary user-friendly name for the
            task instance
          args (List): list of arguments to be supplied to
            the subclass (must match those defined in the
            'init' method)
          kws (Dictionary): dictionary of keyword-value pairs
            to be supplied to the subclass (must match those
            defined in the 'init' method)
        """
        PipelineTask.__init__(self,_name,*args,**kws)
        self._dispatchers = []
        self._result = None

    def add_call(self,name,f,*args,**kwds):
        """
        Add a function call to the task

        Arguments:
           name (str): a user friendly name for the call
           f (function): the function or instance method
             to be invoked in the task
           args (list): arguments to be passed to the
             function when invoked
           kwds (mapping): keyword=value mapping to be
             passed to the function when invoked
        """
        d = Dispatcher()
        cmd = d.dispatch_function_cmd(f,*args,**kwds)
        self.add_cmd(PipelineCommandWrapper(name,
                                            *cmd.command_line))
        self._dispatchers.append(d)

    def finish_task(self):
        """
        Internal: perform actions to finish the task
        """
        result = list()
        try:
            for d in self._dispatchers:
                result.append(d.get_result())
        except Exception as ex:
            logger.critical("'%s': exception collecting results from "
                            "dispatchers (%s)" % (self.name(),self.id()))
            self.report("'%s': exception collecting results from "
                        "dispatchers (%s)" % (self.name(),self.id()))
            self.report("Exception: %s" % ex)
            self._completed = True
            self._exit_code = 1
            result = None
        self._result = result
        PipelineTask.finish_task(self)

    def result(self):
        """
        Return the results of the function invocations
        """
        return self._result

class PipelineCommand(object):
    """
    Base class for constructing program command lines

    This class should be subclassed to implement the 'init'
    and 'cmd' methods.

    The 'init' method should do any preprocessing and
    caching of arguments to be used in the 'cmd' method;
    the 'cmd' method should use these to construct and
    return a 'Command' instance.
    """
    def __init__(self,*args,**kws):
        """
        Create a new PipelineCommand instance

        Arguments:
          args (List): list of arguments to be supplied to
            the subclass (must match those defined in the
            'init' method)
          kws (Dictionary): dictionary of keyword-value pairs
            to be supplied to the subclass (must match those
            defined in the 'init' method)
        """
        # Set internal name
        self._name = self.__class__.__name__
        # Invoke the 'init' method
        self.init(*args,**kws)

    def name(self):
        """
        Return a "sanitized" version of the class name
        """
        return sanitize_name(self._name)

    def make_wrapper_script(self,scripts_dir=None,shell="/bin/bash",
                            envmodules=None,working_dir=None):
        """
        Generate a uniquely-named wrapper script to run the command

        Arguments:
          scripts_dir (str): path of directory to write
            the wrapper scripts to
          shell (str): shell to use (defaults to '/bin/bash')
          envmodules (str): list of environment modules to load
          working_dir (str): explicitly specify the directory
            the script should be executed in

        Returns:
          String: name of the wrapper script.
        """
        # Wrap in a script
        if scripts_dir is None:
            scripts_dir = os.getcwd()
        script_file = os.path.join(scripts_dir,"%s.%s.sh" % (self.name(),
                                                             uuid.uuid4()))
        prologue = ["echo \"#### COMMAND %s\"" % self._name,
                    "echo \"#### HOSTNAME $HOSTNAME\"",
                    "echo \"#### USER $USER\"",
                    "echo \"#### START $(date)\""]
        if envmodules:
            shell += " --login"
            for module in envmodules:
                if module is not None:
                    prologue.append("module load %s" % module)
        if working_dir:
            prologue.append("cd %s" % working_dir)
        prologue.append("echo \"#### CWD $(pwd)\"")
        epilogue = ["exit_code=$?",
                    "echo \"#### END $(date)\"",
                    "echo \"#### EXIT_CODE $exit_code\"",
                    "exit $exit_code"]
        self.cmd().make_wrapper_script(filen=script_file,
                                       shell=shell,
                                       prologue='\n'.join(prologue),
                                       epilogue='\n'.join(epilogue),
                                       quote_spaces=True)
        return script_file

    def init(self):
        """
        Initialise and store parameters

        Must be implemented by the subclass
        """
        raise NotImplementedError("Subclass must implement 'init' method")

    def cmd(self):
        """
        Build the command

        Must be implemented by the subclass and return a
        Command instance
        """
        raise NotImplementedError("Subclass must implement 'cmd' method")

class PipelineCommandWrapper(PipelineCommand):
    """
    Class for constructing program command lines

    This class is based on the PipelineCommand class but
    can be used directly (rather than needing to be
    subclassed).

    For example, to wrap the 'ls' command directly:

    >>> ls_command = PipelineCommandWrapper("List directory",'ls',dirn)

    It is also possible to extend the command line
    using the 'add_args' method, for example:

    >>> ls_command = PipelineCommandWrapper("List directory",'ls')
    >>> ls.command.add_args(dirn)
    """
    def __init__(self,name,*args):
        """
        Create a new PipelineCommandWrapper instance

        Arguments:
          name (str): arbitrary name for the command
          args  (List): initial list of arguments making
            up the command
        """
        PipelineCommand.__init__(self,*args)
        self._name = str(name)
        self._cmd = None
        if args:
            self._cmd = Command(*args)

    def add_args(self,*args):
        """
        Add additional arguments to extend the command being built

        Arguments:
          args  (List): one or more arguments to append to
            the command
        """
        if self._cmd is None:
            self._cmd = Command(*args)
        else:
            self._cmd.add_args(*args)

    def init(self,*args):
        """
        Internal: dummy init which does nothing
        """
        pass

    def cmd(self):
        """
        Internal: implement the 'cmd' method
        """
        return self._cmd

######################################################################
# Parameter-like utility classes for passing values between tasks
######################################################################

class FunctionParam(BaseParam):
    """
    Class for deferred function evaluation as pipeline parameter

    This class wraps a function with a set of
    parameters; the function evaluation is
    deferred until the 'value' property is invoked.

    Any parameters which are PipelineParam-like
    instances will be replaced with their values
    before being passed to the function.
    """
    def __init__(self,f,*args,**kws):
        """
        Create a new FunctionParam instance

        Arguments:
          f (object): function-like object to
            be evaluated
          args (list): positional arguments to
            pass to the function on evaluation
          kws (mapping): keyworded arguments to
            pass to the function on evaluation
        """
        BaseParam.__init__(self)
        self._f = f
        self._args = args
        self._kws = kws
    @property
    def value(self):
        """
        Return value from evaluated function
        """
        args = []
        for arg in self._args:
            try:
                args.append(arg.value)
            except AttributeError:
                args.append(arg)
        kws = {}
        for kw in self._kws:
            try:
                kws[kw] = self._kws[kw].value
            except AttributeError:
                kws[kw] = self._kws[kw]
        return self._f(*args,**kws)

class PathJoinParam(FunctionParam):
    """
    Class for joining file paths as pipeline parameter

    This class implements the pipeline parameter
    equivalent of the `os.path.join` function, taking
    a set of path elements on instantiation (which
    can be strings or PipelineParam-like objects)
    and returning the joined path elements on
    evaluation via the `value` property.

    Example usage:

    >>> pth = PathJoinParam("/path","to","file.txt")
    >>> pth.value
    "/path/to/file.txt"

    >>> base_dir = PipelineParam(value="/path/to/base")
    >>> pth = PathJoinParam(base_dir,"file.txt")
    >>> pth.value
    "/path/to/base/file.txt"
    >>> base_dir.set("/path/to/new/base")
    >>> pth.value
    "/path/to/new/base/file.txt"

    Note that this class doesn't implement a `set`
    method (unlike the standard PipelineParam class)
    so the path elements cannot be changed after
    initialisation.
    """
    def __init__(self,*p):
        """
        Create a new PathJoinParam instance

        Arguments:
          p (iterable): list of path elements
            to join; can be strings or
            PipelineParam-like objects
        """
        FunctionParam.__init__(self,
                               os.path.join,
                               *p)

class PathExistsParam(FunctionParam):
    """
    Class for checking file/directory existance as pipeline parameter

    This class implements the pipeline parameter
    equivalent of the `os.path.exists` function, taking
    a path on instantiation (which can be a string
    or PipelineParam-like object) and returning a
    boolean value indicating whether the path is an
    existing file system object via the `value`
    property.

    Example usage:

    >>> exists = PathExistsParam("/path/to/file.txt")
    >>> exists.value
    True

    >>> f = PipelineParam(value="/path/to/file.txt")
    >>> exists = PathExistsParam(f)
    >>> exists.value
    True
    >>> f.set("/path/to/missing.txt")
    >>> exists.value
    False

    Note that this class doesn't implement a `set`
    method (unlike the standard PipelineParam class)
    so the path elements cannot be changed after
    initialisation.
    """
    def __init__(self,p):
        """
        Create a new PathExistsParam instance

        Arguments:
          p (str): path to check existance of;
            can be a string or a PipelineParam-like
            object
        """
        FunctionParam.__init__(self,
                               os.path.exists,
                               p)

######################################################################
# Dispatching Python functions as separate processes
######################################################################

class Dispatcher(object):
    """
    Class to invoke Python functions in external processes

    Enables a Python function to be run in a separate process
    and collect the results.

    Example usage:

    >>> d = Dispatcher()
    >>> cmd = d.dispatch_function_cmd(bcftbx.utils.list_dirs,os.get_cwd())
    >>> cmd.run_subprocess()
    >>> result = d.get_result()

    The Dispatcher works by pickling the function, arguments
    and keywords to files, and then creating a command which
    can be run as an external process; this unpickles the
    function etc, executes it, and makes the result available
    to the dispatcher on successful completion.

    Currently uses cloudpickle as the pickling module:
    https://pypi.org/project/cloudpickle/
    """
    def __init__(self,working_dir=None,cleanup=True):
        """
        Create a new Dispatcher instance

        Arguments:
          working_dir (str): optional, explicitly specify the
            directory where pickled files and other
            intermediates will be stored
          cleanup (bool): if True then remove the working
            directory on exit
        """
        self._id = str(uuid.uuid4())
        if not working_dir:
            working_dir = "dispatcher.%s" % self._id
        self._working_dir = os.path.abspath(working_dir)
        self._pickled_result_file = None
        self._pickler = cloudpickle
        # Use the current Python interpreter when
        # running the function
        self._python = sys.executable
        if cleanup:
            atexit.register(self._cleanup)

    @property
    def working_dir(self):
        """
        Return the working directory path
        """
        if not os.path.exists(self._working_dir):
            logger.debug("Dispatcher: making working directory: "
                         "%s" % self._working_dir)
            os.mkdir(self._working_dir)
        return self._working_dir

    def _cleanup(self):
        """
        Internal: remove the working directory and its contents
        """
        if os.path.exists(self._working_dir):
            logger.debug("Dispatcher: cleaning up dir %s" %
                         self._working_dir)
            shutil.rmtree(self._working_dir)

    def execute(self,pkl_func_file,pkl_args_file,pkl_kwds_file,
                pkl_result_file=None):
        """
        Internal: execute the function

        Arguments:
          pkl_func_file (str): path to the file with the pickled
            version of the function to be invoked
          pkl_args_file (str): path to the file with the pickled
            version of the function arguments
          pkl_args_file (str): path to the file with the pickled
            version of the keyword arguments
          pkl_result_file (str): optional path to the file where
            the pickled version of the function result will be
            written
        """
        # Unpickle the components
        logger.info("Dispatcher: running 'execute'...")
        try:
            func = self._unpickle_object(pkl_func_file)
        except AttributeError as ex:
            # Note that functions are pickled by 'fully qualified'
            # name reference, not by value i.e. only the function
            # name and the name of the parent module are pickled
            print("Failed to load unpickled function: %s" % ex)
            return 1
        args = self._unpickle_object(pkl_args_file)
        kwds = self._unpickle_object(pkl_kwds_file)
        logger.info("Dispatcher: executing:")
        logger.info("-- function: %s" % func)
        logger.info("-- args    : %s" % (args,))
        logger.info("-- kwds    : %s" % (kwds,))
        # Execute the function
        try:
            result = func(*args,**kwds)
            logger.info("Dispatcher: result: %s" % (result,))
        except Exception as ex:
            print("Dispatched function raised exception: %s" % ex)
            result = ex
        # Pickle the result
        if pkl_result_file:
            self._pickle_object(result,pkl_result_file)
        return 0

    def dispatch_function_cmd(self,func,*args,**kwds):
        """
        Generate a command to execute the function externally

        Arguments:
          func (function): function to be executed
          args (list): arguments to invoke the function with
          kwds (mapping): keyword arguments to invoke the
            function with

        Returns:
          Command: a Command instance that can be used to
            execute the function.
        """
        logger.info("Dispatching function:")
        logger.info("-- function: %s" % func)
        logger.info("-- args    : %s" % (args,))
        logger.info("-- kwds    : %s" % (kwds,))
        # Pickle the components
        pkl_func_file = self._pickle_object(func)
        pkl_args_file = self._pickle_object(args)
        pkl_kwds_file = self._pickle_object(kwds)
        # Filename for results
        self._pickled_result_file = os.path.join(self.working_dir,
                                                 str(uuid.uuid4()))
        # Generate command to run the function
        return Command(self._python,
                       "-c",
                       "from auto_process_ngs.pipeliner import Dispatcher ; "
                       "Dispatcher().execute(\"%s\",\"%s\",\"%s\",\"%s\")" %
                       (pkl_func_file,
                        pkl_args_file,
                        pkl_kwds_file,
                        self._pickled_result_file))

    def get_result(self):
        """Return the result from the function invocation

        Returns:
          Object: the object returned by the function, or None
            if no result was found.
        """
        if os.path.exists(self._pickled_result_file):
            result = self._unpickle_object(self._pickled_result_file)
            if isinstance(result,Exception):
                raise result
            else:
                return result
        else:
            raise Exception("Pickled output not found")

    def _pickle_object(self,obj,pickle_file=None):
        """
        Internal: pickle an object and write to file

        Arguments:
          obj (object): object to be pickled
          pickle_file (str): specified path to output
            file (optional)

        Returns:
          String: path to the file containing the
            pickled object.
        """
        if not pickle_file:
            pickle_file = os.path.join(self.working_dir,
                                       str(uuid.uuid4()))
        with open(pickle_file,'wb') as fp:
            fp.write(self._pickler.dumps(obj))
        return pickle_file

    def _unpickle_object(self,pickle_file):
        """
        Internal: unpickle an object from a file

        Arguments:
          pickle_file (str): path to output file
            with pickled data

        Returns:
          Object: the unpickled object.
        """
        with open(pickle_file,'rb') as fp:
            return self._pickler.loads(fp.read())

######################################################################
# Utility functions
######################################################################

def sanitize_name(s):
    """
    Convert string to lowercase and replace special characters

    Arguments:
      s (str): string to sanitize

    Returns:
      String: sanitized string.
    """
    name = []
    # Convert to lower case and replace special characters
    for c in str(s).lower():
        if c not in ALLOWED_CHARS:
            if len(name) < 2 or name[-2:] != '__':
                name.append('_')
        else:
            name.append(c)
    return ''.join(name)

def collect_files(dirn,pattern):
    """
    Return names of files in a directory which match a glob pattern

    Arguments:
      dirn (str): path to a directory containing the files
      pattern (str): a glob pattern to match

    Returns:
      List: list of matching files
    """
    return sorted(glob.glob(os.path.join(os.path.abspath(dirn),pattern)))

def resolve_parameter(p):
    """
    Resolve the value of an arbitrary parameter

    The supplied "parameter" can be any object; if it
    has a 'value' property then the resolved value will
    whatever this returns, otherwise the supplied
    object is returned as-is.

    Arguments:
      p (object): parameter to resolve

    Returns:
      Object: resolved parameter value.
    """
    try:
        # If arg is a PipelineParam then convert it
        # and store the final value instead of the
        # PipelineParam instance
        return p.value
    except AttributeError:
        return p
