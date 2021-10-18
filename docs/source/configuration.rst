
.. _auto_process_configuration:

*************
Configuration
*************

--------
Overview
--------

The autoprocessor reads its global settings for the local system from a
``auto_process.ini`` file, which it looks for in order in the following
locations:

1. The file specified by the ``AUTO_PROCESS_CONF`` environment
   variable (if set)
2. The current directory
3. The ``config`` subdirectory of the installation directory
4. The installation directory (for legacy installations only)

.. note::

   In previous versions of the package the configuration file was
   called ``settings.ini``, and this will be used as a fallback if
   no ``auto_process.ini`` file is found in the locations above.

To create a ``auto_process.ini`` file for a new installation, use the
command

::

    auto_process.py config --init

To see the current settings, do

::

    auto_process.py config


To change the settings, either use the ``--set`` options, for example

::

    auto_process.py config --set bcl2fastq.nprocessors=4

or simply edit the ``auto_process.ini`` file by hand using a text editor.


.. note::

   If no ``auto_process.ini`` file exists then ``auto_process`` will run
   using the built-in default values.

.. note::

   Many of the configuration options can be over-ridden at run time
   using command line options for the specific ``auto_process``
   subcommands.

.. _basic_configuration:

-------------------
Basic configuration
-------------------

Using the basic ``auto_process`` Fastq generation requires minimal
configuration when running locally; provided that the required
``bcl2fastq`` software is available on the system (see
:ref:`software_dependencies`) it should run without further setup.

Running the QC pipeline requires additional software plus reference data
which are described in more detail in :ref:`reference_data`.

.. _config_sequencer_platforms:

------------------------
Sequencers and platforms
------------------------

Information about sequencers used at the local site can be stored
in ``sequencer`` sections of the configuration file.

A sequencer can be defined by adding a new section of the form
``[sequencer:INSTRUMENT_NAME]``, where ``INSTRUMENT_NAME`` is the
ID name for the instrument. Within each section the following
data items can then be associated with the sequencer:

============= ==============================================
``platform``  **Compulsory** sets the generic platform name
              (one of ``hiseq``, ``hiseq4000``, ``nextseq``,
	      ``miseq``, ``miniseq`` or ``iseq``)
``model``     Text describing the sequencer model (e.g.
              ``HiSeq 4000``)
============= ==============================================

For example: if the local facility has a HiSeq 4000 instrument
with ID ``SN7001250`` then this would be defined in ``auto_process.ini``
as follows:

::

   [sequencer:SN7001250]
   platform = hiseq4000
   model = "HiSeq 4000"

The instrument name can be derived from the name of the directories
produced by the sequencer (see :ref:`run_and_fastq_naming_conventions`).


.. note::

   These sections replace the old ``sequencers`` section used
   to define the sequencer platforms, e.g.

   ::

      [sequencers]
      SN7001250 = hiseq4000

  This section is still supported but is now deprecated.

----------------
Default metadata
----------------

The ``metadata`` section of the configuration file allows defaults
to be specified for metadata items associated with each run.

Currently it is possible to set a default for the ``source``
metadata item, which specifies where the data was received from,
for example:

::

   [metadata]
   default_data_source = "Local sequencing facility"

If no default is set then the values can be updated using the
``metadata`` command (see :ref:`commands_metadata`).

.. _job_runners:

-----------
Job Runners
-----------

*Job runners* are used within ``auto_process`` to tell the pipelines
how to execute commands. There are currently two types of runner available:

* ``SimpleJobRunner`` runs jobs as a subprocess of the current process,
  so they run locally (i.e. on the same hardware that the ``auto_process``
  command was started on)
* ``GEJobRunner`` submits jobs to Grid Engine (GE), which enables it to
  exploit additional resources available on a compute cluster (see
  :ref:`running_on_compute_cluster`)

Job runners can also be configured to specify the number of CPUs
available to commands that are executed using them (see
:ref:`setting_number_of_cpus`).

By default ``auto_process`` is configured to use ``SimpleJobRunner``
for all jobs; the default runner is defined in the settings:

::

   [general]
   default_runner = SimpleJobRunner

This default can be overridden for specific commands and pipeline
stages by explicitly specifying alternative runners in the ``runners``
section of the settings file:

============================= =========================================
Runner name                   Used for
============================= =========================================
``bcl2fastq``                 Running ``bcl2fastq`` in Fastq generation
``stats``                     Running commands to generate statistics
                              after Fastq generation (e.g.
			      ``fastq_statistics.py``)
``rsync``                     Running commands for transferring data
                              (e.g. copying primary data for Fastq
                              generation, archiving etc)
``qc``                        Running computationally intensive QC
                              commands (e.g. ``FastQC``, ``Fastq_screen``,
                              etc)
``star``                      Running pipeline tasks which use ``STAR``
                              (e.g. strandedness, alignment etc)
``cellranger``                Running ``cellranger`` in Fastq generation
                              and QC pipelines
``icell8``                    Default runner for commands in the ICELL8
                              processing pipeline
``icell8_contaminant_filter`` Running the contaminant filtering in the
                              ICELL8 pipeline
``icell8_statistics``         Generating statistics for ICELL8 data
``icell8_report``             Reporting on the ICELL8 pipeline
============================= =========================================

.. _setting_number_of_cpus:

--------------------------------
Setting number of available CPUs
--------------------------------

Job runners allow the number of available CPUs (aka processors or
threads) to be specified, and this information is then used when
running jobs in the ``auto_process`` pipelines.

For ``SimpleJobRunners`` the number of CPUs is specified via the
``nslots`` argument. For example:

::

   [runners]
   qc = SimpleJobRunner(nslots=8)

(Without ``nslots`` the number of CPUs implicitly defaults to 1.)

For ``GEJobRunners`` the number of available CPUs is inferred from the
``-pe smp.pe`` argument (see :ref:`running_on_compute_cluster`).

For some commands the number of available CPUs will be taken implicitly
from this argument unless explicitly overridden by the following settings:

=============== ================================== =====================
Section         Setting                            Overrides runner
=============== ================================== =====================
``bcl2fastq``   ``nprocessors``                    ``bcl2fastq``
``fastq_stats`` ``nprocessors``                    ``stats``
``qc``          ``nprocessors``                    ``qc``
``icell8``      ``nprocessors_contaminant_filter`` ``icell8_contaminant_filter``
``icell8``      ``nprocessors_statistics``         ``icell8_statistics``
``10xgenomics`` ``cellranger_localcores``          ``cellranger`` (*)
=============== ================================== =====================

(*) Used when ``cellranger`` is run with ``--jobmode=local``

.. _running_on_compute_cluster:

----------------------------
Running on a compute cluster
----------------------------

The ``GEJobRunner`` can be used to make ``auto_process`` submit its
computationally intensive jobs to a compute cluster rather than on
the local host; to switch to using ``GEJobRunner``, set the default
runner in the settings:

::

   [general]
   default_runner = GEJobRunner

Additional options for Grid Engine submission can be specified by
enclosing when defining the runner, for example sending all jobs to a
particular queue might use:

::

   default_runner = GEJobRunner(-q ngs.queue)

This default runer can further be overridden for specific commands
and pipeline stages by the settings in the ``runners`` section of the
configuration file (see the previous section :ref:`job_runners`).

For example: to run ``bcl2fastq`` jobs in parallel environment
with 8 cores might look like:

::

   [runners]
   bcl2fastq = GEJobRunner(-pe smp.pe 8)


.. note::

   If you specify multiple processors for the ``bcl2fastq`` runner and are
   using ``GEJobRunner`` then you should ensure that the job runner requests
   a suitable number of cores when submitting jobs.

.. note::

   When running on a cluster the ``auto_process`` driver process should
   run on the cluster login node; it has a small CPU and memory footprint
   which should impact minimally on other users of the system.

.. _limiting_number_of_jobs:

------------------------------------------
Managing concurrent jobs and process loads
------------------------------------------

There are a number of settings available in the ``[general]``
section which allow limits to be set on the resources that
``auto_process`` will try to consume when running jobs and
pipelines:

======================= =============================================
Setting
======================= =============================================
``max_concurrent_jobs`` Maximum number of jobs that ``auto_process``
                        is allowed to run at one time
``max_cores``           Maximum number of cores that ``auto_process``
                        is allowed to use at one time across all
                        jobs
``max_batches``         Dynamically sets batch sizes within pipelines
                        so that number of job batches from each task
                        doesn't exceed this number
======================= =============================================

For example:

::

    [general]
    max_cores = 24

If any of these is set to zero or ``None`` then this means
that resource is not limited by ``auto_process``.

``max_concurrent_jobs`` and ``max_batches`` are useful on
shared cluster systems, to avoid submitting large numbers of
jobs at one time.

``max_cores`` is useful when running on a local workstation,
to avoid exceeding resource limits while ensuring the most
efficient use of the available CPUs.

.. _environment_modules:

-------------------------
Using environment modules
-------------------------

`Environment modules <http://modules.sourceforge.net/>`_ provide a way to
dynamically modify the user's environment. They can be especially useful to
provide access to multiple versions of the same software package, and to
manage conflicts between packages.

The ``[modulefiles]`` section in ``auto_process.ini`` allows specific module
files to be loaded before a specific step, for example::

    [modulefiles]
    make_fastqs = apps/bcl2fastq/2.20

These can be defined for the following stages:

 * ``make_fastqs``
 * ``run_qc``
 * ``publish_qc``
 * ``process_icell8``

(see :ref:`software_dependencies` for details of what software is required
for each of these stages.)

.. note::

   These can be overridden for the ``make_fastqs`` and ``run_qc`` stages
   using the ``--modulefiles`` option.

Environment modules for Fastq generation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the ``make_fastqs`` stage, additional module files can be specified
for individual tasks with the Fastq generation pipeline:

* ``bcl2fastq``
* ``cellranger_mkfastq``
* ``cellranger_atac_mkfastq``
* ``cellranger_arc_mkfastq``

If any of these are defined then they will be loaded for the relevant
tasks in the Fastq generation pipeline.

Environment modules for the QC pipeline
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the ``run_qc`` stage, additional module files can be specified for
individual tasks within the QC pipeline:

 * ``illumina_qc``
 * ``fastq_strand``
 * ``cellranger``
 * ``report_qc``

If any of these are defined then they will be loaded for the relevant
tasks in the QC pipeline.

.. _conda_dependency_resolution:

--------------------------------------------
Using conda to resolve pipeline dependencies
--------------------------------------------

For certain pipelines and tasks it is possible to enable the ``conda``
package management utility to handle setting up appropriate run-time
environments, rather than having to manually install the required
dependencies and specify their locations (e.g. using environment
modules).

To do this by default, set the ``enable_conda`` parameter in the
``[conda]`` section, i.e.::

    [conda]
    enable_conda = true

Note that this requires ``conda`` to be installed and available on the
user's ``PATH`` at run-time.

By default a temporary directory will be used when creating and reusing
``conda`` environments, but this can be overriden by setting the
``env_dir`` parameter, e.g.::

    [conda]
    enable_conda = true
    env_dir = $HOME/conda_envs

.. _required_bcl2fastq_versions:

---------------------------
Required bcl2fastq versions
---------------------------

Different versions of Illumina's ``bcl2fastq`` software can be specified
both as a default and dependent on the sequencer platform, by setting the
appropriate parameters in the ``auto_process.ini`` file.

The ``[bcl2fastq]`` directive specifies the defaults to use for all
platforms in the absence of more specific settings, for example::

    [bcl2fastq]
    default_version = 1.8.4
    nprocessors = 8

These settings can be overriden for specific platforms, by creating optional
directives of the form ``[platform:NAME]`` (where ``NAME`` is the name of the
platform). For example to set the version to use when processing data from a
NextSeq instrument to be specifically ``2.17.1.14``::

    [platform:nextseq]
    bcl2fastq = 2.17.1.14

A range of versions can be specified by prefacing the version number by
one of the operators ``>``, ``>=``, ``<=`` and ``<`` (``==`` can also be
specified explicitly), for example::

    bcl2fastq = >=2.0

Alternatively a comma-separated list can be provided::

    bcl2fastq = >=1.8.3,<2.0

If no bcl2fastq version is explicitly specified then the highest available
version will be used.

.. note::

   This mechanism allows multiple ``bcl2fastq`` versions to be present
   in the environment simultaneously.

.. _data_transfer_destinations:

--------------------------
Data transfer destinations
--------------------------

The ``transfer_data.py`` utility can be used to copy Fastqs and other
data produced by the ``auto_process.py`` pipeline to arbitrary
destinations, typically for sharing with end users of the pipeline.

The utility provides a number of command line options to specify a
destination and the data that are transferred at runtime. However it
is also to define one or more destinations in the configuration file,
with appropriate presets for each destination.

A destination can be defined by adding a new section to the config
file of the form ``[destination:NAME]``, where ``NAME`` is the name
that will be used to refer to the destination when it is specified in
a run of ``transfer_data.py``.

Within each section the following parameters can be set for the
destination:

====================== ==============================================
Parameter              Function
====================== ==============================================
``directory``          **Compulsory** sets the destination directory
                       to copy files to; can be an arbitrary location
                       of the form ``[[USER@]HOST:]DIR``
``subdir``             Subdirectory naming scheme
``readme_template``    Template file to generate ``README`` from
``url``                Base URL to access copied data at
``include_downloader`` Whether to include ``download_fastqs.py``
``include_qc_report``  Whether to include zipped QC reports
``hard_links``         Whether to hard link to Fastqs rather making
                       copies (for local directories on the same file
                       system as the original Fastqs)
====================== ==============================================

For example:

::

    [destination:webserver]
    directory = /mnt/hosted/web
    subdir = random_bin
    readme_template = README.webserver
    url = http://ourdata.com/shared
    hard_links = true

See :ref:`transfer_data` for more information on what these settings do.

-------------------
Bash tab completion
-------------------

The ``auto_process-completion.bash`` file (installed into the
``etc/bash_completion.d`` subdirectory of the installation location) can
used to enable tab completion of auto_process.py commands within ``bash``
shells.

* For a global installation, copy the file to the system's
  ``/etc/bash_completion.d/`` directory, to make it available
  to all users
* For a local installation, source the file when setting up the
  environment for the installation (or source it in your
  ``~/.bashrc`` or similar).
