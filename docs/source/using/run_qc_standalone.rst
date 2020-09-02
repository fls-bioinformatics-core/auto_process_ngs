Running the QC standalone with ``run_qc.py``
============================================

The utility ``run_qc.py`` allows the QC pipeline to be run on an
arbitrary set of Fastqs outside of the ``auto_process`` pipeline.

The general invocation is:

::

   run_qc.py DIR [ DIR ... ]

where ``DIR`` is a directory with the Fastq files to run the QC
on (or which has a ``fastqs`` subdirectory with the Fastqs; use
the ``--fastq_dir`` option to specify a different subdirectory).

Some of the most commonly used options are:

* ``--protocol``: specify the QC protocol
* ``--organism``: specify the organism(s)
* ``--multiqc``: turns on generation of MultiQC reports

(See the documentation in the :ref:`utilities_run_qc` section
for the full range of options.)

Running on different platforms: ``--local``
-------------------------------------------

By default the QC pipeline will run using the settings from the
``auto_process`` configuration file; however it is recommended
to use the ``--local`` option if running the QC on a local
workstation, or within a job submitted to a compute cluster
(for example, if running inside another script).

For example: submitting a QC run as a single job on a Grid
Engine compute cluster might look like:

::

   qsub -b y -V -pe smp.pe 16 'run_qc.py --local /data/Fastqs'

In this mode the pipeline overrides the central configuration
and attempts to adjust parameters for running the QC to suit
the local setup.

It should make reasonable guesses for the number of available
CPUs and memory. However the following options can be used with
``--local`` to override the guesses:

* ``--maxcores``: sets the maximum number of CPUs available;
  the QC will not exceed this number when running jobs. If
  this isn't set explicitly then the pipeline will attempt to
  determine the number of CPUs automatically;
* ``--maxmem``: sets the maximum amount of memory available
  (in Gbs); currently this is only used if ``cellranger`` is
  being run. If this isn't set explicitly then the pipeline will
  attempt to determine the available memory automatically.

For example: submitting a QC run as a single job on a Grid
Engine compute cluster might look like:

::

   qsub -b y -V -pe smp.pe 16 'run_qc.py --local --maxcores=16 --maxmem=64 /data/Fastqs'
