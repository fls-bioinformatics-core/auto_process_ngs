Set up and configuration
========================

Requirements
************

The autoprocessor is written in Python and should work with versions 2.6
and 2.7. It also depends on the `bcftbx` Python module and NGS QC scripts
which are part of the `genomics` github repository.

In addition the following software must be installed:

* bcl2fastq http://support.illumina.com/downloads/bcl2fastq_conversion_software_184.html
* bowtie http://bowtie-bio.sourceforge.net/index.shtml
* fastq_screen http://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/
* fastqc http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

The programs provided by these packages must be found on the `PATH` when
the appropriate autoprocessor commands are issued. :ref:`environment-modules`
can be used to help manage this.

Configuration: settings.ini
***************************

The autoprocessor reads its global settings for the local system from a
`settings.ini` file, which it looks for in order in the following
locations:

1. The current directory;
2. The ``config`` subdirectory of the installation directory;
3. The installation directory (for legacy installations only)

If no `settings.ini` file is found then one can be created using the
command::

    auto_process.py config --init

otherwise the autoprocessor will run using the built-in default values.

To see the current settings, do::

    auto_process.py config

To update the settings use the ``--set`` options, for example::

    auto_process.py config --set bcl2fastq.nprocessors=4

The most important settings are the :ref:`job-runners` and for any
:ref:`environment-modules` that you wish to specify for a particular
processing stage.

.. _job-runners:

Job runner specification
------------------------

Job runners tell the autoprocessor how to run programs. There are
currently only two available:

* ``SimpleJobRunner``: runs programs as a subprocess of the current process
* ``GEJobRunner``: runs programs using Grid Engine (GE)

The ``GEJobRunner`` is recommended when using the autoprocessor on cluster
systems. To specify additional Grid Engine-specific options to use with
the runner, enclose them in parentheses e.g.::

    [runners]
    bcl2fastq = GEJobRunner(-pe smp.pe 8)

.. note::

   If you specify multiple processors for the ``bcl2fastq`` and are using
   the ``GEJobRunner`` then you should ensure that the job runner requests
   a suitable number of cores when submitting jobs.

.. _environment-modules:

Environment modules
-------------------

`Environment modules <http://modules.sourceforge.net/>`_ provide a way to
dynamic modify the user's environment. They can be especially useful to
provide access to multiple versions of the same software package, and to
manage conflicts between packages.

The ``[modulefiles]`` directive allows specific module files to be loaded
before a specific step, for example::

    [modulefiles]
    make_fastqs = apps/bcl2fastq/1.8.4
