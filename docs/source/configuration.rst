
.. _auto_process_configuration:

*************
Configuration
*************

The autoprocessor reads its global settings for the local system from a
``settings.ini`` file, which it looks for in order in the following
locations:

1. The current directory
2. The ``config`` subdirectory of the installation directory
3. The installation directory (for legacy installations only)

If no ``settings.ini`` file is found then one can be created using the
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

   If you specify multiple processors for the ``bcl2fastq`` runner and are
   using ``GEJobRunner`` then you should ensure that the job runner requests
   a suitable number of cores when submitting jobs.

.. _environment-modules:

Environment modules
-------------------

`Environment modules <http://modules.sourceforge.net/>`_ provide a way to
dynamic modify the user's environment. They can be especially useful to
provide access to multiple versions of the same software package, and to
manage conflicts between packages.

The ``[modulefiles]`` directive in ``settings.ini`` allows specific module
files to be loaded before a specific step, for example::

    [modulefiles]
    make_fastqs = apps/bcl2fastq/1.8.4

.. note::

   These can be overridden for the ``make_fastqs`` and ``run_qc`` using
   the ``--modulefiles`` option.

.. _required_bcl2fastq_versions:

Required bcl2fastq versions and other settings
----------------------------------------------

Different versions of Illumina's ``bcl2fastq`` software can be specified
both as a default and dependent on the sequencer platform, by setting the
appropriate parameters in the ``settings.ini`` file.

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

.. warning::

   Previously the ``[bcl2fastq]`` directive allowed the versions to be
   set using platform names specified within that section, for example::

        [bcl2fastq]
        ...
        hiseq = 1.8.4

   This method is now deprecated in favour of the ``[platform:NAME]``
   mechanism.

   If this old method is detected then warnings are issued and the
   software attempts to make an intelligent choice about the versions.
