Set up and configuration
========================

Requirements
************

The autoprocessor is written in Python and should work with versions 2.6
and 2.7. It also depends on the ``genomics-bcftbx`` Python module and NGS
QC scripts which are part of the ``genomics`` github repository.

In addition the following software must be installed:

* ``bcl2fastq``: one or both of:

  * ``bcl2fastq`` 1.8.4: http://support.illumina.com/downloads/bcl2fastq_conversion_software_184.html
  * ``bcl2fastq`` 2.17: https://support.illumina.com/downloads/bcl2fastq-conversion-software-v217.html

* ``bowtie`` http://bowtie-bio.sourceforge.net/index.shtml
* ``fastq_screen`` http://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/
* ``fastqc`` http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

The programs provided by these packages must be found on the ``PATH`` when
the appropriate autoprocessor commands are issued. :ref:`environment-modules`
can be used to help manage this.

..  note::

    If there are multiple ``bcl2fastq`` packages available then see
    :ref:`required_bcl2fastq_versions` for how to specify which version
    is used.

Configuration: settings.ini
***************************

The autoprocessor reads its global settings for the local system from a
``settings.ini`` file, which it looks for in order in the following
locations:

1. The current directory;
2. The ``config`` subdirectory of the installation directory;
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

Required bcl2fastq versions
---------------------------

Different versions of Illumina's ``bcl2fastq`` software can be specified
depending on the sequencer platform, by setting the appropriate parameters
in the ``[bcl2fastq]`` directive.

For example, to specify the version to use when processing data from a
NextSeq instrument to be specifically ``2.17.1.14``::

    [bcl2fastq]
    ...
    nextseq = 2.17.1.14

A range of versions can be specified by prefacing the version number by
one of the operators ``>``, ``>=``, ``<=`` and ``<`` (``==`` can also be
specified explicitly), for example::

    nextseq = >=2.0

Alternatively a comma-separated list can be provided::

    hiseq = >=1.8.3,<2.0

The ``default_version`` sets the required version implicitly in the absence
of an explicit specification, for example::

    [bcl2fastq]
    ...
    default_version = 1.8.4
    hiseq = None

If the ``default_version`` is set to ``None`` then the highest available
version will be used.

.. note::

   This mechanism allows multiple ``bcl2fastq`` versions to be present
   in the environment simultaneously.

Bash tab completion
*******************

The ``auto_process-completion.bash`` file (installed into the
``etc/bash_completion.d`` subdirectory of the installation location)
can used to enable tab completion of ``auto_process.py`` commands
within ``bash`` shells.

* For a global installation, copy the file to the system's
  ``/etc/bash_completion.d/`` directory, to make it available
  to all users
* For a local installation, source the file when setting up the
  environment for the installation (or source it in your ``~/.bashrc``
  or similar).


