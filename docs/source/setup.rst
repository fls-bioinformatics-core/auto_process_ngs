Set up and configuration
========================

Requirements
************

The autoprocessor is written in Python and should work with Python 2.7.

A number of additional Python packages are explicitly required but these
should be installed automatically from the Python Package Index (PyPI)
if the autoprocessor is installed via ``pip``:

 * ``pillow``
 * ``matplotlib``
 * ``pandas``
 * ``nebulizer``

In addition it also depends on the ``genomics-bcftbx`` Python module,
which is not currently available from PyPI; it can be obtained from
the Github repository at
https://github.com/fls-bioinformatics-core/genomics or can be installed
via the ``requirements.txt`` file (see :ref:`auto_process_installation`).

Some of the commands depend on external software packages which must
also be installed. These include

 * ``bcl2fastq``: used for generating Fastq files from raw Illumina
   sequencing data in the ``make_fastqs`` command. Version 2.17+ is
   required for NextSeq data.

   - ``bcl2fastq`` 2.17: https://support.illumina.com/downloads/bcl2fastq-conversion-software-v217.html
   - ``bcl2fastq`` 1.8.4: http://support.illumina.com/downloads/bcl2fastq_conversion_software_184.html

 * QC pipeline:

   - ``bowtie`` http://bowtie-bio.sourceforge.net/index.shtml
   - ``fastq_screen`` http://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/
   - ``fastqc`` http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
   - ``multiqc`` http://multiqc.info/

 * ICell8 processing pipeline:

   - ``cutadapt`` http://cutadapt.readthedocs.io
   - ``bowtie2`` (optional, required if using Bowtie2 genome indexes
     for contaminant filtering - see :ref:`icell8`)
     http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
   - Plus software required for the QC pipeline

Some of this software can be obtained from ``bioconda`` - see
https://bioconda.github.io/ from more details.

The programs provided by these packages must be found on the ``PATH`` when
the appropriate autoprocessor commands are issued. :ref:`environment-modules`
can be used to help manage this.

..  note::

    If there are multiple ``bcl2fastq`` packages available on the path
    at run time then see :ref:`required_bcl2fastq_versions` for how to
    specify which version is used.

.. _auto_process_installation:

Installation
************

The autoprocessing package can obtained from its Github respository at

 * https://github.com/fls-bioinformatics-core/auto_process_ngs

The available versions can be found via
https://github.com/fls-bioinformatics-core/auto_process_ngs/releases

It is recommended to download one of these versions and unpack locally,
e.g.::

    wget https://github.com/fls-bioinformatics-core/auto_process_ngs/archive/0.6.4.tar.gz
    tar zxf 0.6.4.tar.gzhttps://github.com/fls-bioinformatics-core/auto_process_ngs/releases/tag/0.6.4

which will download and unpack version 0.6.4 into a new directory
``auto_process_ngs-0.6.4``.

If the Python dependencies listed above are already available in the
environment then the scripts and utilities can be run from this
directory without further installation.

Alternatively the autoprocessor can be installed into a Python
virtualenv, for example::

    virtualenv venv.auto_process
    . venv.auto_process/bin/activate
    pip install -r auto_process_ngs-0.6.4/requirements.txt
    pip install ./auto_process_ngs

In either case the installation can be configured for local preferences
by creating a configuration file as outlined in the section
:ref:`auto_process_configuration`.

.. _auto_process_configuration:

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


