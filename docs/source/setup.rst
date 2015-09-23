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

The autoprocessor reads settings for the local system from a `settings.ini`
file. An initial `settings.ini` file is created the first time that the
autoprocessor is run.

Job runner specification
------------------------

Job runners tell the autoprocessor how to run programs. There are
currently only two available:

* SimpleJobRunner: runs programs as a subprocess of the current process
* GEJobRunner: runs programs using Grid Engine (GE)

The `GEJobRunner` is recommended when using the autoprocessor on cluster
systems. To specify additional Grid Engine-specific options to use with
the runner, enclose them in parentheses e.g.::

    [runners]
    bcl2fastq = GEJobRunner(-pe smp.pe 8)

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


