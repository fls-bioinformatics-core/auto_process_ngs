auto_process_ngs
================

Scripts and utilities for automatic processing & management of NGS sequencing
data from Illumina sequencers, developed within the Bioinformatics Core Facility
(BCF) at the University of Manchester (UoM).

Full documentation is available at http://auto_process_ngs.readthedocs.org.

Overview
********

This package provides a small set of utilities to generate FASTQ files from
raw ``bcl`` files output by Illumina MiSEQ and HiSEQ sequencers, and to
perform other task such as QC runs and file management.

It also facilitates handling problem situations such as barcoding issues or
incomplete runs.

Installation
************

It is recommended to use::

    pip install .

from within the top-level source directory to install the package.

To use the package without installing it first you will need to add the
directory to your ``PYTHONPATH`` environment.

To install directly from github using ``pip``::

    pip install git+https://github.com/fls-bioinformatics-core/auto_process_ngs.git

**Note**

* For pip 1.5: you will need to specify the ``--process-dependency-links``
  argument to the ``install`` command to pull in the dependencies.
* For pip 1.6: you must first do ``pip install -r requirements.txt`` to
  pull in the dependencies.

Documentation
*************

Documentation based on ``sphinx`` is available under the ``docs`` directory.

To build::

    cd docs
    make html

which creates the documentation in the ``docs/build`` subdirectory.

Running Tests
*************

The tests can be run using::

    python setup.py test

In addition the tests are run via GitHub Actions whenever the repository
is updated:

.. image:: https://github.com/fls-bioinformatics-core/auto_process_ngs/workflows/Python%20CI/badge.svg
   :target: https://github.com/fls-bioinformatics-core/auto_process_ngs/actions?query=workflow%3A%22Python+CI%22

Dependencies
************

The package depends on the ``genomics-bcftbx`` package, available from
https://github.com/fls-bioinformatics-core/genomics

Developmental version
*********************

The developmental branch of the code on github is ``devel``, this can be
installed using::

    pip install git+https://github.com/fls-bioinformatics-core/auto_process_ngs.git@devel

Use the ``-e`` option to install an 'editable' version (see the section on
`"Editable" installs
<https://pip.pypa.io/en/latest/reference/pip_install.html#editable-installs>_`
in the pip documentation).
