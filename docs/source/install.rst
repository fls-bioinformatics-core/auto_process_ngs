
.. _auto_process_installation:

************
Installation
************

The ``auto-process-ngs`` package is available from its GitHub respository at

 * https://github.com/fls-bioinformatics-core/auto_process_ngs

Specific versions can be obtained as ``tar.gz`` archives from:

 * https://github.com/fls-bioinformatics-core/auto_process_ngs/releases

The software is written in Python and works with Python 2.7.

It is recommended to install the package into a Python ``virtualenv``, for
example::

    virtualenv venv.ap
    . venv.ap/bin/activate

To install a specific version, first download and unpack the source code,
e.g.::

    wget https://github.com/fls-bioinformatics-core/auto_process_ngs/archive/0.10.0.tar.gz
    tar zxf 0.10.0.tar.gz

Then install the requirements and then rest of the packages::

    pip install -r auto_process_ngs-0.10.0/requirements.txt
    pip install ./auto_process_ngs-0.10.0

In addition a number of external software packages are required for the
Fastq generation and QC pipelines; these are listed in the section
:ref:`auto_process_requirements`.

The installation may also require some configuration to tailor it to the
local environment (for example, if running on a compute cluster); this is
outlined in more detail in the
:doc:`configuration documentation <configuration>`.
