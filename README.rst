auto_process
============

Scripts and utilities for automatic processing & management of NGS sequencing
data.

Installation
************

It is recommended to use::

    pip install .

from within the top-level source directory to install the package.

To use the package without installing it first you will need to add the
directory to your `PYTHONPATH` environment.

To install directly from github using `pip`::

    pip install git+https://github.com/fls-bioinformatics-core/genomics.git@devel

Documentation
*************

Documentation based on `sphinx` is available under the `docs` directory.

To build::

    cd docs
    make html

which creates the documentation in the `docs/build` subdirectory.

Running Tests
*************

The tests can be run using::

    python setup.py test

Dependencies
************

The package depends on the ``genomics-bcftbx`` package, available from
https://github.com/fls-bioinformatics-core/genomics
