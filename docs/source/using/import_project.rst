Copying analysis projects using ``auto_process import_project``
===============================================================

Under certain circumstances it may be necessary to import one or more
analysis projects from one analysis directory into another (for example,
when the raw sequencing data has been reprocessed to generate new Fastqs
using different parameters).

While the required steps can be performed manually, in these cases
it is recommended that the ``import_project`` command is used as this
automates the process.

The basic usage is:

::

   auto_process.py import_project ANALYSIS_DIR PATH/TO/PROJECT_DIR

(If the ``ANALYSIS_DIR`` is ommitted then the project is imported into
the current analysis directory.)

The following steps are peformed:

* The project directory and its contents are copied to the target
  analysis directory;
* The details of the imported project is added to the ``projects.info``
  file;
* The metadata of the imported project is updated to ensure that any
  stored paths are consistent with the new location;
* QC reports and ZIP archives are automatically regenerated.

.. note::

   If a project directory already exists in the target analysis
   directory with the same name as the one being imported, then the
   ``import_project`` command will stop with a warning.

The command supports the following additional options:

* ``--comment``: allows comment text to be specified which will be
  appended to any existing comments associated with the project.
