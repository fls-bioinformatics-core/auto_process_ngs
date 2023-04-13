Running the QC standalone with ``run_qc.py``
============================================

The utility ``run_qc.py`` allows the QC pipeline to be run on an
arbitrary set of Fastqs outside of the ``auto_process`` pipeline.

The general invocation is:

::

   run_qc.py DIR | FASTQ [ FASTQ ... ]

If ``DIR`` is supplied then it should be a directory containing
the Fastq files to run the QC on:

::

   run_qc.py /mnt/data/project/fastqs.trimmed/

Alternatively a list of Fastq files can be supplied directly,
for example:

::

   run_qc.py /mnt/data/project/fastqs.trimmed/*.fastq

.. note::

   You can use more elaborate shell wildcard patterns to select
   the Fastqs to pass to the QC, for example:

   ::

      run_qc.py /mnt/data/project/fastqs.trimmed/{SC1_*,SC2_*}.trimmed.fastq

   would only select Fastq files matching ``SC1_*.trimmed.fastq``
   and ``SC2_*.trimmed.fastq``.

   See https://tldp.org/LDP/GNU-Linux-Tools-Summary/html/x11655.htm
   for more detail on using shell-style wildcards for pattern
   matching.

Various options are available to control the QC; the following
sections outline the most useful - see the documentation in the
:ref:`utilities_run_qc` section for the full set of options.

Specifying the QC metadata
--------------------------

The following options specify metadata for the QC which will
determine which metrics are run:

* ``--protocol``: specify the QC protocol (see :doc:`run_qc`
  for a complete list, or use the ``--info`` option)
* ``--organism``: specify the organism(s) (e.g. ``human``,
  ``mouse`` etc; use the ``--info`` option to list the
  organisms available in the current configuration)
* ``--library-type``: specify the library type (e.g. ``RNA-seq``,
  ``ChIP-seq`` etc)
* ``--single-cell-platform``: specify the name for the single
  cell platform (e.g. ``10xGenomics Chromium 3'v3``; see
  :doc:`../control_files/projects_info` for a complete list,
  or use the ``--info`` option)

.. note::

   By default ``run_qc.py`` attempts to locate a project
   metadata file associated with the input files, and will
   use the metadata defined there unless overridden by the
   options above; use the ``--ignore-metadata`` option to
   ignore the project metadata even if it is found.

Specifying the outputs
----------------------

By default the QC reports are written to the current directory,
with the QC outputs written to a ``qc`` subdirectory.

The following options can be used to override the defaults:

* ``-o``/``--out_dir``: sets the top-level directory where
  the QC reports and outputs are written
* ``--qc_dir``: sets the location where the QC outputs are
  written; if this is a relative path then it will be a
  subdirectory of the top-level output directory
* ``--filename``: name for the HTML report from the QC
* ``--name``: sets the name for the project (used in the
  QC report title)

Updating existing QC outputs
----------------------------

If the output directory for the QC already exists then by
default ``run_qc.py`` will stop without running the pipeline;
to override this behaviour and update an existing QC output
directory, specify the ``-u`` or ``--update`` option.

Specifying reference data
-------------------------

By default the reference data used in the QC pipeline (for
example the ``STAR`` index or the genome annotation) are
determined from the configuration file according to the
organism associated with the dataset.

If no matching reference data are found then the following
options are available:

* ``--star-index``: sets the path to the ``STAR`` index
  to use
* ``--gtf``: sets the path to the genome annotation in
  GTF format

Note that these options will also override the reference data
which are set in the configuration.

.. note::

   The ``--info`` option can be used to see what default
   reference data has been assigned to each configured
   organism name.

Running 10xGenomics single library analyses
-------------------------------------------

When running the QC pipeline on 10xGenomics data using an
:doc:`appropriate protocol <run_qc>`, additional options are
available for controlling the single library analyses that
are run using the `count` command of the appropriate
10xGenomics software package.

The following options can be used to override the defaults
defined in the configuration:

* ``--cellranger``: explicitly sets the path to the ``cellranger``
  (or other appropriate 10xGenomics package)
* ``--cellranger-reference``: sets the path to the reference
  dataset to use for single library analysis

For single cell RNA-seq additional options are available:

* ``--10x_force_cells``: explicitly specify the number of cells,
  overriding automatic cell detection algorithms (i.e. set the
  ``--force-cells`` option for Cellranger)
* ``--10x_chemistry``: explicitly set the assay configuration,
  overriding the automatic assay detection (i.e. set the
  ``--chemistry`` option for Cellranger)

.. note::

   The full outputs from the single library analyses are
   copied to the ``cellranger_count`` subdirectory of the
   top-level output directory, in addition to the metrics
   and HTML summary files written to the QC directory.

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

Explicitly specifying these parameters for a QC run submitted as
a single job on a Grid Engine compute cluster might look like:

::

   qsub -b y -V -pe smp.pe 16 'run_qc.py --local --maxcores=16 --maxmem=64 /data/Fastqs'

Listing organisms and other information: ``--info``
---------------------------------------------------

The ``--info`` option of ``run_qc.py`` displays various items
from the current configuration, including a list of the
available organisms and the indexes and other reference data
assigned to each.

Other information includes the available QC protocols, single
ell platforms and FastqScreen ``.conf`` files.

Once the information is displayed, ``run_qc.py`` will exit
without performing any further action.
