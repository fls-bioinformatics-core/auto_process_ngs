Projects metadata file: ``projects.info``
=========================================

The ``projects.info`` file is used to define the projects within
the analysis directory, and associates metadata information with
each project.

It is initially created when the :doc:`make_fastqs <../using/make_fastqs>`
command completes, and acts as a control file for subsequent
operations including:

* :doc:`setup_analysis_dirs <../using/setup_analysis_dirs>` (used to
  populate the project directories)
* :doc:`publish_qc <../using/publish_qc>` (only publishes data for projects
  that are listed in ``projects.info``)

``projects.info`` is a tab delimited file with one line for each
project, with the following fields:

===============  =================================================
Field            Data item
===============  =================================================
``Project``      Name of the project
``Samples``      List of sample names associated with the project
``User``         Name(s) of the researcher(s) directly associated
                 with the project
``Library``      Library or application type (e.g. ``RNA-seq``,
                 ``ChIP-seq``)
``SC_Platform``  Single cell platform used to prepare the samples
``Organism``     Organism(s) that the samples in the project
                 came from (e.g. ``Human``, ``Mouse``,
		 ``D. Melanogaster``)
``PI``           Name(s) of the principal investigator(s) (PIs)
                 associated with the project
``Comments``     Free text field for any additional comments
===============  =================================================

The project name and sample names are normally added automatically;
the remaining fields must be populated manually by editing the
file. The following fields are compulsory for
``setup_analysis_dirs``:

* ``User``
* ``Library``
* ``Organism``
* ``PI``

The following conventions are used for data item syntax:

* "Null" values are represented by a full stop (i.e. ``.``)
* Unknown values should be represented with a question mark
  (i.e. ``?``)
* Where there are multiple users, PIs or organisms, they should be
  separated by comma characters (e.g. ``Human,mouse``)
* Projects preceeded by a comment character (i.e. ``#``) are
  ignored

Currently there are no canonical lists of allowed values for libraries
or organism names, however:

* Organism name(s) are used to look up appropriate reference data files
  for the QC pipeline in the configuration file, so names that are used
  should be consistent with the configuration;
* Some combinations of single cell platform and library type are used
  to determine the appropriate :doc:`QC protocol <../using/run_qc>` for
  the data, in which case the single cell platform and library values
  must be a valid combination; :ref:`single_cell_and_spatial` below.

.. note::

   The organism names are converted to lower case and non-alphabetic
   characters are converted to underscores (``_``) when looking up
   names in the configuration: for example, ``Mouse`` is converted
   to ``mouse``, ``D. Melanogaster`` is converted to
   ``d_melanogaster``)

.. _single_cell_and_spatial:

Single cell and spatial data
----------------------------

For single cell and spatial data values in the library type
(``Library``) and single cell platform (``SC_platform``) fields have
special significance.

For more information on the valid options for specific types of data
see:

 * :doc:`../single_cell/10x_single_cell`
 * :doc:`../single_cell/parse`
 * :doc:`../spatial/10x_visium`

.. note::

   Use ``.`` for the single cell platform for projects that don't have
   single cell or spatial data.
