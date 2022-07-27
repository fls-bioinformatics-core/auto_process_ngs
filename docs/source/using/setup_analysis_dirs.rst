Setting up project directories using ``auto_process setup_analysis_dirs``
=========================================================================

Following Fastq generation using the ``make_fastqs`` command, typically
the next step is to create and populate the project directories ready for
subsequent QC and analysis.

This is done using the ``setup_analysis_dirs`` command, for example:

::

   auto_process.py setup_analyis_dirs

This reads the :doc:`projects.info <../control_files/projects_info>`
metadata file (initially created by ``make_fastqs``) and creates a new
subdirectory for each listed project.

Before runing ``setup_analysis_dirs``, the ``projects.info`` file should
be edited to fill in the following information for each project:

* **User**: the name(s) of the user(s) associated with the project
* **PI**: the name(s) of the principal investigator(s) (PIs) associated
  with the project
* **Library**: the library or application type (for example "RNA-seq",
  "ChIP-seq" etc)
* **Organism**: the organism(s) that the samples in the project
  originally came from (for example "Human", "Mouse", "D. Melanogaster"
  etc)
* **SC_Platform**: the single-cell platform used to prepare the samples
  (if appropriate).

See :doc:`../control_files/projects_info` for more information on the
format of the ``projects.info`` file and the allowed values for each
field.

.. note::

   ``setup_analysis_dirs`` checks that at minimum each projects has
   non-null values for the user, PI, library and organism fields. To
   skip this check, specify the ``--ignore-missing-metadata`` option.

.. note::

   In addition to the projects listed in ``projects.info``, if the
   outputs from ``make_fastqs`` included 'undetermined' Fastqs then
   these will be copied to an ``undetermined`` "project".
   
Once the project directories have been created, the next step is to
run the QC pipeline - see :doc:`Running the QC <run_qc>`.

----------------
Additional files
----------------

For certain types of data there are additional files that should
be added to the project directory manually after running
``setup_analysis_dirs``.

Note that for some of these files ``setup_analysis_dirs`` will
generate partially-populated template versions, with the
``.template`` extension:

* :doc:`10x_multiome_libraries.info <../control_files/10x_multiome_libraries_info>`
* :doc:`10x_multi_config.csv <../control_files/10x_multi_config_csv>`

These can be edited and renamed before use in downstream processing
stages (e.g. the QC pipeline).

.. _setup_analysis_dirs-add-identifier:

-----------------------------------------------
Adding an identifier to project directory names
-----------------------------------------------

Sometimes it can be useful to create multiple projects in
parallel from the same ``projects.info`` data (for example,
when a run has been processed in several different ways -
see :ref:`make_fastqs-processing-same-run-multiple-times`).

In this case the ``--id`` option of the ``setup_analysis_dirs``
command can be used to create project directories where
each name is taken from the ``projects.info`` file but has
an identifier (e.g. ``no_trimming``) appended.

For example:

::

   auto_process.py setup_analysis_dirs --id=no_trimming

would produce projects named ``<PROJECT>_no_trimming``.

When multiple ``bcl2fastq`` output directories exist in the
same analysis directory, the ``--id`` option can be paired with
the ``--unaligned-dir`` option to produce sets of projects
derived from specific ``bcl2fastq`` outputs.

For example:

::

   auto_process.py setup_analysis_dirs \
      --unaligned-dir=bcl2fastq_no_trimming --id=no_trimming
