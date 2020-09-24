Setting up project directories using ``auto_process setup_analysis_dirs``
=========================================================================

Following Fastq generation using the ``make_fastqs`` command, typically
the next step is to create and populate the project directories ready for
subsequent QC and analysis.

This is done using the ``setup_analysis_dirs`` command, for example:

::

   auto_process.py setup_analyis_dirs

This reads the :doc:`projects.info <projects_info>` metadata file
(initially created by ``make_fastqs``) and creates a new subdirectory
for each listed project.

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

See :doc:`projects_info` for more information on the format of the
``projects.info`` file and the allowed values for each field.

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
