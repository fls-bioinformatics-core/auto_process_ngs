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

----------------
Additional files
----------------

For certain types of data there are additional files that should
be added to the project directory manually after running
``setup_analysis_dirs``.

Note that for some of these files ``setup_analysis_dirs`` will
generate partially-populated template versions, with the
``.template`` extension. These can be edited and renamed before
use in downstream processing stages (e.g. the QC pipeline).

.. _10x_multiome_libraries_info:

10xGenomics single cell multiome linked samples
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For a run with 10xGenomics single cell multiome gene expression or
ATAC data, an additional file ``10x_multiome_libraries.info`` can
be created which indicates where the complementary ATAC or GEX
data are located.

If this file exists then it will be used by :doc:`run_qc <run_qc>`
to run the ``cellranger-arc`` single library analysis for each sample
listed in the file (if it is absent then the single library analysis
is skipped).

The file should consist of one line for each sample in the project,
with two fields separated by a tab:

::

   <SAMPLE>     [RUN][:PROJECT][/COMPLEMENTARY_SAMPLE]

For example:

::

   AA_multiome_ATAC       201031_NB0123_000011_AHXXXX:AA/AA_multiome_GEX

.. note::

   ``RUN`` can be specified as a full path to an analysis
   directory (e.g. ``/data/201031_NB0123_000011_AHXXXX_analysis``),
   the name of an analysis directory (e.g.
   ``201031_NB0123_000011_AHXXXX_analysis``), the name of a run
   (e.g. ``201031_NB0123_000011_AHXXXX``) or a run reference ID
   (e.g. ``NEXTSEQ_201031#11``).

   If a path is not supplied then the parent directory
   structure will be scanned to try and locate the specified
   analysis directory.

.. note::

   Lines starting with a comment character (``#``) will be
   ignored when setting up the single library analysis.
