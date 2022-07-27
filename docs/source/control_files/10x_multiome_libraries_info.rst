10xGenomics single cell multiome linked samples: ``10x_multiome_libraries.info``
================================================================================

For a run with 10xGenomics single cell multiome gene expression or
ATAC data, an additional file ``10x_multiome_libraries.info`` can
be created which indicates where the complementary ATAC or GEX
data are located.

If this file exists then it will be used by :doc:`run_qc <../using/run_qc>`
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
