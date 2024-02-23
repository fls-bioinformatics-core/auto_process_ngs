10xGenomics CellPlex and Flex datasets: ``10x_multi_config.csv``
================================================================

For 10xGenomics CellPlex (cell multiplexing) and Flex (fixed RNA
profiling) data, multiplexing analyses are run using the
`cellranger`_ ``multi`` command, provided that a
``10x_multi_config.csv`` file is also present in the project
directory.

.. note::

   Template versions of this file will be automatically and
   partially populated by the
   :doc:`setup_analysis_dirs <setup_analysis_dirs>` command
   for the appropriate single cell data type (see
   :doc:`../single_cell/10x_single_cell`).

Depending on the library type, this file should have the format
and content outlined in the following 10x Genomics documentation:

 * **CellPlex**: `cellranger_multi_cellplex`_
 * **Flex**: `cellranger_multi_frp`_
 * **Single cell immune profiling**: `cellranger_multi

.. _cellranger: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger
.. _cellranger_multi_cellplex: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/multi#cellranger-multi
.. _cellranger_multi_frp: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/multi-frp
.. _cellranger_multi_immune_profing: https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-5p-multi
