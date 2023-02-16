10xGenomics CellPlex and Flex datasets: ``10x_multi_config.csv``
================================================================

For 10xGenomics CellPlex (cell multiplexing) and Flex (fixed RNA
profiling) data, multiplexing analyses are run using the
`cellranger`_ ``multi`` command, provided that a
``10x_multi_config.csv`` file is also present in the project
directory.

Depending on the library type, this file should have the format
and content outlined at `cellranger_multi_cellplex`_ (for CellPlex
data) or `cellranger_multi_frp`_ (for Flex data).

.. _cellranger: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger
.. _cellranger_multi_cellplex: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/multi#cellranger-multi
.. _cellranger_multi_frp: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/multi-frp
