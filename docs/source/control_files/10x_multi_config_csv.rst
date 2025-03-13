10x Genomics CellPlex and Flex datasets: ``10x_multi_config[.SAMPLE].csv``
==========================================================================

For certain types of 10x Genomics data, the QC pipeline will
automatically run the `cellranger`_ ``multi`` command to perform
appropriate analyses (for example, multiplexing analyses for
CellPlex data), provided that:

 1. The ``cellranger_multi`` QC module is specified as part of
    the protocol, and
 2. one or more configuration files for running CellRanger ``multi``
    are present in the analysis project directory.

The configuration files should be named either
``10x_multi_config.csv`` (appropriate if there is a single
physical sample in the dataset) or ``10x_multi_config.<SAMPLE>.csv``
(the more general form when there is one or more physical
sample, in which case ``<SAMPLE>`` should be a label for each
physical sample in the dataset).

.. note::

   Currently the physical sample name in the configuration
   file name is arbitrary and doesn't have to match the sample
   names for the Fastq files.
    
A template version of the file (called
``10x_multi_config.csv.template``) is automatically
generated and partially populated by the
:doc:`setup_analysis_dirs <setup_analysis_dirs>` command for the
appropriate single cell data type (see
:doc:`../single_cell/10x_single_cell`). The template file is
ignored by the QC pipeline, so it must be manually copied,
renamed and edited for each physical sample that the ``multi``
command should be run on.

Specifically, at a minimum:

 1. The physical sample name should be added into the file name
    (if there are multiple physical samples);
 2. The ``[libraries]`` section should be edited to remove any
    Fastqs associated with other physical samples, and the
    ``feature_types`` field for each Fastq ID should be updated
    to assign the correct feature type (e.g. ``Gene expression``,
    ``Multiplexing capture`` etc;
 3. The ``samples`` section should be removed or edited as
    appropriate for the type of data being analysed (for example,
    assigning CMOs tp multiplexed sample names for CellPlex
    dataset).

Any additional options (for example, forcing numbers of cells)
can be included if required - details of the format of the
``multi`` configuration file for different applications can be
found via
`Cell Ranger Multi Config CSV <https://www.10xgenomics.com/support/software/cell-ranger/latest/advanced/cr-multi-config-csv-opts>`_

The QC pipeline will attempt to run ``cellranger multi`` for
each of the appropriately named configuration files in the project
directory. Currently this is supported for the following library
types:

 * `CellPlex <https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/multi#cellranger-multicellranger_multi_cellplex>`_
 * `Flex <https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/multi-frp>`_

The following library types are not yet supported:

 * `Single cell immune profiling <https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-5p-multi>`_

.. _cellranger: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger
