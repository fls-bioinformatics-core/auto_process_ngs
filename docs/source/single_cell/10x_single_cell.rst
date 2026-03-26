Processing 10x Genomics single cell data
========================================

Background
----------

10xGenomics provides a range of platforms for various types of single
cell samples including:

* Single cell and single nuclei RNA-seq
* Single cell ATAC-seq
* Single cell multiome gene expression and ATAC data
* Single cell immune profiling data
* Cellplex (cell multiplexing)
* Flex (fixed RNA profiling)

The ``auto_process`` pipeline has integrated support for Fastq generation
and QC (in the ``make_fastqs`` and ``run_qc`` commands) for these data
as outlined in the following sections.

Requirements
------------

The appropriate 10x Genomics software pipeline must already be installed
and available to ``auto_process.py``:

================================= =================
10x Genomics data                 Software required
================================= =================
Single cell/single nuclei RNA-seq CellRanger
Single cell ATAC-seq              CellRanger-ATAC
Single cell multiome              CellRanger-ARC
Single cell immune profiling      Cellranger
CellPlex (cell multiplexing)      CellRanger
Flex (fixed RNA profiling)        CellRanger
================================= =================

In addition:

* Fastq generation for all data requires ``bcl2fastq``;
* QC for single cell multiome GEX also uses CellRanger;
* QC for single cell multiome ATAC also uses Cellranger-ATAC.

Fastq generation
----------------

General
~~~~~~~

Fastq generation for 10x Genomics single cell data should be performed
using the :doc:`make_fastqs <../using/make_fastqs>` command with the
appropriate ``10x_*`` protocol specified via the ``--protocol`` option:

.. include:: ../auto/10x_single_cell_fq_protocols.rst

.. note::

   By default adapter trimming is automatically disabled for all
   ``10x_*`` protocols, by removing any adapter sequences specified
   in the sample sheet.

Sample sheets for 10x Genomics data can specify either Illumina index
sequences (recommended) or 10x Genomics-style index IDs (now deprecated by
10x Genomics); generally the ``10x_*`` protocols can operate with either.

.. note::

   The ``10x_multiome`` protocol (which can be used for "unpooled" 10x
   single cell multiome data) still requires 10x Genomics indexes. However
   it is recommended that even for unpooled data the specific protocol
   is used (i.e. either ``10x_multiome_gex`` or ``10x_multiome_atac``) as
   outlined in the following section.

Choosing Fastq generation protocol for single cell multiome data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are three Fastq generation protocols for single cell
multiome data; which should be used will depend on the specific
configuration of the sequencing run:

 * **Run only has either the GEX or the ATAC component of the
   single cell multiome experiment (unpooled):** use the appropriate
   ``10x_multiome_*`` protocol (i.e. ``10x_multiome_gex`` or
   ``10x_multiome_atac``) according to the type of data in the run.

   For example:

   ::

       auto_process.py make_fastqs --protocol 10x_multiome_gex

 * **Run has both GEX and ATAC components of the single cell
   multiome experiment in different lanes (pooled):** use the
   ``--lanes`` option to assign the appropriate Fastq generation
   protocol to each of the lane subsets.

   For example:

   ::

      auto_process.py make_fastqs \
      --lanes=1:10x_multiome_atac \
      --lanes=2:10x_multiome_gex


.. note::

   While the ``10x_multiome`` protocol is still supported for
   unpooled data with 10x Genomics index IDs, it is now deprecated
   and should not be used.

Analysis project setup and QC
-----------------------------

Once Fastqs have been successfully generated, the ``SC_platform``
and ``Library`` metadata items should be set to the appropriate values
for the 10x Genomics single cell project(s) in the ``projects.info``
control file.

The following values are valid options for 10x Genomics single cell
data within ``auto-process-ngs``:

.. include:: ../auto/10x_single_cell_apps.rst

Library types and extensions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Extensions are optional additions to the base library type which
give additional information about the type of experiment that was
performed.

One or more extensions can be specified as part of the library type,
by appending them using a plus symbol (``+``). For example:

* ``scRNA-seq+CSP`` indicates single cell RNA-seq (gene expression)
  performed with cell surface proteins;
* ``Immune Profiling+CSP+VDJ`` indicates single cell immune profiling
  with both cell surface proteins and V(D)J.

Running the :doc:`setup_analysis_dirs <../using/setup_analysis_dirs>`
command will automatically transfer these values into the single cell
project metadata on creation.

Additionally for certain types of data, ``setup_analysis_dirs`` will
also create template control files for use in subsequent QC runs:

 * **Single cell multiome**: a template
   :doc:`10x_multiome_libraries.info <../control_files/10x_multiome_libraries_info>`
   file, which should be renamed and populated in order to link each
   ATAC (or GEX) sample to the complementary GEX (or ATAC) sample.

 * **CellPlex and Flex**: a template
   :doc:`10x_multi_config.csv <../control_files/10x_multi_config_csv>`
   file, which should be copied, renamed and appropriately edited
   for each physical sample in the dataset.

 * **Single Cell immune profiling**: a template
   :doc:`10x_multi_config.csv <../control_files/10x_multi_config_csv>`
   file, which should be copied, renamed and appropriately edited
   for each physical sample in the dataset.

For the ``10x_multi_config.csv`` the extensions will determine the form
of the template (so that only the necessary elements are included).

The :doc:`run_qc <../using/run_qc>` command
will then determine the appropriate QC protocol to use based on the
metadata values.

Troubleshooting
---------------

Single-library analyses fail for low read counts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It has been observed that when the Fastq files produced by the ``mkfastq``
command have very low read counts then the single-library analyses may
fail, with ``cellranger count`` reporting an error of the form e.g.:

::

    Could not auto-detect Single Cell 3' chemistry. Fraction of barcodes
    on whitelist was at best 0.23%, while we expected at least 10.00% for
    one of the chemistries.

There is currently no workaround for this issue.

Single-library analyses fail to detect chemistry automatically
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default ``cellranger count`` attempts to determine the chemistry used
automatically, however this may fail if a low number of reads map to the
reference genome and give an error of the form:

::

    The chemistry was unable to be automatically determined. This can
    happen if not enough reads originate from the given reference. Please
    verify your choice of reference or explicitly specify the chemistry
    via the --chemistry argument.

If the reference data being used is correct then use the ``--chemistry``
option to specify the appropriate assay configuration - see
https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count

Appendices
----------

.. _10x_sc-immune-profiling-data:

Manual QC steps for single cell immune profiling data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Currently a full automated QC protocol is not available for Chromium
5' single cell immune profiling: specifically, there is no provision
for running Cellranger's ``multi`` pipeline for each sample, or for
automatically integrating the resulting outputs into the QC report.

It is possible to run the ``multi`` pipeline manually for each sample,
using the sample-specific ``10x_multi_config.<SAMPLE>.csv`` files.

For example, a script of the form:

::

   #!/usr/bin/bash
   #$ -N cellranger_multi_PJB01
   #$ -V
   #$ -cwd
   #$ -j y
   #$ -pe smp.pe 16
   #$ -l mem256
   mkdir -p cellranger_multi && cd cellranger_multi
   /PATH/TO/cellranger multi \
   --id PJB01 --csv PATH/TO/10x_multi_config.PJB01.csv \
   --jobmode=local \
   --localcores=16 \
   --localmem=128 \
   --maxjobs=24 \
   --jobinterval=100

could be used to submit a Cellranger ``multi`` job for the ``PJB01``
sample, with the outputs being created in a subdirectory
``cellranger_multi/PJB01`` in the current directory.

To include the outputs in the QC report, copy the relevant files
(specifically the ``web_summary.html`` files for each sample) into
the QC directory and then create an ``extra_outputs.tsv`` which
references these (as described in
:ref:`run_qc_including_external_outputs`).

For example:

::

   cellranger_multi/PJB01/web_summary.html    CellRanger multi output for PJB01

Rerunning ``run_qc`` will force update of the QC report which should
then also link in these additional reports.
