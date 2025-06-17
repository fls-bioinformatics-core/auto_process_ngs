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

If a sample sheet with the appropriate 10x Genomics indexes is provided
then all 10x Genomics single cell data should be processed using the
:doc:`make_fastqs <../using/make_fastqs>` command, with the appropriate
``10x_*`` specified via the ``--protocol`` option:

=================================================== =====================
10x Genomics data                                   Protocol
=================================================== =====================
Single cell/single nuclei RNA-seq                   ``10x_chromium_sc``
Single cell ATAC-seq                                ``10x_atac``
Single cell multiome (unpooled GEX or ATAC)         ``10x_multiome``
Single cell multiome GEX (pooled GEX and ATAC)      ``10x_multiome_gex``
Single cell multiome ATAC (pooled GEX and ATAC)     ``10x_multiome_atac``
Single cell immune profiling                        ``10x_chromium_sc``
CellPlex (cell multiplexing)                        ``10x_chromium_sc``
Flex (fixed RNA profiling)                          ``10x_chromium_sc``
=================================================== =====================

.. note::

   By default adapter trimming is automatically disabled for all
   ``10x_*`` protocols, by removing any adapter sequences specified
   in the sample sheet.

.. note::

   If the sample sheet contains Illumina index sequences then the
   ``standard`` protocol should be used instead (note that in this case
   the defaults used for masking and trimming compared to the defaults
   may differ from those used by the 10x Genomics pipeline).

Choosing Fastq generation protocol for single cell multiome data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are three Fastq generation protocols for single cell
multiome data; which should be used will depend on the specific
configuration of the sequencing run:

 * **Run only has either the GEX or the ATAC component of the single
   cell multiome experiment:** the ``10x_multiome`` protocol is
   preferred as Cellranger-ARC should be able to automatically
   determine which component the data are.

 * **Run has both GEX and ATAC components of the single cell
   multiome experiment in different lanes:** in this situation
   Cellranger-ARC cannot automatically determine which component
   is in which lane, so the ``10x_multiome_gex`` protocol should be
   explicitly specified for the lanes with the GEX data, and the
   ``10x_multiome_atac`` specified for those with the ATAC data,
   via the ``--lanes`` option.

   For example:

   ::

      auto_process.py make_fastqs \
      --lanes=1:10x_multiome_atac \
      --lanes=2:10x_multiome_gex

   See :ref:`10x_multiome-pooled-data` for more information on
   how the multiome protocols are implemented and used.

Analysis project setup and QC
-----------------------------

Once Fastqs have been successfully generated, the ``SC_platform``
and ``Library`` metadata items should be set to the appropriate values
for the 10x Genomics single cell project(s) in the ``projects.info``
control file.

The following values are valid options for 10x Genomics single cell
data within ``auto-process-ngs``:

========================================= ==============================
Single cell platform                      Library types
========================================= ==============================
``10xGenomics Chromium 3'[*]``            ``scRNA-seq``, ``snRNA-seq``,
                                          ``CellPlex``,
                                          ``CellPlex scRNA-seq``, ``Flex``
``10xGenomics Chromium 3'v3``             ``scRNA-seq``, ``snRNA-seq``
                                          ``CellPlex``,
                                          ``CellPlex scRNA-seq``, ``Flex``
``10xGenomics Chromium 3'v2``             ``scRNA-seq``, ``snRNA-seq``
                                          ``CellPlex``,
                                          ``CellPlex scRNA-seq``, ``Flex``
``10xGenomics Chromium GEM-X``            ``Flex``
``10xGenomics Chromium GEM-X 3'v4``       ``scRNA-seq``
``10xGenomics Chromium Next GEM``         ``scRNA-seq``,
                                          ``CellPlex scRNA-seq``, ``Flex``
``10xGenomics Chromium Next GEM 3'v3.1``  ``scRNA-seq``,
                                          ``CellPlex scRNA-seq``
``10xGenomics Chromium 5'``               ``Single Cell Immune Profiling``
``10xGenomics Single Cell ATAC``          ``scATAC-seq``, ``snATAC-seq``
``10xGenomics Single Cell Multiome``      ``ATAC``, ``GEX``
========================================= ==============================

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
   file, which should be copied, renamed and appropriatedly edited
   for each physical sample in the dataset.

 * **Single Cell immune profiling**: a template
   :doc:`10x_multi_config.csv <../control_files/10x_multi_config_csv>`
   file, which should be copied, renamed and appropriatedly edited
   for each physical sample in the dataset.

The :doc:`run_qc <../using/run_qc>` command
will then determine the appropriate QC protocol to use based on the
metadata values.

.. note::

   Currently a full QC pipeline is not implemented for single cell
   immune profiling data: see :ref:`10x_sc-immune-profiling-data`
   for additional manual steps that can be performed for these types
   of data.

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

.. _10x_multiome-pooled-data:

Details for handling pooled single cell multiome ATAC and GEX data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If 10x Genomics single cell multiome ATAC and multiome GEX libraries
are sequenced together in the same run then the standard ``10x_multiome``
protocol of the ``make_fastqs`` command is unable to correctly process
the data.

Pooling the ATAC and GEX components of a single cell multiome experiment
is not officially supported by 10x Genomics, and this limitation is due
to this configuration not being supported by the ``cellranger-arc``
pipeline. However they do provide information on how to handle this
situation in this knowledge base article:

https://kb.10xgenomics.com/hc/en-us/articles/360049373331-Can-Multiome-ATAC-and-Multiome-GEX-libraries-be-sequenced-together-

and the two sub-protocols outlined in that article have been implemented
within ``make_fastqs`` as the ``10x_multiome_atac`` and ``10_multiome_gex``
protocols, which should be used as follows:

 1. Ensure that ATAC and GEX data are assigned to separate projects
    in the input sample sheet
 2. Use the ``--lanes`` option to explicitly specify the appropriate
    sub-protocol for the lanes with the ATAC and GEX samples

For example:

::

   auto_process.py make_fastqs \
      --lanes=1:10x_multiome_atac \
      --lanes=2:10x_multiome_gex

assuming that the ATAC data are in lane 1 and the GEX data in lane 2.

.. warning::

   These protocols should only be used when the single cell
   multiome data has been pooled with other types of data;
   when the single cell multiome data for a single component
   (either GEX or ATAC) comprises the whole sequencing run
   then the ``10x_multiome`` protocol should be used instead.

The ``10x_multiome_atac`` protocol then runs ``cellranger-arc mkfastq``
with the following custom options:

 1. ``--use-bases-mask`` with a bases mask string that has been
    adjusted appropriately to match the template ``Y*,I8n*,Y24,Y*``
 2. ``--filter-single-index`` is explicitly specified

The ``10x_multiome_gex`` protocol runs ``cellranger-arc mkfastq`` with
the following custom options:

 1. ``--use-bases-mask`` with a bases mask string that has been
    adjusted appropriately to match the template
    ``Y28n*,I10,I10n*,Y*``
 2. ``--filter-dual-index`` is explicitly specified
