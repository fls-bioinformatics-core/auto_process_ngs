Generating QC reports using ``auto_process run_qc``
===================================================

The ``run_qc`` command is used to run the QC pipeline on the
Fastqs for each project within the analysis directory, and
generate a report for each one.

The general invocation of the command is:

::

   auto_process.py run_qc *options* [ANALYSIS_DIR]

There is also a standalone utility called ``run_qc.py`` which
can be used to run the QC pipeline on an arbitrary subdirectory.

The QC pipeline protocol used for each project will differ slightly
depending on the nature of the data within that project:

===================== ==========================
QC protocol           Used for
===================== ==========================
``standardPE``        Standard paired-end data i.e. R1/R2 Fastq pairs
``standardSE``        Standard single-end data i.e. R1 Fastqs only
``10x_scATAC``        10xGenomics single cell & single nuclei ATAC-seq
``10x_scRNAseq``      10xGenomics single cell RNA-seq
``10x_snRNAseq``      10xGenomics single nuclei RNA-seq
``10x_Visium``        10xGenomics Visium spatial RNA-seq
``10x_Multiome_ATAC`` 10xGenomics Visium multiome ATAC-seq data
``10x_Multiome_GEX``  10xGenomics Visium multiome gene expression data
``singlecell``        ICELL8 single cell RNA-seq
``ICELL8_scATAC``     ICELL8 single cell ATAC-seq
===================== ==========================

The protocol is determined automatically for each project.

The standard QC pipelines run the following external packages for
each Fastq file in the project:

 * `fastqc`_ (for general quality metrics)
 * `fastq_screen`_ for a set of genome indexes (to confirm that
   sequences are from the expected organisms, and check for any
   contamination)
 * `fastq_strand`_ is run either pair-wise (for paired end data),
   or on R1 (for single-end data), to determine the strandedness
   of the sequence data

For the single-cell and single-nuclei RNA-seq data from ICELL8 and
10xGenomics Chromium platforms, for 10xGenomics Visium RNA-seq, and
for 10xGenomics Multiome GEX data:
data:

 * `fastqc`_ is run for each Fastq file
 * `fastq_screen`_ and `fastq_strand`_ are run on just the R2
   reads

.. _fastqc:  http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
.. _fastq_screen: http://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/
.. _fastq_strand: https://genomics-bcftbx.readthedocs.io/en/latest/reference/qc_pipeline.html#fastq-strand

In addition for 10xGenomics scRNA-seq, snRNA-seq and scATAC data:

 * Single library analysis is run using the ``count`` command of
   either `cellranger`_ or `cellranger_atac`_ (as appropriate)

.. _cellranger: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger
.. _cellranger_atac: https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/what-is-cell-ranger-atac

`multiQC`_ is also run to summarise the QC from all the Fastqs in the
project.

On successful completion of the pipeline for an HTML report is
generated for each project; these are described in
:doc:`QC reports <../output/qc_reports>`. If a QC server has been
specified in the configuration then the reports can be copied
there for sharing using the :doc:`publish_qc command <publish_qc>`.

.. note::

   The QC pipeline can be run outside of the ``auto_process``
   pipeline by using the ``run_qc.py`` utility; see the
   section on :doc:`running the QC standalone <run_qc_standalone>`.

.. _multiqc: http://multiqc.info/

---------------------------
Configuring the QC pipeline
---------------------------

See :ref:`software_dependencies` for details of the additional
software required to run the QC pipeline. Environment modules can be
used to set up the runtime environment for the pipeline (see
:ref:`environment_modules`).

Suitable :ref:`job runners <job_runners>` should be defined,
particularly if running the pipeline on a compute cluster (see
:ref:`running_on_compute_cluster`).

Some of the pipeline stages also require appropriate reference
data to be set up before they can run; see the :ref:`reference_data`
configuration documentation for more details.
