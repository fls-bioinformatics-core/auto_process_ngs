Generating QC reports using ``auto_process run_qc``
===================================================

--------
Overview
--------

The ``run_qc`` command is used to run the QC pipeline on the
Fastqs for each project within the analysis directory, and
generate a report for each one.

The general invocation of the command is:

::

   auto_process.py run_qc *options* [ANALYSIS_DIR]

There is also a standalone utility called ``run_qc.py`` which
can be used to run the QC pipeline on an arbitrary subdirectory.

------------
QC protocols
------------

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
``10x_Multiome_ATAC`` 10xGenomics single cell multiome ATAC-seq data
``10x_Multiome_GEX``  10xGenomics single cell multiome gene expression data
``10x_CellPlex``      10xGenomics CellPlex cell multiplexing data
``ParseEvercode``     Parse Evercode single cell RNA-seq
``singlecell``        ICELL8 single cell RNA-seq
``ICELL8_scATAC``     ICELL8 single cell ATAC-seq
===================== ==========================

The protocol is determined automatically for each project, based
on the metadata.

In turn each protocol defines a set of "QC modules" that are run
by the QC pipeline.

----------
QC modules
----------

The following QC modules are available within the QC pipeline:

============================== ======================
QC module                      Details
============================== ======================
``fastqc``                     Runs `fastqc`_ for general quality metrics
``fastq_screen``               Runs `fastq_screen`_ for a set of genome
                               indexes, to verify sequences correspond to
                               the expected organism (and check for
                               contaminants)
``sequence_lengths``           Examines distribution of sequence lengths
                               and checks for padding and masking
``strandedness``               Runs `fastq_strand`_ on the appropriate
                               reads to indicate the strand-specificity of
                               the sequence data (requires appropriate
			       ``STAR`` indexes)
``picard_insert_size_metrics`` Runs `picard`_ ``CollectInsertSizeMetrics``
                               on paired data to determine insert sizes
``rseqc_genebody_coverage``    Runs `rseqc`_ ``geneBody_coverage.py`` to
                               generate gene body coverage plot for all
			       samples
``qualimap_rnaseq``            Runs `qualimap`_ ``rnaseq`` to generate
                               various metrics including coverage and
			       genomic origin of reads
``cellranger_count``           Single library analysis for each sample using
                               `cellranger`_ ``count``
``cellranger-atac_count``      Single library analysis for each sample using
                               `cellranger_atac`_ ``count``
``cellranger-arc_count``       Single cell multiome analysis using
                               `cellranger_arc`_ ``count`` (requires
                               :doc:`10x_multiome_libraries.info <../control_files/10x_multiome_libraries_info>`)
``cellranger_multi``           Cell multiplexing analysis using
                               `cellranger`_ ``multi`` (requires
                               :doc:`10x_multi_config.csv <../control_files/10x_multi_config_csv>`)
============================== ======================

Appropriate reference data must be available (for example,
``STAR`` indexes or 10x Genomics reference datasets), and
certain additional files (noted in the QC module descriptions
above) may be required for some of the QC modules to run;
typically the modules are skipped when appropriate reference
data is not available.

.. _fastqc:  http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
.. _fastq_screen: http://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/
.. _fastq_strand: https://genomics-bcftbx.readthedocs.io/en/latest/reference/qc_pipeline.html#fastq-strand
.. _picard: https://gatk.broadinstitute.org/hc/en-us/articles/360037055772-CollectInsertSizeMetrics-Picard-
.. _rseqc: http://rseqc.sourceforge.net/#
.. _qualimap: http://qualimap.conesalab.org/doc_html/command_line.html#rna-seq-qc
.. _cellranger: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger
.. _cellranger_atac: https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/what-is-cell-ranger-atac
.. _cellranger_arc: https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/what-is-cell-ranger-arc
.. _multiqc: http://multiqc.info/

On successful completion of the pipeline an HTML report is
generated for each project; these are described in
:doc:`QC reports <../output/qc_reports>`. By default `multiQC`_
is also run as part of the reporting.

If a QC server has been specified in the configuration then the
reports can be copied there for sharing using the
:doc:`publish_qc command <publish_qc>`.

.. note::

   The QC pipeline can be run outside of the ``auto_process``
   pipeline by using the ``run_qc.py`` utility; see the
   section on :doc:`running the QC standalone <run_qc_standalone>`.

------------------
Additional options
------------------

For 10xGenomics datasets, the following options can be used to
override the defaults defined in the configuration:

* ``--cellranger``: explicitly sets the path to the ``cellranger``
  (or other appropriate 10xGenomics package)
* ``--10x_force_cells``: explicitly specify the number of cells,
  overriding automatic cell detection algorithms (i.e. set the
  ``--force-cells`` option for CellRanger)
* ``--10x_extra_projects``: specify additional project directories
  to fetch Fastqs from when running single library analyses (i.e.
  add the Fastq directory paths for each project to the
  ``--fastqs`` option for CellRanger)
