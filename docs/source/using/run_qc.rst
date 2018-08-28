Generating QC reports using ``auto_process setup``
==================================================

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

============== ========================
QC protocol    Used for
============== ========================
``standardPE`` Standard paired-end data i.e. R1/R2 Fastq pairs
``standardSE`` Standard single-end data i.e. R1 Fastqs only
``singlecell`` Data from single-cell platforms (ICELL8, 10xGenomics)
============== ========================

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

For the single-cell data from ICELL8 and 10xGenomics Chromium
platforms:

 * `fastqc`_ is run for each Fastq file
 * `fastq_screen`_ and `fastq_strand`_ are run on just the R2
   reads

.. _fastqc:  http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
.. _fastq_screen: http://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/
.. _fastq_strand: https://genomics-bcftbx.readthedocs.io/en/latest/reference/qc_pipeline.html#fastq-strand

`multiQC`_ is also run to summarise the QC from all the Fastqs in the
project.

On successful completion of the pipeline for an HTML report is
generated for each project; these are described in
:doc:`QC reports <../output/qc_reports>`. If a QC server has been
specified in the configuration then the reports can be copied
there for sharing using the :doc:`publish_qc command <publish_qc>`.

.. _multiqc: http://multiqc.info/
