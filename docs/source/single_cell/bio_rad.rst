Processing Bio-Rad ddSEQ single cell data
=========================================

Background
----------

Bio-Rad provides ddSEQ Single-Cell 3' RNA-seq and ddSEQ Single-Cell
ATAC-seq kits. This document outlines using the ``auto-process-ngs``
pipeline to perform Fastq generation and QC for these data.

Requirements
------------

No additional external software is required for the Fastq generation
or QC.

Fastq generation
----------------

The ``biorad_ddseq`` Fastq generation protocol should be used for
Bio-Rad's ddSEQ Single-Cell RNA-Seq and ATAC-seq data when running
the :doc:`make_fastqs <../using/make_fastqs>` command.


Analysis project setup and QC
-----------------------------

Once Fastqs have been successfully generated, the ``SC_platform``
and ``Library`` metadata items should be set to the appropriate values
for the Parse Evercode project(s) in the ``projects.info`` control file.

The following values are valid options:

.. include:: ../auto/bio_rad_single_cell_apps.rst

Running the :doc:`setup_analysis_dirs <../using/setup_analysis_dirs>`
command will automatically transfer these values into the project
metadata on creation. The :doc:`run_qc <../using/run_qc>` command
will then determine the appropriate QC protocol to use based on the
metadata values.

Appendix: position of DNA insert sequences for QC (ddSEQ ATAC)
--------------------------------------------------------------

The Bio-Rad SureCell ATAC library configuration for R1 reads is assumed
to look like:

* 7bp barcode
* 0-4bp "phase block"
* Fixed sequence ``TATGCATGAC`` (10bp)
* 7bp barcode
* Fixed sequence ``AGTCACTGAG`` (10bp)
* 7bp barcode
* Fixed sequence ``TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG`` (33bp)
* DNA insert (remainder of the sequence, 40-44bp)

The DNA insert in any individual read sequence is therefore expected
to start between positions 75 and 79 (depending on the size of the
variable length phase block sequence).

For the R2 reads it is assumed to look like:

* DNA insert (40bp)
* 6bp TI adapter

Based on this, for the ``BioRad_ddSEQ_ATAC`` QC protocol the mapped
metrics are calculated using R1 positions 79-118 and all of R2.
