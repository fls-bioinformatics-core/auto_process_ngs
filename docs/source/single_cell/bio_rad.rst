Processing Bio-Rad ddSEQ single cell data
=========================================

Background
----------

Bio-Rad provides ddSEQ Single-Cell ATAC-seq and ddSEQ Single-Cell 3'
RNA-seq kits. This document outlines using the ``auto-process-ngs``
pipeline to perform Fastq generation and QC for ddSEQ Single-Cell
ATAC-seq data (there is no support currently for ddSEQ single cell
RNA-seq data).

Requirements
------------

No additional external software is required for the Fastq generation
or QC.

Fastq generation
----------------

The ``biorad_ddseq`` Fastq generation protocol should be used for
Bio-Rad's ddSEQ Single-Cell ATAC-seq data when running the
:doc:`make_fastqs <../using/make_fastqs>` command.


Analysis project setup and QC
-----------------------------

Once Fastqs have been successfully generated, the ``SC_platform``
and ``Library`` metadata items should be set to the appropriate values
for the Parse Evercode project(s) in the ``projects.info`` control file.

The following values are valid options:

===================================== =================================
Single cell platform                  Library types
===================================== =================================
``Bio-Rad ddSEQ Single Cell ATAC``    ``scATAC-seq``, ``snATAC-seq``
===================================== =================================

Running the :doc:`setup_analysis_dirs <../using/setup_analysis_dirs>`
command will automatically transfer these values into the project
metadata on creation. The :doc:`run_qc <../using/run_qc>` command
will then determine the appropriate QC protocol to use based on the
metadata values.
