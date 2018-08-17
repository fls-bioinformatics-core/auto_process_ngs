********
Overview
********

The core utility ``auto_process.py`` implements pipeline stages
as subcommands. For a standard sequencing run the workflow looks
like:

* :doc:`auto_process.py setup <using/setup>` initialises
  a new analysis directory for processing a sequencing
  run.

* :doc:`auto_process.py make_fastqs <using/make_fastqs>`
  performs the demultiplexing of the raw sequencer output
  (BCL files) into Fastqs files. It is a wrapper around
  Illumina's ``bcl2fastq`` (for 10xGenomics single cell
  data it wraps ``cellranger``) with additional
  functionality to fetch data, detect errors, generate
  statistics, and handle special cases.

* :doc:`auto_process.py analyse_barcodes <using/analyse_barcodes>`
  analyses the barcode sequences in the Fastqs produced
  by ``make_fastqs``. This can help for example if
  demultiplexing has failed due to mislabeling of samples.

* :doc:`auto_process.py setup_analysis_dirs <using/setup_analysis_dirs>`
  sets up individual project directories for each of the
  projects within the sequencing run.

* :doc:`auto_process.py run_qc <using/run_qc>` runs the
  QC pipeline (including ``fastqc`` and ``fastq_screen``)
  on each of the project directories and generates a report
  for each one.

* :doc:`auto_process.py publish_qc <using/publish_qc>`
  publishes the QC reports from make_fastqs and run_qc
  to a web-accessible location to be viewed by core
  facility staff and bioinformaticians.

* :doc:`auto_process.py archive <using/archive>` copies
  the analysis directory and contents from the working
  area to the archive storage area, for subsequent
  bioinformatics analysis.

Additional commands are available for reporting:

* :doc:`auto_process.py report <using/report>` generates
  reports on analysis directories in various formats.

===========
Terminology
===========

This documentation refers to various concepts such as
*analysis directories*, *analysis projects*,
*sequencing runs* and so on. We define these as
follows:

* **(Sequencing) run**: the primary sequencing data (i.e.
  Bcl files) produced by a run of a sequencer
* **Analysis directory**: the top-level directory which
  holds the outputs from ``auto_process`` when processing
  data from a sequencing run
* **Analysis project**: a subdirectory of an analysis
  directory which holds the FASTQs and QC outputs for
  a set of samples (aka a *project*) from the sequencing
  run; there can be multiple analysis projects within a
  single analysis directory

Additionally the terms *data source*, *working area*,
*QC server* and *archive storage* refer to components of
the compute infrastructure where the processing takes
place:

* **Data source**: location where the data from the
  sequencing run are located
* **Working area**: location where the processing is
  performed, and which holds the analysis directory
  and contents during the processing
* **QC server**: location where QC reports are published
  to, and which can be accessed via a web server
* **Archive storage**: location where the final outputs
  (i.e. the analysis directory and projects) are copied
  to and stored once processing is completed, for
  subsequent analysis by the bioinformaticians

These can all be on the same filesystem on a single machine;
or one or more parts can be NFS filesystems, or even
filesystems mounted on other machines which are accessed
using ``ssh``.

=============================
Supported sequencer platforms
=============================

The pipeline is currently used for output from the following
Illumina sequencers:

* HISeq 4000
* MISeq
* NextSeq
* MiniSeq

Earlier versions have been used on GAIIx and HISeq 2000/2500.

===============================
Supported single-cell platforms
===============================

The pipeline supports handling data from the Takara Bio SMARTer
ICELL8 and 10xGenomics Chromium single-call RNA-seq platforms:

* :doc:`Handling ICELL8 scRNA-seq data <icell8>`
* :doc:`Handling 10xGenomics Chromium scRNA-seq data <10xgenomics>`
