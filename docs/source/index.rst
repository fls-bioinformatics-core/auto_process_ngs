*************************************************************
auto-process-ngs: automatic processing of NGS sequencing data
*************************************************************

.. toctree::
   :maxdepth: 2
   :caption: Getting started

   requirements
   install

.. _running-pipelines:

.. toctree::
   :maxdepth: 2
   :caption: Running Pipelines

   Analysis setup <using/setup>
   Fastq generation <using/make_fastqs>
   Analysing barcodes <using/analyse_barcodes>
   Setting up projects <using/setup_analysis_dirs>
   Running QC <using/run_qc>
   Publishing QC <using/publish_qc>
   Archiving analyses <using/archive>
   commands
   barcode_analysis
   metadata
   utilities
   managing_data
   problems

.. toctree::
   :maxdepth: 2
   :caption: Outputs
	     
   analysis_dirs

.. _single-cell-docs:

.. toctree::
   :maxdepth: 2
   :caption: Single-cell Platforms

   Takara Bio ICELL8 data <icell8>
   10xGenomics Chromium data  <10xgenomics>

.. _developers-docs:

.. toctree::
   :maxdepth: 2
   :caption: Developer Documentation

   advanced
   developers/index

============
What it does
============

``auto_process_ngs`` provides a set of utilities which automate
the generation and QC of Fastqs from the sequencing output of
various Illumina Next Generation Sequencing (NGS) platforms.
It can also handle data prepared from a number of single-cell
RNA-seq (scRNA-seq) systems. Together these utilities form the
pipeline used for the initial processing, QC and management of
sequencing data within the
`Bioinformatics Core Facility <https://www.bmh.manchester.ac.uk/research/facilities/bioinformatics/>`_
at the `University of Manchester <https://www.manchester.ac.uk/>`_.

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

* :doc:`auto_process.py setup_analysis_dirs <using/setup_analysis_dirs>`
  sets up individual project directories for each of the
  projects within the sequencing run.

* :doc:`auto_process.py run_qc <using/run_qc>` runs the
  QC on the project directories. This runs a set of
  utilities including ``fastqc`` and ``fastq_screen`` to
  and generates a report for each project.

* :doc:`auto_process.py publish_qc <using/publish_qc>`
  publishes the QC reports from make_fastqs and run_qc
  to a web-accessible location to be viewed by core
  facility staff and bioinformaticians.

* :doc:`auto_process.py archive <using/archive>` copies
  the analysis directory and contents from the working
  area to the archive storage area.

Additional commands are available for special cases:

* :doc:`auto_process.py analyse_barcodes <using/analyse_barcodes>`
  analyses the barcode sequences in the Fastqs produced
  by ``make_fastqs``. This can help for example if
  demultiplexing has failed due to mislabeling of samples.

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


