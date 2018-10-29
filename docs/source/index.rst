**************************************************
auto-process-ngs: automated processing of NGS data
**************************************************

``auto_process_ngs`` provides a set of utilities which automate
the processing of sequencing data from Illumina Next Generation
Sequencing (NGS) platforms, specifically:

* Generation of Fastq files from raw Bcl data produced by
  the sequencer
* Dividing Fastqs into projects for subsequent analysis
* Running quality checks (QC) on each project
* Copying final data to an archive location ready for
  appropriate bioinformatic analyses

In addition to standard Illumina sequencing data it can also
handle data prepared using a number of single-cell (SC) RNA-seq
platforms.

Together these utilities form the pipeline used for the initial
processing, QC and management of sequencing data within the
`Bioinformatics Core Facility <https://www.bmh.manchester.ac.uk/research/facilities/bioinformatics/>`_
at the `University of Manchester <https://www.manchester.ac.uk/>`_.

.. toctree::
   :maxdepth: 2
   :caption: Getting started

   overview
   requirements
   install
   configuration

.. _running-pipelines:

.. toctree::
   :maxdepth: 2
   :caption: Pipeline stages

   Starting an analysis <using/setup>
   Fastq generation <using/make_fastqs>
   Analysing barcodes <using/analyse_barcodes>
   Setting up projects <using/setup_analysis_dirs>
   Running QC <using/run_qc>
   Publishing QC <using/publish_qc>
   Archiving analyses <using/archive>
   Troubleshooting <using/troubleshooting>

.. toctree::
   :maxdepth: 2
   :caption: Post-processing

   Reporting analyses <using/report>
   Managing and sharing data <using/managing_data>

.. toctree::
   :maxdepth: 2
   :caption: Outputs
	     
   Analysis and project directories <output/analysis_dirs>
   Processing QC <output/processing_qc>
   QC reports <output/qc_reports>

.. _single-cell-docs:

.. toctree::
   :maxdepth: 2
   :caption: Single-cell Platforms

   Takara Bio ICELL8 data <icell8>
   10xGenomics Chromium data  <10xgenomics>

.. _reference-docs:

.. toctree::
   :maxdepth: 2
   :caption: Reference Documentation

   reference/commands
   reference/utilities

.. _developers-docs:

.. toctree::
   :maxdepth: 2
   :caption: Developer Documentation

   advanced
   developers/index
