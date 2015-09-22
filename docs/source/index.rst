auto_process_ngs: automatic processing of NGS sequencing data
=============================================================

`auto_process_ngs` provides a set of utilities to help with the automated
processing and management of sequencing data from Illumina's HiSEQ and
MiSEQ platforms.

Basic usage
***********

Set up the analysis project from the sequencer output directory `<SOURCE>`
(containing the bcl files)::

    auto_process.py setup <SOURCE>

Move into the resulting analysis directory and generate FASTQ files from
the bcls::

    cd <SOURCE>_analysis
    auto_process.py make_fastqs

Edit the `projects.info` file then set up the project directories and run
the QC pipeline::

    auto_process.py setup_analysis_dirs
    auto_process.py run_qc

Contents
********

.. toctree::
   :maxdepth: 2

   setup
   standard_protocol
   commands
   analysis_dirs
   metadata
   problems
   count_barcodes
   managing_data

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

