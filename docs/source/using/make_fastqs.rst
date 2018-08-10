Fastq generation using ``auto_process make_fastqs``
===================================================

Overview
--------

The ``make_fastqs`` command is the backbone of the ``auto_process``
pipeline. It performs the key step of generating Fastq files from
the raw BCL data produced by the sequencer.

Standard Fastq generation
-------------------------

After creating the analysis directory using the :doc:`setup` command,
the Fastq generation pipeline is performed using a command of the
form:

::

   auto_process.py make_fastqs *options* [ANALYSIS_DIR]

.. _icell8_fastq_generation:

Fastq generation for ICELL8 single-cell data
--------------------------------------------

Initial Fastqs can be generated from ICELL8 single-cell8 data using the
``--protocol=icell8`` option:

::

    auto_process.py make_fastqs --protocol=icell8 ...

Subsequently the read pairs can be processed using the
``process_icell8.py`` utility described in the
:ref:`icell8_qc_and_filtering_protocol` section.

.. note::

   ``--protocol=icell8`` runs the standard bcl2fastq commands with
   with the following settings:

   * Disable adapter trimming and masking by setting
     ``--minimum-trimmed-read-length=21`` and
     ``--mask-short-adapter-reads=0`` (recommended by Wafergen
     specifically for NextSeq data)
   * Updating the bases mask setting so that only the first 21 bases
     of the R1 read are kept.

   This is recommended to stop unintentional trimming of UMI sequences
   (which are mostly random) from the R1, should they happen to match
   part of an adapter sequence.
