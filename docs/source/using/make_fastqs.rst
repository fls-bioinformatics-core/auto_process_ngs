Fastq generation using ``auto_process make_fastqs``
===================================================

Overview
--------

The ``make_fastqs`` command is the backbone of the ``auto_process``
pipeline. It performs the key step of generating Fastq files from
the raw BCL data produced by the sequencer.

The general invocation of the command is:

::

   auto_process.py make_fastqs [--protocol=PROTOCOL] *options* [ANALYSIS_DIR]

The ``--protocol`` option should be chosen according to the type
of data that is being processed:

=================== =====================================
Protocol option     Used for
=================== =====================================
``standard``        Standard Illumina sequencing data
``icell8``          ICELL8 single-cell data
``10x_chromium_sc`` 10xGenomics Chromium single-cell data
=================== =====================================

By default ``make_fastqs`` performs the following steps:

* Fetches the BCL data and copies it to the ``primary_data`` subdirectory
  of the analysis directory
* Checks the sample sheet for errors and possible barcode collisions
* Runs ``bcl2fastq`` (or ``cellranger`` for 10xGenomics Chromium data)
  and verifies that the outputs match what was expected from the input
  sample sheet
* Generates statistics for the Fastq data and produces a report on the
  analysing the numbers of reads assigned to each sample, lane and
  project (``processing.html``).

Standard Fastq generation
-------------------------

After creating the analysis directory using the :doc:`setup <setup>`
command, the Fastq generation pipeline is performed using a command
of the form:

::

   auto_process.py make_fastqs *options* [ANALYSIS_DIR]

.. _icell8_fastq_generation:

Fastq generation for ICELL8 single-cell data
--------------------------------------------

Initial Fastqs can be generated from ICELL8 single-cell8 data using the
``--protocol=icell8`` option:

::

    auto_process.py make_fastqs --protocol=icell8 ...

Subsequently the read pairs must be processed using the
``process_icell8.py`` utility described in the
:ref:`icell8_qc_and_filtering_protocol` section, to post-process
the Fastqs.

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

.. _10x_chromium_sc_fastq_generation:

Fastq generation for 10xGenomics Chromium single-cell data
----------------------------------------------------------

Fastq generation can be performed for 10xGenomics Chromium
single-cell data by using the ``--protocol=10x_chromium_sc``
option:

::

    auto_process.py make_fastqs --protocol=10x_chromium_sc ...

which fetches the data and runs ``cellranger mkfastq``.

.. note::

   This replaces the ``process_10xgenomics.py mkfastq`` command,
   which is now deprecated and likely to be removed in future.

This will generate the Fastqs in the specified output directory
(e.g. ``bcl2fastq``) along with an HTML report derived from the
``cellranger`` JSON QC summary file, statistics for the Fastqs.

Specifying the platform
-----------------------

``make_fastqs`` will try to identify the sequencer platform based on
the instrument name. If the instrument name is not recognised then
it can be explicitly specified using the ``--platform`` option.

Handling runs with mixed data: splitting lanes
----------------------------------------------

The HISeq platform enables runs with samples requiring different
processing protocols appearing within a single run. For example:

* Lanes 1 to 6 have samples with 8bp dual indexes, but lanes 7
  and 8 have 6p single index
* Lanes 1 and 2 have 10xGenomics Chromium or ICELL8 single-cell
  samples, but the remaining have contain standard samples

In cases such as these the recommended procedure is to prepare a
single sample sheet which contains appropriate indexes for each
lane, and split the processing by running ``make_fastqs`` multiple
times for each set of lanes using the ``--lanes`` option, and
specifying the appropriate options in each case.
