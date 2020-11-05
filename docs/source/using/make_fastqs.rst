Fastq generation using ``auto_process make_fastqs``
===================================================

Overview
--------

The ``make_fastqs`` command is the backbone of the ``auto_process``
pipeline. It is run after creating the analysis directory using the
:doc:`setup <setup>` command and performs the key step of generating
Fastq files from the raw BCL data produced by the sequencer.

The general invocation of the command is:

::

   auto_process.py make_fastqs [--protocol=PROTOCOL] *options* [ANALYSIS_DIR]

The ``--protocol`` option should be chosen according to the type
of data that is being processed:

======================== =====================================
Protocol option          Used for
======================== =====================================
``standard``             Standard Illumina sequencing data
                         (default)
``mirna``                miRNA-seq data
``icell8``               ICELL8 single-cell RNA-seq data
``icell8_atac``          ICELL8 single-cell ATAC-seq data
``10x_chromium_sc``      10xGenomics Chromium single cell
                         RNA-seq data
``10x_atac``             10xGenomics Chromium single cell
                         ATAC-seq data
``10x_visium``           10xGenomics Visium spatial RNA-seq
                         data
``10x_multiome``         10xGenomics Multiome single cell
                         ATAC + GEX data
======================== =====================================

The protocols are described in the section :ref:`make_fastqs-protocols`;
information on most useful options can be found in
:ref:`make_fastqs-commonly-used-options`.

By default ``make_fastqs`` performs the following steps:

* Fetches the BCL data and copies it to the ``primary_data`` subdirectory
  of the analysis directory
* Checks the sample sheet for errors and possible barcode collisions
* Runs ``bcl2fastq`` (or ``cellranger`` for 10xGenomics Chromium data)
  and verifies that the outputs match what was expected from the input
  sample sheet
* Generates statistics for the Fastq data and produces a report on the
  analysing the numbers of reads assigned to each sample, lane and
  project (``processing.html``)
* Analyses the barcode index sequences to look for issues with
  demultiplexing (``standard`` protocol only)

Various options are available to skip or control each of these stages;
more detail on the different usage modes can be found in the
subsequent sections:

* :ref:`make_fastqs-adapter-trimming-and-masking`
* :ref:`make_fastqs-mixed-protocols`

The outputs produced on successful completion are described below
in the section :ref:`make_fastqs-outputs`.

.. note::

   Advice on handling unusual cases and problems can be found
   in the section :doc:`troubleshooting`.

Once the Fastqs have been generated, the next step is to set up the
project directories - see
:doc:`Setting up project directories <setup_analysis_dirs>`.

.. _make_fastqs-protocols:

Fastq generation protocols
--------------------------

The pre-defined protocols for Fastq generation are described in
the sections below:

* :ref:`make_fastqs-standard-protocol`
* :ref:`make_fastqs-mirna-protocol`
* :ref:`make_fastqs-icell8-protocol`
* :ref:`make_fastqs-icell8-atac-protocol`
* :ref:`make_fastqs-10x_chromium_sc-protocol`
* :ref:`make_fastqs-10x_atac-protocol`
* :ref:`make_fastqs-10x_visium-protocol`
* :ref:`make_fastqs-10x_multiome-protocol`

.. _make_fastqs-standard-protocol:

Standard data (``--protocol=standard``)
***************************************

The Fastq generation for standard data is performed using a command
of the form:

::

   auto_process.py make_fastqs ...

The outputs produced on successful completion are described below
in the section :ref:`make_fastqs-outputs`; it is recommended to check
the :doc:`processing QC <../output/processing_qc>` and
:doc:`barcode analysis <../output/barcode_analysis>` reports which
will highlight issues with the demultiplexing.

.. _make_fastqs-mirna-protocol:

miRNA-seq data (``--protocol=mirna``)
*************************************

Initial Fastqs can be generated from miRNA-seq data using the
``--protocol=mirna`` option:

::

    auto_process.py make_fastqs --protocol=mirna ...

This adjusts the adapter trimming and masking options as follows:

 * Sets the minimum trimmed read length to 10 bases
 * Turn off short read masking by setting the threshold length
   to zero

Subsequently the Fastq generation is the same as the standard
protocol described in :ref:`make_fastqs-standard-protocol`.

More details about adapter trimming and short read masking can be
found in the section :ref:`make_fastqs-adapter-trimming-and-masking`.

.. _make_fastqs-icell8-protocol:

ICELL8 single-cell RNA-seq data (``--protocol=icell8``)
*******************************************************

Initial Fastqs can be generated from ICELL8 single-cell8 RNA-seq data
using the ``--protocol=icell8`` option:

::

    auto_process.py make_fastqs --protocol=icell8 ...

Subsequently the read pairs must be processed using the
``process_icell8.py`` utility described in the
:ref:`icell8_scRNA-seq_qc_and_filtering_protocol` section, to post-process
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

.. _make_fastqs-icell8-atac-protocol:

ICELL8 single-cell ATAC-seq data (``--protocol=icell8_atac``)
*************************************************************

Initial Fastqs can be generated from ICELL8 single-cell8 ATAC-seq data
using the ``--protocol=icell8_atac`` option:

::

    auto_process.py make_fastqs --protocol=icell8_atac --well-list=WELL_LIST_FILE...

This runs bcl2fastq to perform initial standard demultiplexing based on
the samples defined in the sample sheet, followed by a second round of
demultiplexing into ICELL8 samples based on the contents of the well list
file which must be supplied via the mandatory ``--well-list`` argument.

.. warning::

   This protocol is still under development.

.. _make_fastqs-10x_chromium_sc-protocol:

10xGenomics Chromium single-cell RNA-seq data (``--protocol=10x_chromium_sc``)
******************************************************************************

Fastq generation can be performed for 10xGenomics Chromium
single-cell RNA-seq data by using the ``--protocol=10x_chromium_sc``
option:

::

    auto_process.py make_fastqs --protocol=10x_chromium_sc ...

which fetches the data and runs ``cellranger mkfastq``.

This will generate the Fastqs in the specified output directory
(e.g. ``bcl2fastq``) along with an HTML report derived from the
``cellranger`` JSON QC summary file, statistics for the Fastqs.

.. note::

   ``make_fastqs`` offers various options for controlling the
   behaviour of ``cellranger mkfastqs``, for example setting the
   jobmode (see :ref:`10xgenomics-additional-options`).

.. _make_fastqs-10x_atac-protocol:

10xGenomics single-cell ATAC-seq data (``--protocol=10x_atac``)
***************************************************************

Fastq generation can be performed for 10xGenomics single-cell
ATAC-seq data by using the ``--protocol=10x_atac`` option:

::

    auto_process.py make_fastqs --protocol=10x_atac ...

which fetches the data and runs ``cellranger-atac mkfastq``.

This will generate the Fastqs in the specified output directory
(e.g. ``bcl2fastq``) along with an HTML report derived from the
``cellranger-atac`` JSON QC summary file, statistics for the Fastqs.

.. note::

   ``make_fastqs`` offers various options for controlling the
   behaviour of ``cellranger-atac mkfastqs``, for example setting the
   jobmode (see :ref:`10xgenomics-additional-options`).

.. _make_fastqs-10x_multiome-protocol:

10xGenomics single cell Multiome ATAC + GEX data (``--protocol=10x_multiome``)
******************************************************************************

Fastq generation can be performed for 10xGenomics single cell
multiome ATAC and gene expression (GEX) data by using the
``--protocol=10x_multiome`` option:

::

    auto_process.py make_fastqs --protocol=10x_multiome ...

which fetches the data and runs ``cellranger-atac mkfastq``.

This will generate the Fastqs in the specified output directory
(e.g. ``bcl2fastq``) along with an HTML report derived from the
``cellranger-arc`` JSON QC summary file, statistics for the Fastqs.

By default adapter trimming is automatically disabled by removing
the adapter sequences in the sample sheet.

.. note::

   ``make_fastqs`` offers various options for controlling the
   behaviour of ``cellranger-arc mkfastqs``, for example setting the
   jobmode (see :ref:`10xgenomics-additional-options`).

.. _make_fastqs-10x_visium-protocol:

10xGenomics spatial RNA-seq data (``--protocol=10x_visium``)
************************************************************

Fastq generation can be performed for 10xGenomics spatial RNA-seq
ata by using the ``--protocol=10x_visium`` option:

::

    auto_process.py make_fastqs --protocol=10x_visium ...

which fetches the data and runs ``spaceranger mkfastq``.

This will generate the Fastqs in the specified output directory
(e.g. ``bcl2fastq``) along with an HTML report derived from the
``spaceranger`` JSON QC summary file, statistics for the Fastqs.

.. note::

   ``make_fastqs`` offers various options for controlling the
   behaviour of ``spaceranger mkfastqs``, for example setting the
   jobmode (see :ref:`10xgenomics-additional-options`).

.. _make_fastqs-commonly-used-options:

Commonly used options
---------------------

Some of the most commonly used options are:

* ``--protocol``: specifies the Fastq generation protocol
* ``--output-dir``: specifies the directory to write the output
  Fastqs to (defaults to ``bcl2fastq``)
* ``--sample-sheet``: specifies a non-default sample sheet file
  to use (defaults to ``custom_SampleSheet.csv``; the new sample
  sheet file will become the default for subsequent runs)
* ``--lanes``: allows a subset of lanes to be processed (useful
  for multi-lane sequencers when samples with a mixture
  of processing protocols have been run). Lanes can be specified
  as a range (e.g. ``1-4``), a list (e.g. ``6,8``) or a
  combination (e.g. ``1-4,6,8``). See
  :ref:`make_fastqs-mixed-protocols` for more details.
* ``--use-bases-mask``: allows a custom bases mask string (which
  controls how each cycle of raw data is used) to be specified
  (default is to determine the bases mask automatically; set to
  ``auto`` to restore this behaviour)
* ``--platform``: if the sequencer platform cannot be identified
  from the instrument name it can be explicitly specified using
  this option (see :ref:`config_sequencer_platforms` for how to
  associate sequencers and platforms in the configuration)
* ``--no-barcode-analysis`` skips the barcode analysis for
  standard runs
* ``--no-stats`` skips the generation of statistics and processing
  QC reporting

The full set of options can be found in the
:ref:`'make_fastqs' section of the command reference <commands_make_fastqs>`.

.. _make_fastqs-adapter-trimming-and-masking:

Configuring adapter trimming and masking
----------------------------------------

By default Fastq generation includes adapter trimming and masking of
short reads via ``bcl2fastq``.

Adapter sequences used for trimming are taken from those specified
in the input sample sheet, but these can be overriden by using the
``--adapter`` and ``--adapter-read2`` options to specify different
sequences.

Adapter trimming can be disabled by specifying the
``--no-adapter-trimming`` option (or by setting both adapter
sequences to empty strings).

When adapter trimming is performed two additional operations are
applied:

* **Minium read length** is enforced for reads which are shorter
  than this length after trimming, by padding them with N's
  up to the minimum length
* **Masking of short reads** is performed for reads below a
  masking threshold length, by masking *all* bases in the read
  with N's

Minimum read length defaults to 35 bases but can set explicitly by
using the ``--minimum-trimmed-read-length`` option; the masking
threshold defaults to 22 bases but can be set using the
``--mask-short-adapter-reads`` option. Set this to zero to turn
off masking.

.. warning::

   Setting the minimum read length to zero when using adapter
   trimming can result in read records with zero-length sequences,
   which may cause problems in downstream analyses.

.. _make_fastqs-mixed-protocols:

Fastq generation for runs with mixed protocols and options
----------------------------------------------------------

Multi-lane instruments such as the HiSeq platform provide the
option to run mixtures of samples requiring different processing
protocols in a single sequencing run, for example:

* Samples in some lanes have different barcode index
  characteristics (e.g. different lengths) to those in
  other lanes
* Some lanes contain standard samples whilst others contain
  10xGenomics or ICELL8 single-cell samples

``make_fastqs`` is able to process these in a single run provided
that:

* the sample sheet has the appropriate index sequences for
  each lane (for example, truncating index sequences, or
  inserting the appropriate 10xGenomics indexes); and
* where different protocols or processing options need to
  be specified for groups of lanes, that these are specified
  via multiple ``--lanes`` options.

``make_fastqs`` will process each set of lanes separately
before combining them into a single output directory at the
end.

For example: say we have a HiSeq run with non-standard samples
in lanes 5 and 6, and standard samples in all other lanes.

If the samples in lanes 5 and 6 have different barcode lengths
to those in the other lanes, but should otherwise be treated
the same, then the following command line would be sufficient
to handle this:

::

   auto_process.py make_fastqs \
	    --sample-sheet=SampleSheet.updated.csv

However if the samples in lanes 5 and 6 were 10xGenomics
Chromium single cell data, then it is necessary to explicitly
specify which lanes to group together and how each group should
be handled. This is done using the ``--lanes`` option to
indicate that the ``10x_chromium_sc`` protocol should be used
with lanes 5 and 6, and that the ``standard`` protocol should
be used with the other lanes:

::

   auto_process.py make_fastqs \
            --lanes=1-4,7-8:standard \
	    --lanes=5,6:10x_chromium_sc \
	    --sample-sheet=SampleSheet.updated.csv


.. note::

   If the ``--lanes`` option is used one or more times then
   only those lanes explicitly listed will be processed.
   Lanes that aren't specified will be excluded from the
   processing.

More generally it's possible to set multiple options on a
set of lanes using the lanes option, for example to explicitly
specify the adapter sequences for lane 8:

::

   auto_process.py make_fastqs \
            --lanes=1-7 \
	    --lanes=8:adapter=CTGTCTCTTATACACATCT \
	    --sample-sheet=SampleSheet.updated.csv

The general form of the ``--lanes`` option is:

::

   --lanes=LANES[:protocol][:OPTION=VALUE[:OPTION=VALUE...]]

The available options are:

===================================== ==================================
Option                                Description
===================================== ==================================
``bases_mask=BASES_MASK``             Set bases mask
``trim_adapters=yes|no``              Turn adapter trimming on or off
``adapter=SEQUENCE``                  Set adapter sequence for trimming
``adapter_read2=SEQUENCE``            Set read2 adapter sequence
``minimum_trimmed_read_length=N``     Set minimum trimmed read length
``mask_short_adapter_reads=N``        Set minimum read length below
                                      which sequences are masked
``icell8_well_list=FILE``             Well list file (``icell8`` and
                                      ``icell8_atac`` protocols only)
``icell8_atac_swap_i1_and_i2=yes|no`` Turn I1/I2 swapping on or off
                                      (``icell8_atac`` protocol only)
``icell8_atac_reverse_complement``    Set reverse complementing option
                                      (``icell8_atac`` protocol only)
``analyse_barcodes=yes|no``           Turn barcode analysis on or off
===================================== ==================================

These options will override the defaults and any global values
set by the top-level options.

It is also possible to process subsets of lanes manually, and
then use the ``merge_fastq_dirs``, ``update_fastq_stats`` and
``analyse_barcodes`` commands to combine and analyse the Fastqs.

For example, for the mixture of standard and 10xGenomics samples
previously described this might look like:

::

   # Process lanes 1-4,7-8 (standard samples)
   auto_process.py make_fastqs \
            --lanes=1-4,7-8 \
	    --sample-sheet=SampleSheet.updated.csv \
            --output-dir=bcl2fastq.L123478 \
            --use-bases-mask=auto \
            --no-barcode-analysis \
	    --no-stats

   # Process lanes 5-6 (10xGenomics samples)
   auto_process.py make_fastqs \
            --lanes=5-6 \
	    --sample-sheet=SampleSheet.updated.csv \
	    --protocol=10x_chromium_sc \
            --output-dir=bcl2fastq.L56 \
            --use-bases-mask=auto \
	    --no-stats

   # Combine outputs
   auto_process.py merge_fastq_dirs \
             --primary-unaligned-dir=bcl2fastq.L123478 \
	     --output-dir=bcl2fastq

   # Generate statistics
   auto_process.py update_fastq_stats

   # Analyse barcodes (standard samples only)
   auto_process.py analyse_barcodes --lanes=1-4,7-8

See the appropriate sections of the command reference for
the full set of available options:

* :ref:`commands_merge_fastq_dirs`
* :ref:`commands_update_fastq_stats`
* :ref:`commands_analyse_barcodes`

.. _make_fastqs-outputs:

Outputs
-------

On completion the ``make_fastqs`` command will produce:

* An output directory called ``bcl2fastq`` with the demultiplexed
  Fastq files (see below for more detail)
* A set of tab-delimited files with statistics on each of the
  Fastq files
* An HTML report on the processing QC (see the section on
  :doc:`Processing QC reports <../output/processing_qc>`)
* A :doc:`projects.info <projects_info>` metadata file which
  is used by the :doc:`setup_analysis_dirs <setup_analysis_dirs>`
  command when setting up analysis project directories (see
  :doc:`Setting up project directories <setup_analysis_dirs>`)

For standard runs there will additional outputs:

* A directory called ``barcode_analysis`` which will contain
  reports with analysis of the barcode index sequences (see the
  section on :doc:`Barcode analysis <../output/barcode_analysis>`)

If the run included 10xGenomics Chromium 3'v2 data then there will
be some additional outputs:

* A report in the top-level analysis directory called
  ``cellranger_qc_summary[_LANES].html``, which is an HTML copy
  of the QC summary JSON file produced by ``cellranger mkfastq``
  (nb ``LANES`` will be the subset of lanes from the
  run which contained the Chromium data, if the run consisted
  of a mixture of Chromium and non-Chromium samples, for example:
  ``--lanes=5,6`` results in ``56``).

.. note::

   The processing QC reports can be copied to the QC server using
   the :doc:`publish_qc command <publish_qc>`.

Output Fastq files
******************

Each sample defined in the input sample sheet will produce one
or more output Fastq files, depending on:

* if the run was single- or paired-end,
* whether the sample appeared in more than one lane, and
* whether the ``--no-lane-splitting`` option was specified

By default if samples appear in more than one lane in a sequencing
run then ``make_fastqs`` will generate multiple Fastqs with
each Fastq only containing reads from a single lane, and with
the lane number appearing in the Fastq file name.

However if the ``--no-lane-splitting`` option is specified then
the reads from all lanes that the sample appeared in will be
combined into the same Fastq file.

The default lane splitting behaviour can be controlled via the
configuration options in the ``auto_process.ini`` file (see
:doc:`configuration <../configuration>`).

.. note::

   Lane splitting is always performed for 10xGenomics single cell
   data, regardless of the settings or options supplied to
   ``make_fastqs``.
