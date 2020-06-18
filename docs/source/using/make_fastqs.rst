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
``10x_chromium_sc``      10xGenomics Chromium single-cell
                         RNA-seq data
``10x_atac``             10xGenomics Chromium single-cell
                         ATAC-seq data
======================== =====================================

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

* :ref:`make_fastqs-standard-protocol`
* :ref:`make_fastqs-mirna-protocol`
* :ref:`make_fastqs-icell8-protocol`
* :ref:`make_fastqs-icell8-atac-protocol`
* :ref:`make_fastqs-10x_chromium_sc-protocol`
* :ref:`make_fastqs-10x_atac-protocol`
* :ref:`make_fastqs-adapter-trimming-and-masking`
* :ref:`make_fastqs-mixed-protocols`

Information on other commonly used options can be found
:ref:`below <make_fastqs-commonly-used-options>`.

Advice on handling other unusual cases and problems can be found
here:

* :doc:`Troubleshooting unusual situations <troubleshooting>`

Once the Fastqs have been generated, the next step is to set up the
project directories - see
:doc:`Setting up project directories <setup_analysis_dirs>`.

.. _make_fastqs-commonly-used-options:

Commonly used options
---------------------

Some of the most commonly used options are:

* ``--output-dir``: specifies the directory to write the output
  Fastqs to (defaults to ``bcl2fastq``)
* ``--sample-sheet``: specifies a non-default sample sheet file
  to use (defaults to ``custom_SampleSheet.csv``; the new sample
  sheet file will become the default for subsequent runs)
* ``--lanes``: allows a subset of lanes to be processed (useful
  for multi-lane sequencers when samples with a mixture
  of processing protocols have been run). Lanes can be specified
  as a range (e.g. ``1-4``), a list (e.g. ``6,8``) or a
  combination (e.g. ``1-4,6,8``)
* ``--use-bases-mask``: allows a custom bases mask string (which
  controls how each cycle of raw data is used) to be specified
  (default is to determine the bases mask automatically; set to
  ``auto`` to restore this behaviour)
* ``--platform``: if the sequencer platform cannot be identified
  from the instrument name it can be explicitly specified using
  this option (see :ref:`config_sequencer_platforms` for how to
  associate sequencers and platforms in the configuration)
* ``--no-barcode-analysis`` skips the barcode analysis for
  standard runs (helpful when handling runs requiring multiple
  rounds of processing; see :ref:`make_fastqs-mixed-protocols`)
* ``--no-stats`` skips the generation of statistics and processing
  QC reporting (helpful when handling runs requiring multiple
  rounds of processing; see :ref:`make_fastqs-mixed-protocols`)

The full set of options can be found in the
:ref:`'make_fastqs' section of the command reference <commands_make_fastqs>`.

.. _make_fastqs-standard-protocol:

Standard Fastq generation (``--protocol=standard``)
---------------------------------------------------

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

miRNA-seq Fastq generation (``--protocol=mirna``)
-------------------------------------------------

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

Fastq generation for ICELL8 single-cell RNA-seq data (``--protocol=icell8``)
----------------------------------------------------------------------------

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

Fastq generation for ICELL8 single-cell ATAC-seq data (``--protocol=icell8_atac``)
----------------------------------------------------------------------------------

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

Fastq generation for 10xGenomics Chromium single-cell RNA-seq data (``--protocol=10x_chromium_sc``)
---------------------------------------------------------------------------------------------------

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

Fastq generation for 10xGenomics single-cell ATAC-seq data (``--protocol=10x_atac``)
---------------------------------------------------------------------------------------------------------

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
   than this length after trimming, by padding them with ``N``s
   up to the minimum length
 * **Masking of short reads** is performed for reads below a
   masking threshold length, by masking *all* bases in the read
   with ``N``s

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

Fastq generation for runs with mixed protocols
----------------------------------------------

Multi-lane instruments such as the HISeq platform provide the
option to run mixtures of samples requiring different processing
protocols in a single sequencing run, for example:

* Samples in some lanes have different barcode index
  characteristics (e.g. different lengths) to those in
  other lanes
* Some lanes contain standard samples whilst others contain
  10xGenomics or ICELL8 single-cell samples

In these cases the data cannot be processed in a single
``make_fastqs`` run. Instead the recommended procedure for
handling these situations is:

1. Prepare a single sample sheet with the appropriate indexes
   for each lane (for example truncating index sequences, or
   inserting the appropriate 10xGenomics indexes)
2. Run ``make_fastqs`` multiple times to process each subset of
   lanes on their own using the ``--lanes`` option, specifying the
   appropriate protocol and processing options and writing the
   Fastqs for each to a different output directory using the
   ``--output-dir`` option
3. Combine the outputs from each subset into a single output
   directory using the ``merge_fastq_dirs`` command
4. (Re)generate the statistics and QC report on the merged
   data using the ``update_fastq_stats`` command

For example: say we have a HISeq run with non-standard samples
in lanes 5 and 6, and standard samples in all other lanes. In
this case, after updating the samplesheet the standard samples
would be processed first:

::

   auto_process.py make_fastqs \
            --lanes=1-4,7-8 \
	    --sample-sheet=SampleSheet.updated.csv \
            --output-dir=bcl2fastq.L123478 \
            --no-barcode-analysis \
	    --no-stats

The ``--lanes`` option restricts the lanes to just those with
the standard samples. ``--output-dir`` writes the Fastqs to a
custom output directory. Specifying ``--no-stats`` suppresses
the statistics generation at this stage.

Next process the non-standard samples, for example: if the
samples in lanes 5 and 6 had different barcode lengths:

::

   auto_process.py make_fastqs \
            --lanes=5-6 \
            --output-dir=bcl2fastq.L56 \
            --use-bases-mask=auto \
            --no-barcode-analysis \
	    --no-stats

Alternatively if the data in these lanes were 10xGenomics
Chromium single cell data:

::

   auto_process.py make_fastqs \
            --lanes=5-6 \
	    --protocol=10x_chromium_sc \
            --output-dir=bcl2fastq.L56 \
            --use-bases-mask=auto \
	    --no-stats

The outputs from each subset of lanes can be merged into a
single output directory using the ``merge_fastq_dirs`` command.
For example:

::

   auto_process.py merge_fastq_dirs \
             --primary-unaligned-dir=bcl2fastq.L123478 \
	     --output-dir=bcl2fastq

To generate the statistics and processing QC report for the
merged data use the ``update_fastq_stats`` command:

::

   auto_process.py update_fastq_stats

To perform the barcode analysis for the merged data use the
``analyse_barcodes`` command:

::

   auto_process.py analyse_barcodes

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
* A ``projects.info`` metadata file which is used for setting up
  analysis project directories (see
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
