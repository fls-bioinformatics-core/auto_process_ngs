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

By default ``make_fastqs`` performs the following steps:

* Fetches the BCL data and copies it to the ``primary_data`` subdirectory
  of the analysis directory
* Checks the sample sheet for errors and possible barcode collisions
* Runs BCL to Fastq conversion (either ``bcl2fastq`` or ``bcl-convert``,
  or the appropriate pipeline for  for 10x Genomics single cell or spatial
  data e.g. ``cellranger``) and verifies that the outputs match what was
  expected from the input sample sheet
* Generates statistics for the Fastq data and produces a report on the
  analysing the numbers of reads assigned to each sample, lane and
  project (``processing_qc.html``)
* Analyses the barcode index sequences to identify possible demultiplexing
  issues

Various options are available to skip or control each of these stages;
more detail on the different usage modes can be found in the
subsequent sections:

* :ref:`make_fastqs-adapter-trimming-and-masking`
* :ref:`make_fastqs-mixed-protocols`
* :ref:`make_fastqs-bcl-converter`

Information on the different Fastq generation protocols can be found in
:ref:`make_fastqs-protocols`, and some of the other useful options can be
for the found in :ref:`make_fastqs-commonly-used-options`.

The outputs produced on successful completion are described below
in the section :ref:`make_fastqs-outputs`; it is recommended to check
the :doc:`processing QC <../output/processing_qc>` and
:doc:`barcode analysis <../output/barcode_analysis>` reports which
will highlight issues with the demultiplexing.

Once the Fastqs have been generated, the next step is to set up the
project directories - see
:doc:`Setting up project directories <setup_analysis_dirs>`.

.. _make_fastqs-protocols:

Fastq generation protocols
--------------------------

The ``--protocol`` option should be set according to the type of data
that are being processed:

======================== =====================================
Protocol option          Used for
======================== =====================================
``standard``             Standard Illumina sequencing data
                         (default)
``mirna``                miRNA-seq data
``10x_chromium_sc``      10xGenomics Chromium single cell
                         RNA-seq, CellPlex and Flex data
``10x_atac``             10xGenomics Chromium single cell
                         ATAC-seq data
``10x_visium``           10xGenomics Visium spatial RNA-seq
                         data
``10x_multiome``         10xGenomics Multiome single cell
                         ATAC or GEX data only in single
                         run
``10x_multiome_atac``    10xGenomics Multiome single cell
                         ATAC data (run with pooled GEX
                         and ATAC data)
``10x_multiome_gex``     10xGenomics Multiome single cell
                         GEX data only (run with pooled GEX
                         and ATAC data)
``icell8``               ICELL8 single-cell RNA-seq data
``icell8_atac``          ICELL8 single-cell ATAC-seq data
======================== =====================================

Typically the ``standard`` protocol is sufficient for most types of
data (RNA-seq, ATAC-seq, ChIP-seq, metagenomics etc).

.. note::

   For data where the sequences are expected to be very short (such
   as miRNA-seq data), the ``mirna`` protocol should be used instead -
   this is the same as the ``standard`` protocol but adjusts the
   adapter trimming and masking options as follows:

   * Sets the minimum trimmed read length to 10 bases
   * Turns off short read masking by setting the threshold length
     to zero

   More details about adapter trimming and short read masking can be
   found in the section :ref:`make_fastqs-adapter-trimming-and-masking`.

For other types of data (typically single cell and spatial), refer to
the appropriate section of the documentation for more details of which
Fastq generation protocols should be used:

* :doc:`10x Genomics single cell data <../single_cell/10x_single_cell>`
* :doc:`10x Genomics spatial data <../spatial/10x_visium>`
* :doc:`Parse Evercode single cell data <../single_cell/parse>`
* :doc:`ICELL single cell data <../single_cell/icell8>`

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
  :ref:`make_fastqs-mixed-protocols` for more details
* ``--bcl-converter``: allows the Illumina Fastq generation
  software to be specified, see :ref:`make_fastqs-bcl-converter`
  for more details
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
:ref:`'make_fastqs' <commands_make_fastqs>` section of the command
reference.

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
``tenx_filter_single_index=yes|no``   Set ``--filter-single-index``
                                      option for ``cellranger``
                                      or ``cellranger-arc``
``tenx_filter_dual_index=yes|no``     Set ``--filter-dual-index``
                                      option for ``cellranger``
                                      or ``cellranger-arc``
``spaceranger_rc_i2_override=BOOL``   Set ``--rc-i2-override`` option
                                      for ``spaceranger`` (can be
                                      either ``true`` or ``false``)
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

.. _make_fastqs-processing-same-run-multiple-times:

Processing a single run multiple times
--------------------------------------

Sometimes it is necessary to process a single run multiple times,
(for example, to try different parameter sets) while keeping the
outputs from each processing attempt in the same analysis
directory.

The ``--id`` option of the ``make_fastqs`` command can be used to
facilitate this, by allowing an identifier (e.g. ``no_trimming``)
to be supplied which will then be appended to the outputs from the
Fastq generation (including the output directories holding the
generated Fastqs, the barcode analysis directories, and the
statistics and processing report files).

For example:

::

   auto_process.py make_fastqs --id=no_trimming --no-adapter-trimming

would produce ``bcl2fastq_no_trimming``, ``barcodes_no_trimming``,
``statistics_no_trimming.info`` and so on.

.. note::

   The ``--id`` option of the ``setup_analysis_dirs`` command
   can be used to create projects which carry the same identifier,
   see :ref:`setup_analysis_dirs-add-identifier`.

.. note::

   A simpler alternative is to set up a completely new parallel
   analysis directory for reprocessing, and expliciting assigning
   a unique analysis number to distinguish it from other analysis
   attempts.

   This can be done via the ``setup`` command using the ``-n``
   option (see :ref:`setup_specifying_analysis_run_number`), or by
   setting the ``analysis_number`` metadata item within an existing
   analysis directory.

.. _make_fastqs-bcl-converter:

Specifying Illumina BCL conversion software
-------------------------------------------

For the ``standard`` and ``mirna`` Fastq generation protocols,
it is possible to use either the ``bcl2fastq`` or ``bcl-convert``
software packages to convert raw BCL data into Fastq files.

The ``--bcl-converter`` command line option can be used to
specify both the BCL converter software and optionally also
restrict to a range (or single version), for example:

::

   auto_process.py make_fastqs --bcl-converter 'bcl-convert>=3.7'

Default BCL conversion software can be specified in the config
file, both generally and on a per-platform basis (see
:ref:`specifying_bcl_conversion_software`).

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
* A :doc:`projects.info <../control_files/projects_info>` metadata
  file which is used by the :doc:`setup_analysis_dirs <setup_analysis_dirs>`
  command when setting up analysis project directories (see
  :doc:`Setting up project directories <setup_analysis_dirs>`)

For standard runs there will additional outputs:

* A directory called ``barcode_analysis`` which will contain
  reports with analysis of the barcode index sequences (see the
  section on :doc:`Barcode analysis <../output/barcode_analysis>`)

If the run included 10xGenomics Chromium 3' data then there will
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
