Analyse index sequences in demultiplexed Fastqs using ``auto_process analyse_barcodes``
=======================================================================================

Overview
--------

The demultiplexing stage of processing assigns reads to different samples
by matching the index sequence (aka "barcode index") of each read against
the reference index sequence for the samples, as defined in the input
``SampleSheet.csv`` file.

Reads which cannot be matched against any sample are assigned to the
"undetermined" category. Typically the percentage of these undetermined reads
should be low compared to the total number of reads.

However in cases where demultiplexing hasn't worked well, the number of
undetermined reads maybe unusually high (typically corresponding with an
unexpectedly low number of reads assigned to one or more samples). This can
happen artifically because of errors in the sample sheet, or because of
errors in the sequencing.

The ``analyse_barcodes`` command (and underlying ``analyse_barcodes.py``
utility) can be used to help diagnose and understand problems with barcode
assignment in these situations.

Performing the analysis
-----------------------

Barcode analysis can be run after the :doc:`make_fastqs <make_fastqs>`
stage. Its simplest invocation is:

::

    auto_process.py analyse_barcodes

This counts the reads associated with each barcode sequence for each
FASTQ file output by ``bcl2fastq``, and reports the "top" (i.e. most
frequent) barcode sequences across the run. For example:

::

    Barcodes which cover less than 0.1% of reads have been excluded
    Reported barcodes cover 93.9% of the data
    
    #Rank	Index	Sample	N_seqs	N_reads	%reads	(%Total_reads)
        1	AAGAT	KL05_L1	1	17623152	9.5%	(9.5%)
        2	CAATG	KL08_L1	1	12784251	6.9%	(16.4%)
        3	TAGTA	KL06_L1	1	12361719	6.7%	(23.1%)
        4	GTATT	KL15_L1	1	10011071	5.4%	(28.5%)
    ...

Where an index sequence corresponds to a sample defined in the sample
sheet, the name of the sample is shown. We expect that the top ranked
barcodes should all correspond to samples, so "gaps" in the list
indicate that there may be a problem, e.g.:

::

    ...
       16	CAACT		1	4170435	2.3%	(76.1%)
       17	CTATT		1	3930029	2.1%	(78.2%)
       18	ACACG	KL13_L1	1	3692276	2.0%	(80.2%)
       19	AGTTC		1	3514693	1.9%	(82.1%)
    ...

The analysis is performed in the ``barcode_analysis`` subdirectory.
The counts are cached in files in the ``barcode_analysis/counts/``
directory, allowing the the command can be run more efficiency on
subsequent occasions using different parameters if necessary.

Tuning the reporting
--------------------

The ``analyse_barcodes`` command supports a number of options to allow
"tuning" of the analysis and reporting:

 * ``--cutoff``: by default barcodes are excluded from the analysis if
   they are associated with fewer than 0.1% of all the reads counted.
   Use this option to adjust the read limit to report more or fewer
   barcodes.

 * ``--mismatches``: by default only exact matches to barcode sequences
   are considered (i.e. zero mismatches). If the number of mismatches
   (as specified by this option) are more than zero then similar barcodes
   will be grouped together.

 * ``--lanes``: by default all barcodes and counts are combined and
   analysed together. Use this option to restrict analysis to a subset
   of lanes.

 * ``--sample-sheet``: by default the sample names and corresponding
   reference barcodes are taked from the sample sheet used in the
   ``make_fastqs`` stage; this option can be used to specify a different
   sample sheet file to use.

.. note::

   When ``--cutoff`` is used with ``--mismatches`` then the read cutoff
   is applied **after** the grouping.

Running analyse_barcodes.py directly
------------------------------------

The ``analyse_barcodes.py`` utility can be run directly, either using
the cached counts produced by the ``analyse_barcodes`` command, or from
scratch, to perform more nuanced analyses of the barcode sequences if
required.

``analyse_barcodes.py`` operates in two stages:

 1. Counting the frequency of each index sequence in the target FASTQs

 2. Analysis of these raw counts to determine groups of similar barcode
    sequences (optionally), and report the most numerous barcodes (or
    groups) matched against reference sequences from a sample sheet.

The first stage can be very time consuming so the counts can be output to an
intermediate ``counts`` file using the ``-o`` option::

    analyse_barcodes.py ... -o SAMPLE.counts SAMPLE.fq

.. note::

   To suppress analysis and reporting when generating counts
   use the ``--no-report`` option.

Multiple analyses can be performed using the cached counts, which are
reloaded into the program using the ``-c`` option::

    analyse_barcodes.py ... -c SAMPLE.counts

Multiple counts files can be combined via the ``-c`` option::

    analyse_barcodes.py ... -c SAMPLE_1.counts SAMPLE_2.counts ...

.. note::

   The ``analyse_barcodes`` command generates counts files for each
   FASTQ file, in the ``barcode_analysis/counts/`` directory, using
   the naming convention of ``FASTQ.counts``.

By default the results of the analysis are written to stdout; use
the ``-r`` option to specify an output file instead.

Analysing undetermined barcodes only
------------------------------------

Currently this can be done by running the ``analyse_barcodes.py`` utility
directly on the cached counts for just the "undetermined" FASTQ files,
for example::

    analyse_barcodes.py -c barcode_analysis/Undetermined*.counts
