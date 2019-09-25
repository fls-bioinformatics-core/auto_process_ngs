Barcode Analysis Reports
========================

The demultiplexing stage of processing assigns reads to different
samples by matching the index sequence (aka "barcode index") of each
read against the reference index sequence for the samples, as defined
in the input ``SampleSheet.csv`` file. Reads which cannot be matched
against any sample are assigned to the "undetermined" category.

Typically the percentage of these undetermined reads should be low
compared to the total number of reads. However in cases where
demultiplexing hasn't worked well, the number of undetermined reads
may be unusually high (typically corresponding with an unexpectedly
low number of reads assigned to one or more samples). This can
happen artifically because of errors in the sample sheet, or because
of errors in the sequencing.

In these cases the barcode analysis report can be useful for diagnosing
and understanding various problems with the demultiplexing, by
identifying the following issues:

 * **Underrepresented samples**: samples which have extremely low
   or zero numbers of reads assigned to them;
 * **Overrepresented unassigned index sequences**: barcodes which
   have high numbers of associated reads but which aren't assigned
   to a sample in the sample sheet.

.. note::

   By default the analysis discounts *individual* barcode sequences
   with extremely low associated reads (less than 0.0001% of the
   total number of reads) before grouping similar sequences according
   to allowed mismatches.

   This is done to remove large numbers of spurious sequences which
   would otherwise cause the analysis to take much longer to complete.

   However it also means that the report read counts are no longer
   accurate, so they should be treated as approximations and are
   likely to underreport the number of reads compared with those
   reported in the per-sample statistics.

Barcode analysis is normally run automatically as part of the
:doc:`auto_process make_fastqs <../using/make_fastqs>` command (for
standard sequencing runs). It can also be run as a separate operation
using the
:ref:`auto_process.py analyse_barcodes <commands_analyse_barcodes>`
command.

The analysis counts the reads associated with each barcode sequence
in each lane across all the FASTQs generated for those lanes. It then
and reports the "top" (i.e. most frequent) barcode sequences across the
run, for example:

::

    Barcodes which cover less than 0.1% of reads have been excluded
    Reported barcodes cover 93.9% of the data
    
    #Rank	Index	Sample	N_seqs	N_reads	%reads	(%Total_reads)
        1	AAGAT	KL05_L1	1	17623152	9.5%	(9.5%)
        2	CAATG	KL08_L1	1	12784251	6.9%	(16.4%)
        3	TAGTA	KL06_L1	1	12361719	6.7%	(23.1%)
        4	GTATT	KL15_L1	1	10011071	5.4%	(28.5%)
    ...

If an index sequence corresponds to a sample defined in the sample
sheet, then the name of the sample is also shown.

When demultiplexing has been successful, the top ranked barcodes should
all correspond to sample names. Any gaps or samples missing from the top
barcodes indicate that there may be a problems.

For example:

::

    #Rank	Index	Sample	N_seqs	N_reads	%reads	(%Total_reads)
        1	GCTACGCT+CTAAGCCT	SM9	1	2894178	13.3%	(13.3%)
        2	CAGAGAGG+AGAGTAGA	SM8	1	2454006	11.3%	(24.6%)
        3	TCCTGAGC+CTCTCTAT	SM4	1	2241825	10.3%	(34.9%)
        4	AGGCAGAA+CTCTCTAT	SM3	1	2206475	10.1%	(45.0%)
        5	CGAGGCTG+CTAAGCCT		1	2141506	9.8%	(54.9%)
        6	CTCTCTAC+AGAGTAGA	SM7	1	2094131	9.6%	(64.5%)
        7	TAGGCATG+TATCCTCT	SM6	1	1943974	8.9%	(73.4%)
        8	GGACTCCT+TATCCTCT	SM5	1	1868127	8.6%	(82.0%)
        9	CGTACTAG+TAGATCGC	SM2	1	1659786	7.6%	(89.7%)
       10	TAAGGCGA+TAGATCGC		1	942645	4.3%	(94.0%)
       11	CGGCAGAA+CTCTCTAT		1	61005	0.3%	(94.3%)
       12	CAGAGAGG+CGAGTAGA		1	38768	0.2%	(94.5%)
       13	CTCTCTAC+CGAGTAGA		1	26927	0.1%	(94.6%)

Here the 5th ranked barcode sequence doesn't have an associated
sample name. This overrepresented sequence is highlighted after
the list:

::
    
    The following unassigned barcodes are overrepresented compared to the assigned barcodes:
    #Index	N_reads	%reads
    CGAGGCTG+CTAAGCCT	2141506	9.84%

In addition the analysis reports any underrepresented and
missing samples, for example:

::

    The following samples are underrepresented:
    
    	#Sample	Index	N_reads	%reads
    	SM1	CAGAGAGG+AGAGTAGT	9971	0.05%
    	SM10	CGAGGCTG+AAGCCTCT		<0.05%

.. note::

   Note that if weeding of the initial data has been performed (which
   is the default) then the counts are only approximate and so it's
   not possible to say if underrepresented samples with zero counts are
   really missing. Therefore they are reported as having read fractions
   below the cutoff threshold, as in the example above.

The analysis is performed in the ``barcode_analysis`` subdirectory
by default, where the report is written as text, XLS and HTML files.

The raw counts are also cached in files in the
``barcode_analysis/counts/`` subdirectory. This means that the analysis
can be rerun quickly with different parameters (for example a different
sample sheet or cut-off) if desired.

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
