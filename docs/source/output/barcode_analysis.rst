Barcode analysis reports
========================

Index sequence (aka barcode) analysis is run as part of the
:doc:`auto_process make_fastqs <../using/make_fastqs>` command (for
standard sequencing runs); it can also be run as a separate operation
using the
:doc:`auto_process.py analyse_barcodes <commands_analyse_barcodes>`
command. It can be useful for diagnosing and understanding various
problems with the demultiplexing, for example if there are errors in
the sample sheet resulting in low or zero reads being assigned to
one or more samples.

Barcode anlysis can identify the following issues:

 * **Underrepresented and missing samples**: samples which have
   extremely low (underrepresented) or zero (missing) numbers of
   reads assigned to them;
 * **Overrepresented unassigned index sequences**: barcodes which
   have high numbers of associated reads but which aren't assigned
   to a sample in the sample sheet.

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
        1	GCTACGCTCTAAGCCT	SM9	1	2894178	13.3%	(13.3%)
        2	CAGAGAGGAGAGTAGA	SM8	1	2454006	11.3%	(24.6%)
        3	TCCTGAGCCTCTCTAT	SM4	1	2241825	10.3%	(34.9%)
        4	AGGCAGAACTCTCTAT	SM3	1	2206475	10.1%	(45.0%)
        5	CGAGGCTGCTAAGCCT		1	2141506	9.8%	(54.9%)
        6	CTCTCTACAGAGTAGA	SM7	1	2094131	9.6%	(64.5%)
        7	TAGGCATGTATCCTCT	SM6	1	1943974	8.9%	(73.4%)
        8	GGACTCCTTATCCTCT	SM5	1	1868127	8.6%	(82.0%)
        9	CGTACTAGTAGATCGC	SM2	1	1659786	7.6%	(89.7%)
       10	TAAGGCGATAGATCGC		1	942645	4.3%	(94.0%)
       11	CGGCAGAACTCTCTAT		1	61005	0.3%	(94.3%)
       12	CAGAGAGGCGAGTAGA		1	38768	0.2%	(94.5%)
       13	CTCTCTACCGAGTAGA		1	26927	0.1%	(94.6%)

Here the 5th ranked barcode sequence doesn't have an associated
sample name. This overrepresented sequence is highlighted after
the list:

::
    
    The following unassigned barcodes are overrepresented compared to the assigned barcodes:
    #Index	N_reads	%reads
    CGAGGCTGCTAAGCCT	2141506	9.84%

In addition the analysis reports any underrepresented and
missing samples, for example:

::

    The following samples are underrepresented:
    
    	#Sample	Index	N_reads	%reads
    	SM1	CAGAGAGGAGAGTAGT	9971	0.05%
    
    The following samples had no counts:
    
    	#Sample	Index
    	SM10	CGAGGCTGAAGCCTCT

The analysis is performed in the ``barcode_analysis`` subdirectory
by default, where the report is written as text, XLS and HTML files.

The raw counts are also cached in files in the
``barcode_analysis/counts/`` subdirectory. This means that the analysis
can be rerun quickly with different parameters (for example a different
sample sheet or cut-off) if desired.
