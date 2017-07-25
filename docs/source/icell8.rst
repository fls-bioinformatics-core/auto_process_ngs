Processing Wafergen ICell8 single-cell data
===========================================

Background
----------

The Wafergen ICell8 system prepares single cell (SC) samples which
are then sequenced as part of an Illumina sequencing run.

Initial processing of the sequencing data produces a single FASTQ file
R1/R2 pair, where the R1 reads hold an inline barcode and a UMI for
the read pair:

 * **Inline barcode:** bases 1-11 in R1
 * **UMI:** bases 12-21 in R1

FASTQs can be generated using the ``auto_process.py make_fastqs``
command, with the following settings suggested by Wafergen specifically
for NextSeq data:

 * Disable adapter trimming and masking by using
   ``--minimum-trimmed-read-length=21 --mask-short-adapter-reads=0``

This is recommended to stop unintentional trimming of UMI sequences
(which are mostly random) from the R1, should they happen to match
part of an adapter sequence.

Subsequently the read pairs can be processed using the utility script
``process_icell8.py`` to perform initial filtering and QC as described
below. The utility also splits the read pairs into FASTQs by both
barcode and by sample.

QC and filtering protocol
-------------------------

The following steps are performed:

 * **Initial quality screen:** apply quality filters to the barcode and
   UMI sequences and reject read pairs that fail to meet the criteria.

   - Barcode bases must have Q >= 10
   - UMI bases must have Q >= 30

   NB this filtering can be skipped data by specifying the
   ``--no-quality-filter`` option. This is recommended for NextSeq
   data 


 * **Barcode filtering:** barcodes are checked against the list of
   expected barcodes in the input "well list" file; read pairs that
   don't have an exact match are rejected.


 * **Read trimming:** ``cutadapt`` is used to perform the following
   trimming and filtering operations on the R2 read:

   - Remove sequencing primers
   - Remove poly-A/T and poly-N sequences
   - Apply quality filter of Q <= 25
   - Remove short reads (<= 20 bases) post-trimming

   NB if an R2 read fails any of the filters then the read pair is
   rejected.


 * **Poly-G region estimation:** ``cutadapt`` is used to look for
   R2 reads which contain poly-G regions. These reads are counted but
   no other action is taken; the step simply estimates the size of
   the effect in the data.

   NB this step is performed after the barcode filtering.


 * **Contamination screen:** ``fastq_screen`` is run to check the
   R2 reads against a set of mammalian and contaminant organisms, and
   exclude any read pairs where there is an exclusive match to the
   contaminants.

   The screen files must be explicitly supplied to the utility using
   the ``-m``/``--mammalian`` and ``-c``/``--contaminants`` options.

Reorganisation by barcode and sample
------------------------------------

At the end of the QC and filter pipeline the read pairs are
reorganised in two different ways:

 * **Reorganisation by barcode:** the read pairs are sorted into
   individual FASTQs according to their inline barcodes. This set of
   FASTQs forms the final outputs of the pipeline. Note that the
   number of files is equal to the number of barcodes (~1000).

 * **Reorganisation by sample:** the read pairs are sorted into FASTQs
   according to the sample name associated with the barcodes in the
   "well list" file.

This results in two fastq directories: ``fastqs.barcodes`` and
``fastqs.samples``. Note that the read pairs themselves are the same
in each set.

Reports
-------

The standard QC procedure is run on each set of FASTQS (barcodes and
samples) and reports are generated for each.

In addition the ``stats`` directory contains a summary of the read
and UMI counts after each stage (in TSV and XLSX format).

There is also a final summary report ``icell8_processing.html``.
