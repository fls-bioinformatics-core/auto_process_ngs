QC protocol specification
=========================

--------
Overview
--------

Different QC protocol definitions within ``auto_process``
are used to run the appropriate metrics for different
types of data within the QC pipeline, by specifying
which metrics should be used on which reads within
a dataset.

Each protocol is defined by the following attributes:

=================== ============================
Attribute           Description
=================== ============================
Name                Short name for the protocol
Description         Longer free-text description 
Sequence data reads Reads with biological data
Index data reads    Reads with index data
QC modules          List of QC modules and arguments
=================== ============================

-----------------------
Protocol specifications
-----------------------

The attributes for a protocol can be specified using a
string of the form:

::

   NAME:DESCRIPTION:seq_reads=SEQ_READS:index_reads=INDEX_READS:qc_modules=QC_MODULES

where:

=============== ================== ======================
Field           Description        Example
=============== ================== ======================
``NAME``        Protocol name      ``standardPE``
``DESCRIPTION`` Description text   ``Standard paired-end data``
``SEQ_READS``   List of reads      ``[r1,r2]``
``INDEX_READS`` List of reads      ``[r1]``
``QC_MODULES``  List of QC modules ``[fastqc,fastq_screen]``
=============== ================== ======================

For example:

::

   basicQC:'Basic QC protocol':seq_reads=[r1,r2]:index_reads=[]:qc_modules=[sequence_lengths,fastqc]

The name and description must always be the first two fields
of a protocol specification. The ``seq_reads``, ``index_reads``
and ``qc_modules`` fields can appear in any order and can be
omitted if empty. The same read ID should not appear in both
``seq_reads`` and ``index_reads``.

-------------------------------
Protocol names and descriptions
-------------------------------

Names cannot contain spaces or colons; descriptions can
contain any characters (including spaces).

Optionally, descriptions can be enclosed in matching quotes
(either single or double); these should be used if the
description includes special characters such as colons or
quotes.

-----------------------------
Sequence data and index reads
-----------------------------

QC protocols distinguish between "sequence data" and "index"
reads as follows:

* Sequence data reads are those reads containing biologically
  significant data, and are used for mapped metrics;
* Index reads contain non-biologically significant data (e.g.
  feature barcodes). Non-mapped metrics are run on both
  sequence data and index data reads.

Reads can optionally also specify subsequences to use for
QC, for example for 10x Genomics Flex data only the first
50 bases of R1 contains the biologically significant data.
Subsequences can be specified by ranges attached to reads,
for example in this case ``seq_reads=[r2:1-50]``, and tell
the pipeline to only use these portions of the the reads
when running the mapped metrics.

----------
QC modules
----------

QC modules specific which software will be run and which
metrics will be produced.

The following QC modules are available:

============================== =============================== ==========
QC module                      Software                        Metrics
============================== =============================== ==========
``fastqc``                     FastQC                          General stats and quality information
``fastq_screen``               FastqQScreen                    Screens against panels of organisms
``picard_insert_size_metrics`` Picard CollectInsertSizeMetrics Insert sizes
``qualimap_rnaseq``            Qualimap 'rnaseq'               Coverage and genomic origin of reads
``rseqc_genebody_coverage``    RSeQC geneBody_coverage.py      Coverage
``sequence_lengths``           Built-in                        Sequence length and masking stats
``strandedness``               fastq_strand.py                 Strandedness information
``cellranger_count``           CellRanger 'count'              Single library analysis
``cellranger-atac_count``      CellRanger-ATAC 'count'         Single library analysis
``cellranger-arc_count``       Cellranger-ARC 'count'          Single library analysis
``cellranger_multi``           CellRanger 'multi'              Multiplexing analysis
============================== =============================== ==========

Some QC modules allow additional optional arguments to be
specified, to modify either the metric or pipeline behaviour,
or how the metric validation is performed.

Arguments are specified using the general syntax:

::

   QC_MODULE(ARG1=VALUE1;ARG2=VALUE2;...)

For example:

::

   cellranger_count(set_cell_count=false;set_metadata=false)

The following options are recognised:

========================================== =================== ===========
Option                                     QC modules          Description
========================================== =================== ===========
``chemistry=CHEMISTRY``                    cellranger[*]_count Explicitly set chemistry when running CellRanger 'count' (e.g. ``ARC-v1``)
``library=LIBRARY_TYPE``                   cellranger[*]_count Explicitly set library type in pipeline tasks (e.g. ``snRNA-seq``)
``cellranger_version=VERSION``             cellranger[*]_count Set expected CellRanger version for validation (default is the version identified by the pipeline; use ``*`` to match any version)
``cellranger_refdata=REFDATA``             cellranger[*]_count Set expected CellRanger reference dataset for validation (default is the reference identified by the pipeline; use ``*`` to match any reference dataset)                             
``set_cell_count=true|false``              cellranger[*]_count Whether to use outputs to set the cell count (default is ``true``)
``set_metadata=true|false``                cellranger[*]_count Whether to set metadata from this module (default is ``true``)
``cellranger_use_multi_config=true|false`` cellranger_count    If ``true`` then use CellRanger 'multi' config file to identify samples with biological data (default is ``false``)
========================================== =================== ===========

