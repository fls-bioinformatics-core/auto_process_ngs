QC pipeline
===========

.. _qc_pipeline_overview:

Overview
--------

The standard QC pipeline is run on each project in a sequencing
run, and in turn runs the following external packages for each
Fastq file in the project:

 * ``fastqc`` (for general quality metrics)
 * ``fastq_screen`` for a set of genome indexes (to confirm that
   sequences are from the expected organisms, and check for any
   contamination)

In addition:

 * ``fastq_strand`` is run either pair-wise (for paired end data),
   or on R1 (for single-end data) or R2 (for single-cell data from
   ICELL8 and 10xGenomics Chromium platforms), to determine the
   strandedness of the sequence data
 * ``multiQC`` is run to summarise the QC from all the Fastqs in
   the project

A HTML report is also generated for the QC outputs, and is described
in the sections :ref:`qc_pipeline_qc_reports` below.

.. _qc_pipeline_qc_reports:

QC reports
----------

An example of the top of a QC report index page is shown below:

.. image:: images/qc/qc_report_full.png

The first part consists of two tables: the first summarises any metadata
associated with the project, and the second attempts to summarise the
key QC results for each Fastq, each Fastq pair (for paired end data),
and each sample.

For example:

.. image:: images/qc/qc_report_summary.png

The purpose of this table is to help pick up on trends and identify any
outliers within the dataset as a whole; hence it uses small plots to
convey a general sense of the data.

Quality boxplots
----------------

The summary table includes a small version of the sequence quality
boxplot from ``fastqc``, for example:

.. image:: images/qc/uboxplot.png

Fastqc summary plots
--------------------

The output from ``fastqc`` includes a summary table with a set of
metrics and an indication of whether the Fastq has passed, failed
or triggered a warning for each.

The summary table includes a small plot which gives an impression of
the overall state of the metrics for each Fastq file, for example:

.. image:: images/qc/fastqc_uplot.png

Fastq_screen summary plots
--------------------------

The summary table includes a small plot which summarises the
outputs from ``fastq_screen``, for example:

.. image:: images/qc/fastq_screen_uplot.png

Strandedness
------------

``fastq_strand.py`` runs ``STAR`` to get the number of reads which
map to the forward and reverse strands; it then calculates a
pseudo-percentage ("pseudo" because it can exceed 100%) for foward
and reverse.

The summary table reports the pseudo-percentages as a barplot with
a pair of barplots, where the top bar represents the forward
pseudo-percentage and the bottom bar the reverse value.

Some examples:

Likely forward strand:

.. image:: images/qc/strandedness_forward.png


Likely reverse strand:

.. image:: images/qc/strandedness_reverse.png


Likely unstranded:

.. image:: images/qc/strandedness_no_strand.png
