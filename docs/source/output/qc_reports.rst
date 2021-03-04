QC Reports
==========

The :doc:`auto_process run_qc <../using/run_qc>` command outputs an
HTML report for the QC for each of the projects in the analysis
directory, which enables an assessment of the quality of the Fastq
files for each sample in the project and tries to highlight aspects
which might pose problems with the data in subsequent analyses.

An example of the top of a QC report index page is shown below:

.. image:: ../images/qc/qc_report_full.png
   :align: center

The report consists of:

* :ref:`qc_report_project_metadata`
* :ref:`qc_report_qc_summary_table`
* :ref:`qc_report_qc_outputs_per_fastq`

.. _qc_report_project_metadata:

************************
Project metadata summary
************************

The project metadata table summarises information associated with the
project, including the user, PI, library type, organisms and QC
protocol.

.. _qc_report_qc_summary_table:

****************
QC summary table
****************

The QC summary table summarises the key QC results for each sample
and Fastq or Fastq pair (for paired end data).

For example:

.. image:: ../images/qc/qc_report_summary.png
   :align: center

The summary includes numbers of reads or read pairs and the sequence
length ranges for each Fastq; other metrics are represented by small
plots:

The following data are shown:

* :ref:`qc_report_quality_boxplots`
* :ref:`qc_report_fastqc_summary_plots`
* :ref:`qc_report_fastq_screen_summary_plots`
* :ref:`qc_report_strandedness`

One purpose of this table is to help pick up on trends and identify
any outliers within the dataset as a whole; hence the main function
of these plots are to convey a general sense of the data.

Note that not all outputs might appear, depending on the
:doc:`QC protocol <../using/run_qc>` that was used.

The sample and Fastq names in the table link through to the
full QC outputs for the sample or Fastqs in question; other items
(e.g. the quality boxplots) link to the relevant parts of the full
QC outputs section (see :ref:`qc_report_qc_outputs_per_fastq`).

An additional summary table may appear after this one with details
of outputs from 10xGenomics single library analyses (see
:ref:`qc_report_single_library_analyses`).

.. note::

   In earlier versions of the QC reports, links to single library
   analyses were appended directly to the main summary table, and
   no separate sigle library analyses table was present.

.. _qc_report_quality_boxplots:

Quality boxplots
----------------

The summary table includes a small version of the sequence quality
boxplot from ``fastqc``, for example:

.. image:: ../images/qc/uboxplot.png
   :align: center

A larger version of the plot is presented in the
:ref:`qc_report_qc_outputs_per_fastq` section.

.. _qc_report_fastqc_summary_plots:

Fastqc summary plots
--------------------

The output from ``fastqc`` includes a summary table with a set of
metrics and an indication of whether the Fastq has passed, failed
or triggered a warning for each.

The summary table includes a small plot which gives an impression of
the overall state of the metrics for each Fastq file, for example:

.. image:: ../images/qc/fastqc_uplot.png
   :align: center

Each bar in the plot represents one of the ``fastqc`` metrics,
(for example "Basic statistics", "Per base sequence quality", and
so on); the colour (red, amber, green) and position (left, centre,
right) indicate the status of the metric as determined by
``fastqc``.

The data are presented in more detail in a table in the
:ref:`qc_report_qc_outputs_per_fastq` section.

.. _qc_report_fastq_screen_summary_plots:

Fastq_screen summary plots
--------------------------

The summary table includes a small plot which represents the
outputs from ``fastq_screen``, for example:

.. image:: ../images/qc/fastq_screen_uplot.png
   :align: center

The three boxes represent (from left to right) the model organisms,
other organisms and rRNA plots produced by ``fastq_screen``. The
full plots and links to the raw data for each screen can be found
in the :ref:`qc_report_qc_outputs_per_fastq` section.

.. _qc_report_strandedness:

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

.. table:: Example strandedness plots
   :widths: auto

   ======================= ===============
   Strandedness            Example
   ======================= ===============
   Likely forward stranded .. image:: ../images/qc/strandedness_forward.png
   Likely reverse stranded .. image:: ../images/qc/strandedness_reverse.png
   Likely unstranded       .. image:: ../images/qc/strandedness_no_strand.png
   ======================= ===============

More detailed information about the strandedness statistics
is given in the :ref:`qc_report_qc_outputs_per_fastq` section.

.. _qc_report_single_library_analyses:

Single library analyses
-----------------------

For 10xGenomics datasets single library analyses may also have
been performed for each sample using the ``count`` command of the
appropriate 10xGenomics pipeline (e.g. ``cellranger`` for scRNA-seq
data, ``cellranger-atac`` for scATAC-seq data etc).

In these cases an additional summary table will appear in the report
with appropriate metrics for each sample along with links to the HTML
reports from the ``count`` command. For example, for an scRNA-seq
dataset:

.. image:: ../images/qc/qc_report_single_library_summary.png
   :align: center

The reported metrics will depend on the pipeline and type of data.

Details of the contents of the linked ``web_summary.html`` report can
be found in the appropriate documentation for the 10xGenomics pipeline:

 * ``cellranger``: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/summary
 * ``cellranger-atac``: https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/summary
 * ``cellranger-arc``: https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/output/summary

.. note::

   The full set of outputs can be found under the ``cellranger_count``
   subdirectory of the project directory, when single library
   analysis has been performed.
   
.. _qc_report_qc_outputs_per_fastq:

*************************
Full QC outputs per Fastq
*************************

After the summary table, the full QC outputs for each Fastq or
Fastq pair are grouped by sample, for example:

.. image:: ../images/qc/qc_outputs_per_fastq.png
   :align: center

For each Fastq the subsections consist of:

* ``fastqc`` outputs including the sequence quality boxplot
  and a table of the quality metrics with links to the full
  report:

  .. image:: ../images/qc/fastqc_full.png

* ``fastq_screen`` outputs for each screen, for example:

  .. image:: ../images/qc/fastq_screen_full.png

* ``fastq_strand`` data:
  
  .. image:: ../images/qc/strandedness_full.png
