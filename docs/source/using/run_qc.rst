Generating QC reports using ``auto_process run_qc``
===================================================

--------
Overview
--------

The ``run_qc`` command is used to run the QC pipeline on the
Fastqs for each project within the analysis directory, and
generate a report for each one.

The general invocation of the command is:

::

   auto_process.py run_qc *options* [ANALYSIS_DIR]

There is also a standalone utility called ``run_qc.py`` which
can be used to run the QC pipeline on an arbitrary subdirectory.

------------
QC protocols
------------

QC protocols define sets of "QC modules" to be run by the QC
pipeline (see below), as well as "sequence data reads" (reads that
contain biological data) and "index data reads" (reads that
contain index data).

The QC pipeline protocol used for each project will differ slightly
depending on the nature of the data within that project:

.. include:: ../auto/qc_protocols.rst

By default the protocol is determined automatically for each project,
based on the metadata associated with that project.

.. note::

   The automatic protocol assignment can be overridden on a
   per-project basis using the ``--protocol`` option:

   ::

      auto_process.py run_qc --protocol PROJECT=QCPROTOCOL ...

   The ``--protocol`` option can be specified multiple times
   (once for each project where the protocol should be set
   explicitly)

----------
QC modules
----------

The following QC modules are available within the QC pipeline:

============================== ======================
QC module                      Details
============================== ======================
``fastqc``                     Runs `fastqc`_ for general quality metrics
``fastq_screen``               Runs `fastq_screen`_ for a set of genome
                               indexes, to verify sequences correspond to
                               the expected organism (and check for
                               contaminants)
``sequence_lengths``           Examines distribution of sequence lengths
                               and checks for padding and masking
``strandedness``               Runs `fastq_strand`_ on the appropriate
                               reads to indicate the strand-specificity of
                               the sequence data (requires appropriate
			       ``STAR`` indexes)
``picard_insert_size_metrics`` Runs `picard`_ ``CollectInsertSizeMetrics``
                               on paired data to determine insert sizes
``rseqc_genebody_coverage``    Runs `rseqc`_ ``geneBody_coverage.py`` to
                               generate gene body coverage plot for all
			       samples
``rseqc_infer_experimemt``     Run `rseqc`_ ``infer_experiment.py`` to
                               determine strand-specificity
``qualimap_rnaseq``            Runs `qualimap`_ ``rnaseq`` to generate
                               various metrics including coverage and
			       genomic origin of reads
``cellranger_count``           Single library analysis for each sample using
                               `cellranger`_ ``count``
``cellranger-atac_count``      Single library analysis for each sample using
                               `cellranger_atac`_ ``count``
``cellranger-arc_count``       Single cell multiome analysis using
                               `cellranger_arc`_ ``count`` (requires
                               :doc:`10x_multiome_libraries.info <../control_files/10x_multiome_libraries_info>`)
``cellranger_multi``           Cell multiplexing and fixed RNA profiling
                               analyses using `cellranger`_ ``multi``
                               (requires
                               :doc:`10x_multi_config[.SAMPLE].csv <../control_files/10x_multi_config_csv>`)
============================== ======================

Appropriate reference data must be available (for example,
``STAR`` indexes or 10x Genomics reference datasets), and
certain additional files (noted in the QC module descriptions
above) may be required for some of the QC modules to run;
typically the modules are skipped when appropriate reference
data is not available.

.. _fastqc:  http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
.. _fastq_screen: http://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/
.. _fastq_strand: https://genomics-bcftbx.readthedocs.io/en/latest/reference/qc_pipeline.html#fastq-strand
.. _picard: https://gatk.broadinstitute.org/hc/en-us/articles/360037055772-CollectInsertSizeMetrics-Picard-
.. _rseqc: http://rseqc.sourceforge.net/#
.. _qualimap: http://qualimap.conesalab.org/doc_html/command_line.html#rna-seq-qc
.. _cellranger: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger
.. _cellranger_atac: https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/what-is-cell-ranger-atac
.. _cellranger_arc: https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/what-is-cell-ranger-arc
.. _multiqc: http://multiqc.info/

On successful completion of the pipeline an HTML report is
generated for each project; these are described in
:doc:`QC reports <../output/qc_reports>`. By default `multiQC`_
is also run as part of the reporting.

If a QC server has been specified in the configuration then the
reports can be copied there for sharing using the
:doc:`publish_qc command <publish_qc>`.

.. note::

   The QC pipeline can be run outside of the ``auto_process``
   pipeline by using the ``run_qc.py`` utility; see the
   section on :doc:`running the QC standalone <run_qc_standalone>`.

---------------------------------------------
Biological versus non-biological data samples
---------------------------------------------

For some types of dataset (e.g. 10x Genomics CellPlex data), not
all samples in the dataset contain biological data (for example,
CellPlex datasets also have "multiplexing capture" samples which
contain feature barcodes).

Biological and non-biological data samples may be identified
implicitly (for example, by using the library type information
in the ``10x_multi_config.csv`` file for CellPlex datasets).
Alternatively samples with biological data can be explicitly
defined in the ``Biological samples`` field of the ``README.info``
metadata file in the analysis project directory, as a
comma-separated list of sample names. For example:

::

   Biological samples    SMPL1,SMPL2

Samples in the project which are not in this list are treated as
containing non-biological data; if no samples are listed then all
samples are assumed to contain biological data.

When biological and non-biological samples are differentiated,
the pipeline will only run and report "mapped" metrics (for
example screens, strandedness, gene body coverage etc) for the
biological samples; these metrics will be omitted for non-biological
samples (even if they have been specified as part of the QC
protocol).

--------------------------------------
QC metrics using subsequences in reads
--------------------------------------

Some metrics can be applied to subsequences within reads (rather
than the whole sequence), if a range of bases is defined within
the protocol.

Specifically, where subsequences are specified for sequence data
reads (i.e. reads containing biological data) then the "mapped"
QC modules (for example, FastqScreen or strandedness) will only
use those subsequences.

(See the :doc:`QC protocol specification <../reference/qc_protocol_specification.rst>`
documentation for the subsequence specification syntax.)

-------------------------------------------
Per-lane QC metrics for undetermined Fastqs
-------------------------------------------

By default when handling Fastqs for the ``undetermined``
project, the pipeline runs in mode whereby it generates
copies of the input Fastqs for each lane that appears in the
read headers of each Fastq, and then run the QC on those
per-lane Fastqs (rather the originals, which are not changed).

This results in per-lane QC metrics, which can be useful for
diagnostic purposes when handling Fastqs which have been merged
across multiple lanes (for example, to determine whether
contanimation is confined to a single lane).

The behaviour is controlled by the ``split_undetermined_fastqs``
setting in the ``qc`` section of the configuration file (see
:ref:`qc_pipeline_configuration`).

.. _run_qc_including_external_outputs:

-----------------------------------------
Including external (non-pipeline) outputs
-----------------------------------------

It is possible to include links within the QC reports to
additional output files produced outside of the pipeline,
by specifying them within an ``extra_outputs.tsv`` file.

If this file is present in the QC directory being reported
then any additional external files which it specifies
will be linked from the "Extra outputs" section in the
report. The specified files will also be included in the
ZIP archive.

.. note::

   The additional files must also be within the QC
   directory.

Each line in the TSV file specifies an output file to
include, with the minimal specification being:

::

   FILE_PATH     DESCRIPTION

For example:

::

   extra_metrics/metrics.html    Manually generated QC metrics

If the output has additional associated files or directories
(for example, the main output is an index file that then
links to other files) then these can be included as a
comma-separated list via an optional third field in the TSV
file:

::

   FILE_PATH     DESCRIPTION     PATH1[,PATH2[,...]]

For example:

::

   more_metrics/index.html    More QC metrics    more_metrics/results

The additional files will then be included in the ZIP archive.

.. note::

   If any additional "files" are actually subdirectories
   then the contents of those directories will also be
   included automatically.

------------------
Additional options
------------------

For 10xGenomics datasets, the following options can be used to
override the defaults defined in the configuration:

* ``--cellranger``: explicitly sets the path to the ``cellranger``
  (or other appropriate 10xGenomics package)
* ``--10x_force_cells``: explicitly specify the number of cells,
  overriding automatic cell detection algorithms (i.e. set the
  ``--force-cells`` option for CellRanger)
* ``--10x_extra_projects``: specify additional project directories
  to fetch Fastqs from when running single library analyses (i.e.
  add the Fastq directory paths for each project to the
  ``--fastqs`` option for CellRanger)
