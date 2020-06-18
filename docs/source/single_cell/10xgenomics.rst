Processing 10xGenomics Chromium SC 3'v2 single-cell data
========================================================

Background
----------

The 10xGenomics Chromium SC 3'v2 system prepares single cell (SC) samples
which are then sequenced as part of an Illumina sequencing run. 10xGenomics
provide the ``cellranger`` and ``cellranger-atac`` software packages to
perform Fastq generation and subsequent analyses:

* ``cellranger`` is used for single cell RNA-seq data
* ``cellranger-atac`` is used for single cell ATAC-seq data

The auto-process package currently provides a utility script called
``process_10xgenomics.py`` which wraps a subset of the ``cellranger``
and ``cellranger-atac`` commands, whilst also providing a degree of
integration with the ``auto_process`` pipeline.

Processing protocol for 10xGenomics Chromium data
-------------------------------------------------

The recommended steps are:

1. Generate the Fastqs as described in
   :ref:`make_fastqs-10x_chromium_sc-protocol` or
   :ref:`make_fastqs-10x_atac-protocol`
2. Set up analysis directories and run initial QC as per the standard
   protocol
3. Perform initial single library analysis by running the
   ``process_10xgenomics.py`` utility, as described in
   :ref:`10xgenomics-initial-single-library-analysis`

.. _10xgenomics-initial-single-library-analysis:

Perform initial single-library analysis
---------------------------------------

Initial single-library analysis can be performed by using
``process_10xgenomics.py count`` (for scRNA-seq data) or
``process_10xgenomics.py count-atac`` (for scATAC-seq data) commands.
These are wrappers which run ``cellranger count`` or
``cellranger-atac count`` on all the samples in the Chromium-based
projects.

.. _10xgenomics-count-options:

Single-library analysis for scRNA-seq data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The general invocation for running the single-library analysis on
single cell RNA-seq data is:

::

       process_10xgenomics.py count \
           -t /path/to/refdata-cellranger-mm10-1.2.0
	   PROJECT1 [ PROJECT2 ... ]

where ``PROJECT1`` etc represent the projects with Chromium
RNA-seq datasets.

The ``count`` command supports the following options::

    -u : specify the name of the output directory from 'mkfastq'
    -t : specify the path to the directory with the appropriate
         10xGenomics transcriptome data
    -c : specify the assay configuration (aka chemistry)

in addition to the options outlined in the section
:ref:`10xgenomics-additional-options`.

.. note::

   If a project metadata defines an organism name which matches one
   of the entries in the ``10xgenomics_transcriptomes`` section of
   the configuration file then the ``-t`` option isn't required;
   instead the matching reference data will be used automatically
   for the single-library analysis.

See https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count
for more details on the single-library analysis for scRNA-seq.

.. _10xgenomics-count-atac-options:

Single-library analysis for sATAC-seq data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The general invocation for running the single-library analysis on
single cell ATAC-seq data is:

::

       process_10xgenomics.py count-atac \
           -r /path/to/refdata-cellranger-atac-mm10-1.0.1
	   PROJECT1 [ PROJECT2 ... ]

where ``PROJECT1`` etc represent the projects with Chromium
ATAC-seq datasets.

The ``count-atac`` command supports the following options::

    -u : specify the name of the output directory from 'mkfastq'
    -r : specify the path to the directory with the appropriate
         10xGenomics ATAC genome reference data

in addition to the options outlined in the section
:ref:`10xgenomics-additional-options`.

.. note::

   If a project metadata defines an organism name which matches one
   of the entries in the ``10xgenomics_atac_genome_references``
   section of the configuration file then the ``-r`` option isn't
   required; instead the matching reference data will be used
   automatically for the single-library analysis.

See https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/using/count
for more details on the single-library analysis for scATAC-seq.

.. _10xgenomics-additional-options:

Options for controlling cellranger
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``process_10xgenomics.py`` has a number of additional options for
controlling how the ``cellranger`` pipeline is run::

    --jobmode JOB_MODE : job mode to run cellranger in
    --jobinterval JOB_INTERVAL : how often jobs are submitted (in ms)
    --maxjobs MAX_JOBS : maxiumum number of concurrent jobs to run

The default ``JOB_MODE`` is ``local``, which runs the ``cellranger``
pipelines on the local system. In this case the following additional
options can be used to control the resources used by the pipeline::

    --localcores LOCAL_CORES : maximum number of cores the pipeline
                               can request
    --localmem LOCAL_MEM     : maximum memory the pipeline can
                               request (in Gbs)

If ``cellranger`` is configured to use additional job submission
systems (e.g. Grid Engine) then ``JOB_MODE`` can specify one of these
(e.g. ``sge``). In this case the following additional options can
be used::

    --mempercore MEM_PER_CORE : memory assumed per core (in Gbs)

Note that all the above options  map onto the equivalent ``cellranger``
options; there are also the following general non-``cellranger`` options::

   --modulefiles MODULEFILES : comma-separated list of environment
                               modules to load before executing commands

.. _10xgenomics-outputs:

Outputs and reports
*******************

After running the ``process_10xgenomics.py counts`` command, the project
directory will contain the following output directories:

 ========================== =================================================
 **Directory**              **Description and contents**
 -------------------------- -------------------------------------------------
 ``fastqs``                 FASTQs from ``cellranger mkfastq``/``bcl2fastq``
 ``qc``                     The standard QC outputs
 ``cellranger_fastq_path``  Bcl2fastq-like directory with links to FASTQs
 ``cellranger_count``       Single-library analyses from ``cellranger count``
 ========================== =================================================

The ``cellranger_count`` directories each further contain one
subdirectory for each sample, within which there is the ``outs``
directory produced by ``cellranger_count``.

.. note::

   By default these ``outs`` directories only contain the
   ``web_summary.html`` files; to collect all the outputs from
   ``cellranger count`` (i.e. the ``.cloupe``, ``BAM``, and gene
   matrix files required for subsequent analyses), use the
   ``--all-outputs`` option.

The ``cellranger_fastq_path`` directory is a facsimile of the bcl2fastq
output directory produced by ``cellranger mkfastq``, which can be supplied
as the input to one of the ``cellranger`` analysis commands if desired.

The directory will also contain:

 * The report from ``cellranger count`` (``cellranger_count_report.html``)
   which links to the ``web_summary.html`` file for each sample
 * A ZIP archive file with the report plus the summaries for each sample,
   for viewing elsewhere
 * A ``README.info`` file

Troubleshooting
***************

Single-library analyses fail for low read counts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It has been observed that when the Fastq files produced by the ``mkfastq``
command have very low read counts then the single-library analyses may
fail, with ``cellranger count`` reporting an error of the form e.g.::

    Could not auto-detect Single Cell 3' chemistry. Fraction of barcodes
    on whitelist was at best 0.23%, while we expected at least 10.00% for
    one of the chemistries.

There is currently no workaround for this issue.

Single-library analyses fail to detect chemistry automatically
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default ``cellranger count`` attempts to determine the chemistry used
automatically, however this may fail if a low number of reads map to the
reference genome and give an error of the form::

    The chemistry was unable to be automatically determined. This can
    happen if not enough reads originate from the given reference. Please
    verify your choice of reference or explicitly specify the chemistry
    via the --chemistry argument.

If the reference data being used is correct then use the ``--chemistry``
option to specify the appropriate assay configuration - see
https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count
