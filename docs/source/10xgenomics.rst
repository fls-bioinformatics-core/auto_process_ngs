Processing 10xGenomics Chromium SC 3'v2 single-cell data
========================================================

Background
----------

The 10xGenomics Chromium SC 3'v2 system prepares single cell (SC) samples
which are then sequenced as part of an Illumina sequencing run. 10xGenomics
provide the ``cellranger`` software package to perform Fastq generation
and subsequent analyses.

The auto-process package currently provides a utility script called
``process_10xgenomics.py`` to help with processing the 10xGenomics Chromium
data. This sits alongside the standard ``auto_process`` commands, and
wraps a subset of the ``cellranger`` commands whilst also providing a
degree of integration with the ``auto_process`` protocols.

The stages are:

1. :ref:`10xgenomics-fastq-generation`
2. :ref:`10xgenomics-set-up-analysis-dirs`
3. :ref:`10xgenomics-initial-single-library-analysis`

Each of these is covered in detail in the subsequent sections.

.. _10xgenomics-fastq-generation:

Fastq generation
----------------

The first step is to generate Fastq files from the raw sequencing data.
There are currently two possible scenarios:

1. The sequencing run only contains Chromium samples, or
2. The sequencing run contains a mixture of Chromium and non-Chromium
   samples.

The Fastq protocols for these two scenarios are slightly different and
are described in the following sections.

Sequencing run only contains Chromium samples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this case the initial step after setp is to fetch the primary data
for the run, by running the command::

    auto_process.py make_fastqs --only-fetch-primary-data

Once the data has been retrieved, the Fastq generation can be performed
using ``process_10xgenomics.py`` to run ``cellranger mkfastq``, e.g.::

    process_10xgenomics.py mkfastq \
        -s custom_SampleSheet.10xgenomics.csv \
        -r primary_data/170426_K00311_0033_AHJCY7BBXX \
        -o bcl2fastq

This will generate the Fastqs in the specified output directory (e.g.
``bcl2fastq``) along with an HTML report derived from the ``cellranger``
JSON QC summary file.

Finally, to generate statistics for the Fastqs run::

    auto_process.py update_fastq_statistics

Sequencing run contains both Chromium and non-Chromium samples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this case it is necessary to split the processing to handle the
Chromium and non-Chromium datasets separately. See the section
:ref:`problem-split-processing` for how to split the processing into
batches.

It is recommended to generate Fastqs for the lanes containing
'standard' (i.e. non-Chromium) datasets first, using the
``auto_process.py make_fastqs --lanes=...`` command.

When this has completed, the Fastqs can be generated for the lanes
with the Chromium datasets using the ``process_10xgenomics.py``
utility and specifying the lanes which contain these data via the
``-l/--lanes`` option. E.g.::

    process_10xgenomics.py mkfastq \
        -s custom_SampleSheet.10xgenomics.csv \
        -r primary_data/170426_K00311_0033_AHJCY7BBXX \
	-l 5,6 \
        -o bcl2fastq.chromium

This will generate the Fastqs in the specified output directory (e.g.
``bcl2fastq.chromium``) along with an HTML report derived from the
``cellranger`` JSON QC summary file. It will also update the existing
``projects.info`` file with entries for the Chromium projects.

Finally, to update the statistics with the Chromium Fastqs, run e.g.::

    auto_process.py update_fastq_stats \
         --unaligned-dir bcl2fastq.chromium \
         --add

``--add`` is necessary to combine the statistics for the Chromium data
with those from the standard datasets.

.. _10xgenomics-set-up-analysis-dirs:

Set up project analysis directories and run QC
----------------------------------------------

Regardless of the scenario described previously, the ``projects.info``
file should contain initial entries for both the Chromium and
non-Chromium datasets (if any).

.. note::

   The ``projects.info`` can also be updated using the
   ``update_projects`` command::

       process_10xgenomics.py update_projects \
           -u bcl2fastq.chromium

   which will add entries for the Chromium projects and samples
   to the ``projects.info`` file.

Once the file has been edited to include additional data on the user,
PI etc, the ``setup_analysis_dirs`` command can be used to create
analysis directories for each of the projects.

For runs which consist solely of Chromium datasets, do e.g.::

    auto_process.py setup_analysis_dirs \
        --unaligned_dir=bcl2fastq.chromium

For runs with a mixture of Chromium and non-Chromium datasets it's
necessary to run this command twice: once to set up projects for the
non-Chromium samples and once to set up to the Chromium projects.
E.g.::

    # Set up non-Chromium-based projects
    auto_process.py setup_analysis_dirs

    # Set up Chromium-based projects
    auto_process.py setup_analysis_dirs \
        --unaligned_dir=bcl2fastq.chromium \
        --undetermined=undetermined.chromium

.. note::

   The ``--undetermined`` option is required to also create a 'project'
   specifically for the undetermined reads from the Chromium data.

Once the projects are set up, the standard QC pipeline can be run
using the ``run_qc`` command::

       auto_process.py run_qc

.. _10xgenomics-initial-single-library-analysis:

Perform initial single-library analysis
---------------------------------------

Initial single-library analysis can be performed by using
``process_10xgenomics.py`` to run ``cellranger count`` on all the
samples in the Chromium-based projects, e.g.::

       process_10xgenomics.py count \
           -t .../refdata-cellranger-mm10-1.2.0
	   PROJECT1 PROJECT2

(where ``PROJECT1`` and ``PROJECT2`` represent the projects with Chromium
datasets).

   This generates the 'count' output in a temporary directory and copies
   the ``web_summary.html`` files into the project directories on
   completion.

   It also creates an index file and a ZIP archive with the HTML summary
   reports from ``cellranger count``. These reports are copied by the
   ``publish_qc`` command, along with the standard QC reports.

.. note::

   See https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count
   for more details on the single-library analysis.

.. _10xgenomics-requirements:

Requirements
------------

The ``process_10xgenomics.py`` utility requires that the ``cellranger``
package be installed and available in the environment (for example, a
suitable environment module can be specified via the ``--modulefiles``
option.)

See https://support.10xgenomics.com/single-cell/software/pipelines/latest/what-is-cell-ranger
for more information on installing and running ``cellranger``.

.. _10xgenomics-outputs:

Outputs and reports
-------------------

The ``process_10xgenomics.py mkfastq`` command generates a report in the
top-level analysis directory called ``cellranger_qc_summary[_LANES].html``,
which is an HTML copy of the QC summary JSON file produced by
``cellranger mkfastq`` (nb ``LANES`` will be the subset of lanes from the
run which contained the Chromium data, if the run consisted of a mixture of
Chromium and non-Chromium samples; , for example: ``--lanes=5,6`` results
in ``56``).

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

Appendix: known issues
----------------------

It has been observed that when the Fastq files produced by the ``mkfastq``
command have very low read counts then the single-library analyses may
fail, with ``cellranger count`` reporting an error of the form e.g.::

    Could not auto-detect Single Cell 3' chemistry. Fraction of barcodes
    on whitelist was at best 0.23%, while we expected at least 10.00% for
    one of the chemistries.

There is currently no workaround for this issue.

Appendix: options for 'process_10xgenomics.py'
----------------------------------------------

.. _10xgenomics-mkfastq-options:

'mkfastq' options
~~~~~~~~~~~~~~~~~

The ``mkfastq`` command supports the following options::

    -s : specify the sample sheet to use
    -r : specify the location of the primary data for the run
    -l : optionally, specify the lane numbers
    -o : specify the output directory

See also :ref:`10xgenomics-additional-options`.

.. _10xgenomics-count-options:

'count' options
~~~~~~~~~~~~~~~

 The ``count`` command supports the following options::

    -u : specify the name of the output directory from 'mkfastq'
    -t : specify the path to the directory with the appropriate
         10xGenomics transcriptome data

See also :ref:`10xgenomics-additional-options`.

.. _10xgenomics-additional-options:

Additional options for 'process_10xgenomics.py'
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``process_10xgenomics.py`` has a number of additional options for
controlling how the ``cellranger`` pipeline is run::

    --jobmode JOB_MODE : job mode to run cellranger in
    --mempercore MEM_PER_CORE : memory assumed per core (in Gbs)
    --maxjobs MAX_JOBS : maxiumum number of concurrent jobs to run
    --jobinterval JOB_INTERVAL : how often jobs are submitted (in ms)

These map onto the equivalent ``cellranger`` options.

There are also the following general options::

   --modulefiles MODULEFILES : comma-separated list of environment
                               modules to load before executing commands


