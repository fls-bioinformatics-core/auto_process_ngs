10xGenomics Chromium SC 3'v2
============================

The utility script ``process_10xgenomics.py`` has been added to help with
processing single cell (SC) data from 10xGenomics Chromium SC 3'v2 system,
alongside the standard ``auto_process`` commands.

(``process_10xgenomics.py`` is essentially a wrapper for a subset of the
commands provided by the 10xGenomics' ``cellranger`` pipeline.)

The suggested protocol is:

0. **Generate Fastqs for 'standard' datasets**

   Assuming that the sequencing run includes a mixture of 'standard' (i.e.
   non-Chromium data) alongside the Chromium data, it's recommended that
   the Fastqs are generated for these standard datasets using the
  ``auto_process.py make_fastqs`` command.

   See :ref:`problem-split-processing` for how to split the
   processing into batches.

1. **Generate Fastqs for the Chromium data**

   Use ``process_10xgenomics.py`` to run ``cellranger mkfastq`` e.g.::

       process_10xgenomics.py mkfastq \
           -s custom_SampleSheet.10xgenomics.csv \
           -r primary_data/170426_K00311_0033_AHJCY7BBXX \
           -l 5,6 \
           -o bcl2fastq.10xgenomics

   This also generates an HTML report called ``cellranger_qc_summary_56.html``
   from the ``cellranger`` JSON QC summary file, which will be included by
   the standard ``publish_qc`` command.

   The relevant options are::

       -s : specify the sample sheet to use
       -r : specify the location of the primary data for the run
       -l : optionally, specify the lane numbers
       -o : specify the output directory

   (See also :ref:`10xgenomics-additional-options`.)

2. **Generate statistics**

   Use the ``update_fastq_stats`` command::

       auto_process.py update_fastq_stats \
           --unaligned_dir bcl2fastq.10xgenomics \
           --add

   ``--add`` is necessary to combine the statistics for the Chromium data
   with those from the standard datasets.

3. **Set up analysis directories**

  Use the ``setup_analysis_dirs`` command::

      auto_process.py setup_analysis_dirs \
          --unaligned_dir=bcl2fastq.10xgenomics \
          --undetermined=undetermined.10xgenomics

  This will create analysis directories for projects listed in
   ``projects.info`` which also appear in the referenced 'unaligned'
   directory.

  The ``--undetermined`` option is required to also create a 'project'
  for the undetermined reads from the Chromium data.

4. **Run 'cellranger count' on all samples**

   Use ``process_10xgenomics.py`` to run ``cellranger count`` on all
   samples::

       process_10xgenomics.py count \
           -u bcl2fastq.10xgenomics \
           -t .../refdata-cellranger-mm10-1.2.0

   This generates the 'count' output in a temporary directory and copies
   the relevant files into the project directories on completion.

   It also creates an index file and a ZIP archive with the HTML summary
   reports from ``cellranger count``. These reports are copied by the
   ``publish_qc`` command, along with the standard QC reports.

   The relevant options are::

       -u : specify the name of the output directory from 'mkfastq'
       -t : specify the path to the directory with the appropriate
            10xGenomics transcriptome data

   (See also :ref:`10xgenomics-additional-options`.)

4. **Run the standard QC pipeline**

   Use the ``run_qc`` command::

      auto_process.py run_qc [ --projects=10xGenomics ]

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

Outputs
-------

 * The ``process_10xgenomics.py mkfastq`` command produces an HTML copy
   of the QC summary JSON file produced by ``cellranger mkfastq``,
   called ``cellranger_qc_summary[_LANES].html`` (where ``LANES`` is a
   list of the lane numbers specified when running the command, for
   example: ``--lanes=5,6`` results in ``56``).

 * The ``process_10xgenomics.py count`` command produces a subdirectory
   called ``cellranger_count`` in analysis directories where there is
   Chromium data.

   This contains one subdirectory for each sample, within which there is
   the ``outs`` directory produced by ``cellranger_count``. These ``outs``
   directories contain the ``.cloupe``, ``BAM`` and gene matrix files
   required for subsequent analyses.

   There is also a ``cellranger_count_report.html`` file which links to
   the ``web_summary.html`` file for each sample, and a ZIP archive file
   which contains this index file plus the summaries, for viewing
   elsewhere.

.. _10xgenomics-additional-options:

Additional options for 'process_10xgenomics.py'
-----------------------------------------------

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


