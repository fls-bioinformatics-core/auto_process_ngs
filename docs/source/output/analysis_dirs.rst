Analysis and Project Directories
================================

********************
Analysis directories
********************

The top-level analysis directory is created by the
:doc:`setup command <../using/setup>`, and typically will have
a name based on the sequencing run name with the suffix
``_analysis``, for example:

::

   180817_M00123_0001_000000000-BV1X2_analysis

The analysis directory will contain the following files and
directories produced by the ``auto_process`` commands:

.. table::
   :widths: auto

   ========================== ================================== ============
   **File or Directory**      **Description and contents**       **Stage**
   -------------------------- ---------------------------------- ------------
   auto_process.info          Parameter file used for processing :doc:`setup <../using/setup>`
   metadata.info              Metadata for the run               :doc:`setup <../using/setup>`
   logs/                      Log files from each processing     :doc:`setup <../using/setup>`
                              command
   ScriptCode/                Directory for custom user scripts  :doc:`setup <../using/setup>`
   SampleSheet.orig.csv       Copy of the original sample sheet  :doc:`setup <../using/setup>`
                              file
   custom_SampleSheet.csv     Updated version of sample sheet    :doc:`setup <../using/setup>`
                              used for processing
   primary_data/              Raw sequencing data                :doc:`make_fastqs <../using/make_fastqs>`
   processing_qc.html         Processing QC report               :doc:`make_fastqs <../using/make_fastqs>`
   per_lane_stats.info        Per-lane statistics                :doc:`make_fastqs <../using/make_fastqs>`
   per_lane_sample_stats.info Per-lane statistics for samples    :doc:`make_fastqs <../using/make_fastqs>`
   statistics_full.info       Per-Fastq statistics               :doc:`make_fastqs <../using/make_fastqs>`
   statistics.info            Per-Fastq statistics               :doc:`make_fastqs <../using/make_fastqs>`
   projects.info              Metadata for all projects          :doc:`make_fastqs <../using/make_fastqs>`
   <bcl2fastq>/               Output from ``bcl2fastq``          :doc:`make_fastqs <../using/make_fastqs>`
                              (can be set explicitly using the
                              ``--output-dir`` option)
   <PROJECT>/                 Project directory (one for each
                              project defined in
                              ``projects.info``)
   undetermined/              Project directory for undetermined
                              reads
   README.txt                 Text file with user notes on the   readme
                              run (e.g. unusual processing
                              steps)
   ========================== ================================== ============

---------------------------
Analysis directory metadata
---------------------------

The following additional information is associated with each analysis:

.. table::
   :widths: auto

   ====================== ========================================
   **Item**               **Description**
   ---------------------- ----------------------------------------
   ``run_number``         Locally-assigned identifier which
                          can differ from the instrument run
   ``source``             Source of the sequencing data, for
                          example the name of the facility,
	                  institution or service that
		          provided it
   ``platform``           The sequencing platform (e.g. ``miseq``)
   ``bcl2fastq_software`` Location and version of the package
                          used to perform the Fastq generation
   ====================== ========================================

This metadata is not required for processing, but should be set before
the QC is published and the analysis is completed.

The metadata for an analysis directory can be inspected using the
``metadata`` command::

    auto_process.py metadata

and updated using the ``--set`` option::

    auto_process.py metadata --set run_number=88

*******************
Project directories
*******************

Project directories are created within the analysis directory by
the :doc:`setup command <../using/setup>` command, based on the
contents of the ``projects.info`` file.

Each project directory will contain the following files and
directories:

   ========================== ==================================
   **Directory**              **Description and contents**
   -------------------------- ----------------------------------
   README.info                Project metadata
   fastqs/                    Fastq files
   ScriptCode/                Directory for custom user scripts
   qc/                        QC pipeline outputs
   qc_report.html             :doc:`QC report <qc_reports>`
   qc_report.PROJECT.RUN.zip  ZIP file containing all QC outputs
   multiqc_report.html        ``multiqc`` outputs
   multiqc_report_data/       Data associated with ``multiqc``
   ========================== ==================================

**********************************
``undetermined`` project directory
**********************************

This is a special project that is created for storage and QC of
the reads which couldn't be assigned to any samples by the
Fastq generation.
