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
   <PROJECT>/                 Project directory (one for each    :doc:`setup_analysis_dirs <../using/setup_analysis_dirs>`
                              project defined in
                              ``projects.info``)
   undetermined/              Project directory for undetermined :doc:`setup_analysis_dirs <../using/setup_analysis_dirs>`
                              reads
   README.txt                 Text file with user notes on the   :ref:`commands_readme`
                              run (e.g. unusual processing
                              steps)
   ========================== ================================== ============

---------------------------
Analysis directory metadata
---------------------------

Each analysis has additional data items associated with it which are
stored in the ``metadata.info`` file.

The most commonly used metadata items are listed in the table below:

.. table::
   :widths: auto

   ======================= ========================================
   **Item**                **Description**
   ----------------------- ----------------------------------------
   ``run_number``          Facility-assigned identifier which
                           can differ from the instrument run
                           number
   ``source``              Source of the sequencing data, for
                           example the name of the facility,
	                   institution or service that
		           provided it
   ``platform``            The sequencing platform (e.g. ``miseq``)
   ``bcl2fastq_software``  Location and version of the package
                           used to perform the Fastq generation
   ``cellranger_software`` Location and version of the package
                           for handling 10xGenomics single cell
			   data
   ``assay``               The library prep kit (taken from the
                           sample sheet by default)
   ``sequencer_model``     The model of the sequencing instrument
   ======================= ========================================

The full set of metadata items and values for can be viewed using the
``metadata`` command:

::

    auto_process.py metadata [ANALYSIS_DIR]

This metadata is not required for processing, but should be set before
the QC is published and the analysis is completed. Items can be set or
updated using the ``--set`` option of the ``metadata`` command, for
example:

::

    auto_process.py metadata --set run_number=88

----------------
Run reference ID
----------------

Each analysis directory has a "run reference ID" which is generated
automatically from the associated metadata (specifically the platform,
run datestamp, instrument run number and facility run number).

The general form of the reference ID is:

::

    PLATFORM_DATESTAMP[/INSTRUMENT_RUN_NUMBER]#FACILITY_RUN_NUMBER

The instrument run number is included if it differs from the facility
run number (or if the facility run number is not supplied).

For example:

::

    HISEQ4000_181029/88#72

is a run from a HiSeq 4000 instrument with datestamp ``181029`` and
instrument run number ``88``; the facility assigned run ID ``72`` to
the run as its local identifier.

::

    MISEQ_180912#3

is a run a MiSeq instrument with datestamp ``180912``, where both the
instrument and facility run numbers are ``3``.

*******************
Project directories
*******************

Project directories are created within the analysis directory by
the :doc:`setup command <../using/setup>` command, based on the
contents of the ``projects.info`` file.

Each project directory will contain the following files and
directories:

   ========================== ======================================
   **File or Directory**      **Description and contents**
   -------------------------- --------------------------------------
   README.info                Project metadata
   fastqs/                    Fastq files
   ScriptCode/                Directory for custom user scripts
   qc/                        QC pipeline outputs
   qc_report.html             :doc:`QC report <qc_reports>`
   qc_report.PROJECT.RUN.zip  ZIP file containing all QC outputs
                              and reports
   multiqc_report.html        ``multiqc`` outputs
   multiqc_report_data/       Data associated with ``multiqc``
   cellranger_count/          Full outputs from ``cellranger count``
                              single library analyses
                              (10xGenomics projects only)
   ========================== ======================================

--------------------------
Project directory metadata
--------------------------

Each analysis project has additional data items associated with it
which are stored in project's ``README.info`` file.

The most commonly used metadata items are listed in the table below:

.. table::
   :widths: auto

   ======================== =========================================
   **Item**                 **Description**
   ------------------------ -----------------------------------------
   ``Run``                  Parent run name
   ``Platform``             Sequencing platform (e.g. ``miseq``)
   ``Sequencer model``      The model of the sequencing instrument
   ``User``                 Name of the user(s)
   ``PI``                   Name of PI(s)
   ``Organism``             Organism name(s)
   ``Library type``         The type of experiment (e.g. ``RNA-seq``)
   ``Library prep kit``     The kit used to prepare the libraries
                            (e.g. ``TruSeq Stranded mRNA``)
   ``Single cell platform`` Single cell platform, if applicable
   ``Number of cells``      Number of cells (single cell only)
   ``ICELL8 well list``     Well list file (ICELL8 only)
   ``Paired_end``           Whether the data are single- or paired-
                            end
   ``Primary fastqs``       Subdirectory holding the 'primary' set of
                            Fastq files for the project
   ``Samples``              Number and list of sample names
   ``Comments``             Any additional comments about the project
   ======================== =========================================

Typically most of the values are populated at setup time from the
contents of the ``projects.info`` file
(see :doc:`Setting up analysis directories <../using/setup_analysis_dirs>`),
with the others being set automatically (for example after running
single cell analyses).

-----------------------------------
Multiple Fastq sets within projects
-----------------------------------

Normally each project will only have one set of Fastq files
associated with it, and these will be in the ``fastqs``
subdirectory of the project directory.

However some analyses may have more than one sets of
associated Fastqs, and in these cases there will be multiple
subdirectories (each of which contains one of these sets).

For example, ICELL8 single cell projects typically have two
or three sets of Fastqs:

* ``fastqs.samples`` are the Fastqs after filtering and QC,
  with the reads assigned to samples according to the
  well list file
* ``fastqs.barcodes`` are the Fastqs after filtering and
  QC, with the reads assigned to barcodes according to the
  well list file
* (if present) ``fastqs`` are the original Fastq files
  produced by the BCL to Fastq conversion, without any
  additional filtering or QC

The project metadata file includes the item ``Primary fastqs``
which indicates which of the Fastq sets is the principal
one.

**********************************
``undetermined`` project directory
**********************************

This is a special project that is created for storage and QC of
the reads which couldn't be assigned to any samples by the
Fastq generation.
