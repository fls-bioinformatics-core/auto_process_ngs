Processing Takara Bio ICELL8 single-cell data
=============================================

Background
----------

The Wafergen ICell8 system prepares single cell (SC) samples which
are then sequenced as part of an Illumina sequencing run.

Initial processing of the sequencing data produces a single FASTQ file
R1/R2 pair, where the R1 reads hold an inline barcode and a unique
molecular identifier (UMI) for the read pair:

 * **Inline barcode:** bases 1-11 in R1
 * **UMI:** bases 12-25 in R1

This is illustrated in the figure below:

.. image:: images/icell8_read1.png 

.. note::

   Each inline barcode corresponds to a specific well in the
   experiment, which in turn corresponds to a specific single
   cell.

.. note::

   Older versions of the pipeline (autoprocess version 0.8.0
   and earlier) assumed that the UMI length was 10 bases.

The corresponding R2 read contains the actual sequence data
corresponding to the cell and UMI referenced by its R1 partner.

.. _icell8_processing_protocol:

Processing protocol for ICell8 data
-----------------------------------

The recommended steps are:

1. Generate initial FASTQs as described in
   :ref:`icell8_fastq_generation`
2. Set up analysis directories and run initial QC as per the standard
   protocol
3. Perform ICell8-specific filtering and additional QC on the reads
   by running the ``process_icell8.py`` utility, as described in
   :ref:`icell8_qc_and_filtering_protocol`
4. Manually update the sample name information in the ``project.info``
   and ``README.info`` files as described in
   :ref:`icell8_updating_sample_lists`

..  _icell8_qc_and_filtering_protocol:

QC and filtering protocol
-------------------------

The ``process_icell8.py`` utility script performs initial filtering
and QC according to the protocol described below. The utility also splits
the read pairs into FASTQs by both barcode and by sample.

In addition to the FASTQS produced from the  :ref:`icell8_fastq_generation`
step, the processing also requires an input :ref:`icell8_well_list_file`.
This lists the valid ICell8 barcodes and maps those barcodes to their
parent sample.

The following steps are performed:

 * **Optional initial quality screen:** this applies quality filtering
   to the barcode and UMI sequences, and rejects read pairs that fail to
   meet the following criteria:

   - Barcode bases must have Q >= 10
   - UMI bases must have Q >= 30

   By default this filtering is not used (this is recommended for
   NextSeq data). The quality filtering can be turned on by specifying
   the ``--quality-filter`` option.


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

   If the screen files aren't defined in the ``settings.ini`` file
   then they must be explicitly supplied to the utility using
   the ``-m``/``--mammalian`` and ``-c``/``--contaminants`` options.

.. note::

   The contaminant filtering step can be turned off by specifying
   the ``--no-contaminant-filter`` option, for example if analysing
   data from a non-mammalian organism.

Reorganisation by barcode and sample
------------------------------------

At the end of the QC and filter pipeline the read pairs are
reorganised in two different ways:

 * **Reorganisation by barcode:** the filtered read pairs are
   sorted into individual FASTQs according to their inline barcodes.
   This set of FASTQs forms the final outputs of the pipeline. Note
   that each barcode corresponds to a single cell, and the number of
   R1/R2 file pairs is equal to the number of barcodes/cells (~1000).

 * **Reorganisation by sample:** the read pairs are sorted into FASTQs
   according to the sample name associated with the barcodes/cells in
   the "well list" file. Essentially these group all the single cells
   from each sample, so the number of R1/R2 file pairs corresponds to
   the number of samples.

The information on valid barcodes and the relationship of barcode to
sample are taken from the :ref:`icell8_well_list_file`.

Each set of FASTQs are stored in their own directories:
``fastqs.barcodes`` and ``fastqs.samples``. Note that the read pairs
themselves are the same in each set.

The standard QC procedure is run on each set of FASTQS (barcodes and
samples) and QC reports are generated for each.

Outputs and reports
-------------------

The pipeline directory will contain the following output
directories:

 ========================== ===============================================
 **Directory**              **Description and contents**
 -------------------------- -----------------------------------------------
 ``fastqs``                 Initial FASTQs from ``bcl2fastq``
 ``fastqs.barcodes``        FASTQs with reads sorted by ICell8 barcode
                            (i.e. cell), plus QC outputs.
                            The FASTQs will be named according to the
                            convention ``NAME.BARCODE.r[1|2].fastq.gz``.
 ``fastqs.samples``         FASTQs with reads sorted by ICell8 sample
                            name (as defined in the input well list file),
                            plus QC outputs.
                            The FASTQs will be named according to the
                            convention ``SAMPLE.r[1|2].fastq.gz``.
 ``qc``                     QC for the initial FASTQs
 ``qc.barcodes``            QC for the FASTQs in ``fastqs.barcodes``
 ``qc.samples``             QC for the FASTQs in ``fastqs.samples``
 ``stats``                  Summary of the read and UMI counts after each
                            processing stage, in TSV (``icell8_stats.tsv``)
                            and XLSX format (``icell8_stats.xlsx``)
 ``logs``                   Logs from the pipeline execution
 ``scripts``                Scripts generated as part of the pipeline
                            execution.
 ``icell8_processing_data`` Data and plots for the final summary report
                            (see below)
 ========================== ===============================================

The directory will also contain:

 * A copy of the :ref:`icell8_well_list_file` (name preserved)
 * A final summary report ``icell8_processing.html``
 * A ``README.info`` file (nb only if the directory was set up as
   an autoprocess project)

The final report summarises information on the following:

 * Numbers of reads assigned to barcodes
 * Overall numbers of reads filtered after each stage
 * Initial and final read count distributions against barcodes
 * Number of reads assigned and filtered at each stage by sample
 * Poly-G region counts and distribution

.. _icell8_well_list_file:

Well list file
--------------

The well list file is a tab-delimited file output from the ICell8 which
amongst other things lists the valid ICell8 barcodes for the experiment
and the mapping of barcodes to samples.

(As noted elsewhere, each barcode corresponds to a well which in turn
corresponds to a single cell.)

.. _icell8_pipeline_configuration:

Appendix: configuring the ICell8 processing pipeline
----------------------------------------------------

The running of the pipeline can be configured via command line options,
or by setting the appropriate parameters options in the ``settings.ini``
configuration file.

Reference data and quality filtering
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 * **Mammalian genome panel**: ``fastq_screen`` conf file with the
   indices for "mammalian" genomes, to use in the contamination
   filtering step.

   Set using the ``-m`` option on the command line, or via
   ``[icell8] mammalian_conf_file`` in the configuration file.

 * **Contaminant genome panel**: ``fastq_screen`` conf file with the
   indices for "contaminant" genomes, to use in the contamination
   filtering step.

   Set using the ``-c`` option on the command line, or via
   ``[icell8] contaminant_conf_file`` in the configuration file.

   To turn off the contaminant filtering, specify the
   ``--no-contaminant-filter`` option.

 * **Quality filtering of barcode and UMI sequences**: by default
   read pairs are *not* removed if the associated barcode or UMI
   sequences don't meet the appropriate quality criteria.

   To turn on quality filtering, specify the
   ``-q``/``--quality_filter`` option (nb there is no equivalent
   parameter in the configuration file).

Runtime environment
~~~~~~~~~~~~~~~~~~~

 * **Environment modules**: specify a list of environment modules
   (separated with commas) to load before running the pipeline.

   Set using the ``--modulefiles`` option on the command line, or
   via ``[modulefiles] process_icell8`` in the configuration file.

 * **Job runners and processors**: specify job runners and number
   of processors to use for specific classes of tasks in the pipeline.
   See :ref:`job_runners_and_processors` for more details.

 * **Aligner**: explicitly specify the aligner (currently either
   ``bowtie`` or ``bowtie2``) to use for contamination filtering.

   Set using the ``-a`` option on the command line, or via
   ``[icell8] aligner`` in the configuration file. (NB if this is
   not set then an appropriate aligner will be selected
   automatically from those available in the execution
   environment.)

FASTQ batching
~~~~~~~~~~~~~~

 * **Read batch size**: number of reads to assign to each "batch"
   when splitting FASTQs for processing.

   Batching the reads enables many of the pipeline tasks to run
   in parallel, if the execution environment allows it (e.g. if
   running on a compute cluster).

   Set using the ``-s`` option on the command, or via
   ``[icell8] batch_size``.

Job control
~~~~~~~~~~~

 * **Maximum number of concurrent jobs**: limits the number of
   processes that the pipeline will attempt to run at any one
   time.

   The default is taken from the ``max_concurrent_jobs``
   parameter in the configuration file; it can be set at run
   time using the ``-j``/``--max-jobs`` command line option.

..  _job_runners_and_processors:

Job runners and processors
~~~~~~~~~~~~~~~~~~~~~~~~~~

Job runners and numbers of processors can be explicitly defined
for different "stages" of the pipeline, where a stage is
essentially a class of tasks).

For the ICell8 processing pipeline the stages are:

 ================== ========================================
 **Name**           **Description**
 ------------------ ----------------------------------------
 contaminant_filter Tasks for filtering "contaminated" reads
 qc                 Tasks for performing QC on the FASTQs
 statistics         Tasks for generating statistics
 ================== ========================================

Use the ``-n``/``--nprocessors`` and ``-r``/``--runners`` options
to specify the number of cores that can be used, and an appropriate
runner (if necessary) for each of these stages.

Via the command line e.g.::

    process_icell.py ... -r statistics='GEJobRunner(-pe smp.pe 4)' -n 4

Via the configuration file::

    [icell8]
    nprocessors_statistics = 4

    [runners]
    icell8_statistics = GEJobRunner(-pe smp.pe 4)

.. _icell8_updating_sample_lists:

Appendix: manually updating sample lists
----------------------------------------

Currently the processing pipeline implemented in ``process_icell8.py``
doesn't automatically update the sample lists in ``projects.info``
and the ``README.info`` in the ICell8 project directories.

To do this manually requires extracting the sample names and editing
the files to update them with the correct data.

The sample names can be extracted from the well list file using the
command::

    tail -n +2 WellList.TXT | cut -f5 | sed 's/ /_/g' | sort -u | paste -s -d","

which produces a list suitable for ``projects.info`` (e.g.
``Pos_Ctrl,SC004,SC005,SC006``).

To get them in a format suitable for the ``README.info`` file::

    tail -n +2 WellList.TXT | cut -f5 | sed 's/ /_/g' | sort -u | paste -s -d"," | sed 's/,/, /g'

(e.g. ``Pos_Ctrl, SC004, SC005, SC006, SC007``).

The number of samples can be obtained by::

    tail -n +2 WellList.TXT | cut -f5 | sort -u | wc -l
