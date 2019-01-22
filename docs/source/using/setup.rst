Starting an analysis using ``auto_process setup``
=================================================

The ``setup`` command is used to start a new analysis: it
initialises a new analysis directory for processing a
sequencing run. Most subsequent ``auto_process`` commands
(for example ``make_fastqs``) are normally issued from within
the analysis directory.

The simplest invocation of the command is:

::

   auto_process.py setup DATA_DIR

where ``DATA_DIR`` specifies the location of the data source
(i.e. the top-level directory containing the output from the
sequencer run which is to be processed). For example:

::

   auto_process.py setup /mnt/data/seqruns/180817_M00123_0001_000000000-BV1X2

``DATA_DIR`` can be a local or a remote directory; see
:ref:`setup_remote_data_dir`.

By default the name of the new analysis directory will consist
of the basename of ``DATA_DIR`` with the suffix ``_analysis``
appended, for example the command above will produce an analysis
directory called

::

   180817_M00123_0001_000000000-BV1X2_analysis

To create an analysis directory with a different name, specify
it using the ``--analysis-dir`` option:

::

   auto_process.py setup DATA_DIR --analysis-dir ANALYSIS_DIR

See :doc:`Analysis Directories <../output/analysis_dirs>` for
details of the analysis directory structure.

The ``setup`` command will also report the expected outputs
based on the sample sheet associated with the sequencing run,
for example:

::

   Predicted projects:
   ===================
   - JohnBleakley
   - LauraBridges
   - MarcusDreng
   - StevenYound

   JohnBleakley (24 samples)
   -------------------------
   JB1	S11	TGCGGCGT-TACCGAGG	L3,4
   JB2	S12	CATAATAC-CGTTAGAA	L3,4
   JB3	S13	GATCTATC-AGCCTCAT	L3,4
   JB4	S14	AGCTCGCT-GATTCTGC	L3,4
   JB5	S15	CGGAACTG-TCGTAGTG	L3,4
   JB6	S16	TAAGGTCA-CTACGACA	L3,4
   ...

It will also flag up any potential issues (for example if
two project names are very similar then this might indicate
a typo in a project name in the sample sheet).

.. note::

   The output prediction can also be generated using
   ``auto_process samplesheet`` command.

.. _setup_specifying_sample_sheet:

********************************
Specifying the sample sheet file
********************************

The ``setup`` command will try to locate the sample sheet
within the source data and will use this by default.

However if the sequencing run doesn't include a sample
sheet file (for example, NextSeq runs), or if you want to
use an alternative sample sheet, then the ``--sample-sheet``
option can be used to explicitly specify the sample sheet
file.

For example:

::

   auto_process.py setup \
      --sample-sheet /mnt/data/samplesheets/SampleSheet_180817.csv \
      /mnt/data/seqruns/180817_M00123_0001_000000000-BV1X2

The sample sheet can also be on a remote system, for example:

::

   auto_process.py setup \
      --sample-sheet pjb@kellerman.man.ac.uk:samplesheets/SampleSheet_180817.csv \
      /mnt/data/seqruns/180817_M00123_0001_000000000-BV1X2

or it can be a URL:

::

   auto_process.py setup \
      --sample-sheet https://example.com/samplesheets/SampleSheet_180817.csv \
      /mnt/data/seqruns/180817_M00123_0001_000000000-BV1X2

.. _setup_remote_data_dir:

********************************************
Specifying a remote sequencing run directory
********************************************

For data on a remote system which is accessible via ``ssh``,
the ``DATA_DIR`` can be specified using the general syntax

::

   [[USER@]HOST:]DATA_DIR

For example:

::

   pjb@kellerman.man.ac.uk:/mnt/data/seqruns/180817_M00123_0001_000000000-BV1X2

(The ``--sample-sheet`` option accepts the same syntax.)

.. note::

   It is recommended that either passwordless ``ssh`` access
   is configured, or that ``ssh-agent`` is used for the
   current session, to suppress multiple password prompts
   each time the remote system is accessed.

.. _setup_import_fastqs:

************************************
Setup from existing bcl2fastq output
************************************

A new analysis directory can be created from an existing
``bcl2fastq`` output directory using the ``--fastq-dir``
option, which should be used to specify the subdirectory
of the ``DATA_DIR`` which contains the output Fastq files.

For example:

::

   auto_process.py setup \
      --fastq-dir bcl2fastq2 \
      /mnt/data/seqruns/180817_M00123_0001_000000000-BV1X2

where ``bcl2fastq2`` is the output directory from the
BCL-to-Fastq conversion software, within the run data
directory ``180817_M00123_0001_000000000-BV1X2``.
