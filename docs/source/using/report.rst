Reporting analysis directory contents using ``auto_process report``
===================================================================

Generates reports about the contents of an analysis directory.

General invocation of the command is:

::

   auto_process.py report *reporting_option* [ANALYSIS_DIR]

where *reporting_option* determines the data that are reported and
the format that is used.

The following options are available:

=================== =====================================
Reporting option    Description
=================== =====================================
``--logging``       Single line summary of all projects
                    in the analysis directory
``--projects``      One tab-delimited line per project
``--summary``       Longer format report suitable for
                    bioinformaticians
=================== =====================================

These options are described in more detail below.

``--logging`` report
--------------------

Produces a single line summary of the projects in the analysis
directory.

For example for a run with a single project:

::

    Paired end: 'AB': Abby Brown, Mouse RNA-seq (PI: Carl Dover) (4 samples)

For runs with multiple projects, the details of the additional
projects are appended with semi-colons as separators.

``--projects`` report
---------------------

Outputs a tab-delimited line for each project in the analysis
directory, which could be inserted into a spreadsheet.

For example for a run with a single project:

::

    MISEQ_150729#88  88   GTCF       Abby Brown   Carl Dover  RNA-seq   Mouse   MISEQ   4    yes     AB1-4


``--summary`` report
--------------------

Outputs a longer format report which summarises all the data and
projects in the analysis directory.

For example:

::

    MISEQ run #88 datestamped 150729
    ================================
    Run name : 150729_M00789_0088_000000000-ABCD1
    Reference: MISEQ_150729#88
    Platform : MISEQ
    Directory: /runs/2015/miseq/150729_M00789_0088_000000000-ABCD1_analysis
    Endedness: Paired end
    Bcl2fastq: bcl2fastq 2.17.1.14
    Assay    : TruSeq HT

    1 project:
    - 'AB': Abby Brown Mouse RNA-seq 4 samples (PI Carl Dover)

    Additional notes/comments:
    - AB: 1% PhiX spike in

Customing data that are reported
--------------------------------

The data that are reported can be customised by using the ``--fields``
option of the ``report`` command.

The available fields are:

========================= ========================
Field name                Associated value
========================= ========================
``run_id``                Run ID (e.g. ``MISEQ_150729#88``)
``run_number``            Facility run number (e.g. ``88``)
``sequencer_platform``    Sequencer platform for the run
                          (e.g. ``MISEQ``)
``platform``              Alias for ``sequencer_platform``
``source``                Source of the sequencing data
``data_source``           Alias for ``source``
``user``                  User associated with project
``PI``                    Principle investigator associated
                          with the project
``pi``                    Alias for ``PI``
``application``           Application associated with the
                          project (e.g. ``RNA-seq``)
``library_type``          Alias for ``application``
``single_cell_platform``  Name of the single-cell platform
                          (e.g. ``ICELL8``), or empty field
			  if not a single-cell project
``organism``              Organism(s) associated with the
                          project
``no_of_samples``         Number of samples in the project
``#samples``              Alias for ``no_of_samples``
``no_of_cells``           Number of cells in the project,
                          for single-cell projects
``#cells``                Alias for ``no_of_cells``
``paired_end``            ``yes`` if the run was paired-end,
                          ``no`` if it was single-end
``sample_names``          Comma-separated list of sample
                          names associated with the project
``samples``               Alias for ``sample_names``
``null``                  Writes an empty field
========================= ========================

.. note::

   Custom fields are only available for the ``--projects``
   reporting mode.
