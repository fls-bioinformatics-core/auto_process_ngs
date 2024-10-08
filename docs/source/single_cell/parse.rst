Processing Parse Evercode single cell data
==========================================

Background
----------

Parse Biosciences provides the Evercode platform for single cell
sequencing. This document outlines using the ``auto-process-ngs``
pipeline to perform Fastq generation and QC for Parse Evercode data.

Requirements
------------

No additional external software is required for the Fastq generation
or QC.

Fastq generation
----------------

The ``parse_evercode`` Fastq generation protocol should be used
for Parse Evercode samples when running the
:doc:`make_fastqs <../using/make_fastqs>` command.


Analysis project setup and QC
-----------------------------

Once Fastqs have been successfully generated, the ``SC_platform``
and ``Library`` metadata items should be set to the appropriate values
for the Parse Evercode project(s) in the ``projects.info`` control file.

The following values are valid options:

===================================== =================================
Single cell platform                  Library types
===================================== =================================
``Parse Evercode``                    ``scRNA-seq``, ``snRNA-seq``,
                                      ``TCR``, ``TCR scRNA-seq``,
                                      ``WT``, ``WT scRNA-seq``
===================================== =================================

Running the :doc:`setup_analysis_dirs <../using/setup_analysis_dirs>`
command will automatically transfer these values into the project
metadata on creation. The :doc:`run_qc <../using/run_qc>` command
will then determine the appropriate QC protocol to use based on the
metadata values.
