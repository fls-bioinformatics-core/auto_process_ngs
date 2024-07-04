Processing 10x Genomics Visium spatial transcriptomics data
===========================================================

Background
----------

10x Genomics provides the Visium platform for spatial transcriptomics.
This document outlines using the ``auto-process-ngs`` pipeline to
perform Fastq generation and QC for various types of Visium data.

Requirements
------------

10x Genomics' SpaceRanger software is used for Fastq generation for
Visium data. It should be downloaded and installed before running the
pipeline. No additional external software is required for the QC.

Fastq generation
----------------

If a sample sheet with the appropriate 10x Genomics indexes is provided
then all Visium data should be processed using the ``10x_visium`` protocol
with the :doc:`make_fastqs <../using/make_fastqs>` command, for example:

::

   auto_process.py make_fastqs --protocol=10x_visium

This protocol wraps the ``spaceranger mkfastq`` command.

.. note::

   If the sample sheet contains Illumina index sequences then the
   ``standard`` protocol should be used instead (note that in this case
   the defaults used for masking and trimming compared to the defaults
   may differ from those used by SpaceRanger).

Analysis project setup and QC
-----------------------------

Once Fastqs have been successfully generated, the ``SC_platform``
and ``Library`` metadata items should be set to the appropriate values
for the Visium project(s) in the ``projects.info`` control file.

The following values are valid options for spatial data:

===================================== ==============================
Single cell platform                  Library types
===================================== ==============================
``10xGenomics Visium``                ``Spatial RNA-seq``,
                                      ``FFPE Spatial RNA-seq``,
                                      ``Fresh Frozen RNA-seq``
``10xGenomics CytAssist Visium``      ``FFPE Spatial RNA-seq``
                                      ``FFPE Spatial GEX``,
                                      ``FFPE Spatial PEX``,
                                      ``HD Spatial GEX``
===================================== ==============================

.. note::

   Although 10xGenomics Visium data are not single cell data,
   for historical reasons the platform information is currently
   being stored in the single cell platform metadata field.

Running the :doc:`setup_analysis_dirs <../using/setup_analysis_dirs>`
command will automatically transfer these values into the Visium
project metadata on creation.

The :doc:`run_qc <../using/run_qc>` command
will then determine the appropriate QC protocol to use based on the
metadata values.

Image data
----------

The ``setup_analysis_dirs`` command also creates an empty
``Visium_images`` subdirectory within the project, which should
either be populated manually with the corresponding image files,
or else be deleted.

.. note::

   The presence of an empty images directory will cause the final
   archiving to stop with an error (as a reminder to copy in the
   associated images).
