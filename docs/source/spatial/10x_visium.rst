Processing 10x Genomics Visium spatial transcriptomics data
===========================================================

Background
----------

10x Genomics provides the Visium platform for spatial transcriptomics.
This document outlines using the ``auto-process-ngs`` pipeline to
perform Fastq generation and QC for various types of Visium data.

Requirements
------------

10x Genomics' SpaceRanger software is required for Fastq generation for
Visium data for sample sheets where 10x Genomics index codes have been
specified; in this case, it should be downloaded and installed before
running the pipeline.

If the sample sheet uses standard Illumina index sequences then
SpaceRanger is not needed.

SpaceRanger is not required for the QC.

Fastq generation
----------------

The following Fastq generation protocols are available for Visium data:

.. include:: ../auto/10x_spatial_fq_protocols.rst

The appropriate protocol should be selected from the table above according
to the spatial platform, tissue preparation method and library type, and
supplied to the ``make_fastqs`` command - for example:

::

   auto_process.py make_fastqs --protocol=10x_visium

Depending on the format of the index sequences in the input sample sheet,
this will either use the ``spaceranger mkfastq`` command (if the indexes
are 10x Genomics index codes) or ``bcl2fastq`` or ``bcl-convert`` directly
(if the indexes are standard Illumina sequences).

Analysis project setup and QC
-----------------------------

Once Fastqs have been successfully generated, the ``SC_platform``
and ``Library`` metadata items should be set to the appropriate values
for the Visium project(s) in the ``projects.info`` control file.

The following values are valid options for spatial data:

.. include:: ../auto/10x_visium_apps.rst

.. note::

   Although 10xGenomics Visium data are not single cell data,
   for historical reasons the platform information is currently
   being stored in the single cell platform metadata field.

.. note::

   The platform/library combinations above are based on the
   information at https://www.10xgenomics.com/platforms/visium

   Where ``GEX`` and ``PEX`` appear in the library types, the
   long form versions ``Gene Expression`` and ``Protein Expression``
   can also be used.

   Some legacy values are also still recognised, for example
   ``10xGenomics CytAssist Visium`` (as an alternative to
   ``10xGenomics Visium (CytAssist)``), and ``spatial RNA-seq``
   (as an alternative to ``spatial Gene Expression`` or
   ``spatial GEX``).

Running the :doc:`setup_analysis_dirs <../using/setup_analysis_dirs>`
command will automatically transfer these values into the Visium
project metadata on creation.

The :doc:`run_qc <../using/run_qc>` command
will then determine the appropriate QC protocol to use based on the
metadata values.

Currently the following protocols are defined:

============================== ==================================
Protocol name                  Used for
============================== ==================================
``10x_Visium_GEX``             Spatial Gene Expression and HD
                               Spatial Gene Expression data with
                               50bp insert in R2
``10x_Visium_PEX``             Spatial Gene Expression data with
                               50bp insert in R2
``10x_Visium_GEX_90bp_insert`` Spatial Gene Expression data with
                               90bp insert in R2 (fresh frozen
                               tissues only)
``10x_Visium_legacy``          Spatial Gene Expression data where
                               the R2 insert size is unknown (so
                               is not accounted for) (NB this
                               protocol is deprecated and is only
                               maintained for backwards
                               compatibility)
============================== ==================================

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
