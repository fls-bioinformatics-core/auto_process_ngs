Autoprocessing commands
***********************

* :ref:`core-commands`
* :ref:`metadata-management`
* :ref:`global-configuration`

.. _core-commands:

Core commands
=============

setup
-----

Creates a new top-level analysis directory for processing data from
a sequencing run::

   auto_process.py setup DATA_DIR [ANALYSIS_DIR]

params
------

Query and update project-specific parameters and settings::

   auto_process.py params [--set PARAMTER=VALUE ...] [ANALYSIS_DIR]

make_fastqs
-----------

Performs Fastq generation from bcl files::

   auto_process.py make_fastqs [ANALYSIS_DIR]

setup_analysis_dirs
-------------------

Creates subdirectories populated with Fastq files for each project::

   auto_process.py setup_analyis_dirs [ANALYSIS_DIR]

run_qc
------

Run the QC scripts on the Fastqs files for each project::

   auto_process.py run_qc [ANALYSIS_DIR]

publish_qc
----------

Copy reports from the QC runs to a webserver or other location::

   auto_process.py publish_qc [ANALYSIS_DIR]

archive
-------

Copy the final data to an 'archive' location::

   auto_process.py archive [ANALYSIS_DIR]

.. _metadata-management:

Metadata management
===================

metadata
--------

Query and update metadata associated with a project::

   auto_process.py metadata [--set ITEM=VALUE ...] [ANALYSIS_DIR]

.. _global-configuration:

Global configuration
====================

config
------

Query and set global configuration options that are used for all
projects::

    auto_process.py config [--set SECTION.PARAM=VALUE ...]
