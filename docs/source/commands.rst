=================
Command reference
=================

* :ref:`core-commands`
* :ref:`special-cases`
* :ref:`metadata-management`
* :ref:`global-configuration`

.. _core-commands:

*************
Core commands
*************

* :doc:`setup <using/setup>`
* :doc:`make_fastqs <using/make_fastqs>`
* :doc:`analyse_barcodes <using/analyse_barcodes>`
* :doc:`setup_analysis_dirs <using/setup_analysis_dirs>`
* :doc:`run_qc <using/run_qc>`
* :doc:`publish_qc <using/publish_qc>`
* :doc:`archive <using/archive>`

params
------

Query and update project-specific parameters and settings::

   auto_process.py params [--set PARAMTER=VALUE ...] [ANALYSIS_DIR]

.. _special-cases:

*************
Special cases
*************

import_project
--------------

Import a project directory from a separate location into the current
analysis::

   auto_process.py import_project [ANALYSIS_DIR] PROJECT_DIR

.. _metadata-management:

*******************
Metadata management
*******************

metadata
--------

Query and update metadata associated with a project::

   auto_process.py metadata [--set ITEM=VALUE ...] [ANALYSIS_DIR]

.. _global-configuration:

********************
Global configuration
********************

config
------

Query and set global configuration options that are used for all
projects::

    auto_process.py config [--set SECTION.PARAM=VALUE ...]
