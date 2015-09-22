Metadata
========

The following additional information is associated with each analysis:

* ``run_number``: a locally-assigned identifier (which can differ
  from the run number encoded in the run name)
* ``source``: where the data originated from, for example the
  name of the facility, institution or service that provided it
* ``platform``: indicates the sequencing platform (e.g. ``miseq``)
* ``bcl2fastq_software``: information on the location and version of
  the package that was used to perform the BCL-to-FASTQ conversion.

This metadata is not required for processing, but should be set before
the QC is published and the analysis is completed.

The metadata for an analysis directory can be inspected using the
``metadata`` command::

    auto_process.py metadata

and updated using the ``--set`` option::

    auto_process.py metadata --set run_number=88

