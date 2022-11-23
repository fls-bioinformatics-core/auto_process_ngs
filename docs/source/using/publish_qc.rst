Making QC reports available using ``auto_process publish_qc``
=============================================================

--------
Overview
--------

Once QC reports have been generated from :doc:`run_qc <run_qc>`
they can be copied to a webserver or other location using the
``publish_qc`` command.

If a QC server has been specified in the configuration file
then the reports can be published there using

::

   auto_process.py publish_qc

If no server was configured (or to copy to an alternative
destination) use the ``--qc_dir`` option, for example:

::

   auto_process.py publish_qc --qc_dir=/mnt/qc_reports

.. note::

   The destination specified by ``--qc_dir`` can be either
   a local or a remote directory; to specify a remote
   destination, use the ``[USER@]HOST:DIR`` syntax.

If a base URL is also specified (either in the configuration,
or via the ``--url`` option) then the URL of the published
HTML index will also be reported.

Repeated runs of ``publish_qc`` will update the copies on the
server.

---------------------------------
Using a hierarchy for publication
---------------------------------

By default the QC will be copied to a directory with the
same name as the analysis directory, directly underneath
the target directory - for example:

::

   /mnt/qc_reports/180817_M00123_0001_000000000-BV1X2_analysis

However this can quickly lead to an extremely cluttered
QC server directory if reports are not removed periodically.

To mitigate this, specifying ``--use-hierarchy=yes`` inserts
two additional directory levels to the final destination,
to create a ``YEAR/PLATFORM`` hierarchy. For the previous
example this would result in something like:

::

   /mnt/qc_reports/2018/miseq/180817_M00123_0001_000000000-BV1X2_analysis

---------------------------------------------
Selecting subsets of projects for publication
---------------------------------------------

By default the QC for all the projects in the analysis
directory will be copied; to only publish the QC for a subset,
specify their names using the ``--projects`` option.

--------------------------------
Handling projects with failed QC
--------------------------------

Reports will only be published if the QC is verified for all the
specified projects associated with the run; otherwise the publish
command will stop.

However there are a number of options available to handle this
situation:

1. Use the ``--ignore-missing-qc`` to only publish projects which
   pass the verification, skipping those with missing QC, or

2. Specify the ``--force`` command to ignore the failed QC and
   force generation and publication of the missing QC.

For #2, the reports for the projects with failed QC will contain a
warning message in the header, and the links to the reports will be
marked with a warning icon in the index page. The
``--suppress-warnings`` options can be used along with ``--force``
to suppress these warnings if necessary.

.. note::

   ``--suppress-warnings`` will only suppress warnings in new
   QC reports; to make this option apply to all reports, combine
   it with the ``--regenerate-reports`` option.

----------------------------------
Updating QC reports on publication
----------------------------------

Sometimes it is useful to regenerate the QC reports for all projects,
in which case the ``--regenerate-reports`` option can be specified
to force all the reports to be regenerated prior to publication.
