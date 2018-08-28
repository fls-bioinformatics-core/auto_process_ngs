Making QC reports available using ``auto_process publish_qc``
=============================================================

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

By default the QC for all the projects in the analysis
directory will be copied; to only publish the QC for a subset,
specify their names using the ``--projects`` option.

The QC reports will only be published if the QC is verified
for all the specified projects; to force copying even if
some projects have missing QC, use the ``--ignore-missing-qc``
option.

By default the QC will be copied to a directory with the
same name as the analysis directory, directly underneath
the target directory - for example:

::

   /mnt/qc_reports/180817_M00123_0001_000000000-BV1X2_analysis

However this can quickly lead to an extremely cluttered
QC server directory if reports are not removed periodically.

To mitigate this, specifying ``--use-hierarchy=yes`` inserts
two additional directory levels to the final destination,
to create a ``YEAR/PLATFORM`` hierarchary. For the previous
example this would result in something like:

::

   /mnt/qc_reports/2018/miseq/180817_M00123_0001_000000000-BV1X2_analysis

Repeated runs of ``publish_qc`` will update the copies on the
server.
