Moving data to the archive location using ``auto_process archive``
==================================================================

Once the initial processing and QC have been completed and approved,
the data can be copied to its final location using the ``archive``
command.

.. note::

   In this context the "archive" location is where the data should
   be stored and accessed by bioinformaticians for subsequent
   analyses.

Archiving is a two stage process:

1. The default command:

   ::

       auto_process.py archive

   copies the final data (a subset of the data in the working
   directory) to a "staging" directory in the archive location.
   Multiple invocations synchronise the staging directory with
   the current state of the working directory.

2. The command:

   ::

      auto_process.py archive --final

   updates the staging area and copies the data to the final
   location; it also performs actions such as setting the
   group and permissions on the final data.

   Once ``--final`` has been used subsequent ``archive``
   commands cannot be used.

The staging directory name is the name of the analysis directory
with a double underscore prepended and with the suffix
``.pending``; for example:

::

   __180817_M00123_0001_000000000-BV1X2_analysis.pending

The final directory has the same name as the analysis directory.

By default ``archive`` inserts two additional directory levels
to the final destination, to create a ``YEAR/PLATFORM``
hierarchary. For example, if the archive location was
``/mnt/archive/`` then the full path to the staging directory
would look like

::

   /mnt/archive/2018/miseq/__180817_M00123_0001_000000000-BV1X2_analysis.pending

and the final location would be

::

   /mnt/archive/2018/miseq/180817_M00123_0001_000000000-BV1X2_analysis

---------------------
Archiving failed runs
---------------------

By default it is not possible to archive an analysis directory
which doesn't have any project directories.

However in some cases it might be desirable to archive an incomplete
analysis directory, for example if the original run had failed.

In this case the ``--force`` option of the ``archive`` command
can be used to force archiving of the analysis directory, provided
that a ``bcl2fastq`` output subdirectory (or similar) also exists.
