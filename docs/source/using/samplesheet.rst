Sample sheet manipulation using ``auto_process samplesheet``
============================================================

The ``samplesheet`` command provides functionality for interrogating
and manipulating the Illumina sample sheet control file.

When the command is invoked without arguments then the contents
of the sample sheet are displayed.

The ``-p``/``--predict`` argument shows the predicted outputs from
the sample sheet.

The following arguments enable key data in the sample sheet to be
updated in batch:

* ``--set-project``: set the sample project for some or all lines
  in the file;
* ``--set-sample-id``: set the sample ID for some or all lines
  in the file (typically to synchronise the IDs with the sample
  names);
* ``--set-sample-name``: set the sample project for some or all
  lines in the file (typically to synchronise the names with the
  sample IDs).

The same syntax is used for the ``--set-*`` functions to specify
the new values and to optionally select subsets of lines to
apply the updates to:

::

   [LANES:][COL=PATTERN:]NEW_VALUE

where ``LANES`` specifies one or more lane numbers (for example:
``1``, ``1,2,3``, ``1-3``, ``1,3-5`` etc), and ``COL`` specifies the
name of a column in the sample sheet where lines will be selected
only if the value in that column matches the glob-style
``PATTERN`` (for example: ``Sample_Name=ITS*``).

The ``-e``/``--edit`` argument brings up the sample sheet in an
editor to allow changes to be made manually.
