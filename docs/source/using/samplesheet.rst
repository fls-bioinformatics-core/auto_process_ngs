Sample sheet manipulation using ``auto_process samplesheet``
============================================================

The ``samplesheet`` command provides functionality for interrogating
and manipulating the Illumina sample sheet control file.

--------------------
Default sample sheet
--------------------

The default sample sheet file is set in the parameters of the
current analysis directory; all operations will be performed
using this sample sheet.

When the ``samplesheet`` command is invoked without arguments
then the path and contents of the default sample sheet are
displayed.

The default sample sheet file can be updated with the ``--use``
option, for example:

::

   auto_process.py samplesheet --use new_samplesheet.csv


------------------
Predicting outputs
------------------

The ``-p``/``--predict`` argument shows the predicted outputs from
the sample sheet.

------------------------------------------
Updating project and sample names in batch
------------------------------------------

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

---------------------------------
Manually editing the sample sheet
---------------------------------

The ``-e``/``--edit`` argument brings up the sample sheet in an
editor to allow changes to be made manually.

---------------------------------
Importing a new sample sheet file
---------------------------------

The ``-i``/``--import`` argument enables content from a new file to
be imported into the analysis directory and used as the new default
sample sheet.

The imported file can be a local or remote file, or a URL, for
example:

::

   auto_process.py samplesheet -i user@remote.com:/Data/new_samplesheet.csv

.. note::

   Importing differs from the ``--use`` argument in that a copy
   of the imported file is made within the analysis directory
   (whereas ``--use`` simply changes the location that the default
   sample sheet parameter points to).
