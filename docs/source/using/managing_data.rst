Managing and sharing data
=========================

These additional tasks are not really part of the automated processing, but
the utilities for performing them are currently part of the same package so
they are outlined here.

 * :ref:`manage_fastqs`
 * :ref:`transfer_data`
 * :ref:`download_fastqs`
 * :ref:`update_project_metadata`
 * :ref:`audit_projects`

There are no specific utilities for exporting Fastqs to a Galaxy data
library, but suggestions on how to do this can be found in the section
:ref:`exporting_to_galaxy`.

.. _manage_fastqs:

``manage_fastqs.py``: managing and copy Fastq files
***************************************************

The ``manage_fastqs.py`` utility can be used to explore and copy data from
an analysis directory to another location on a per-project basis.

For example, to get a list of projects within a run::

    manage_fastqs.py ANALYSIS_DIR

.. note::

   If a project contains multiple sets of FASTQs then these
   will be listed under the project name, with the default
   "primary" project marked with an asterisk.

To see a list of the FASTQ files associated with a particular project::

    manage_fastqs.py ANALYSIS_DIR PROJECTNAME

To copy the FASTQs to a local or remote directory, use the ``copy`` command::

    manage_fastqs.py ANALYSIS_DIR PROJECTNAME copy /path/to/local/dir
    manage_fastqs.py ANALYSIS_DIR PROJECTNAME copy me@remote.org:/path/to/remote/dir

This will also generate an MD5 checksum file for the transferred files; to
generate the MD5 sums on their own, use the ``md5`` command::

    manage_fastqs.py ANALYSIS_DIR PROJECTNAME md5

Finally to make a zip file containing the Fastqs, use the ``zip`` command::

    manage_fastqs.py ANALYSIS_DIR PROJECTNAME zip

.. note::

    The ``zip`` option works best if the Fastqs are relatively small.

Working with multiple Fastq sets
--------------------------------

If a project has more than one Fastq set associated with it then by
default the operations described above will use the "primary" set
(typically, the set of Fastq files in the ``fastqs`` subdirectory
of the project).

To operate on an alternative set, use the ``--fastq_dir`` option to
switch e.g.::

    manage_fastqs.py ANALYSIS_DIR PROJECTNAME --fastq_dir=ALT_FASTQS_DIR

Handling subsets of files
-------------------------

Use the ``--filter`` option to work with a subset of files - this allows a
'glob'-style pattern to be specified so that only files with matching names
will be included.

For example to only copy ``R1`` files::

    manage_fastqs.py ANALYSIS_DIR PROJECTNAME copy /path/to/local/dir --filter *_R1_*

.. _transfer_data:

``transfer_data.py``: copying data for transfer to end users
************************************************************

The ``transfer_data.py`` utility can be used to copy data from analysis
projects to different destinations, typically to transfer copies of
data to end users.

A destination is defined as a local or remote directory where files
will be copied, for example in its most basic mode:

::

    transfer_data.py /mnt/data/shared ANALYSIS_DIR PROJECT

will copy Fastq files from the ``PROJECT`` project in ``ANALYSIS_DIR``
to the local directory ``/mnt/data/shared``.

Destinations can also be defined in the configuration file (see
:ref:`data_transfer_destinations`) and then referred to by their
name when copying the Fastqs.

For example:

::

    transfer_data.py webserver ANALYSIS_DIR PROJECT

where ``webserver`` is a pre-defined destination.

Schemes for dymanic subdirectory specification
----------------------------------------------

By default the data are copied directly to the specified directory.
However it is possible to specify a scheme for dynamic subdirectory
assignment, which can be useful for example if copying to a
webserver.

The scheme can be specified via either the ``--subdir`` command line
option or the ``subdir`` parameter in the configuration file.

The following schemes are available:

==============  ==========================================
Scheme name     Behaviour
==============  ==========================================
``random_bin``  Locates an empty pre-existing subdirectory
                (aka 'bin') at random
``run_id``      Creates a new subdirectory named
                ``PLATFORM_DATESTAMP.RUN_NUMBER-PROJECT``
                (must not already exist)
==============  ==========================================

Generating a README file from a template
----------------------------------------

It is possible to generate a ``README`` for the copied data by
specifying a template file via either the ``--readme`` command line
option or the ``readme_template`` parameter in the configuration
file.

The template should be a plain text file but it can also contain
placeholders for 'template variables' which will be substituted with
the appropriate values when the ``README`` file is generated:

================  =================================
Placeholder       Value
================  =================================
``%PLATFORM%``    Run platform (uppercase)
``%RUN_NUMBER%``  Run number
``%DATESTAMP%``   Run datestamp
``%PROJECT%``     Name of project being copied
``%WEBURL%``      Base URL for the webserver
``%BIN%``         Name of the subdirectory, if any
``%DIR%``         Directory data were copied to
``%TODAY%``       Today's date
================  =================================

Including downloader and QC reports
-----------------------------------

By default only Fastqs are copied by ``transfer_data.py``, however it
is possible to include additional files:

 * A standalone downloader script (see :ref:`download_fastqs`)
   (specify the ``--include_downloader`` or set the
   ``include_downloader`` parameter in the configuration);
 * The zipped QC reports for the project (specify the
   ``--include_qc_report`` or set the ``include_qc_report``
   parameter)

.. _download_fastqs:

``download_fastqs.py``: fetch Fastqs from a webserver in batch
**************************************************************

Fastq files pushed to a webserver using ``manage_fastqs.py`` can be retrieved
in batch using the ``download_fastqs.py`` utility::

     download_fastqs.py http://example.com/fastqs/

This fetches the checksum file from the URL and then uses that to get a
list of Fastq files to download. Once the files are downloaded it runs
the Linux ``md5sum`` program to verify the integrity of the downloads.

.. note::

   This utility is stand-alone so it can be sent to end users and
   used independently of other components of the autoprocess package.

.. _update_project_metadata:

``update_project_metadata.py``: manage metadata associated with a project
*************************************************************************

The projects within a run each have a file called ``README.info`` which is
used to hold metadata about that project (for example, user, PI, organism,
library type and so on).

Use the ``update_project_metadata.py`` utility to check and update the
metadata associated with a project, for example to update the PI::

    update_project_metadata.py ANALYSIS_DIR PROJECT -u PI="Andrew Jones"

.. note::

    Project directories created using very old versions of ``auto_process``,
    or predating the automated processing system, might not have metadata
    files. To create one use::

        update_project_metadata.py ANALYSIS_DIR PROJECT -i

    before using ``-u`` to populate the fields.

.. _audit_projects:

``audit_projects.py``: auditing disk usage for multiple runs
************************************************************

Collections of runs that are copied to an 'archive' location via the
``archive`` function of ``auto_process.py`` will form a directory structure
of the form::

    ARCHIVE_DIR/
      |
      +--- 2015/
            |
            +--- hiseq/
                  |
                  +--- 150429_HISEQ_XXYYY_12345BB_analysis/
                  |
                  +--- 150408_HISEQ_XXYYY_67890CC_analysis/
                  |
                  .

Within each run dir there will be one or more project directories.

The projects can be audited according to PI and disk usage using the
``audit_projects.py`` utility, for example::

    audit_projects.py ARCHIVE_DIR/2015/hiseq/

Multiple directories can be specified, e.g.::

    audit_projects.py ARCHIVE_DIR/2015/hiseq/ ARCHIVE_DIR/2014/hiseq/

This will print out a summary of usage for each PI, e.g.::

    Summary (PI, # of projects, total usage):
    =========================================
    Peter Brooks	12	3.7T
    Trevor Smith	8	2.3T
    Donald Raymond	6	2.2T
    ...
    Total usage	164	22.3T

plus a breakdown of the usage for each of the projects belonging to each
PI, for example::

    Breakdown by PI/project:
    ========================
    Peter Brooks:
	150121_HISEQ001_0123_ABCD123XX:	SteveAustin	128.1G
	150306_HISEQ001_0234_ABCD123XX:	MartinLouis	159.7G
	150415_HISEQ001_0345_ABCD123XX:	MartinLouis	72.8G
        ...

There is also a summary of the amount of space used for storing the
'undetermined' read data, for each run.

.. note::

   The disk usage for each file is calculated by using Python's ``os.lstat``
   function to get the number of 512-byte blocks per file. The total usage
   is then the sum of all the files and directories.

   However these values can differ from the sizes returned by the Linux
   ``du`` program, for various reasons including using a different block
   size (e.g. ``du`` uses 1024-byte blocks). So the returned values should
   not be treated as absolutes.

.. _exporting_to_galaxy:

Exporting Fastqs to a data library in a local Galaxy instance
*************************************************************

Upload of Fastq files from a run into a data library on a Galaxy instance
can be performed using the ``nebulizer`` utility.

.. note::

   You will need access to an admin account on the target Galaxy
   server to create and add to the data libraries.

The ``create_library`` and ``create_library_folder`` commands can be used
to make the target data library and folder, if these don't already exist -
for example:

::

    nebulizer create_library MyGalaxy "MISEQ_190626#26" \
        --description "Data from MISEQ run 26 datestamp 190626"
    nebulizer create_library_folder MyGalaxy "MISEQ_190626#26/Fastqs"

would create a data library called *MISEQ_190626#26* on the *MyGalaxy*
instance, and a new folder called *Fastqs* within that library.

Then the ``add_library_datasets`` command can be used to upload Fastqs
to the library.

To upload files from the local system to the server:

::

    nebulizer add_library_datasets MyGalaxy /path/to/fastqs/PB_S1_R1_001.fastq.gz ...

If the files are on the same system as the Galaxy server then the
``--server`` option can be used, for example:

::

    nebulizer add_library_datasets mygalaxy --server Data_Library/Fastqs /path/to/fastqs/on/server/PB_S1_R1_001.fastq.gz ...

It is possible in this case to get Galaxy to create links to the Fastqs
(rather than making copies) which can potentially save time and disk
space, by including the ``--link`` option:

::

    nebulizer add_library_datasets mygalaxy --server --link Data_Library/Fastqs /path/to/fastqs/on/server/PB_S1_R1_001.fastq

.. warning::

   Making links only seems to work for uncompressed Fastq files.

For information on ``nebulizer`` see
https://nebulizer.readthedocs.io/en/latest/
