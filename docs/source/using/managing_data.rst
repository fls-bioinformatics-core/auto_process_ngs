Managing and sharing data
=========================

These additional tasks are not really part of the automated processing, but
the utilities for performing them are currently part of the same package so
they are outlined here.

 * :ref:`manage_fastqs`
 * :ref:`download_fastqs`
 * :ref:`export_to_galaxy`
 * :ref:`update_project_metadata`
 * :ref:`audit_projects`

.. _manage_fastqs:

``manage_fastqs.py``: copying data for transfer to end users
************************************************************

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

    manage_fastqs.py RUN_DIR PROJECTNAME copy /path/to/local/dir --filter *_R1_*

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

.. _export_to_galaxy:

``export_to_galaxy.py``: export Fastqs to a data library in a local Galaxy instance
***********************************************************************************

The ``export_to_galaxy.py`` utility will upload Fastq files from a run into
a data library on a local Galaxy instance::

    export_to_galaxy.py http://localhost:8080 "MISeq Data" /mnt/galaxy-data ANALYSIS_DIR

The script performs the following actions:

 * Copy uncompressed versions of the Fastq files from ``ANALYSIS_DIR`` to a
   ``RUN/PROJECT`` directory structure under the target directory (e.g.
   ``/mnt/galaxy-data``)
 * Make a new folder for the run in the target data library (e.g.
   ``MISeq Data``) with subfolders for each project linked to the
   uncompressed Fastq files created in the previous step.

.. note::

   The destination directory must be mounted locally from both where the
   ``export_to_galaxy.py`` utility is being run, and the Galaxy server
   itself.

.. note::

   The target data library must already exist on the Galaxy server.

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
