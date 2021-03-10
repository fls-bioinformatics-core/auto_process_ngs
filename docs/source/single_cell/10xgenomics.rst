Processing 10xGenomics single cell data
=======================================

Background
----------

10xGenomics provides a range of platforms for various types of single
cell samples including gene expression, ATAC-seq and spatial
transcriptomics.

The ``auto_process`` pipeline has integrated support for Fastq generation
and QC (in the ``make_fastqs`` and ``run_qc`` commands) for the following
platforms:

* Chromium SC 3'v2/v3 (single cell/single nuclei RNA-seq)
* Single cell/single nuclei ATAC-seq
* Visuim spatial transcriptomics
* Single cell multiome gene expression/ATAC

.. warning::

   The ``process_10xgenomics.py`` utility present in previous versions
   of ``auto_process`` is no longer supported and has now been removed
   from the package.

Fastq generation
----------------

The :doc:`make_fastqs <../using/make_fastqs>` command can be used
to generate Fastqs from sequencing data generated from 10xGenomics
samples, by specifying the appropriate protocol via the
``--protocol`` option.

See :ref:`make_fastqs-protocols` for more details.

Single library analyses
-----------------------

Single library analyses are performed automatically as part of the
QC pipeline, provided that the appropriate 10xGenomics software and
references are available.

See the documentation on the :doc:`run_qc command <../using/run_qc>`
and the :doc:`standalone run_qc.py utility <../using/run_qc_standalone>`
for more details.

Troubleshooting
---------------

Single-library analyses fail for low read counts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It has been observed that when the Fastq files produced by the ``mkfastq``
command have very low read counts then the single-library analyses may
fail, with ``cellranger count`` reporting an error of the form e.g.::

    Could not auto-detect Single Cell 3' chemistry. Fraction of barcodes
    on whitelist was at best 0.23%, while we expected at least 10.00% for
    one of the chemistries.

There is currently no workaround for this issue.

Single-library analyses fail to detect chemistry automatically
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default ``cellranger count`` attempts to determine the chemistry used
automatically, however this may fail if a low number of reads map to the
reference genome and give an error of the form::

    The chemistry was unable to be automatically determined. This can
    happen if not enough reads originate from the given reference. Please
    verify your choice of reference or explicitly specify the chemistry
    via the --chemistry argument.

If the reference data being used is correct then use the ``--chemistry``
option to specify the appropriate assay configuration - see
https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count
