============
Requirements
============

.. _supported_python_versions:

*************************
Supported Python versions
*************************

The package consists predominantly of code written in Python, and the
following versions are supported:

* Python 2.7
* Python 3.6
* Python 3.7
* Python 3.8

However as Python 2.7 is now end-of-life, support for this version with
be withdrawn in the near future.

.. _software_dependencies:

*********************
Software dependencies
*********************

The following ``auto_process`` subcommands depend on additional
third-party software packages which must be installed separately:

=================== ================== ===================
Pipeline stage      Software packages  Notes
=================== ================== ===================
make_fastqs         `bcl2fastq 2.17`_  2.17+ recommended
make_fastqs         `cellranger`_      10xGenomics Chromium single-cell RNA-seq data only
make_fastqs         `cellranger-atac`_ 10xGenomics Chromium single-cell ATAC-seq data only
make_fastqs         `cellranger-arc`_  10xGenomics Multiome ATAC + GEX data
make_fastqs         `spaceranger`_     10xGenomics Visium spatial RNA-seq data only
run_qc (*)          `fastqc`_
run_qc (*)          `fastq_screen`_
run_qc (*)          `bowtie`_          Required by fastq_screen
run_qc (*)          `STAR`_            Required for strandedness determination
run_qc              `cellranger`_      10xGenomics Chromium single-cell RNA-seq data only
run_qc              `cellranger-atac`_ 10xGenomics Chromium single-cell ATAC-seq data only
run_qc (*)          `multiqc`_
process_icell8      `cutadapt`_
process_icell8      `fastq_screen`_
process_icell8      `bowtie2`_         Required by fastq_screen
=================== ================== ===================

.. _bcl2fastq 2.17: https://support.illumina.com/downloads/bcl2fastq-conversion-software-v217.html
.. _bcl2fastq1.8.4: http://support.illumina.com/downloads/bcl2fastq_conversion_software_184.html
.. _cellranger: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger
.. _cellranger-atac: https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/what-is-cell-ranger-atac
.. _cellranger-arc: https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/what-is-cell-ranger-arc
.. _spaceranger: https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/what-is-space-ranger
.. _fastqc:  http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
.. _fastq_screen: http://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/
.. _bowtie: http://bowtie-bio.sourceforge.net/index.shtml
.. _bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
.. _STAR: https://github.com/alexdobin/STAR
.. _multiqc: http://multiqc.info/
.. _cutadapt: http://cutadapt.readthedocs.io

(*) indicates packages that only need to be installed if
:ref:`conda_dependency_resolution` hasn't been enabled in the
configuration (or by an appropriate command line option); otherwise
the programs provided by these packages must be available on the
``PATH`` when the appropriate autoprocessor commands are issued.
:ref:`environment_modules` can be used to help manage this.
Alternatively many of these packages can be obtained from the
`bioconda project <https://bioconda.github.io/>`_.

..  note::

    Fastq generation requires Illumina's ``bcl2fastq`` software.
    The recommended version is 2.17+ (earlier versions should work
    but note that they cannot handle NextSeq data); if there are
    multiple ``bcl2fastq`` packages available on the path at run
    time then see :ref:`required_bcl2fastq_versions` for how to
    specify which version is used.

.. _reference_data:

**************
Reference data
**************

The following ``auto_process`` stages require additional reference
data:

* :ref:`auto_process_reference_data_run_qc`
* :ref:`auto_process_reference_data_icell8`
  
.. _auto_process_reference_data_run_qc:

------
run_qc
------

The QC pipeline uses the ``illumina_qc.sh`` script from
`genomics-bcftbx <https://genomics-bcftbx.readthedocs.io/>`_,
which requires a set of ``fastq_screen`` conf files and
underlying ``bowtie`` indexes to be created - these are
described here:

* https://genomics-bcftbx.readthedocs.io/en/latest/config.html#set-up-reference-data

Other reference data are required for the calculation of
various QC metrics:

* Strandedness determination requires ``STAR`` indexes for
  each organism of interest;
* Single library analyses of 10xGenomics single cell data
  require the appropriate compatible reference datasets for
  ``cellranger[-atac|-arc] count``:
  - **scRNA-seq data**: transcriptome reference data set
  - **snRNA-seq data**: "pre-mRNA" reference data set (which
    includes both intronic and exonic information)
  - **sc/snATAC-seq**: Cell Ranger ATAC compatible genome
    reference
  - **single cell multiome GEX+ATAC data**: ``cellranger-arc``
    compatible reference package

These can be defined in ``[organism:...]`` sections of the
``auto_process.ini`` file, for example:

::

   [organism: human]
   star_index = /data/genomeIndexes/hg38/STAR/
   cellranger_reference = /data/10x/refdata-cellranger-GRCh38-1.2.0
   cellranger_premrna_reference = /data/10x/refdata-cellranger-GRCh38-1.2.0_premrna
   cellranger_atac_reference = /data/10x/refdata-cellranger-atac-GRCh38-1.0.1
   cellranger_arc_reference = /data/10x/refdata-cellranger-arc-GRCh38-2020-A
   
   [organism: mouse]
   star_index = /data/genomeIndexes/mm10/STAR/
   cellranger_reference = /data/10x/refdata-cellranger-mm10-1.2.0
   cellranger_atac_reference = /data/10x/refdata-cellranger-atac-mm10-1.0.1
   cellranger_arc_reference = /data/10x/refdata-cellranger-arc-mm10-2020-A

.. note::

   Alternatively reference data sets can be specified at run-time
   for single cell and single nuclei RNA-seq using the
   ``--10x_transcriptome`` and ``--10x_premrna_reference``
   command line options respectively with the ``run_qc`` command
   and the ``run_qc.py`` utility.

10xGenomics provide a number of reference data sets for scRNA-seq,
ATAC-seq and single cell multiome data, which can be downloaded via:

* https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation
* https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/installation
* https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/installation

There are also instructions for constructing reference data for
novel organisms that are not supported by 10xGenomics.

Pre-mRNA references are currently not available, but the documentation
explains how to generate a custom reference package for these data:

* https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references#premrna

.. note::

   The ``[organism:...]`` sections supersede the old
   ``fastq_strand_indexes`` and ``10xgenomics...`` sections
   of the ``auto_process.ini`` file; the old sections are
   still recognised for now but are deprecated and likely to
   be dropped in future.
  
.. _auto_process_reference_data_icell8:

--------------------------------------
process_icell8 (contaminant filtering)
--------------------------------------

The contaminant filtering stage of ``process_icell8`` needs
two ``fastq_screen`` conf files to be set up, one containing
``bowtie`` indexes for "mammalian" genomes (typically human
and mouse) and another containing indexes for "contaminant"
genomes (yeast, E.coli, UniVec7, PhiX, mycoplasma, and
adapter sequences).

These can be defined in the ``icell8`` section of the
``auto_process.ini`` file, for example::

  [icell8]
  mammalian_conf_file = /data/icell8/mammalian_genomes.conf
  contaminants_conf_file = /data/icell8/contaminant_genomes.conf

or else must be specified using the relevant command line
options.
