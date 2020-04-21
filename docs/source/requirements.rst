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
run_qc              `fastqc`_
run_qc              `fastq_screen`_
run_qc              `bowtie`_          Required by fastq_screen
run_qc              `STAR`_            Required for strandedness determination
run_qc              `cellranger`_      10xGenomics Chromium single-cell RNA-seq data only
run_qc              `cellranger-atac`_ 10xGenomics Chromium single-cell ATAC-seq data only
run_qc              `multiqc`_
process_icell8      `cutadapt`_
process_icell8      `fastq_screen`_
process_icell8      `bowtie2`_         Required by fastq_screen
=================== ================== ===================

.. _bcl2fastq 2.17: https://support.illumina.com/downloads/bcl2fastq-conversion-software-v217.html
.. _bcl2fastq1.8.4: http://support.illumina.com/downloads/bcl2fastq_conversion_software_184.html
.. _cellranger: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger
.. _cellranger-atac: https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/what-is-cell-ranger-atac
.. _fastqc:  http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
.. _fastq_screen: http://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/
.. _bowtie: http://bowtie-bio.sourceforge.net/index.shtml
.. _bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
.. _STAR: https://github.com/alexdobin/STAR
.. _multiqc: http://multiqc.info/
.. _cutadapt: http://cutadapt.readthedocs.io

These programs provided by these packages must be found on the
``PATH`` when the appropriate autoprocessor commands are issued.
:ref:`environment-modules` can be used to help manage this.
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

In addition the strandedness determination requires ``STAR``
indexes for each organism of interest. These can then be
defined in the ``fastq_strand_indexes`` section of the
``auto_process.ini`` file, for example::

  [fastq_strand_indexes]
  human = /data/genomeIndexes/hg38/STAR/
  mouse = /data/genomeIndexes/mm10/STAR/

For 10xGenomics single cell data, the single library analysis
requires appropriate compatible reference data for
``cellranger[-atac] count``:

* **scRNA-seq data**: transcriptome reference data set
* **snRNA-seq data**: "pre-mRNA" reference data set (which
  includes both intronic and exonic information)
* **sc/snATAC-seq**: 

The reference data sets can be assigned to different organisms
in the ``10xgenomics_transcriptomes``,
``10xgenomics_premrna_references`` and
``10xgenomics_atac_genome_references``
sections of the ``auto_process.ini`` file.

For example:

::

   [10xgenomics_transcriptomes]
   human = /data/cellranger/refdata-cellranger-GRCh38-1.2.0
   mouse = /data/cellranger/refdata-cellranger-mm10-1.2.0
   
   [10xgenomics_premrna_references]
   human = /data/cellranger/refdata-cellranger-GRCh38-1.2.0_premrna
   mouse = /data/cellranger/refdata-cellranger-mm10-1.2.0_premrnaferences``

   [10xgenomics_atac_genome_references]
   human = /data/cellranger/refdata-cellranger-atac-GRCh38-1.0.1
   mouse = /data/cellranger/refdata-cellranger-atac-mm10-1.0.1

Alternatively reference data sets can be specified at run-time
using the ``--10x_transcriptome`` and ``--10x_premrna_reference``
command line options of ``run_qc`` and the ``run_qc.py`` utility.

10xGenomics provide a number of reference data sets for scRNA-seq
and ATAC-seq data, which can be downloaded via:

* https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation
* https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/installation

There are also instructions for constructing reference data for
novel organisms that are not supported by 10xGenomics.

Pre-mRNA references are currently not available, but the documentation
explains how to generate a custom reference package for these data:

* https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references#premrna
  
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
