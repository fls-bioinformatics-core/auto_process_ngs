============
Requirements
============

*********************
Software dependencies
*********************

Many of the functions of ``auto-process-ngs`` depend on additional
third-party software packages which currenty must be installed
separately.

The programs provided by these packages must be found on the
``PATH`` when the appropriate autoprocessor commands are issued.
:ref:`environment-modules` can be used to help manage this.
Alternatively many of these packages can be obtained from the
`bioconda project <https://bioconda.github.io/>`_.

-----------
make_fastqs
-----------

Fastq generation requires Illumina's ``bcl2fastq`` software.
The recommended version is 2.17+ (earlier versions should work
but note that they cannot handle NextSeq data):

=============== ==============================================================================
bcl2fastq 2.17  https://support.illumina.com/downloads/bcl2fastq-conversion-software-v217.html
bcl2fastq 1.8.4 http://support.illumina.com/downloads/bcl2fastq_conversion_software_184.html
=============== ==============================================================================

10xGenomics Chromium single-cell data also requires the
``cellranger`` pipeline (latest version recommended):

========== =========================================================================================================
cellranger https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger
========== =========================================================================================================

..  note::

    If there are multiple ``bcl2fastq`` packages available on the path
    at run time then see :ref:`required_bcl2fastq_versions` for how to
    specify which version is used.

------
run_qc
------

The standard QC pipeline requires the following external
software:

============ ===============================================================
bowtie       http://bowtie-bio.sourceforge.net/index.shtml
fastq_screen http://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/
fastqc       http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
multiqc      http://multiqc.info/
============ ===============================================================

-----------------
process_icell8.py
-----------------

The ICELL8 processing pipeline requires the following software
(in addition to the packages specified for ``run_qc``:

* ``cutadapt`` http://cutadapt.readthedocs.io
* ``bowtie2`` (optional, required if using Bowtie2 genome indexes
  for contaminant filtering - see :doc:`icell8`)
  http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

*****************************
Supported sequencer platforms
*****************************

The pipeline is currently used for output from the following
Illumina sequencers:

* HISeq 4000
* MISeq
* NextSeq
* MiniSeq

Earlier versions have been used on GAIIx and HISeq 2000/2500.

*******************************
Supported single-cell platforms
*******************************

The pipeline supports handling data from the Takara Bio SMARTer
ICELL8 and 10xGenomics Chromium single-call RNA-seq platforms:

* :doc:`Handling ICELL8 scRNA-seq data <icell8>`
* :doc:`Handling 10xGenomics Chromium scRNA-seq data <10xgenomics>`
