============
Requirements
============

*********************
Software dependencies
*********************

Many of the functions of ``auto_process`` depend on additional
third-party software packages which must be installed separately:

=================== ================= ===================
Pipeline stage      Software packages Notes
=================== ================= ===================
make_fastqs         `bcl2fastq 2.17`_ 2.17+ recommended
make_fastqs         `cellranger`_     10xGenomics Chromium single-cell data only
run_qc              `fastqc`_
run_qc              `fastq_screen`_
run_qc              `bowtie`_         Required by fastq_screen
run_qc              `STAR`_           Required for strandedness determination
run_qc              `multiqc`_
process_icell8      `cutadapt`_
process_icell8      `fastq_screen`_
process_icell8      `bowtie2`_        Required by fastq_screen
process_10xgenomics `cellranger`_
=================== ================= ===================

.. _bcl2fastq 2.17: https://support.illumina.com/downloads/bcl2fastq-conversion-software-v217.html
.. _bcl2fastq1.8.4: http://support.illumina.com/downloads/bcl2fastq_conversion_software_184.html
.. _cellranger: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger
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
