============
Requirements
============

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
run_qc              `multiqc`_
process_icell8      `cutadapt`_
process_icell8      `fastq_screen`_
process_icell8      `bowtie2`_         Required by fastq_screen
process_10xgenomics `cellranger`_
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
* :ref:`auto_process_reference_data_10xgenomics`
  
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
  
.. _auto_process_reference_data_10xgenomics:

----------------------------------------------------------------
process_10xgenomics (single library analysis for scRNA-seq data)
----------------------------------------------------------------

The single library analysis step of ``process_10xgenomics`` for
single cell RNA-seq wraps ``cellranger count`` and requires a
compatible ``cellranger`` transcriptome reference data set for the
organism in question to be provided.

10xGenomics provide a number of reference data sets which can
be downloaded via:

* https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation

(There are also instructions for constructing reference data
for novel organisms that are not supported.)

These can be defined in the ``10xgenomics_transcriptomes``
section of the ``auto_process.ini`` file, for example::

  [10xgenomics_transcriptomes]
  human = /data/cellranger/refdata-cellranger-GRCh38-1.2.0
  mouse = /data/cellranger/refdata-cellranger-mm10-1.2.0

or else must be specified using the relevant command line
option.

.. _auto_process_reference_data_10xgenomics_snrna_seq:

----------------------------------------------------------------
process_10xgenomics (single library analysis for snRNA-seq data)
----------------------------------------------------------------

When dealing with single-nuclei RNA-seq (snRNA-seq) 10xGenomics
data, it is recommended that ``cellranger count`` is run with a
compatible ``cellranger`` "pre-mRNA" reference package (which
includes both intronic and exonic information) instead of the
standard transcriptome reference used for scRNA-seq.

10xGenomics don't provide pre-mRNA references, but the
documentation explains how to generate a custom reference
package for these data:

* https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references#premrna

These can be made available within ``auto-process`` by adding
definitions into the ``10xgenomics_premrna_references``
section of the ``auto_process.ini`` file, for example::

  [10xgenomics_premrna_references]
  human = /data/cellranger/refdata-cellranger-GRCh38-1.2.0_premrna
  mouse = /data/cellranger/refdata-cellranger-mm10-1.2.0_premrna

.. _auto_process_reference_data_10xgenomics_atac:

-----------------------------------------------------------------
process_10xgenomics (single library analysis for scATAC-seq data)
-----------------------------------------------------------------

The single library analysis step of ``process_10xgenomics`` for
single cell ATAC-seq data wraps ``cellranger-atac count`` and
requires a compatible ``cellranger-atac`` ATAC genome reference
data set for the organism in question to be provided.

10xGenomics provide a number of reference data sets which can
be downloaded via:

* https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/installation

(There are also instructions for constructing reference data
for novel organisms that are not supported.)

These can be defined in the ``10xgenomics_atac_genome_references``
section of the ``auto_process.ini`` file, for example::

  [10xgenomics_atac_genome_references]
  human = /data/cellranger/refdata-cellranger-atac-GRCh38-1.0.1
  mouse = /data/cellranger/refdata-cellranger-atac-mm10-1.0.1

or else must be specified using the relevant command line
option.
