============
Requirements
============

.. _supported_python_versions:

*************************
Supported Python versions
*************************

The package consists predominantly of code written in Python, and the
following versions are supported:

* Python 3.6
* Python 3.7
* Python 3.8

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
make_fastqs         `bcl-convert`_     Alternative to ``bcl2fastq``
make_fastqs         `cellranger`_      10xGenomics Chromium single-cell RNA-seq data only
make_fastqs         `cellranger-atac`_ 10xGenomics Chromium single-cell ATAC-seq data only
make_fastqs         `cellranger-arc`_  10xGenomics Multiome ATAC + GEX data
make_fastqs         `spaceranger`_     10xGenomics Visium spatial RNA-seq data only
run_qc (*)          `fastqc`_
run_qc (*)          `fastq_screen`_
run_qc (*)          `bowtie`_          Required by fastq_screen
run_qc (*)          `STAR`_            Required for strandedness and alignment
run_qc (*)          `picard`_          Required for insert size metrics
run_qc (*)          `rseqc`_           Required for gene body coverage
run_qc (*)          `qualimap`_        Required for per-Fastq genomic origin of reads etc
run_qc              `cellranger`_      10xGenomics Chromium single-cell RNA-seq data only
run_qc              `cellranger-atac`_ 10xGenomics single-cell ATAC-seq data only
run_qc              `cellranger-arc`_  10xGenomics Multiome ATAC + GEX data
run_qc (*)          `multiqc`_
process_icell8      `cutadapt`_
process_icell8      `fastq_screen`_
process_icell8      `bowtie2`_         Required by fastq_screen
=================== ================== ===================

.. _bcl2fastq 2.17: https://support.illumina.com/downloads/bcl2fastq-conversion-software-v217.html
.. _bcl2fastq1.8.4: http://support.illumina.com/downloads/bcl2fastq_conversion_software_184.html
.. _bcl-convert: https://support.illumina.com/sequencing/sequencing_software/bcl-convert.html
.. _cellranger: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger
.. _cellranger-atac: https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/what-is-cell-ranger-atac
.. _cellranger-arc: https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/what-is-cell-ranger-arc
.. _spaceranger: https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/what-is-space-ranger
.. _fastqc:  http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
.. _fastq_screen: http://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/
.. _bowtie: http://bowtie-bio.sourceforge.net/index.shtml
.. _bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
.. _STAR: https://github.com/alexdobin/STAR
.. _picard: https://gatk.broadinstitute.org/hc/en-us/articles/360037055772-CollectInsertSizeMetrics-Picard-
.. _rseqc: http://rseqc.sourceforge.net/#
.. _qualimap: http://qualimap.conesalab.org/doc_html/command_line.html#rna-seq-qc
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

Reference data required for the calculation of various QC metrics
are described in the sections below, along with the configuration
settings needed to make them available to the QC pipeline.

FastqScreen
^^^^^^^^^^^

FastqScreen requires one or more ``conf`` files (each of which
defines a specific "screen") along with the underlying ``bowtie``
indexes for each organism which are included in the screens.

Indexes can be created manually, or by using the ``build_index.py``
utility (see :ref:`build_indexes`); the ``conf`` files must be
created manually (see the
`FastqScreen documentation <https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/_build/html/index.html#configuration>`_);
each screen can then be added to the configuration using a
``screen`` section which references the corresponding ``conf``
file, e.g.:

::

   [screen:model_organisms]
   conf_file = /data/model_organisms.conf

The screens to use in the pipeline must be set using the
``fastq_screens`` parameter in the ``qc`` section, e.g.:

::

   [qc]
   fastq_screens = model_organisms,other_organisms,rRNA
   ...

.. note::

   This replaces the old ``qc.setup`` script that was used
   to define the location of a set of standard screen ``conf``
   files, used in earlier versions of the pipeline. Note
   that ``qc.setup`` is not longer needed (and will be ignored
   if present).

Strandedness
^^^^^^^^^^^^

Strandedness determination requires ``STAR`` indexes for each
organism of interest. These can be defined using appropriate
settings in ``[organism:...]`` sections of the ``auto_process.ini``
file, for example:

::

   [organism: human]
   star_index = /data/genomeIndexes/hg38/STAR/

   [organism: mouse]
   star_index = /data/genomeIndexes/mm10/STAR/

Indexes can be created manually, or by using the
``build_index.py`` utility (see :ref:`build_indexes`).

.. note::

   The ``[organism:...]`` sections supersede the old
   ``fastq_strand_indexes`` section of the ``auto_process.ini``
   file; the older section is still recognised for now but is
   deprecated and likely to be dropped in future.

Insert size metrics (Picard)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Picard's ``CollectInsertSizeMetrics`` needs a STAR index for
each organism of interest (in order to generate a BAM file from
the sequences). This should be specfied in the ``[organism:...]``
sections of the ``auto_process.ini`` configuration file, for example:

::

   [organism: human]
   star_index = /data/genomeIndexes/hg38/STAR/

STAR indexes can be created manually, or by using the
``build_index.py`` utility (see :ref:`build_indexes`).

RSeQC gene body coverage
^^^^^^^^^^^^^^^^^^^^^^^^

RSeQC ``geneBody_coverage.py`` needs both a STAR index (in order
to generate a BAM file from the sequences) and gene annotation in
BED format, for each organism of interest. These should be specfied
in the ``[organism:...]`` sections of the ``auto_process.ini``
configuration file, for example:

::

   [organism: human]
   star_index = /data/genomeIndexes/hg38/STAR/
   annotation_bed = /data/genomeIndexes/hg38/hg38.HouseKeepingGenes.bed

.. note::

   STAR indexes can be created manually, or by using the
   ``build_index.py`` utility (see :ref:`build_indexes`). Suitable
   gene model files for human and mouse can be downloaded from
   the RSeQC webpages at
   http://rseqc.sourceforge.net/#download-gene-models-update-on-12-14-2021


Qualimap RNA-seq metrics
^^^^^^^^^^^^^^^^^^^^^^^^

Qualimap's ``rnaseq`` command a STAR index (in order to generate a BAM
file from the sequences) and gene annotation in GTF format, for each
organism of interest. The pipeline also requires annotation in BED
format, in order to run RSeQC's ``infer_experiment.py`` command to
determine strand specificity (which is needed as input to Qualimap).

All these should be specfied in the ``[organism:...]`` sections of the
``auto_process.ini`` configuration file, for example:

::

   [organism: human]
   star_index = /data/genomeIndexes/hg38/STAR/
   annotation_bed = /data/genomeIndexes/hg38/hg38.HouseKeepingGenes.bed
   annotation_gtf = /data/genomeIndexes/hg38/gencode.v40.annotation.gtf

STAR indexes can be created manually, or by using the ``build_index.py``
utility (see :ref:`build_indexes`).

Single cell analyses
^^^^^^^^^^^^^^^^^^^^

Single library analyses of 10xGenomics single cell data require
the appropriate compatible reference datasets for
``cellranger[-atac|-arc] count``:

* **scRNA-seq data**: transcriptome reference data set
* **snRNA-seq data**: "pre-mRNA" reference data set (which
  includes both intronic and exonic information)
* **sc/snATAC-seq**: Cell Ranger ATAC compatible genome
  reference
* **single cell multiome GEX+ATAC data**: ``cellranger-arc``
  compatible reference package

These can all be defined using appropriate settings in
``[organism:...]`` sections of the ``auto_process.ini`` file,
for example:

::

   [organism: human]
   cellranger_reference = /data/10x/refdata-cellranger-GRCh38-1.2.0
   cellranger_premrna_reference = /data/10x/refdata-cellranger-GRCh38-1.2.0_premrna
   cellranger_atac_reference = /data/10x/refdata-cellranger-atac-GRCh38-1.0.1
   cellranger_arc_reference = /data/10x/refdata-cellranger-arc-GRCh38-2020-A

   [organism: mouse]
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
   ``10xgenomics...`` sections of the ``auto_process.ini`` file;
   the old sections are still recognised for now but are
   deprecated and likely to be dropped in future.

Annotation data
^^^^^^^^^^^^^^^

Annotation data in BED and GTF formats can be specified for
organisms of interest via the ``annotation_bed`` and ``annotation_gtf``
settings respectively in ``[organism:...]`` sections of the
``auto_process.ini`` file.

For example:

::

   [organism: human]
   annotation_bed = /data/genomeIndexes/hg38/annotation/hg38_NCBI_RefSeq_All.bed
   annotation_gtf = /data/genomeIndexes/hg38/annotation/hg38_NCBI_RefSeq_All.gtf

   [organism: mouse]
   annotation_bed = /data/genomeIndexes/mm10/annotation/gencode.vM25.annotation.bed
   annotation_gtf = /data/genomeIndexes/mm10/annotation/gencode.vM25.annotation.gtf
  
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

.. _build_indexes:

*****************************
Building indexes for aligners
*****************************

The :ref:`_utilities_build_index.py` utility can be used to
build indexes for ``bowtie``, ``bowtie2`` and ``STAR`` from
the appropriate data files (which must be obtained
separately).

For example: to build indexes for ``hg38`` using STAR version
2.7.7a:

::

   build_index.py star -V 2.7.7a \
       -o hg38_STAR_2.7.7a_gencode40 \
       /mnt/genome_data/hg38/hg38.fa \
       /mnt/genome_data/hg38/hg38.gencode.v40.annotation.gtf

.. note::

   If :ref:`conda_dependency_resolution` isn't enabled then
   the required aligner must be accessible on the ``PATH``,
   and the requested aligner version will be ignored.
