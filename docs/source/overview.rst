Overview
========

The ``auto-process-ngs`` package provides utilities for initial processing
of sequencing data produced by Illumina's range of sequencers, including
the HISeq, MISeq, NextSeq and MiniSeq platforms.

It also has functionality to deal specifically with data produced using
the ICELL8 and 10xGenomics Chromium single-cell platforms.

The package performs the following operations:

- Generation of Fastq files from the raw data produced by the sequencer
- Dividing Fastqs into 'projects' for subsequent analysis
- Performing initial QC on each project
